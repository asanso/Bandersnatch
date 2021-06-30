import os
from sage.all import *
from sage.rings.integer_ring import Z as ZZ
from sage.arith.misc import fundamental_discriminant
from sage.schemes.elliptic_curves.cm import hilbert_class_polynomial
from sage.schemes.elliptic_curves.constructor import EllipticCurve

p = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
assert ZZ(p).is_prime()
Fp = GF(p)
FpX = Fp['X']
X = FpX.gen()


DD = []

for D in range(1,100):
    D = fundamental_discriminant(-D)
    if not(D in DD):
        DD.append(D)
DD.sort()
DD.reverse()

for D in DD:
    H = FpX(hilbert_class_polynomial(D))
    
    for (j,_) in H.roots():
        E = EllipticCurve_from_j(j)
        if not(E.is_supersingular()):
            t = E.trace_of_frobenius()
            PossibleOrders = [p+1-t, p+1+t]
            if j == 1728 :
                y = ((t**2 - 4 * p)//(-4)).sqrt()
                PossibleOrders.append(p+1+2*y)
                PossibleOrders.append(p+1-2*y)
            if j == 0 :
                y = ((t**2 - 4 * p)//(-3)).sqrt()
                PossibleOrders.append(p+1-(-3*y+t)/2)
                PossibleOrders.append(p+1-(-3*y-t)/2)
                PossibleOrders.append(p+1-(3*y-t)/2)
                PossibleOrders.append(p+1-(3*y+t)/2)
            foo=True
            for OrderE in PossibleOrders:
                if foo:
                    print('${}$'.format(D), end=' & ')
                    foo=False
                else:
                    print('', end=' & ')
                FactorsOrderE = OrderE.factor()
                print('${}$-bit'.format(
                    FactorsOrderE[-1][0].nbits()//2
                ), end=' & $')
                for (f,m) in FactorsOrderE[:-1]:
                    if m == 1:
                        print('{} '.format(f), end='')
                    else:
                        print(str(f) + '^{' + str(m) + '} ', end='')
                    print(' \cdot ', end='')
                # last one
                for (f,m) in [FactorsOrderE[-1]]:
                    if m == 1:
                        print('p_{' + str(f.nbits()) + '}', end='')
                    else:
                        print('p_{' + str(f.nbits()) + '}^{'+str(m)+'}', end='')
                print('$\\\\')
        break
