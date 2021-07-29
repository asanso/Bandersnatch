import os
from sage.all import *
from sage.rings.integer_ring import Z as ZZ
from sage.arith.misc import fundamental_discriminant
from sage.schemes.elliptic_curves.cm import hilbert_class_polynomial
from sage.schemes.elliptic_curves.constructor import EllipticCurve

f = open('./' + os.path.dirname(__file__) +'/small-disc-curves.txt', 'w')
p = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
assert ZZ(p).is_prime()
Fp = GF(p)
FpX = Fp['X']
X = FpX.gen()

f.write('p = {}\n\n\n'.format(hex(p)))

discs = []

for D in range(1,100):
    D = fundamental_discriminant(-D)
    if not(D in discs):
        discs.append(D)
        H = FpX(hilbert_class_polynomial(D))

        for (j,_) in H.roots():
            E = EllipticCurve_from_j(j)
            if not(E.is_supersingular()):
                print('D = {} = {}'.format(D, D.factor()))
                f.write('D = {} = {}\n\n'.format(D, D.factor()))
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
                for OrderE in PossibleOrders:
                    FactorsOrderE = OrderE.factor()
                    print('Curve order:\n {}'.format(FactorsOrderE))
                    f.write('Curve order:\n {}\n'.format(FactorsOrderE))
                    print('Curve security:\n {}-bit'.format(
                        FactorsOrderE[-1][0].nbits()//2
                    ))
                    f.write('Curve security:\n {}-bit\n\n'.format(
                        (FactorsOrderE[-1][0].nbits()//2)
                    ))
                print()
                f.write('\n')
            break
        
    
