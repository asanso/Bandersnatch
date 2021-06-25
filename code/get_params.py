import os
from sage.all import *
from sage.schemes.elliptic_curves.cm import hilbert_class_polynomial
from sage.schemes.elliptic_curves.constructor import EllipticCurve

D = -8
H = hilbert_class_polynomial(D)
p = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
Fp = GF(p)
FpX = Fp['X']
X = FpX.gen()
j = FpX(H).roots()[0][0]

# Weierstrass curve: y² = x³ + A * x + b
a=20856223409030707214153620029242927260100058688354534505372746703669794261922
b=4118395616383396936705917269457205707251888519744679123386326773690169742488
E = EllipticCurve(Fp,[a,b])
assert E.j_invariant() == j
r = E.order()//4
assert r.is_prime()
assert (E.quadratic_twist().order()//(2**7*3**3)).is_prime()

Fr = GF(r)
L = -Fr(-2).sqrt()

alpha = E.division_polynomial(2).roots()[0][0]
P = E.lift_x(alpha)
phi0,phi1 = E.isogeny(P)
E2 = E.isogeny_codomain(P)
u = list(set((X**4 - E.a4()/E2.a4()).roots()) & set((X**6 - E.a6()/E2.a6()).roots()))[0][0]
assert u**4 == E.a4()/E2.a4()
assert u**6 == E.a6()/E2.a6()
u2 = u**2
u3 = u**3
r0,r1,_ = phi0.numerator()(X,1).list()
s0 = phi0.denominator()(X,1).list()[0]
t0,t1,_ = phi1.numerator()(X,1).list()

f = open(os.path.dirname(__file__) +'/params-W.py', 'w')
f.write('import sage.all\n')
f.write('from sage.rings.finite_rings.finite_field_constructor import FiniteField\n')
f.write('# Weierstrass parameters\n')
f.write('p={}\n'.format(p))
f.write('Fp=FiniteField(p)\n')
f.write('a=Fp({})\n'.format(a))
f.write('b=Fp({})\n'.format(b))
f.write('r0=Fp({})\n'.format(r0))
f.write('r1=Fp({})\n'.format(r1))
f.write('s0=Fp({})\n'.format(s0))
f.write('t0=Fp({})\n'.format(t0))
f.write('t1=Fp({})\n'.format(t1))
f.write('u2=Fp({})\n'.format(u2))
f.write('u3=Fp({})\n\n'.format(u3))
f.close()

# Montgomery curve: By² = x³ + A * x² + x
s = 1/(X**2-(3*alpha**2+a)).roots()[1][0]
A = 3*alpha*s
B = 5*s # we multiply by a non-square to be on the right twist
E_m = EllipticCurve([0,A,0,1,0]) # up to the twist
assert E_m.j_invariant() == j
assert E_m.quadratic_twist().order() == E.order()
x2 = E_m.division_polynomial(2).roots()[1][0]
# A+2*x2 == 2
c = A+2
u = 1/Fp(-2).sqrt()

f = open(os.path.dirname(__file__) +'/params-M.py', 'w')
f.write('import sage.all\n')
f.write('from sage.rings.finite_rings.finite_field_constructor import FiniteField\n')
f.write('# Montgomery parameters\n')
f.write('p={}\n'.format(p))
f.write('Fp=FiniteField(p)\n')
f.write('B=Fp({})\n'.format(B))
f.write('A=Fp({})\n'.format(A))
f.write('c=Fp({})\n'.format(c))
f.write('u=Fp({})\n'.format(u))
f.close()

# Twised Edwards curve: a * x² + y² = 1 + d * x² * y²
a = (A+2)/B
d = (A-2)/B
# We find nice a and d.
# We cannot put it in a=-1 form because it reaches the twist.
# a=-5 is possible.
d = -5*d/a
a = Fp(-5)
M = Matrix([[ZZ(d),1],[p,0]])
N = M.LLL()
d1,d2=N[0]
assert Fp(d1/d2) == d

'''
Twisted Edwards point (x,y)
maps to a Montgomery point (fx,fy)
fx = (1+y)/(1-y)
fy = fx/x.
The endomorphism map on the montgomery curve is
(psix,psiy) where
psix = -((xx-1)**2/xx + A +2) / Aplus2x2
psiy =  u*yy*(1-1/xx**2) / Aplus2x2
and it goes back to the Twisted Edwards point (ffx,ffy)
ffx = psix/psiy
ffy = (psix-1)/(psix+1).
'''

Fpxy = Fp['x,y']
G = FractionField(Fpxy)
x,y = G.gens()
# TE to M
fx = (1+y)/(1-y)
fy = fx/x
# psi
psix = -((fx-1)**2/fx + A +2) / 2
psiy =  u*fy*(1-1/fx**2) / 2
# M to TE
ffx = psix/psiy
ffy = (psix-1)/(psix+1)

fxn,fxd = ffx(1,X).numerator(), ffx(1,X).denominator()
fyn,fyd = ffy(1,X).numerator(), ffy(1,X).denominator()
a1 = fxn.leading_coefficient()
[(a2,_),(a3,_)] = fxn.roots()
b1 = a2
[(b2,_),(b3,_)] = fyn.roots()
[(c1,_),(c2,_)] = fyd.roots()

f = open(os.path.dirname(__file__) +'/params-TE.py', 'w')
f.write('import sage.all\n')
f.write('from sage.rings.finite_rings.finite_field_constructor import FiniteField\n')
f.write('# Twisted Edwards parameters\n')
f.write('p={}\n'.format(p))
f.write('Fp=FiniteField(p)\n')
f.write('# curve equation coefficients\n')
f.write('a=Fp({})\n'.format(a))
f.write('d_num=Fp({})\n'.format(d1))
f.write('d_den=Fp({})\n'.format(d2))
f.write('d=Fp({})\n'.format(d))
f.write('# endomorphism coefficients\n')
f.write('a1=Fp({})\n'.format(a1))
f.write('a2=Fp({})\n'.format(a2))
f.write('a3=Fp({})\n'.format(a3))
f.write('b1=Fp({})\n'.format(b1))
f.write('b2=Fp({})\n'.format(b2))
f.write('b3=Fp({})\n'.format(b3))
f.write('c1=Fp({})\n'.format(c1))
f.write('c2=Fp({})\n'.format(c2))
f.write('# endomorphism eigenvalue L on the subgroup of size r\n')
f.write('r={}\n'.format(r))
f.write('cofactor = 4\n')
f.write('curve_order=r*cofactor\n')
f.write('L={}\n'.format(L))
f.close()
