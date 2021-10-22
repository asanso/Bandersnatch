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
E = EllipticCurve(j=j)
if not((E.order()//4).is_prime()):
    # we want its twist, defined with sqrt5 (of Fp²)
    assert Fp(5).is_square() == False
    E=EllipticCurve([E.a4() * 5**2, E.a6() * 5**3])
a = E.a4()
b = E.a6()
r = E.order()//4
assert r.is_prime()
assert (E.quadratic_twist().order()//(2**7*3**3)).is_prime()

Fr = GF(r)
L = -Fr(-2).sqrt()

alpha = E.division_polynomial(2).roots()[0][0]
P = E.lift_x(alpha)
phi0,phi1 = E.isogeny(P)
E2 = E.isogeny_codomain(P)
# choice of u leads to psi or -psi
u = list(set((X**4 - E.a4()/E2.a4()).roots()) & set((X**6 - E.a6()/E2.a6()).roots()))[1][0]
assert u**4 == E.a4()/E2.a4()
assert u**6 == E.a6()/E2.a6()
u2 = u**2
u3 = u**3
r0,r1,_ = phi0.numerator()(X,1).list()
s0 = phi0.denominator()(X,1).list()[0]
t0,t1,_ = phi1.numerator()(X,1).list()

# generator
x_rand = 1
def next_small(x) :
    # y = -x
    # if y > 0 :
    #     y+=1
    # return y
    return x+1

while (x_rand**3 + a*x_rand + b).is_square() == False:
    x_rand = next_small(x_rand)
P_rand = E.lift_x(Fp(x_rand))
G = 4*P_rand
while G.is_zero():
    x_rand = next_small(x_rand)
    while (x_rand**3 + a*x_rand + b).is_square() == False:
        x_rand = next_small(x_rand)
    P_rand = E.lift_x(Fp(x_rand))
    G = 4*P_rand
assert not(G.is_zero()) and (r*G).is_zero()
xG_W = G[0]
yG_W = G[1]
assert G[2] == 1

#######
# check
psiG = E(u2*(xG_W**2 + 44800*xG_W + 2257920000)/(xG_W + 44800),u3*yG_W*(xG_W**2+2*44800*xG_W + t0)/(xG_W+44800)**2)
assert psiG == int(L)*G
#######

f = open('./' + os.path.dirname(__file__) +'/params-W.py', 'w')
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
f.write('xG_W=Fp(0x{})\n'.format(ZZ(xG_W).hex()))
f.write('yG_W=Fp(0x{})\n'.format(ZZ(yG_W).hex()))
f.close()

f = open('./' + os.path.dirname(__file__) +'/../paper/params-W.tex', 'w')
f.write('\\begin{verbatim}\n')
f.write('xW={}\n'.format(ZZ(xG_W).hex()))
f.write('yW={}.\n'.format(ZZ(yG_W).hex()))
f.write('\\end{verbatim}\n')
f.close()

# Montgomery curve: By² = x³ + A * x² + x
s = 1/(X**2-(3*alpha**2+a)).roots()[1][0]
A = 3*alpha*s
B = s
E_m = EllipticCurve([0,A,0,1,0]) # up to the twist
assert E_m.j_invariant() == j
assert E_m.quadratic_twist().order() == E.order()
# because here it is up to the twist
x2 = E_m.division_polynomial(2).roots()[1][0]
# A+2*x2 == 2
c = A+2
u = 1/Fp(-2).sqrt()
xG_M = s * (xG_W - alpha)
yG_M = s * yG_W
assert B*yG_M**2 == xG_M**3 + A*xG_M**2 + xG_M

f = open('./' + os.path.dirname(__file__) +'/params-M.py', 'w')
f.write('import sage.all\n')
f.write('from sage.rings.finite_rings.finite_field_constructor import FiniteField\n')
f.write('# Montgomery parameters\n')
f.write('p={}\n'.format(p))
f.write('Fp=FiniteField(p)\n')
f.write('B=Fp({})\n'.format(B))
f.write('A=Fp({})\n'.format(A))
f.write('c=Fp({})\n'.format(c))
f.write('u=Fp({})\n'.format(u))
f.write('xG_M=Fp(0x{})\n'.format(ZZ(xG_M).hex()))
f.write('yG_M=Fp(0x{})\n'.format(ZZ(yG_M).hex()))
f.write('zG_M=Fp(0x1)\n')
f.close()

f = open('./' + os.path.dirname(__file__) +'/../paper/params-M.tex', 'w')
f.write('\\begin{verbatim}\n')
f.write('xM={}.\n'.format(ZZ(xG_M).hex()))
f.write('\\end{verbatim}\n')
f.close()

# Twised Edwards curve: a * x² + y² = 1 + d * x² * y²
a_TE_0 = (A+2)/B
d_TE_0 = (A-2)/B
xG_TE_0 = xG_M / yG_M
yG_TE_0 = (xG_M-1) / (xG_M+1)
# We find nice a and d.
# We cannot put it in a=-1 form because it reaches the twist.
# a=-5 is possible.
d_TE = -5*d_TE_0/a_TE_0
a_TE = Fp(-5)
xG_TE = xG_TE_0 * (Fp(-a_TE_0/5).sqrt())
yG_TE = yG_TE_0
assert a_TE*xG_TE**2 + yG_TE**2 == 1 + d_TE* xG_TE**2 * yG_TE**2
M = Matrix([[ZZ(d_TE),1],[p,0]])
N = M.LLL()
d1,d2=N[0]
assert Fp(d1/d2) == d_TE

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
a = fxn.list()
b = fyn.list()
c = fyd.list()
assert a[1] == 0 and b[1] == 0 and c[1] == 0 and b[2] == -c[0]

f = open('./' + os.path.dirname(__file__) +'/params-TE.py', 'w')
f.write('import sage.all\n')
f.write('from sage.rings.finite_rings.finite_field_constructor import FiniteField\n')
f.write('# Twisted Edwards parameters\n')
f.write('p={}\n'.format(p))
f.write('Fp=FiniteField(p)\n')
f.write('# curve equation coefficients\n')
f.write('a=Fp({})\n'.format(a_TE))
f.write('d_num=Fp({})\n'.format(d1))
f.write('d_den=Fp({})\n'.format(d2))
f.write('d=Fp({})\n'.format(d_TE))
f.write('# endomorphism coefficients\n')
f.write('b=Fp({})\n'.format(b[2]))
f.write('c=Fp({})\n'.format(a[0]/a_TE))
f.write('# endomorphism eigenvalue L on the subgroup of size r\n')
f.write('r={}\n'.format(r))
f.write('cofactor = 4\n')
f.write('curve_order=r*cofactor\n')
f.write('L={}\n'.format(L))
f.write('xG_TE=Fp(0x{})\n'.format(ZZ(xG_TE).hex()))
f.write('yG_TE=Fp(0x{})\n'.format(ZZ(yG_TE).hex()))
f.write('zG_TE=Fp(0x1)\n')
f.close()

f = open('./' + os.path.dirname(__file__) +'/../paper/params-TE.tex', 'w')
f.write('\\begin{verbatim}\n')
f.write('xTE={}\n'.format(ZZ(xG_TE).hex()))
f.write('yTE={}.\n'.format(ZZ(yG_TE).hex()))
f.write('\\end{verbatim}\n')
f.close()

