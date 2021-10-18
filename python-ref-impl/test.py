import random
from curve import Curve
from bandersnatch import Bandersnatch, BandersnatchPoint

###
# tests for the Curve class on the Jubjub curve
###

p = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
a = -1
d = 0x2a9318e74bfa2b48f5fd9207e6bd7fd4292d7f6d37579d2601065fd6d6343eb1
r = 0x0e7db4ea6533afa906673b0101343b00a6682093ccc81082d0970e5ed6f72cb7
E = Curve(p, a, d, r, 8)

P = E.random_point()
Q = E.random_point()

# add
R = P.add(Q)
assert R.on_curve()
assert R == Q.add(P)

# double
R = P.double()
assert R.on_curve()

# scalar mul
n = 5
R = P.scalar_mul(n)
assert R.on_curve()
assert P.scalar_mul(E.cofactor * E.r).is_zero()

# multi scalar mul
for [k1,k2] in [[12,13], [-12,13], [12,-13], [-12,-13]]:
    P1 = E.random_point()
    P2 = E.random_point()
    assert P1.scalar_mul(k1).add(P2.scalar_mul(k2)) == P1.multi_scalar_mul(k1,P2,k2)

# clear cofactor
assert P.clear_cofactor().scalar_mul(E.r).is_zero()


###
# Test for Bandersnatch
###

#
# code from the params-TE.py file
#
from sage.rings.finite_rings.finite_field_constructor import FiniteField
# Twisted Edwards parameters
p=52435875175126190479447740508185965837690552500527637822603658699938581184513
Fp=FiniteField(p)
# curve equation coefficients
a=Fp(52435875175126190479447740508185965837690552500527637822603658699938581184508)
d_num=Fp(138827208126141220649022263972958607803)
d_den=Fp(171449701953573178309673572579671231137)
d=Fp(45022363124591815672509500913686876175488063829319466900776701791074614335719)
# endomorphism coefficients
b=Fp(37446463827641770816307242315180085052603635617490163568005256780843403514036)
c=Fp(49199877423542878313146170939139662862850515542392585932876811575731455068989)
# endomorphism eigenvalue L on the subgroup of size r
r=13108968793781547619861935127046491459309155893440570251786403306729687672801
cofactor = 4
curve_order=r*cofactor
L=8913659658109529928382530854484400854125314752504019737736543920008458395397
xG_TE=Fp(0x29c132cc2c0b34c5743711777bbe42f32b79c022ad998465e1e71866a252ae18)
yG_TE=Fp(0x2a6c669eda123e0f157d8b50badcd586358cad81eee464605e3167b6cc974166)
zG_TE=Fp(0x1)
#
#
#

C = Bandersnatch(p, a, d, L, r, cofactor, b, c)

P = C.random_point()
Q = C.random_point()

G = C.generator()
assert G.scalar_mul(r).is_zero() and not(G.is_zero())

# psi
R = P.psi()
assert R.on_curve()
assert R.psi() == P.double().neg()

# psi is homogeneous
P = C.random_point()
lambda_P = P
lambda_P.x *= 12
lambda_P.y *= 12
lambda_P.z *= 12
R = P.psi()
lambda_R = lambda_P.psi()
assert lambda_R.normalize() == R.normalize()

# L value
P_r = C.point_of_order_r()
assert P_r.psi() == P_r.scalar_mul(C.L)

# GLV
n = random.randint(0,r)
# n = 4257185345094557079734489188109952172285839137338142340240392707284963971010
R1 = P_r.scalar_mul(n)
R2 = P_r.fast_scalar_mul(n)
assert R1 == R2

