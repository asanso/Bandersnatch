import random
from sage.misc.persist import load
from curve import Curve
from bandersnatch import Bandersnatch


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

# clear cofactor
assert P.clear_cofactor().scalar_mul(E.r).is_zero()


###
# Test for Bandersnatch
###

load('params-TE.py')

C = Bandersnatch(p, a, d, L, r, cofactor, a1,a2,a3, b1,b2,b3, c1,c2)

P = C.random_point()
Q = C.random_point()

# psi
R = P.psi()
assert R.on_curve()
assert R.psi() == P.double().neg()

# L value
P_r = C.point_of_order_r()
assert P_r.psi() == P_r.scalar_mul(C.L)

# GLV
n = random.randint(0,r)
R1 = P_r.scalar_mul(n)
R2 = P_r.fast_scalar_mul(n)
assert R1 == R2
