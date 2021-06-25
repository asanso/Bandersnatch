import random
from sage.misc.persist import load
from bandersnatch import Bandersnatch

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
