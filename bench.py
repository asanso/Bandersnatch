import random
from sage.misc.misc import cputime
from sage.misc.persist import load

# Jubjub curve
from te_curve import TE_Curve
p_jubjub = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
a_jubjub = -1
d_jubjub = 0x2a9318e74bfa2b48f5fd9207e6bd7fd4292d7f6d37579d2601065fd6d6343eb1
r_jubjub = 0x0e7db4ea6533afa906673b0101343b00a6682093ccc81082d0970e5ed6f72cb7
J = TE_Curve(p_jubjub, a_jubjub, d_jubjub, r_jubjub, 8)
Q = J.point_of_order_r()
m = random.randint(0,r_jubjub)

t = cputime()
for i in range(1000):
    R1 = Q.scalar_mul(m)
t = cputime(t)
print("Jubjub:\t\t\t{:5.3f}ms".format(t))

# Bandersnatch curve
from bandersnatch import Bandersnatch
load('params-TE.py')
C = Bandersnatch(p, a, d,
                 L, r, cofactor,
                 a1,a2,a3, b1,b2,b3, c1,c2)
P = C.point_of_order_r()
n = random.randint(0,r)

t2 = cputime()
for i in range(1000):
    R2 = P.fast_scalar_mul(n)
t2 = cputime(t2)
print('Bandersnatch:\t\t{:5.3f}ms ({:2.0f}% faster)'.format(t2, 100*(t-t2)/t2))
