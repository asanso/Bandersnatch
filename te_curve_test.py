from te_curve import TE_Curve

p = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
a = -1
d = 0x2a9318e74bfa2b48f5fd9207e6bd7fd4292d7f6d37579d2601065fd6d6343eb1
r = 0x0e7db4ea6533afa906673b0101343b00a6682093ccc81082d0970e5ed6f72cb7
E = TE_Curve(p, a, d, r, 8)

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
