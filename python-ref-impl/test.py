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

# clear cofactor
assert P.clear_cofactor().scalar_mul(E.r).is_zero()


###
# Test for Bandersnatch
###

# this code can be recomputed using `sage get_params.py`
import sage.all
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
a1=Fp(16179988757916560824577558193084210236647645729299773892093730683504906651604)
a2=Fp(37446463827641770816307242315180085052603635617490163568005256780843403514036)
a3=Fp(14989411347484419663140498193005880785086916883037474254598401919095177670477)
b1=Fp(37446463827641770816307242315180085052603635617490163568005256780843403514036)
b2=Fp(36553259151239542273674161596529768046449890757310263666255995151154432137034)
b3=Fp(15882616023886648205773578911656197791240661743217374156347663548784149047479)
c1=Fp(42910309089382041158038545419309140955400939872179826051492616687477682993077)
c2=Fp(9525566085744149321409195088876824882289612628347811771111042012460898191436)
# endomorphism eigenvalue L on the subgroup of size r
r=13108968793781547619861935127046491459309155893440570251786403306729687672801
cofactor = 4
curve_order=r*cofactor
L=8913659658109529928382530854484400854125314752504019737736543920008458395397

C = Bandersnatch(p, a, d, L, r, cofactor, a1,a2,a3, b1,b2,b3, c1,c2)

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
