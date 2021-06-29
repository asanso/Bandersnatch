import random
from sage.misc.misc import cputime

# Jubjub curve
from curve import Curve, Point
p_jubjub = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
a_jubjub = -1
d_jubjub = 0x2a9318e74bfa2b48f5fd9207e6bd7fd4292d7f6d37579d2601065fd6d6343eb1
r_jubjub = 0x0e7db4ea6533afa906673b0101343b00a6682093ccc81082d0970e5ed6f72cb7
J = Curve(p_jubjub, a_jubjub, d_jubjub, r_jubjub, 8)
Q = J.point_of_order_r()
m = random.randint(0,r_jubjub)
# fixed Q and m
Q = Point(8265506417937309428145637863657901096326803541555843920545271411571019173188,43078099614076354221581797221645845941231628185943136834944943003176856371463,2895621973748525387977156146651243021079099683508751012131548696447119063105,J)
m = 2363201799752215928069851892839815906013532842010879554933277259503426635572

t = cputime()
for i in range(1000):
    R1 = Q.scalar_mul(m)
t = cputime(t)
print("Jubjub:\t\t\t{:5.3f}ms".format(t))

# Bandersnatch curve
from bandersnatch import Bandersnatch, BandersnatchPoint

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

C = Bandersnatch(p, a, d,
                 L, r, cofactor,
                 a1,a2,a3, b1,b2,b3, c1,c2)
P = C.point_of_order_r()
n = random.randint(0,r)
# fixed P and n
P = BandersnatchPoint(13738737789055671334382939318077718462576533426798874551591468520593954805549,11575885077368931610486103676191793534029821920164915325066801506752632626968,14458123306641001284399433086015669988340559992755622870694102351476334505845,C)
n = 4257185345094557079734489188109952172285839137338142340240392707284963971010

t2 = cputime()
for i in range(1000):
    R2 = P.fast_scalar_mul(n)
t2 = cputime(t2)
print('Bandersnatch:\t\t{:5.3f}ms ({:2.0f}% faster)'.format(t2, 100*(t-t2)/t2))
