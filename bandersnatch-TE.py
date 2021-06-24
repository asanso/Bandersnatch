import random
import sage.all
from sage.misc.misc import cputime
from sage.rings.integer_ring import Z as ZZ
from sage.matrix.constructor import Matrix
from sage.rings.finite_rings.finite_field_constructor import FiniteField
from sage.modules.free_module_element import free_module_element as vector

# Twisted Edwards Bandersnatch curve

p=52435875175126190479447740508185965837690552500527637822603658699938581184513
Fp = FiniteField(p)

a=Fp(-5)
d=Fp(138827208126141220649022263972958607803)/Fp(171449701953573178309673572579671231137)

def random_point(a,d):
    x = Fp.random_element()
    y2 = (1-a*x**2) / (1-d*x**2)
    while not(y2.is_square()):
        x = Fp.random_element()
        y2 = (1-a*x**2) / (1-d*x**2)
    y = y2.sqrt()
    
    z = Fp.random_element()
    x *= z
    y *= z
    return [x,y,z]

def on_curve(P):
    X,Y,Z = P
    return a*X**2*Z**2 + Y**2*Z**2 == Z**4 + d * X**2 * Y**2

def neg(P):
    return [-P[0], P[1], P[2]]

def eq(P,Q):
    return P[0] * Q[2] == Q[0]*P[2] and P[1] * Q[2] == Q[1]*P[2]

def is_zero(P):
    return P[0] == 0 and P[1] == P[2]

def normalize(P):
    return [P[0]/P[2], P[1]/P[2], 1]

# endomorphism

a1=Fp(16179988757916560824577558193084210236647645729299773892093730683504906651604)
a2=Fp(37446463827641770816307242315180085052603635617490163568005256780843403514036)
a3=Fp(14989411347484419663140498193005880785086916883037474254598401919095177670477)
b1=Fp(37446463827641770816307242315180085052603635617490163568005256780843403514036)
b2=Fp(36553259151239542273674161596529768046449890757310263666255995151154432137034)
b3=Fp(15882616023886648205773578911656197791240661743217374156347663548784149047479)
c1=Fp(42910309089382041158038545419309140955400939872179826051492616687477682993077)
c2=Fp(9525566085744149321409195088876824882289612628347811771111042012460898191436)

def psi(P):
    x,y,z = P
    z2y = z**2*y
    fy = a1 * (y+a2*z) * (y+a3*z)
    gy = b1 * (y+b2*z) * (y+b3*z)
    hy = (y+c1*z)*(y+c2*z)
    return [x*fy*hy,gy*z2y, z2y*hy]

P = random_point(a,d)
R = psi(P)
assert on_curve(R)


# addition

def add(P, Q):
    x,y,z = P
    xx,yy,zz = Q
    A = z*zz
    B = A**2
    C = x*xx
    D = y*yy
    E = d*C*D
    F = B-E
    G = B+E
    X = A*F*((x+y) * (xx+yy) - C - D)
    Y = A*G*(D-a*C)
    Z = F*G
    return [X,Y,Z]

P = random_point(a,d)
Q = random_point(a,d)
R = add(P,Q)
assert on_curve(R)


# doubling

def double(P):
    x,y,z = P
    B = (x+y)**2
    C = x**2
    D = y**2
    E = a*C
    F = E+D
    H = z**2
    J = F-2*H
    X = (B-C-D)*J
    Y = F*(E-D)
    Z = F*J
    return [X,Y,Z]

P = random_point(a,d)
R = double(P)
assert on_curve(R)

# checking that the endomorphism is sqrt(-2)
P = random_point(a,d)
R = psi(psi(P))
assert eq(R, neg(double(P)))


# scalar multiplication
def scalar_mul(n, P):
    if n<0:
        n = -n
        P = neg(P)
    R = P
    for b in ZZ(n).bits()[-2::-1]:
        R = double(R)
        if b == 1:
            R = add(R,P)
    return R

# multi scalar multiplication
def multi_scalar_mul(k1,P,k2,Q):
    if k1<0:
        k1=-k1
        P = neg(P)
    if k2<0:
        k2=-k2
        Q = neg(Q)
    PplusQ = add(P,Q)
    bits_k1 = ZZ(k1).bits()
    bits_k2 = ZZ(k2).bits()
    while len(bits_k1) < len(bits_k2):
        bits_k1.append(0)
    while len(bits_k2) < len(bits_k1):
        bits_k2.append(0)
    R = [0,1,1]
    for i in range(len(bits_k1)-1,-1,-1):
        R = double(R)
        if bits_k1[i] == 1 and bits_k2[i] == 0:
            R = add(R, P)
        if bits_k1[i] == 0 and bits_k2[i] == 1:
            R = add(R, Q)
        if bits_k1[i] == 1 and bits_k2[i] == 1:
            R = add(R, PplusQ)
    return R

P = random_point(a,d)
n = 5
R = scalar_mul(n,P)
assert on_curve(R)

# GLV scalar multiplication
r = 13108968793781547619861935127046491459309155893440570251786403306729687672801
cofactor = 4
curve_order = r * cofactor

P = random_point(a,d)
assert is_zero(scalar_mul(curve_order, P))

def clear_cofactor(P):
    R = scalar_mul(cofactor,P)
    while is_zero(R):
        R = scalar_mul(cofactor,P)
    return R

# L is the eigen value of psi on the r-torsion subgroup of the curve
L=8913659658109529928382530854484400854125314752504019737736543920008458395397
P = clear_cofactor(random_point(a,d))
assert eq(scalar_mul(L,P), psi(P))
M = Matrix([[-L,1], [r,0]])
N = M.LLL()
N_inv = N**-1

def fast_scalar_mul(n,P):
    '''
    Compute [n] P = [k1]P + [k2] psi(P) for a point P in the r-torsion subgroup
    '''
    # endomorphism precomputation
    psiP = psi(P)
    # scalar decomposition
    alpha = vector([n,0]) * N_inv
    v = vector([int(alpha[0]), int(alpha[1])]) * N
    k1 = n-v[0]
    k2 = -v[1]
    return multi_scalar_mul(k1,P,k2,psiP)

# random scalar
n = random.randint(0,r)
# random point of order r
P = clear_cofactor(random_point(a,d))
assert is_zero(scalar_mul(r,P))

t = cputime()
for i in range(1000):
    R1 = scalar_mul(n,P)
print('Scalar multiplication: {}ms'.format(round(cputime(t),3)))

t = cputime()
for i in range(1000):
    R2 = fast_scalar_mul(n,P)
print('GLV scalar multiplication: {}ms'.format(round(cputime(t),3)))

assert eq(R1,R2)
