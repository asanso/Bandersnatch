import sage.all
from sage.rings.integer_ring import Z as ZZ
from sage.matrix.constructor import Matrix
from sage.rings.finite_rings.finite_field_constructor import FiniteField
from sage.modules.free_module_element import free_module_element as vector


def base64FieldInt(a, hex_flag=False):
    a = ZZ(a)
    a1 = a % (2**64)
    a = a >> 64
    a2 = a % (2**64)
    a = a >> 64
    a3 = a % (2**64)
    a = a >> 64
    a4 = a % (2**64)
    if hex:
        return '[' + str(hex(a1)) + ',' + str(hex(a2)) + ',' + str(hex(a3)) \
            + ',' + str(hex(a4)) + ']'
    else:
        return '[' + str(a1) + ',' + str(a2) + ',' + str(a3) + ',' + str(a4) + ']'
    
class Curve():
    
    def __init__(self, p, a, d, r, cofactor):
        self.p = p
        self.Fp = FiniteField(p)
        self.a = self.Fp(a)
        self.d = self.Fp(d)
        self.cofactor = cofactor
        self.r = r

    def __str__(self):
        a = ZZ(self.a)
        d = ZZ(self.d)
        p = ZZ(self.p)
        if len(str(a-p)) < 5 :
            a = a-p
        if len(str(d-p)) < 5 :
            d = d-p
        return "Elliptic curve with equation {} * x² + y² = 1 + {} * x² * y² defined over GF({})".format(a,d,p)

    def random_point(self):
        x = self.Fp.random_element()
        y2 = (1-self.a*x**2) / (1-self.d*x**2)
        while not(y2.is_square()):
            x = self.Fp.random_element()
            y2 = (1-self.a*x**2) / (1-self.d*x**2)
        y = y2.sqrt()
        z = self.Fp.random_element()
        x *= z
        y *= z
        return Point(x, y, z, self)
    
    def point_of_order_r(self):
        P = self.random_point().clear_cofactor()
        while P.is_zero():
            P = self.random_point().clear_cofactor()
        assert P.scalar_mul(self.r).is_zero()
        return P

class Point():
    
    def __init__(self, x, y, z, curve):
        self.x = x
        self.y = y
        self.z = z
        self.curve = curve

    def __str__(self):
        return 'x:' + base64FieldInt(self.x, True) + '\n' + \
            'y:' + base64FieldInt(self.y, True) + '\n' + \
            'z:' + base64FieldInt(self.z, True)
        
    def __eq__(self, other):
        return self.x * other.z == \
            other.x * self.z and self.y * other.z == other.y*self.z

    def on_curve(self):
        X,Y,Z = self.x, self.y, self.z
        a = self.curve.a
        d = self.curve.d
        return a*X**2*Z**2 + Y**2*Z**2 == Z**4 + d * X**2 * Y**2

    def neg(self):
        return Point(-self.x, self.y, self.z, self.curve)
            
    def is_zero(self):
        return self.x.is_zero() and self.y == self.z
    
    def normalize(self):
        return Point(self.x/self.z, self.y/self.z, 1, self.curve)

    def add(self, other):
        x,y,z = self.x, self.y, self.z
        xx,yy,zz = other.x, other.y, other.z
        A = z*zz
        B = A**2
        C = x*xx
        D = y*yy
        E = self.curve.d*C*D
        F = B-E
        G = B+E
        X = A*F*((x+y) * (xx+yy) - C - D)
        Y = A*G*(D-self.curve.a*C)
        Z = F*G
        return Point(X, Y, Z, self.curve)

    def double(self):
        x,y,z = self.x, self.y, self.z
        B = (x+y)**2
        C = x**2
        D = y**2
        E = self.curve.a*C
        F = E+D
        H = z**2
        J = F-2*H
        X = (B-C-D)*J
        Y = F*(E-D)
        Z = F*J
        return Point(X, Y, Z, self.curve)

    def scalar_mul(self, n):
        if n<0:
            n = -n
            P = self.neg()
        else:
            P = self
        R = P
        for b in ZZ(n).bits()[-2::-1]:
            R = R.double()
            if b == 1:
                R = R.add(P)
        return R

    def multi_scalar_mul(self, k1, other, k2):
        P = self
        if k1<0:
            k1=-k1
            P = P.neg()
        if k2<0:
            k2=-k2
            other = other.neg()
        PplusOther = P.add(other)
        bits_k1 = ZZ(k1).bits()
        bits_k2 = ZZ(k2).bits()
        while len(bits_k1) < len(bits_k2):
            bits_k1.append(0)
        while len(bits_k2) < len(bits_k1):
            bits_k2.append(0)
        R = Point(0, 1, 1, self.curve)
        for i in range(len(bits_k1)-1,-1,-1):
            R = R.double()
            if bits_k1[i] == 1 and bits_k2[i] == 0:
                R = R.add(self)
            if bits_k1[i] == 0 and bits_k2[i] == 1:
                R = R.add(other)
            if bits_k1[i] == 1 and bits_k2[i] == 1:
                R = R.add(PplusOther)
        return R

    def clear_cofactor(self):
        return self.scalar_mul(self.curve.cofactor)
