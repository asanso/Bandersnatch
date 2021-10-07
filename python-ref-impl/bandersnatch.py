import sage.all
from sage.rings.integer_ring import Z as ZZ
from sage.matrix.constructor import Matrix
from sage.rings.finite_rings.finite_field_constructor import FiniteField
from sage.modules.free_module_element import free_module_element as vector
from curve import Curve
from curve import Point

class Bandersnatch(Curve):
    
    def __init__(self,
                 p, a, d,
                 L, r, cofactor,
                 a1,a2, b1,b2, c1):
        super().__init__(p, a, d, r, cofactor)
        self.D = -8
        self.L = L
        self.cofactor = cofactor
        self.r = r
        M = Matrix([[-L,1], [r,0]])
        self.N = M.LLL()
        self.N_inv = self.N**-1
        self.a1 = a1
        self.a2 = a2
        self.b1 = b1
        self.b2 = b2
        self.c1 = c1

    def random_point(self):
        P = super().random_point()
        return BandersnatchPoint(P.x, P.y, P.z, self)

    def point_of_order_r(self):
        P = super().point_of_order_r()
        return BandersnatchPoint(P.x, P.y, P.z, self)

    def generator(self):
        p=52435875175126190479447740508185965837690552500527637822603658699938581184513
        Fp=FiniteField(p)
        x=Fp(0x29c132cc2c0b34c5743711777bbe42f32b79c022ad998465e1e71866a252ae18)
        y=Fp(0x2a6c669eda123e0f157d8b50badcd586358cad81eee464605e3167b6cc974166)
        z=Fp(0x1)
        return BandersnatchPoint(x, y, z, self)

    
class BandersnatchPoint(Point):
    
    def __init__(self, x, y, z, curve):
        super().__init__(x, y, z, curve)
        
    def psi(self):
        x,y,z = self.x, self.y, self.z
        a1 = self.curve.a1
        a2 = self.curve.a2
        b1 = self.curve.b1
        b2 = self.curve.b2
        c1 = self.curve.c1
        z2 = z**2
        y2 = y**2
        fy = a1*y2+a2*z2
        gy = b1*y2+b2*z2
        hy = y2+c1*z2
        return BandersnatchPoint(x*fy*hy, gy*z2y, hy*z2*y, self.curve)

    def fast_scalar_mul(self, n):
        psiP = self.psi()
        beta = vector([n,0]) * self.curve.N_inv
        b = vector([int(beta[0]), int(beta[1])]) * self.curve.N
        k1 = n-b[0]
        k2 = -b[1]
        return self.multi_scalar_mul(k1, psiP, k2)
