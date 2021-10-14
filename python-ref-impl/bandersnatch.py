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
                 b,c):
        super().__init__(p, a, d, r, cofactor)
        self.D = -8
        self.L = L
        self.cofactor = cofactor
        self.r = r
        M = Matrix([[-L,1], [r,0]])
        self.N = M.LLL()
        self.N_inv = self.N**-1
        self.b = b
        self.c = c
        
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
        # computed in 9 multiplications and 3 additions.
        x,y,z = self.x, self.y, self.z
        b = self.curve.b
        c = self.curve.c
        xy = x*y
        y2 = y**2
        z2 = z**2
        bz2 = b*z2
        fy = c * (z2-y2)
        gy = b * (y2 + bz2)
        hy = y2 - bz2
        return BandersnatchPoint(fy*hy, gy*xy, hy * xy, self.curve)
    
    def fast_scalar_mul(self, n):
        psiP = self.psi()
        beta = vector([n,0]) * self.curve.N_inv
        b = vector([int(beta[0]), int(beta[1])]) * self.curve.N
        k1 = n-b[0]
        k2 = -b[1]
        return self.multi_scalar_mul(k1, psiP, k2)
