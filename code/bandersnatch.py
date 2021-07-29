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
                 a1,a2,a3, b1,b2,b3, c1,c2):
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
        self.a3 = a3
        self.b1 = b1
        self.b2 = b2
        self.b3 = b3
        self.c1 = c1
        self.c2 = c2

    def random_point(self):
        P = super().random_point()
        return BandersnatchPoint(P.x, P.y, P.z, self)

    def point_of_order_r(self):
        P = super().point_of_order_r()
        return BandersnatchPoint(P.x, P.y, P.z, self)

    def generator(self):
        p=52435875175126190479447740508185965837690552500527637822603658699938581184513
        Fp=FiniteField(p)
        x=Fp(0x14acd3bb6ae0c7074ac851ebf149c7b735178562dd178209297d705d8c618427)
        y=Fp(0x3350e8e5c888b3387da2c9332cb6b3bdee7a8872268184c136945494d48e6888)
        z=Fp(0x1)
        return BandersnatchPoint(x, y, z, self)

    
class BandersnatchPoint(Point):
    
    def __init__(self, x, y, z, curve):
        super().__init__(x, y, z, curve)
        
    def psi(self):
        x,y,z = self.x, self.y, self.z
        a1 = self.curve.a1
        a2 = self.curve.a2
        a3 = self.curve.a3
        b1 = self.curve.b1
        b2 = self.curve.b2
        b3 = self.curve.b3
        c1 = self.curve.c1
        c2 = self.curve.c2
        z2y = z**2*y
        fy = a1 * (y+a2*z) * (y+a3*z)
        gy = b1 * (y+b2*z) * (y+b3*z)
        hy = (y+c1*z)*(y+c2*z)
        return BandersnatchPoint(x*fy*hy, gy*z2y, z2y*hy, self.curve)

    def fast_scalar_mul(self, n):
        psiP = self.psi()
        beta = vector([n,0]) * self.curve.N_inv
        b = vector([int(beta[0]), int(beta[1])]) * self.curve.N
        k1 = n-b[0]
        k2 = -b[1]
        return self.multi_scalar_mul(k1, psiP, k2)
