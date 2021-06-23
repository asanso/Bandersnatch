from sage.all import *
from sage.rings.finite_rings.finite_field_constructor import FiniteField

# Montgomery Bandersnatch curve in xz-only coordinates

p=52435875175126190479447740508185965837690552500527637822603658699938581184513
Fp = FiniteField(p)

B=Fp(33322015990175592254575859866067371076354504870126610304022993361853120164277)
A=Fp(22457052480157351153166744122174204267516718734452689313406854861748225843561)

c=Fp(22457052480157351153166744122174204267516718734452689313406854861748225843563)

def psi(x,z):
    xz = x*z
    return [-(x-z)**2 - c * xz, 2*xz]

# test
x = Fp.random_element()
y2 = Fp((x**3+A*x**2+x)/B)
while not(y2.is_square()):
    x = Fp.random_element()
    y2 = Fp((x**3+A*x**2+x)/B)
y = y2.sqrt()

z = Fp.random_element()
x *= z

X,Z = psi(x,z)
assert (((X/Z)**3+A*(X/Z)**2+X/Z)/B).is_square()
