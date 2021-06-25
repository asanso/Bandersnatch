# to be used together with the `params-M.py` file

def psi(x,z):
    xz = x*z
    return [-(x-z)**2 - c * xz, 2*xz]

def full_psi(x,y,z):
    xz = x*z
    return [x*(-(x-z)**2 - c * xz), u*y*(x**2-z**2) ,2*x**2*z]

# test
x = Fp.random_element()
y2 = Fp((x**3+A*x**2+x)/B)
while not(y2.is_square()):
    x = Fp.random_element()
    y2 = Fp((x**3+A*x**2+x)/B)
y = y2.sqrt()

z = Fp.random_element()
x *= z
y *= z

X,Y,Z = full_psi(x,y,z)
assert B*Y**2*Z == X**3 + A*X**2*Z + X*Z**2

X2,Z2 = psi(x,z)
assert X/Z == X2/Z2

