# to be used together with the `params-TE.py` file

    
# endomorphism
def psi(x,y,z):
    z2y = z**2*y
    fy = a1 * (y+a2*z) * (y+a3*z)
    gy = b1 * (y+b2*z) * (y+b3*z)
    hy = (y+c1*z)*(y+c2*z)
    return [x*fy*hy,gy*z2y, z2y*hy]

x = Fp.random_element()
y2 = (1-a*x**2) / (1-d*x**2)
while not(y2.is_square()):
    x = Fp.random_element()
    y2 = (1-a*x**2) / (1-d*x**2)
y = y2.sqrt()

z = Fp.random_element()
x *= z
y *= z

X,Y,Z = psi(x,y,z)
assert a*X**2*Z**2 + Y**2*Z**2 == Z**4 + d * X**2 * Y**2
