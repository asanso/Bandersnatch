# to be used together with the `params-W.py` file
def psi(x,y):
    rx = (x+r1)*x+r0
    sx = x+s0
    sx2 = sx**2
    tx = (x+t1)*x+t0
    return [u2*rx/sx, u3*y*tx/sx2]

# test
x = Fp.random_element()
y2 = x**3+a*x+b
while not(y2.is_square()):
    x = Fp.random_element()
    y2 = x**3+a*x+b
y = y2.sqrt()

X,Y = psi(x,y)
assert Y**2 == X**3+a*X+b
