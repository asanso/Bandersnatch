from sage.rings.finite_rings.finite_field_constructor import FiniteField

# Weierstrass Bandersnatch curve in affine coordinates

p=52435875175126190479447740508185965837690552500527637822603658699938581184513
Fp = FiniteField(p)

a=Fp(20856223409030707214153620029242927260100058688354534505372746703669794261922)
b=Fp(4118395616383396936705917269457205707251888519744679123386326773690169742488)

r0=Fp(8460616024632051863286924185728629979016185787198334425817815457773555916652)
r1=Fp(14417655457111367516429032275726150553899552124793034239520187109589559348972)
s0=Fp(14417655457111367516429032275726150553899552124793034239520187109589559348972)
t0=Fp(4886139905610459846240090702495259539852707412592144821865093693573891696429)
t1=Fp(28835310914222735032858064551452301107799104249586068479040374219179118697944)
u2=Fp(26217937587563095239723870254092982918845276250263818911301829349969290592256)
u3=Fp(44345880796167910065426388998607034559978594870533879487421667123803446877207)

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
