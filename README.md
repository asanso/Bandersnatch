# Bandersnatch

Bandersnatch is an elliptic curve defined over the field `GF(p)` where
`p` is the BLS12-381 curve subgroup order. It can be represented in
Twisted Edwards form `a*x**2 + y**2 = 1 + d * x**2 * y**2` where:
```python3
p = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
a = -5
d = 0x6389c12633c267cbc66e3bf86be3b6d8cb66677177e54f92b369f2f5188d58e7
```

## How we obtained Bandersnatch

We looked for small discriminant curves defined over the prime field
of characteristic `p` above.
The command 
```shell
$ make curvesearch
```
outputs the possible ordinary curves for `D < 389` (see
`small-disc-curves.txt` for the output).
Bandersnatch is the one of these curves that fits very well for
cryptographic applications: it has a subgroup of 253-bit order `r`.
```python3
# Curver order
 2^2 * 13108968793781547619861935127046491459309155893440570251786403306729687672801

# Quadratic twist order
 2^7 * 3^3 * 15172417585395309745210573063711216967055694857434315578142854216712503379
```

## Efficient Scalar multiplication

Bandersnatch has `j=8000` and `-D=-8`, meaning that it has a fast endomorphism
`psi = sqrt(-2)` leading to an efficient GLV scalar multiplication
algorithm on the subgroup of order `r`.
We provide three models in order to describe our curve, together with
the coefficients for computing this endomorphism.
Three files (`python-ref-impl/params-W.py`, `python-ref-impl/params-M.py` and
`python-ref-impl/params-TE.py`) can be generated using
```shell
make getparams
```
* In affine Weierstrass coordinates.<br>
The file `python-ref-impl/params-W.py` contains the curve parameters `p, a, b` such
that the Weierstrass equation of the curve over `GF(p)` is `y**2 =
x**3 + a*x + b`. It also includes the coefficients
`r0,r1,s0,t0,t1,u2,u3` such that
```python3
psi_W(x,y) = ( u2*((x+r1)*x+r0)/(x+s0) , u3*y*((x+t1)*x+t0)/(x+s0)**2 )
```
* In projective `x-z` Montgomery coordinates.<br>
The file `python-ref-impl/params-M.py` includes the parameters of the Montgomery curve
`B*y**2 = x**3 + A*x**2 + x` over `GF(p)`, together with the
coefficient `c` for the endomorphism in projective `x-z` coordinates:
```python3
psi_M(x,z) = ( -(x-z)**2 - c * x * z , 2 * x * z )
```
* In projective Twisted Edwards coordinates.<br>
The file `python-ref-impl/params-TE.py` includes the parameters of the curve in the
form `a*x**2 + y**2 = 1 + d * x**2 * y**2`, and the endomorphism
coefficients `a1,a2,a3,b1,b2,b3,c1,c2` such that
```python3
psi_TE(x,y,z) = ( x * a1 * (y+a2*z) * (y+a3*z) * (y+c1*z) * (y+c2*z) ,
b1 * (y+b2*z) * (y+b3*z) * z**2 * y , (y+c1*z) * (y+c2*z) * z**2 * y )
```

## Comparison with Jubjub

We implemented a non-optimized elliptic curve group law arithmetic in
`python-ref-impl/curve.py` and `python-ref-impl/bandersnatch.py`.
Our benchmarks lead to a Bandersnatch scalar multiplication ~35%
faster than on Jubjub.
This estimation can be reproducible using
```shell
$ make bench
```
and our code can be tested using
```shell
$ make test
```
We also provide Rust benchmarks, where we obtained a 42%
improvement. TODO coming soon.

## Technical details
We are writing a paper for more technical details. It is still a work
in progress, available using
```shell
make pdf
evince paper/bandersnatch.pdf
```

## Reference implementation
We provide two reference implementations
- [python](https://github.com/asanso/Bandersnatch/tree/main/python-ref-impl)
- [rust](https://github.com/zhenfeizhang/bandersnatch)
