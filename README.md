# Bandersnatch

Bandersnatch is an elliptic curve defined over the field `GF(p)` where
`p` is the BLS12-381 curve subgroup order.
```python3
p = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
a = -5
d = 0x6389c12633c267cbc66e3bf86be3b6d8cb66677177e54f92b369f2f5188d58e7
```

## Finding Bandersnatch
Bandersnatch is a curve of small discriminant `-D = -4*2` and `j=8000`.
The command outputs the file `small-disc-curves.txt`:
```shell
$ sage small-disc-curves.py
```
It computes all the possible ordinary curves with `D < 389`. Only
Bandersnatch has a cryptographically secure structure:
```python3
Curver order:
 2^2 * 13108968793781547619861935127046491459309155893440570251786403306729687672801

Quadratic twist order:
 2^7 * 3^3 * 15172417585395309745210573063711216967055694857434315578142854216712503379
```

## Group law in different models

The three files `bandersnatch-W.py`, `bandersnatch-M.py` and
`bandersnatch-TE.py` give the formulas of the endomorphism in
different models:
* in affine Weierstrass coordinates,
* in projective Montgomery coordinates,
* in Twisted Edwards projective coordinates.


## Accelerated scalar multiplication

The file `bandersnatch-TE.py` estimates the cost of a scalar
multiplication with and without the GLV method. *We obtain an algorithm
30% faster*:
```shell
$ sage bandersnatch-TE.py
Scalar multiplication: 1.831ms
GLV scalar multiplication: 1.253ms

```
