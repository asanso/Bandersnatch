# Bandersnatch

Bandersnatch is an elliptic curve defined over the field `GF(p)` where
`p` is the BLS12-381 curve subgroup order.

The following code output is in the file `small-disc-curves.txt`:
```
sage small-disc-curves.py
```

The three files `bandersnatch-W.py`, `bandersnatch-M.py` and
`bandersnatch-TE.py` give the formulas of the endomorphism in
different models:
* in affine Weierstrass coordinates,
* in projective Montgomery coordinates,
* in Twisted Edwards projective coordinates.
