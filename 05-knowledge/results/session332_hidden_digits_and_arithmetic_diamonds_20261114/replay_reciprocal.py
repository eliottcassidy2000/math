"""Finite exact frontier for the generalized reciprocal cycle equation."""
from math import isqrt
from sympy import divisors

hits=[]
for j in range(1,62):
    q=1<<j
    for d in divisors(q-1):
        p=isqrt(q+d)
        H=(q-1)//d
        if p>=3 and p%2 and H%2 and p*p==q+d:
            hits.append((p,j,H))
if hits!=[(3,3,7)]:
    raise RuntimeError(("unexpected reciprocal loop",hits))
print("exact reciprocal loop catalogue through exponent 61:",hits)
