"""
The even-run condition <=> the consecutive-value blocks all have ODD size <=> compositions of n into ODD parts.
CLOSED FORM:  G(n) = Sum_{compositions of n into odd parts} k! * 2^{#parts>=3},  k=#parts.
Then Var(H) = (n!/4^{n-1})(G(n)-n!).
"""
from functools import lru_cache
import math
def comps_odd(n):
    # yield compositions of n into odd parts
    if n==0: yield (); return
    for first in range(1,n+1,2):
        for rest in comps_odd(n-first):
            yield (first,)+rest
def G_comp(n):
    tot=0
    for c in comps_odd(n):
        k=len(c); big=sum(1 for p in c if p>=3)
        tot+=math.factorial(k)*2**big
    return tot
print("verify composition-into-odd-parts closed form:")
for n in range(3,9):
    print(f"  n={n}: G(n)={G_comp(n)}")
print(f"\nG(n) sequence n=1..12:", [G_comp(n) for n in range(1,13)])
from fractions import Fraction as F
print("\nVar(H) = (n!/4^(n-1))(G(n)-n!):")
for n in range(3,9):
    v=F(math.factorial(n),4**(n-1))*(G_comp(n)-math.factorial(n))
    print(f"  n={n}: Var(H) = {v} = {float(v):.4f}")
