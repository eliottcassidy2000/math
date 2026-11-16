# A smaller finite frontier for the prime odd wheel layers

For an odd prime q, if period 2^a q covers with a≤5, then at least q
distinct primes of odd order part q divide 2^(32q)−1. Each is of
the form p_j=1+2q k_j; order the k_j increasingly. For q≥52000,
none of these primes is 3,5,7,11. Thus the allowable k_j occupy
at most φ(1155)=480 classes modulo 1155=3·5·7·11, because 2q
is invertible modulo 1155. A block bound gives

    k_j ≥ 1155 floor((j−1)/480)+1 ≥ (77/32)j−1157.

Also k_j≥j. Therefore k_j/j has the following certified lower
bounds for j starting at the listed thresholds:

    j    1500  3000  6000  12000  24000
    bound 8/5     2  11/5  23/10  47/20.

At q=52000 the product of the stepwise ratios exceeds
(11/5)^52000, by the integer check below; all further factors
contribute at least 47/20>11/5. The product of the required primes
would therefore be greater than

    [(2q²/e)(11/5)]^q > [(8/5)q²]^q > 2^(32q),

using q!≥(q/e)^q and e<11/4. This is impossible.
Thus a=6 is the exact least dyadic exponent for every prime
q≥52000, with the ambient-period convention in the
[wheel proof](wheel/WHEEL_FOLLOWUP.md). This is a cutoff for the
remaining factorization work, not a classification of the primes
below 52000.

Exact arithmetic for the two numerical inequalities:

```python
from math import prod
groups=[(8,5,1500),(2,1,3000),(11,5,6000),(23,10,12000),(47,20,28000)]
assert prod(n**k for n,d,k in groups)*5**52000 > prod(d**k for n,d,k in groups)*11**52000
assert 8*52000**2 > 5*2**32
```

A counterexample to the unqualified nesting claim even with a nontrivial
odd part: q=8,589,934,609 is prime, as checked by trial division to its
integer square root. It exceeds 2^33. At ambient period 64q, W(64)
already covers. But an odd-only cover would require q distinct labels
1 mod 2q, whose product is at least (2q²/e)^q>2^(64q).
This proves the dyadic-only alternative is not just a convention for
the degenerate m=1 case.
