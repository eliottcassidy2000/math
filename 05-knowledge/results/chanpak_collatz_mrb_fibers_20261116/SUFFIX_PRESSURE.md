# Sharper MRB pressure by keeping the live suffix

The [grouped-tail pressure](MRB_PRESSURE.md) weights *every* isolated block or block group. The same idea is substantially stronger if we retain the state known to be present at the true seam. For the established 168-block bank this gives **o(X^(23/30))** for the full language, and **o(X^(25/36))** for the mixed binary/decimal controller with 111 distinct actions. Ordinal-sum closure remains power-sparse at these same exponents.

## Conditional growth without an LP

In the nonnegative literal frontier transfer let f read H, and let P be the transfer of the common suffix S=(0,3,5,2,4). The physical seed ends in S, and so does every completed return block. Thus every boundary raw state is Pz for a nonnegative integer vector z. This is stronger than knowing the rank or the abstract cone image.

For each physical block w with raw transfer B_w define
```
λ_w = min_{j:(fP)_j>0} (f B_w P)_j/(fP)_j .
```
Numerator coordinates with denominator zero are nonnegative and need no comparison. Therefore H after appending w is at least λ_w times its old value, on **every** physical boundary state. This is a genuine absolute growth comparison, not the count of w on an isolated prefix. The exact rational minimum and maximum across all 168 constants are 2439/5 and 39436/25. The marked suffix is doing real work: the worst isolated block had count only 55.

For the invertible-at-five bank retain the last listed physical representative of each exact integer action on the rank-four return lattice. Replacing a complete word by representatives does not change the semantic state at any block boundary; a partial stop is evaluated from that same state. Use each representative's own λ_w in the lower growth bound; no equality of raw transfer matrices is assumed.

## Directed integer test for the whole tail

For a positive rational s, if Σ_w λ_w^(-s)<1, then after m blocks
```
Σ_{words of m blocks} H(word)^(-s) ≤ H(seed)^(-s) [Σ_w λ_w^(-s)]^m.
```
Summing in m (also allowing a finite menu of partial stops) gives finite Dirichlet mass. For a rational λ_w=a/b and s=p/q the replay uses the least integer n with n^q a^p≥10^(12q) b^p as an upward numerator. Summation gives:

| Bank | s | Exact upper numerator over 10^12 |
|---|---:|---:|
| All 168 | 23/30 | 991486913547 |
| Mixed 111 | 25/36 | 996591617658 |

Both are strictly less than 1. Finite mass implies o(X^s), not just O(X^s): the contribution from counts above an arbitrary fixed cutoff bounds the normalized counting limsup by a tail of the convergent series. The multiset Euler product of [the grouped proof](MRB_PRESSURE.md) gives the same o-result for ordinal products. For the mixed bank this leaves runs of **ω(X^(11/36)) omitted consecutive odds** in dyadic bands. If all odd classes modulo an even M are reached, the worst least witness is **ω(M^(36/25))**. Neither claim extends to other return languages.

Run `python check_suffix_pressure.py` with the earlier sibling lattice packets and generated automaton. Every λ is evaluated by exact Fraction comparisons, all raw integer matrix products are checked for overflow, and every fractional power bound used in the sums is replayed with integer exponent comparisons. A constructive next test is to compare pressure after retaining the *actual* replacement seam for each changing core, rather than giving all proposed moves the generous cost of isolation.
