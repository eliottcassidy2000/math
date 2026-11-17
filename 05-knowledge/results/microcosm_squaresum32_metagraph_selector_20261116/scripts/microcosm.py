#!/usr/bin/env python3
"""microcosm.py (M-1): the first 32 odd numbers, their least orders n(h) from the census to order 9 (counts.txt), the split
3 | 1 | 6 | 1 | 21 around the holes, triangular block sizes, block sums, the tripling chains 7·3^k and 3^k, and the
sums of the first 32 odd numbers."""
from collections import defaultdict
first = {}
for line in open('counts.txt'):
    n, h = map(int, line.split())
    if h not in first or n < first[h]: first[h] = n
odd = list(range(1, 64, 2)); print("the 32 odd numbers below 64 with least order n(h) (- = hole):", [(h, first.get(h, '-')) for h in odd])
blocks = [[1, 3, 5], [7], list(range(9, 20, 2)), [21], list(range(23, 64, 2))]
tri = lambda m: int((8 * m + 1) ** 0.5) ** 2 == 8 * m + 1
print("block sizes:", [len(b) for b in blocks], "triangular:", [tri(len(b)) for b in blocks], "; block sums:", [sum(b) for b in blocks], "triangular:", [tri(sum(b)) for b in blocks])
print("least orders by block:", [[first.get(h, '-') for h in b] for b in blocks])
print("sum of the first 32 odd numbers:", sum(odd), "= 32^2 =", 32 ** 2)
S = set(first)
print("tripling chain of 7:", [(7 * 3 ** k, 7 * 3 ** k in S, first.get(7 * 3 ** k, '-')) for k in range(5)], "; of 1:", [(3 ** k, 3 ** k in S, first.get(3 ** k, '-')) for k in range(7)])
print("multiples of 7 below 200 that are counts (orders <= 9):", [m for m in range(7, 200, 14) if m in S], "; not counts:", [m for m in range(7, 200, 14) if m not in S])
print("odd numbers in [65, 2881] not realised at orders <= 9 (first 30):", [h for h in range(65, 2882, 2) if h not in S][:30], "; count:", sum(1 for h in range(65, 2882, 2) if h not in S))
