#!/usr/bin/env python3
"""Rank distribution of a uniformly random alternating m x m matrix over F_2 (exact), for comparison with the census
of beta_R = (mu_R / 2 mod 2).  N(m, 2r) = prod_{i=1..r} 2^(2i-2)/(2^(2i)-1) * prod_{i=0..2r-1} (2^(m-i) - 1)."""
from fractions import Fraction as Fr
def dist(m):
    out = {}
    for r in range(m//2 + 1):
        num = Fr(1)
        for i in range(1, r + 1): num *= Fr(2**(2*i - 2), 2**(2*i) - 1)
        for i in range(2*r): num *= (2**(m - i) - 1)
        out[2*r] = num / 2**(m*(m - 1)//2)
    assert sum(out.values()) == 1
    return out
for m in range(4, 13):
    d = dist(m)
    print(m, {k: round(float(v), 5) for k, v in d.items() if v > 1e-6}, " E[2^(-rank/2)] =", round(float(sum(v*Fr(1, 2**(k//2)) for k, v in d.items())), 6))
