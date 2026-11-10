#!/usr/bin/env python3
"""Section 4.2: the one-apex words 1^(d-1) 0 b 1^d 0 (read from the source of the chain to its sink; the apex beats
the vertices marked 1).  With v the integer whose binary digits are the word (first letter most significant) and
g(0) = 1, g(2v+1) = g(v), g(2v) = g(v) + v + (v mod 2) (session25; THM-4885 section 0), the counts are
2^(2d+1) - 3*2^(d-1) - 1 (b = 0) and 2^(2d+1) - 3*2^(d-1) + 1 (b = 1)."""
from functools import lru_cache
@lru_cache(None)
def g(v):
    if v == 0: return 1
    if v & 1: return g(v >> 1)
    w = v >> 1
    return g(w) + w + (w & 1)
for d in range(2, 12):
    lo = int("1"*(d - 1) + "00" + "1"*d + "0", 2); hi = int("1"*(d - 1) + "01" + "1"*d + "0", 2)
    base = 2**(2*d + 1) - 3*2**(d - 1)
    assert g(lo) == base - 1 and g(hi) == base + 1
    print(d, g(lo), g(hi))
print("all d = 2..11 agree with 2^(2d+1) - 3*2^(d-1) -+ 1")
