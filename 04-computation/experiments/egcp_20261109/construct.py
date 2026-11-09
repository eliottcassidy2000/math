"""End-to-end witnesses for THM-5040: build(n, x) returns an explicit even graph on exactly n vertices with d = x,
for every x in O(n), following the proof step by step (census base n <= 9 from base_witnesses.txt; padding; cones;
E_2-multi-joins; the leaf stable recursion SR; the jump/frontier family).  The main program brute-forces d and
evenness of the returned graphs.

usage: python construct.py NMAX SAMPLE   (checks all delicate values and SAMPLE random others for 10 <= n <= NMAX)
"""
import sys, random
from functools import lru_cache
from egcp_lib import (M, O, Tj, from_hex, pad, cone, E2join, SR, frontier, cocktail, dcount, is_even,
                      has_leaf, no_universal)

W = {}; Lb = {}
for line in open('base_witnesses.txt'):
    if line.startswith('#'): continue
    tag, n, d, *e = line.split(); g = from_hex(int(n), e[0] if e else '')
    if tag == 'W': W[int(d)] = g
    else: Lb[(int(n), int(d))] = g

@lru_cache(maxsize=None)
def leaf_top(m, j):
    """leaf graph without universal vertex on m vertices, d = 2^(m-1) - j, j in Tj(m)  (Theorem T)"""
    if m <= 9: return Lb[(m, 2**(m - 1) - j)]
    if j in Tj(m - 2): return SR(leaf_top(m - 2, j))
    return frontier(m, j)

@lru_cache(maxsize=None)
def leaf_up(m, c):
    """leaf graph without universal vertex on m vertices, d = 2^(m-1) + c, c in {1,3,5,7,11,13}  (Lemma 10)"""
    if m <= 9: return Lb[(m, 2**(m - 1) + c)]
    return SR(leaf_up(m - 2, c))

@lru_cache(maxsize=None)
def build(n, x):
    assert x in O(n), (n, x)
    if n <= 9: return pad(W[x], n)
    if x in O(n - 1): return pad(build(n - 1, x), n)
    if x == M(n):                                           # Lemma 6
        return cocktail(n // 2) if n % 2 == 0 else cone(cocktail((n - 1) // 2))
    if x < 2**(n - 1):                                      # Case 1
        j = 2**(n - 1) - x
        return leaf_top(n, j)                               # (x = M(n-1) was caught by O(n-1))
    D = 2**n - x
    if D <= 2**(n - 2) + 1:                                 # Case 2a
        g = 1 if n % 2 else 2
        while 2**g + n - g - 1 < D: g += 2
        dG = 2**g - D + n - g
        if dG in (3, 7, 13) and g + 2 <= n - 2:
            g += 2; dG = 2**g - D + n - g
        if dG not in (3, 7, 13):
            return E2join(build(g, dG), (n - g) // 2)
        return cone(leaf_up(n - 1, x - 3 * 2**(n - 2)))    # x = 3*2^(n-2) + {1,5,11}
    y = 2**(n - 1) - D                                      # Case 2b
    if y in (3, 7, 13): return leaf_up(n, y)
    if y in O(n - 2): return cone(build(n - 1, y))          # build(n-1, y) is padded: no universal vertex
    return cone(leaf_top(n - 1, 2**(n - 2) - y))

def check(g, n, x):
    return g[0] == n and dcount(g) == x and is_even(g)

if __name__ == "__main__":
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    per = int(sys.argv[2]) if len(sys.argv) > 2 else 30
    rng = random.Random(89)
    bad = total = 0
    for n in range(10, nmax + 1):
        new = sorted(O(n) - O(n - 1))
        delicate = [x for x in new if x < 2**(n - 1) + 2**(n - 3) or x > 2**n - 4 * n
                    or abs(x - 3 * 2**(n - 2)) < 40]
        sample = sorted(set(delicate) | set(rng.sample(new, min(per, len(new)))))
        for x in sample:
            g = build(n, x); total += 1
            if not check(g, n, x):
                bad += 1; print("FAIL n=%d x=%d" % (n, x))
        # the leaf families themselves
        for j in sorted(Tj(n)):
            g = leaf_top(n, j); total += 1
            if not (check(g, n, 2**(n - 1) - j) and has_leaf(g) and no_universal(g)):
                bad += 1; print("FAIL leaf_top", n, j)
        print("n=%d: %d witnesses for the %d values needing n vertices, plus T(%d); failures %d"
              % (n, len(sample), len(new), n, bad), flush=True)
    print("TOTAL", total, "BAD", bad)
