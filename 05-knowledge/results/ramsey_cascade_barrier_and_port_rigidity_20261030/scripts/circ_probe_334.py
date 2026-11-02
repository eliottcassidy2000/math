#!/usr/bin/env python3
"""circ_probe_334.py -- all circulant (3,3,4;29) colourings: colour the 14 differences {1..14} of Z/29 with colours
0,1 (triangle-free) and 2 (K4-free); a cyclic colouring has edge counts 29*|D_c| and degrees 2|D_c|.  Exhaustive
DFS over 3^14 assignments with incremental monochromatic-triangle checks and a final K4 check; colours 0 and 1
are interchangeable, so we fix colour(1) in {0, 2} ... no: report all, then quotient by the swap and by the
multiplier group (Z/29)^x acting on differences.  Standard library only."""
import itertools, sys
p = 29; D = list(range(1, 15))
def dist(x): x %= p; return min(x, p - x)
def tri_free(S):
    """no triangle: no a,b in S (with repetition) such that dist(a+b) or dist(a-b) in S"""
    Sset = set(S)
    for a in S:
        for b in S:
            if dist(a + b) in Sset or dist(a - b) in Sset: return False
    return True
def k4_free(S):
    Sset = set(S)
    # vertices 0 < a < b < c with all six differences in S
    for a in range(1, p):
        if dist(a) not in Sset: continue
        for b in range(a + 1, p):
            if dist(b) not in Sset or dist(b - a) not in Sset: continue
            for c in range(b + 1, p):
                if dist(c) in Sset and dist(c - a) in Sset and dist(c - b) in Sset: return False
    return True
sols = []
def dfs(i, cls):
    if i == len(D):
        S2 = [d for d in D if cls[d] == 2]
        if k4_free(S2): sols.append(dict(cls))
        return
    d = D[i]
    for c in (0, 1, 2):
        cls[d] = c
        if c < 2:
            S = [x for x in D[:i+1] if cls[x] == c]
            if not tri_free(S): continue
        dfs(i + 1, cls)
    del cls[d]
dfs(0, {})
print('circulant (3,3,4;29) colourings (labelled by difference classes):', len(sols))
from collections import Counter
splits = Counter()
for s in sols:
    sizes = tuple(sum(1 for d in D if s[d] == c) for c in range(3))
    splits[sizes] += 1
print('class-size splits (|D0|,|D1|,|D2|) -> edge counts 29*size, degrees 2*size:')
for k, v in sorted(splits.items()): print('  ', k, 'edges', tuple(29 * x for x in k), 'degrees', tuple(2 * x for x in k), 'count', v)
# quotient by multipliers and by the colour swap
units = [u for u in range(1, p) if True]
canon = set()
for s in sols:
    best = None
    for u in range(1, p):
        for swap in (False, True):
            img = tuple(((1 - s[dist(u * d)]) if (swap and s[dist(u * d)] < 2) else s[dist(u * d)]) for d in D)
            if best is None or img < best: best = img
    canon.add(best)
print('up to multipliers and the colour-0/1 swap:', len(canon), 'classes')
for c in sorted(canon): print('  ', ''.join(map(str, c)), 'split', tuple(c.count(x) for x in range(3)))
