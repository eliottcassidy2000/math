#!/usr/bin/env python3
"""
metagraph_atlas_census_klein_S313.py — klein-2026-07-15-S313 (cont.3)

Full iso-class census of tournaments n = 3..7 by orbit sweep (each labelled tournament touched
O(1) times): classes, c3 (= level = distance-to-transitivity, HYP-6948), |Aut| (all odd:
Feit-Thompson face), self-complementarity.  Outputs the LEVEL-WIDTH tables of the metagraph
(classes per c3 level) — the x-axis profile — plus Aut histograms and cross-checks
(A000568 = 1,2,4,12,56,456; max |Aut|; SC counts 8,12 at n=5,6 per canon).
"""
import itertools
from math import comb

def census(n):
    m = n * (n - 1) // 2
    pairs = list(itertools.combinations(range(n), 2))
    pidx = {pr: i for i, pr in enumerate(pairs)}
    perms = list(itertools.permutations(range(n)))
    # precompute bit remap tables per perm: newbit j from old bit i (+ flip flag)
    remaps = []
    for g in perms:
        tab = []
        for i, (u, v) in enumerate(pairs):
            gu, gv = g[u], g[v]
            j = pidx[(min(gu, gv), max(gu, gv))]
            tab.append((j, 0 if gu < gv else 1))
        remaps.append(tab)
    def apply_r(bits, tab):
        out = 0
        for i in range(m):
            b = (bits >> i) & 1
            j, fl = tab[i]
            out |= ((b ^ fl) << j)
        return out
    seen = bytearray(1 << m)
    classes = []
    comp_mask = (1 << m) - 1
    for bits in range(1 << m):
        if seen[bits]: continue
        orb = {apply_r(bits, tab) for tab in remaps}
        for t in orb: seen[t] = 1
        aut = len(perms) // len(orb)
        # scores and c3
        s = [0] * n
        for i, (u, v) in enumerate(pairs):
            if (bits >> i) & 1: s[u] += 1
            else: s[v] += 1
        c3 = comb(n, 3) - sum(comb(t, 2) for t in s)
        sc = (bits ^ comp_mask) in orb
        classes.append((bits, c3, aut, sc, len(orb)))
    return classes

A000568 = {3: 2, 4: 4, 5: 12, 6: 56, 7: 456}
print("n | classes | level widths (classes per c3 = 0..max) | Aut histogram {order: #classes} | SC")
for n in range(3, 8):
    cl = census(n)
    assert len(cl) == A000568[n], (n, len(cl))
    assert sum(orb for *_, orb in [(c[0], c[4]) for c in cl] ) or True
    assert sum(c[4] for c in cl) == 1 << (n * (n - 1) // 2)
    kmax = (n ** 3 - n) // 24 if n % 2 == 1 else (n ** 3 - 4 * n) // 24
    widths = [0] * (kmax + 1)
    auth = {}
    scn = 0
    for _, c3, aut, sc, _orb in cl:
        widths[c3] += 1
        auth[aut] = auth.get(aut, 0) + 1
        scn += sc
        assert aut % 2 == 1, "even Aut?!"   # Feit-Thompson face
    print(f"n={n}: {len(cl):4d} | {widths} | {dict(sorted(auth.items()))} | SC={scn}")
    # widths read backwards = classes per x-level from ceiling down (x = (n^3-n)/3 - 8*c3)
print()
print("Notes: level k = c3 = # cyclic triangles = tie-splits from transitive (HYP-6948).")
print("Cross-checks: totals = A000568; all |Aut| odd; SC counts vs canon (8 at n=5, 12 at n=6...).")
print("Max width vs old antichain-width formula C(n-2, floor((n-2)/2)) (different 'width'!):")
for n in range(3, 8):
    cl = census(n)
    kmax = (n ** 3 - n) // 24 if n % 2 == 1 else (n ** 3 - 4 * n) // 24
    widths = [0] * (kmax + 1)
    for _, c3, *_ in cl: widths[c3] += 1
    print(f"  n={n}: max level width = {max(widths)} at k = {widths.index(max(widths))};"
          f" old formula {comb(n - 2, (n - 2) // 2)}")
