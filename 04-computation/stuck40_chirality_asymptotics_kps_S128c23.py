#!/usr/bin/env python3
"""stuck40_chirality_asymptotics_kps_S128c23.py -- kind-pasteur S128 cont.23.
(A) classify the 40 stuck n=5 tournaments (classes, scores); (B) the A5/2I chirality census at n=5
(every Aut is odd-order => even => class = 2 chiral A5-orbits; refine SC by even/odd op-witness);
(C) false-peak asymptotics: exhaustive n=7."""
import sys
from math import comb
from itertools import permutations
from collections import Counter, defaultdict
sys.stdout.reconfigure(line_buffering=True)
def canon_and_scores(B, n, perms):
    pairs = [(u, v) for u in range(n) for v in range(u + 1, n)]
    pidx = {p: k for k, p in enumerate(pairs)}
    best = None
    for pm in perms:
        x = 0
        for (i, j) in pairs:
            pi, pj = pm[i], pm[j]
            lo, hi = (pi, pj) if pi < pj else (pj, pi)
            bit = B[i][j] if pm[i] == lo else (not B[i][j])
            x |= (1 if bit else 0) << pidx[(lo, hi)]
        if best is None or x < best:
            best = x
    return best
# (A)+(B)
n = 5
perms5 = list(permutations(range(n)))
evens5 = [p for p in perms5 if sum(1 for i in range(5) for j in range(i+1,5) if p[i]>p[j]) % 2 == 0]
pairs = [(u, v) for u in range(n) for v in range(u + 1, n)]
m = comb(n, 2)
c3max = 5
stuck = Counter()
classes_all = {}
for mask in range(1 << m):
    B = [[False]*n for _ in range(n)]
    for k, (u, v) in enumerate(pairs):
        if (mask >> k) & 1: B[u][v] = True
        else: B[v][u] = True
    d = [sum(r) for r in B]
    c3 = comb(n,3) - sum(comb(x,2) for x in d)
    key = canon_and_scores(B, n, perms5)
    classes_all.setdefault(key, (tuple(sorted(d)), c3))
    if c3 < c3max and not any(B[u][v] and d[u]==d[v]+2 for u in range(n) for v in range(n) if u!=v):
        stuck[(key, tuple(sorted(d)), c3)] += 1
print("(A) the 40 stuck n=5 tournaments by class:")
for (key, sc, c3), cnt in stuck.items():
    print("   class canon=%d scores=%s c3=%d : %d labeled members" % (key, sc, c3, cnt))
print("(B) chirality census n=5 (op-witness parity per class):")
chir = Counter()
for key, (sc, c3) in sorted(classes_all.items()):
    # reconstruct a rep from canon key
    B = [[False]*n for _ in range(n)]
    pidx = {p: k for k, p in enumerate(pairs)}
    for (u, v) in pairs:
        if (key >> pidx[(u,v)]) & 1: B[u][v] = True
        else: B[v][u] = True
    even_w = odd_w = False
    for pm in perms5:
        ok = all(B[i][j] == (not B[[pm.index(x) for x in range(n)][0]][0] if False else B[pm.index(0)][0]) for i in [] ) # placeholder
    # direct: pi is op-witness iff for all u!=v: B[u][v] == B_op[pi(u)][pi(v)] = B[pi(v)][pi(u)]
    for pm in perms5:
        if all((B[u][v] == B[pm[v]][pm[u]]) for u in range(n) for v in range(n) if u != v):
            par = sum(1 for i in range(5) for j in range(i+1,5) if pm[i]>pm[j]) % 2
            if par == 0: even_w = True
            else: odd_w = True
    tag = ("SC-even" if even_w else "") + ("+" if even_w and odd_w else "") + ("SC-odd" if odd_w else "") or "NS"
    chir[tag if (even_w or odd_w) else "NS"] += 1
print("   ", dict(chir), " (12 classes total; SC split by op-witness parity = the 2I/chirality refinement)")
# (C)
print("(C) false peaks exhaustive n=7:")
n7 = 7
pairs7 = [(u, v) for u in range(n7) for v in range(u + 1, n7)]
m7 = comb(n7, 2)
c3max7 = 14
cnt7 = 0; tot_below = 0
for mask in range(1 << m7):
    d = [0]*n7
    for k, (u, v) in enumerate(pairs7):
        if (mask >> k) & 1: d[u] += 1
        else: d[v] += 1
    c3 = comb(n7,3) - sum(comb(x,2) for x in d)
    if c3 == c3max7: continue
    tot_below += 1
    up = False
    for k, (u, v) in enumerate(pairs7):
        if (mask >> k) & 1:
            if d[u] == d[v] + 2: up = True; break
        else:
            if d[v] == d[u] + 2: up = True; break
    if not up: cnt7 += 1
print("   n=7: false peaks %d of %d below-max (fraction %.5f); sequence 0, 40, 368, %d; fractions %.4f, %.4f, %.5f"
      % (cnt7, tot_below, cnt7/tot_below, cnt7, 40/1000, 368/30128, cnt7/tot_below))
print("DONE")
