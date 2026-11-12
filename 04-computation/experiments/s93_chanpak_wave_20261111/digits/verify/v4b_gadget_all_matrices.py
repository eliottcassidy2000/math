#!/usr/bin/env python3
"""Verifier (lane digits), check V4b: EVERY gadget matrix of gadget_matrices_{3,4,6}.json against own direct counts.
Six fixed pairs (B, O) whose states span Q^4; for each gadget entry (X, P, G): state of (B+, P) by direct count = G s.
A 4 x 4 matrix is determined by its action on four independent states, so each stored matrix is certified."""
import json, os, sys, time
import numpy as np
from v4_gadgets import pair_state_direct, closing
from vlin import Space

here = os.path.dirname(os.path.abspath(__file__))
t0 = time.time()
# bases: (n, adjacency rows, O)
cand = []
for n in range(1, 4):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for bits in range(1 << len(pairs)):
        a = [0] * n
        for t, (i, j) in enumerate(pairs):
            if (bits >> t) & 1:
                a[i] |= 1 << j
            else:
                a[j] |= 1 << i
        for O in range(1, 1 << n):
            cand.append((n, np.array(a, np.int64), O))
s = np.zeros(4, np.int64)
V = Space(4)
bases = []
for (n, a, O) in cand:
    pair_state_direct(n, a, O, s)
    st = [int(x) for x in s]
    if V.add(st) or len(bases) < 6:
        bases.append((n, a, O, st))
    if V.dim() == 4 and len(bases) >= 6:
        break
print(f"{len(bases)} test pairs (B, O), B of order <= 3; rank of their states: {V.dim()}; states: {[b[3] for b in bases]}")
s2 = np.zeros(4, np.int64)
for name in ("gadget_matrices_3.json", "gadget_matrices_4.json", "gadget_matrices_6.json"):
    gad = json.load(open(f"{here}/../{name}"))
    bad = 0
    for g in gad:
        adjX = np.array(g["adj"], np.int64)
        for (n, a, O, st) in bases:
            nn, aa = closing(n, a, O, g["n"], adjX)
            pair_state_direct(nn, aa, g["P"] << n, s2)
            pred = [sum(g["G"][r][q] * st[q] for q in range(4)) for r in range(4)]
            if pred != [int(x) for x in s2]:
                bad += 1
    print(f"{name}: {len(gad)} matrices x {len(bases)} bases, failures {bad}   [{time.time() - t0:.1f} s]", flush=True)
