#!/usr/bin/env python3
"""Verifier (lane digits), check V2b: EXACT (no LP) reason why the 52 translation class spaces C(N; 1) on which the
offset is not identically zero contain no actual state: on each of them g = c B with c > 1 or g = B = 0 fails B >= g >= 1."""
import json, sys
from collections import Counter
from fractions import Fraction
from vlin import class_space, int_eigs
from v2_onemark_pairs import ann_basis

ops = json.load(open(sys.path[0] + "/../onemark_matrices_5.json"))
L = [1, 0, 0, 0]
res = Counter()
ratios = Counter()
for o in ops:
    G = o["N"]
    if 1 not in int_eigs(G):
        continue
    V, lam = class_space(L, [G], [1])
    if V.dim() == 4 or V.contains(lam[0]):
        continue
    basis = ann_basis(V, 4)
    proj = [(b[2], b[3]) for b in basis if (b[2], b[3]) != (0, 0)]
    # rank of the projection to (B, g)
    rank2 = any(p[0] * q[1] - p[1] * q[0] != 0 for p in proj for q in proj)
    if rank2:
        res["projection to (B, g) is 2-dimensional: NOT excluded by B >= g >= 1"] += 1
    elif not proj:
        res["B = g = 0 on C: excluded"] += 1
    else:
        Bc, gc = proj[0]
        if gc == 0:
            res["g = 0 on C: excluded"] += 1
        else:
            r = Fraction(Bc, gc)
            ratios[r] += 1
            res["B = r g on C with r < 1: excluded" if r < 1 else "B = r g, r >= 1: NOT excluded"] += 1
print("translation class spaces C(N; 1) with offset not identically 0:", sum(res.values()))
for k, v in sorted(res.items()):
    print("  ", k, ":", v)
print("   ratios B/g on C:", sorted((str(r), c) for r, c in ratios.items()))
