#!/usr/bin/env python3
"""Verifier (lane digits), check V1b: EVERY matrix of onemark_matrices_5.json against own direct counts on eight
marked bases of orders 3..5 (not the lane's interpolation bases) whose states span Q^4."""
import json, os, random, time
import numpy as np
from v1_onemark_matrices import direct_state
from vlin import Space

here = os.path.dirname(os.path.abspath(__file__))
t0 = time.time()
rng = random.Random(4)
s = np.zeros(4, np.int64)
bases = []
V = Space(4)
while len(bases) < 8 or V.dim() < 4:
    n = rng.choice([3, 4, 5])
    adj = np.zeros(n, np.int64)
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < 0.5:
                adj[i] |= 1 << j
            else:
                adj[j] |= 1 << i
    z = rng.randrange(n)
    direct_state(n, adj, z, s)
    st = [int(x) for x in s]
    if st in [b[3] for b in bases]:
        continue
    V.add(st)
    bases.append((n, adj, z, st))
print(f"{len(bases)} marked bases of orders 3..5, rank of their states {V.dim()}: {[b[3] for b in bases]}")
ops = json.load(open(here + "/../onemark_matrices_5.json"))
bad = 0
s2 = np.zeros(4, np.int64)
for o in ops:
    k = o["k"]; W = o["W"]; m = o["m"]
    for (n, adj, z, st) in bases:
        N = n + k
        a = np.zeros(N, np.int64)
        xmask = ((1 << k) - 1) << n
        for i in range(n):
            a[i] = int(adj[i])
            if i != z:
                a[i] |= xmask
            else:
                a[i] |= (((1 << k) - 1) & ~W) << n
        for i in range(k):
            a[n + i] = int(o["adjX"][i]) << n
            if (W >> i) & 1:
                a[n + i] |= 1 << z
        direct_state(N, a, z if m < 0 else n + m, s2)
        if [sum(o["N"][r][q] * st[q] for q in range(4)) for r in range(4)] != [int(x) for x in s2]:
            bad += 1
print(f"{len(ops)} stored matrices x {len(bases)} bases: failures {bad}   [{time.time() - t0:.1f} s]")
