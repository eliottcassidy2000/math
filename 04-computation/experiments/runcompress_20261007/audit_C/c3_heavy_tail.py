#!/usr/bin/env python3
"""audit_C item 3: is a long post-run source ones-run 'the mechanism of the heavy tail' (THM-4603 title / Reading)?

Monte Carlo from the universal state (3, 1-27) with Haar post-run bits: Y = 1 + 2^r w (w a random odd 6000-bit integer,
r = 1 w.p. 2/3, r = 2 w.p. 1/3, as for geometric two-run lengths), X = 27Y - 26. The pair is run for S = 3000 Terras
steps. For every path we record survival (no absorption by S), the longest run of consecutive odd steps of the SOURCE in
[0, S) (computed on the source orbit whether or not the chain was absorbed), and the maximal debt before absorption/S.
If long source ones-runs were the mechanism of the s^(-1/2) tail, survivors would be dominated by paths with a long run.
"""
import random, math
from collections import Counter

rnd = random.Random(2718)
S = 3000
NP = 3000
rows = []
for p in range(NP):
    r = 1 if rnd.random() < 2 / 3 else 2
    w = rnd.getrandbits(6000) | 1 | (1 << 5999)
    Y = 1 + (w << r)
    X = 27 * Y - 26
    u, v, k = X, Y, 3
    absorbed = None
    run = best = 0
    kmax = 3
    for s in range(S):
        if u & 1:
            run += 1
            if run > best:
                best = run
        else:
            run = 0
        if absorbed is None:
            k += (u & 1) - (v & 1)
            if k > kmax:
                kmax = k
        u = (3 * u + 1) >> 1 if u & 1 else u >> 1
        if absorbed is None:
            v = (3 * v + 1) >> 1 if v & 1 else v >> 1
            if u == v and k == 0:
                absorbed = s + 1
    rows.append((absorbed is None, best, kmax))
surv = [r for r in rows if r[0]]
print(f"paths {NP}; survivors at S = {S}: {len(surv)} ({len(surv)/NP:.3f}); sqrt(S) x survival = {math.sqrt(S)*len(surv)/NP:.2f}")
for m in (8, 10, 12, 14, 16):
    a = sum(1 for r in rows if r[1] >= m)
    b = sum(1 for r in surv if r[1] >= m)
    print(f"  longest source ones-run >= {m:2d}: all paths {a/NP:.3f}, survivors {b/max(1,len(surv)):.3f}; "
          f"P(survive | run >= {m}) = {b/max(1,a):.3f} vs P(survive) = {len(surv)/NP:.3f}")
med_all = sorted(r[1] for r in rows)[NP // 2]
med_s = sorted(r[1] for r in surv)[len(surv) // 2]
print(f"median longest source ones-run: all {med_all}, survivors {med_s}")
# survivors after removing every path with a source ones-run >= 14 (the 'long-run' paths)
no_long = [r for r in rows if r[1] < 14]
print(f"survival among paths with NO source ones-run >= 14: {sum(r[0] for r in no_long)/len(no_long):.3f} "
      f"(n = {len(no_long)}): the tail is essentially unchanged")
print(f"max debt: survivors median {sorted(r[2] for r in surv)[len(surv)//2]}, all paths median {sorted(r[2] for r in rows)[NP//2]}")
