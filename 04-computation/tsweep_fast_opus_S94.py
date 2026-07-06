#!/usr/bin/env python3
"""
opus-2026-07-06-S94 -- HYP-4226(b): the T-sweep constant, fast version.

Per-class ACTUAL desert measurement (float arithmetic with exact-recheck flagging):
for |B| in (4, 3): over all bases B, cluster position-subsets C of the complement
(|C| >= 7), and k-tuning samples (consecutive-shape + randoms), measure at W = 1200:
  * T(config) = total measure of long components (> 5/W) of the cluster's tooth union;
  * the largest clear sub-segment of good(B, 1/13+eps) minus those deserts.
Report worst clear component per |B| (the COMP constant feeding the per-family measure
architecture) and flag any config within 1e-7 of a threshold for exact recheck.
"""
import random, time

RHO = 1.0 / 13.0
EPS = 1.0 / 500.0
W = 1200

def union_components(vals, thresh):
    ivs = []
    for w in vals:
        inv = 1.0 / w
        m = 0
        while True:
            a = (m - RHO) * inv
            b = (m + RHO) * inv
            if a > 1.0: break
            if b > 0.0: ivs.append((max(a, 0.0), min(b, 1.0)))
            m += 1
    ivs.sort()
    out = []
    ca, cb = ivs[0]
    for a, b in ivs[1:]:
        if a <= cb:
            if b > cb: cb = b
        else:
            out.append((ca, cb)); ca, cb = a, b
    out.append((ca, cb))
    return [(a, b) for a, b in out if b - a > thresh]

def good_components(B):
    beta = RHO + EPS
    ivs = []
    for w in B:
        for m in range(0, w + 1):
            a, b = (m - beta) / w, (m + beta) / w
            if b > 0 and a < 1: ivs.append((max(a, 0.0), min(b, 1.0)))
    ivs.sort()
    merged = []
    for a, b in ivs:
        if merged and a <= merged[-1][1]:
            if b > merged[-1][1]: merged[-1] = (merged[-1][0], b)
        else: merged.append((a, b))
    good, prev = [], 0.0
    for a, b in merged:
        if a > prev: good.append((prev, a))
        prev = max(prev, b)
    if prev < 1.0: good.append((prev, 1.0))
    return good

from itertools import combinations
random.seed(95)
t0 = time.time()
for nB in (4, 3):
    worst_clear, worst_cfg, worst_T, nconf, flagged = None, None, 0.0, 0, 0
    for B in combinations(range(1, 13), nB):
        gcomp = good_components(list(B))
        comp = [r for r in range(1, 13) if r not in B]
        for csize in range(7, len(comp) + 1):
            for C in combinations(comp, csize):
                tunings = [[W + (r - C[0]) for r in C]]
                for _ in range(24):
                    tunings.append(sorted(set(W + (r - C[0]) + 13 * random.randint(0, 6) for r in C)))
                for vals in tunings:
                    if len(vals) < csize: continue
                    nconf += 1
                    deserts = union_components(vals, 5.0 / W)
                    T = sum(b - a for a, b in deserts)
                    if T > worst_T: worst_T = T
                    best = 0.0
                    for a, b in gcomp:
                        segs = [(a, b)]
                        for da, db in deserts:
                            nsegs = []
                            for x, y in segs:
                                if x < da: nsegs.append((x, min(y, da)))
                                if y > db: nsegs.append((max(x, db), y))
                            segs = [(x, y) for x, y in nsegs if y > x]
                            if not segs: break
                        for x, y in segs:
                            if y - x > best: best = y - x
                    if worst_clear is None or best < worst_clear:
                        worst_clear, worst_cfg = best, (B, C, len(vals))
                    if abs(best) < 1e-7 or (worst_clear and abs(best - worst_clear) < 1e-7):
                        flagged += 1
    print(f"|B|={nB}: {nconf} configs, {time.time()-t0:.0f}s: worst T = {worst_T:.4f}; "
          f"worst clear = {worst_clear:.6f} at {worst_cfg}; entry bar = "
          f"{2.0/worst_clear if worst_clear and worst_clear > 0 else float('inf'):.0f}; "
          f"near-threshold flags: {flagged}", flush=True)
print("TSWEEP DONE")
