#!/usr/bin/env python3
"""r6_allshapes_finebranch_enum_deathstar_S58.py

Exhaustive fine-killer enumeration for the r=6 all-shapes tail L-bound:
  For core P (7-subset of {1..12}) and 5 killers K in [K0,332], compute the largest
  component L of S(P) \\ union(danger(k)). Verify min L > 1/2331 (=> the 'exactly one
  killer >= 333' tail is dodgeable by the sharp horn, for all-fine-killer shapes).

SCOPE / HONESTY:
- This is a PARTIAL verification: it enumerates only configs with all 5 killers in
  [K0,332]. A violation could in principle use a coarse killer (the fine restriction is
  empirically motivated by min-L search, NOT proven). A COMPLETE verification is the full
  [92,332] range = kind-pasteur's r=6 finite horn ~3.64e12 (~140 days).
- Cost: ~400s/core in pure Python (widest-arc prune) => ~88h for all 792 cores. This is a
  cron/overnight/cluster job, or an optimized-C reimplementation. Use --cores to slice.

Usage: python3 THIS.py [K0] [core_start] [core_end]
"""
import itertools, sys, time
K0 = int(sys.argv[1]) if len(sys.argv) > 1 else 283
CS = int(sys.argv[2]) if len(sys.argv) > 2 else 0
CE = int(sys.argv[3]) if len(sys.argv) > 3 else 792
THR = 1.0/2331
KR = list(range(K0, 333))

def csafe(P):
    S = [(0.0, 1.0)]
    for v in P:
        w = 1.0/(14*v); arcs = []
        for j in range(v):
            c = j/v; lo = (c-w) % 1; hi = (c+w) % 1
            if lo < hi: arcs.append((lo, hi))
            else: arcs.append((lo, 1.0)); arcs.append((0.0, hi))
        for clo, chi in sorted(arcs):
            nn = []
            for a, b in S:
                if chi <= a or clo >= b: nn.append((a, b)); continue
                if clo > a: nn.append((a, clo))
                if chi < b: nn.append((chi, b))
            S = nn
    return S

def arc_gap(a, bb, K):
    cuts = []
    for k in K:
        w = 1.0/(14*k); jlo = int(a*k); jhi = int(bb*k)+1
        for j in range(jlo, jhi+1):
            c = j/k; lo = c-w
            if lo >= bb: continue
            hi = c+w
            if hi <= a: continue
            cuts.append((a if lo < a else lo, bb if hi > bb else hi))
    cuts.sort(); cur = a; best = 0.0
    for lo, hi in cuts:
        if lo > cur and lo-cur > best: best = lo-cur
        if hi > cur: cur = hi
    if bb-cur > best: best = bb-cur
    return best

cores = [list(c) for c in itertools.combinations(range(1, 13), 7)]
gmin = 1.0; gcfg = None; viol = 0; t0 = time.time()
for ci in range(CS, min(CE, len(cores))):
    P = cores[ci]; S = csafe(P)
    arcs = sorted(S, key=lambda x: -(x[1]-x[0])); A0 = arcs[0]
    for combo in itertools.combinations(KR, 5):
        g0 = arc_gap(A0[0], A0[1], combo)
        if g0 > THR: continue                     # widest-arc gap big => L>THR, not a violation
        L = g0
        for A in arcs[1:]:
            g = arc_gap(A[0], A[1], combo)
            if g > L: L = g
        if L <= THR: viol += 1
        if L < gmin: gmin = L; gcfg = (P, combo)
    print("core %d/%d %s  elapsed %.0fs  min L=%.7f (ratio %.3f) viol=%d"
          % (ci, len(cores), P, time.time()-t0, gmin, gmin/THR, viol), flush=True)
print("DONE cores[%d,%d) killers[%d,332]: min L=%.7f (1/2331=%.7f ratio %.3f) at %s ; VIOLATIONS=%d -> %s"
      % (CS, CE, K0, gmin, THR, gmin/THR, gcfg, viol,
         "L-BOUND HOLDS (fine branch)" if viol == 0 else "VIOLATION FOUND"), flush=True)
