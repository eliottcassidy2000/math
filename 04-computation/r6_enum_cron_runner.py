#!/usr/bin/env python3
"""r6_enum_cron_runner.py -- resumable, time-budgeted chunk runner for the r=6 all-shapes
fine-branch L-bound enumeration (death-star-S58). Designed to be driven by cron.

Each invocation:
  - reads progress from  05-knowledge/results/r6_finebranch_enum.progress  (JSON: next core, gmin, viol)
  - processes cores from `next` onward until a wall-clock BUDGET is exceeded (at a core boundary)
  - appends per-core results to  05-knowledge/results/r6_finebranch_enum.out
  - rewrites the progress file
When `next` reaches 792 it writes a COMPLETE line and does nothing further.

SCOPE (unchanged, honest): PARTIAL verification -- 5 killers in [K0,332] only; a violation could
in principle use a coarse killer. Complete verification is the full [92,332] range (~3.64e12).
"""
import itertools, json, os, sys, time

REPO = "/home/claude/math"
OUT  = os.path.join(REPO, "05-knowledge/results/r6_finebranch_enum.out")
PROG = os.path.join(REPO, "05-knowledge/results/r6_finebranch_enum.progress")
K0   = 283
BUDGET = int(sys.argv[1]) if len(sys.argv) > 1 else 1500   # seconds per invocation
THR  = 1.0/2331
KR   = list(range(K0, 333))
NCORES = 792

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

# load progress (robust: any missing/empty/corrupt file => fresh start)
st = {"next": 0, "gmin": 1.0, "viol": 0, "worst": None}
try:
    with open(PROG) as _f:
        _d = json.load(_f)
    if isinstance(_d, dict) and "next" in _d:
        st = _d
except Exception:
    pass

if st["next"] >= NCORES:
    with open(OUT, "a") as f:
        f.write("COMPLETE: all %d cores done. global min L=%.7f (ratio %.3f) worst=%s total_viol=%d -> %s\n"
                % (NCORES, st["gmin"], st["gmin"]/THR, st["worst"], st["viol"],
                   "L-BOUND HOLDS (fine branch)" if st["viol"] == 0 else "VIOLATION"))
    print("COMPLETE"); sys.exit(0)

t0 = time.time(); ci = st["next"]; lines = []
while ci < NCORES and time.time()-t0 < BUDGET:
    P = cores[ci]; S = csafe(P)
    arcs = sorted(S, key=lambda x: -(x[1]-x[0])); A0 = arcs[0]
    cmin = 1.0
    for combo in itertools.combinations(KR, 5):
        g0 = arc_gap(A0[0], A0[1], combo)
        if g0 > THR: continue
        L = g0
        for A in arcs[1:]:
            g = arc_gap(A[0], A[1], combo)
            if g > L: L = g
        if L <= THR: st["viol"] += 1
        if L < cmin: cmin = L
        if L < st["gmin"]: st["gmin"] = L; st["worst"] = [P, list(combo)]
    lines.append("core %d %s  min_widest_pass_L=%.7f  running_min=%.7f (ratio %.3f)  viol=%d"
                 % (ci, P, cmin, st["gmin"], st["gmin"]/THR, st["viol"]))
    ci += 1

st["next"] = ci
with open(OUT, "a") as f:
    for ln in lines: f.write(ln + "\n")
    f.write("-- chunk end: reached core %d/%d in %.0fs, %d cores this run --\n"
            % (ci, NCORES, time.time()-t0, len(lines)))
json.dump(st, open(PROG, "w"))
print("processed to core %d/%d (%d cores, %.0fs); running min L ratio %.3f, viol=%d"
      % (ci, NCORES, len(lines), time.time()-t0, st["gmin"]/THR, st["viol"]))
