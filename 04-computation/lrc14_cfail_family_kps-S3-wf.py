"""
Phase 3: map the C(S)-failure FAMILY in S3 and confirm M(S)>=1/14 on all of it.

Phase 2 found exactly one C-fail in the structured boundary sweep:
  S* = [1,2,3,5,7,8,9,10,11,12,13,38,42], best_ratio 429/532, M=2/23.
This phase:
  (1) widens the boundary sweep (more drops, larger 2-large range) to collect ALL C-fails,
  (2) computes exact M on every C-fail to confirm none breaks LRC (M<1/14),
  (3) reports the min best-criterion ratio and min M over the C-fail family,
  (4) extends to 3-large clusters (k=3) to check C-fails are not confined to k=2.
Exact Fraction. A C-fail with M<1/14 would be an LRC(14) counterexample.
"""
import sys, time
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations

def log(*a): print(*a); sys.stdout.flush()
H = F(1, 14)

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r

def safe_components(A, h=H):
    iv = []
    for u in A:
        for j in range(0, u):
            c = F(j, u); a = (c - h / u) % 1; b = (c + h / u) % 1
            if a < b: iv.append((a, b))
            else: iv.append((a, F(1))); iv.append((F(0), b))
    iv.sort(); merged = []
    for a, b in iv:
        if merged and a <= merged[-1][1]: merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else: merged.append((a, b))
    safe = []; prev = F(0)
    for a, b in merged:
        if a > prev: safe.append((prev, a))
        prev = max(prev, b)
    if prev < 1: safe.append((prev, F(1)))
    return safe

def Wwidth(A):
    sc = safe_components(A)
    if not sc: return F(0)
    ws = [b - a for a, b in sc]
    if sc[0][0] == 0 and sc[-1][1] == 1 and len(sc) > 1: ws.append((sc[0][1]) + (1 - sc[-1][0]))
    return max(ws)

def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2 * k + 1, 2 * v) <= F(1, 2): C.add(F(2 * k + 1, 2 * v)); k += 1
    for i in range(len(S)):
        for j in range(i + 1, len(S)):
            for d in (S[i] + S[j], S[j] - S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2): C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C

def Mval(S):
    b = F(0)
    for t in cand(S):
        v = min(nrm(x * t) for x in S)
        if v > b: b = v
    return b

def is_covering(S): return all(any(v % q == 0 for v in S) for q in range(2, 15))
def primitive(S): return reduce(gcd, S) == 1
def best_criterion(S):
    br = F(0); bv = None
    for v in S:
        rest = [x for x in S if x != v]
        r = Wwidth(rest) * 7 * v
        if r > br: br = r; bv = v
    return br, bv
def case_of(S):
    S = sorted(S); Vmin, Vmax = S[0], S[-1]
    k = sum(1 for v in S if v > 13)
    if k <= 1: return 'S1'
    if Vmax < 13 * Vmin: return 'S2'
    return 'S3'

t0 = time.time()
cfails = []   # (S, best_ratio, M)
seen = set()
nchk = 0

log("Phase 3: k=2 boundary sweep (P = 11 small from {1..13}, two large a<b).")
# P = {1..13} minus two drops (one may be from a 'covering-supplied' position).
for drops in combinations(range(2, 14), 2):  # keep 1 -> Vmin=1
    P = [x for x in range(1, 14) if x not in drops]  # 11 elements incl 1
    if len(P) != 11: continue
    for a in range(14, 60):
        for b in range(a + 1, a + 36):
            S = tuple(sorted(P + [a, b]))
            if len(S) != 13 or S in seen: continue
            seen.add(S)
            Sl = list(S)
            if not primitive(Sl) or not is_covering(Sl): continue
            if case_of(Sl) != 'S3': continue
            nchk += 1
            br, bv = best_criterion(Sl)
            if br <= 1:
                m = Mval(Sl)
                cfails.append((Sl, br, bv, m))
    if time.time() - t0 > 320:
        log("  (time cap hit during k=2 sweep)"); break

log(f"k=2 sweep: checked {nchk} S3 sets in {time.time()-t0:.1f}s")
log(f"C(S) FAILURES (k=2): {len(cfails)}")
cfails.sort(key=lambda z: z[1])
for S, br, bv, m in cfails:
    flag = "  *** M<1/14 LRC-COUNTEREXAMPLE ***" if m < F(1, 14) else ""
    log(f"  {S}  best_ratio {br} ~{float(br):.5f} via {bv}  M={m} ~{float(m):.5f}{flag}")

# k=3 spot check on the same small-part skeleton near S*
log("\nPhase 3b: k=3 spot check (P=10 small, three large near 38..50).")
n3chk = 0; cf3 = []
for drops in combinations(range(2, 14), 3):
    P = [x for x in range(1, 14) if x not in drops]  # 10 small incl 1
    for triple in combinations(range(14, 52), 3):
        S = tuple(sorted(P + list(triple)))
        if len(S) != 13 or S in seen: continue
        seen.add(S)
        Sl = list(S)
        if not primitive(Sl) or not is_covering(Sl): continue
        if case_of(Sl) != 'S3': continue
        n3chk += 1
        br, bv = best_criterion(Sl)
        if br <= 1:
            cf3.append((Sl, br, bv, Mval(Sl)))
        if time.time() - t0 > 520: break
    if time.time() - t0 > 520: break
log(f"k=3 spot: checked {n3chk} S3 sets. C-fails: {len(cf3)}")
for S, br, bv, m in cf3[:20]:
    flag = "  *** M<1/14 ***" if m < F(1, 14) else ""
    log(f"  {S}  best_ratio ~{float(br):.5f} via {bv}  M={m} ~{float(m):.5f}{flag}")

all_cf = cfails + cf3
log("\n=== SUMMARY over C(S)-failure family ===")
log("total C-fails found:", len(all_cf))
if all_cf:
    minr = min(z[1] for z in all_cf); minm = min(z[3] for z in all_cf)
    log("min best-criterion ratio over C-fails:", minr, float(minr))
    log("min exact M over C-fails:", minm, float(minm), ">=1/14?", minm >= F(1, 14))
    lrc_breaks = [z for z in all_cf if z[3] < F(1, 14)]
    log("LRC(14) counterexamples (M<1/14):", len(lrc_breaks))
    for S, br, bv, m in lrc_breaks:
        log("   BREAK:", S, "M=", m)
log("\nDONE %.1fs" % (time.time() - t0))
