"""
Phase 2: targeted adversarial drive at the WEAK SPOTS found in phase 1.

Phase 1 found: 0 C(S) failures over 9150 S3 sets; worst best-ratio = 207/176 ~ 1.176
at [1,2,3,5,7,8,9,10,11,12,13,23,28] (a near-S2/S3 boundary set where Vmax fails).

Targets here:
  (T1) Drive the BEST criterion ratio as close to 1 (or below) as possible.
       Focus on S2/S3 boundary (Vmax ~ 13*Vmin) and degenerate clusters (k=2).
  (T2) Exhaustively classify the boundary family {1..13 minus drops} U {two large} with
       Vmax just at/above 13*Vmin, computing best_criterion + exact M.
  (T3) Confirm: for EVERY sampled S3 set with smallest best-ratios, M(S) >= 1/14 exactly.
  (T4) Verify the snippet's G2 boundary set [1,2,3,5,7,8,9,10,11,12,13,27,28] exactly:
       via-Vmax(28) ratio, via-27 ratio, and M.
Exact Fraction. One C(S) failure with M<1/14 => counterexample.
"""
import sys, random, time
from fractions import Fraction as F
from math import gcd
from functools import reduce

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

# ---------------------------------------------------------------------------
log("=" * 70); log("T4: exact audit of the G2 boundary set [1,2,3,5,7,8,9,10,11,12,13,27,28]")
log("=" * 70)
G2 = [1, 2, 3, 5, 7, 8, 9, 10, 11, 12, 13, 27, 28]
log("covering:", is_covering(G2), " primitive:", primitive(G2), " case:", case_of(G2))
for v in [28, 27]:
    rest = [x for x in G2 if x != v]
    r = Wwidth(rest) * 7 * v
    log(f"  via-{v}: W(S\\{v})*7*{v} = {r} ~ {float(r)}")
br, bv = best_criterion(G2)
log("  best_criterion:", br, float(br), "via", bv)
log("  M(G2) =", Mval(G2), float(Mval(G2)), " >=1/14?", Mval(G2) >= F(1, 14))

log("\n" + "=" * 70); log("T2: boundary family sweep. P = drop-d core, two large near 13*Vmin")
log("=" * 70)
# Vmin will be 1 (P contains 1). 13*Vmin=13, so S3 needs Vmax>=13. Use clusters at 14..60.
worst = []  # (ratio, M, S)
t0 = time.time()
random.seed(2024)
seen = set()
nchk = 0
for d in range(2, 14):  # drop d from {1..13} (keep 1 so Vmin=1); P has 12 elts -> need to remove some for 2 large
    P11 = [x for x in range(1, 14) if x != d]  # 12 elements, includes 1
    # we need room for 2 large: drop one more small to make P size 11
    for d2 in range(1, 14):
        if d2 == d: continue
        if d2 == 1: continue  # keep Vmin=1
        P = sorted(set(P11) - {d2})  # 11 small
        for a in range(14, 56):
            for b in range(a + 1, a + 30):
                S = sorted(set(P) | {a, b})
                if len(S) != 13: continue
                if tuple(S) in seen: continue
                seen.add(tuple(S))
                if not primitive(S) or not is_covering(S): continue
                if case_of(S) != 'S3': continue
                br, bv = best_criterion(S)
                nchk += 1
                if br <= F(13, 10):  # near-threshold: record and get exact M
                    m = Mval(S)
                    worst.append((br, m, S, bv))
        if time.time() - t0 > 300: break
    if time.time() - t0 > 300: break

worst.sort(key=lambda x: x[0])
log(f"checked {nchk} boundary S3 sets in {time.time()-t0:.1f}s")
log("smallest best-criterion ratios (near threshold 1):")
anyfail = False
for br, m, S, bv in worst[:30]:
    flag = "  <<< C-FAIL" if br <= 1 else ""
    mflag = "  <<< M<1/14 COUNTEREXAMPLE" if m < F(1, 14) else ""
    if br <= 1: anyfail = True
    log(f"  ratio {float(br):.4f} ({br}) via {bv}  M={float(m):.5f} ({m}){flag}{mflag}")
log("Any C(S) failure (best_ratio<=1) in boundary sweep:", anyfail)
mn = min((m for _, m, _, _ in worst), default=None)
if mn is not None:
    log("min exact M among near-threshold boundary sets:", mn, float(mn), ">=1/14?", mn >= F(1,14))

log("\n" + "=" * 70); log("T1: random drive to minimize best-criterion ratio across all S3 regimes")
log("=" * 70)
random.seed(555)
gmin_ratio = None; gmin_set = None; gmin_M = None
nrand = 0
while time.time() - t0 < 470:
    # adversarial: low-W small part (drop-6-like cores) + small-Vmax tight cluster (k=2 or 3)
    Vmin = random.choice([1, 1, 1, 2, 3])
    # small part avoiding to lower W
    psize = random.randint(8, 11)
    P = sorted(random.sample(range(Vmin, 14), min(psize, 14 - Vmin)))
    if Vmin not in P: P = sorted(set(P) | {Vmin})
    Lsize = 13 - len(set(P))
    if Lsize < 2: continue
    lo = 13 * Vmin
    V0 = lo + random.randint(0, 25)
    spread = random.choice([14, 15, 16, 18, 21, 24, 28])
    if spread + 1 < Lsize: continue
    offs = sorted(random.sample(range(0, spread + 1), Lsize))
    L = [V0 + o for o in offs]
    S = sorted(set(P) | set(L))
    pool = [x for x in range(1, 14) if x not in S]
    while len(S) < 13 and pool: S = sorted(set(S) | {pool.pop()})
    if len(S) != 13 or not primitive(S) or not is_covering(S): continue
    if case_of(S) != 'S3': continue
    nrand += 1
    br, bv = best_criterion(S)
    if gmin_ratio is None or br < gmin_ratio:
        gmin_ratio = br; gmin_set = S; gmin_M = Mval(S)
        log(f"  new min best-ratio {float(br):.4f} ({br}) via {bv}  M={float(gmin_M):.5f}  {S}")

log(f"\nT1 random: {nrand} S3 sets. Global min best-criterion ratio = {gmin_ratio} ~ {float(gmin_ratio) if gmin_ratio else None}")
log("  at", gmin_set, " M =", gmin_M, float(gmin_M) if gmin_M else None, ">=1/14?", (gmin_M>=F(1,14)) if gmin_M else None)
log("\nDONE total %.1fs" % (time.time() - t0))
