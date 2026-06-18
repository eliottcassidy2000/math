"""
Adversarial verification of the 'scale-invariance-AP-and-k3' closure for LRC(14).

Claim under test (AP-FAMILY THEOREM):
  Let S = {1,2,...,12,m}, m>=14 integer.
  (i)  S covering  <=>  182 | m.
  (ii) For every covering S (m=182k), M(S) >= 2/27 > 1/14, closed by two witnesses:
        Witness 1: tau = 2/27. Works unless m mod 27 in {0,13,14}.
        Witness 2: tau = (m/13)/(m+1) = 14k/(182k+1). Level = 14k/(182k+1), min 28/365 at k=2.

We re-derive every load-bearing claim EXACTLY with fractions.Fraction.
We also stress the k=3 'VERIFIED' claims and the a3/a4 'min M' claims.
"""
from fractions import Fraction as F
from math import gcd

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r

def safe_components(A, h=F(1, 14)):
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
    if sc[0][0] == 0 and sc[-1][1] == 1 and len(sc) > 1: ws.append(sc[0][1] + (1 - sc[-1][0]))
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

def is_cov(S):
    return all(any(v % q == 0 for v in S) for q in range(2, 15))

def level_at(S, tau):
    """Exact lower bound on M(S) given by global witness tau: min_v ||v*tau||."""
    return min(nrm(F(v) * tau) for v in S)

print("="*70)
print("STEP 0: Covering reduction for S={1,...,12,m}")
print("="*70)
# Claim: S={1..12,m} covering <=> 182|m
base = list(range(1, 13))
violations = []
for m in range(14, 6000):
    S = base + [m]
    cov = is_cov(S)
    claim = (m % 182 == 0)
    if cov != claim:
        violations.append((m, cov, claim))
print(f"Checked m in [14,6000). Covering<=>182|m violations: {len(violations)}")
if violations:
    print("  VIOLATIONS:", violations[:20])
# Also confirm m=182 IS covering and m=181,183 are not
for m in [181, 182, 183, 364]:
    S = base + [m]
    print(f"  m={m}: covering={is_cov(S)}  (182|m: {m%182==0})")

print()
print("="*70)
print("STEP 1: WITNESS 1  tau = 2/27")
print("="*70)
tau1 = F(2, 27)
# (a) small-part level min_{v=1..12} ||2v/27||
small_levels = {v: nrm(F(2*v, 27)) for v in range(1, 13)}
print("  ||2v/27|| for v=1..12:")
for v in range(1, 13):
    print(f"    v={v:2d}: {small_levels[v]}")
small_min = min(small_levels.values())
print(f"  min over v=1..12 = {small_min}  (claim: 2/27 = {F(2,27)})  ok={small_min==F(2,27)}")
# (b) runner m residue table mod 27 : ||2r/27|| >= 2/27 except r in {0,13,14}
print("  Residue table ||2r/27|| for r=0..26 (flag if < 2/27):")
bad_res = []
for r in range(27):
    lvl = nrm(F(2*r, 27))
    if lvl < F(2, 27):
        bad_res.append((r, lvl))
print(f"    residues with ||2r/27|| < 2/27: {bad_res}")
print(f"    claim says bad set = {{0,13,14}}; computed = {sorted(r for r,_ in bad_res)}")

print()
print("="*70)
print("STEP 2: WITNESS 2  tau = (m/13)/(m+1) = 14k/(182k+1)")
print("="*70)
# Verify closed form level = 14k/(182k+1) for covering m=182k, AND that it's a valid global witness
# Check both the algebraic claim and direct min_v ||v*tau||.
w2_fail = []
for k in range(1, 600):
    m = 182 * k
    q = m // 13  # = 14k
    tau2 = F(q, m + 1)
    lvl = level_at(base + [m], tau2)
    claimed = F(14 * k, 182 * k + 1)
    if lvl != claimed:
        w2_fail.append((k, lvl, claimed))
print(f"  k in [1,600): level == 14k/(182k+1) mismatches: {len(w2_fail)}")
if w2_fail:
    print("    first few:", w2_fail[:5])
# min over k of the level for the residue classes needing witness 2
print(f"  level at k=2 (m=364): {F(14*2,182*2+1)}  (claim 28/365 = {F(28,365)})")
print(f"  28/365 - 1/14 = {F(28,365)-F(1,14)}  (claim 27/5110 = {F(27,5110)})")
print(f"  2/27 - 1/14 = {F(2,27)-F(1,14)}  (claim 1/378 = {F(1,378)})")

print()
print("="*70)
print("STEP 3: COMBINE — which covering m need witness 2, and does W2 cover them?")
print("="*70)
# Witness 1 fails iff m mod 27 in {0,13,14}. m=182k, 182 mod 27 = ?
print(f"  182 mod 27 = {182 % 27}")
# m mod 27 = (182k) mod 27 = (20k) mod 27
# W1 fails iff 20k mod 27 in {0,13,14}
need_w2 = []
for k in range(0, 27):
    r = (182 * k) % 27
    if r in (0, 13, 14):
        need_w2.append((k, r))
print(f"  k mod 27 classes where m mod 27 in {{0,13,14}} (need W2): {need_w2}")
print(f"  claim says k mod 27 in {{0,2,25}}")
# Now the CRITICAL adversarial check: for EVERY covering m=182k up to large bound,
# at least ONE witness gives level > 1/14, and record the min over all k.
print()
print("  Full sweep: for each covering m=182k, take best witness level, find global min.")
worst = None
counterex = []
for k in range(1, 2000):
    m = 182 * k
    S = base + [m]
    l1 = level_at(S, F(2, 27)) if (m % 27) not in (0, 13, 14) else F(0)
    l2 = level_at(S, F(14*k, 182*k+1))
    best = max(l1, l2)
    if best <= F(1, 14):
        counterex.append((k, m, best))
    if worst is None or best < worst[2]:
        worst = (k, m, best)
print(f"  k in [1,2000): covering m closed by best witness > 1/14 always? {len(counterex)==0}")
print(f"  worst (min best-witness level): k={worst[0]} m={worst[1]} level={worst[2]} (={float(worst[2]):.6f})")
print(f"  1/14 = {float(F(1,14)):.6f}; floor exceeds 1/14: {worst[2] > F(1,14)}")
if counterex:
    print("  !!! COUNTEREXAMPLES (best witness <= 1/14):", counterex[:10])

print()
print("="*70)
print("STEP 3b: Cross-check true M(S) via Mval for covering m=182k, small k")
print("="*70)
mfail = []
for k in range(1, 30):
    m = 182 * k
    S = base + [m]
    Mv = Mval(S)
    if Mv < F(1, 14):
        mfail.append((k, m, Mv))
    if k <= 8:
        print(f"  k={k} m={m}: M(S)={Mv} (={float(Mv):.6f}), >=1/14: {Mv>=F(1,14)}")
print(f"  covering m=182k (k<30) with true M < 1/14: {len(mfail)}")
if mfail:
    print("    ", mfail)
