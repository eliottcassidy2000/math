"""
Adversarial verification of the "bounded-speed-reduction for residual case S3 of LRC(14)" claim.

The claim is NEGATIVE (angle cannot close S3) + PARTIAL. We adversarially verify and HUNT for:
  (a) an S3 covering 13-set where the criterion C(S) FAILS  (no v with W(S\v)*7v>1),
  (b) an S3 covering 13-set with M(S) < 1/14  (a genuine LRC(14) counterexample),
  (c) a failure of the claimed dominant-large sub-regime closure
      (12 non-max runners <=13 and Vmax>=53 => via-Vmax fires C(S)),
  (d) sanity on R1 (scaling invariance) and R2 (W not bounded below).

Default skeptical. Exact Fraction throughout. One genuine violation => holds=false.
"""
import sys, random, time
from fractions import Fraction as F
from math import gcd
from functools import reduce

def log(*a):
    print(*a); sys.stdout.flush()

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
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
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
    if sc[0][0] == 0 and sc[-1][1] == 1 and len(sc) > 1:
        ws.append((sc[0][1]) + (1 - sc[-1][0]))
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

def is_covering(S):
    return all(any(v % q == 0 for v in S) for q in range(2, 15))

def primitive(S):
    return reduce(gcd, S) == 1

def criterion_ratio_via(S, v):
    rest = [x for x in S if x != v]
    return Wwidth(rest) * 7 * v

def best_criterion(S):
    br = F(0); bv = None
    for v in S:
        r = criterion_ratio_via(S, v)
        if r > br: br = r; bv = v
    return br, bv

def case_of(S):
    S = sorted(S); Vmin, Vmax = S[0], S[-1]
    k = sum(1 for v in S if v > 13)
    if k <= 1: return 'S1'
    if Vmax < 13 * Vmin: return 'S2'
    return 'S3'

# ---------------------------------------------------------------------------
log("=" * 70); log("R1: scaling invariance of W(S\\{Vmax})*7*Vmax under S->g*S")
log("=" * 70)
random.seed(1); r1_ok = True
for _ in range(8):
    base = sorted(random.sample(range(1, 60), 13))
    if not primitive(base): continue
    g = random.randint(2, 7)
    rS = criterion_ratio_via(base, max(base))
    rgS = criterion_ratio_via([g * x for x in base], max([g * x for x in base]))
    if rS != rgS:
        r1_ok = False; log("R1 VIOLATION", base, g, rS, rgS)
log("R1 holds on 8 pairs:", r1_ok)

log("=" * 70); log("R2: W({1..12}) scaling")
log("=" * 70)
Wbase = Wwidth(list(range(1, 13)))
log("W({1..12}) =", Wbase)
for t in (1, 5, 20, 50):
    A = [t * x for x in range(1, 13)]; w = Wwidth(A)
    log(f"  t={t}: W={w}  W*t={w*t}")

log("=" * 70); log("Dominant-large threshold: drop-one cores of {1..13}")
log("=" * 70)
minw = None; minc = None
for drop in range(1, 14):
    sub = [x for x in range(1, 14) if x != drop]
    w = Wwidth(sub)
    cov = all(any(v % q == 0 for v in sub) for q in range(2, 14))
    if minw is None or w < minw: minw = w; minc = drop
    log(f"  drop {drop}: W={w} ~{float(w):.6f} covers2..13={cov}")
log("min-W drop-one core: drop", minc, "W=", minw, "=5/1848?", minw == F(5, 1848))
log("=> dominant-large threshold Vmax >", 1 / (7 * minw), "~", float(1 / (7 * minw)))

# ---------------------------------------------------------------------------
log("=" * 70); log("ADVERSARIAL HUNT over S3 covering 13-sets")
log("=" * 70)

worst_ratio = None; worst_ratio_set = None
fail_C = []; domlarge_fail = []
suspicious = []   # sets with via-Vmax<=1, to inspect best_criterion and M
nS3 = 0; t0 = time.time()

def consider(S):
    """Fast path: check via-Vmax only; if it fails, do full best_criterion."""
    global worst_ratio, worst_ratio_set, nS3
    S = sorted(set(S))
    if len(S) != 13 or not primitive(S) or not is_covering(S): return
    if case_of(S) != 'S3': return
    nS3 += 1
    Vmax = max(S)
    rmax = criterion_ratio_via(S, Vmax)
    non_max = [v for v in S if v != Vmax]
    domlarge = all(v <= 13 for v in non_max) and Vmax >= 53
    if domlarge and rmax <= 1:
        domlarge_fail.append((S[:], rmax))
    if rmax > 1:
        # C(S) holds via Vmax; track ratio
        if worst_ratio is None or rmax < worst_ratio:
            # only update worst if it's the true best for this set (Vmax is usually best, but
            # to be safe for the 'worst best-ratio' metric we compute best only when rmax small)
            pass
        return
    # via-Vmax failed: compute full best_criterion
    br, bv = best_criterion(S)
    if worst_ratio is None or br < worst_ratio:
        worst_ratio = br; worst_ratio_set = S[:]
    if br <= 1:
        fail_C.append((S[:], br, bv))
    suspicious.append(S[:])

# ---- Strategy A: deterministic dominant-large stress. ----
# 12 non-max all <=13 (a drop-one core, Vmin=1), Vmax a multiple supplying covering.
# Need covering: {1..13}\{drop} already covers most q; the single large Vmax must cover
# whatever q is lost, and we need k>=2 (>=2 speeds >13) for S3... but drop-one core has
# only 1 large speed => k=1 => that's S1 not S3! So pure dominant-large with 12 small + 1
# large is S1 (already proved). For S3 we need >=2 large. So the dominant-large *sub-regime*
# of the claim (non-max<=13, Vmax>=53) can ONLY occur with exactly... let's check:
log("\nStrategy A: can 'non-max all<=13 & Vmax>=53' even be S3? (needs k>=2 large)")
log("  k=#{v>13}. If 12 non-max are all <=13 then only Vmax>13 => k=1 => case S1, NOT S3.")
log("  So the claimed dominant-large sub-regime is DISJOINT from S3. Verifying by search...")

# ---- Strategy B: genuine S3 sets (>=2 large), broad random. ----
random.seed(7); buildsB = 0
TLIMIT = 480
for _ in range(60000):
    if time.time() - t0 > TLIMIT: break
    psize = random.randint(2, 11)
    P = sorted(random.sample(range(1, 14), psize))
    V0 = random.randint(15, 250)
    spread = random.choice([14, 18, 21, 28, 35, 45, 60])
    Lsize = 13 - len(set(P))
    if Lsize < 2: continue
    offs = sorted(random.sample(range(0, spread + 1), min(Lsize, spread + 1)))
    L = [V0 + o for o in offs]
    S = sorted(set(P) | set(L))
    pool = [x for x in range(1, 14) if x not in S]
    while len(S) < 13 and pool: S = sorted(set(S) | {pool.pop()})
    if len(S) != 13: continue
    consider(S); buildsB += 1

log(f"\nStrategy B done: {buildsB} builds, {nS3} valid S3 sets, {time.time()-t0:.1f}s")

# ---- Strategy C: adversarial near-tight S3 (make the via-Vmax ratio small). ----
# Use small parts that are drop-k cores of {1..13} (low W) and tight clusters with small Vmax
# (just above 13*Vmin). Vmin=1 typically. Vmax just over 13.
random.seed(123); buildsC = 0
for _ in range(60000):
    if time.time() - t0 > TLIMIT: break
    # small part: a subset of {1..13} that still helps covering
    psize = random.randint(6, 11)
    P = sorted(random.sample(range(1, 14), psize))
    Vmin = min(P)
    # cluster just above 13*Vmin with small spread (tight)
    base = max(13 * Vmin + 1, 14)
    V0 = base + random.randint(0, 40)
    spread = random.choice([14, 15, 16, 18, 21, 25])
    Lsize = 13 - len(set(P))
    if Lsize < 2: continue
    offs = sorted(random.sample(range(0, spread + 1), min(Lsize, spread + 1)))
    L = [V0 + o for o in offs]
    S = sorted(set(P) | set(L))
    pool = [x for x in range(1, 14) if x not in S]
    while len(S) < 13 and pool: S = sorted(set(S) | {pool.pop()})
    if len(S) != 13: continue
    consider(S); buildsC += 1

log(f"Strategy C done: {buildsC} builds, {nS3} valid S3 total, {time.time()-t0:.1f}s")

log("\n--- HUNT RESULTS ---")
log("S3 covering sets examined:", nS3)
log("C(S) failures (best_ratio<=1):", len(fail_C))
for s, br, bv in fail_C[:30]:
    log("   ", s, "best_ratio", br, float(br), "via", bv)
log("dominant-large claim failures (non-max<=13, Vmax>=53, via-Vmax<=1):", len(domlarge_fail))
for s, r in domlarge_fail[:10]:
    log("   ", s, float(r))
log("# sets where via-Vmax failed (>0 means Vmax not always the witness):", len(suspicious))
for s in suspicious[:15]:
    br, bv = best_criterion(s)
    log("    via-Vmax failed:", s, "-> best via", bv, "ratio", float(br))
if worst_ratio is not None:
    log("worst best-criterion ratio (among via-Vmax failures):", worst_ratio, float(worst_ratio))
    log("  at", worst_ratio_set)

# ---- Exact-M check on the most suspicious sets + worst ratio set ----
log("\n--- EXACT M(S) on suspicious / worst sets (target: M<1/14 = counterexample) ---")
checkset = list({tuple(s) for s in suspicious})
if worst_ratio_set: checkset.append(tuple(worst_ratio_set))
fail_M = []; minM = None; minM_set = None
for tup in checkset[:60]:
    S = list(tup)
    m = Mval(S)
    if minM is None or m < minM: minM = m; minM_set = S
    if m < F(1, 14): fail_M.append((S, m))
log("checked", min(len(checkset), 60), "sets")
log("M<1/14 failures:", len(fail_M))
for s, m in fail_M[:10]:
    log("    COUNTEREXAMPLE?", s, "M=", m, float(m))
if minM is not None:
    log("min exact M among checked:", minM, float(minM), "at", minM_set)
log("\nDONE total time %.1fs" % (time.time() - t0))
