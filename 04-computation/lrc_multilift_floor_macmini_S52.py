#!/usr/bin/env python3
"""
mac-mini-2026-07-05-S52 -- HYP-4103 part 2: the MULTI-LIFT FLOOR.

Part 1 (lrc_multilift_leg_macmini_S52.py) closed RIGIDITY for l=2 on its full
structural domain and confirmed the ladder M({1..12}\\{r} u {14r}) = 14/(13(r+1)).
Its scan floor hit 2/25 at the +13 BLOCK LIFT {1..12}\\{4,6} u {17,19} -- BELOW
the single-lift floor 14/169.  This script pins the floor:

E1  The HEIGHT-1 HYPERCUBE: all 4096 lifts W_C = ({1..12}\\C) u (C+13), C subset
    of {1..12} (k_r in {0,1}).  EXACT M for every C.  Floor + bottom table +
    species anatomy.  (Inside kps-S1's k<=2 rigidity space -- their sweep says
    AP-unique tight; here we get the exact FLOOR the assembly's beta needs.)

E2  l=2 FLOOR-LEVEL closure at beta = 2/25: structural domain grows to
    w_b <= 258 (fee at 2/25: 1/w_a + 1/w_b < 17/2200), w_a <= 24*w_b (window
    at 2/25 over the 11-base).  Scan at >= 2/25; EXACT M for every failure =
    the complete sub-2/25 stratum of the double-lift space.

E3  The ladder witness LAW: a_r = -14*(13-r)^{-1} mod 13(r+1), margin 14/(13(r+1)),
    for r = 7..12 -- verify + per-runner distance anatomy (equioscillation).
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
import sys, time

sys.path.insert(0, '04-computation')
from lonely_profile import profile

T0 = time.time()
def log(msg=""):
    print(msg, flush=True)

ONE13 = F(1, 13)
BETA_STAR = F(14, 169)
BETA_BLOCK = F(2, 25)
AP = list(range(1, 13))

def M_exact(S):
    for cap in (11, 8, 6, 4, 3, 2):
        p = profile(sorted(S), F(1, cap))
        m = p.M()
        if m is not None:
            return m
    return None

def dist_q(x, q):
    x %= q
    return min(x, q - x)

def sieve_missing(W):
    return [m for m in range(2, 13) if not any(v % m == 0 for v in W)]

def is_primitive(W):
    return reduce(gcd, W) == 1

# ============================================================================
log("=" * 78)
log("E1 -- the height-1 hypercube: W_C = ({1..12}\\C) u (C+13), all 4096 C")
log("     (scan-first at 2/25; EXACT on scan-failures = the definitive sub-2/25")
log("      stratum; exact also on |C|<=2, consecutive blocks, even blocks, sample)")
log("=" * 78)
t0 = time.time()
E1_QLIST = ([25, 50] + [13 * u for u in range(2, 21)] +
            [q for q in range(8, 41) if q % 13])
def scan25(W):
    for q in E1_QLIST:
        if any(v % q == 0 for v in W):
            continue
        for a in range(1, q // 2 + 1):
            if gcd(a, q) != 1:
                continue
            m = min(dist_q(a * v, q) for v in W)
            if m * 25 >= 2 * q:
                return (q, a, F(m, q))
    return None

fails = []            # scan-failures -> exact (everything possibly < 2/25)
cert25 = []           # sets whose best cert is exactly 2/25 (floor attainers?)
n_prim = 0
masks_exact = set()
for mask in range(1, 4096):
    C = [r for r in range(1, 13) if mask & (1 << (r - 1))]
    W = sorted([v for v in AP if v not in C] + [c + 13 for c in C])
    if not is_primitive(W):
        continue
    n_prim += 1
    hit = scan25(W)
    if hit is None:
        fails.append((tuple(C), tuple(W)))
        masks_exact.add(mask)
    elif hit[2] == BETA_BLOCK:
        cert25.append((tuple(C), hit))
log(f"primitive height-1 lifts: {n_prim} / 4095; scan(>=2/25) failures: {len(fails)}; "
    f"exact-2/25 certs: {len(cert25)}   [{time.time()-t0:.0f}s]")

log("\nEXACT M for every scan-failure (the complete candidate sub-2/25 stratum):")
sub25 = []
for C, W in fails:
    M = M_exact(W)
    tag = "  << BELOW 2/25" if M < BETA_BLOCK else ""
    if M <= ONE13:
        tag = "  !! RIGIDITY VIOLATION"
    if M < BETA_BLOCK:
        sub25.append((M, C, W))
    log(f"   M = {str(M):>10}  C = {list(C)}  W = {list(W)}{tag}")
sub25.sort()
log(f"\nheight-1 sets with M < 2/25: {len(sub25)}")
log(f"HEIGHT-1 FLOOR = {sub25[0][0] if sub25 else 'not below 2/25 (>= 2/25 everywhere off-AP)'}"
    f"{' at C = ' + str(list(sub25[0][1])) if sub25 else ''}")

log("\nexact M on structured slices (anatomy):")
def exact_row(C):
    W = sorted([v for v in AP if v not in C] + [c + 13 for c in C])
    if not is_primitive(W):
        return None
    return M_exact(W), tuple(C), tuple(W)
slices = []
for r in range(1, 13):
    slices.append((r,))
for r in range(1, 13):
    for s in range(r + 1, 13):
        slices.append((r, s))
for a in range(1, 13):
    for b in range(a + 2, 13):
        slices.append(tuple(range(a, b + 1)))       # consecutive blocks len >= 3
import random
random.seed(52)
for _ in range(150):
    mask = random.randrange(1, 4096)
    slices.append(tuple(r for r in range(1, 13) if mask & (1 << (r - 1))))
seen = set()
rows = []
for C in slices:
    if C in seen:
        continue
    seen.add(C)
    row = exact_row(list(C))
    if row:
        rows.append(row)
rows.sort()
log(f"{'M':>10} {'C (lifted coords)':>30}")
for M, C, W in rows[:25]:
    consec = (list(C) == list(range(C[0], C[0] + len(C))))
    log(f"{str(M):>10} {str(list(C)):>30}  {'CONSEC' if consec and len(C)>1 else ''}")
log(f"...top 3 loosest: {[(str(M), list(C)) for M, C, W in rows[-3:]]}")
below_bstar = [r for r in rows if r[0] < BETA_STAR]
log(f"structured-slice sets with M < 14/169: {len(below_bstar)}: "
    f"{[(str(M), list(C)) for M, C, W in below_bstar[:12]]}")
log(f"[t = {time.time()-T0:.0f}s]")

# ============================================================================
log("\n" + "=" * 78)
log("E2 -- l=2 floor-level closure at beta = 2/25 (complete sub-2/25 stratum)")
log("=" * 78)
d2 = (F(1, 11) - BETA_BLOCK) / 12
fee = d2 * (1 - 4 * BETA_BLOCK) / BETA_BLOCK
log(f"level 2/25: delta_2 = {d2}; fee: 1/w_a + 1/w_b < {fee}; both >= {int(2/fee)+1} closes")
log(f"one-killer window at 2/25 over 11-base: w_a > 24*w_b closes; domain w_b <= {int(2/fee)}, w_a <= 24*w_b")
t0 = time.time()
stats = dict(total=0, sieved=0, scanned=0, exact=0)
sub = []                       # the complete list of double lifts with M < 2/25
QLIST = ([13 * u for u in range(2, 21)] + [q for q in range(8, 41) if q % 13] + [25, 50])
wbmax = int(2 / fee)           # 258
for s in range(1, 13):
    for k in range(1, (wbmax - s) // 13 + 1):
        w_b = s + 13 * k
        if w_b > wbmax:
            break
        for r in range(1, 13):
            if r == s:
                continue
            base10 = [v for v in AP if v != r and v != s]
            jlo = max(1, (w_b - r + 12) // 13)
            jhi = (24 * w_b - r) // 13
            for j in range(jlo, jhi + 1):
                w_a = r + 13 * j
                if w_a < w_b:
                    continue
                stats['total'] += 1
                W = sorted(base10 + [w_a, w_b])
                if sieve_missing(W):
                    stats['sieved'] += 1
                    continue
                found = None
                for q in QLIST + [w_a + w_b, abs(w_a - w_b), w_b + 1, w_a + 1]:
                    if q < 2 or any(v % q == 0 for v in W):
                        continue
                    for a in range(1, q // 2 + 1):
                        if gcd(a, q) != 1:
                            continue
                        m = min(dist_q(a * v, q) for v in W)
                        if m * 25 >= 2 * q:           # margin >= 2/25
                            found = (q, a, F(m, q))
                            break
                    if found:
                        break
                if found:
                    stats['scanned'] += 1
                    continue
                stats['exact'] += 1
                M = M_exact(W)
                if M < BETA_BLOCK:
                    sub.append((M, tuple(W), r, j, s, k))
    log(f"  s={s}: {stats}  t={time.time()-t0:.0f}s")
log(f"\nl=2 floor-level sweep: {stats}")
sub.sort()
log(f"COMPLETE sub-2/25 stratum of double lifts (structural domain): {len(sub)} sets")
for M, W, r, j, s, k in sub[:25]:
    log(f"   M = {str(M):>10}  W = {list(W)}  (r={r},j={j},s={s},k={k})")
if sub:
    log(f"\nl=2 TRUE FLOOR = {sub[0][0]} at {list(sub[0][1])}")
else:
    log("\nl=2 floor >= 2/25 everywhere on the domain")
log(f"[t = {time.time()-T0:.0f}s]")

# ============================================================================
log("\n" + "=" * 78)
log("E3 -- the ladder witness law: a_r = -14*(13-r)^{-1} mod 13(r+1)")
log("=" * 78)
for r in range(7, 13):
    q = 13 * (r + 1)
    inv = pow(13 - r, -1, q)
    a = (-14 * inv) % q
    if a > q // 2:
        a = q - a              # mirror to the canonical half
    W = sorted([v for v in AP if v != r] + [14 * r])
    dists = sorted((dist_q(a * v, q), v) for v in W)
    mu = F(dists[0][0], q)
    binders = [v for d, v in dists if d == dists[0][0]]
    ok = (mu == F(14, q))
    log(f"r={r:>2}: q={q}, a={a}, margin={mu} (=14/q: {ok}); binders (dist=14): {binders}")
log(f"\nDONE  [total t = {time.time()-T0:.0f}s]")
