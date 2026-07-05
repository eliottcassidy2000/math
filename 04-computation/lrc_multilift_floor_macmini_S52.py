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
log("E1 -- the height-1 hypercube: W_C = ({1..12}\\C) u (C+13), all 4096 C, EXACT")
log("=" * 78)
t0 = time.time()
rows = []
n_prim = 0
for mask in range(1, 4096):            # C = emptyset is the AP itself (tight)
    C = [r for r in range(1, 13) if mask & (1 << (r - 1))]
    W = sorted([v for v in AP if v not in C] + [c + 13 for c in C])
    if not is_primitive(W):
        continue
    n_prim += 1
    M = M_exact(W)
    rows.append((M, tuple(C), tuple(W)))
rows.sort()
log(f"primitive height-1 lifts: {n_prim} / 4095 nonzero C   [{time.time()-t0:.0f}s]")
viol = [r for r in rows if r[0] <= ONE13]
log(f"rigidity violations (M <= 1/13): {len(viol)}")
log(f"\nHEIGHT-1 FLOOR = {rows[0][0]}   at C = {list(rows[0][1])}   W = {list(rows[0][2])}")
log("\nbottom 20 of the hypercube:")
log(f"{'M':>10} {'C (lifted coords)':>28}  W")
for M, C, W in rows[:20]:
    log(f"{str(M):>10} {str(list(C)):>28}  {list(W)}")
log("\ntop 5 (loosest):")
for M, C, W in rows[-5:]:
    log(f"{str(M):>10} {str(list(C)):>28}  {list(W)}")
# species anatomy: are the low rows the CONSECUTIVE/EVEN blocks?
log("\nanatomy of the bottom 50: |C|, min/max of C, consecutive?, all-even?")
from collections import Counter
sig = Counter()
for M, C, W in rows[:50]:
    consec = (list(C) == list(range(C[0], C[0] + len(C))))
    alleven = all(c % 2 == 0 for c in C)
    sig[(len(C), consec, alleven)] += 1
for (l, consec, alleven), cnt in sorted(sig.items()):
    log(f"  |C|={l:>2}  consecutive={str(consec):>5}  all-even={str(alleven):>5}  x{cnt}")
below_bstar = [r for r in rows if r[0] < BETA_STAR]
log(f"\nheight-1 lifts with M < 14/169: {len(below_bstar)}")
for M, C, W in below_bstar[:15]:
    log(f"   M = {str(M):>8}  C = {list(C)}")
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
