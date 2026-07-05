#!/usr/bin/env python3
"""
mac-mini-2026-07-05-S52 -- HYP-4103: the MULTI-LIFT LEG of hdich.

Context: residue pinning (opus-S75) => tight-from-above primitive 12-families
with no 13-multiple are lifts W = {r + 13 k_r : r = 1..12} of the AP {1..12}.
Rigidity (= TightLooseDichotomy's tight side) needs: k != 0 => M(W) > 1/13,
quantitatively >= beta for the assembly's loose margin.

Single lifts (one k_r > 0): S51 swept k <= 12; opus-S77 kernel-checked all 144
(floor 14/169 at {1..11,168}, MISTAKE-104 attribution fix).  kps-S1 swept ALL
lift vectors with k_r <= 2 (531k adversaries, AP-unique).  THIS SESSION:

PART A  The deep-well LADDER closed form: the r >= 7 sieve-survivors at k = r
        are the 14r-lifts {1..12}\\{r} u {14r}; conjecture M = 14/(13(r+1)),
        r = 7..12.  Exact M + witness structure + certificates.
        Also re-verify the MISTAKE-104 pair: {1..11,25} = 1/12, {1..11,168} = 14/169.

PART B  SINGLE-LIFT FLOOR CLOSURE to the *floor-level* structural cutoff.
        S51/S77's cutoff 144 = 12B is the RIGIDITY cutoff (window certifies
        margin > 1/13 for killer > 144).  The FLOOR claim ("no single lift
        below beta* = 14/169") needs the beta*-level window: base {1..12}\\{r}
        margin 1/12 (LRC(12) citation), B = 12, window half-width
        delta = (1/12 - beta*)/12 = 1/24336; killer w escapes at level beta*
        iff its closed bad arc (length 2*beta*/w) is shorter than the window
        (2*delta), i.e. w > beta*/delta = 14*144 = 2016 (one-tooth containment,
        IntervalEscape/no_cover_of_long).  So sweep k <= 155 (killer <= 2016);
        killers > 2016 are window-certified >= beta*.  MISTAKE-104 discipline:
        sweep to the STRUCTURAL cutoff of the question being asked.

PART C  DOUBLE-LIFT (l=2) RIGIDITY SWEEP, complete on its structural domain:
        W(r,j,s,k) = {1..12}\\{r,s} u {r+13j, s+13k}, w_a = r+13j >= w_b = s+13k.
        Rigidity-level closures (both via LRC citations + one-tooth/measure):
          (i)  one-killer window over the 11-base incl. w_b (margin 1/12,
               B = w_b): w_a > 12*w_b  => strictly loose.  [opus-S76 shape]
          (ii) two-killer FEE over the 10-base {1..12}\\{r,s} (margin 1/11 =
               LRC(11) citation, B = 12): window L = 2*delta_2,
               delta_2 = (1/11 - 1/13)/12 = 1/858; per-killer bad measure
               inside the window <= 2L/13 + 2/(13w); clear point exists iff
               (2/13)(1/w_a + 1/w_b) < L*(13-4)/13, i.e.
               1/w_a + 1/w_b < 9*delta_2 = 3/286;  both >= 191 suffices.
               [= klein-S134 lonely_of_window_multi fee, at the n=13 level]
        Residual finite domain: w_b <= 190 AND w_b <= w_a <= 12*w_b.
        Sweep it EXHAUSTIVELY: sieve filter -> witness scan (rational-point
        certificates, kps-S2 atom shape) -> exact M escalation.  ZERO sets may
        end <= 1/13.  Also record the l=2 FLOOR over the domain and the
        window-certified level w_a/(12(w_a+w_b)) for the tail.

PART D  The FEE/TOWER THRESHOLD TABLE at the rigidity level for l = 2..6
        (l lifted coords, all-big closure): all w_i >= T_l with
        T_l = 156(13-l)/(13-2l): 191, 223, 281, 416, 1092.  The 2l < 13 wall
        (l <= 6) is the same ceiling as klein-S134's "<= 6 tops".
        Probes: l=3 exhaustive at k <= 4; l=4 at k <= 3 (beyond kps's k <= 2).

Certificates: every witness emitted as (t = a/q, margin mu/q) -- decidable
integer facts in kps-S2's rational_point_margin atom shape.
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
BETA_STAR = F(14, 169)          # the (conjectured) lift floor
AP = list(range(1, 13))

def M_exact(S):
    """Exact M via the profile library, with shrinking rmax caps."""
    for cap in (11, 8, 6, 4, 3, 2):
        p = profile(sorted(S), F(1, cap))
        m = p.M()
        if m is not None:
            return m
    return None

def margin_at(S, t):
    return min(abs(v * t - round(v * t)) for v in S)

def sieve_missing(W):
    """Moduli m in 2..12 with no multiple in W (=> M >= 1/m >= 1/12 by sieve)."""
    return [m for m in range(2, 13) if not any(v % m == 0 for v in W)]

def is_primitive(W):
    return reduce(gcd, W) == 1

# ----------------------------------------------------------------------------
# witness scan machinery: rational points t = a/q, exact integer margin test
# ----------------------------------------------------------------------------
def dist_q(x, q):
    x %= q
    return min(x, q - x)

def scan_witness(W, betanum, betaden, q_list):
    """Find (q, a, margin) with margin = min_v dist(a v mod q)/q >= beta.
    Exact integer compare: dist*betaden >= betanum*q.  Returns first hit."""
    for q in q_list:
        # skip q if some element is 0 mod q (margin 0 for all a)
        if any(v % q == 0 for v in W):
            continue
        for a in range(1, q // 2 + 1):
            if gcd(a, q) != 1:
                continue
            m = min(dist_q(a * v, q) for v in W)
            if m * betaden >= betanum * q:
                return q, a, F(m, q)
    return None

# the theory-driven denominator library: 13s (lift-separating) first, then
# small non-13 moduli, then per-set pair sums added by callers
Q_LIB = [13 * s for s in range(2, 25)] + [q for q in range(8, 41) if q % 13 != 0]

def scan_or_exact(W, betanum, betaden, extra_q=()):
    qs = list(extra_q) + Q_LIB
    hit = scan_witness(W, betanum, betaden, qs)
    if hit:
        return hit, None
    return None, M_exact(W)

# ============================================================================
log("=" * 78)
log("PART A -- the deep-well ladder {1..12}\\{r} u {14r} (attribution-corrected)")
log("=" * 78)
log(f"{'r':>3} {'family tail':>12} {'M exact':>12} {'14/(13(r+1))':>14} {'match':>6}  witness (t, margin)")
ladder_ok = True
for r in range(1, 13):
    W = sorted([v for v in AP if v != r] + [14 * r])
    M = M_exact(W)
    pred = F(14, 13 * (r + 1))
    match = (M == pred)
    if r >= 7 and not match:
        ladder_ok = False
    # exact witness: scan the 13s library densely for t attaining >= M
    hit = scan_witness(W, M.numerator, M.denominator,
                       [13 * (r + 1), 13 * (r + 2)] + Q_LIB + [sum(sorted(W)[:2]), W[-1] + 1, W[-1] - 1])
    wtxt = f"t={hit[1]}/{hit[0]}, mu={hit[2]}" if hit else "(argmax not on scanned grid)"
    log(f"{r:>3} {'+{'+str(14*r)+'}':>12} {str(M):>12} {str(pred):>14} {'YES' if match else 'no':>6}  {wtxt}")
log(f"\nladder closed form M = 14/(13(r+1)) for r = 7..12: {'CONFIRMED' if ladder_ok else 'REFUTED'}")

log("\nMISTAKE-104 re-verification:")
for W, expect in ((sorted(list(range(1, 12)) + [25]), F(1, 12)),
                  (sorted(list(range(1, 12)) + [168]), F(14, 169))):
    M = M_exact(W)
    log(f"  M({W}) = {M} (expected {expect}) {'OK' if M == expect else '** MISMATCH **'}")

# ============================================================================
log("\n" + "=" * 78)
log("PART B -- single-lift FLOOR closure to the beta*-level cutoff (k <= 155)")
log("=" * 78)
delta = (F(1, 12) - BETA_STAR) / 12
cutoff = BETA_STAR / delta
log(f"beta* = {BETA_STAR}; window half-width delta = {delta}; killer cutoff beta*/delta = {cutoff}")
assert cutoff == 2016
worst = None
below = []          # any single lift with M < beta* (would correct the floor)
viol = 0
n_sieved = n_surv = n_scan = n_exact = 0
for r in range(1, 13):
    base = [v for v in AP if v != r]
    for k in range(1, 156):
        w = r + 13 * k
        if w > 2016:
            break
        W = sorted(base + [w])
        if not is_primitive(W):
            continue
        if sieve_missing(W):
            n_sieved += 1
            continue                     # M >= 1/12 > beta*, closed
        n_surv += 1
        hit, M = scan_or_exact(W, BETA_STAR.numerator, BETA_STAR.denominator,
                               extra_q=(13 * (r + 1), w + 1, w - 1, w + r))
        if hit:
            n_scan += 1
            continue                     # certificate: margin >= beta*
        n_exact += 1
        if M <= ONE13:
            viol += 1
            log(f"  !! RIGIDITY VIOLATION r={r} k={k}: M={M}")
        if M < BETA_STAR:
            below.append((r, k, w, M))
        if worst is None or M < worst[0]:
            worst = (M, r, k, w)
log(f"swept: sieve-closed {n_sieved}, survivors {n_surv} (witness-certified {n_scan}, exact {n_exact})")
log(f"rigidity violations: {viol}")
if below:
    log("SINGLE LIFTS BELOW beta* (floor correction!):")
    for r, k, w, M in sorted(below, key=lambda x: x[3]):
        log(f"  r={r} k={k} ({{1..12}}\\{{{r}}} u {{{w}}}): M = {M}")
else:
    log(f"NO single lift below beta* = 14/169 in the full structural domain (k <= 155);")
    log(f"killers > 2016 window-certified >= beta*.  SINGLE-LIFT FLOOR = 14/169 CLOSED.")
if worst:
    log(f"(worst exact-M case seen: M = {worst[0]} at r={worst[1]}, k={worst[2]}, killer {worst[3]})")
log(f"[t = {time.time()-T0:.1f}s]")

# ============================================================================
log("\n" + "=" * 78)
log("PART C -- double-lift (l=2) RIGIDITY sweep, complete on structural domain")
log("=" * 78)
delta2 = (F(1, 11) - ONE13) / 12
fee2 = 9 * delta2
log(f"two-killer fee: 1/w_a + 1/w_b < 9*delta_2 = {fee2} (= 3/286); both >= 191 closes")
log(f"domain: w_b <= 190, w_b <= w_a <= 12*w_b (outside: window/fee-certified loose)")

t0 = time.time()
BSCAN_N, BSCAN_D = 1, 13     # rigidity scan level: margin STRICTLY > 1/13 wanted;
                             # integer test uses >= so scan at 1/13 then verify strict
stats = dict(total=0, nonprim=0, sieved=0, scanned=0, exact=0, viol=0)
floor2 = None                # min exact/cert margin seen (over survivors)
floor2_exact = None          # min exact M computed
hard = []                    # sets needing exact M
cert_sample = []             # (family, q, a, margin) certificate examples
# precompute base margins per (r,s) for the scan library to speed the killer test
for s in range(1, 13):
    for k in range(1, 15):
        w_b = s + 13 * k
        if w_b > 190:
            break
        for r in range(1, 13):
            if r == s:
                continue
            base10 = [v for v in AP if v != r and v != s]
            jlo = max(1, (w_b - r + 12) // 13)
            jhi = (12 * w_b - r) // 13
            for j in range(jlo, jhi + 1):
                w_a = r + 13 * j
                if w_a < w_b:            # enforce ordering (strict dedupe below)
                    continue
                if w_a == w_b:
                    continue             # impossible (r != s mod 13) but cheap guard
                # dedupe the symmetric representation: keep only (big, small) with
                # (w_a, r) lexicographically the big side; pairs are unordered
                if (w_a == w_b and r > s):
                    continue
                stats['total'] += 1
                W = sorted(base10 + [w_a, w_b])
                if not is_primitive(W):
                    stats['nonprim'] += 1
                    continue
                if sieve_missing(W):
                    stats['sieved'] += 1
                    continue
                # witness scan at level > 1/13: use margin >= 1/13 + tie-break:
                # integer test dist*13 > q  <=> margin > 1/13 strictly
                found = None
                for q in ([13 * u for u in (2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)] +
                          [q for q in range(8, 41) if q % 13] +
                          [w_b + 1, w_b - 1, w_a + 1, w_a - 1, w_a + w_b, abs(w_a - w_b)]):
                    if q < 2 or any(v % q == 0 for v in W):
                        continue
                    for a in range(1, q // 2 + 1):
                        if gcd(a, q) != 1:
                            continue
                        m = min(dist_q(a * v, q) for v in W)
                        if m * 13 > q:               # STRICT > 1/13
                            found = (q, a, F(m, q))
                            break
                    if found:
                        break
                if found:
                    stats['scanned'] += 1
                    if floor2 is None or found[2] < floor2[0]:
                        floor2 = (found[2], tuple(W), found[0], found[1])
                    if len(cert_sample) < 5:
                        cert_sample.append((tuple(W), found))
                    continue
                stats['exact'] += 1
                M = M_exact(W)
                hard.append((tuple(W), M, r, j, s, k))
                if M <= ONE13:
                    stats['viol'] += 1
                    log(f"  !! RIGIDITY VIOLATION: W={W} M={M} (r={r},j={j},s={s},k={k})")
                if floor2_exact is None or M < floor2_exact[0]:
                    floor2_exact = (M, tuple(W), r, j, s, k)
        # progress
    log(f"  s={s} done  [total {stats['total']}, sieved {stats['sieved']}, "
        f"scanned {stats['scanned']}, exact {stats['exact']}]  t={time.time()-t0:.0f}s")

log(f"\nl=2 sweep complete: {stats}")
log(f"RIGIDITY: {'ZERO violations -- every double lift in the structural domain is strictly loose' if stats['viol']==0 else 'VIOLATIONS FOUND'}")
if floor2_exact:
    log(f"hard-case exact floor: M = {floor2_exact[0]} at {floor2_exact[1]} (r,j,s,k = {floor2_exact[2:]})")
if hard:
    hard.sort(key=lambda x: x[1])
    log(f"hard cases (needed exact M): {len(hard)}; 10 smallest:")
    for W, M, r, j, s, k in hard[:10]:
        log(f"   M = {str(M):>10}  W = {list(W)}   (r={r},j={j},s={s},k={k})")
log(f"scan-certificate floor (over witness-certified sets): "
    f"{floor2[0] if floor2 else None} at {floor2[1] if floor2 else None}")
log("sample certificates (kps atom shape: t=a/q, margin mu/q):")
for W, (q, a, mu) in cert_sample:
    log(f"   W={list(W)}: t={a}/{q}, mu={mu}")
log(f"[t = {time.time()-T0:.1f}s]")

# ============================================================================
log("\n" + "=" * 78)
log("PART D -- fee/tower thresholds (l = 2..6) + l=3 probe beyond kps k<=2")
log("=" * 78)
log(f"{'l':>3} {'base margin':>12} {'delta_l':>12} {'fee S_l':>14} {'all-big T_l':>12}")
for l in range(2, 7):
    d_l = (F(1, 13 - l) - ONE13) / 12
    S_l = d_l * (13 - 2 * l)
    T_l = int(l / S_l) + 1
    log(f"{l:>3} {'1/'+str(13-l):>12} {str(d_l):>12} {str(S_l):>14} {T_l:>12}")
log("(2l < 13 wall: l <= 6 -- same ceiling as klein-S134's <=6 tops.  l >= 7: fee dies;")
log(" those strata have <= 5 unlifted elements; kps-S1 swept all l at k <= 2; OPEN beyond.)")

log("\nl=3 probe: all C(12,3) coordinate triples, k <= 4 each (beyond kps k <= 2):")
t0 = time.time()
p3 = dict(total=0, sieved=0, scanned=0, exact=0, viol=0)
floor3 = None
hard3 = []
for trip in combinations(range(1, 13), 3):
    base9 = [v for v in AP if v not in trip]
    for k1 in range(1, 5):
        for k2 in range(1, 5):
            for k3 in range(1, 5):
                lifts = [trip[0] + 13 * k1, trip[1] + 13 * k2, trip[2] + 13 * k3]
                W = sorted(base9 + lifts)
                p3['total'] += 1
                if not is_primitive(W):
                    continue
                if sieve_missing(W):
                    p3['sieved'] += 1
                    continue
                found = None
                for q in ([13 * u for u in range(2, 15)] +
                          [q for q in range(8, 41) if q % 13]):
                    if any(v % q == 0 for v in W):
                        continue
                    for a in range(1, q // 2 + 1):
                        if gcd(a, q) != 1:
                            continue
                        m = min(dist_q(a * v, q) for v in W)
                        if m * 13 > q:
                            found = (q, a, F(m, q))
                            break
                    if found:
                        break
                if found:
                    p3['scanned'] += 1
                    if floor3 is None or found[2] < floor3[0]:
                        floor3 = (found[2], tuple(W))
                    continue
                p3['exact'] += 1
                M = M_exact(W)
                hard3.append((tuple(W), M))
                if M <= ONE13:
                    p3['viol'] += 1
                    log(f"  !! VIOLATION: {W} M={M}")
log(f"l=3 probe: {p3}  t={time.time()-t0:.0f}s")
if hard3:
    hard3.sort(key=lambda x: x[1])
    log("l=3 hard cases (exact M), smallest 8:")
    for W, M in hard3[:8]:
        log(f"   M = {str(M):>10}  W = {list(W)}")
log(f"l=3 scan floor: {floor3[0] if floor3 else None} at {floor3[1] if floor3 else None}")

log(f"\nDONE  [total t = {time.time()-T0:.1f}s]")
