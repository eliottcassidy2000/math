#!/usr/bin/env python3
"""
lrc14_mixed_slope_j8plus_census_deathstar_S16.py  (death-star-2026-07-12-S16)

THE MIXED-SLOPE j>=8 STRATUM = PERTURBED TIGHT DILATES.  Reduction (this session):
a family in stratum (3) [some admissible scale L > 91B exists, but every one has
j >= 8 impure] has, at such a scale, ALL 13 elements within L/91 of multiples of L.
Its lift multiset K = (round(v_i/L)):
  - some lift 0            -> tiny-speed runner, far/scale-peel territory (flagged);
  - <= 12 distinct lifts   -> M(K) >= 1/13 (LRC(<=13)), atom => M(V) >= 1/13 - B/L, LOOSE;
  - 13 distinct, M(K) >= 2/27 -> atom => LOOSE once B/L < 2/27 - 1/14 = 1/378;
  - 13 distinct, M(K) < 2/27, K covering (primitivized) -> klein ILP covering-min
    14/183 > 2/27 for speeds <= 182 => LOOSE (k_max > 182 => lift far-element, peel);
  - 13 distinct, M(K) < 2/27, non-covering: K in the WALL WINDOW; the tight locus
    (verified to speed 60, THM-612 GAP-A) = {AP {1..13}, GW {1..11,13,24}}; window
    inhabitants above 1/14 (mediants, M=3/41 etc.) clear by the atom per-family.
SO the stratum's genuinely-tight core = L*AP + b and L*GW + b with j >= 8 nonzero
small offsets b_i, mixed slopes.  This script censuses THOSE.

THE OBJECT.  Profile W = {(k_i, b_i)}: reach2(W) = sup_{s,u} min_i ||k_i s + b_i u||
= sup_s min( pure(s), gammaF(s) ),  pure(s) = min over b_i=0 of ||k_i s||,
gammaF(s) = max_u min over b_i!=0 of ||k_i s + b_i u||.  Every evaluated rational s
gives a CERTIFIED lower bound (exact Fractions).  M(V) >= reach2(W) - B/(2L).

THE LOCAL LP (predictor).  At the AP tight point (s,u) = (1/14, 0) all margins are
exactly 1/14; to first order margin_i gains  sign_i*(i*sigma + b_i*mu),  sign_i = +1
for i <= 6 (position i/14 below 1/2, moving up helps), -1 for i >= 8; i = 7 is free
(position 1/2, margin 1/2).  Feasibility of a uniform gain delta > 0  <=>  the slope
sets A = {b_i/i : i <= 6} and B = {b_i/i : i >= 8} (pure runners contribute slope 0
to their group) are STRICTLY SEPARATED (max B < min A or max A < min B).
Prediction: separated patterns are loose with MACRO margin; doubly-sign-mixed
patterns (0 in the closed hull of both A and B) are the near-tight candidates.

KEY STRUCTURAL FACTS (found while designing; drive the experiment design):
  - j=13 (no pure) at B=1 is TRIVIAL: at s=1/2 all +-1 centers collapse onto {0,1/2},
    gammaF(1/2) = 1/4 -> reach2 >= 1/4. Small-q collapse handles B<=3 similarly.
  - With pure runners, the witness s = a/q DIES iff q divides some pure lift
    (pure margin 0). A CLEAN q in (j, 13] gives gammaF >= 1/q >= 1/13 (j centers
    cannot fill q slots: a >= 2-slot vacancy gap) AND pure margin >= 1/q... > 1/14.
  - Blocking every q in {9..13} requires pure >= {9,10,11,12,13} — exactly the 5
    available pure slots at j=8: THE UNIQUE MAXIMAL BLOCKER pure = {9..13}. It
    falls back to q=7 (pure margin 1/7, gammaF = 1/14 boundary when the 8 impure
    centers fill all 7 slots) — where the Part-6b s-MOTION strictness applies.

EXPERIMENTS
  A. j=13 sanity (AP+GW, B=1..3): confirm the small-q collapse caps everything.
  B. SUPPORT CENSUS (B=1): all pure sets of size 1..5, random + adversarial
     slot-filling signs, early-exit cap 1/12; sub-cap corners collected.
  C. Corner deep-dive: exhaustive signs on the sub-cap supports, full refinement;
     the record-min pattern + its witness s* (expect s near 1/7-motion).
  D. Same census for the GW lift.
  E. End-to-end witness transport at L = 360360: primitive/DC check + exact
     certified margin of t = (m + s*)/L on the integer family.

Pull-frequency note: the runner syncs with origin between experiment blocks (the
harness session does the actual git pull; this script just prints block boundaries).
"""

from fractions import Fraction
from itertools import combinations
from math import gcd
import random

F = Fraction

def distZ(x: F) -> F:
    f = x - (x.numerator // x.denominator)
    return min(f, 1 - f)

# ---------------------------------------------------------------- gammaF
def gammaF_exact(phases, speeds):
    """max_u min_i ||phi_i + b_i u||, exact; fast path when all |b_i| = 1."""
    n = len(speeds)
    if n == 0:
        return F(1, 2), F(0)
    if all(abs(b) == 1 for b in speeds):
        pts = sorted(set(((-speeds[i] * phases[i]) % 1) for i in range(n)))
        m = len(pts)
        if m == 1:
            return F(1, 2), (pts[0] + F(1, 2)) % 1
        best, argu = F(-1), None
        for t in range(m):
            g = (pts[(t + 1) % m] - pts[t]) % 1
            if g > best:
                best, argu = g, (pts[t] + g / 2) % 1
        return best / 2, argu
    cands = set()
    for i in range(n):
        b, ph = speeds[i], phases[i]
        lo, hi = min(ph, ph + b), max(ph, ph + b)
        m0 = (lo - F(1, 2)).__floor__()
        m1 = (hi - F(1, 2)).__ceil__()
        for m in range(m0, m1 + 1):
            u = F(m + F(1, 2) - ph, b)
            if 0 <= u < 1:
                cands.add(u)
    for i in range(n):
        for j in range(i + 1, n):
            for eps in (1, -1):
                d = speeds[i] - eps * speeds[j]
                if d == 0:
                    continue
                rhs0 = eps * phases[j] - phases[i]
                for m in range(-abs(d) - 2, abs(d) + 3):
                    u = F(rhs0 + m, d)
                    if 0 <= u < 1:
                        cands.add(u)
    best, argu = F(-1), None
    for u in cands:
        v = min(distZ(phases[i] + speeds[i] * u) for i in range(n))
        if v > best:
            best, argu = v, u
    return best, argu

def profile_value_at_s(s, lifts, bs):
    """min(pure(s), gammaF(s)) for profile {(lifts_i, bs_i)}; exact."""
    pure = [lifts[i] for i in range(len(lifts)) if bs[i] == 0]
    imp_ph = [lifts[i] * s for i in range(len(lifts)) if bs[i] != 0]
    imp_b = [bs[i] for i in range(len(lifts)) if bs[i] != 0]
    pv = min((distZ(k * s) for k in pure), default=F(1, 2))
    gv, gu = gammaF_exact(imp_ph, imp_b)
    return min(pv, gv), gu

# s-grid: Farey up to Q + mediant refinement around the incumbent
def build_s_grid(Q=40):
    seen = set()
    out = []
    for q in range(2, Q + 1):
        for a in range(1, q):
            if gcd(a, q) == 1:
                s = F(a, q)
                if s not in seen:
                    seen.add(s)
                    out.append(s)
    return out

S_GRID = build_s_grid(40)
S_GRID_SMALL = build_s_grid(22)

def reach2_lower(lifts, bs, refine_rounds=2, top=6, grid=None, cap=None):
    """certified lower bound for reach2: exact values on grid + mediant refinement.
    With cap: return early once the certified bound reaches cap (pattern safely loose)."""
    vals = []
    for s in (grid if grid is not None else S_GRID):
        v, _ = profile_value_at_s(s, lifts, bs)
        if cap is not None and v >= cap:
            return v, s
        vals.append((v, s))
    vals.sort(reverse=True)
    incumbent = vals[0]
    pool = [sv[1] for sv in vals[:top]]
    for _ in range(refine_rounds):
        newpool = []
        for s in pool:
            p, q = s.numerator, s.denominator
            for (a2, b2) in ((p * 2 - 1, q * 2), (p * 2 + 1, q * 2),
                             (p * 3 - 1, q * 3), (p * 3 + 1, q * 3)):
                if 0 < a2 < b2:
                    s2 = F(a2, b2)
                    v2, _ = profile_value_at_s(s2, lifts, bs)
                    if v2 > incumbent[0]:
                        incumbent = (v2, s2)
                    newpool.append(s2)
        pool = [s for _, s in sorted(((profile_value_at_s(s, lifts, bs)[0], s)
                                      for s in newpool), reverse=True)[:top]]
    return incumbent  # (value, s*)

# ---------------------------------------------------------------- predictor
def slope_separated(lifts, bs):
    """Local LP feasibility at the AP tight point: slope sets of low (i<=6) vs
    high (i>=8) groups strictly separated (pure runners contribute slope 0)."""
    A = [F(bs[i], lifts[i]) for i in range(len(lifts)) if lifts[i] <= 6]
    Bg = [F(bs[i], lifts[i]) for i in range(len(lifts)) if lifts[i] >= 8]
    if not A or not Bg:
        return True
    return max(Bg) < min(A) or max(A) < min(Bg)

# ---------------------------------------------------------------- experiments
AP = list(range(1, 14))
GW = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 24]
WALL = F(1, 14)

def first_clean_q(pure_vals, qmax=16):
    """smallest q >= 2 with q dividing no pure lift (a 'clean' witness base)."""
    out = []
    for q in range(2, qmax + 1):
        if all(p % q != 0 for p in pure_vals):
            out.append(q)
    return out

def greedy_fill_signs(lifts, impure_idx, q, a=1):
    """adversarial signs: try to OCCUPY as many u-slots mod q as possible
    (centers -b_i k_i a mod q), to kill the base-q witness."""
    used = set()
    bs = {}
    for i in impure_idx:
        r = (lifts[i] * a) % q
        c_plus, c_minus = (-r) % q, r % q
        if c_plus not in used:
            bs[i] = 1; used.add(c_plus)
        elif c_minus not in used:
            bs[i] = -1; used.add(c_minus)
        else:
            bs[i] = 1
    return [bs.get(i, 0) for i in range(len(lifts))]

def exp_A():
    print("=" * 76)
    print("EXP A — j=13 (no pure): the small-q collapse makes B<=3 trivially loose")
    print("=" * 76)
    rng = random.Random(7)
    for lifts, name in ((AP, "AP"), (GW, "GW")):
        for Bmax in (1, 2, 3):
            worst = (F(10), None)
            for _ in range(40):
                bs = [rng.choice([1, -1]) * rng.randrange(1, Bmax + 1) for _ in range(13)]
                v, s = reach2_lower(lifts, bs, refine_rounds=0, top=1, cap=F(1, 12))
                if v < worst[0]:
                    worst = (v, s)
            print(f"  {name} lift, B={Bmax}: min over 40 random patterns "
                  f"R >= {worst[0]} (~{float(worst[0]):.5f})  [cap 1/12]")

def exp_B(lifts, name):
    print("=" * 76)
    print(f"EXP B — {name} lift, B=1: SUPPORT census (all pure sets size 1..5),")
    print("        random + adversarial slot-filling signs, cap 1/12")
    print("=" * 76)
    rng = random.Random(16)
    n = len(lifts)
    subcap = []
    tested = 0
    CAP = F(1, 12)
    for psz in range(1, 6):
        j = n - psz
        rec = (F(10), None, None)
        for pure_idx in combinations(range(n), psz):
            pure_vals = [lifts[i] for i in pure_idx]
            impure_idx = [i for i in range(n) if i not in pure_idx]
            cq = first_clean_q(pure_vals)
            trials = []
            for q in cq[:3]:
                for a in range(1, min(q, 4)):
                    if gcd(a, q) == 1:
                        trials.append(greedy_fill_signs(lifts, impure_idx, q, a))
            for _ in range(8):
                bs = [0] * n
                for i in impure_idx:
                    bs[i] = rng.choice([1, -1])
                trials.append(bs)
            for bs in trials:
                tested += 1
                v, s_at = reach2_lower(lifts, bs, refine_rounds=0, top=1, cap=CAP)
                if v < rec[0]:
                    rec = (v, tuple(bs), s_at)
                if v < CAP:
                    subcap.append((v, tuple(bs), s_at, tuple(pure_vals)))
        print(f"  |pure|={psz} (j={j}): min R >= {rec[0]} (~{float(rec[0]):.5f}) "
              f" b={rec[1]}")
    subcap.sort(key=lambda r: r[0])
    print(f"  tested {tested} patterns; sub-cap (<1/12): {len(subcap)}")
    for v, bs, s_at, pv in subcap[:8]:
        print(f"    R>= {v} (~{float(v):.5f}) pure={pv} s={s_at}")
    return subcap

def exp_C(subcap, lifts, name):
    print("=" * 76)
    print(f"EXP C — {name}: deep-dive the sub-cap corners (full refine, no cap)")
    print("=" * 76)
    if not subcap:
        # the known corner: pure = {9..13} (AP indices 8..12)
        corners = [tuple(range(8, 13))] if lifts == AP else []
        print("  (no sub-cap patterns from EXP B; deep-diving the theory corner)")
    else:
        corners = []
        seen = set()
        for v, bs, s_at, pv in subcap:
            key = tuple(sorted(i for i in range(len(lifts)) if bs[i] == 0))
            if key not in seen:
                seen.add(key)
                corners.append(key)
            if len(corners) >= 3:
                break
    rec_all = (F(10), None, None)
    for pure_idx in corners:
        impure_idx = [i for i in range(len(lifts)) if i not in pure_idx]
        j = len(impure_idx)
        rec = (F(10), None, None)
        if j <= 9:
            masks = range(2 ** (j - 1))  # fix last sign + (negation symmetry)
        else:
            masks = None
        pats = []
        if masks is not None:
            for mask in masks:
                bs = [0] * len(lifts)
                for t2, i in enumerate(impure_idx[:-1]):
                    bs[i] = 1 if (mask >> t2) & 1 else -1
                bs[impure_idx[-1]] = 1
                pats.append(bs)
        else:
            rng = random.Random(99)
            for _ in range(256):
                bs = [0] * len(lifts)
                for i in impure_idx:
                    bs[i] = rng.choice([1, -1])
                pats.append(bs)
        for bs in pats:
            v, s_at = reach2_lower(lifts, bs, refine_rounds=2, top=4)
            if v < rec[0]:
                rec = (v, tuple(bs), s_at)
        pure_vals = tuple(lifts[i] for i in pure_idx)
        print(f"  pure={pure_vals} (j={j}): EXHAUSTIVE({len(pats)}) min R >= {rec[0]} "
              f"(~{float(rec[0]):.5f}) at s={rec[2]}  b={rec[1]}")
        print(f"    margin over wall: {rec[0] - WALL} (~{float(rec[0] - WALL):.5f})")
        if rec[0] < rec_all[0]:
            rec_all = rec
    return rec_all

def exp_E(rec, lifts):
    print("=" * 76)
    print("EXP E — end-to-end transported witness at L = 360360")
    print("=" * 76)
    if rec[1] is None:
        print("  (no record pattern)"); return
    v_prof, bs, s_at = rec
    L = 360360
    v = [L * lifts[i] + bs[i] for i in range(len(lifts))]
    g = 0
    for x in v:
        g = gcd(g, x)
    dc = all(any(x % d == 0 for x in v) for d in range(2, 15))
    pure = [lifts[i] for i in range(len(lifts)) if bs[i] == 0]
    imp_ph = [lifts[i] * s_at for i in range(len(lifts)) if bs[i] != 0]
    imp_b = [bs[i] for i in range(len(lifts)) if bs[i] != 0]
    gv, ustar = gammaF_exact(imp_ph, imp_b)
    m = round(L * ustar - s_at)
    t = F(m + s_at, L)
    marg = min(distZ(x * t) for x in v)
    print(f"  worst profile b={bs}  (R>= {v_prof} ~{float(v_prof):.5f})")
    print(f"  V = 360360*k + b: gcd={g} DC={dc}")
    print(f"  witness t=({m}+{s_at})/{L}: EXACT margin = {marg} (~{float(marg):.6f})"
          f"  beats 1/14? {marg > WALL}")

if __name__ == "__main__":
    print("[BLOCK A start]")
    exp_A()
    print("[BLOCK A end — sync point]")
    sub_ap = exp_B(AP, "AP")
    print("[BLOCK B end — sync point]")
    rec_ap = exp_C(sub_ap, AP, "AP")
    print("[BLOCK C end — sync point]")
    sub_gw = exp_B(GW, "GW")
    rec_gw = exp_C(sub_gw, GW, "GW")
    print("[BLOCK D end — sync point]")
    exp_E(rec_ap if rec_ap[0] <= rec_gw[0] else rec_gw,
          AP if rec_ap[0] <= rec_gw[0] else GW)
    print("[ALL BLOCKS done]")
