#!/usr/bin/env python3
"""
lrc14_uescape_j7_boundary_deathstar_S15.py  (death-star-2026-07-12-S15)

THE j=7 BOUNDARY OF THE U-ESCAPE (THM-721 Part 5 residual).

THM-721 Part 2 gives reach2(W) >= min(M_P, 1/(2j)) for a compressed profile
W = ((k_i,b_i)), j = #impure.  At j = 7 the union bound yields exactly 1/14 —
the wall, not strict.  This probe quantifies the two candidate STRICTNESS
mechanisms for j = 7:

  gamma(s) := max_u min_{i in F} ||k_i s + b_i u||       (the u-escape value at s)

  (M1) s-MOTION / wrap-gap defect: the impure forbidden u-arc systems can tile
       exactly at an isolated good s* (e.g. the near-dilate adversary at s=1/7:
       phases {7s..13s} = staggered 1/7-AP, gamma = 1/14 EXACTLY).  Moving s by
       eps inside the pure good set changes the circular gaps at rates =
       consecutive lift differences; the rates sum to 0, so unless ALL rates
       vanish the maxgap grows linearly: gamma = 1/14 + Omega(eps).
       All rates vanish  <=>  all impure lifts k_i equal.
  (M2) equal-lift corner: k_i = k common => common offset c = k s, and
       gamma = max_u min_i ||c + b_i u|| with 7 DISTINCT integer speeds b_i =
       a common-offset shifted LRC(8) slice.  Exact tiling by 7 arc systems of
       distinct |b_i| is DMNR-obstructed; +-b pairs (equal |b|) are reflections.
       Probe: inf_c gamma_b(c) over b-sets — is it > 1/14 always?

PROBES
  A. mechanism check on the near-dilate adversary profile (exact Fractions):
     gamma(1/7) = 1/14 exactly; gamma(1/7 + eps) = 1/14 + eps/2; the gap-rate law.
  B. adversarial profile search (j=7): minimize sup over good-s grid of gamma(s)
     over structured + random profiles (AP lifts, equal lifts, mixed, +-b).
     Record: does anything PIN sup_s gamma at 1/14?
  C. equal-lift corner exact scan: for 7-subsets of {+-1..+-B}, inf over a c-grid
     (with local refinement) of gamma_b(c); reflected pairs vs distinct moduli.
  D. stratum inhabitation: for record profiles, scan small L for primitive,
     divisor-complete, compressed-at-scale-L realizations (j=7 at that scale).

All arc/crossing arithmetic in exact Fractions.
"""

from fractions import Fraction
from itertools import combinations
import random, sys

F = Fraction

def distZ(x: Fraction) -> Fraction:
    f = x - (x.numerator // x.denominator)  # frac in [0,1)
    return min(f, 1 - f)

# ---------------------------------------------------------------- gamma(s)
def gamma_exact(phases, speeds):
    """max_u min_i ||phi_i + b_i u|| over u in [0,1), exact.
    Candidates: peaks (phi_i + b_i u = 1/2 mod 1) and pairwise crossings
    phi_i + b_i u = +-(phi_j + b_j u) mod 1."""
    n = len(speeds)
    if all(abs(b) == 1 for b in speeds):
        # bad centers -sign(b)*phi; gamma = (max circular gap)/2
        pts = sorted((-speeds[i] * phases[i]) % 1 for i in range(n))
        mg = max(((pts[(t + 1) % n] - pts[t]) % 1 for t in range(n)))
        # argu = midpoint of the max gap
        for t in range(n):
            if (pts[(t + 1) % n] - pts[t]) % 1 == mg:
                return mg / 2, (pts[t] + mg / 2) % 1
    cands = set()
    for i in range(n):
        b, ph = speeds[i], phases[i]
        lo = min(ph, ph + b)   # range of ph + b*u over u in [0,1]
        hi = max(ph, ph + b)
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
    if not cands:
        cands = {F(0)}
    best = F(-1)
    argu = None
    for u in cands:
        v = min(distZ(phases[i] + speeds[i] * u) for i in range(n))
        if v > best:
            best, argu = v, u
    return best, argu

def gamma_at_s(s, lifts, bs):
    phases = [k * s for k in lifts]
    return gamma_exact(phases, bs)

# ---------------------------------------------------------------- good set
def pure_margin(s, pure_lifts):
    return min(distZ(k * s) for k in pure_lifts)

def good_s_grid(pure_lifts, floor, qmax):
    """rationals a/q, q<=qmax, with pure margin >= floor (dedup)."""
    seen, out = set(), []
    for q in range(2, qmax + 1):
        for a in range(1, q):
            s = F(a, q)
            if s in seen:
                continue
            seen.add(s)
            if pure_margin(s, pure_lifts) >= floor:
                out.append(s)
    return out

# ---------------------------------------------------------------- probes
def probe_A():
    print("=" * 72)
    print("PROBE A — mechanism check: near-dilate adversary, impure lifts 7..13, b=1")
    print("=" * 72)
    pure = [1, 2, 3, 4, 5, 6]
    lifts = [7, 8, 9, 10, 11, 12, 13]
    bs = [1] * 7
    s0 = F(1, 7)
    g0, u0 = gamma_at_s(s0, lifts, bs)
    print(f"  s = 1/7 : pure margin = {pure_margin(s0, pure)} (= 1/7), "
          f"gamma = {g0} (= 1/14? {g0 == F(1,14)})  argu = {u0}")
    print("  moving s = 1/7 + eps (inside 1/13-good set needs eps <= 1/91):")
    for den in (1000, 300, 182, 91):
        eps = F(1, 7 * den)
        s = s0 + eps
        g, _ = gamma_at_s(s, lifts, bs)
        pm = pure_margin(s, pure)
        pred = F(1, 14) + eps / 2
        print(f"    eps = 1/{7*den:5d}: gamma = {g}  predicted 1/14+eps/2 = {pred} "
          f"match={g == pred}  pure_margin = {pm} (>=1/13 {pm >= F(1,13)})")
    # gap-rate law: gaps of {k s mod 1} for k in lifts at s=1/7+eps
    eps = F(1, 700)
    s = s0 + eps
    pts = sorted(F((k * s) - (k * s).numerator // (k * s).denominator) for k in lifts)
    pts = sorted((k * s) % 1 for k in lifts)
    gaps = [(pts[(t + 1) % 7] - pts[t]) % 1 for t in range(7)]
    print(f"  gaps at eps=1/700: {gaps}")
    print(f"  (six gaps 1/7+eps, one wrap gap 1/7-6eps — rates = lift differences)")

_GRID_CACHE = {}
def sup_gamma_over_good(pure, lifts, bs, floor=F(1, 13), qmax=48):
    key = (tuple(pure), floor, qmax)
    if key not in _GRID_CACHE:
        _GRID_CACHE[key] = good_s_grid(pure, floor, qmax)
    grid = _GRID_CACHE[key]
    if not grid:
        return None, None, 0
    best, args = F(-1), None
    for s in grid:
        g, _ = gamma_at_s(s, lifts, bs)
        if g > best:
            best, args = g, s
    return best, args, len(grid)

def probe_B():
    print("=" * 72)
    print("PROBE B — adversarial j=7 profiles: minimize sup_{good s} gamma(s)")
    print("=" * 72)
    pure = [1, 2, 3, 4, 5, 6]          # 6 distinct pure lifts, M_P >= 1/7
    rng = random.Random(721)
    profiles = []
    profiles.append(("AP lifts 7..13, b=1", [7, 8, 9, 10, 11, 12, 13], [1] * 7))
    profiles.append(("AP lifts 7..13, b=+-1 alt", [7, 8, 9, 10, 11, 12, 13],
                     [1, -1, 1, -1, 1, -1, 1]))
    profiles.append(("AP lifts 7..13, b=1..7", [7, 8, 9, 10, 11, 12, 13],
                     [1, 2, 3, 4, 5, 6, 7]))
    profiles.append(("equal lifts k=7, b=+-{1,2,3}+4", [7] * 7,
                     [1, -1, 2, -2, 3, -3, 4]))
    profiles.append(("equal lifts k=7, b=1..7", [7] * 7, [1, 2, 3, 4, 5, 6, 7]))
    profiles.append(("two blocks k=7,11, b mixed", [7, 7, 7, 11, 11, 11, 11],
                     [1, -1, 2, 1, -1, 2, -2]))
    profiles.append(("lifts=pure clones 1..6+7, b=1", [1, 2, 3, 4, 5, 6, 7], [1] * 7))
    profiles.append(("wide lifts 2,3,5,8,13,21,34 b=+-1", [2, 3, 5, 8, 13, 21, 34],
                     [1, -1, 1, -1, 1, -1, 1]))
    for _ in range(20):
        lifts = sorted(rng.randrange(0, 15) for _ in range(7))
        bs = [rng.choice([1, -1, 2, -2, 3, -3]) for _ in range(7)]
        profiles.append((f"rand k={lifts} b={bs}", lifts, bs))
    rec = (F(10), None)
    for name, lifts, bs in profiles:
        sup, args, ng = sup_gamma_over_good(pure, lifts, bs)
        tag = ""
        if sup is not None and sup < rec[0]:
            rec = (sup, name)
            tag = "  <-- record"
        if sup is not None and (sup <= F(1, 14) or tag):
            print(f"  sup_s gamma = {sup} (~{float(sup):.5f}) at s={args} "
                  f"[{ng} good s] : {name}{tag}")
    print(f"  RECORD min-over-profiles of sup_good-s gamma = {rec[0]} "
          f"(~{float(rec[0]):.5f})  1/14={float(F(1,14)):.5f}  profile: {rec[1]}")
    print(f"  strict? {rec[0] > F(1, 14)}")

def probe_C():
    print("=" * 72)
    print("PROBE C — equal-lift corner: inf_c gamma_b(c), 7 distinct b in +-{1..B}")
    print("=" * 72)
    def inf_gamma_c(bs, q1=97):
        # coarse grid + structured candidates (c = a/(14*lcm-ish)) + local refine
        best = (F(10), None)
        cand = [F(a, q1) for a in range(q1)]
        cand += [F(a, 56) for a in range(56)]
        for c in cand:
            g, _ = gamma_exact([c] * len(bs), bs)
            if g < best[0]:
                best = (g, c)
        c0 = best[1]
        for scale in (F(1, 10**3), F(1, 10**4)):
            improved = True
            while improved:
                improved = False
                for dc in (-scale, scale):
                    c = (c0 + dc) % 1
                    g, _ = gamma_exact([c] * len(bs), bs)
                    if g < best[0]:
                        best, c0, improved = (g, c), c, True
        return best
    universe4 = [1, -1, 2, -2, 3, -3, 4, -4]
    print("  B=4 (all 8 choose 7):")
    rec = (F(10), None)
    for bs in combinations(universe4, 7):
        g, c = inf_gamma_c(list(bs))
        flag = "  <-- record" if g < rec[0] else ""
        if flag:
            rec = (g, bs)
        print(f"    b={bs}: inf_c gamma ~ {g} (~{float(g):.6f}) at c={c}{flag}")
    print(f"  B=4 record: {rec[0]} (~{float(rec[0]):.6f}) at b={rec[1]}; "
          f"1/14 ~ {float(F(1,14)):.6f}; strictly above? {rec[0] > F(1,14)}")
    print("  B=5 sample (8 random 7-subsets of +-{1..5}):")
    rng = random.Random(183)
    universe5 = [1, -1, 2, -2, 3, -3, 4, -4, 5, -5]
    rec5 = (F(10), None)
    for _ in range(8):
        bs = tuple(sorted(rng.sample(universe5, 7)))
        g, c = inf_gamma_c(list(bs))
        print(f"    b={bs}: inf_c gamma ~ {g} (~{float(g):.6f})")
        if g < rec5[0]:
            rec5 = (g, bs)
    print(f"  B=5 sampled record: {rec5[0]} (~{float(rec5[0]):.6f}) at b={rec5[1]}; "
          f"strictly above 1/14? {rec5[0] > F(1,14)}")

def is_divisor_complete(v, nmax=14):
    return all(any(x % d == 0 for x in v) for d in range(2, nmax + 1))

from math import gcd
def probe_D():
    print("=" * 72)
    print("PROBE D — stratum inhabitation: DC + primitive + compressed j=7 at scale L")
    print("=" * 72)
    pure = [1, 2, 3, 4, 5, 6]
    lifts = [7, 8, 9, 10, 11, 12, 13]
    bs = [1] * 7
    found = 0
    for L in list(range(14, 200)) + [2 * 3 * 4 * 7, 420, 840]:
        v = sorted([k * L for k in pure] + [lifts[i] * L + bs[i] for i in range(7)])
        if len(set(v)) < 13:
            continue
        g = 0
        for x in v:
            g = gcd(g, x)
        if g != 1:
            continue
        if is_divisor_complete(v):
            found += 1
            if found <= 5:
                print(f"    L={L}: v={v}  primitive DC, j=7 at scale L (B=1)")
    print(f"  DC+primitive j=7-at-scale-L realizations found: {found}")
    print("  (inhabited => the j=7 leg is not vacuous; its gamma strictness is PROBE A/B)")

if __name__ == "__main__":
    probe_A()
    probe_B()
    probe_C()
    probe_D()
    print("=" * 72)
    print("SUMMARY hooks: M1 (s-motion) rates law; M2 (equal-lift) inf_c record;")
    print("adversarial record vs 1/14; DC inhabitation count.")
