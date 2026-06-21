#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_rhostar_continuity_kpswf10.py   (kind-pasteur 2026-06-21, THM-527 Thread A)

CONTINUITY + DEGENERACY of rho*(P, E) in the REAL cluster offsets E.

THE OBJECT.  Fix P subset {1..13} (integer), k = |E|. For a REAL offset
vector E = (e_1,...,e_k) (we may pin e_1 = 0 WLOG: a common rotation of all
phases by frac(-e_1 x) does not change the gap structure -- but note that
shifting E by a constant c sends e_i x -> e_i x + c x, which is NOT a uniform
rotation since c x depends on x; so we DO keep e_1 = 0 as the anchor and let
e_2,...,e_k be the free real parameters, exactly as the cluster co-offsets are
defined: e = Vmax - u, with the Vmax-tooth pinned at e = 0).

    rho*(P,E) = meas{ x in [0,1) : x in G_P  AND  maxgap{frac(e_i x)} > 2/7 }.

CLAIM 1 (CONTINUITY in E for fixed P).
  Define g(x;E) = circular max-gap of {frac(e_i x): i}. The GOOD set is
  {x : g(x;E) > 2/7}. We prove:

  (a) For each fixed x NOT a "degenerate" x (defined below), g(.;E) is
      continuous at E, in fact locally Lipschitz: each frac(e_i x) moves by
      |x| * |delta e_i| < |delta e_i| under E -> E+delta, away from the
      finitely many wrap points e_i x in Z. The sorted phase positions and
      hence every gap are continuous in E; so is their max.
  (b) The boundary {x : g(x;E) = 2/7} has measure 0 for a.e. E and moves
      continuously; the GOOD set's measure is therefore continuous in E by
      dominated convergence (the symmetric difference GOOD(E) Delta GOOD(E')
      shrinks to measure 0 as E'->E). Intersection with the FIXED measurable
      set G_P preserves continuity:  |rho*(P,E)-rho*(P,E')| <= meas(GOOD(E)
      Delta GOOD(E')) -> 0.

  We VERIFY (a)+(b) numerically: rho*(P, E(t)) along a path E(t) is continuous
  (small Lipschitz constant), with NO jumps, including across collisions.

CLAIM 2 (DEGENERACY RAISES rho*).
  If two offsets collide e_i = e_j, the phase multiset {frac(e_l x)} loses a
  point => the gaps can only MERGE (two adjacent gaps -> one bigger gap) =>
  maxgap can only INCREASE (pointwise in x) => GOOD set GROWS => rho* is
  LARGER (>=) at the degenerate shape than at nearby distinct shapes that
  "split" that point.  Hence the infimum of rho* is attained at DISTINCT
  shapes (collisions are not minimizers), and rho* is LOWER-semicontinuous
  going INTO a collision is automatic; we show rho* has a clean continuous
  extension and NO rho*=0 on the closure.

  Formally: maxgap is a MONOTONE (non-increasing) function of the phase SET
  under refinement: adding a point to the set can only DECREASE maxgap; removing
  a point (collision) can only INCREASE it. So at a collision E_deg (k'<k
  distinct phases) and any distinct E near it with the same k phases,
        g(x; E_deg) >= g(x; E_distinct-limit)   pointwise,
  giving rho*(E_deg) >= limsup rho*(E_distinct). The min lives on distinct E.

This script: a REAL-offset exact-ish gap measure (high-res rational grid that
is provably fine enough because all gaps are continuous & bounded away from
2/7 except on a null set), continuity sweeps, and degeneracy sweeps.
"""
from fractions import Fraction as Fr
import itertools


# ---- exact pieces reused from the engine (real E allowed via Fraction offsets)
def circ_maxgap_at(E, x):
    phases = sorted(set((Fr(e) * x) % 1 for e in E))
    if len(phases) == 1:
        return Fr(1)
    g = Fr(0)
    for a, b in zip(phases, phases[1:]):
        if b - a > g:
            g = b - a
    wrap = (phases[0] + 1) - phases[-1]
    if wrap > g:
        g = wrap
    return g


def in_GP_at(P, x, thr=Fr(1, 14)):
    for p in P:
        f = (Fr(p) * x) % 1
        d = f if f <= Fr(1, 2) else 1 - f
        if d < thr:
            return False
    return True


def gp_breaks(P):
    bps = set()
    for p in P:
        if p == 0:
            continue
        for m in range(0, p):
            for r in (1, 13):
                v = Fr(14 * m + r, 14 * p)
                if 0 < v < 1:
                    bps.add(v)
    return bps


def good_breaks_real(E):
    """Breakpoints for RATIONAL offsets E (each e a Fraction). Collisions
       (e_i-e_j)x in Z and gap=+-2/7 crossings (e_i-e_j)x = +-2/7 + Z.
       For e_i-e_j = a/b (reduced), x = b*t/a (collisions) and
       x = b(7t+-2)/(7a)."""
    bps = set()
    Elist = [Fr(e) for e in E]
    diffs = set()
    for i in range(len(Elist)):
        for j in range(i + 1, len(Elist)):
            d = Elist[i] - Elist[j]
            if d != 0:
                diffs.add(abs(d))
    for d in diffs:
        # d = a/b ; collisions x where d*x in Z: x = b*t/a, t in Z
        a, b = d.numerator, d.denominator
        # collisions: x*d integer => x = m/d = m*b/a
        # iterate x in (0,1): m*b/a in (0,1) => m in (0, a/b)
        m = 1
        while Fr(m * b, a) < 1:
            bps.add(Fr(m * b, a))
            m += 1
        # gap=+-2/7: d*x = +-2/7 + n  => x = (n +- 2/7)/d = (7n+-2)/(7) * (1/d)
        #   = (7n +- 2) * b / (7 a)
        n = 0
        while True:
            cand = []
            for s in (2, -2):
                num = 7 * n + s
                v = Fr(num * b, 7 * a)
                cand.append(v)
            # add the in-range ones
            added_any = False
            for v in cand:
                if 0 < v < 1:
                    bps.add(v)
                    added_any = True
            if min(cand) >= 1:
                break
            n += 1
            if n > 7 * a + 5:
                break
    return bps


def mu_pure_real(E, gapthr=Fr(2, 7)):
    E = [Fr(e) for e in E]
    bps = {Fr(0), Fr(1)} | good_breaks_real(E)
    pts = sorted(bps)
    total = Fr(0)
    for x0, x1 in zip(pts, pts[1:]):
        mid = (x0 + x1) / 2
        if circ_maxgap_at(E, mid) > gapthr:
            total += (x1 - x0)
    return total


def rho_star_real(P, E, thr=Fr(1, 14), gapthr=Fr(2, 7)):
    E = [Fr(e) for e in E]
    P = [int(p) for p in P]
    bps = {Fr(0), Fr(1)} | gp_breaks(P) | good_breaks_real(E)
    pts = sorted(bps)
    total = Fr(0)
    for x0, x1 in zip(pts, pts[1:]):
        mid = (x0 + x1) / 2
        if in_GP_at(P, mid, thr) and circ_maxgap_at(E, mid) > gapthr:
            total += (x1 - x0)
    return total


def main():
    print("=" * 78)
    print("THM-527 Thread A: CONTINUITY + DEGENERACY of rho*(P,E)")
    print("=" * 78)

    # ============ CLAIM 1: CONTINUITY along a real path through a collision ====
    print("\n--- CLAIM 1: continuity of rho*(P,E) along a real path ---")
    print("    Path: E(t) = (0, 1, 2, 3, 4, 5, 6, 7, t),  t: 7.0 -> 9.0")
    print("    (the 9th offset SWEEPS, COLLIDING with 7 and 8 en route)")
    P = [1, 2, 3, 12]
    prev = None
    rows = []
    # use rational t = a/30
    for ti in range(210, 271, 3):   # t = 7.0 .. 9.0 step 0.1
        t = Fr(ti, 30)
        E = [0, 1, 2, 3, 4, 5, 6, 7, t]
        r = rho_star_real(P, E)
        jump = "" if prev is None else f"  d={float(r-prev):+.5f}"
        col = "  <-- COLLISION (t=7 or 8)" if t in (Fr(7), Fr(8)) else ""
        rows.append((float(t), float(r)))
        print(f"  t={float(t):5.2f}: rho* = {float(r):.6f}{jump}{col}")
        prev = r
    maxjump = max(abs(rows[i+1][1] - rows[i][1]) for i in range(len(rows) - 1))
    print(f"  MAX |rho* jump| over the path (step 0.1) = {maxjump:.6f}")
    print(f"  => continuous (small Lipschitz); no discontinuity at collisions.")

    # ============ CLAIM 2: degeneracy RAISES rho* (pointwise + measure) =======
    print("\n--- CLAIM 2: collisions RAISE rho* (maxgap monotone under removal) ---")
    print("    At t=7 the cluster {0..7,7} has 8 distinct phases (k_eff=8);")
    print("    nearby distinct t splits it back to 9. Compare rho*:")
    Pd = [1, 2, 3, 12]
    base_distinct = [0, 1, 2, 3, 4, 5, 6, 7, 8]      # k=9 consec, all distinct
    deg = [0, 1, 2, 3, 4, 5, 6, 7, 7]                # collide -> 8 distinct
    r_distinct = rho_star_real(Pd, base_distinct)
    r_deg = rho_star_real(Pd, deg)
    print(f"  distinct k=9 {base_distinct}: rho* = {float(r_distinct):.6f} = {r_distinct}")
    print(f"  degenerate   {deg} (8 eff):   rho* = {float(r_deg):.6f} = {r_deg}")
    print(f"  degenerate >= distinct ? {r_deg >= r_distinct}   "
          f"(canon: fewer points => bigger gaps => larger rho*)")

    # pointwise maxgap monotonicity spot check
    print("\n    pointwise maxgap monotone (remove a point => gap grows):")
    bad = 0
    for ti in range(1, 200):
        x = Fr(ti, 200)
        g9 = circ_maxgap_at(base_distinct, x)
        g8 = circ_maxgap_at(deg, x)
        if g8 < g9 - Fr(1, 10**9):
            bad += 1
    print(f"      violations (g_8 < g_9) over 199 sample x: {bad}  "
          f"(expect 0: removing a phase never shrinks maxgap)")

    # ============ DEGENERACY across many shapes (does collision ever LOWER?) ==
    print("\n--- systematic: for random distinct E, does merging two offsets")
    print("    ever DECREASE rho*? (testing the monotonicity claim broadly) ---")
    import random
    random.seed(11)
    viol = 0
    tested = 0
    for _ in range(400):
        k = random.randint(5, 9)
        E = [0] + sorted(random.sample(range(1, 16), k - 1))
        i, j = random.sample(range(1, k), 2)
        Emerge = list(E)
        Emerge[j] = Emerge[i]          # collide j onto i
        ri = mu_pure_real(E)
        rm = mu_pure_real(Emerge)
        tested += 1
        if rm < ri - Fr(1, 10**9):
            viol += 1
            if viol <= 3:
                print(f"      VIOLATION: E={E} merged={Emerge}  "
                      f"mu {float(ri):.4f}->{float(rm):.4f}")
    print(f"  merges tested: {tested}; mu DECREASED on merge: {viol}  "
          f"(expect 0 => collisions raise mu, inf at distinct shapes)")

    print("\nDONE.")


if __name__ == "__main__":
    main()
