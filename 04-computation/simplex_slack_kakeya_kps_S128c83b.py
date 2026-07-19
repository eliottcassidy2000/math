#!/usr/bin/env python3
"""simplex_slack_kakeya_kps_S128c83b.py  -- kind-pasteur-2026-07-18-S128 (cont.83)

Builds on death-star-S58b's corrected maximiser bound (HYP-7735/7736) and on the
Kakeya/X-ray reading of THM-1154.

death-star's B (verified == bad_from_g) is, in the ordering region g_{s0}<=g_{s1}<=g_{s2}:

    g_{s0} <= 2/7 ,  g_{s1}-g_{s0} <= 2/7 ,  g_{s2}-g_{s1} <= 2/7 ,  g_{s2} >= 5/7 .

THIS SCRIPT TESTS FOUR THINGS.

(A) THE SLACK IDENTITY.  Put
        f1 = 2/7 - g_{s0},  f2 = 2/7 - (g_{s1}-g_{s0}),
        f3 = 2/7 - (g_{s2}-g_{s1}),  f4 = g_{s2} - 5/7 .
    Then f1+f2+f3+f4 = 6/7 - g_{s2} + g_{s2} - 5/7 = 1/7 IDENTICALLY.
    So each piece of B is a 3-SIMPLEX {f>=0, sum f = 1/7}, not a box, and its
    incentre (all slacks 1/28) is the permutation of (1/4,1/2,3/4).

(B) |B| IN CLOSED FORM.  The change of variables (g_{s0},g_{s1},g_{s2}) -> (f1,f2,f3)
    is unimodular, so each piece has volume (1/7)^3/6 and, over six orderings,
        |B| = (1/7)^3 = 1/343
    exactly -- PROVIDED each simplex sits inside its ordering region (checked here
    by verifying all four vertices satisfy the ordering).  Cross-checked by grid.

(C) THE GENERAL SINGLE-RUN BOUND, with NO centre-hit hypothesis.
    f4 decreases at rate c = d_{s2} and 0 <= f4 <= 1/7 on a run, so
        run <= (1/7)/c .
    death-star proved run = (1/28)(1/c + 1/max(a,b-a,c-b)) AT A CENTRE HIT; this is
    the same bound for EVERY run.  Checked exactly on every cell-run.

(D) THE RESIDUE, REDUCED TO ONE COMBINATORIAL INEQUALITY.
    Summing (C) over cell-runs and grouping by top coordinate i (N_i = #cell-runs
    whose max coordinate is i):
        sojourn <= (1/7) * SUM_i N_i/d_i .
    So the whole Kakeya sup bound sojourn <= 2/21 would follow from
        W(d) := SUM_i N_i/d_i  <=  2/3 ,
    which is TIGHT at (1,2,3) (2 runs, both topped by the rate-3 coordinate).
    This would close death-star's >=2-mirror-pair residue as well, since it never
    splits on the number of mirror pairs.  Tested exhaustively below.
"""
import sys
from fractions import Fraction as F
from math import gcd

N = int(sys.argv[1]) if len(sys.argv) > 1 else 20


# ---------- (A) slack identity + (B) simplex vertices ----------------------------
def slack_identity_check():
    """f1+f2+f3+f4 == 1/7 identically, on a lattice of ordered triples."""
    bad = 0
    S = 12
    for i in range(S + 1):
        for j in range(i, S + 1):
            for k in range(j, S + 1):
                g0, g1, g2 = F(i, S), F(j, S), F(k, S)
                f1 = F(2, 7) - g0
                f2 = F(2, 7) - (g1 - g0)
                f3 = F(2, 7) - (g2 - g1)
                f4 = g2 - F(5, 7)
                if f1 + f2 + f3 + f4 != F(1, 7):
                    bad += 1
    return bad


def simplex_vertices():
    """The four vertices of {f>=0, sum f = 1/7}, back in g-coordinates."""
    verts = []
    for which in range(4):
        f = [F(0)] * 4
        f[which] = F(1, 7)
        g0 = F(2, 7) - f[0]
        g1 = g0 + F(2, 7) - f[1]
        g2 = g1 + F(2, 7) - f[2]
        verts.append((g0, g1, g2, g2 - F(5, 7) == f[3]))
    return verts


def corner_simplex_derivation():
    """|B| = (1/7)^3 = 1/343, exactly.  Gap coordinates turn each piece into a cube corner.

    In the ordering region g_{s0}<=g_{s1}<=g_{s2} put the GAP coordinates
        x = g_{s0},   y = g_{s1}-g_{s0},   z = g_{s2}-g_{s1} .
    This map is unimodular (triangular, unit diagonal), so it preserves volume.
    The four constraints read
        x <= 2/7,  y <= 2/7,  z <= 2/7,  x+y+z >= 5/7,     x,y,z >= 0.
    Now set u = 2/7-x, v = 2/7-y, w = 2/7-z >= 0.  Then
        x+y+z >= 5/7  <=>  6/7 - (u+v+w) >= 5/7  <=>  u+v+w <= 1/7,
    and u,v,w <= 2/7 is then automatic since 1/7 < 2/7.  So each piece is the
    SIMPLEX {u,v,w >= 0, u+v+w <= 1/7}, of volume (1/7)^3/6.

    Two boundary checks, both strict, so there is no clipping or wraparound:
        x = 2/7 - u >= 2/7 - 1/7 = 1/7 > 0,
        g_{s2} = x+y+z <= 6/7 < 1.
    Six ordering regions => |B| = 6 * (1/7)^3/6 = (1/7)^3 = 1/343.
    """
    leg = F(1, 7)
    per_piece = leg ** 3 / 6
    return per_piece, 6 * per_piece


def weyl_measure_B(M):
    """Low-discrepancy (Kronecker) estimate of |B| -- converges far better than a
    midpoint grid, whose error on a polytope is O(surface * h)."""
    a1, a2, a3 = 0.8191725133961644, 0.6710436067037893, 0.5497004779019702
    x = y = z = 0.0
    cnt = 0
    for _ in range(M):
        x += a1
        y += a2
        z += a3
        if x >= 1.0:
            x -= 1.0
        if y >= 1.0:
            y -= 1.0
        if z >= 1.0:
            z -= 1.0
        lo, mid, hi = sorted((x, y, z))
        if lo <= 2 / 7 and mid - lo <= 2 / 7 and hi - mid <= 2 / 7 and hi >= 5 / 7:
            cnt += 1
    return cnt / M


# ---------- cell decomposition (death-star's constraints, + top coordinate) -------
def cell_runs(DD):
    """Maximal bad sub-intervals inside each ordering/lattice cell.
    Returns (lo, hi, top_index).  Ordering is constant on a cell, so `top` is well
    defined and the f4-rate on that run is exactly DD[top]."""
    bps = {F(0), F(1)}
    for d in DD:
        for k in range(d + 1):
            bps.add(F(k, d))
    for i in range(3):
        for j in range(i + 1, 3):
            dd = abs(DD[i] - DD[j])
            if dd > 0:
                for k in range(dd + 1):
                    bps.add(F(k, dd))
    bps = sorted(bps)
    out = []
    for idx in range(len(bps) - 1):
        lo, hi = bps[idx], bps[idx + 1]
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        gmid = [(-DD[i] * mid) % 1 for i in range(3)]
        a = [gmid[i] + DD[i] * mid for i in range(3)]
        order = sorted(range(3), key=lambda i: gmid[i])
        s0, s1, s2 = order
        cons = [("ge", a[s0] - F(2, 7), DD[s0]), ("le", a[s2] - F(5, 7), DD[s2])]
        feas = True
        for (x, y) in [(s1, s0), (s2, s1)]:
            c = a[x] - a[y]
            dc = DD[x] - DD[y]
            if dc > 0:
                cons.append(("ge", c - F(2, 7), dc))
            elif dc < 0:
                cons.append(("le", c - F(2, 7), dc))
            elif not (c <= F(2, 7)):
                feas = False
        if not feas:
            continue
        ulo, uhi = lo, hi
        for typ, const, dd in cons:
            b = const / dd
            if typ == "ge":
                ulo = max(ulo, b)
            else:
                uhi = min(uhi, b)
        if uhi > ulo:
            out.append((ulo, uhi, s2))
    return out


def merge(ivs):
    ivs = sorted((lo, hi) for lo, hi, _ in ivs)
    m = []
    for lo, hi in ivs:
        if m and lo <= m[-1][1]:
            m[-1] = (m[-1][0], max(m[-1][1], hi))
        else:
            m.append((lo, hi))
    return m


def merge_with_exit_top(ivs):
    """Merge cell-runs, and tag each merged run with the top coordinate of its LAST cell.

    Along a merged run the identity of the max coordinate can only switch to a SLOWER
    coordinate (a switch means g_i = g_j with both >= 5/7, and after it the larger one
    is the slower one), so the exit rate is the smallest top rate on the run.  That is
    the rate that matters for the f4 = 0 exit, so it is the right one to charge."""
    ivs = sorted(ivs)
    m = []
    for lo, hi, top in ivs:
        if m and lo <= m[-1][1]:
            if hi >= m[-1][1]:
                m[-1] = (m[-1][0], hi, top)
            else:
                m[-1] = (m[-1][0], m[-1][1], m[-1][2])
        else:
            m.append((lo, hi, top))
    return m


def is123orbit(d):
    s = sorted(d)
    return s[1] == 2 * s[0] and s[2] == 3 * s[0]


# ---------- main -----------------------------------------------------------------
print("=" * 78)
print("(A) SLACK IDENTITY  f1+f2+f3+f4 = 1/7")
print("=" * 78)
print("  violations on the ordered 12^3 lattice : %d   -> identity %s"
      % (slack_identity_check(), "HOLDS" if slack_identity_check() == 0 else "FAILS"))

print()
print("=" * 78)
print("(B) B IS SIX SIMPLICES ;  |B| = (1/7)^3 = 1/343")
print("=" * 78)
vs = simplex_vertices()
allin = True
for (g0, g1, g2, ok4) in vs:
    inorder = (0 <= g0 <= g1 <= g2 <= 1)
    allin = allin and inorder and ok4
    print("  vertex g = (%s, %s, %s)   ordered=%s  f4 consistent=%s"
          % (g0, g1, g2, inorder, ok4))
print("  all four vertices inside their ordering region : %s" % allin)
print(corner_simplex_derivation.__doc__)
per_piece, total = corner_simplex_derivation()
print("  => each piece = 3-simplex of volume (1/7)^3/6 = %s" % per_piece)
print("  => |B| = 6 * that = %s = %.7f" % (total, float(total)))
for M in (200000, 1000000, 4000000):
    est = weyl_measure_B(M)
    print("     Kronecker %8d pts : %.7f   (exact %.7f, rel err %+.2f%%)"
          % (M, est, float(total), 100 * (est / float(total) - 1)))
print("  incentre: all four slacks equal 1/28 at g = (1/4,1/2,3/4):",
      [F(2, 7) - F(1, 4), F(2, 7) - F(1, 4), F(2, 7) - F(1, 4), F(3, 4) - F(5, 7)])

print()
print("=" * 78)
print("(C)+(D) run <= (1/7)/d_top  ;  and  W(d) = SUM_i N_i/d_i  <=  2/3 ?")
print("=" * 78)
viol_run = 0
maxW = F(0)
argW = None
maxsoj = F(0)
argsoj = None
worst_by_pairs = {}
Wties = []
tested = 0
stats = {"viol_merged": 0, "viol_E": 0, "maxWp": F(0), "argWp": None,
         "Wpties": [], "maxWp_multi": F(0), "argWp_multi": None,
         "mirror_missing": 0, "viol_F": 0, "viol_Fsum": 0, "maxWpp": F(0),
         "argWpp": None, "Wppties": [], "maxWpp_multi": F(0), "argWpp_multi": None,
         "rho_by_P": {}, "viol_G": 0, "ex_G": None}
for d2 in range(1, N + 1):
    for d3 in range(1, N + 1):
        for d4 in range(1, N + 1):
            DD = [d2, d3, d4]
            if gcd(gcd(d2, d3), d4) != 1:
                continue
            tested += 1
            runs = cell_runs(DD)
            if not runs:
                continue
            # (C) per-cell-run bound
            for lo, hi, top in runs:
                if hi - lo > F(1, 7) / DD[top]:
                    viol_run += 1
            # (D) the weight W
            Nc = [0, 0, 0]
            for _, _, top in runs:
                Nc[top] += 1
            W = sum(F(Nc[i], DD[i]) for i in range(3))
            soj = sum(hi - lo for lo, hi in merge(runs))
            if W > maxW:
                maxW = W
                argW = tuple(DD)
                Wties = [tuple(DD)]
            elif W == maxW:
                Wties.append(tuple(DD))
            if soj > maxsoj:
                maxsoj = soj
                argsoj = tuple(DD)
            npairs = len(merge(runs))
            key = "1-pair (<=2 runs)" if npairs <= 2 else ">=2-pairs (>=4 runs)"
            if key not in worst_by_pairs or soj > worst_by_pairs[key][0]:
                worst_by_pairs[key] = (soj, tuple(DD), W)
            # ---- (E) charge each MERGED run to the top rate of its LAST cell ----
            mr = merge_with_exit_top(runs)
            Mc = [0, 0, 0]
            for lo, hi, top in mr:
                Mc[top] += 1
                if hi - lo > F(1, 7) / DD[top]:
                    stats["viol_merged"] += 1
            Wp = sum(F(Mc[i], DD[i]) for i in range(3))
            if soj > F(1, 7) * Wp:
                stats["viol_E"] += 1
            if Wp > stats["maxWp"]:
                stats["maxWp"] = Wp
                stats["argWp"] = tuple(DD)
                stats["Wpties"] = [tuple(DD)]
            elif Wp == stats["maxWp"]:
                stats["Wpties"].append(tuple(DD))
            if npairs > 2 and Wp > stats["maxWp_multi"]:
                stats["maxWp_multi"] = Wp
                stats["argWp_multi"] = tuple(DD)
            # ---- (F) mirror symmetry: charge each run to the FASTER of the two ----
            # mirror of the run [lo,hi] is [1-hi,1-lo], a run of EQUAL length whose
            # top coordinate is this run's BOTTOM coordinate.  So
            #     run <= (1/7)/max(d_exit(r), d_exit(r*)).
            byleft = {}
            for lo, hi, top in mr:
                byleft[lo] = (hi, top)
            Wpp = F(0)
            for lo, hi, top in mr:
                mlo = 1 - hi
                if mlo in byleft and byleft[mlo][0] == 1 - lo:
                    rate = max(DD[top], DD[byleft[mlo][1]])
                else:
                    stats["mirror_missing"] += 1
                    rate = DD[top]
                if hi - lo > F(1, 7) / rate:
                    stats["viol_F"] += 1
                Wpp += F(1, rate)
            if soj > F(1, 7) * Wpp:
                stats["viol_Fsum"] += 1
            if Wpp > stats["maxWpp"]:
                stats["maxWpp"] = Wpp
                stats["argWpp"] = tuple(DD)
                stats["Wppties"] = [tuple(DD)]
            elif Wpp == stats["maxWpp"]:
                stats["Wppties"].append(tuple(DD))
            if npairs > 2 and Wpp > stats["maxWpp_multi"]:
                stats["maxWpp_multi"] = Wpp
                stats["argWpp_multi"] = tuple(DD)
            # ---- (G) is rho >= 3P, P = #mirror pairs ?  (would prove W'' <= 2/3) ----
            P = len(mr) // 2 if len(mr) % 2 == 0 else (len(mr) + 1) // 2
            rhos = []
            for lo, hi, top in mr:
                mlo = 1 - hi
                if mlo in byleft and byleft[mlo][0] == 1 - lo:
                    rhos.append(max(DD[top], DD[byleft[mlo][1]]))
                else:
                    rhos.append(DD[top])
            if rhos:
                mr_ = min(rhos)
                rec = stats["rho_by_P"].setdefault(P, [10 ** 9, None])
                if mr_ < rec[0]:
                    rec[0] = mr_
                    rec[1] = tuple(DD)
                if mr_ < 3 * P:
                    stats["viol_G"] += 1
                    if stats["ex_G"] is None:
                        stats["ex_G"] = (tuple(DD), P, mr_)

print("  primitive directions tested (1..%d)^3 : %d" % (N, tested))
print("  (C) cell-runs violating run <= (1/7)/d_top : %d  -> %s"
      % (viol_run, "BOUND HOLDS" if viol_run == 0 else "BOUND BROKEN"))
print("  (D) MAX W(d) = %s = %.6f at %s   (2/3 = %.6f)  -> W<=2/3 %s"
      % (maxW, float(maxW), argW, 2 / 3, "HOLDS" if maxW <= F(2, 3) else "BROKEN"))
print("      directions attaining max W : %d ; all (1,2,3)-orbit : %s"
      % (len(Wties), all(is123orbit(t) for t in Wties)))
if not all(is123orbit(t) for t in Wties):
    print("      non-orbit maximisers (first 10):",
          [t for t in Wties if not is123orbit(t)][:10])
print("  MAX sojourn = %s = %.7f at %s  (2/21 = %.7f)"
      % (maxsoj, float(maxsoj), argsoj, 2 / 21))
for k in sorted(worst_by_pairs):
    soj, dd, W = worst_by_pairs[k]
    print("    %-22s max sojourn %s = %.6f at %s ; W = %s = %.4f ; (1/7)W = %.6f"
          % (k, soj, float(soj), dd, W, float(W), float(W) / 7))
print()
print("=" * 78)
print("(E) THE FIX: charge each MERGED run to the top rate of its LAST cell")
print("=" * 78)
print("  merged runs violating run <= (1/7)/d_exit : %d  -> %s"
      % (stats["viol_merged"], "HOLDS" if stats["viol_merged"] == 0 else "BROKEN"))
print("  directions violating sojourn <= (1/7)*W'  : %d  -> %s"
      % (stats["viol_E"], "HOLDS" if stats["viol_E"] == 0 else "BROKEN"))
print("  MAX W'(d) = %s = %.6f at %s   (2/3 = %.6f)  -> W'<=2/3 %s"
      % (stats["maxWp"], float(stats["maxWp"]), stats["argWp"], 2 / 3,
         "HOLDS" if stats["maxWp"] <= F(2, 3) else "BROKEN"))
print("      attaining : %d ; all (1,2,3)-orbit : %s"
      % (len(stats["Wpties"]), all(is123orbit(t) for t in stats["Wpties"])))
if not all(is123orbit(t) for t in stats["Wpties"]):
    print("      non-orbit maximisers (first 10):",
          [t for t in stats["Wpties"] if not is123orbit(t)][:10])
print("  MAX W' over >=2-mirror-pair directions = %s = %.6f at %s"
      % (stats["maxWp_multi"], float(stats["maxWp_multi"]), stats["argWp_multi"]))
print()
print("=" * 78)
print("(F) MIRROR-CHARGED:  run <= (1/7)/max(d_exit(r), d_exit(r*))")
print("=" * 78)
print("  merged runs whose mirror is not itself a merged run : %d  -> mirror symmetry %s"
      % (stats["mirror_missing"], "OK" if stats["mirror_missing"] == 0 else "FAILS"))
print("  runs violating run <= (1/7)/max(d_exit,d_exit*) : %d  -> %s"
      % (stats["viol_F"], "HOLDS" if stats["viol_F"] == 0 else "BROKEN"))
print("  directions violating sojourn <= (1/7)*W''       : %d  -> %s"
      % (stats["viol_Fsum"], "HOLDS" if stats["viol_Fsum"] == 0 else "BROKEN"))
print("  MAX W''(d) = %s = %.6f at %s   (2/3 = %.6f)  -> W''<=2/3 %s"
      % (stats["maxWpp"], float(stats["maxWpp"]), stats["argWpp"], 2 / 3,
         "HOLDS" if stats["maxWpp"] <= F(2, 3) else "BROKEN"))
print("      attaining : %d ; all (1,2,3)-orbit : %s"
      % (len(stats["Wppties"]), all(is123orbit(t) for t in stats["Wppties"])))
if not all(is123orbit(t) for t in stats["Wppties"]):
    print("      non-orbit maximisers (first 10):",
          [t for t in stats["Wppties"] if not is123orbit(t)][:10])
print("  MAX W'' over >=2-mirror-pair directions = %s = %.6f at %s"
      % (stats["maxWpp_multi"], float(stats["maxWpp_multi"]), stats["argWpp_multi"]))
print()
print("=" * 78)
print("(G) WHY W'' <= 2/3 ?  test  rho >= 3P  (P = # mirror pairs)")
print("=" * 78)
print("  P  |  min rho observed  |  witness           |  3P  |  rho>=3P")
for P in sorted(stats["rho_by_P"]):
    lo, wit = stats["rho_by_P"][P]
    print("  %-3d|  %-17d|  %-18s|  %-4d|  %s"
          % (P, lo, wit, 3 * P, "yes" if lo >= 3 * P else "NO"))
print("  directions violating rho >= 3P : %d %s"
      % (stats["viol_G"], "" if stats["ex_G"] is None else
         "(first: d=%s P=%d rho=%d)" % stats["ex_G"]))
print("  NOTE rho >= 3P  =>  sum over pairs of 1/rho <= P/(3P) = 1/3  =>  W'' <= 2/3.")
print()
print("  IMPLICATION: sojourn <= (1/7)*W'' always (C+F).  If W'' <= 2/3 is a theorem,")
print("  then sojourn <= 2/21 for EVERY direction -- no split on mirror pairs, so")
print("  death-star's >=2-pair residue closes too.")
