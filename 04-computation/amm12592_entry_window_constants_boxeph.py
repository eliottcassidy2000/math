"""LANE F1 -- extended-window constants for the ENTRY reduction.

Exact-rational (int/Fraction) certification, dyadic R = 2^9 .. 2^40,
eps in {1/32, 1/16}:

  (W1) a-priori F4 with the WIDE window: for every post-feed certifying row
       i0 <= i_pf + W, W = floor(c_w R), the Theorem S-cone-fc clocks fit
       the budget GIVEN F1^F2^F3 alone (a_t <= C(D0-1,t+1) from F3; exact
       staircase floors; K1c/K2c closed-form majorants; T_t <= T_2 for
       t >= 2 by the proved monotonicity lemma).  Largest c_w on the grid
       {1/8, 3/16, 1/4, 5/16, 3/8} passing at ALL scales is reported;
       the ENTRY window uses c_w = 1/8.
  (W2) drain deadline for ANY i0 >= i_pf: exact cell-0 clock invariance
       i0 + ceil(a0(i0)/2) = i_pf + ceil(a0(i_pf)/2)  <=  (R-2)/2 + (R-2)/2,
       certified via i_pf_ub <= (R-2)/2 (S3 route, per-R exact check).
  (W3) death-free window: junk front grows <= 2/row (T6a), death needs
       junk at cell d >= d_fe >= Dlb; threshold m_win = Dlb - 2W - 2 such
       that  m(i_pf) <= m_win  =>  no death in [i_pf, i_pf + W] whatever
       happens.  m_win/R reported (~0.49 at eps = 1/32, vs measured
       m(i_pf) ~ 0.93 sqrt R).

g sandwich: fresh integer Lucas/Fibonacci comparisons at M = 2^20
(627035/2^20 < g < 627036/2^20), re-derived here, no floats anywhere.
Output -> 05-knowledge/results/amm12592_entry_window_constants_boxeph.{out,json}
"""
import json, os, sys
from fractions import Fraction

RESULTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "05-knowledge", "results")


def fib_luc(n):
    if n == 0:
        return 0, 2
    F, L = fib_luc(n // 2)
    F2 = F * L
    L2 = L * L - (2 if n // 2 % 2 == 0 else -2)
    if n % 2 == 0:
        return F2, L2
    return (F2 + L2) // 2, (5 * F2 + L2) // 2


def g_gt(p, q, Fq, Lq):
    rhs = 2 * 5 ** p - Lq
    if rhs < 0:
        return True
    return 5 * Fq * Fq > rhs * rhs


def window_ok(cw_num, cw_den, g_lo, g_hi, eps_num, eps_den, kmin=9, kmax=40):
    """True iff (W1)+(W2)+(W3 threshold positive) verify for all dyadic
    R = 2^kmin..2^kmax at this eps and W = floor(cw*R).  Returns
    (ok, table) with per-R exact ints in the table."""
    table = []
    ok_all = True
    for kk in range(kmin, kmax + 1):
        R = 1 << kk
        D0 = -((-R * eps_num) // eps_den)
        W = (R * cw_num) // cw_den
        ipf_ub = int((R * (1 - g_lo) - D0) / (1 + g_lo)) + 2
        i_ub = ipf_ub + W
        Dlb = int((2 * g_lo * R + D0) / (1 + g_hi)) - 1
        Dub = int(g_hi * (R + i_ub)) + D0 + 1
        need1 = Dub * (Dub - 1) // 2
        need2 = Dub * (Dub - 1) * (Dub - 2) // 6

        def k1ok(K):
            return (g_lo + 2) * K * K - 4 * K >= need1

        def k2ok(K):
            return (g_lo * K * K / 2 - K) * (2 * Dlb - 3) >= need2

        def bs(okf):
            lo, hi = 1, 8 * R
            if not okf(hi):
                return None
            while lo < hi:
                mid = (lo + hi) // 2
                if okf(mid):
                    hi = mid
                else:
                    lo = mid + 1
            return lo

        K1c, K2c = bs(k1ok), bs(k2ok)
        m_win = Dlb - 2 * W - 2
        ok = (K1c is not None and K2c is not None and
              i_ub + max(K1c, K2c) <= R - 2 and       # (W1)
              ipf_ub <= (R - 2) // 2 and              # (W2)
              m_win > 0)                              # (W3)
        ok_all = ok_all and ok
        table.append({"R_log2": kk, "W": W, "i_ub": i_ub, "K1c": K1c,
                      "K2c": K2c, "budget_slack": (R - 2) - (i_ub +
                      max(K1c or 0, K2c or 0)), "m_win": m_win, "ok": ok})
    return ok_all, table


def main():
    sys.setrecursionlimit(100)
    M = 1 << 20
    F, L = fib_luc(2 * M)
    p_lo, p_hi = 627035, 627036
    assert g_gt(p_lo, M, F, L) and not g_gt(p_hi, M, F, L), "g sandwich FAIL"
    g_lo, g_hi = Fraction(p_lo, M), Fraction(p_hi, M)
    lines = [f"[PASS] g sandwich {p_lo}/2^20 < g < {p_hi}/2^20"]
    out = {"g_lo": [p_lo, M], "g_hi": [p_hi, M]}
    grid = [(1, 8), (3, 16), (1, 4), (5, 16), (3, 8)]
    for eps_num, eps_den in ((1, 32), (1, 16)):
        best = None
        for cn, cd in grid:
            ok, tab = window_ok(cn, cd, g_lo, g_hi, eps_num, eps_den)
            if ok:
                best = (cn, cd, tab)
        okW8, tabW8 = window_ok(1, 8, g_lo, g_hi, eps_num, eps_den)
        assert okW8, f"W = R/8 FAILED at eps={eps_num}/{eps_den}"
        worst = min(t["budget_slack"] for t in tabW8)
        mw_frac_min = min(Fraction(t["m_win"], 1 << t["R_log2"])
                          for t in tabW8)
        lines.append(
            f"[PASS] eps={eps_num}/{eps_den}: W = floor(R/8) a-priori-F4 "
            f"verified for ALL dyadic 2^9..2^40; min budget slack {worst}; "
            f"min m_win/R = {float(mw_frac_min):.4f}; largest grid c_w "
            f"passing = {best[0]}/{best[1]}")
        out[f"eps_{eps_num}_{eps_den}"] = {
            "W8_ok": okW8, "table_W8": tabW8,
            "largest_cw": [best[0], best[1]],
            "min_m_win_over_R": [mw_frac_min.numerator,
                                 mw_frac_min.denominator]}
        print(lines[-1], flush=True)
    json.dump(out, open(os.path.join(
        RESULTS, "amm12592_entry_window_constants_boxeph.json"), "w"))
    with open(os.path.join(
            RESULTS, "amm12592_entry_window_constants_boxeph.out"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("[done] JSON + out written.", flush=True)


if __name__ == "__main__":
    main()
