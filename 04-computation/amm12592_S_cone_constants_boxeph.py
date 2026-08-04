"""LANE E1 -- exact rational constants for the S invariant-cone theorem.

All exact (int / Fraction).  Independent re-certification of the g sandwich
(fresh Lucas/Fibonacci fast doubling; integer comparisons 5^a <=> phi^(2M)),
then exact brackets for:
  eps* = 2(1-g-g^2)/(3+2g)      (Theorem B threshold; decreasing in g)
  1 + g + eps* = (5+3g)/(3+2g)  (NEW closed form; derivative -1/(3+2g)^2)
  i_feed/(R) bound (1-g)/(1+g) < 1/2   (deadline never binding)
  d_feedend lower bound (2g R + D0)/(1+g) - 1
  K1 = ceil((2 Lambda - 1)/(1 - g))    (cell-1 extinction row budget)
  R_b (budget threshold) and R_1 (ENTRY(mu*sqrt R) threshold), Lambda=2^11.
Output: amm12592_S_cone_constants_boxeph.json + stdout report.
"""
import json, os, sys
from fractions import Fraction
from math import isqrt

RESULTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "05-knowledge", "results")

def fib_luc(n):
    """(F_n, L_n) by fast doubling, exact ints."""
    if n == 0:
        return 0, 2
    F, L = fib_luc(n // 2)
    # F_{2k} = F_k L_k ; L_{2k} = L_k^2 - 2(-1)^k
    F2 = F * L
    L2 = L * L - (2 if n // 2 % 2 == 0 else -2)
    if n % 2 == 0:
        return F2, L2
    # step +1: F_{m+1} = (F_m + L_m)/2, L_{m+1} = (5F_m + L_m)/2
    return (F2 + L2) // 2, (5 * F2 + L2) // 2

def g_gt(p, q, Fq, Lq):
    """Exact test  g = log_5 phi^2 > p/q  <=>  phi^(2q) > 5^p, phi^(2q) =
    (L_{2q} + F_{2q} sqrt5)/2, via integer squaring."""
    rhs = 2 * 5 ** p - Lq
    if rhs < 0:
        return True
    return 5 * Fq * Fq > rhs * rhs

def main():
    sys.setrecursionlimit(100)
    M = 1 << 20
    F, L = fib_luc(2 * M)        # F_{2^21}, L_{2^21}
    p_lo, p_hi = 627035, 627036
    a = g_gt(p_lo, M, F, L)
    b = g_gt(p_hi, M, F, L)
    assert a and not b, ("g sandwich FAILED", a, b)
    g_lo = Fraction(p_lo, M)
    g_hi = Fraction(p_hi, M)
    print(f"[PASS] fresh g sandwich: {p_lo}/2^20 < g < {p_hi}/2^20 "
          f"({float(g_lo):.10f}, {float(g_hi):.10f})")

    def epsstar(g):
        return 2 * (1 - g - g * g) / (3 + 2 * g)

    # eps* decreasing in g: f'(g) = -2(5+6g+2g^2)/(3+2g)^2 < 0 (symbolic).
    es_lo, es_hi = epsstar(g_hi), epsstar(g_lo)
    assert es_lo < es_hi
    print(f"[PASS] eps* in ({float(es_lo):.10f}, {float(es_hi):.10f}); "
          f"eps* < 1/32: {es_hi < Fraction(1,32)}")

    # closed form check: 1 + g + eps*(g) == (5+3g)/(3+2g) identically
    for gt in (g_lo, g_hi, Fraction(3, 5), Fraction(1, 2)):
        assert 1 + gt + epsstar(gt) == Fraction(5 + 3 * gt, 1) / (3 + 2 * gt)
    # decreasing in g with derivative -1/(3+2g)^2  =>  bracket:
    cstar_lo = (5 + 3 * g_hi) / (3 + 2 * g_hi)
    cstar_hi = (5 + 3 * g_lo) / (3 + 2 * g_lo)
    assert cstar_lo < cstar_hi
    print(f"[PASS] 1 + g + eps* = (5+3g)/(3+2g) in "
          f"({float(cstar_lo):.10f}, {float(cstar_hi):.10f})")
    print(f"       width = {float(cstar_hi - cstar_lo):.3e}")

    # deadline: (1-g)/(1+g) < 1/2  <=>  g > 1/3
    assert g_lo > Fraction(1, 3)
    ifeed_frac_ub = (1 - g_lo) / (1 + g_lo)
    print(f"[PASS] i_feed/R <= (1-g)/(1+g) < {float(ifeed_frac_ub):.6f} < 1/2")

    LAM = 1 << 11
    # K1 = ceil((2*LAM - 1)/(1 - g_hi))  (cell-1: needs 2(Lambda-1)+1 delta0
    # rows; #delta0 rows in k consecutive >= k(1-g) - 1)
    num = (2 * LAM - 1)
    K1 = -((-num * (1 + 0)) // 1)  # placeholder
    K1 = (num * g_hi.denominator + (g_hi.denominator - g_hi.numerator) - 1) \
        // (g_hi.denominator - g_hi.numerator)
    K1 = int(K1)
    print(f"K1 (cell-1 extinction budget) = {K1} rows  (2Lam = {2*LAM})")

    out = {"g_lo": [p_lo, M], "g_hi": [p_hi, M],
           "epsstar": [[es_lo.numerator, es_lo.denominator],
                       [es_hi.numerator, es_hi.denominator]],
           "one_plus_g_plus_epsstar_eq": "(5+3g)/(3+2g)",
           "cstar_bracket": [[cstar_lo.numerator, cstar_lo.denominator],
                             [cstar_hi.numerator, cstar_hi.denominator]],
           "cstar_bracket_float": [float(cstar_lo), float(cstar_hi)],
           "K1": K1, "Lambda": LAM}

    # R_b: least dyadic R with ifeed_ub(R) + 2 + layer(33) + K1 <= R - 2 and
    # (for both eps) window certificates (Thm B: R >= 162 suffices).
    for eps_num, eps_den in ((1, 32), (1, 16)):
        eps = Fraction(eps_num, eps_den)
        Rb = None
        for k in range(7, 36):
            R = 1 << k
            D0 = -((-R * eps_num) // eps_den)
            ifeed_ub = (R * (1 - g_lo) - D0) / (1 + g_lo)   # rational UB
            if ifeed_ub + 2 + 33 + K1 <= R - 2 and R >= 162:
                Rb = R
                break
        # R_1: least dyadic R with 6*LAM*(MU*ceil(sqrt R)+2) <= d_fe lower
        # bound (2 g_lo R + D0)/(1+g_hi) - 1
        MU = 2
        R1 = None
        for k in range(7, 40):
            R = 1 << k
            D0 = -((-R * eps_num) // eps_den)
            s = isqrt(R)
            cs = s if s * s == R else s + 1
            dfe_lb = (2 * g_lo * R + D0) / (1 + g_hi) - 1
            if 6 * LAM * (MU * cs + 2) <= dfe_lb:
                R1 = R
                break
        print(f"eps = {eps_num}/{eps_den}:  R_b = {Rb} = 2^{Rb.bit_length()-1}"
              f"   R_1(Lambda=2^11, mu=2) = {R1} = 2^{R1.bit_length()-1}")
        out[f"Rb_eps_{eps_num}_{eps_den}"] = Rb
        out[f"R1_eps_{eps_num}_{eps_den}"] = R1

    # ---- Corollary (a-priori F4): from F3 alone, a_t <= C(D0-1, t+1) for
    # every t >= 1 and a_0 <= D0-1 (F2).  Exact clock bounds:
    #   K1c: smallest K with sum_k 2*[s_k + 2(k-1)] >= C(Dub-1,2)*... using
    #        sum s_k >= g_lo K^2/2 - K  (s_k = floor(g_lo k)):
    #        (g_lo+2)K^2 - 4K >= Dub*(Dub-1)/2       [a_1 <= C(Dub-1,2)]
    #   K2c: (g_lo K^2/2 - K)*(2 Dlb - 3) >= C(Dub-1,3) [a_2 <= C(Dub-1,3)]
    #   and T_t <= T_2 for all t >= 2 (need decreasing, rate increasing).
    # Verify i_ub + 2 + max(K1c, K2c) <= R - 2 for dyadic R (both eps).
    print("\nCorollary a-priori-F4 verification (exact rationals):")
    for eps_num, eps_den in ((1, 32), (1, 16)):
        R0 = None
        bad = []
        for kk in range(7, 41):
            R = 1 << kk
            D0 = -((-R * eps_num) // eps_den)
            i_ub = int((R * (1 - g_lo) - D0) / (1 + g_lo)) + 2 + 64
            Dlb = int((2 * g_lo * R + D0) / (1 + g_hi)) - 1
            Dub = int(g_hi * (R + i_ub)) + D0 + 1
            need1 = Dub * (Dub - 1) // 2
            need2 = Dub * (Dub - 1) * (Dub - 2) // 6
            # binary search exact K's
            def k1ok(K):
                return (g_lo + 2) * K * K - 4 * K >= need1
            def k2ok(K):
                return (g_lo * K * K / 2 - K) * (2 * Dlb - 3) >= need2
            def bs(ok):
                lo, hi = 1, 4 * R
                if not ok(hi):
                    return None
                while lo < hi:
                    mid = (lo + hi) // 2
                    if ok(mid):
                        hi = mid
                    else:
                        lo = mid + 1
                return lo
            K1c, K2c = bs(k1ok), bs(k2ok)
            ok = (K1c is not None and K2c is not None and
                  i_ub + max(K1c, K2c) <= R - 2 and
                  2 * i_ub <= R - 2)
            if ok and R0 is None:
                R0 = R
            if not ok:
                bad.append(R)
                R0 = None      # require ok for ALL R >= R0
        # find final stable R0: last bad + one doubling
        R0 = (max(bad) * 2) if bad else (1 << 7)
        print(f"  eps={eps_num}/{eps_den}: a-priori F4 holds for all dyadic "
              f"R >= {R0} = 2^{R0.bit_length()-1}; failures at {bad if bad else 'none'}")
        out[f"apriori_F4_R0_eps_{eps_num}_{eps_den}"] = R0
    json.dump(out, open(os.path.join(
        RESULTS, "amm12592_S_cone_constants_boxeph.json"), "w"), indent=1)
    print("constants JSON written.")

if __name__ == "__main__":
    main()
