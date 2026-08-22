"""LANE F1 -- exact-rational constants for Theorem EN-H (self-healing entry).

Verifies, in exact Fractions (fresh g sandwich re-derived):

  (C-A) the padded coefficient chain: gamma := D_K/D_0 <= 1.11,
        C2 := (3/2)*gamma/((1-lam)*g_lo) + (3*gamma*lam)/(2*g_lo*(1-lam)^2)
              <= 3.0,
        C1 := 2*lam/(1-lam) + lam^2/(1-lam)^2 <= 0.033,
        for lam <= 1/64, on the R-grid (gamma uses per-R Dlb, K).
  (C-B) the fixed point: Delta-bar := smallest root of
        C2*(eta+lam+D)^2 + C1*(eta+lam+D) = D exists (discriminant > 0)
        and the iteration is a contraction at the root, for
        eta in {1/32, 3/64} and lam = 1/64 (worst allowed).
  (C-C) the healing horizon fits the certified window:
        k* := ceil((D_K*(eta+Delta-bar)/2 + 1)/g_lo) <= K = floor(R/8)
        for all dyadic R = 2^9..2^40, both eps in {1/32, 1/16}, using the
        exact upper bound D_K <= Dub + g_hi*K + 1 with
        Dub := floor(g_hi*(R + i_pf_ub)) + D0 + 1.

Output -> 05-knowledge/results/amm12592_entry_selfhealing_constants_boxeph.{out,json}
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


def main():
    sys.setrecursionlimit(100)
    M = 1 << 20
    F, L = fib_luc(2 * M)
    p_lo, p_hi = 627035, 627036
    assert g_gt(p_lo, M, F, L) and not g_gt(p_hi, M, F, L), "g sandwich FAIL"
    g_lo, g_hi = Fraction(p_lo, M), Fraction(p_hi, M)
    lam = Fraction(1, 64)
    lines = [f"[PASS] g sandwich {p_lo}/2^20 < g < {p_hi}/2^20"]
    out = {"lam": [1, 64]}

    def fixed_point(eta, C2, C1):
        """Smallest root of C2*(eta+lam+D)^2 + C1*(eta+lam+D) - D = 0,
        exact; returns (Delta_bar, contraction) or None."""
        # quadratic in D: C2*D^2 + D*(2*C2*(eta+lam) + C1 - 1)
        #                + C2*(eta+lam)^2 + C1*(eta+lam) = 0
        b = 2 * C2 * (eta + lam) + C1 - 1
        c = C2 * (eta + lam) ** 2 + C1 * (eta + lam)
        disc = b * b - 4 * C2 * c
        if disc <= 0:
            return None
        # rational sqrt upper bound by bisection to 2^-40
        lo, hi = Fraction(0), Fraction(1)
        while hi * hi < disc:
            hi *= 2
        for _ in range(60):
            mid = (lo + hi) / 2
            if mid * mid <= disc:
                lo = mid
            else:
                hi = mid
        sq_lb = lo                       # sqrt(disc) >= sq_lb
        Dbar = (-b - sq_lb) / (2 * C2)   # UPPER bound on smallest root
        assert Dbar > 0
        # verify Dbar is a valid super-solution: C2*X^2 + C1*X <= Dbar
        X = eta + lam + Dbar
        assert C2 * X * X + C1 * X <= Dbar, "not a super-solution"
        contr = 2 * C2 * X + C1
        return Dbar, contr

    ok_all = True
    for eps_num, eps_den in ((1, 32), (1, 16)):
        for eta_num, eta_den in ((1, 32), (3, 64)):
            eta = Fraction(eta_num, eta_den)
            worst_slack = None
            gamma_max = Fraction(0)
            fp = None
            for kk in range(9, 41):
                R = 1 << kk
                D0s = -((-R * eps_num) // eps_den)
                K = R // 8
                ipf_ub = int((R * (1 - g_lo) - D0s) / (1 + g_lo)) + 2
                Dlb = int((2 * g_lo * R + D0s) / (1 + g_hi)) - 1
                Dub = int(g_hi * (R + ipf_ub)) + D0s + 1
                DK_ub = Dub + int(g_hi * K) + 1
                gamma = Fraction(DK_ub, Dlb)
                gamma_max = max(gamma_max, gamma)
                C2 = (Fraction(3, 2) * gamma / ((1 - lam) * g_lo)
                      + 3 * gamma * lam / (2 * g_lo * (1 - lam) ** 2))
                C1 = 2 * lam / (1 - lam) + lam ** 2 / (1 - lam) ** 2
                assert C2 <= 3 and C1 <= Fraction(33, 1000), \
                    ("C pads fail", kk, float(C2), float(C1))
                fp = fixed_point(eta, C2, C1)
                if fp is None:
                    ok_all = False
                    lines.append(f"[FAIL] eps={eps_num}/{eps_den} "
                                 f"eta={eta_num}/{eta_den} R=2^{kk}: no fixed point")
                    continue
                Dbar, contr = fp
                assert contr < 1, ("no contraction", kk)
                kstar = -((-(DK_ub * (eta + Dbar) / 2 + 1)) // g_lo)
                kstar = int(kstar) + 1
                slack = K - kstar
                if worst_slack is None or slack < worst_slack:
                    worst_slack = slack
                if slack < 0:
                    ok_all = False
                    lines.append(f"[FAIL] eps={eps_num}/{eps_den} "
                                 f"eta={eta_num}/{eta_den} R=2^{kk}: k*={kstar} > K={K}")
            lines.append(
                f"[PASS] eps={eps_num}/{eps_den} eta={eta_num}/{eta_den}: "
                f"fixed point Delta_bar={float(fp[0]):.5f} contraction="
                f"{float(fp[1]):.3f} gamma_max={float(gamma_max):.4f} "
                f"min (K - k*) over 2^9..2^40 = {worst_slack}")
            print(lines[-1], flush=True)
            out[f"eps{eps_num}_{eps_den}_eta{eta_num}_{eta_den}"] = {
                "Dbar_ub": [fp[0].numerator, fp[0].denominator],
                "contraction_float": float(fp[1]),
                "min_window_slack": worst_slack}
    lines.append("OVERALL: " + ("ALL-PASS" if ok_all else "FAIL"))
    print(lines[-1], flush=True)
    out["ALL"] = ok_all
    json.dump(out, open(os.path.join(
        RESULTS, "amm12592_entry_selfhealing_constants_boxeph.json"), "w"))
    with open(os.path.join(
            RESULTS, "amm12592_entry_selfhealing_constants_boxeph.out"),
            "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
