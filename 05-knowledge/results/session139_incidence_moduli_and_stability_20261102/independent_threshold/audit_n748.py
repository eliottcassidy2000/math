#!/usr/bin/env python3
"""Independent exact controls for Session139's n >= 748 analytic theorem.

Every check stays active under python -O.  --threshold 747 is a deliberately
false extension, rejected by the shared-wedge-budget endpoint certificate.
Finite checks support (but do not replace) the uniform arguments in AUDIT.md.
"""
import argparse
from fractions import Fraction as F
import json
from math import factorial
import sys


class CheckFailure(Exception):
    pass


def run(n):
    evidence = {}

    def check(name, condition, value=None):
        if not condition:
            raise CheckFailure(f"{name}: {value}")
        evidence[name] = True if value is None else str(value)

    c, eta, delta = F(249999, 1000000), F(1, 55), F(1, 3500)
    k_upper = F(8, n - 1)
    denominator = c - 4 * eta - delta
    check("light_coefficient_positive", denominator > 0, denominator)
    C = 4 / denominator
    h = eta - C * eta**2
    check("order_domain", n >= 225, n)
    check("kappa_positive_domain", n // 2 >= 2)
    exp4partial = sum((F(4)**j / factorial(j) for j in range(8)), F(0))
    check("exp4_partial_identity", exp4partial == F(16319, 315), exp4partial)
    check("exp4_partial_exceeds_50", exp4partial > 50, exp4partial - 50)
    check("small_variance_excluded", F(49, 200) > k_upper, F(49, 200) - k_upper)
    v_lower = F(49 * (n - 1), 1600)
    check("variance_bootstrap_exceeds_4", v_lower > 4, v_lower - 4)
    g_lower = (1 - F(1, 50**4)) / 4
    check("nonlinear_c_bootstrap", g_lower > c, g_lower - c)
    check("variance_upper_uniform_domain", n - 1 >= 224)
    check("variance_upper_endpoint", F(n, 7) >= F(n, 8) + F(4 * n, n - 1),
          F(n, 7) - F(n, 8) - F(4 * n, n - 1))
    slope = c / 8 - eta / 2 - F(1, 60)
    surplus = c * (n - 1) / 8 - F(4 * n, n - 1) - eta * n / 2 - F(n, 60)
    check("heavy_mass_uniform_slope", slope > 0, slope)
    check("heavy_mass_endpoint_surplus", surplus > 0, surplus)
    check("exponential_tail_endpoint", (2000 * n)**15 * 3**n < 4**n)
    check("exponential_tail_uniform_derivative", F(1, 60) > F(1, n),
          F(1, 60) - F(1, n))
    check("light_error_delta", 4 * F(1, 2000) * F(1, 7) == delta, delta)
    check("shared_budget_interior_minimum", 0 < C * eta < F(1, 2), F(1, 2) - C * eta)
    check("shared_budget_strict_contradiction", h > k_upper, h - k_upper)
    check("no_light_closure", eta > k_upper, eta - k_upper)
    # Polynomial coefficients are checked, not sample values of x or D.
    x0 = 1 - 2 * C * eta
    left = (1 / (4 * C) - h, eta - 1 / (2 * C), 1 / (4 * C))
    right = (x0**2 / (4 * C), -x0 / (2 * C), 1 / (4 * C))
    check("shared_budget_square_identity", left == right)
    nk = n * k_upper
    original = (F(0), F(-1), 1 / (2 * nk))
    completed = (nk**2 / (2 * nk) - nk / 2,
                 -2 * nk / (2 * nk), 1 / (2 * nk))
    check("heavy_mass_square_identity", original == completed)
    # Explicitly verify both parity coefficient bounds at the endpoint pair.
    for order in (n, n + 1):
        r = order // 2
        kappa = F(4, r) * (1 - F(3, 4)**r - (F(1, 4)**r if order % 2 == 0 else 0))
        check(f"parity_coefficient_{order}", 0 < kappa <= F(8, order - 1)
              and kappa * r <= 4)
    return {"status": "PASS", "threshold": n,
            "constants": {"c": str(c), "eta": str(eta), "delta": str(delta),
                          "C": str(C), "h": str(h)}, "active_checks": evidence}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threshold", type=int, default=748)
    args = parser.parse_args()
    try:
        result = run(args.threshold)
    except CheckFailure as exc:
        print(json.dumps({"status": "REJECTED", "reason": str(exc)}, sort_keys=True))
        return 1
    print(json.dumps(result, sort_keys=True, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
