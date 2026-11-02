#!/usr/bin/env python3
"""Independent exact controls for the light-deficit n >= 526 argument.

All checks remain active under python -O. --threshold 525 is a deliberately
false extension of this certificate and must be rejected by its final margin.
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

    r = n // 2
    c, d = F(249999, 1000000), F(487, 2000)
    t, lam = d / 8, d / 64
    k_upper, k_rough = F(4, r), F(8, n - 1)
    check("order_domain", n >= 225 and r > 18)
    check("small_variance_excluded", F(49, 200) > k_rough, F(49, 200) - k_rough)
    e4 = sum((F(4)**j / factorial(j) for j in range(8)), F(0))
    check("exp4_partial_exceeds_50", e4 > 50, e4)
    check("variance_bootstrap_exceeds_4", F(49 * (n - 1), 1600) > 4,
          F(49 * (n - 1), 1600) - 4)
    check("nonlinear_c_bootstrap", (1 - F(1, 50**4)) / 4 > c,
          (1 - F(1, 50**4)) / 4 - c)
    rh_lower = c * (n - 1) / 2 - F(16 * n, n - 1) - 2 * t * n
    check("heavy_count_endpoint_exceeds_17", rh_lower > 17, rh_lower - 17)
    check("heavy_count_uniform_slope", c / 2 - 2 * t > 0, c / 2 - 2 * t)
    A = F(r) + F(33, 2) + F(8, r)
    check("variance_upper_even_identity", F(2 * r, 2) + 2 * (2 * r) * F(4, r) == r + 16)
    check("variance_upper_odd_identity", F(2 * r + 1, 2) + 2 * (2 * r + 1) * F(4, r) == A)
    tail_left = A * F(3, 4)**(r - 1)
    tail_right = A * F(3, 4)**18 / (r - 18)
    check("tail_budget_positive", c - d > 0, c - d)
    check("tail_left_endpoint", tail_left < c - d)
    check("tail_right_endpoint", tail_right < c - d, c - d - tail_right)
    check("tail_left_uniform_derivative", F(1, r) < F(1, 4), F(1, 4) - F(1, r))
    check("tail_right_uniform_derivative", 144 <= 16 * r)
    check("tail_right_derivative_identity",
          (1 - F(8, r*r)) * (r - 18) - A == -F(69, 2) - F(16, r) + F(144, r*r))
    check("light_deficit_positive", d - 4 * t > 0, d - 4 * t)
    check("local_square_constant", lam * d == t*t)
    check("local_square_linear_w", 4 * lam == t / 2)
    check("local_square_remainder_nonnegative", t > 0 and 2*t-t >= 0)
    check("global_budget_identity", t - 4 * lam == d / 16, t - 4 * lam)
    check("final_strict_margin", d / 16 > k_upper, d / 16 - k_upper)
    check("no_light_closure", t > k_upper, t - k_upper)
    # Check the local bivariate polynomial identity coefficient by coefficient.
    # Indices are exponents (S,w), not sample evaluation points.
    phi_minus_tS_plus_deficit = {(2, 0): F(1, 4), (1, 1): F(1, 2),
        (1, 0): -t, (0, 1): -4*lam, (0, 0): lam*d}
    square_plus_remainder = {(2, 0): F(1, 4), (1, 1): F(1, 2),
        (1, 0): -t, (0, 1): -t + t/2, (0, 0): t*t}
    check("local_square_polynomial_identity", phi_minus_tS_plus_deficit == square_plus_remainder)
    for order in (n, n + 1):
        R = order // 2
        kappa = F(4, R) * (1 - F(3, 4)**R - (F(1, 4)**R if order % 2 == 0 else 0))
        check(f"parity_bound_{order}", 0 < kappa <= F(4, R))
    return {"status": "PASS", "threshold": n, "R": r,
            "constants": {"c": str(c), "d": str(d), "t": str(t), "lambda": str(lam)},
            "active_checks": evidence}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threshold", type=int, default=526)
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
