#!/usr/bin/env python3
"""Exact controls for odd-square shape measures; all checks survive -O."""
import argparse
from fractions import Fraction
import json
from math import asin, atan, gcd, isqrt, pi, sqrt
from pathlib import Path


CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def ratio(a, b):
    value = Fraction(a, b)
    return {"numerator": value.numerator, "denominator": value.denominator}


def arithmetic_sieve(limit):
    phi = list(range(limit + 1))
    mu = [1] * (limit + 1)
    omega = [0] * (limit + 1)
    for p in range(2, limit + 1):
        if omega[p] == 0:
            for n in range(p, limit + 1, p):
                phi[n] -= phi[n] // p
                mu[n] = -mu[n]
                omega[n] += 1
            for n in range(p * p, limit + 1, p * p):
                mu[n] = 0
    return phi, mu, omega


def euclid_triples(limit):
    triples = []
    for m in range(2, isqrt(limit) + 1):
        top = min(m - 1, isqrt(limit - m * m))
        for n in range(1 if m % 2 == 0 else 2, top + 1, 2):
            if gcd(m, n) == 1:
                A, B, c = m * m - n * n, 2 * m * n, m * m + n * n
                triples.append((A, B, c, m, n))
    return sorted(triples, key=lambda t: (t[2], t[0], t[1]))


def odd_root_triples(limit):
    found = set()
    for u in range(3, isqrt(2 * limit - 1) + 1, 2):
        top = min(u - 2, isqrt(2 * limit - u * u))
        for v in range(1, top + 1, 2):
            if gcd(u, v) == 1:
                found.add((u * v, (u * u - v * v) // 2, (u * u + v * v) // 2))
    return found


def direct_triangles(limit):
    found = set()
    for c in range(2, limit + 1):
        for a in range(1, isqrt(c * c // 2) + 1):
            b = isqrt(c * c - a * a)
            if a * a + b * b == c * c and gcd(a, b) == 1:
                found.add((a, b, c) if a % 2 else (b, a, c))
    return found


def fibre_controls(phi, omega):
    records = []
    selected = {3, 5, 9, 15, 105, 315, 999, 1001}
    aggregate = 0
    aggregate_angle = 0
    for u in range(3, 1002, 2):
        values = [v for v in range(1, u, 2) if gcd(u, v) == 1]
        N = len(values)
        require(2 * N == phi[u], f"fibre count u={u}")
        require(N > 0, "nonempty fibre")
        max_num = max(max(abs(j * u - v * N), abs((j - 1) * u - v * N))
                      for j, v in enumerate(values, 1))
        require(2 * max_num <= (2 ** omega[u]) * u, f"fibre discrepancy u={u}")
        require(4 * u * 4 ** omega[u] <= 15 * phi[u] ** 2, f"uniform discrepancy bound u={u}")
        angle_count = sum(3 * v <= u or 2 * v >= u for v in values)
        direct_angle = 0
        for v in values:
            a, b = u * v, (u * u - v * v) // 2
            small, large = sorted((a, b))
            direct_angle += 4 * small <= 3 * large
        require(angle_count == direct_angle, f"folded angle preimages u={u}")
        # Difference from G(atan(3/4))=5/6 bounded by two CDF errors.
        require(abs(6 * angle_count - 5 * N) <= 6 * 2 ** omega[u], f"folded fibre error u={u}")
        aggregate += N
        aggregate_angle += angle_count
        if u in selected:
            records.append({"u": u, "size": N,
                            "slope_Kolmogorov_discrepancy": ratio(max_num, N * u),
                            "proved_bound": ratio(2 ** omega[u], phi[u]),
                            "angle_at_atan_3_over_4_count": angle_count,
                            "angle_at_atan_3_over_4_fraction": ratio(angle_count, N)})
    require(4 * 15 * 4 ** omega[15] == 15 * phi[15] ** 2, "u15 bound equality")
    return {"universe": "all odd u=3..1001, every admissible odd v", "fibres": 500,
            "selected": records, "aggregate_outer_height_1001": {
                "triples": aggregate, "angle_at_atan_3_over_4_count": aggregate_angle,
                "asymptotic_count_diagnostic": 1001 ** 2 / pi ** 2,
                "limit_angle_CDF": ratio(5, 6)}}


def parity_controls(mu):
    rows = []
    for M in [10, 100, 1000]:
        coprime = opposite = oddodd = 0
        for a in range(1, M + 1):
            for b in range(1, M + 1):
                if gcd(a, b) == 1:
                    coprime += 1
                    opposite += (a - b) % 2 == 1
                    oddodd += a % 2 == 1 and b % 2 == 1
        mob_all = sum(mu[d] * (M // d) ** 2 for d in range(1, M + 1))
        mob_opposite = sum(mu[d] * 2 * ((M // d + 1) // 2) * ((M // d) // 2)
                           for d in range(1, M + 1, 2))
        mob_odd = sum(mu[d] * ((M // d + 1) // 2) ** 2 for d in range(1, M + 1, 2))
        require((coprime, opposite, oddodd) == (mob_all, mob_opposite, mob_odd), "independent parity Mobius count")
        require(coprime == opposite + oddodd, "coprime parity partition")
        rows.append({"M": M, "pairs": M * M, "coprime": coprime,
                     "coprime_opposite_parity": opposite, "coprime_odd_odd": oddodd,
                     "opposite_parity_denominator": 2 * ((M + 1) // 2) * (M // 2),
                     "odd_odd_denominator": ((M + 1) // 2) ** 2})
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path(__file__).with_suffix(".json"))
    args = parser.parse_args()
    phi, mu, omega = arithmetic_sieve(1001)
    triples = euclid_triples(1_000_000)
    small_euclid = {t[:3] for t in triples if t[2] <= 10_000}
    require(small_euclid == odd_root_triples(10_000), "independent odd-root enumeration")
    require({t[:3] for t in triples if t[2] <= 300} == direct_triangles(300), "independent integer-leg enumeration")
    require(sum(t[2] <= 1000 for t in triples) == 158, "inherited c1000 positive control")
    require(triples[0][:3] == (3, 4, 5), "first primitive triple")
    require((9, 12, 15) not in small_euclid, "nonprimitive hostile")
    for A, B, c, m, n in triples:
        require(A * A + B * B == c * c, "right triangle")
        require(gcd(A, B) == 1 and A % 2 == 1 and B % 2 == 0, "primitive labels")
        require(c + B == (m + n) ** 2 and c - B == (m - n) ** 2, "odd-square pair")
        if c <= 1000:
            a, b = sorted((A, B))
            e, d, l = Fraction(a * a, c * c), Fraction(b * b, c * c), Fraction(a * b, c * c)
            r = e / d
            require(e + d == 1 and l * l == e * d, "exact altitude projection identities")
            require(r == Fraction(a, b) ** 2 and l * l == r / (1 + r) ** 2, "exact shape map")
            require((e - Fraction(1, 2)) ** 2 + l * l == Fraction(1, 4), "unit-diameter circle")
    require(Fraction(3 * 4, 5 * 5) == Fraction(12, 25), "3-4-5 normalized altitude")
    require(Fraction(3 * 4, 5) ** 2 != 2, "sqrt2 altitude hostile")
    heights = []
    for X in [100, 1000, 10_000, 100_000, 1_000_000]:
        current = [t for t in triples if t[2] <= X]
        heights.append({"hypotenuse_height": X, "primitive_triples": len(current),
                        "X_over_2pi_diagnostic": X / (2 * pi),
                        "count_divided_by_X": ratio(len(current), X)})
    shape_cdfs = []
    for numerator in range(1, 9):
        count = sum(8 * min(t[0], t[1]) <= numerator * max(t[0], t[1]) for t in triples)
        shape_cdfs.append({"tan_theta_cutoff": ratio(numerator, 8), "count": count,
                           "finite_fraction": ratio(count, len(triples)),
                           "angle_uniform_limit_diagnostic": 4 * atan(numerator / 8) / pi})
    epsilon_rows = []
    for eps in map(Fraction, ["1/10", "1/4", "2/5", "9/20", "19/40", "49/100", "1/2"]):
        num, den = eps.numerator, eps.denominator
        both = ratio_only = altitude_only = 0
        for A, B, c, _, _ in triples:
            a, b = sorted((A, B))
            r_pass = den * a * a > num * b * b
            l_pass = den * a * b > num * c * c
            ratio_only += r_pass
            altitude_only += l_pass
            both += r_pass and l_pass
        binding = "empty" if eps >= Fraction(1, 2) else "ratio" if eps * (1 + eps) ** 2 <= 1 else "altitude"
        require(both == (0 if binding == "empty" else ratio_only if binding == "ratio" else altitude_only), "epsilon nesting and binding coordinate")
        limit = 0.0 if binding == "empty" else 1 - 4 * max(atan(sqrt(float(eps))), asin(2 * float(eps)) / 2) / pi
        epsilon_rows.append({"epsilon": ratio(num, den), "binding": binding,
                             "ratio_count": ratio_only, "altitude_count": altitude_only,
                             "joint_count": both, "limit_fraction_diagnostic": limit})
    fibres = fibre_controls(phi, omega)
    parities = parity_controls(mu)
    result = {
        "status": "FINITE-EXACT PASS",
        "scope": {"max_hypotenuse": 1_000_000, "all_triples_enumerated": len(triples),
                  "independent_odd_root_max_hypotenuse": 10_000,
                  "independent_integer_leg_max_hypotenuse": 300,
                  "exact_Fraction_shape_max_hypotenuse": 1000,
                  "all_fixed_outer_fibres": "odd u=3..1001",
                  "asymptotic_decimals": "diagnostics only; limits proved in note"},
        "heights": heights, "shape_cdfs": shape_cdfs, "epsilon_rows": epsilon_rows,
        "fibres": fibres, "parity_boxes": parities,
        "hostiles": {"3_4_5_altitude": ratio(12, 5), "3_4_5_normalized_altitude": ratio(12, 25),
                     "3_4_5_angle_outer_fibre_limit_CDF": ratio(5, 6),
                     "3_4_5_angle_hypotenuse_limit_CDF_diagnostic": 4 * atan(3 / 4) / pi,
                     "sampling_measures_differ": True},
        "checks": CHECKS
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": result["status"], "checks": CHECKS,
                      "triples": len(triples), "output": str(args.output)}))


if __name__ == "__main__":
    main()
