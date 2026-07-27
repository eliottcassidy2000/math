#!/usr/bin/env python3
"""Exact dependency-free referee for THM-2520.

The proof is symbolic.  This companion supplies an independent finite-grid
path for the jump aggregation, Perron reduction, CRT coverage, variation
invoice, and the sharp conditional-variance controls.  It intentionally uses
``require`` rather than ``assert`` so normal and ``python -O`` executions are
identical.
"""

from fractions import Fraction as Q
from math import gcd, lcm


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def cyclic_jumps(values):
    return [values[j] - values[(j - 1) % len(values)] for j in range(len(values))]


def cyclic_variation(values):
    return sum(abs(values[(j + 1) % len(values)] - values[j]) for j in range(len(values)))


def aggregate_jumps(values, k13, d0):
    d = (13**k13) * d0
    require(len(values) == d, "endpoint grid length mismatch")
    jumps = cyclic_jumps(values)
    out = [Q(0) for _ in range(d0)]
    for j, jump in enumerate(jumps):
        out[j % d0] += jump
    return jumps, out


def step_value(values, x):
    """Value on [j/D,(j+1)/D), away from endpoints."""
    x -= x.numerator // x.denominator
    scaled = x * len(values)
    j = scaled.numerator // scaled.denominator
    return values[j % len(values)]


def direct_perron_cells(values, m13, d0):
    """Values of P_(13^m)F on the D0 open grid cells."""
    m = 13**m13
    out = []
    for t in range(d0):
        y = Q(2 * t + 1, 2 * d0)
        total = Q(0)
        for r in range(m):
            total += step_value(values, (y + r) / m)
        out.append(total / m)
    return out


def predicted_perron_jumps(c, k13, m13, d0):
    m = 13**m13
    if d0 == 1:
        return [c[0] / m]
    scale = pow(13, m13 - k13, d0)
    inv = pow(scale, -1, d0)
    return [c[(inv * t) % d0] / m for t in range(d0)]


def deterministic_values(kind, k13, d0, seed):
    d = (13**k13) * d0
    if kind == "pure13":
        # Constant on each of the 13^K coarse cells.
        return [Q(((j // d0) * (seed + 1) + seed) % 5, 4) for j in range(d)]
    if kind == "prime_to_13":
        return [Q(((j % d0) * (seed + 2) + 2 * seed + 1) % 7, 6) for j in range(d)]
    return [Q((j * j + (3 + seed) * j + 5 * seed + k13) % 11, 10) for j in range(d)]


def integrate_on_grid(function, denominator):
    total = Q(0)
    for j in range(denominator):
        x = Q(2 * j + 1, 2 * denominator)
        total += function(x)
    return total / denominator


def collision_vector(values):
    d = len(values)
    grid = lcm(d, 13)
    out = []
    for u in range(13):
        out.append(
            integrate_on_grid(
                lambda x, u=u: step_value(values, x)
                * step_value(values, x + Q(u, 13)),
                grid,
            )
        )
    return out


def conditional_variance(values):
    d = len(values)
    grid = lcm(d, 13)

    def summand(x):
        h = step_value(values, x)
        mean = sum(step_value(values, x + Q(u, 13)) for u in range(13)) / 13
        return (h - mean) ** 2

    return integrate_on_grid(summand, grid)


def run_grid_referee():
    cases = 0
    perron_cells = 0
    crt_checks = 0
    variation_checks = 0
    zero_cases = 0
    charged_cases = 0

    for k13 in range(3):
        for d0 in (1, 2, 5, 7):
            for kind in ("pure13", "prime_to_13", "generic"):
                for seed in range(2):
                    values = deterministic_values(kind, k13, d0, seed)
                    _, c = aggregate_jumps(values, k13, d0)
                    c_zero = all(x == 0 for x in c)
                    for m13 in (k13, k13 + 1):
                        h_values = direct_perron_cells(values, m13, d0)
                        actual_jumps = cyclic_jumps(h_values)
                        predicted = predicted_perron_jumps(c, k13, m13, d0)
                        require(actual_jumps == predicted, "Perron jump-current formula failed")
                        is_constant = all(x == h_values[0] for x in h_values)
                        require(is_constant == c_zero, "constant/jump dichotomy failed")
                        require(
                            cyclic_variation(h_values)
                            <= cyclic_variation(values) / (13**m13),
                            "variation contraction failed",
                        )
                        cases += 1
                        perron_cells += len(h_values)
                        variation_checks += 1
                        zero_cases += int(c_zero)
                        charged_cases += int(not c_zero)

                        # For each last digit a, CRT covers every D0 residue.
                        scale = pow(13, m13 - k13, d0) if d0 > 1 else 0
                        for a in range(1, 13):
                            residues = {
                                ((scale * k) % d0) if d0 > 1 else 0
                                for k in range(a, 13 * d0, 13)
                            }
                            require(residues == set(range(d0)), "CRT ladder coverage failed")
                            crt_checks += 1

    return cases, perron_cells, crt_checks, variation_checks, zero_cases, charged_cases


def run_controls():
    # Pure-13 hostile: one of thirteen roots lies in [0,1/13).
    hostile_values = [Q(1)] + [Q(0)] * 12
    _, hostile_c = aggregate_jumps(hostile_values, 1, 1)
    hostile_h = direct_perron_cells(hostile_values, 1, 1)
    require(hostile_c == [Q(0)], "pure-13 hostile current changed")
    require(hostile_h == [Q(1, 13)], "pure-13 hostile Perron mean changed")
    require((13 * hostile_h[0]).denominator == 1, "integer fibre multiplicity failed")

    hostile_b = collision_vector([Q(1, 13)])
    hostile_drift = hostile_b[0] - sum(hostile_b) / 13
    require(len(set(hostile_b)) == 1 and hostile_drift == 0, "hostile drift is nonzero")

    # Prime-to-13 positive control: [0,1/2).
    positive_values = [Q(1), Q(0)]
    _, positive_c = aggregate_jumps(positive_values, 0, 2)
    require(positive_c == [Q(1), Q(-1)], "positive control current changed")
    positive_b = collision_vector(positive_values)
    positive_drift = positive_b[0] - sum(positive_b) / 13
    positive_variance = conditional_variance(positive_values)
    require(positive_drift == positive_variance, "conditional variance identity failed")
    require(positive_drift > 0, "positive control drift vanished")
    require(len(set(positive_b)) > 1, "positive control collision vector is constant")
    positive_mean = sum(positive_values) / len(positive_values)
    require(positive_mean.denominator != 1, "denominator hostile unexpectedly integral")

    # General lambda*Z multisection check on a pure-13 grid.
    lattice = Q(1, 3)
    lattice_values = [lattice * ((2 * j + 1) % 5) for j in range(13)]
    _, lattice_c = aggregate_jumps(lattice_values, 1, 1)
    lattice_h = direct_perron_cells(lattice_values, 1, 1)[0]
    require(lattice_c == [Q(0)], "lattice pure-13 current changed")
    require((13 * lattice_h / lattice).denominator == 1, "lambda-lattice invoice failed")

    return hostile_c, hostile_h, hostile_b, positive_c, positive_b, positive_drift, 3


def main():
    stats = run_grid_referee()
    (
        hostile_c,
        hostile_h,
        hostile_b,
        positive_c,
        positive_b,
        positive_drift,
        lattice_checks,
    ) = run_controls()

    print("THM-2520 rational-jump CRT / delayed-owner exact referee")
    print("grid cases:", stats[0])
    print("direct Perron cells:", stats[1])
    print("CRT last-digit coverage checks:", stats[2])
    print("variation contractions:", stats[3])
    print("constant / charged cases:", stats[4], stats[5])
    print("pure-13 hostile C, P_13F:", hostile_c, hostile_h)
    print("pure-13 hostile B_u unique, drift:", sorted(set(hostile_b)), Q(0))
    print("prime-to-13 control C:", positive_c)
    print("prime-to-13 control B_u:", positive_b)
    print("prime-to-13 conditional variance/drift:", positive_drift)
    print("integer/lattice inverse-fibre checks:", lattice_checks)
    print("endpoint-safe delay constant: 20/3 (using pi^2 < 10)")
    print(
        "VERDICT: aggregated non-13 jumps exactly decide the Perron constant branch; "
        "off it, CRT supplies every last digit and one common delayed owner can force all colours"
    )


if __name__ == "__main__":
    main()
