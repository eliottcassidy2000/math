#!/usr/bin/env python3
"""Import-free independent referee for THM-4434's coefficient box.

This implementation does not import any repository code.  It generates
oriented coefficient vectors directly from a signed Cartesian cube and finds
both error-slice and speed-polytope vertices as plane/cube-edge intersections.
That is independent of the candidate's polygon-clipping construction.
"""

from collections import defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from math import gcd
import sys


sys.stdout.reconfigure(newline="\n")
R = F(3, 14)
TARGET = F(15, 98)
EXPECTED_DIGEST = "3be81c2522a1df6a146e50634754620103f9d2d8840d17f34c9e9a4954e849f7"
EXCLUDED = {(0, 1, 1), (1, 1, 2)}
CHECKS = 0


def gate(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def gcd3(v):
    return gcd(gcd(abs(v[0]), abs(v[1])), abs(v[2]))


def dot(u, v):
    return sum(x * y for x, y in zip(u, v))


def cross(u, v):
    return (
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    )


def canonical_sign(v):
    for x in v:
        if x:
            return x > 0
    return False


def cube_edge_section(coefficients, level, low, high):
    """All vertices of {x in [low,high]^3: coefficients dot x=level}."""
    vertices = set()
    for solved in range(3):
        if coefficients[solved] == 0:
            continue
        other = [i for i in range(3) if i != solved]
        for left, right in product((low, high), repeat=2):
            point = [None, None, None]
            point[other[0]], point[other[1]] = left, right
            point[solved] = (
                level
                - coefficients[other[0]] * left
                - coefficients[other[1]] * right
            ) / coefficients[solved]
            if low <= point[solved] <= high:
                vertices.add(tuple(point))
    gate(bool(vertices), ("empty plane/cube section", coefficients, level, low, high))
    for point in vertices:
        gate(
            all(low <= x <= high for x in point) and dot(coefficients, point) == level,
            ("bad plane/cube vertex", coefficients, level, point),
        )
    return tuple(sorted(vertices))


def defects(pattern):
    radius_numerator = 3 * sum(pattern)
    unit = all(x % 3 for x in pattern)
    return tuple(
        d
        for d in range(-(radius_numerator // 14) - 1, radius_numerator // 14 + 2)
        if 14 * abs(d) < radius_numerator and ((d % 3 == 0) == unit)
    )


def oriented_vectors(maximum=18):
    """Generate sectors without sorting magnitudes or choosing an isolated sign."""
    grouped = defaultdict(list)
    for v in product(range(-maximum, maximum + 1), repeat=3):
        if not canonical_sign(v) or gcd3(v) != 1:
            continue
        if sum(x != 0 for x in v) < 2 or sum(abs(x) for x in v) % 2:
            continue
        if sum(x % 3 == 0 for x in v) > 1:
            continue
        if not (any(x > 0 for x in v) and any(x < 0 for x in v)):
            continue
        pattern = tuple(sorted(abs(x) for x in v))
        if pattern in EXCLUDED:
            continue
        grouped[pattern].append(v)
    return {pattern: tuple(sorted(vectors)) for pattern, vectors in grouped.items()}


def width_on_slice(w, v, slice_vertices):
    pivot = next(i for i, x in enumerate(v) if x)
    values = tuple(F(cross(w, e)[pivot], v[pivot]) for e in slice_vertices)
    return max(values) - min(values)


def strict_integer_levels(pattern):
    radius_numerator = 3 * sum(pattern)
    return tuple(
        d
        for d in range(-(radius_numerator // 14) - 1, radius_numerator // 14 + 2)
        if 14 * abs(d) < radius_numerator
    )


def determinant(first, second):
    return first[0] * second[1] - first[1] * second[0]


def compile_pattern(pattern, vectors):
    D = defects(pattern)
    all_levels = strict_integer_levels(pattern)
    rho = F(2, 3) if all(x % 3 for x in pattern) else F(1, 3)
    best = F(-1)
    witness = None
    for v in vectors:
        speed_vertices = cube_edge_section(v, F(0), F(0), F(1))
        all_sections = {
            d: cube_edge_section(v, F(d), -R, R)
            for d in all_levels
        }
        for w in speed_vertices:
            widths = {d: width_on_slice(w, v, all_sections[d]) for d in all_levels}
            value = rho * sum(widths[d] for d in D)

            # Independent audit of the analytic zonotope/lattice-rule chain.
            # Endpoint levels are intentionally absent: f is defined as zero
            # there to match the strict physical defect condition.
            integral = F(9, 49) * sum(w)
            f_zero = widths[0]
            sum_one = sum(widths.values())
            sum_three = sum(widths[d] for d in all_levels if d % 3 == 0)
            gate(all(widths[d] == widths[-d] for d in all_levels), ("non-even slice word", v, w))
            gate(abs(sum_one - integral) <= f_zero, ("h=1 lattice rule", v, w, sum_one, integral, f_zero))
            gate(abs(sum_three - integral / 3) <= f_zero, ("h=3 lattice rule", v, w, sum_three, integral, f_zero))
            if rho == F(2, 3):
                gate(value == F(2, 3) * sum_three, ("unit load formula", v, w))
            else:
                gate(value == F(1, 3) * (sum_one - sum_three), ("one-zero load formula", v, w))
            gate(value <= F(2, 9) * integral + F(2, 3) * f_zero, ("combined lattice error", v, w))

            maximum_index = max(range(3), key=lambda i: abs(v[i]))
            maximum_coefficient = abs(v[maximum_index])
            gate(
                f_zero <= F(6, 7) * max(w) / maximum_coefficient,
                ("central slice bound", v, w, f_zero),
            )
            generators = []
            for coordinate in range(3):
                basis = tuple(F(int(i == coordinate)) for i in range(3))
                generators.append((F(v[coordinate]), F(cross(w, basis)[maximum_index], v[maximum_index])))
            determinant_sum = sum(
                abs(determinant(generators[i], generators[j]))
                for i in range(3)
                for j in range(i + 1, 3)
            )
            gate(determinant_sum == sum(w), ("zonotope determinant sum", v, w, determinant_sum))

            if value > best:
                best = value
                witness = (v, w)
    intercept = (F(4, 3) if rho == F(2, 3) else F(1)) * len(D)
    return D, best, intercept, witness


def support_two_formula(pattern, D):
    zero, p, q = pattern
    gate(zero == 0 and p and q, ("not support two", pattern))
    return (
        2 * R * p * len(D)
        + sum(min(2 * p * R, R * (p + q) - abs(d)) for d in D)
    ) / (3 * p * q)


def main():
    grouped = oriented_vectors()
    patterns = tuple(sorted(grouped))
    gate(len(patterns) == 308, ("pattern count", len(patterns)))
    support_two = tuple(p for p in patterns if p[0] == 0)
    gate(len(support_two) == 15, ("support-two count", len(support_two)))
    gate(len(patterns) - len(support_two) == 293, "full-support count")

    semantic = sha256()
    results = {}
    maximum = F(0)
    equalities = []
    for pattern in patterns:
        D, slope, intercept, witness = compile_pattern(pattern, grouped[pattern])
        gate(slope <= TARGET, ("coefficient-box target", pattern, slope, witness))
        if pattern[0] == 0:
            gate(slope == support_two_formula(pattern, D), ("rectangle mismatch", pattern, slope))
        results[pattern] = (D, slope, intercept, witness)
        maximum = max(maximum, slope)
        if slope == TARGET:
            equalities.append(pattern)
        semantic.update(repr((pattern, D, slope, intercept)).encode())

    norm_four_vectors = tuple(
        v
        for v in product(range(-2, 3), repeat=3)
        if tuple(sorted(abs(x) for x in v)) == (1, 1, 2)
        and canonical_sign(v)
        and any(x > 0 for x in v)
        and any(x < 0 for x in v)
    )
    _, norm_four_slope, _, _ = compile_pattern((1, 1, 2), norm_four_vectors)

    gate(maximum == TARGET and equalities == [(1, 7, 8)], ("maximum/equality", maximum, equalities))
    gate(norm_four_slope == F(2, 7), ("norm-four hostile", norm_four_slope))
    gate(semantic.hexdigest() == EXPECTED_DIGEST, ("semantic digest", semantic.hexdigest()))

    threshold_gap = F(2, 11) - TARGET
    gate(threshold_gap == F(31, 1078), ("count threshold gap", threshold_gap))
    gate(F(6, 49) + F(4, 7 * 19) == F(142, 931), "M=19 slope value")
    gate(TARGET - F(142, 931) == F(1, 1862), "M=19 strict margin")
    gate(F(2, 7) / threshold_gap == F(308, 31), "relation-norm coefficient")
    gate(F(4, 3) / threshold_gap == F(4312, 93), "relation-norm intercept")
    gate(F(308, 31) * 56 + F(4312, 93) == F(56056, 93) < 603, "S<=56 cutoff")
    g58 = F(3, 16) * 58 * 58 - F(308, 31) * 58 - F(4312, 93)
    gate(g58 == F(3023, 372) > 0, ("quadratic cutoff at S=58", g58))
    gate(F(3, 8) * 58 - F(308, 31) > 0, "quadratic monotonicity from S=58")

    print("THM-4434 COEFFICIENT-BOX IMPORT-FREE CUBE-EDGE REFEREE")
    print("arithmetic=exact Fraction; repository_imports=0; optimizable_assertions=0")
    print(f"universe={len(patterns)}; full_support={len(patterns)-len(support_two)}; support_two={len(support_two)}")
    print(f"oriented_sectors={sum(len(vectors) for vectors in grouped.values())}")
    print(f"maximum={maximum}; equality_patterns={equalities}; norm_four_hostile={norm_four_slope}")
    print("actual_zero=" + ",".join(f"{p}:{results[p][1]}" for p in support_two))
    print(f"semantic_sha256={semantic.hexdigest()}")
    print("analytic_constants=M19:142/931:margin1/1862; threshold:308S/31+4312/93; S56:56056/93; g58:3023/372")
    print(f"checks={CHECKS}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
