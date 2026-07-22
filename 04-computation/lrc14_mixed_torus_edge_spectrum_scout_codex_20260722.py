#!/usr/bin/env python3
"""Exact referee/scout for THM-2099's mixed-torus positive-edge spectrum.

It computes the
Haar mass of

    |u.X|<1/14, |v.X|<1/14, |g.X|>1/7

exactly from the primitive relation among three rank-two characters.  It then
tests the proposed strict rank-eight maximum-spanning-tree gate on structured
and seeded-random coefficient sets.  The exact dyadic row refutes that gate
and is the proof-facing output; the bounded relation/random scans are labelled
diagnostics.  Runtime checks do not use ``assert``.

Tournament Analysis: the pair observable is the exact restricted overlap.
Orienting an edge toward the larger endpoint after sorting by incident weight
is only a deterministic tie scheduler; it discards threshold colors and the
graphic-matroid optimum.  The script therefore reports the weight histogram,
low-edge graph degrees/components, MST, and tie count, not tournament scores
as if they were proof invariants.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from math import gcd
from random import Random


GUARD = (1, 0)
BUDGET = F(5, 49)
GENERIC = F(5, 343)
NONFACTOR_FLOOR = F(1, 98)
FIRST_POSITIVE = F(1, 392)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def fmt(value: F) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def det(a: tuple[int, int], b: tuple[int, int]) -> int:
    return a[0] * b[1] - a[1] * b[0]


def positive_part_square(value: int) -> int:
    return max(value, 0) ** 2


def convolution_cdf(s: int, b: int, c: int) -> F:
    """P(b U+c V <= s), for independent U,V uniform on [-1,1]."""
    b, c = abs(b), abs(c)
    require(b > 0 and c > 0, "convolution coefficients must be nonzero")
    t = s + b + c
    numerator = (
        positive_part_square(t)
        - positive_part_square(t - 2 * b)
        - positive_part_square(t - 2 * c)
        + positive_part_square(t - 2 * b - 2 * c)
    )
    return F(numerator, 8 * b * c)


def convolution_interval(left: int, right: int, b: int, c: int) -> F:
    """P(left <= b U+c V <= right); endpoints have measure zero."""
    if right <= left:
        return F(0)
    return convolution_cdf(right, b, c) - convolution_cdf(left, b, c)


def independent_outside_fraction(a: int, b: int, c: int) -> F:
    """Outside-guard fraction conditional on the two terminal danger bands.

    The primitive relation is ``a*g=b*u+c*v`` with ``a>0``.  The roots of
    this equation on each fiber are averaged exactly, and the convolution of
    the two terminal coordinates is integrated by its piecewise-quadratic CDF.
    """
    require(a > 0 and b != 0 and c != 0, "invalid primitive relation")
    b, c = abs(b), abs(c)
    support = b + c
    period = 14 * a
    radius = 2 * a
    inside_total = F(0)
    margin = (support + period + radius) // period + 2
    for k in range(a):
        inside = F(0)
        for m in range(-margin, margin + 2):
            center = period * m - 14 * k
            left = max(-support, center - radius)
            right = min(support, center + radius)
            inside += convolution_interval(left, right, b, c)
        require(F(0) <= inside <= 1, "fiber inside mass escaped [0,1]")
        inside_total += inside
    return 1 - inside_total / a


def fold(residue: int, modulus: int) -> int:
    residue %= modulus
    return residue * (modulus - residue)


def rho_fold(a: int, b: int) -> F:
    """Haar measure of the two radius-1/14 integer combs."""
    a, b = abs(a), abs(b)
    require(a > 0 and b > 0, "dependent speeds must be nonzero")
    if a > b:
        a, b = b, a
    common = gcd(a, b)
    modulus = 14 * common
    numerator = 4 * a * b + fold(a + b, modulus) - fold(b - a, modulus)
    return F(numerator, 196 * a * b)


def primitive_multiple(vector: tuple[int, int]) -> tuple[tuple[int, int], int]:
    common = gcd(abs(vector[0]), abs(vector[1]))
    require(common > 0, "zero character is not allowed")
    primitive = (vector[0] // common, vector[1] // common)
    if primitive[0] < 0 or (primitive[0] == 0 and primitive[1] < 0):
        primitive = (-primitive[0], -primitive[1])
        common = -common
    return primitive, common


def edge_data(
    u: tuple[int, int], v: tuple[int, int], g: tuple[int, int] = GUARD
) -> tuple[F, tuple[int, int, int] | tuple[str, int, int]]:
    """Return exact restricted overlap and its primitive relation label."""
    duv = det(u, v)
    if duv == 0:
        primitive, first = primitive_multiple(u)
        require(det(primitive, g) != 0, "terminal must be transverse to guard")
        if primitive[0] * v[1] == primitive[1] * v[0]:
            second = v[0] // primitive[0] if primitive[0] else v[1] // primitive[1]
        else:
            raise RuntimeError("dependent-character decomposition failed")
        return F(5, 7) * rho_fold(first, second), ("dep", abs(first), abs(second))

    dvg = det(v, g)
    dgu = det(g, u)
    common = gcd(abs(duv), gcd(abs(dvg), abs(dgu)))
    a = duv // common
    b = -dvg // common
    c = -dgu // common
    if a < 0:
        a, b, c = -a, -b, -c
    require(a * g[0] == b * u[0] + c * v[0], "relation x-coordinate failed")
    require(a * g[1] == b * u[1] + c * v[1], "relation y-coordinate failed")
    fraction = independent_outside_fraction(a, b, c)
    return fraction / 49, (a, abs(b), abs(c))


def maximum_spanning_tree(weights: list[list[F]]) -> tuple[F, tuple[tuple[int, int, F], ...]]:
    n = len(weights)
    edges = sorted(
        ((weights[i][j], i, j) for i in range(n) for j in range(i + 1, n)),
        reverse=True,
    )
    parent = list(range(n))

    def root(vertex: int) -> int:
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    chosen: list[tuple[int, int, F]] = []
    for weight, i, j in edges:
        ri, rj = root(i), root(j)
        if ri != rj:
            parent[ri] = rj
            chosen.append((i, j, weight))
            if len(chosen) == n - 1:
                break
    require(len(chosen) == n - 1, "complete weighted graph was disconnected")
    return sum((edge[2] for edge in chosen), F(0)), tuple(chosen)


def fingerprint(vectors: tuple[tuple[int, int], ...]) -> dict[str, object]:
    n = len(vectors)
    weights = [[F(0) for _ in range(n)] for _ in range(n)]
    relations: Counter[object] = Counter()
    for i in range(n):
        for j in range(i + 1, n):
            weight, relation = edge_data(vectors[i], vectors[j])
            weights[i][j] = weights[j][i] = weight
            relations[relation] += 1
    tree_weight, tree = maximum_spanning_tree(weights)
    histogram = Counter(weights[i][j] for i in range(n) for j in range(i + 1, n))
    low_degrees = tuple(
        sum(weights[i][j] < NONFACTOR_FLOOR for j in range(n) if i != j)
        for i in range(n)
    )
    tie_count = sum(count * (count - 1) // 2 for count in histogram.values())
    return {
        "tau": tree_weight,
        "tree": tree,
        "histogram": histogram,
        "relations": relations,
        "low_degrees": low_degrees,
        "tie_count": tie_count,
    }


def print_case(name: str, vectors: tuple[tuple[int, int], ...]) -> F:
    data = fingerprint(vectors)
    tau = data["tau"]
    print(f"{name}: vectors={vectors}")
    print(
        "  tau=", fmt(tau), " delta=", fmt(tau - BUDGET),
        " low-degrees=", data["low_degrees"], " tie-pairs=", data["tie_count"],
        sep="",
    )
    print(
        "  weights=",
        ", ".join(f"{fmt(weight)}:{count}" for weight, count in sorted(data["histogram"].items())),
        sep="",
    )
    print(
        "  tree=",
        tuple((i, j, fmt(weight)) for i, j, weight in data["tree"]),
        sep="",
    )
    return tau


def relation_spectrum() -> None:
    minimum = (F(1), None)
    nonfactor_minimum = (F(1), None)
    below_half: list[tuple[F, tuple[int, int, int]]] = []
    ranges = [(1, 1, 40), (2, 8, 16)]
    for a_lo, a_hi, bc_hi in ranges:
        for a in range(a_lo, a_hi + 1):
            for b in range(1, bc_hi + 1):
                for c in range(1, bc_hi + 1):
                    if gcd(a, gcd(b, c)) != 1:
                        continue
                    fraction = independent_outside_fraction(a, b, c)
                    if 0 < fraction < minimum[0]:
                        minimum = (fraction, (a, b, c))
                    if a > 1 and fraction < nonfactor_minimum[0]:
                        nonfactor_minimum = (fraction, (a, b, c))
                    if a == 1 and fraction < F(1, 2):
                        below_half.append((fraction, (a, b, c)))
    require(minimum == (F(1, 8), (1, 1, 2)), "first positive layer changed")
    require(nonfactor_minimum[0] >= F(1, 2), "nonfactor floor failed on scan")
    print("diagnostic relation scan: a=1,b,c<=40 and 2<=a<=8,b,c<=16")
    print("  minimum positive outside fraction:", fmt(minimum[0]), minimum[1])
    print("  minimum nonfactor outside fraction:", fmt(nonfactor_minimum[0]), nonfactor_minimum[1])
    print("  distinct integral-factor fractions below 1/2:")
    spectrum = Counter(fraction for fraction, _ in below_half)
    print("   ", ", ".join(f"{fmt(value)}:{count}" for value, count in sorted(spectrum.items())))


def structured_cases() -> None:
    cases = {
        "dyadic": tuple((1, 2**k) for k in range(8)),
        "triadic": tuple((1, 3**k) for k in range(8)),
        "consecutive": tuple((1, k) for k in range(1, 9)),
        "fibonacci": ((1, 1), (1, 2), (1, 3), (1, 5), (1, 8), (1, 13), (1, 21), (1, 34)),
        "alternating-x": tuple(((-1) ** k, k + 1) for k in range(8)),
    }
    values = {name: print_case(name, vectors) for name, vectors in cases.items()}
    require(values["dyadic"] == F(70541889, 691400192), "dyadic MST changed")
    require(values["dyadic"] < BUDGET, "dyadic row no longer refutes strict tree gate")
    require(values["triadic"] > BUDGET, "triadic control should cross the tree gate")


def dyadic_exact_gate() -> None:
    vectors = tuple((1, 2**k) for k in range(8))
    expected = {
        1: F(1, 392),
        2: F(2, 147),
        3: F(5, 343),
        4: F(227, 15680),
        5: F(177, 12152),
        6: F(5, 343),
        7: F(46441, 3186176),
    }
    for gap, weight in expected.items():
        for start in range(8 - gap):
            actual, relation = edge_data(vectors[start], vectors[start + gap])
            require(actual == weight, f"dyadic gap-{gap} mass changed")
            require(
                relation == (2**gap - 1, 2**gap, 1),
                f"dyadic gap-{gap} relation changed",
            )
    require(expected[3] == expected[6] == GENERIC, "generic dyadic layers changed")
    require(expected[7] < expected[3], "gap-7 ordering changed")
    require(expected[5] < expected[7], "gap-5 ordering changed")
    require(expected[4] < expected[5], "gap-4 ordering changed")
    expected_tree = 5 * GENERIC + expected[7] + expected[5]
    require(expected_tree == F(70541889, 691400192), "Kruskal sum changed")
    print("exact dyadic gate:")
    print(
        "  gap weights:",
        ", ".join(f"r={gap}:{fmt(weight)}" for gap, weight in expected.items()),
    )
    print("  five generic edges form three mod-3 components; Kruskal adds gaps 7 and 5")
    print("  tau=", fmt(expected_tree), " < 5/49 by ", fmt(BUDGET - expected_tree), sep="")
    x, y = F(1, 3), F(0)
    require(abs(x) > F(1, 7), "dyadic witness guard check failed")
    require(
        all(abs((x + 2**k * y) % 1) in (F(1, 3), F(2, 3)) for k in range(8)),
        "dyadic witness terminal check failed",
    )
    print("  safe point X=(1/3,0): guard and all eight terminals equal 1/3 mod 1")


def dyadic_residue_law() -> None:
    ys = tuple(2**k for k in range(8))
    constants = {
        0: F(0),
        1: -F(9151, 691400192),
        2: F(213, 24304),
        3: -F(1031, 41818560),
        4: F(1, 49),
        5: F(3, 49),
        6: F(33, 392),
    }
    for a in range(1, 36):
        vectors = tuple((a, y) for y in ys)
        delta = fingerprint(vectors)["tau"] - BUDGET
        require(a * delta == constants[a % 7], f"dyadic residue law failed at a={a}")
    print("dyadic affine-pencil residue law: a*(tau-5/49), indexed by a mod 7")
    print("  " + ", ".join(f"{residue}:{fmt(value)}" for residue, value in constants.items()))
    print("  negative exactly for residues 1,3; zero exactly for residue 0")


def random_scout(samples: int = 250) -> None:
    rng = Random(2099)
    best: tuple[F, tuple[tuple[int, int], ...] | None] = (F(10), None)
    closed = 0
    for _ in range(samples):
        ys = sorted(rng.sample(range(1, 25), 8))
        vectors = tuple((rng.randint(-8, 8), y) for y in ys)
        data = fingerprint(vectors)
        tau = data["tau"]
        if tau < best[0]:
            best = (tau, vectors)
        if tau <= BUDGET:
            closed += 1
    require(best[1] is not None, "random scout produced no rows")
    print(f"seeded random scout: samples={samples}, tau<=5/49 rows={closed}")
    print_case("random-minimum", best[1])


def main() -> None:
    require(independent_outside_fraction(1, 1, 1) == 0, "zero layer failed")
    require(independent_outside_fraction(1, 2, 1) == F(1, 8), "first layer failed")
    require(independent_outside_fraction(1, 2, 2) == F(1, 4), "balanced layer failed")
    require(independent_outside_fraction(1, 3, 1) == F(1, 3), "3+1 layer failed")
    require(edge_data((0, 1), (1, -2))[0] == FIRST_POSITIVE, "vector relation check failed")
    print("THM-2099 MIXED-TORUS POSITIVE-EDGE SCOUT")
    print("budget=", fmt(BUDGET), " generic-edge=", fmt(GENERIC), sep="")
    relation_spectrum()
    dyadic_exact_gate()
    dyadic_residue_law()
    structured_cases()
    random_scout()
    print("PASS -- exact formula/dyadic no-go verified; bounded scans remain diagnostics")


if __name__ == "__main__":
    main()
