#!/usr/bin/env python3
"""Exact referee for THM-2114's quarter-torsion/toothpick lemmas.

The proof-facing checks use only integer and ``Fraction`` arithmetic.  They
verify the moving-root count, the exact hostile guard-pencil spectrum and its
maximum spanning tree, and the quarter-torsion parity witness.  A bounded
primitive-lattice scan checks (but does not prove) the all-height C4-free
toothpick argument supplied in the theorem.

Tournament Analysis: the pair observable is the restricted collision mass.
Orienting its complete graph by weight is only a scheduler; the proof carrier
is the nested undirected threshold graph plus its maximum graphic-matroid
bases and finite-ring kernel sidecar.  Consequently this referee reports
exact weight layers and maximum-tree/equality data, not label-dependent
tournament scores.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from itertools import product
from math import gcd


BUDGET = F(5, 49)
PAIR_BASELINE = F(5, 343)
K = (0, 7, 14, 15, 16, 17, 24, 31)
EQUALITY_VECTORS = (
    (0, 7),
    (-4, 10),
    (11, 2),
    (-18, 10),
    (1, 8),
    (-19, 9),
    (-3, 18),
    (9, 2),
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def fmt(value: F) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def circle_norm(value: F) -> F:
    residue = value % 1
    return min(residue, 1 - residue)


def moving_root_count(m: int, s: F) -> int:
    return sum(circle_norm((s + k) / m) < F(1, 7) for k in range(m))


def audit_moving_root_law() -> int:
    checks = 0
    samples = tuple(F(j, 211) for j in range(211))
    for m in range(1, 101):
        q, r = divmod(2 * m, 7)
        center = F(q, 2)
        radius = F(r, 14)
        for s in samples:
            if circle_norm(s - center) == radius:
                continue
            expected = q + int(r > 0 and circle_norm(s - center) < radius)
            require(moving_root_count(m, s) == expected, "moving-root law failed")
            checks += 1

        ceiling = (2 * m + 6) // 7
        if m >= 2:
            floor = F(m - ceiling, 49 * m)
            require(floor >= F(1, 98), "finite-monodromy floor failed")
    return checks


def pencil_fraction(distance: int) -> F:
    """Conditional guard-safe fraction for f_k=(k,1), f_l=(l,1)."""
    require(distance > 0, "pencil distance must be positive")
    q, r = divmod(distance, 7)
    if r == 0:
        return F(5, 7)
    return F(5 * q + r - 1, distance)


def positive_part_square(value: int) -> int:
    return max(value, 0) ** 2


def convolution_cdf(s: int, a: int, b: int) -> F:
    """P(a U+b V <= s) for independent U,V uniform on [-1,1]."""
    a, b = abs(a), abs(b)
    t = s + a + b
    numerator = (
        positive_part_square(t)
        - positive_part_square(t - 2 * a)
        - positive_part_square(t - 2 * b)
        + positive_part_square(t - 2 * a - 2 * b)
    )
    return F(numerator, 8 * a * b)


def convolution_interval(left: int, right: int, a: int, b: int) -> F:
    if right <= left:
        return F(0)
    return convolution_cdf(right, a, b) - convolution_cdf(left, a, b)


def relation_outside_fraction(m: int, a: int, b: int) -> F:
    """Conditional outside fraction for m*g=a*f+b*f'."""
    require(m > 0 and a != 0 and b != 0, "invalid primitive relation")
    a, b = abs(a), abs(b)
    support = a + b
    period = 14 * m
    radius = 2 * m
    inside_total = F(0)
    margin = (support + period + radius) // period + 2
    for k in range(m):
        inside = F(0)
        for shift in range(-margin, margin + 2):
            center = period * shift - 14 * k
            left = max(-support, center - radius)
            right = min(support, center + radius)
            inside += convolution_interval(left, right, a, b)
        require(F(0) <= inside <= 1, "fiber probability escaped [0,1]")
        inside_total += inside
    return 1 - inside_total / m


def determinant(a: tuple[int, int], b: tuple[int, int]) -> int:
    return a[0] * b[1] - a[1] * b[0]


def vector_edge_fraction(
    first: tuple[int, int], second: tuple[int, int]
) -> tuple[F, tuple[int, int, int]]:
    """Return 49*w and the primitive relation against g=(1,0)."""
    det_pair = determinant(first, second)
    require(det_pair != 0, "equality example unexpectedly has a dependent pair")
    det_second_guard = -second[1]
    det_guard_first = first[1]
    common = gcd(abs(det_pair), gcd(abs(det_second_guard), abs(det_guard_first)))
    m = det_pair // common
    a = -det_second_guard // common
    b = -det_guard_first // common
    if m < 0:
        m, a, b = -m, -a, -b
    require(
        (m, 0) == (a * first[0] + b * second[0], a * first[1] + b * second[1]),
        "primitive relation failed",
    )
    return relation_outside_fraction(m, a, b), (m, a, b)


def maximum_spanning_tree(
    labels: tuple[int, ...],
) -> tuple[F, tuple[tuple[int, int, F], ...]]:
    edges = sorted(
        (
            (pencil_fraction(abs(a - b)), a, b)
            for i, a in enumerate(labels)
            for b in labels[i + 1 :]
        ),
        reverse=True,
    )
    parent = {label: label for label in labels}

    def root(label: int) -> int:
        while parent[label] != label:
            parent[label] = parent[parent[label]]
            label = parent[label]
        return label

    chosen: list[tuple[int, int, F]] = []
    for weight, a, b in edges:
        ra, rb = root(a), root(b)
        if ra != rb:
            parent[ra] = rb
            chosen.append((a, b, weight))
            if len(chosen) == len(labels) - 1:
                break
    require(len(chosen) == len(labels) - 1, "Kruskal did not span")
    return sum((edge[2] for edge in chosen), F(0)), tuple(chosen)


def audit_hostile_pencil() -> tuple[Counter[F], F]:
    histogram = Counter(
        pencil_fraction(abs(a - b))
        for i, a in enumerate(K)
        for b in K[i + 1 :]
    )
    expected = Counter(
        {
            F(5, 7): 6,
            F(22, 31): 1,
            F(17, 24): 2,
            F(12, 17): 3,
            F(7, 10): 2,
            F(11, 16): 2,
            F(2, 3): 5,
            F(5, 8): 2,
            F(1, 2): 2,
            F(0): 3,
        }
    )
    require(histogram == expected, "hostile-pencil edge histogram changed")
    normalized_tau, tree = maximum_spanning_tree(K)
    require(normalized_tau == F(8579, 1736), "hostile-pencil tree changed")
    require(
        Counter(weight for _, _, weight in tree)
        == Counter({F(5, 7): 4, F(22, 31): 1, F(11, 16): 2}),
        "hostile-pencil Kruskal layers changed",
    )
    tau = normalized_tau / 49
    require(tau == BUDGET - F(101, 85064), "hostile-pencil deficit changed")

    # Direction d=(1,2): g.d=1 and f_k.d=k+2.
    speeds = tuple(k + 2 for k in K)
    require(speeds == (2, 9, 16, 17, 18, 19, 26, 33), "specialization changed")
    require(len(set(speeds)) == 8 and min(speeds) > 1, "specialization invalid")
    return histogram, tau


def audit_quarter_torsion() -> int:
    """Audit g=c p and f_i=b_i p+a_i f at (p,f)=(1/2,1/4)."""
    checks = 0
    p_value, f_value = F(1, 2), F(1, 4)
    for content in range(1, 32, 2):
        require(circle_norm(content * p_value) == F(1, 2), "guard parity failed")
        for transverse in range(-31, 32, 2):
            for horizontal in range(-20, 21):
                value = horizontal * p_value + transverse * f_value
                require(circle_norm(value) == F(1, 4), "quarter witness failed")
                checks += 1
    return checks


def affine_rank_one(points: tuple[tuple[int, int], ...]) -> bool:
    differences = tuple((x - points[0][0], y - points[0][1]) for x, y in points[1:])
    pivot = next((v for v in differences if v != (0, 0)), None)
    return pivot is None or all(determinant(pivot, v) == 0 for v in differences)


def audit_equality_counterexample() -> tuple[Counter[F], int]:
    weights: dict[tuple[int, int], F] = {}
    for i in range(8):
        for j in range(i + 1, 8):
            weights[i, j] = vector_edge_fraction(EQUALITY_VECTORS[i], EQUALITY_VECTORS[j])[0]

    special = {
        (1, 2): F(42, 59),
        (2, 3): F(52, 73),
        (2, 4): F(61, 86),
        (2, 5): F(2347, 3288),
        (2, 6): F(655, 918),
        (2, 7): F(1, 2),
    }
    for edge, weight in weights.items():
        expected = special.get(edge, F(5, 7))
        require(weight == expected, f"equality example edge {edge} changed")
    histogram = Counter(weights.values())
    require(histogram == Counter({F(5, 7): 22, **{v: 1 for v in special.values()}}), "histogram changed")

    # The equality graph is K_7 on S plus the pendant equality edge 0--2.
    s = {0, 1, 3, 4, 5, 6, 7}
    equality_edges = {edge for edge, weight in weights.items() if weight == F(5, 7)}
    expected_edges = {
        tuple(sorted((i, j))) for i in s for j in s if i < j
    } | {(0, 2)}
    require(equality_edges == expected_edges, "equality graph changed")

    # Fix the first sign: global sign does not change affine rank.
    affine_gauges = 0
    for tail_signs in product((-1, 1), repeat=7):
        signs = (1,) + tail_signs
        signed = tuple(
            (sign * vector[0], sign * vector[1])
            for sign, vector in zip(signs, EQUALITY_VECTORS)
        )
        affine_gauges += int(affine_rank_one(signed))
    require(affine_gauges == 0, "THM-2103 counterexample became a signed pencil")

    direction = (1, 3)
    speeds = tuple(x * direction[0] + y * direction[1] for x, y in EQUALITY_VECTORS)
    require(speeds == (21, 26, 17, 12, 25, 8, 51, 15), "positive specialization changed")
    require(len(set(speeds)) == 8 and min(speeds) > 1, "positive specialization invalid")

    # Missing projective kernel direction (1,-1) modulo five.
    needle = (1, -1)
    characters = ((1, 0),) + EQUALITY_VECTORS
    require(all((x * needle[0] + y * needle[1]) % 5 for x, y in characters), "mod-5 needle blocked")
    torus_point = (F(1, 5), -F(1, 5))
    distances = tuple(circle_norm(x * torus_point[0] + y * torus_point[1]) for x, y in characters)
    require(distances[0] > F(1, 7), "mod-5 guard witness failed")
    require(all(value > F(1, 14) for value in distances[1:]), "mod-5 terminal witness failed")
    return histogram, affine_gauges


def low_neighbor(vertex: tuple[int, int], other: tuple[int, int]) -> bool:
    x, d = vertex
    y, e = other
    if d == e:
        return abs(x - y) == 1
    if e == 2 * d:
        return abs(y - 2 * x) == 1
    if d == 2 * e:
        return abs(x - 2 * y) == 1
    return False


def audit_toothpick_window() -> tuple[int, int, int]:
    vertices = tuple(
        (x, d)
        for d in range(1, 17)
        for x in range(-48, 49)
        if gcd(abs(x), d) == 1
    )
    adjacency = {
        vertex: {other for other in vertices if other != vertex and low_neighbor(vertex, other)}
        for vertex in vertices
    }
    max_degree = max(map(len, adjacency.values()))
    max_codegree = max(
        len(adjacency[a] & adjacency[b])
        for i, a in enumerate(vertices)
        for b in vertices[i + 1 :]
    )
    require(max_degree <= 6, "toothpick degree exceeded six")
    require(max_codegree <= 1, "bounded toothpick window contains a C4")
    return len(vertices), max_degree, max_codegree


def main() -> None:
    root_checks = audit_moving_root_law()
    histogram, tau = audit_hostile_pencil()
    quarter_checks = audit_quarter_torsion()
    equality_histogram, affine_gauges = audit_equality_counterexample()
    vertex_count, max_degree, max_codegree = audit_toothpick_window()

    print("THM-2114 QUARTER-TORSION / TOOTHPICK REFEREE")
    print(f"moving-root law: {root_checks} exact samples; m=1..100: PASS")
    print("hostile guard pencil K=", K, sep="")
    print(
        "  normalized layers=",
        ", ".join(f"{fmt(weight)}:{count}" for weight, count in sorted(histogram.items(), reverse=True)),
        sep="",
    )
    print("  normalized tau=8579/1736; actual tau=", fmt(tau), sep="")
    print("  cap deficit=101/85064; direction speeds=(2,9,16,17,18,19,26,33)")
    print(f"quarter-torsion odd-parity checks: {quarter_checks}: PASS")
    print(
        "THM-2103 refuter: normalized weights=",
        ", ".join(
            f"{fmt(weight)}:{count}"
            for weight, count in sorted(equality_histogram.items(), reverse=True)
        ),
        sep="",
    )
    print(f"  tau=5/49 exactly; signed affine-rank-one gauges={affine_gauges}/128")
    print("  mod-5 Kakeya needle X=(1,-1)/5: strict mixed escape")
    print(
        f"diagnostic primitive toothpick window: vertices={vertex_count}, "
        f"max-degree={max_degree}, max-codegree={max_codegree}"
    )
    print("proof carrier: threshold/max-basis matroid + finite-ring kernels; tournament is scheduler only")
    print("PASS -- exact finite checks agree with THM-2114")


if __name__ == "__main__":
    main()
