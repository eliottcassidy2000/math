#!/usr/bin/env python3
"""Exact audit for THM-774, the folded-diamond two-sheet obstruction.

The proof in THM-774 is elementary and uniform.  This companion performs
three independent finite checks around its fragile interfaces:

* atom-by-atom equivalence with THM-769's eligibility/parity predicate;
* a piecewise-linear measure computation versus the closed formula;
* the complete finite exceptional table in the sharp universal bound.

All geometric arithmetic is rational.  The sampled ranges are audits of the
uniform proof, not substitutes for it.
"""

from fractions import Fraction
from math import gcd


EPSILON = Fraction(2, 13)
SHARP_CAP = Fraction(8, 117)


def circle_norm(z):
    """Distance from the rational z to the nearest integer."""
    r = z.numerator % z.denominator
    return Fraction(min(r, z.denominator - r), z.denominator)


def half_distance(n, tau):
    """Distance from n*tau to Z+1/2."""
    return Fraction(1, 2) - circle_norm(n * tau)


def nearest_integer(z):
    """Unique nearest integer, or None at a half-integer tie."""
    q, r = divmod(z.numerator, z.denominator)
    if 2 * r < z.denominator:
        return q
    if 2 * r > z.denominator:
        return q + 1
    return None


def thm769_two_sheet_predicate(x, y, tau):
    """Closed eligibility plus opposite nearest-integer parity."""
    if circle_norm(x * tau) > EPSILON:
        return False
    if circle_norm(y * tau) > EPSILON:
        return False
    nx = nearest_integer(x * tau)
    ny = nearest_integer(y * tau)
    assert nx is not None and ny is not None
    return (nx - ny) % 2 == 1


def diamond_predicate(a, b, tau):
    return circle_norm(a * tau) + circle_norm(b * tau) >= 1 - EPSILON


def two_lift_clearance(x, y, tau):
    return max(
        min(
            circle_norm(x * (tau + j) / 2),
            circle_norm(y * (tau + j) / 2),
        )
        for j in (0, 1)
    )


def tent_breakpoints(n):
    return {Fraction(k, 2 * n) for k in range(2 * n + 1)}


def danger_boundaries(n):
    points = set()
    for k in range(n):
        points.add((Fraction(k) + EPSILON) / n)
        points.add((Fraction(k + 1) - EPSILON) / n)
    return points


def diamond_boundaries(a, b):
    """All roots h_a+h_b=epsilon, found on exact affine atoms."""
    base = sorted(tent_breakpoints(a) | tent_breakpoints(b))
    roots = set()
    for left, right in zip(base, base[1:]):
        fl = half_distance(a, left) + half_distance(b, left)
        fr = half_distance(a, right) + half_distance(b, right)
        if fl == EPSILON:
            roots.add(left)
        if fr == EPSILON:
            roots.add(right)
        if (fl - EPSILON) * (fr - EPSILON) < 0:
            ratio = (EPSILON - fl) / (fr - fl)
            roots.add(left + ratio * (right - left))
    return roots


def audit_identity(a, b):
    x, y = a + b, a - b
    points = (
        {Fraction(0), Fraction(1)}
        | tent_breakpoints(a)
        | tent_breakpoints(b)
        | tent_breakpoints(x)
        | tent_breakpoints(y)
        | danger_boundaries(x)
        | danger_boundaries(y)
        | diamond_boundaries(a, b)
    )
    ordered = sorted(points)
    probes = set(ordered[:-1])
    probes.update((left + right) / 2 for left, right in zip(ordered, ordered[1:]))
    for tau in probes:
        lhs = thm769_two_sheet_predicate(x, y, tau)
        rhs = diamond_predicate(a, b, tau)
        assert lhs == rhs, (a, b, tau, lhs, rhs)
        q_value = circle_norm(a * tau) + circle_norm(b * tau)
        assert two_lift_clearance(x, y, tau) == (1 - q_value) / 2
    return len(probes)


def piecewise_diamond_measure(a, b):
    """Independent integration of h_a+h_b<=epsilon on affine atoms."""
    points = sorted(tent_breakpoints(a) | tent_breakpoints(b))
    total = Fraction()
    for left, right in zip(points, points[1:]):
        fl = half_distance(a, left) + half_distance(b, left)
        fr = half_distance(a, right) + half_distance(b, right)
        width = right - left
        if fl <= EPSILON and fr <= EPSILON:
            total += width
        elif fl <= EPSILON < fr:
            total += width * (EPSILON - fl) / (fr - fl)
        elif fr <= EPSILON < fl:
            total += width * (EPSILON - fr) / (fl - fr)
    return total


def qualifying_tooth_count(a):
    return (4 * a + 13) // 26


def closed_formula(a, b):
    """Exact measure for reduced coprime, opposite-parity a>b."""
    assert a > b > 0 and gcd(a, b) == 1 and (a - b) % 2 == 1
    r, s = a - b, a + b
    n = qualifying_tooth_count(a)
    tooth_sum = sum(
        (
            min(
                Fraction(4 * r, 13),
                Fraction(4 * a, 13) - (2 * j - 1),
            )
            for j in range(1, n + 1)
        ),
        Fraction(),
    )
    return Fraction(2, r * s) * tooth_sum


def safe_measure(core):
    """Exact measure of {tau:min_u ||u*tau||>1/13}."""
    delta = Fraction(1, 13)
    points = {Fraction(0), Fraction(1)}
    for u in core:
        for k in range(u):
            points.add((Fraction(k) + delta) / u)
            points.add((Fraction(k + 1) - delta) / u)
    ordered = sorted(points)
    total = Fraction()
    for left, right in zip(ordered, ordered[1:]):
        tau = (left + right) / 2
        if min(circle_norm(u * tau) for u in core) > delta:
            total += right - left
    return total


def total_order_tournament(vertices, key):
    """Tournament from a scalar gauge and an explicit lexicographic tie path."""
    order = sorted(vertices, key=lambda v: (key(v), v))
    rank = {v: i for i, v in enumerate(order)}
    edges = set()
    scores = {v: 0 for v in vertices}
    for i, u in enumerate(vertices):
        for v in vertices[i + 1 :]:
            edge = (u, v) if rank[u] < rank[v] else (v, u)
            edges.add(edge)
            scores[edge[0]] += 1
    return {
        "order": tuple(order),
        "edges": edges,
        "score_histogram": tuple(sorted(scores.values())),
        "directed_cycles": 0,
        "scc_sizes": (1,) * len(vertices),
        "hamiltonian_paths": 1,
    }


def main():
    print("THM-774 folded-diamond two-sheet obstruction certificate")
    print("closed exception danger / strict core-loose convention")
    print()

    identity_cases = 0
    identity_probes = 0
    for a in range(2, 81):
        for b in range(1, a):
            if (a - b) % 2 != 1:
                continue
            identity_cases += 1
            identity_probes += audit_identity(a, b)
    print("predicate identity audit")
    print("  opposite-parity (a,b) cases through a=80:", identity_cases)
    print("  exact atom/end-point probes:", identity_probes)
    print("  failures: 0")
    print("  clearance-fold B=(1-Q)/2 failures: 0")
    print()

    formula_cases = 0
    for a in range(2, 161):
        for b in range(1, a):
            if gcd(a, b) != 1 or (a - b) % 2 != 1:
                continue
            formula_cases += 1
            direct = piecewise_diamond_measure(a, b)
            formula = closed_formula(a, b)
            assert direct == formula, (a, b, direct, formula)
            assert formula <= SHARP_CAP, (a, b, formula)
    print("measure-formula audit")
    print("  reduced coprime cases through a=160:", formula_cases)
    print("  piecewise-linear/formula mismatches: 0")
    print("  cap violations: 0")
    print()

    # The proof uses two elementary bounds.  If s=a+b>=9n, the first gives
    # 8n/(13s)<=8/117.  Otherwise, the second works when 4a>=23n.  Since
    # n<=2a/13+1/2, the latter is automatic for a>=25; enumerate a<25.
    exceptional_a = tuple(
        (a, qualifying_tooth_count(a))
        for a in range(1, 25)
        if 4 * a < 23 * qualifying_tooth_count(a)
    )
    assert exceptional_a == ((4, 1), (5, 1), (10, 2), (11, 2), (17, 3))
    expected_maxima = {
        4: (Fraction(6, 91), 3),
        5: (Fraction(8, 117), 4),
        10: (Fraction(8, 169), 3),
        11: (Fraction(16, 273), 10),
        17: (Fraction(2, 39), 16),
    }
    print("sharp-cap proof audit")
    print("  automatic for a>=25 from n<=2a/13+1/2")
    print("  finite exceptional (a,n):", exceptional_a)
    exceptional_vertices = []
    equality_pairs = []
    for a, n in exceptional_a:
        rows = []
        for b in range(1, a):
            if gcd(a, b) != 1 or (a - b) % 2 != 1:
                continue
            measure = closed_formula(a, b)
            rows.append((measure, b))
            exceptional_vertices.append((a, b))
            if measure == SHARP_CAP:
                equality_pairs.append((a, b))
        maximum = max(rows)
        assert maximum == expected_maxima[a]
        print(f"  a={a:2d}, n={n}: max={maximum[0]} at b={maximum[1]}")
    assert equality_pairs == [(5, 4)]
    print("  sharp reduced equality pair: (a,b)=(5,4), hence x:y=9:1")
    print()

    print("guardrail: measure is a filter, not a closure")
    deletion_measures = []
    for dropped in range(1, 12):
        core = tuple(u for u in range(1, 12) if u != dropped)
        measure = safe_measure(core)
        deletion_measures.append((dropped, measure))
        relation = "below" if measure < SHARP_CAP else "above"
        print(f"  U={{1,...,11}}\\{{{dropped:2d}}}: |G_U|={measure} ({relation} cap)")
    assert dict(deletion_measures)[2] == Fraction(569, 10010)
    assert dict(deletion_measures)[4] == Fraction(499, 10010)
    assert dict(deletion_measures)[10] == Fraction(746, 15015)
    print("  conclusion: several loose cores survive the scalar cap")
    print()

    # Tournament Analysis.  We deliberately use finite proof obligations as
    # vertices.  Runner vertices lose the half-frequency intersection; safe
    # components lose the tooth label.  Even this obligation tournament is
    # telemetry: its scalar gauges discard component locations and cannot test
    # containment G_U subset H_{a,b}.
    metric = total_order_tournament(
        exceptional_vertices, lambda v: closed_formula(*v)
    )
    coarse = total_order_tournament(
        exceptional_vertices,
        lambda v: Fraction(8 * qualifying_tooth_count(v[0]), 13 * (v[0] + v[1])),
    )
    flips = len(metric["edges"] ^ coarse["edges"]) // 2
    print("Tournament Analysis (vertices = finite exceptional proof obligations)")
    print("  pairwise observable: relative obstruction size")
    print("  gauge A: exact folded-diamond measure")
    print("  gauge B: coarse first-bound value 8n/(13s)")
    print("  tie Hamiltonian path: lexicographic (a,b)")
    print("  score histogram:", metric["score_histogram"])
    print("  directed cycles:", metric["directed_cycles"], coarse["directed_cycles"])
    print("  SCC sizes:", metric["scc_sizes"])
    print("  Hamiltonian-path counts:", metric["hamiltonian_paths"], coarse["hamiltonian_paths"])
    print("  edge flips between gauges:", flips)
    print("  challenged vertex sets: runners, half-frequencies, safe components")
    print("  preserved predicate: pairwise size ranking only")
    print("  destroyed data: tooth locations and the containment G_U subset H_{a,b}")
    print()
    print("FINAL: PASS")


if __name__ == "__main__":
    main()
