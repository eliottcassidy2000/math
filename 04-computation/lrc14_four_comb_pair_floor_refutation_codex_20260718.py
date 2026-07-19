#!/usr/bin/env python3
"""Exact referee for THM-1191, the four-comb pair-floor refutation.

The proof-facing arithmetic is Fraction-exact.  The script independently
checks the folded and capped-trapezoid overlap formulas, proves the closed
formula on the 13-adic ladder, verifies the exact counterexample, and checks
the incidence constants in the six-wall averaging guardrail.  The finite
q-scan is telemetry only and is explicitly not used as a universal floor.
"""

from fractions import Fraction as F
from itertools import combinations, permutations
from math import gcd


def fold14(r: int) -> int:
    r %= 14
    return r * (14 - r)


def pair_overlap(a: int, b: int) -> F:
    """Haar measure of D_a intersect D_b, D_s={t: ||st||<1/14}."""
    assert a > 0 and b > 0
    g = gcd(a, b)
    a //= g
    b //= g
    if a > b:
        a, b = b, a
    return F(4 * a * b + fold14(a + b) - fold14(b - a), 196 * a * b)


def pair_overlap_trapezoid(a: int, b: int) -> F:
    """Independent capped-trapezoid evaluation of the same overlap."""
    assert a > 0 and b > 0
    g = gcd(a, b)
    width = F(a + b, 14 * a * b)
    cap = F(1, 7 * max(a, b))
    lattice_step = F(g, a * b)
    total = F(0)
    j = 0
    while j * lattice_step < width:
        piece = min(width - j * lattice_step, cap)
        total += piece if j == 0 else 2 * piece
        j += 1
    return g * total


def pair_mass(speeds: tuple[int, ...]) -> F:
    return sum((pair_overlap(a, b) for a, b in combinations(speeds, 2)), F(0))


def ladder_overlap(gap: int) -> F:
    """Closed form rho(1,13^gap)."""
    q = 13**gap
    sign = -1 if gap % 2 else 1
    return F(1, 49) * (1 + sign * F(6, q))


def strongly_connected_components(
    vertices: tuple[int, ...], edges: set[tuple[int, int]]
) -> list[tuple[int, ...]]:
    """Tiny exact SCC routine for the four-vertex fingerprint."""
    unseen = set(vertices)
    components: list[tuple[int, ...]] = []
    while unseen:
        seed = min(unseen)

        def reachable(start: int, reverse: bool = False) -> set[int]:
            found = {start}
            changed = True
            while changed:
                changed = False
                for u, v in edges:
                    if reverse:
                        u, v = v, u
                    if u in found and v not in found:
                        found.add(v)
                        changed = True
            return found

        component = reachable(seed).intersection(reachable(seed, reverse=True))
        components.append(tuple(sorted(component)))
        unseen.difference_update(component)
    return sorted(components)


def main() -> None:
    print("THM-1191 exact referee: four-comb pair-floor refutation")
    print("=" * 76)

    print("pair-overlap formula crosscheck")
    crosschecks = 0
    for a in range(1, 65):
        for b in range(a + 1, 97):
            assert pair_overlap(a, b) == pair_overlap_trapezoid(a, b)
            crosschecks += 1
    print(f"folded formula == capped trapezoid: {crosschecks} pairs")

    print("\nthirteen-adic self-similar ladder")
    for gap in range(1, 9):
        q = 13**gap
        assert q % 14 == (13 if gap % 2 else 1)
        assert pair_overlap(1, q) == ladder_overlap(gap)
    for common_scale in (1, 2, 7, 31, 169):
        speeds = tuple(common_scale * 13**i for i in range(4))
        assert pair_mass(speeds) == F(210, 2197)
    print("rho(1,13^d)=1/49*(1-6/13^d), d odd")
    print("rho(1,13^d)=1/49*(1+6/13^d), d even")
    print("common-scale checks: u=1,2,7,31,169")

    speeds = (1, 13, 169, 2197)
    rows = tuple((a, b, pair_overlap(a, b)) for a, b in combinations(speeds, 2))
    expected = (
        (1, 13, F(1, 91)),
        (1, 169, F(25, 1183)),
        (1, 2197, F(313, 15379)),
        (13, 169, F(1, 91)),
        (13, 2197, F(25, 1183)),
        (169, 2197, F(1, 91)),
    )
    assert rows == expected
    total = sum((row[2] for row in rows), F(0))
    old_candidate = F(13, 135)
    assert total == F(210, 2197)
    assert old_candidate - total == F(211, 296595) > 0
    print("pair rows:", ", ".join(f"({a},{b}):{rho}" for a, b, rho in rows))
    print(f"R4={total} < 13/135 by {old_candidate - total}")

    print("\nfinite integer-ratio telemetry (not a universal proof)")
    geometric_rows = []
    for q in range(2, 201):
        value = pair_mass((1, q, q * q, q * q * q))
        geometric_rows.append((value, q))
    geometric_rows.sort()
    assert geometric_rows[0] == (F(210, 2197), 13)
    assert sum(value < old_candidate for value, _ in geometric_rows) == 1
    print("q=2..200 consecutive geometric chains: unique old-floor refuter q=13")
    print("scope: telemetry only; no corrected universal R4 floor is claimed")

    print("\nsix-wall four-subset averaging guardrail")
    vertices = tuple(range(6))
    four_sets = tuple(combinations(vertices, 4))
    assert len(four_sets) == 15
    edge_multiplicity = {
        edge: sum(set(edge).issubset(four_set) for four_set in four_sets)
        for edge in combinations(vertices, 2)
    }
    assert set(edge_multiplicity.values()) == {6}
    route_cap = F(5, 2) * total
    target_r6 = F(13, 54)
    assert route_cap == F(525, 2197)
    assert target_r6 - route_cap == F(211, 118638) > 0
    print("sum_[|A|=4] R4(A)=6 R6; number of four-subsets=15")
    print("any uniform four-floor L gives R6 >= (5/2)L")
    print(f"counterexample forces L<=210/2197, so black-box route <= {route_cap}")
    print(f"13/54 exceeds that route ceiling by {target_r6 - route_cap}")

    # If B is the six-pair floor, a gap of length 6/(7a) contributes
    # (6B)/(7a) before endpoint debts in the post-chi inequality.
    local_route_cap = F(6, 7) * route_cap
    local_target = F(6, 7) * target_r6
    assert local_route_cap == F(450, 2197)
    assert local_target == F(13, 63)
    assert local_target - local_route_cap == F(211, 138411)
    gcd_route_cap = F(196, 13) * local_route_cap
    gcd_target = F(196, 13) * local_target
    assert gcd_route_cap == F(88200, 28561)
    assert gcd_target == F(28, 9)
    assert gcd_target - gcd_route_cap == F(5908, 257049)
    print(f"post-chi local constant: route <= {local_route_cap} < target {local_target}")
    print(f"post-chi gcd RHS constant: route <= {gcd_route_cap} < target {gcd_target}")
    print("the drift coefficient 72/13 is unchanged; only the constant term drops")

    print("\ntournament and alternate-carrier audit")
    exponent_vertices = (0, 1, 2, 3)
    edges: set[tuple[int, int]] = set()
    natural_edges = set(combinations(exponent_vertices, 2))
    for i, j in combinations(exponent_vertices, 2):
        correction = pair_overlap(13**i, 13**j) - F(1, 49)
        assert correction != 0
        # Gauge: low overlap follows increasing exponent; high overlap reverses.
        edges.add((i, j) if correction < 0 else (j, i))
    assert edges == {(0, 1), (2, 0), (0, 3), (1, 2), (3, 1), (2, 3)}
    scores = tuple(sorted(sum(u == v for u, _ in edges) for v in exponent_vertices))
    assert scores == (1, 1, 2, 2)
    directed_triangles = []
    for triple in combinations(exponent_vertices, 3):
        if all(sum((v, w) in edges for w in triple if w != v) == 1 for v in triple):
            directed_triangles.append(triple)
    assert directed_triangles == [(0, 1, 2), (1, 2, 3)]
    components = strongly_connected_components(exponent_vertices, edges)
    assert components == [(0, 1, 2, 3)]
    paths = tuple(
        p
        for p in permutations(exponent_vertices)
        if all((p[i], p[i + 1]) in edges for i in range(3))
    )
    assert len(paths) == 5
    flips = len(edges.symmetric_difference(natural_edges)) // 2
    assert flips == 2
    print("vertices challenged: 13-adic exponent levels, not runners or arcs")
    print("observable: signed correction rho(13^i,13^j)-1/49")
    print("gauge: low correction points up the exponent line; high reverses")
    print("tie path: natural exponent order (unused here: there are no ties)")
    print("fingerprint: scores 1,1,2,2; cycles 2; SCC sizes [4]; paths 5; flips 2")
    print("preserved: parity sign pattern; destroyed: 13^-gap magnitude and exact R4")
    print("faithful carrier: ordered exponent stalk + parity colour + gap weight")
    print("challenged carriers: runners, gaps, residues, valuation levels, obligations")
    print("=" * 76)
    print("done")


if __name__ == "__main__":
    main()
