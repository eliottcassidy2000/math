#!/usr/bin/env python3
"""
Unit-distance n=22 Mathieu/M21 residue scout.

Prompt: consider the 22-point unit-distance frontier and the fact that M24
produces M23 and M22 easily, but M21 is isomorphic to L3(4).

This script turns that group-theoretic observation into a finite carrier.

Facts used:
  * M22 is the two-point stabilizer in M24 / one-point stabilizer in M23.
  * A point stabilizer in M22 is M21 ~= PSL(3,4), acting on the 21 points of
    PG(2,4).
  * The Steiner system S(3,6,22) has 77 hexads.  Hexads through a fixed point,
    with that point removed, are the 21 lines of PG(2,4).

Unit-distance translation:
  A hypothetical 61-edge 22-point unit-distance graph has a degree-4 or degree-5
  deletion ear over a 57-edge or 56-edge 21-core.  The Mathieu residue suggests
  classifying that ear's neighbor set inside PG(2,4):

      degree 5: line-hexad, near-line, or scattered 5-set
      degree 4: punctured-line, 3-collinear, or 4-arc

The output is a design ledger, not a proof of u(22).  Its job is to preserve
the side-channel that raw edge counts forget.

Tournament Analysis:
  * Vertices: proof carriers / ear types, not points.
  * Pairwise observable: design coherence, ear relevance, circle-cap sharpness,
    geometry retention, and scalar-forgetfulness penalty.
  * Switch: majority vote over those proof criteria.
  * Tie Hamiltonian path: listed carrier order.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations
from math import comb


def gf_add(a: int, b: int) -> int:
    return a ^ b


def gf_mul(a: int, b: int) -> int:
    # Elements are c0 + c1*x over F2, encoded as c0 + 2*c1.
    a0, a1 = a & 1, (a >> 1) & 1
    b0, b1 = b & 1, (b >> 1) & 1
    c0 = a0 & b0
    c1 = (a0 & b1) ^ (a1 & b0)
    c2 = a1 & b1
    # x^2 = x + 1 for the primitive polynomial x^2 + x + 1.
    return (c0 ^ c2) | ((c1 ^ c2) << 1)


def gf_inv(a: int) -> int:
    if a == 0:
        raise ZeroDivisionError
    for b in (1, 2, 3):
        if gf_mul(a, b) == 1:
            return b
    raise AssertionError(a)


def gf_dot(a: tuple[int, int, int], b: tuple[int, int, int]) -> int:
    out = 0
    for x, y in zip(a, b):
        out = gf_add(out, gf_mul(x, y))
    return out


def normalize(v: tuple[int, int, int]) -> tuple[int, int, int]:
    if v == (0, 0, 0):
        raise ValueError("zero projective vector")
    pivot = next(x for x in v if x)
    inv = gf_inv(pivot)
    return tuple(gf_mul(inv, x) for x in v)


def pg24_points() -> tuple[tuple[int, int, int], ...]:
    pts = {
        normalize((a, b, c))
        for a in range(4)
        for b in range(4)
        for c in range(4)
        if (a, b, c) != (0, 0, 0)
    }
    return tuple(sorted(pts))


def pg24_lines(points: tuple[tuple[int, int, int], ...]):
    covectors = pg24_points()
    lines = []
    for cov in covectors:
        line = tuple(i for i, p in enumerate(points) if gf_dot(cov, p) == 0)
        lines.append(line)
    return tuple(sorted(set(lines)))


POINTS = pg24_points()
LINES = pg24_lines(POINTS)


def pair_line_index() -> dict[tuple[int, int], int]:
    out = {}
    for idx, line in enumerate(LINES):
        for a, b in combinations(line, 2):
            out[tuple(sorted((a, b)))] = idx
    return out


PAIR_LINE = pair_line_index()


def validate_plane() -> None:
    assert len(POINTS) == 21
    assert len(LINES) == 21
    assert all(len(line) == 5 for line in LINES)
    assert len(PAIR_LINE) == comb(21, 2)
    point_degrees = Counter(p for line in LINES for p in line)
    assert sorted(point_degrees.values()) == [5] * 21


def line_intersection_hist(subset: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    counts = Counter(len(set(subset) & set(line)) for line in LINES)
    return tuple(sorted(counts.items()))


def secant_count(subset: tuple[int, ...]) -> int:
    return sum(1 for line in LINES if len(set(subset) & set(line)) >= 2)


def max_line_intersection(subset: tuple[int, ...]) -> int:
    return max(len(set(subset) & set(line)) for line in LINES)


def classify_subset(subset: tuple[int, ...]) -> str:
    k = len(subset)
    m = max_line_intersection(subset)
    if m == k:
        return f"line_{k}"
    if k == 4 and m == 3:
        return "three_collinear_4"
    if k == 4 and m == 2:
        return "arc_4_no_three"
    if k == 5 and m == 4:
        return "near_line_4_plus_1"
    if k == 5 and m == 3:
        return "three_collinear_5"
    if k == 5 and m == 2:
        return "arc_5_no_three"
    return f"other_{k}_{m}"


def subset_type_table(k: int):
    rows = []
    by_type: dict[str, list[tuple[int, tuple[tuple[int, int], ...]]]] = {}
    for subset in combinations(range(21), k):
        typ = classify_subset(subset)
        by_type.setdefault(typ, []).append((secant_count(subset), line_intersection_hist(subset)))
    for typ, vals in sorted(by_type.items()):
        secants = Counter(x for x, _ in vals)
        hists = Counter(hist for _, hist in vals)
        rows.append(
            {
                "type": typ,
                "count": len(vals),
                "secant_hist": dict(sorted(secants.items())),
                "line_profile_count": len(hists),
                "sample_profile": hists.most_common(1)[0][0],
            }
        )
    return rows


def max_unit_chords_on_unit_circle(k: int) -> int:
    # Unit chords among points on a unit circle occur at angular separation 60
    # degrees.  Extremal finite subsets live on a 6-cycle.
    best = 0
    vertices = range(6)
    edges = {tuple(sorted((i, (i + 1) % 6))) for i in vertices}
    for subset in combinations(vertices, k):
        s = set(subset)
        best = max(best, sum(1 for a, b in edges if a in s and b in s))
    return best


def print_design_header() -> None:
    print("=" * 78)
    print("Unit-distance n=22 Mathieu/M21 residue scout")
    print("=" * 78)
    print("External facts verified from standard references:")
    print("  M24 point stabilizer -> M23; two-point stabilizer -> M22.")
    print("  M22 point stabilizer is M21 ~= PSL(3,4) on PG(2,4)'s 21 points.")
    print("  S(3,6,22) has 77 hexads; through a fixed point are 21 residual")
    print("  5-sets, exactly the lines of PG(2,4).")
    print()
    print("Small unit-distance frontier context:")
    print("  Published small-n data give u(21)=57 and 60 <= u(22) <= 61.")
    print("  A 61-edge 22-point graph has a degree-4 or degree-5 deletion ear:")
    print("    degree 4 -> 57-edge 21-core")
    print("    degree 5 -> 56-edge 21-core")
    print()


def print_plane_audit() -> None:
    print("(A) PG(2,4) residue audit")
    print(f"  points={len(POINTS)} lines={len(LINES)} line_size={len(LINES[0])}")
    print(f"  point_line_degree={Counter(p for line in LINES for p in line).most_common(1)[0][1]}")
    print(f"  pairs covered by unique line={len(PAIR_LINE)}/{comb(21, 2)}")
    print()
    print("  Witt design arithmetic:")
    print(f"    S(3,6,22) hexads = C(22,3)/C(6,3) = {comb(22, 3) // comb(6, 3)}")
    print(f"    hexads through a fixed point = C(21,2)/C(5,2) = {comb(21, 2) // comb(5, 2)}")
    print("    deleting the fixed point turns those 21 hexads into PG(2,4) lines")
    print()


def print_ear_tables() -> None:
    print("(B) Mathieu residue types for degree-4 and degree-5 ears")
    for k in (4, 5):
        print(f"  k={k} ear-neighbor subsets of the 21-point residue")
        for row in subset_type_table(k):
            print(
                f"    {row['type']}: count={row['count']} "
                f"secants={row['secant_hist']} profile={row['sample_profile']}"
            )
        print(
            f"    unit-circle cap: at most {max_unit_chords_on_unit_circle(k)} "
            f"unit chords among the {k} neighbors of one extension point"
        )
        print()


def print_reduction() -> None:
    print("(C) Reduction suggested by the M22 -> M21 break")
    print("  The 22-point problem should be split by the deleted point p:")
    print("    1. degree-5 ear over a 56-core")
    print("       line_5: p plus a PG(2,4) line is a coherent Witt hexad")
    print("       near_line/scattered: the ear spends many residue secants")
    print("    2. degree-4 ear over a 57-core")
    print("       line_4: p plus four of a PG line is a punctured hexad")
    print("       three_collinear/arc: the ear is projective-scattered")
    print()
    print("  Why this is not decorative numerology:")
    print("    M24/M23/M22 are sporadic highly transitive layers; deleting once")
    print("    from M22 exposes the classical projective plane PG(2,4).")
    print("    That is exactly the same shape as the n=22 unit-distance frontier:")
    print("    deciding 61 is not another raw edge-count step, but a one-point")
    print("    extension problem over a 21-point core.")
    print()
    print("  Candidate proof split:")
    print("    coherent ears: check line_5 / line_4 against unit-circle chord caps")
    print("      and Moser cap-endpoint ledgers;")
    print("    scattered ears: use many-secant residue profiles as an obstruction")
    print("      library, analogous to totally-unfaithful graph side channels.")
    print()


@dataclass(frozen=True)
class Carrier:
    name: str
    design_coherence: int
    ear_relevance: int
    circle_cap_sharpness: int
    geometry_retention: int
    computational_burden: int
    scalar_forgetfulness: int


CARRIERS = (
    Carrier("line_5_hexad_ear", 5, 5, 5, 4, 2, 1),
    Carrier("line_4_punctured_hexad_ear", 5, 5, 4, 4, 2, 1),
    Carrier("near_line_4_plus_1", 4, 4, 4, 3, 3, 2),
    Carrier("three_collinear_scattered_ear", 3, 4, 3, 3, 4, 2),
    Carrier("arc_scattered_ear", 2, 4, 2, 3, 5, 2),
    Carrier("raw_21_core_extension", 1, 5, 1, 2, 4, 4),
    Carrier("Moser_27_quantum", 2, 4, 3, 5, 2, 2),
    Carrier("graph_only_F_free_coimage", 0, 3, 0, 1, 5, 5),
)


def votes(a: Carrier, b: Carrier) -> tuple[int, int]:
    criteria = [
        (a.design_coherence > b.design_coherence, b.design_coherence > a.design_coherence),
        (a.ear_relevance > b.ear_relevance, b.ear_relevance > a.ear_relevance),
        (a.circle_cap_sharpness > b.circle_cap_sharpness, b.circle_cap_sharpness > a.circle_cap_sharpness),
        (a.geometry_retention > b.geometry_retention, b.geometry_retention > a.geometry_retention),
        (a.computational_burden < b.computational_burden, b.computational_burden < a.computational_burden),
        (a.scalar_forgetfulness < b.scalar_forgetfulness, b.scalar_forgetfulness < a.scalar_forgetfulness),
    ]
    av = sum(1 for x, y in criteria if x and not y)
    bv = sum(1 for x, y in criteria if y and not x)
    return av, bv


def tournament_adj(items: tuple[Carrier, ...]) -> list[list[int]]:
    n = len(items)
    adj = [[0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        iv, jv = votes(items[i], items[j])
        if iv > jv or (iv == jv and i < j):
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return adj


def directed_3cycles(adj: list[list[int]]) -> int:
    total = 0
    for i, j, k in combinations(range(len(adj)), 3):
        if adj[i][j] and adj[j][k] and adj[k][i]:
            total += 1
        if adj[i][k] and adj[k][j] and adj[j][i]:
            total += 1
    return total


def scc_sizes(adj: list[list[int]]) -> list[int]:
    n = len(adj)

    def reach(start: int, reverse: bool = False) -> set[int]:
        seen = {start}
        stack = [start]
        while stack:
            v = stack.pop()
            for u in range(n):
                edge = adj[u][v] if reverse else adj[v][u]
                if edge and u not in seen:
                    seen.add(u)
                    stack.append(u)
        return seen

    left = set(range(n))
    out = []
    while left:
        root = next(iter(left))
        comp = reach(root) & reach(root, reverse=True)
        out.append(len(comp))
        left -= comp
    return sorted(out, reverse=True)


def hamiltonian_paths(adj: list[list[int]]) -> int:
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for v in range(n):
            cur = dp[mask][v]
            if not cur:
                continue
            for u in range(n):
                if not (mask & (1 << u)) and adj[v][u]:
                    dp[mask | (1 << u)][u] += cur
    return sum(dp[-1])


def print_tournament_analysis() -> None:
    print("(D) Tournament Analysis over proof carriers")
    adj = tournament_adj(CARRIERS)
    scores = [sum(row) for row in adj]
    order = sorted(zip(CARRIERS, scores), key=lambda item: (-item[1], item[0].name))
    print("  vertices: proof carriers / ear types, not points")
    print("  pairwise observable: design coherence, ear relevance, circle cap,")
    print("    geometry retention, burden, scalar-forgetfulness penalty")
    print("  tie Hamiltonian path: listed carrier order")
    for carrier, score in order:
        print(
            f"    score={score} {carrier.name} "
            f"features=({carrier.design_coherence},{carrier.ear_relevance},"
            f"{carrier.circle_cap_sharpness},{carrier.geometry_retention},"
            f"{carrier.computational_burden},{carrier.scalar_forgetfulness})"
        )
    print(f"  score_hist={dict(sorted(Counter(scores).items()))}")
    print(f"  directed_3_cycles={directed_3cycles(adj)}")
    print(f"  scc_sizes={scc_sizes(adj)}")
    print(f"  hamiltonian_paths={hamiltonian_paths(adj)}")
    print()


def main() -> None:
    validate_plane()
    print_design_header()
    print_plane_audit()
    print_ear_tables()
    print_reduction()
    print_tournament_analysis()
    print("SYNTHESIS")
    print("  1. The Mathieu observation suggests a concrete deletion residue:")
    print("     n=22 -> n=21 exposes PG(2,4), not a featureless 21-core.")
    print("  2. A possible 61-edge graph is a degree-4/5 ear problem.  The ear")
    print("     can be coherent (line/punctured-line) or scattered in PG(2,4).")
    print("  3. Coherent ears are small enough to attack with unit-circle chord")
    print("     caps and the Moser cap-endpoint ledger; scattered ears should feed")
    print("     an obstruction library via their many-secant residue profiles.")
    print("  4. This is the unit-distance analogue of the repo's no-leak lesson:")
    print("     raw graph quotients forget the side channel where the proof lives.")


if __name__ == "__main__":
    main()
