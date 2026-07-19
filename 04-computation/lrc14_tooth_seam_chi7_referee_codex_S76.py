#!/usr/bin/env python3
"""Exact referee for THM-1156's tooth-seam/chi7 bipartition.

At radius 1/14 the m-th open danger tooth of speed a has endpoints

    L_a(m) = (14m-1)/(14a),   R_a(m) = (14m+1)/(14a).

The script checks, without floating point, that an oriented exact seam
R_a(m)=L_b(n) exists precisely when the gcd-reduced speeds A=a/g and
B=b/g satisfy 14 | A+B.  It verifies the complete seam-phase formula,
the 2- and 7-adic layer invariants, the chi_7 sign flip, and the resulting
triangle-free/Fano obstruction.  A deliberately tie-completed tournament
is reported only as a lossy diagnostic.

No assertion is used, so ``python -O`` replays every check.
"""

from collections import Counter, deque
from fractions import Fraction as F
from itertools import combinations
from math import gcd


BOUND = 196
TOURNAMENT_BOUND = 14


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def valuation(value: int, prime: int) -> int:
    exponent = 0
    while value % prime == 0:
        value //= prime
        exponent += 1
    return exponent


def seven_free_part(value: int) -> int:
    while value % 7 == 0:
        value //= 7
    return value


def chi7(value: int) -> int:
    residue = value % 7
    require(residue != 0, "chi7 called on a multiple of seven")
    return 1 if residue in (1, 2, 4) else -1


def seam_criterion(a: int, b: int) -> bool:
    common = gcd(a, b)
    return (a // common + b // common) % 14 == 0


def right_endpoints(speed: int) -> frozenset[F]:
    return frozenset(F(14 * tooth + 1, 14 * speed) % 1 for tooth in range(speed))


def left_endpoints(speed: int) -> frozenset[F]:
    return frozenset(F(14 * tooth - 1, 14 * speed) % 1 for tooth in range(speed))


def oriented_edge(a: int, b: int) -> tuple[int, int]:
    """Tie-complete the seam semigraph to a tournament.

    Genuine seam edges point from chi7=+1 to chi7=-1.  Every nonedge is a
    tie and is broken by numerical order.  The completion is intentionally
    not used by the theorem.
    """

    require(a < b, "oriented_edge expects a<b")
    if not seam_criterion(a, b):
        return a, b
    color_a = chi7(seven_free_part(a))
    color_b = chi7(seven_free_part(b))
    require(color_a == -color_b, "seam edge failed the chi7 switch")
    return (a, b) if color_a == 1 else (b, a)


def strongly_connected_component_sizes(
    vertices: tuple[int, ...], edges: set[tuple[int, int]]
) -> tuple[int, ...]:
    def reachable(source: int, reverse: bool = False) -> set[int]:
        seen = {source}
        queue = deque([source])
        while queue:
            vertex = queue.popleft()
            for candidate in vertices:
                edge = (candidate, vertex) if reverse else (vertex, candidate)
                if edge in edges and candidate not in seen:
                    seen.add(candidate)
                    queue.append(candidate)
        return seen

    unseen = set(vertices)
    sizes: list[int] = []
    while unseen:
        source = min(unseen)
        component = reachable(source) & reachable(source, reverse=True)
        sizes.append(len(component))
        unseen -= component
    return tuple(sorted(sizes, reverse=True))


def hamiltonian_path_count(vertices: tuple[int, ...], edges: set[tuple[int, int]]) -> int:
    index = {vertex: position for position, vertex in enumerate(vertices)}
    size = 1 << len(vertices)
    dynamic = [[0] * len(vertices) for _ in range(size)]
    for vertex in vertices:
        dynamic[1 << index[vertex]][index[vertex]] = 1
    for mask in range(1, size):
        for last in vertices:
            last_index = index[last]
            count = dynamic[mask][last_index]
            if count == 0:
                continue
            for nxt in vertices:
                next_index = index[nxt]
                if mask & (1 << next_index) == 0 and (last, nxt) in edges:
                    dynamic[mask | (1 << next_index)][next_index] += count
    return sum(dynamic[-1])


def main() -> None:
    rights = {speed: right_endpoints(speed) for speed in range(1, BOUND + 1)}
    lefts = {speed: left_endpoints(speed) for speed in range(1, BOUND + 1)}

    seam_edges: set[tuple[int, int]] = set()
    seam_phase_events = 0
    reduced_direction_pairs: set[tuple[int, int]] = set()
    layer_counts: Counter[tuple[int, int]] = Counter()
    for a, b in combinations(range(1, BOUND + 1), 2):
        common = gcd(a, b)
        reduced_a, reduced_b = a // common, b // common
        actual = rights[a] & lefts[b]
        expected = seam_criterion(a, b)
        require(bool(actual) == expected, f"seam criterion failed at {(a, b)}")
        require(len(actual) == (common if expected else 0),
                f"seam multiplicity failed at {(a, b)}")
        if not expected:
            continue

        seam_edges.add((a, b))
        seam_phase_events += len(actual)
        reduced_direction_pairs.add((reduced_a, reduced_b))
        require(gcd(reduced_a, 14) == gcd(reduced_b, 14) == 1,
                f"reduced seam directions are not units at {(a, b)}")
        require(valuation(a, 2) == valuation(b, 2),
                f"2-adic seam layer changed at {(a, b)}")
        require(valuation(a, 7) == valuation(b, 7),
                f"7-adic seam layer changed at {(a, b)}")
        require(chi7(seven_free_part(a)) == -chi7(seven_free_part(b)),
                f"chi7 did not flip at {(a, b)}")

        inverse = pow(reduced_a, -1, 14)
        constructed = frozenset(
            F(inverse + 14 * index, 14 * common) % 1
            for index in range(common)
        )
        require(actual == constructed, f"seam phase formula failed at {(a, b)}")
        layer_counts[(valuation(a, 2), valuation(a, 7))] += 1

    triangle_count = 0
    maximum_triple_edges = 0
    for a, b, c in combinations(range(1, BOUND + 1), 3):
        edge_count = sum(
            pair in seam_edges for pair in ((a, b), (a, c), (b, c))
        )
        maximum_triple_edges = max(maximum_triple_edges, edge_count)
        triangle_count += edge_count == 3
    require(triangle_count == 0, "exact-seam graph contains a triangle")
    require(maximum_triple_edges == 2, "unexpected triple seam-edge maximum")

    # A standard labelled Fano plane has seven 3-point lines.  Triangle
    # freeness is label-independent, hence every line misses a seam edge.
    fano_lines = (
        (1, 2, 3), (1, 4, 5), (1, 6, 7), (2, 4, 6),
        (2, 5, 7), (3, 4, 7), (3, 5, 6),
    )
    require(len(fano_lines) == 7 and all(len(line) == 3 for line in fano_lines),
            "Fano incidence changed")

    vertices = tuple(range(1, TOURNAMENT_BOUND + 1))
    tournament_edges = {
        oriented_edge(a, b) for a, b in combinations(vertices, 2)
    }
    require(len(tournament_edges) == len(vertices) * (len(vertices) - 1) // 2,
            "tie completion is not a tournament")
    scores = tuple(
        sorted(sum((vertex, other) in tournament_edges for other in vertices if other != vertex)
               for vertex in vertices)
    )
    directed_triangles = sum(
        (a, b) in tournament_edges
        and (b, c) in tournament_edges
        and (c, a) in tournament_edges
        or (a, c) in tournament_edges
        and (c, b) in tournament_edges
        and (b, a) in tournament_edges
        for a, b, c in combinations(vertices, 3)
    )
    scc_sizes = strongly_connected_component_sizes(vertices, tournament_edges)
    path_count = hamiltonian_path_count(vertices, tournament_edges)
    require(path_count > 0, "tournament lost every Hamiltonian path")

    degree_histogram = Counter()
    for vertex in range(1, BOUND + 1):
        degree_histogram[sum(vertex in edge for edge in seam_edges)] += 1

    print("THM-1156 exact tooth-seam / chi7 referee")
    print(f"speed audit: 1 <= a < b <= {BOUND}")
    print("tooth convention: open ((14m-1)/(14a),(14m+1)/(14a))")
    print("exact algebra: R_a(m)=L_b(n) iff 14*(a*n-b*m)=a+b")
    print("projective criterion: a seam exists iff 14 divides a/g+b/g")
    print(f"seam pairs: {len(seam_edges)}")
    print(f"seam phase events: {seam_phase_events}")
    print(f"distinct reduced direction pairs: {len(reduced_direction_pairs)}")
    print("phase formula: t=(A^(-1) mod 14 + 14j)/(14g), 0<=j<g")
    print("layer law: nu_2(a)=nu_2(b), nu_7(a)=nu_7(b)")
    print("quadratic switch: chi7(a/7^nu7(a))=-chi7(b/7^nu7(b))")
    print(f"occupied (nu2,nu7) layers: {len(layer_counts)}")
    print("seam-graph degree histogram: " + str(tuple(sorted(degree_histogram.items()))))
    print(f"seam triangles: {triangle_count}; maximum edges on any triple: {maximum_triple_edges}")
    print("Fano consequence: each of seven 3-point lines contains a nonseam pair")
    print("Tournament Analysis (lossy tie completion on speeds 1..14):")
    print("  observable: exact seam; switch: chi7 +/-; nonseams tie-broken by speed")
    print(f"  score sequence: {scores}")
    print(f"  directed 3-cycles: {directed_triangles}")
    print(f"  SCC sizes: {scc_sizes}")
    print(f"  Hamiltonian-path count: {path_count}")
    print("  edge flips from numerical order: " +
          str(sum((b, a) in tournament_edges for a, b in combinations(vertices, 2))))
    print("faithful vertices: tooth-boundary seam events, not runners alone")
    print("faithful data: (nu2,nu7)-layer + chi7 side + rational phase + tooth addresses")
    print("destroyed by tournament: nonedge gap/overlap sign, metric size, repeated chronology, coverage")
    print("VERDICT: exact seams form a layered bipartite graph; no odd seam cycle exists")
    print("SCOPE: structural zero-overlap obstruction only; no localized eta>0 is claimed")


if __name__ == "__main__":
    main()
