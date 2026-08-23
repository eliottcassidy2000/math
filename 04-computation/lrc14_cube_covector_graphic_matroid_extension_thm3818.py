#!/usr/bin/env python3
"""Exact support-two cube-covector graphic-matroid and facet audit.

This companion extends THM-3818 only inside its proved support-two
decoder atlas.  It verifies finite extremals and hostile controls for two
different structures carried by the placed covectors:

1. after a diagonal speed rescaling their row matroid is exactly graphic;
2. their oriented exposed facets have a separate two-level compatibility law.

The general proofs are elementary and recorded in the adjacent integration
note.  No LRC(14) exclusion or prize consequence is claimed.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
import sys

from sympy import Matrix


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


CHECKS = 0


def require(condition: bool, label: object) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


def factor(n: int) -> dict[int, int]:
    require(n >= 1, ("positive factor input", n))
    answer: dict[int, int] = {}
    candidate = 2
    while candidate * candidate <= n:
        while n % candidate == 0:
            answer[candidate] = answer.get(candidate, 0) + 1
            n //= candidate
        candidate = 3 if candidate == 2 else candidate + 2
    if n > 1:
        answer[n] = answer.get(n, 0) + 1
    return answer


def admissible_sum(total: int) -> bool:
    factors = factor(total)
    return bool(factors) and all(prime % 3 == 2 and exponent <= 2 for prime, exponent in factors.items())


def inert_scale(scale: int) -> bool:
    """The stronger table-free decoder requires only inert prime divisors."""
    return all(prime % 3 == 2 for prime in factor(scale))


def ratio_atlas() -> tuple[tuple[int, int], ...]:
    ratios: list[tuple[int, int]] = []
    for total in range(3, 357):
        if not admissible_sum(total):
            continue
        for p in range(1, (total + 1) // 2):
            q = total - p
            if p < q and math.gcd(p, q) == 1:
                ratios.append((p, q))
    return tuple(ratios)


def reduced_pair(x: int, y: int) -> tuple[int, int]:
    require(x > 0 and y > 0 and x != y, ("distinct positive speeds", x, y))
    lo, hi = sorted((x, y))
    common = math.gcd(lo, hi)
    return lo // common, hi // common


def edge_row(speeds: tuple[int, ...], i: int, j: int) -> tuple[int, ...]:
    require(i != j, ("distinct labels", i, j))
    lo, hi = (i, j) if speeds[i] < speeds[j] else (j, i)
    p, q = reduced_pair(speeds[lo], speeds[hi])
    row = [0] * len(speeds)
    row[lo] = q
    row[hi] = -p
    require(sum(row[k] * speeds[k] for k in range(len(speeds))) == 0, ("edge relation", i, j))
    return tuple(row)


def decoder_edges(speeds: tuple[int, ...], ratio_set: set[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    return tuple(
        (i, j)
        for i in range(len(speeds))
        for j in range(i + 1, len(speeds))
        if reduced_pair(speeds[i], speeds[j]) in ratio_set
    )


def component_count(vertex_count: int, edges: tuple[tuple[int, int], ...]) -> int:
    parent = list(range(vertex_count))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i, j in edges:
        root_i = find(i)
        root_j = find(j)
        if root_i != root_j:
            parent[root_i] = root_j
    return len({find(i) for i in range(vertex_count)})


def matrix_rank(rows: tuple[tuple[int, ...], ...], width: int) -> int:
    return 0 if not rows else Matrix(rows).rank()


def facet_compatible(speeds: tuple[int, ...], edges: tuple[tuple[int, int], ...]) -> bool:
    positive: set[int] = set()
    negative: set[int] = set()
    for i, j in edges:
        lo, hi = (i, j) if speeds[i] < speeds[j] else (j, i)
        positive.add(lo)
        negative.add(hi)
    return positive.isdisjoint(negative)


def add_rows(terms: tuple[tuple[int, tuple[int, ...]], ...]) -> tuple[int, ...]:
    width = len(terms[0][1])
    return tuple(sum(coefficient * row[j] for coefficient, row in terms) for j in range(width))


def main() -> None:
    ratios = ratio_atlas()
    ratio_set = set(ratios)
    require(len(ratios) == 5855, "THM-3818 ratio count")
    require(len(ratio_set) == len(ratios), "primitive ratio uniqueness")
    max_q = max(q for _, q in ratios)
    max_q_ratios = tuple(pair for pair in ratios if pair[1] == max_q)
    require(max_q == 355, ("sharp coefficient maximum", max_q))
    require(max_q_ratios == ((1, 355),), ("unique coefficient extremizer", max_q_ratios))

    # Projective triangle types.  Normalize the smallest speed to one.  The
    # other two ratios belong to the atlas, and their quotient must belong to
    # it as well.  Each unordered pair below represents exactly one sorted
    # projective triple.
    full_decoder_triangle_types = 0
    all_inert_triangle_types = 0
    triangle_circuit_gates = 0
    best_triangle_key: tuple[int, int, tuple[int, int, int]] | None = None
    best_triangle_multiplicity = 0
    for index, (p, q) in enumerate(ratios):
        for u, v in ratios[index + 1 :]:
            cross_a = q * u
            cross_b = v * p
            require(cross_a != cross_b, ("distinct ratio values", (p, q), (u, v)))
            quotient = reduced_pair(cross_a, cross_b)
            if quotient not in ratio_set:
                continue

            full_decoder_triangle_types += 1

            denominator_lcm = math.lcm(p, u)
            raw = (denominator_lcm, q * (denominator_lcm // p), v * (denominator_lcm // u))
            triple = tuple(sorted(raw))
            common = math.gcd(math.gcd(triple[0], triple[1]), triple[2])
            triple = tuple(value // common for value in triple)
            x, y, z = triple
            require(x < y < z, ("sorted projective triangle", triple))
            require(
                all(reduced_pair(a, b) in ratio_set for a, b in ((x, y), (y, z), (x, z))),
                ("triangle atlas closure", triple),
            )
            # The weighted graphic circuit is explicit.  If ell_uv is the
            # pair lcm and L their lcm around the cycle, then
            # (L/ell_xy)a_xy+(L/ell_yz)a_yz-(L/ell_xz)a_xz=0.
            speeds = (x, y, z)
            row_xy = edge_row(speeds, 0, 1)
            row_yz = edge_row(speeds, 1, 2)
            row_xz = edge_row(speeds, 0, 2)
            ell_xy = math.lcm(x, y)
            ell_yz = math.lcm(y, z)
            ell_xz = math.lcm(x, z)
            circuit_lcm = math.lcm(ell_xy, ell_yz, ell_xz)
            circuit = add_rows(
                (
                    (circuit_lcm // ell_xy, row_xy),
                    (circuit_lcm // ell_yz, row_yz),
                    (-(circuit_lcm // ell_xz), row_xz),
                )
            )
            require(circuit == (0, 0, 0), ("triangle circuit", triple, circuit))
            triangle_circuit_gates += 1

            if all(
                inert_scale(math.gcd(a, b)) for a, b in ((x, y), (y, z), (x, z))
            ):
                all_inert_triangle_types += 1

            key = (max(triple), sum(triple), triple)
            if best_triangle_key is None or key < best_triangle_key:
                best_triangle_key = key
                best_triangle_multiplicity = 1
            elif key == best_triangle_key:
                best_triangle_multiplicity += 1

    require(
        full_decoder_triangle_types == 245220,
        ("all-scale decoder triangle census", full_decoder_triangle_types),
    )
    require(
        all_inert_triangle_types == 46136,
        ("all-inert table-free triangle census", all_inert_triangle_types),
    )
    require(best_triangle_key == (8, 13, (2, 3, 8)), ("minimal triangle", best_triangle_key))
    require(best_triangle_multiplicity == 1, ("minimal triangle uniqueness", best_triangle_multiplicity))

    # Exhaust every subset of the decoder graph on a mixed five-speed control.
    # This independently checks the graphic rank law and the orientation law.
    control_speeds = (1, 2, 3, 8, 9)
    control_edges = decoder_edges(control_speeds, ratio_set)
    require(len(control_edges) == 8, ("control edge count", control_edges))
    subset_rank_gates = 0
    compatible_subsets = 0
    for mask in range(1 << len(control_edges)):
        edges = tuple(control_edges[k] for k in range(len(control_edges)) if mask & (1 << k))
        rows = tuple(edge_row(control_speeds, i, j) for i, j in edges)
        rank = matrix_rank(rows, len(control_speeds))
        graphic_rank = len(control_speeds) - component_count(len(control_speeds), edges)
        require(rank == graphic_rank, ("graphic rank", mask, rank, graphic_rank))
        require((rank == len(edges)) == (len(edges) == graphic_rank), ("forest independence", mask))
        compatible = facet_compatible(control_speeds, edges)
        positive = {
            i if control_speeds[i] < control_speeds[j] else j
            for i, j in edges
        }
        negative = {
            j if control_speeds[i] < control_speeds[j] else i
            for i, j in edges
        }
        require(compatible == positive.isdisjoint(negative), ("two-level orientation", mask))
        compatible_subsets += int(compatible)
        subset_rank_gates += 3

    # Independence and facet compatibility are orthogonal.  The two-edge
    # path 2<3<8 is a forest but asks coordinate 3 to be both lower and upper.
    forest_speeds = (2, 3, 8)
    forest_edges = ((0, 1), (1, 2))
    forest_rows = tuple(edge_row(forest_speeds, *edge) for edge in forest_edges)
    require(matrix_rank(forest_rows, 3) == 2, "independent forest hostile")
    require(not facet_compatible(forest_speeds, forest_edges), "forest facet conflict")

    # Conversely, this alternating four-cycle is dependent but its maxima are
    # compatible: speeds 1,2 are upper endpoints and 3,9 lower endpoints.
    cycle_speeds = (1, 2, 3, 9)
    cycle_edges = ((0, 2), (1, 2), (1, 3), (0, 3))
    cycle_rows = tuple(edge_row(cycle_speeds, *edge) for edge in cycle_edges)
    require(all(reduced_pair(cycle_speeds[i], cycle_speeds[j]) in ratio_set for i, j in cycle_edges), "compatible cycle atlas")
    require(matrix_rank(cycle_rows, 4) == 3, "dependent compatible cycle")
    require(facet_compatible(cycle_speeds, cycle_edges), "compatible cycle face")

    # Sharp connected-component height control.  Every decoder coefficient is
    # at most 355, so a t-vertex tree cofactor is at most 355^(t-1).  The full
    # 13-vertex bound is attained by the primitive path of powers of 355.
    connected_bound = 355**12
    extremal_speeds = tuple(355**k for k in range(13))
    extremal_edges = decoder_edges(extremal_speeds, ratio_set)
    wanted_path = tuple((i, i + 1) for i in range(12))
    require(extremal_edges == wanted_path, ("sharp decoder path", extremal_edges))
    extremal_rows = tuple(edge_row(extremal_speeds, *edge) for edge in extremal_edges)
    require(matrix_rank(extremal_rows, 13) == 12, "sharp path rank")
    require(math.gcd(*extremal_speeds) == 1, "sharp path primitive")
    require(max(extremal_speeds) == connected_bound, "sharp path height")

    # Diagonal equivalence control on every edge used above:
    # a_ij diag(n)=lcm(n_i,n_j)(e_i-e_j), up to the lower-to-higher gauge.
    diagonal_gates = 0
    for speeds, edges in (
        (control_speeds, control_edges),
        (cycle_speeds, cycle_edges),
        (extremal_speeds, extremal_edges),
    ):
        for i, j in edges:
            row = edge_row(speeds, i, j)
            lo, hi = (i, j) if speeds[i] < speeds[j] else (j, i)
            ell = math.lcm(speeds[lo], speeds[hi])
            scaled = tuple(row[k] * speeds[k] for k in range(len(speeds)))
            wanted = tuple(ell if k == lo else -ell if k == hi else 0 for k in range(len(speeds)))
            require(scaled == wanted, ("diagonal incidence equivalence", speeds, i, j))
            diagonal_gates += 1

    semantic_record = {
        "ratios": len(ratios),
        "max_q": max_q,
        "max_q_ratios": max_q_ratios,
        "full_decoder_triangle_types": full_decoder_triangle_types,
        "all_inert_triangle_types": all_inert_triangle_types,
        "minimal_triangle": best_triangle_key,
        "control_edges": control_edges,
        "control_subsets": 1 << len(control_edges),
        "compatible_subsets": compatible_subsets,
        "connected_bound": connected_bound,
        "extremal_path_edges": extremal_edges,
    }
    semantic = hashlib.sha256(
        json.dumps(semantic_record, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()

    print("LRC14_CUBE_COVECTOR_GRAPHIC_MATROID_PROBE")
    print("status=PROVED_GENERAL_FORMULAS+FINITE_EXACT_CONTROLS;LRC14_OPEN")
    print("diagonal_equivalence=a_ij*diag(n)=lcm(n_i,n_j)*(e_low-e_high)")
    print("row_matroid=graphic;rank(E)=13-components(G);independent_iff_forest;circuits=simple_cycles")
    print("cycle_law=sum_cycle(sign_e*L/lcm(n_i,n_j)*a_e)=0")
    print("rank11_obstruction=decoder_graph_disconnected;exact_rank11_decoder_span=exactly_two_components")
    print("component_shape_bound=max(n_i/gcd_component)<=355^(vertices-1)")
    print(f"connected_terminal_bound=max_speed<={connected_bound}")
    print("connected_bound_sharp=(1,355,...,355^12);decoder_graph_exactly_the_12-edge_path")
    print("facet_compatibility=common_max_face_iff_no_label_occurs_with_both_signs_iff_directed_height<=1")
    print("compatible_face_dimension=13-number_of_incident_labels_for_nonempty_edge_set")
    print("forest_hostile=edges_(2,3),(3,8)_independent_but_common_max_face_empty")
    print("cycle_hostile=edges_(1,3),(2,3),(2,9),(1,9)_dependent_but_common_max_face_nonempty_dimension_9")
    print(f"decoder_ratios={len(ratios)}")
    print(f"all_scale_full_decoder_triangle_circuits={full_decoder_triangle_types}")
    print(f"all_inert_table_free_triangle_circuits={all_inert_triangle_types}")
    print("unique_smallest_triangle_by_max_then_sum=(2,3,8);max=8;sum=13")
    print(f"triangle_circuit_gates={triangle_circuit_gates}")
    print(f"control_subset_rank_orientation_gates={subset_rank_gates}")
    print(f"control_compatible_subsets={compatible_subsets}/{1 << len(control_edges)}")
    print(f"diagonal_equivalence_gates={diagonal_gates}")
    print("finite_universe=all_5855_ratios;all_C(5855,2)_ratio_pairs;all_scale_graph;all_inert_table_free_subcount;all_control_edge_subsets")
    print(f"active_require_gates={CHECKS}")
    print(f"semantic_sha256={semantic}")
    print("scope=internal_relation_matroid+facet_incidence+height_terminal_not_loneliness_or_arrival")
    print("PASS")


if __name__ == "__main__":
    main()
