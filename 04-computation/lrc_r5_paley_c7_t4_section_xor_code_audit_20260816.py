#!/usr/bin/env python3
"""Exact section-code audit for the folded Paley C7 carrier.

Folding the Paley tournament by r~-r gives a symmetric weighted four-block
quotient, not a tournament.  Choosing one representative in each nonzero
block gives eight gauge-dependent T4 sections.  This script classifies their
six edge bits, tournament types, cyclic triangles, and skew determinants using
only integer/Boolean arithmetic.

The resulting affine F2^3 code is a typing sidecar for the folded-C7/K4
transporter obstruction.  It is not an ancestry map, response channel,
physical current, H1 class, or LRC(14) certificate.
"""

from __future__ import annotations

from hashlib import sha256
from itertools import combinations, product
from pathlib import Path


EXPECTED_SEMANTIC_SHA256: str | None = "85640c02220b246781bab2e447a6a2f9fec7ca912d27c56377fbcd943bf041bc"

RESIDUES = frozenset((1, 2, 4))
ORBITS = ((0,), (1, 6), (2, 5), (3, 4))
BASE_REPRESENTATIVES = (1, 2, 3)
EDGES = tuple(combinations(range(4), 2))
MATCHINGS = (
    ((0, 1), (2, 3)),
    ((0, 2), (1, 3)),
    ((0, 3), (1, 2)),
)


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


def arc(source: int, target: int) -> bool:
    return (target - source) % 7 in RESIDUES


def folded_arc_counts() -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple(sum(arc(source, target) for source in left for target in right) for right in ORBITS)
        for left in ORBITS
    )


def edge_code(bits: tuple[int, int, int]) -> tuple[int, ...]:
    signs = tuple(1 if bit == 0 else -1 for bit in bits)
    vertices = (0,) + tuple((sign * base) % 7 for sign, base in zip(signs, BASE_REPRESENTATIVES))
    # Bit zero means the lower-labelled endpoint points to the higher-labelled
    # endpoint; bit one means the reverse orientation.
    return tuple(0 if arc(vertices[i], vertices[j]) else 1 for i, j in EDGES)


def xor_rank(rows: tuple[tuple[int, ...], ...]) -> int:
    work = [list(row) for row in rows]
    rank = 0
    columns = len(work[0]) if work else 0
    for column in range(columns):
        pivot = next((index for index in range(rank, len(work)) if work[index][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        for index in range(len(work)):
            if index != rank and work[index][column]:
                work[index] = [left ^ right for left, right in zip(work[index], work[rank])]
        rank += 1
    return rank


def boundary(edge_vector: tuple[int, ...]) -> tuple[int, ...]:
    vertices = [0, 0, 0, 0]
    for value, (left, right) in zip(edge_vector, EDGES):
        if value:
            vertices[left] ^= 1
            vertices[right] ^= 1
    return tuple(vertices)


def tournament_row(bits: tuple[int, int, int]) -> tuple[object, ...]:
    code = edge_code(bits)
    expected = (bits[0], bits[1], 1 ^ bits[2], bits[0], bits[2], bits[1])
    require(code == expected, (bits, code, expected))

    outdegrees = [0, 0, 0, 0]
    skew = [[0] * 4 for _ in range(4)]
    for bit, (left, right) in zip(code, EDGES):
        if bit == 0:
            outdegrees[left] += 1
            skew[left][right] = 1
            skew[right][left] = -1
        else:
            outdegrees[right] += 1
            skew[left][right] = -1
            skew[right][left] = 1

    cyclic_triangles = 0
    for vertices in combinations(range(4), 3):
        local = {vertex: 0 for vertex in vertices}
        for left, right in combinations(vertices, 2):
            edge_index = EDGES.index((left, right))
            if code[edge_index] == 0:
                local[left] += 1
            else:
                local[right] += 1
        cyclic_triangles += tuple(sorted(local.values())) == (1, 1, 1)

    pfaffian = (
        skew[0][1] * skew[2][3]
        - skew[0][2] * skew[1][3]
        + skew[0][3] * skew[1][2]
    )
    determinant = pfaffian * pfaffian
    x, y, z = bits
    exceptional = z ^ (x & y) ^ (x & z) ^ (y & z)
    require(determinant == 1 + 8 * exceptional, (bits, determinant, exceptional))
    require(cyclic_triangles == 2 - exceptional, (bits, cyclic_triangles, exceptional))

    score = tuple(sorted(outdegrees))
    if score == (1, 1, 2, 2):
        kind = "strong"
    elif score == (1, 1, 1, 3):
        kind = "cycle_plus_source"
    elif score == (0, 2, 2, 2):
        kind = "cycle_plus_sink"
    elif score == (0, 1, 2, 3):
        kind = "transitive"
    else:
        raise RuntimeError((bits, score, "unknown T4 type"))

    return bits, code, tuple(outdegrees), cyclic_triangles, pfaffian, determinant, exceptional, kind


def main() -> None:
    counts = folded_arc_counts()
    expected_counts = (
        (0, 1, 1, 1),
        (1, 1, 2, 2),
        (1, 2, 1, 2),
        (1, 2, 2, 1),
    )
    require(counts == expected_counts, counts)
    require(
        all(counts[left][right] == counts[right][left] for left in range(4) for right in range(4)),
        "fold must lose orientation",
    )

    rows = tuple(tournament_row(bits) for bits in product((0, 1), repeat=3))
    codes = tuple(row[1] for row in rows)
    require(len(set(codes)) == 8, "section code injectivity")

    base = codes[0]
    translated = tuple(tuple(value ^ offset for value, offset in zip(code, base)) for code in codes)
    generators = (
        tuple(left ^ right for left, right in zip(edge_code((1, 0, 0)), base)),
        tuple(left ^ right for left, right in zip(edge_code((0, 1, 0)), base)),
        tuple(left ^ right for left, right in zip(edge_code((0, 0, 1)), base)),
    )
    require(xor_rank(generators) == 3 and len(set(translated)) == 8, (generators, translated))
    require(
        all(
            tuple(left ^ right for left, right in zip(one, two)) in set(translated)
            for one in translated
            for two in translated
        ),
        "translated codes form F2^3",
    )
    require(
        generators == (
            (1, 0, 0, 1, 0, 0),
            (0, 1, 0, 0, 0, 1),
            (0, 0, 1, 0, 1, 0),
        ),
        generators,
    )

    gauge_space = frozenset(translated)
    cycle_space = frozenset(vector for vector in product((0, 1), repeat=6) if boundary(vector) == (0, 0, 0, 0))
    reduced_vertices = frozenset(vector for vector in product((0, 1), repeat=4) if sum(vector) % 2 == 0)
    require(len(cycle_space) == 8 and xor_rank(tuple(cycle_space)) == 3, "K4 cycle rank")
    require(gauge_space.intersection(cycle_space) == {(0, 0, 0, 0, 0, 0)}, "section gauge is not cycle data")
    gauge_boundaries = frozenset(boundary(vector) for vector in gauge_space)
    require(gauge_boundaries == reduced_vertices, (gauge_boundaries, reduced_vertices))
    direct_sums = {
        tuple(left ^ right for left, right in zip(gauge, cycle))
        for gauge in gauge_space
        for cycle in cycle_space
    }
    require(len(direct_sums) == 64, "C1(K4)=cycle_space direct_sum section_gauge")

    # Global negation adds 111 to section bits and reverses all six arcs.
    for bits, code, *_ in rows:
        complement_bits = tuple(bit ^ 1 for bit in bits)
        complement_code = edge_code(complement_bits)
        require(complement_code == tuple(bit ^ 1 for bit in code), (bits, "global reversal"))

    kinds = tuple(row[-1] for row in rows)
    require(kinds.count("strong") == 6, kinds)
    require(kinds.count("cycle_plus_source") == 1, kinds)
    require(kinds.count("cycle_plus_sink") == 1, kinds)
    require(kinds.count("transitive") == 0, kinds)
    exceptional_bits = tuple(row[0] for row in rows if row[-2])
    require(exceptional_bits == ((0, 0, 1), (1, 1, 0)), exceptional_bits)

    # Multiplication by two cycles the three nonzero fold blocks.  The Walsh
    # side has three perfect matchings, but the section gauge generators are
    # three disjoint two-edge paths, so no canonical carrier identification
    # is supplied by the common cardinality three.
    times_two = tuple(next(index for index, orbit in enumerate(ORBITS) if (2 * ORBITS[source][0]) % 7 in orbit) for source in (1, 2, 3))
    require(times_two == (2, 3, 1), times_two)
    matching_edge_indices = tuple(tuple(EDGES.index(edge) for edge in matching) for matching in MATCHINGS)
    path_edge_indices = tuple(tuple(index for index, value in enumerate(generator) if value) for generator in generators)
    require(path_edge_indices == ((0, 3), (1, 5), (2, 4)), path_edge_indices)
    require(set(path_edge_indices) != set(matching_edge_indices), "path gauge is not Walsh matching data")

    proof = (
        "folded_quotient_is_symmetric_weighted_looped_K4_support_not_tournament",
        "section_edge_code=(x,y,1_XOR_z,x,z,y)",
        "eight_sections_form_affine_F2^3_inside_six_edge_bits",
        "global_negation_reverses_all_six_arcs",
        "six_strong_plus_one_source_cycle_plus_one_sink_cycle;no_transitive_section",
        "exceptional=z_XOR_xy_XOR_xz_XOR_yz",
        "skew_determinant=1+8*exceptional",
        "C1_K4=cycle_space_direct_sum_section_gauge;boundary_isomorphism_on_gauge",
        "invertible_tournament_left_mixing_cannot_raise_source_rank_three",
    )
    scope = (
        "orientation_requires_noncanonical_section_gauge",
        "quadratic_XOR_type_is_not_a_linear_Walsh_response_channel",
        "no_ancestry_Boolean_realization_current_H1_bispectrum_or_LRC14_claim",
    )
    semantic_surface = (
        counts,
        rows,
        generators,
        tuple(sorted(cycle_space)),
        tuple(sorted(gauge_boundaries)),
        times_two,
        matching_edge_indices,
        path_edge_indices,
        proof,
        scope,
    )
    semantic_sha256 = sha256(repr(semantic_surface).encode("utf-8")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(
            semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
            (semantic_sha256, EXPECTED_SEMANTIC_SHA256),
        )

    source = Path(__file__).resolve()
    print("PALEY C7 FOLD: T4 SECTION XOR CODE")
    print("status=FINITE_EXACT_STRUCTURE_SIDECAR;scope=folded_transporter_typing")
    print(f"fold=(orbits={ORBITS},arc_counts={counts},off_diagonal=bidirected_K4)")
    print(f"edge_order={EDGES};section_formula=(x,y,1_XOR_z,x,z,y)")
    print(f"affine_code=(base={base},generators={generators},dimension=3,size=8)")
    print(f"section_rows={rows}")
    print(f"type_census=(strong={kinds.count('strong')},source_cycle={kinds.count('cycle_plus_source')},sink_cycle={kinds.count('cycle_plus_sink')},transitive={kinds.count('transitive')})")
    print(f"xor_law=(exceptional_bits={exceptional_bits},E=z_XOR_xy_XOR_xz_XOR_yz,det_skew=1+8E,cyclic_triangles=2-E)")
    print(f"three_grammars=(times2_orbit_cycle={times_two},walsh_matching_edge_indices={matching_edge_indices},section_path_edge_indices={path_edge_indices})")
    print("chain_typing=(dim_C1=6,dim_cycle=3,dim_section_gauge=3,gauge_intersection_cycle=0,boundary_on_gauge=isomorphism,C1=Z1_direct_sum_gauge)")
    print("rank_boundary=all_section_skew_matrices_are_invertible_left_mixers_but_rank(S*X*C)<=rank(X)=3")
    print(f"proof={proof}")
    print(f"scope={scope}")
    print(f"semantic_sha256={semantic_sha256}")
    print(f"script_sha256_lf={lf_sha256(source)}")
    print("reproducibility=no_assert_truth_gates;integer_boolean_only;no_randomness;no_elapsed_fields")
    print("commands=python -B 04-computation/lrc_r5_paley_c7_t4_section_xor_code_audit_20260816.py;python -B -O 04-computation/lrc_r5_paley_c7_t4_section_xor_code_audit_20260816.py")
    print("PASS")


if __name__ == "__main__":
    main()
