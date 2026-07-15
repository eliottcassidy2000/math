#!/usr/bin/env python3
"""Exact Fano/chi_7 bridge and obstruction for the THM-741 flood tail.

This script deliberately changes vertices twice.

1.  A THM-741 body is viewed through its uncovered divisibility obligations
    q in {2,...,14}, rather than through its runners.  The three already-added
    speeds v1,v2,v3 give each covered obligation a nonzero three-bit incidence
    mask.  Those seven masks are the points of PG(2,2).
2.  A THM-841 Farey-gap ladder is viewed through individual endpoint-needle
    atoms (prefix/suffix intervals), rather than through speeds.  Three such
    atoms trace a weakly monotone path in the three-cube.

The computation certifies four exact statements:

* the 21 flood bodies are exactly {8,...,14} plus an edge of K_7, and the
  edges partition into the seven triples belonging to the Fano lines;
* Boolean B_3 inclusion-exclusion has sign word ++-+--+, exactly the extended
  Legendre chi_7 word, and its three negative masks form one Fano line;
* none of the 167 nonidentity Fano automorphisms preserves even one of the
  exact THM-741 flood-body observables r(E), m(E), or V1(E), so this Fano
  organization is not a symmetry quotient of the numerical flood work;
* r(E) is a seven-point additive potential, whereas m(E) and V1(E) have full
  21-dimensional edge orbit spans and survive outside the point/Fano rank-13
  marginal space;
* three one-sided endpoint needles realize at most four masks and at most two
  points of the negative Fano line, including coincident switching times.

Tournament Analysis uses the 21 flood bodies (equivalently K_7 edges) as
vertices.  Its pairwise observable is the exact THM-741 threshold difference
V1(F)-V1(E); the proof-work gauge points from smaller to larger V1, with
lexicographic edge order as an unused tie switch and as the reference gauge.
This keeps predicted root-search extent but destroys the lower-tree interval
geometry.  The obligation-mask quotient instead keeps the modular lcm filter
and destroys exact interval geometry; the endpoint-needle quotient keeps local
switch order and destroys speed identity and cross-gap gluing.

Scope guardrail: the last statement is about endpoint-needle atoms.  A speed
can contribute one needle from each end of a cyclic Farey gap, so the statement
must not be promoted to arbitrary triples of whole speed predicates.  Likewise,
the script organizes and diagnoses the j=4 flood tail; it does not close it.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
from fractions import Fraction
from itertools import combinations, permutations, product
from pathlib import Path


HERE = Path(__file__).resolve().parent
THM741_PATH = HERE / "lrc14_thm741_2002_body_j4_tree_kps_S128c5.py"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_thm741():
    spec = importlib.util.spec_from_file_location("thm741_flood_dependency", THM741_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def uncovered_obligations(body: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(q for q in range(2, 15) if not any(w % q == 0 for w in body))


def fano_lines() -> tuple[tuple[int, int, int], ...]:
    """PG(2,2) on the nonzero binary vectors labelled 1,...,7."""
    return tuple(
        sorted(
            {
                tuple(sorted((a, b, a ^ b)))
                for a, b in combinations(range(1, 8), 2)
            }
        )
    )


def gl32_maps() -> tuple[tuple[int, ...], ...]:
    """All invertible linear maps F_2^3 -> F_2^3 as point permutations."""
    maps = []
    for x, y, z in product(range(1, 8), repeat=3):
        span = {x, y, z, x ^ y, x ^ z, y ^ z, x ^ y ^ z}
        if len(span) != 7:
            continue
        image = [0]
        for a in range(1, 8):
            image.append(
                (x if a & 1 else 0)
                ^ (y if a & 2 else 0)
                ^ (z if a & 4 else 0)
            )
        maps.append(tuple(image))
    return tuple(maps)


def ordered_event_blocks() -> tuple[tuple[tuple[int, ...], ...], ...]:
    """All weak switch orders of three needles (duplicates intentionally removed)."""
    paths = set()
    for order in permutations(range(3)):
        # Compositions of three encode distinct or coincident switch times.
        for sizes in ((3,), (1, 2), (2, 1), (1, 1, 1)):
            blocks = []
            cursor = 0
            for size in sizes:
                blocks.append(tuple(order[cursor : cursor + size]))
                cursor += size
            paths.add(tuple(blocks))
    return tuple(sorted(paths))


def needle_path(initial_mask: int, blocks: tuple[tuple[int, ...], ...]) -> tuple[int, ...]:
    masks = [initial_mask]
    mask = initial_mask
    for block in blocks:
        for bit in block:
            mask ^= 1 << bit
        masks.append(mask)
    return tuple(masks)


def rational_rank(rows: list[list[Fraction | int]]) -> int:
    """Exact row rank over Q."""
    if not rows:
        return 0
    matrix = [[Fraction(value) for value in row] for row in rows]
    height, width = len(matrix), len(matrix[0])
    pivot_row = 0
    for column in range(width):
        pivot = next(
            (row for row in range(pivot_row, height) if matrix[row][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        scale = matrix[pivot_row][column]
        matrix[pivot_row] = [value / scale for value in matrix[pivot_row]]
        for row in range(height):
            if row == pivot_row or not matrix[row][column]:
                continue
            factor = matrix[row][column]
            matrix[row] = [
                left - factor * right
                for left, right in zip(matrix[row], matrix[pivot_row], strict=True)
            ]
        pivot_row += 1
        if pivot_row == height:
            break
    return pivot_row


def main() -> None:
    thm741 = load_thm741()

    all_bodies = tuple(combinations(range(1, 15), 9))
    floods = tuple(E for E in all_bodies if not uncovered_obligations(E))
    expected_floods = tuple(
        tuple(sorted((*range(8, 15), a, b)))
        for a, b in combinations(range(1, 8), 2)
    )
    assert set(floods) == set(expected_floods)
    assert len(floods) == 21

    lines = fano_lines()
    assert len(lines) == 7
    assert all(len(line) == 3 and line[0] ^ line[1] ^ line[2] == 0 for line in lines)
    edge_to_line = {}
    for line in lines:
        for edge in combinations(line, 2):
            assert frozenset(edge) not in edge_to_line
            edge_to_line[frozenset(edge)] = line
    assert len(edge_to_line) == 21

    # Inclusion-exclusion on nonempty B_3 masks: odd support +, even support -.
    boolean_sign = {mask: 1 if mask.bit_count() % 2 else -1 for mask in range(1, 8)}
    chi7_extended = {
        mask: (1 if mask == 7 or mask in (1, 2, 4) else -1)
        for mask in range(1, 8)
    }
    assert boolean_sign == chi7_extended
    negative_line = tuple(mask for mask in range(1, 8) if boolean_sign[mask] < 0)
    positive_cap = tuple(mask for mask in range(1, 8) if boolean_sign[mask] > 0)
    assert negative_line == (3, 5, 6)
    assert negative_line in lines
    assert positive_cap == (1, 2, 4, 7)

    # Exact THM-741 body data, with one edge/flood-body per row.
    rows = []
    weights = {}
    for a, b in combinations(range(1, 8), 2):
        body = tuple(sorted((*range(8, 15), a, b)))
        _, r_E, m_E = thm741.good_norm(body)
        threshold = 3 * m_E / (thm741.S2 * r_E)
        V1 = thm741.minV(4, threshold.numerator, threshold.denominator)
        line = edge_to_line[frozenset((a, b))]
        weight = (r_E, m_E, V1)
        weights[frozenset((a, b))] = weight
        rows.append(
            {
                "edge": [a, b],
                "fano_line": list(line),
                "r": r_E,
                "m": str(m_E),
                "V1": V1,
            }
        )

    maps = gl32_maps()
    assert len(maps) == 168
    stabilizers = {}
    for coordinate, name in enumerate(("r", "m", "V1")):
        preserving = []
        for point_map in maps:
            if all(
                weights[edge][coordinate]
                == weights[
                    frozenset(point_map[p] for p in edge)
                ][coordinate]
                for edge in weights
            ):
                preserving.append((point_map[1], point_map[2], point_map[4]))
        stabilizers[name] = preserving
        assert preserving == [(1, 2, 4)]

    # The three exact observables live in different K7 edge modules.
    edges = tuple(frozenset(edge) for edge in combinations(range(1, 8), 2))
    point_rows = [
        [int(vertex in edge) for edge in edges]
        for vertex in range(1, 8)
    ]
    fano_rows = [
        [int(edge.issubset(set(line))) for edge in edges]
        for line in lines
    ]
    point_rank = rational_rank(point_rows)
    fano_rank = rational_rank(fano_rows)
    point_fano_rank = rational_rank(point_rows + fano_rows)
    assert (point_rank, fano_rank, point_fano_rank) == (7, 7, 13)

    orbit_ranks = {}
    centered_orbit_ranks = {}
    in_point_fano_span = {}
    for coordinate, name in enumerate(("r", "m", "V1")):
        base = [weights[edge][coordinate] for edge in edges]
        orbit = [
            [
                weights[frozenset(point_map[point] for point in edge)][coordinate]
                for edge in edges
            ]
            for point_map in maps
        ]
        centered = [
            [left - right for left, right in zip(row, base, strict=True)]
            for row in orbit
        ]
        orbit_ranks[name] = rational_rank(orbit)
        centered_orbit_ranks[name] = rational_rank(centered)
        in_point_fano_span[name] = (
            rational_rank(point_rows + fano_rows + [base]) == point_fano_rank
        )
    assert orbit_ranks == {"r": 7, "m": 21, "V1": 21}
    assert centered_orbit_ranks == {"r": 6, "m": 20, "V1": 20}
    assert in_point_fano_span == {"r": True, "m": False, "V1": False}

    point_potential = {1: 16, 2: 16, 3: 14, 4: 14, 5: 10, 6: 14, 7: 6}
    assert all(
        weights[edge][0] == sum(point_potential[point] for point in edge)
        for edge in edges
    )
    m_curl = (
        weights[frozenset((1, 2))][1]
        + weights[frozenset((6, 7))][1]
        - weights[frozenset((1, 6))][1]
        - weights[frozenset((2, 7))][1]
    )
    v1_curl = (
        weights[frozenset((1, 3))][2]
        + weights[frozenset((5, 6))][2]
        - weights[frozenset((1, 5))][2]
        - weights[frozenset((3, 6))][2]
    )
    assert m_curl == Fraction(1, 21)
    assert v1_curl == 48

    # Tournament Analysis on proof carriers (the 21 flood bodies/K7 edges).
    # All V1 values are distinct, so the gauge is a transitive total order.
    assert len({row["V1"] for row in rows}) == len(rows)
    v1_hamiltonian_path = tuple(
        tuple(row["edge"]) for row in sorted(rows, key=lambda row: row["V1"])
    )
    lex_gauge_flips = sum(
        1
        for left_index, left in enumerate(rows)
        for right in rows[left_index + 1 :]
        if left["V1"] > right["V1"]
    )
    tournament = {
        "vertices": len(rows),
        "score_histogram": {score: 1 for score in range(len(rows))},
        "directed_3cycles": 0,
        "scc_sizes": [1] * len(rows),
        "edge_flips_vs_lex_gauge": lex_gauge_flips,
        "hamiltonian_path_count": 1,
        "hamiltonian_path": v1_hamiltonian_path,
    }
    assert lex_gauge_flips == 187

    # Each endpoint needle toggles exactly once.  Weak orders include all
    # possible coincidences, and constant needles only delete a switch/event.
    negative_set = set(negative_line)
    paths = {
        needle_path(initial, blocks)
        for initial in range(8)
        for blocks in ordered_event_blocks()
    }
    max_masks = max(len(set(path)) for path in paths)
    max_negative = max(len(set(path) & negative_set) for path in paths)
    assert max_masks == 4
    assert max_negative == 2
    assert all(not set(range(1, 8)).issubset(path) for path in paths)
    assert all(not negative_set.issubset(path) for path in paths)

    # The modular-obligation carrier vanishes exactly on the flood tail.
    assert all(uncovered_obligations(E) == () for E in floods)
    assert thm741.lcm(()) == 1

    payload = {
        "flood_count": len(floods),
        "flood_form": "{8,...,14} union one edge of K7",
        "fano_lines": lines,
        "chi7_word": "++-+--+",
        "negative_line": negative_line,
        "positive_cap": positive_cap,
        "rows": rows,
        "gl32_count": len(maps),
        "stabilizer_sizes": {name: len(value) for name, value in stabilizers.items()},
        "edge_module": {
            "orbit_ranks": orbit_ranks,
            "centered_orbit_ranks": centered_orbit_ranks,
            "point_rank": point_rank,
            "fano_rank": fano_rank,
            "point_plus_fano_rank": point_fano_rank,
            "invisible_dimension": len(edges) - point_fano_rank,
            "in_point_plus_fano_span": in_point_fano_span,
            "r_point_potential": point_potential,
            "m_curl": str(m_curl),
            "V1_curl": v1_curl,
        },
        "tournament": tournament,
        "weak_needle_paths": len(paths),
        "max_masks": max_masks,
        "max_negative_line_points": max_negative,
        "flood_obligation_count": 0,
        "flood_lcm_filter": thm741.lcm(()),
    }
    payload_bytes = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()

    print("FANO/CHI7/J=4 FLOOD/FAREY-NEEDLE EXACT AUDIT")
    print("method: exact Fraction interval unions + exhaustive finite incidence/group audit")
    print(f"bodies audited: {len(all_bodies)} ; flood bodies: {len(floods)}")
    print("flood characterization: {8,...,14} union {a,b}, 1<=a<b<=7")
    print(f"fano lines: {lines}")
    print("edge partition: 21 K7 edges = 7 Fano lines x 3 edges")
    print(f"B3/chi7 signs: {''.join('+' if boolean_sign[i] > 0 else '-' for i in range(1, 8))}")
    print(f"negative Fano line: {negative_line} ; positive complement/cap: {positive_cap}")
    print("flood rows (edge, line, r(E), m(E), V1(E)):")
    for row in rows:
        print(
            "  %s  line=%s  r=%s  m=%s  V1=%s"
            % (tuple(row["edge"]), tuple(row["fano_line"]), row["r"], row["m"], row["V1"])
        )
    print(f"GL(3,2) maps audited: {len(maps)}")
    print(
        "Fano-automorphism stabilizers: "
        + ", ".join(f"{name}={len(value)} (identity only)" for name, value in stabilizers.items())
    )
    print(
        f"edge-module orbit ranks r/m/V1={orbit_ranks['r']}/{orbit_ranks['m']}/{orbit_ranks['V1']} ; "
        f"centered={centered_orbit_ranks['r']}/{centered_orbit_ranks['m']}/{centered_orbit_ranks['V1']}"
    )
    print(
        f"point/Fano/combined incidence ranks={point_rank}/{fano_rank}/{point_fano_rank}; "
        f"invisible edge dimension={len(edges)-point_fano_rank}"
    )
    print(f"r edge potential x={tuple(point_potential.values())}; r_ab=x_a+x_b")
    print(
        f"m,V1 outside point+Fano span; curls m12+m67-m16-m27={m_curl}, "
        f"V1(13)+V1(56)-V1(15)-V1(36)={v1_curl}"
    )
    print("Tournament Analysis vertices: 21 flood bodies/K7 edges (not runners)")
    print("pair observable: V1(F)-V1(E); switch: E->F iff V1(E)<V1(F); tie: lexicographic")
    print(
        "fingerprint: score_hist={0..20:1}, directed_3cycles=0, SCC_sizes=21x1, "
        f"edge_flips_vs_lex={lex_gauge_flips}, Hamiltonian_paths=1"
    )
    print(f"tie Hamiltonian path (increasing V1): {v1_hamiltonian_path}")
    print(
        f"three endpoint needles: {len(paths)} distinct weak-order mask paths; "
        f"max masks={max_masks}; max negative-line points={max_negative}"
    )
    print("flood modular carrier: Q(E)=empty, Qb=empty, lcm(Qb)=1 on all 21 bodies")
    print("scope: no Fano symmetry quotient; no local 3-needle realization of the full chi7 carrier")
    print("scope: this diagnoses but does not close the exact-geometry flood sweeps")
    print("quotients: body keeps (r,m,V1); obligation keeps lcm filter; needle keeps local switch order")
    print("destroyed: lower-tree geometry; exact intervals; speed identity/cross-gap gluing, respectively")
    print(f"source_sha256={sha256(Path(__file__).resolve())}")
    print(f"thm741_dependency_sha256={sha256(THM741_PATH)}")
    print(f"payload_sha256={hashlib.sha256(payload_bytes).hexdigest()}")
    print("VERDICT: ALL EXACT ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
