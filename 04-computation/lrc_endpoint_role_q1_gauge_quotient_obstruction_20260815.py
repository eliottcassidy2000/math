#!/usr/bin/env python3
"""Exact q1-gauge obstruction for the THM-2334 -> role-carrier transplant.

This standalone sidecar types the cheapest coordinate-labelled use of the
THM-2334 target aggregate.  For a scalar word w and retained coordinate i,
the q1-gauged residue

    r_i = e_i - (w_i / w_q1) e_q1

lies in K_w = ker(w mod 13).  The only direct way to feed that coordinate to
the published target response A_w is

    P_i = A_w([r_i] in K_w/L_w).

For both THM-3479 tuples, [r_H] = [r_q5].  The role contract puts H and q5 on
the unique bridge of the two-K4 carrier, so every edge-gradient determinant
built from P has zero bridge and vanishes identically.  This is a quotient
obstruction, not a computation of a physical LRC current.

Among the six named single-row deletions, deleting the H-labelled owner-packet
row is uniquely able to separate H from q5.  Its dual adds the single
character e_H, has 13^3=2197 twists, and is dimension-minimal.  This does not
delete the Boolean guard factor.  The script verifies the refined class
geometry, and pins the later exact U_full bank realizing all 72 graph-factor
products.  U_clock and physical/common-ancestry realization remain open.
"""

from __future__ import annotations

import ast
import hashlib
from itertools import permutations
from pathlib import Path


P = 13
ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "04-computation/lrc_endpoint_role_q1_gauge_quotient_obstruction_20260815.py"
OUTPUT = "05-knowledge/results/lrc_endpoint_role_q1_gauge_quotient_obstruction_20260815.out"

PINS = (
    (
        "THM-2334-PROVED",
        ROOT / "01-canon/theorems/THM-2334-relation-residue-current-and-character-twist-pushforward.md",
        "c80f7bb3d31274a02f046fa6cea3b36cd56c62be611936ca32ee723881cd3899",
        "PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED",
    ),
    (
        "THM-3479-PROVED",
        ROOT / "01-canon/theorems/THM-3479-literal-half-twist-relation-current-two-transplant-certificate.md",
        "025998551e3cdf3c6e4db5c0a4f208dd32f6845970fd4729d4a276035e0fdfeb",
        "PROVED STRUCTURAL + VERIFIED-EXACT + FINITE-EXACT + INDEPENDENTLY AUDITED",
    ),
    (
        "THM-3479-COMPANION",
        ROOT / "04-computation/lrc_half_twist_relation_current_bridge_thm3479.py",
        "ad2a620cdc238f28e3384698b2c612f38cdf2566bd56b76d1cbabcc03107ec0b",
        None,
    ),
    (
        "ROLE-CHART-SIDECAR",
        ROOT / "04-computation/lrc_relation_role_chart_weighted_closure_probe_20260815.py",
        "207c65ca235ea5647e346027d424264e8abbcf27c5f574b5901cca13611d7e03",
        None,
    ),
    (
        "U_FULL-REFINED-BANK-SCRIPT",
        ROOT / "04-computation/lrc14_guard_deleted_refined_endpoint_role_probe_20260816.py",
        "ee2105742abee578a9c41ff7ec954a07ada324fccc2c643429e7ac6e6e6f8fc2",
        None,
    ),
    (
        "U_FULL-REFINED-BANK-OUTPUT",
        ROOT / "05-knowledge/results/lrc14_guard_deleted_refined_endpoint_role_probe_20260816.out",
        "10a98351cc59615a5b6d2b8f555e0936d1a39566d9906127edc2b0fbc3918e73",
        "STATUS=PASS",
    ),
)

CUR_LABELS = ("H", "q1", "q2", "q3", "q4", "q5", "c1", "c2", "c3")
CUR_INDEX = {label: index for index, label in enumerate(CUR_LABELS)}
ROLE_LABELS = ("c1", "c2", "c3", "H", "q2", "q3", "q4", "q5")

RELATION_REL_ORDER = (-27, -27, -27, 20110798, -41, -27, -27, -27, 38)
RELATION_CURRENT = (20110798, -41, -27, -27, -27, 38, -27, -27, -27)

U_FULL = (1, 183, 27, 131, 53, 313, 13, 2197, 742586)
U_CLOCK = (5, 661549, 655231, 658533, 661445, 291, 65, 2197, 742586)
WORDS = (("U_full", U_FULL), ("U_clock", U_CLOCK))

OWNER_PACKET_LABELS = ("c1", "H", "q1", "q2", "q3", "q4")
OMITTED_UNIT = CUR_INDEX["q5"]
TARGET_A = CUR_INDEX["c2"]
TARGET_B = CUR_INDEX["c3"]

# These span the published two-dimensional target-character group modulo <w>.
V1 = (0, -1, 0, 0, 0, 0, 0, 1, 0)  # -e_q1 + e_c2
V2 = (0, 0, -1, 0, 0, 0, 0, 0, 1)  # -e_q2 + e_c3
E_H = (1, 0, 0, 0, 0, 0, 0, 0, 0)

# Vertices u1,...,u8 are encoded by 0,...,7.
EDGES = (
    (0, 3), (0, 4), (0, 5),
    (1, 2), (1, 4), (1, 7),
    (2, 4), (2, 7),
    (3, 4), (3, 5),
    (4, 5), (4, 6), (4, 7),
)
HUB = 4
LEAF = 6
BRIDGE = (HUB, LEAF)
WINGS = ((0, 3, 5), (1, 2, 7))
BLOCKERS = ("c1", "c2", "c3")
MIDDLE_UNITS = ("q2", "q3", "q4")


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    data = path.read_bytes().replace(b"\r\n", b"\n")
    return hashlib.sha256(data).hexdigest()


def dot(left: tuple[int, ...], right: tuple[int, ...]) -> int:
    return sum(a * b for a, b in zip(left, right))


def mod_vector(vector: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(value % P for value in vector)


def rref_mod(rows: tuple[tuple[int, ...], ...]) -> tuple[
    tuple[tuple[int, ...], ...], tuple[int, ...]
]:
    if not rows:
        return (), ()
    width = len(rows[0])
    matrix = [list(mod_vector(row)) for row in rows if any(value % P for value in row)]
    pivots: list[int] = []
    pivot_row = 0
    for column in range(width):
        source = next(
            (row for row in range(pivot_row, len(matrix)) if matrix[row][column]),
            None,
        )
        if source is None:
            continue
        matrix[pivot_row], matrix[source] = matrix[source], matrix[pivot_row]
        inverse = pow(matrix[pivot_row][column], -1, P)
        matrix[pivot_row] = [value * inverse % P for value in matrix[pivot_row]]
        for row in range(len(matrix)):
            if row == pivot_row:
                continue
            factor = matrix[row][column]
            if factor:
                matrix[row] = [
                    (left - factor * right) % P
                    for left, right in zip(matrix[row], matrix[pivot_row])
                ]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    nonzero = tuple(tuple(row) for row in matrix if any(row))
    return nonzero, tuple(pivots)


def rank_mod(rows: tuple[tuple[int, ...], ...]) -> int:
    return len(rref_mod(rows)[0])


def owner_packet_rows(word: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    rows = []
    for label in OWNER_PACKET_LABELS:
        index = CUR_INDEX[label]
        row = [0] * 9
        row[OMITTED_UNIT] += word[index]
        row[index] -= word[OMITTED_UNIT]
        if label == "q1":
            row[OMITTED_UNIT] += word[TARGET_A]
            row[TARGET_A] -= word[OMITTED_UNIT]
        if label == "q2":
            row[OMITTED_UNIT] += word[TARGET_B]
            row[TARGET_B] -= word[OMITTED_UNIT]
        rows.append(mod_vector(tuple(row)))
    return tuple(rows)


def q1_gauged_residue(word: tuple[int, ...], label: str) -> tuple[int, ...]:
    q1 = CUR_INDEX["q1"]
    index = CUR_INDEX[label]
    ratio = word[index] * pow(word[q1], -1, P) % P
    residue = [0] * 9
    residue[index] = 1
    residue[q1] = (residue[q1] - ratio) % P
    answer = tuple(residue)
    require(dot(answer, word) % P == 0, (label, "q1 gauge leaves K", answer))
    return answer


def quotient_coordinate(
    residue: tuple[int, ...], dual_basis: tuple[tuple[int, ...], ...]
) -> tuple[int, ...]:
    return tuple(dot(character, residue) % P for character in dual_basis)


def role_charts() -> tuple[tuple[str, ...], ...]:
    charts = []
    for swap in (0, 1):
        blocker_wing = WINGS[swap]
        unit_wing = WINGS[1 - swap]
        for blocker_order in permutations(BLOCKERS):
            for unit_order in permutations(MIDDLE_UNITS):
                chart = {HUB: "H", LEAF: "q5"}
                chart.update(zip(blocker_wing, blocker_order))
                chart.update(zip(unit_wing, unit_order))
                charts.append(tuple(chart[vertex] for vertex in range(8)))
    answer = tuple(sorted(set(charts)))
    require(len(answer) == 72, ("role chart count", len(answer)))
    return answer


def component_count(edges: tuple[tuple[int, int], ...]) -> int:
    adjacency = {vertex: set() for vertex in range(8)}
    for left, right in edges:
        adjacency[left].add(right)
        adjacency[right].add(left)
    unseen = set(range(8))
    components = 0
    while unseen:
        components += 1
        stack = [min(unseen)]
        unseen.remove(stack[0])
        while stack:
            vertex = stack.pop()
            for neighbor in adjacency[vertex]:
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    stack.append(neighbor)
    return components


def determinant(matrix: list[list[int]]) -> int:
    if not matrix:
        return 1
    total = 0
    for column, value in enumerate(matrix[0]):
        if not value:
            continue
        minor = [row[:column] + row[column + 1 :] for row in matrix[1:]]
        total += (-1 if column % 2 else 1) * value * determinant(minor)
    return total


def weighted_tree_determinant(
    chart: tuple[str, ...], responses: dict[str, int]
) -> int:
    potentials = tuple(responses[label] for label in chart)
    laplacian = [[0] * 8 for _ in range(8)]
    for left, right in EDGES:
        weight = potentials[left] - potentials[right]
        laplacian[left][left] += weight
        laplacian[right][right] += weight
        laplacian[left][right] -= weight
        laplacian[right][left] -= weight
    reduced = [row[:-1] for row in laplacian[:-1]]
    return determinant(reduced)


def synthetic_response(
    classes: dict[str, tuple[int, ...]], base: int
) -> dict[str, int]:
    distinct = tuple(sorted(set(classes.values())))
    values = {coordinate: base**index for index, coordinate in enumerate(distinct)}
    return {label: values[coordinate] for label, coordinate in classes.items()}


def security_gate(source: Path) -> tuple[int, tuple[str, ...]]:
    tree = ast.parse(source.read_text(encoding="utf-8"))
    forbidden_nodes = (ast.Assert,)
    forbidden_calls = {"eval", "exec", "compile", "__import__"}
    bad = []
    for node in ast.walk(tree):
        if isinstance(node, forbidden_nodes):
            bad.append(type(node).__name__)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
            if node.func.id in forbidden_calls:
                bad.append(node.func.id)
    require(not bad, ("security gate", bad))
    return len(tuple(ast.walk(tree))), tuple(bad)


def main() -> None:
    pin_rows = []
    for label, path, expected_hash, status_text in PINS:
        actual_hash = lf_sha256(path)
        require(actual_hash == expected_hash, (label, "hash drift", actual_hash))
        if status_text is not None:
            text = path.read_text(encoding="utf-8")
            require(status_text in text, (label, "status drift"))
        pin_rows.append((label, actual_hash))

    require(RELATION_REL_ORDER == (-27, -27, -27, 20110798, -41, -27, -27, -27, 38),
            "relation-order drift")
    charts = role_charts()
    require(component_count(EDGES) == 1, "carrier is disconnected")
    graph_bridges = tuple(
        edge for index, edge in enumerate(EDGES)
        if component_count(EDGES[:index] + EDGES[index + 1 :]) == 2
    )
    require(graph_bridges == (BRIDGE,), ("carrier bridge census", graph_bridges))

    packets = []
    coarse_profiles = []
    refined_profiles = []
    fixed_relation_classes = []
    removed_row_census = []
    synthetic_census = []

    for name, word in WORDS:
        require(dot(RELATION_CURRENT, word) == 0, (name, "fixed relation"))
        packet = owner_packet_rows(word)
        require(rank_mod(packet) == 6, (name, "packet rank"))
        require(all(dot(row, word) % P == 0 for row in packet),
                (name, "packet outside K"))
        require(all(dot(row, V1) % P == 0 and dot(row, V2) % P == 0
                    for row in packet), (name, "coarse dual basis"))
        require(rank_mod((mod_vector(word), mod_vector(V1), mod_vector(V2))) == 3,
                (name, "coarse dual independence"))

        residues = {label: q1_gauged_residue(word, label) for label in ROLE_LABELS}
        bridge_difference = tuple(
            (left - right) % P
            for left, right in zip(residues["H"], residues["q5"])
        )
        guard_row = packet[OWNER_PACKET_LABELS.index("H")]
        guard_multiplier = (-word[CUR_INDEX["q5"]]) % P
        require(guard_row == tuple(guard_multiplier * value % P
                                   for value in bridge_difference),
                (name, "guard row is not the bridge difference"))
        coarse = {
            label: quotient_coordinate(residue, (V1, V2))
            for label, residue in residues.items()
        }
        require(coarse["H"] == coarse["q5"], (name, "bridge did not collapse"))
        require(coarse == {
            "c1": (0, 0), "c2": (1, 0), "c3": (0, 1),
            "H": (1, 0), "q2": (1, 12), "q3": (1, 0),
            "q4": (1, 0), "q5": (1, 0),
        }, (name, "coarse profile", coarse))

        address_mod = mod_vector(RELATION_CURRENT)
        require(dot(address_mod, word) % P == 0, (name, "address outside K"))
        relation_class = quotient_coordinate(address_mod, (V1, V2))
        require(relation_class == (1, 0), (name, "relation target class", relation_class))

        packet_without_guard = tuple(
            row for label, row in zip(OWNER_PACKET_LABELS, packet) if label != "H"
        )
        require(rank_mod(packet_without_guard) == 5, (name, "refined packet rank"))
        require(all(dot(row, E_H) % P == 0 for row in packet_without_guard),
                (name, "e_H does not annihilate refined packet"))
        require(rank_mod((mod_vector(word), mod_vector(V1), mod_vector(V2), E_H)) == 4,
                (name, "refined dual independence"))
        refined = {
            label: quotient_coordinate(residue, (V1, V2, E_H))
            for label, residue in residues.items()
        }
        require(refined["H"] != refined["q5"],
                (name, "minimal refinement still collapses bridge"))
        require(refined == {
            "c1": (0, 0, 0), "c2": (1, 0, 0), "c3": (0, 1, 0),
            "H": (1, 0, 1), "q2": (1, 12, 0), "q3": (1, 0, 0),
            "q4": (1, 0, 0), "q5": (1, 0, 0),
        }, (name, "refined profile", refined))

        removal_results = []
        for removed_index, removed_label in enumerate(OWNER_PACKET_LABELS):
            subpacket = packet[:removed_index] + packet[removed_index + 1 :]
            require(rank_mod(subpacket) == 5, (name, "subpacket rank", removed_label))
            separated = rank_mod(subpacket + (bridge_difference,)) == 6
            removal_results.append((removed_label, separated))
        require(removal_results == [
            ("c1", False), ("H", True), ("q1", False),
            ("q2", False), ("q3", False), ("q4", False),
        ], (name, "minimal row removal census", removal_results))

        coarse_response = synthetic_response(coarse, 2)
        coarse_determinants = tuple(
            weighted_tree_determinant(chart, coarse_response) for chart in charts
        )
        require(set(coarse_determinants) == {0},
                (name, "coarse determinant escaped zero"))

        refined_control = None
        for base in range(2, 101):
            response = synthetic_response(refined, base)
            determinants = tuple(
                weighted_tree_determinant(chart, response) for chart in charts
            )
            if all(determinants):
                refined_control = (
                    base,
                    min(abs(value) for value in determinants),
                    len(set(determinants)),
                    sum(value > 0 for value in determinants),
                    sum(value < 0 for value in determinants),
                )
                break
        require(refined_control is not None, (name, "no refined positive control"))

        packets.append((name, rank_mod(packet), rank_mod(packet_without_guard)))
        coarse_profiles.append((name, tuple((label, coarse[label]) for label in ROLE_LABELS)))
        refined_profiles.append((name, tuple((label, refined[label]) for label in ROLE_LABELS)))
        fixed_relation_classes.append((name, relation_class,
                                       quotient_coordinate(address_mod, (V1, V2, E_H))))
        removed_row_census.append((name, tuple(removal_results)))
        synthetic_census.append((name, len(coarse_determinants), len(charts),
                                 refined_control))

    require(coarse_profiles[0][1] == coarse_profiles[1][1],
            "tuple-dependent coarse profile")
    require(refined_profiles[0][1] == refined_profiles[1][1],
            "tuple-dependent refined profile")

    consequence = (
        tuple(packets), tuple(coarse_profiles), tuple(refined_profiles),
        tuple(fixed_relation_classes), tuple(removed_row_census),
        tuple(synthetic_census), len(charts), BRIDGE, 13**2, 13**3,
    )
    semantic_hash = hashlib.sha256(repr(consequence).encode("utf-8")).hexdigest()
    security = security_gate(ROOT / SCRIPT)

    print("THM-2334 endpoint-to-role q1-gauge quotient obstruction")
    print(f"script={SCRIPT}")
    print(f"output={OUTPUT}")
    print(f"pins={tuple(pin_rows)}")
    print("STATUS=FINITE-EXACT QUOTIENT/BACK-MAP COMPONENT OF PROVED THM-3479; INDEPENDENTLY AUDITED")
    print("GAMMA_TYPE=gamma_w: Ghat_w=L_w^perp/<w> -> cyclotomic endpoint-current scalars")
    print("A_TYPE=A_w: G_w=K_w/L_w -> cyclotomic unrestricted relation-residue aggregates; A is the normalized inverse DFT of gamma")
    print("RELATION_ADDRESS_TYPE=a in Lambda(w)=ker(w:Z^9->Z); C(a;X,m) is one grouped exact-address coefficient and is not computed here")
    print("COORDINATE_RESPONSE_CANDIDATE=P_i=A_w([e_i-(w_i/w_q1)e_q1]); gamma(e_i) would be ill-typed")
    print(f"packet_ranks=(full,remove_guard)={tuple(packets)}")
    print(f"coarse_coordinate_profiles={tuple(coarse_profiles)}")
    print(f"fixed_relation_classes=(coarse,refined)={tuple(fixed_relation_classes)}")
    print(f"role_charts={len(charts)} forced_bridge={BRIDGE} labels=(H,q5)")
    print("COARSE_VERDICT=[r_H]=[r_q5], hence P_H=P_q5 for every A_w and all 72 weighted tree determinants vanish identically")
    print(f"coarse_synthetic_zero_determinants={tuple((name, zero, total) for name, zero, total, _ in synthetic_census)}")
    print(f"single_row_removal_separation_census={tuple(removed_row_census)}")
    print("MINIMAL_REFINEMENT=among the six named one-row deletions, only deleting the H-labelled owner-packet row separates the bridge; the Boolean guard_safe factors remain; add character e_H; quotient dimension 3; twist bank size 2197")
    print(f"refined_coordinate_profiles={tuple(refined_profiles)}")
    print(f"refined_synthetic_positive_controls={tuple((name, control) for name, _, _, control in synthetic_census)}")
    print("COMPLETED_U_FULL_REFINEMENT=the pinned 2197-twist bank, normalized inverse DFT at five distinct role classes, coarse fibre-sum back-map, and graph census give nonzero bridge and both K4 factors in all 72 labelled charts")
    print("MISSING_EXACT_DATA=the U_clock refined bank and a proof that the U_full edge differences are a lawful physical common-ancestry current")
    print("PHYSICAL_CURRENT_TYPE=an actual Boolean/common-ancestry endpoint or THM-2512 response observable; neither gamma, A, nor a synthetic graph gradient supplies it")
    print("NONCONSEQUENCES=no grouped C(a;X,m) nonvanishing, all-91-unit B(q), ancestry, bispectrum, scalar-row exclusion, or LRC(14)")
    print(f"semantic_sha256={semantic_hash}")
    print(f"security_ast_nodes_and_forbidden={security}")
    print("VERDICT=the canonical coarse THM-2334 target response cannot power the role determinant; one extra guard character gives the dimension-minimal named repair, realized exactly for U_full but not upgraded to a physical current")


if __name__ == "__main__":
    main()
