#!/usr/bin/env python3
"""Derived exact P7/arc audit for the r=5 pointed-six carrier.

This probe deliberately avoids rebuilding the five-coordinate tensor of
``b1baa781a``.  It does three smaller things.

1. It reconstructs the pointed relation table through the independently
   audited pointed-root implementation and identifies its six independent
   rows with both the edges of the alternating incidence path P7 and the arcs
   of the bidirected owner-state path A<->B<->C<->D.
2. It reopens only the exact two-digit source profiles.  This proves the
   arc-reversal-even rank-three coefficient space away from multiples of
   three, the even-plus-middle-current rank-four space at r0=6, and the
   exceptional tent-section/source-quotient ledger.
3. It combines those source facts with the already pinned diagonal-bundle
   theorem and the independently audited conditional ranks.  No output tensor
   is recomputed here; the r0=3,9 completion is a derived consequence of those
   exact parents.

Everything is static finite response data.  In particular, no chronology,
local system, cocycle, physical current, row exclusion, or LRC(14) claim is
made.
"""

from __future__ import annotations

from hashlib import sha256
import importlib
import json
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

POINT = importlib.import_module(
    "lrc_r5_ufull_owner_node_pointed_six_state_root_difference_"
    "independent_audit_20260816"
)
BUNDLE = importlib.import_module(
    "lrc_r5_two_current_digits_pointed_root_difference_carrier_"
    "transition_probe_20260816"
)
TWO = importlib.import_module(
    "lrc_r5_ufull_owner_node_boolean_square_two_digit_current_"
    "ancestry_independent_audit_20260816"
)

P = 13
PRIME = POINT.PRIME
POINTS = POINT.POINTS
STATE_PATH = (0, 1, 3, 2)
STATE_NAMES = ("A", "B", "C", "D")
ROOTS = (0, 6, 12)
P7_VERTICES = (
    ("state", 0), ("root", 0), ("state", 1), ("root", 6),
    ("state", 3), ("root", 12), ("state", 2),
)
ARC_REVERSAL = (1, 0, 3, 2, 5, 4)
PATH_REFLECTION = (5, 4, 3, 2, 1, 0)
NONROOT_DIGITS = (1, 2, 3, 4, 5, 7, 8, 9, 10, 11)
TENT_SECTION = (12, 12, 9, 3, 0, 0)

EXPECTED_COMMON_RELATION_ROWSPACE = (
    "6e9083f15408f6d2d85fb3f2747ba0bd1f987e83ce4b836cb7298aaccc84e0c4"
)
EXPECTED_BOOLEAN_PARENT_ROWSPACE = (
    "1b4fef00a23dcb79dc52ace662bae2f858ce3da27b6ef19b902ae40f5a79e755"
)
EXPECTED_RECORD_SHA256 = (
    "c1e17b60d62fdafe6a2b99fd5df6d097db1d327590ab0995bda182466d609ec4"
)
EXPECTED_SEMANTIC_SHA256 = (
    "54f2fd0702cfee4b3af954ca3719b9966cb06e9227a94bdc8b14651c452d5feb"
)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def digest_json(value: object) -> str:
    return sha256(
        json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    ).hexdigest()


def rank(rows) -> int:
    return POINT.R.rank_mod(
        tuple(tuple(value % PRIME for value in row) for row in rows)
    )


def canonical_row_basis(rows):
    matrix = [list(value % PRIME for value in row) for row in rows]
    if not matrix:
        return ()
    columns = len(matrix[0])
    require(all(len(row) == columns for row in matrix), "ragged matrix")
    pivot_row = 0
    for column in range(columns):
        pivot = next(
            (row for row in range(pivot_row, len(matrix)) if matrix[row][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        inverse = pow(matrix[pivot_row][column], -1, PRIME)
        matrix[pivot_row] = [value * inverse % PRIME for value in matrix[pivot_row]]
        for row in range(len(matrix)):
            if row == pivot_row or not matrix[row][column]:
                continue
            factor = matrix[row][column]
            matrix[row] = [
                (left - factor * right) % PRIME
                for left, right in zip(matrix[row], matrix[pivot_row])
            ]
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    return tuple(tuple(row) for row in matrix[:pivot_row])


def rowspace_digest(rows) -> str:
    return digest_json(canonical_row_basis(rows))


def standard_basis(size: int, index: int):
    return tuple(int(column == index) for column in range(size))


def add(left, right):
    return tuple((a + b) % PRIME for a, b in zip(left, right))


def subtract(left, right):
    return tuple((a - b) % PRIME for a, b in zip(left, right))


def linear_image(coefficients, rows):
    return tuple(
        sum(coefficient * row[column]
            for coefficient, row in zip(coefficients, rows)) % PRIME
        for column in range(len(rows[0]))
    )


def module_split(rows, permutation):
    reflected = tuple(rows[index] for index in permutation)
    plus = tuple(add(row, image) for row, image in zip(rows, reflected))
    minus = tuple(subtract(row, image) for row, image in zip(rows, reflected))
    return rank(rows), rank(plus), rank(minus)


def reconstruct_pointed_tensor():
    """Replay the clean-room pointed construction, not the four-way candidate."""

    R = POINT.R
    R.split_field_certificate()
    source_u, source_v, source_boundaries, profile_digest, _total, _types = (
        R.source_profiles()
    )
    source_grid = R.SRC.T_DEN
    hierarchy = POINT.SQ.source_hierarchy(
        source_u, source_v, source_boundaries, source_grid
    )
    require(hierarchy is not None, "source hierarchy")
    word, endpoint_grid = R.endpoint_word_and_grid()
    source_u_scaled = tuple(
        R.scale_profile(profile, R.JOINT_COORDINATE, source_grid)
        for profile in source_u
    )
    source_v_scaled = tuple(
        R.scale_profile(profile, R.JOINT_COORDINATE, source_grid)
        for profile in source_v
    )
    boundaries = R.fixed_boundaries(source_boundaries, source_grid)
    harmonic = R.HarmonicPrimitive(word, endpoint_grid)
    danger = R.danger_arcs()
    gamma, work_counts, point_counts = POINT.build_banks(
        word, endpoint_grid, source_u_scaled, source_v_scaled,
        boundaries, harmonic, danger,
    )
    require(work_counts == POINT.EXPECTED_WORK_COUNTS,
            ("pointed work counts", work_counts))
    require(point_counts == POINT.EXPECTED_POINT_SEGMENTS,
            ("pointed segment counts", point_counts))
    require(tuple(digest_json(bank) for bank in gamma)
            == POINT.EXPECTED_DIGESTS[0], "pointed gamma drift")
    zeta = pow(R.JOINT_ROOT, R.JOINT_ORDER // P, PRIME)
    tensor = POINT.inverse_tensor(gamma[0], zeta)
    parent = POINT.parent_marginal(tensor)
    require(digest_json(tensor) == POINT.EXPECTED_DIGESTS[1][0],
            "pointed tensor drift")
    require(digest_json(parent) == POINT.EXPECTED_PARENT_TENSOR[0],
            "pointed parent drift")
    return tensor, parent, profile_digest


def p7_arc_record(relation_rows, parent):
    p7_edges = tuple(
        (P7_VERTICES[index], P7_VERTICES[index + 1]) for index in range(6)
    )
    require(tuple((edge[0][1], edge[1][1])
                  if edge[0][0] == "state" else (edge[1][1], edge[0][1])
                  for edge in p7_edges) == POINTS,
            ("P7 edge labels", p7_edges, POINTS))

    root_neighbours = {0: (0, 1), 6: (1, 3), 12: (3, 2)}
    arcs = []
    for state, root in POINTS:
        neighbours = root_neighbours[root]
        require(state in neighbours, ("root incidence", state, root, neighbours))
        other = neighbours[1] if state == neighbours[0] else neighbours[0]
        arcs.append((state, other))
    arcs = tuple(arcs)
    expected_arcs = ((0, 1), (1, 0), (1, 3), (3, 1), (3, 2), (2, 3))
    require(arcs == expected_arcs, ("directed arc order", arcs))

    adjacency = tuple(tuple(int((left, right) in arcs)
                            for right in STATE_PATH) for left in STATE_PATH)
    pair_status = []
    for left_index, left in enumerate(STATE_PATH):
        for right in STATE_PATH[left_index + 1:]:
            forward = (left, right) in arcs
            backward = (right, left) in arcs
            pair_status.append("both" if forward and backward else
                               "single" if forward or backward else "missing")
    status_counts = tuple(pair_status.count(name)
                          for name in ("both", "single", "missing"))
    require(status_counts == (3, 0, 3), ("partial tournament status", status_counts))

    edge_basis = tuple(standard_basis(6, edge) for edge in range(6))
    state_tail = tuple(
        tuple(int(tail == state) for tail, _head in arcs) for state in STATE_PATH
    )
    root_even = tuple(
        add(edge_basis[left], edge_basis[right])
        for left, right in ((0, 1), (2, 3), (4, 5))
    )
    root_odd = tuple(
        subtract(edge_basis[left], edge_basis[right])
        for left, right in ((0, 1), (2, 3), (4, 5))
    )

    state_images = tuple(linear_image(row, relation_rows) for row in state_tail)
    parent_relation = tuple(
        tuple(sum(parent[state][difference][relation]
                  for difference in range(P)) % PRIME
              for relation in range(P))
        for state in STATE_PATH
    )
    require(state_images == parent_relation,
            ("tail-star/Boolean-parent mismatch", state_images, parent_relation))

    root_images = tuple(linear_image(row, relation_rows) for row in root_even)
    relation_rank = rank(relation_rows)
    require(relation_rank == 6, ("edge response rank", relation_rank))
    require(rank(state_images) == 4, "tail-star rank")
    require(rank(root_images) == 3, "root-even rank")
    relation_digest = rowspace_digest(relation_rows)
    parent_digest = rowspace_digest(state_images)
    require(relation_digest == EXPECTED_COMMON_RELATION_ROWSPACE,
            ("common relation rowspace", relation_digest))
    require(parent_digest == EXPECTED_BOOLEAN_PARENT_ROWSPACE,
            ("Boolean parent rowspace", parent_digest))

    # Standard P7 homology: six alternating-incidence edges, seven vertices.
    p7_boundary = tuple(
        tuple((-1 if vertex == edge else 1 if vertex == edge + 1 else 0)
              for edge in range(6))
        for vertex in range(7)
    )
    p7_boundary_rank = rank(p7_boundary)
    require(p7_boundary_rank == 6, ("P7 boundary rank", p7_boundary_rank))

    # The arc-reversal-odd quotient is the ordinary oriented C1(P4).
    # Its boundary is injective because P4 is a tree.  The middle current has
    # nonzero divergence, so it cannot be an H1 class.
    arc_boundary = tuple(
        tuple((-1 if state == tail else 1 if state == head else 0)
              for tail, head in arcs)
        for state in STATE_PATH
    )
    odd_boundaries = tuple(
        tuple(sum(coefficient * arc_boundary[state][arc] for arc, coefficient
                  in enumerate(row)) % PRIME for state in range(4))
        for row in root_odd
    )
    require(rank(odd_boundaries) == 3, ("P4 current boundary", odd_boundaries))
    require(any(odd_boundaries[1]), "middle current boundary")

    require(tuple(ARC_REVERSAL[PATH_REFLECTION[index]] for index in range(6))
            == tuple(PATH_REFLECTION[ARC_REVERSAL[index]] for index in range(6)),
            "arc reversal/path reflection do not commute")
    representation = (
        module_split(edge_basis, ARC_REVERSAL),
        module_split(edge_basis, PATH_REFLECTION),
        module_split(state_tail,
                     tuple(STATE_PATH.index(state ^ 2) for state in STATE_PATH)),
        module_split(root_even, (2, 1, 0)),
    )
    require(representation == ((6, 3, 3), (6, 3, 3),
                               (4, 2, 2), (3, 2, 1)),
            ("arc representations", representation))

    return (
        p7_edges,
        arcs,
        adjacency,
        status_counts,
        relation_rank,
        rank(state_images),
        rank(root_images),
        relation_digest,
        parent_digest,
        rowspace_digest(root_images),
        representation,
        (p7_boundary_rank, 6 - p7_boundary_rank),
        (rank(odd_boundaries), tuple(odd_boundaries[1])),
    ), edge_basis, root_even, root_odd, state_tail


def source_profile_vector(cells, point: int, r0: int, r1: int):
    state, root = POINTS[point]
    return tuple(
        cell[4][root][r0][r1]
        if cell[0] == state and cell[1][root] else 0
        for cell in cells
    )


def block_vector(point: int, row, cell_count: int):
    return ((0,) * (point * cell_count) + tuple(row)
            + (0,) * ((5 - point) * cell_count))


def reflect_block_vector(row, cell_count: int):
    return tuple(
        row[PATH_REFLECTION[point] * cell_count + (cell_count - 1 - cell)]
        for point in range(6) for cell in range(cell_count)
    )


def tent_source_record(cells, source_ratio_record):
    expected_groups = tuple(
        ((point, POINTS[point], TENT_SECTION[point]),
         tuple((r1, 3 if point in (0, 5) else 2)
               for r1 in NONROOT_DIGITS))
        for point in range(6)
    )
    require(source_ratio_record[0][2] == expected_groups,
            ("exception graph", source_ratio_record[0][2]))
    require(source_ratio_record[0][1] == 60, "exception count")
    require(source_ratio_record[1][1] == 0, "support hostile exceptions")

    drops = tuple(TENT_SECTION[index] - TENT_SECTION[index + 1]
                  for index in range(5))
    line_laplacian = tuple(
        ((1 if index in (0, 5) else 2) * TENT_SECTION[index]
         - (TENT_SECTION[index - 1] if index else 0)
         - (TENT_SECTION[index + 1] if index < 5 else 0))
        for index in range(6)
    )
    require(drops == (0, 3, 6, 3, 0), ("tent drops", drops))
    require(line_laplacian == (0, 3, 3, -3, -3, 0),
            ("line-graph Laplacian", line_laplacian))
    require(all(TENT_SECTION[PATH_REFLECTION[edge]]
                == 12 - TENT_SECTION[edge] for edge in range(6)),
            "affine tent reflection")

    cell_count = len(cells)
    children = []
    parents = []
    stalk_records = []
    for point, r0 in enumerate(TENT_SECTION):
        rows = tuple(source_profile_vector(cells, point, r0, r1)
                     for r1 in range(P))
        parent = tuple(sum(row[cell] for row in rows) % PRIME
                       for cell in range(cell_count))
        nonroot = tuple(rows[r1] for r1 in NONROOT_DIGITS)
        raw_rank = rank(nonroot)
        quotient_rank = rank((parent,) + nonroot) - 1
        all_rank = rank(rows)
        stalk_records.append((raw_rank, quotient_rank, all_rank,
                              digest_json((parent, nonroot))))
        parents.append(block_vector(point, parent, cell_count))
        children.extend(block_vector(point, row, cell_count) for row in nonroot)

    stalk_ranks = tuple(record[:3] for record in stalk_records)
    require(stalk_ranks == ((3, 2, 3), (2, 1, 2), (2, 1, 2),
                            (2, 1, 2), (2, 1, 2), (3, 2, 3)),
            ("defect stalk ranks", stalk_ranks))
    require(rank(children) == 14 and rank(parents) == 6
            and rank(tuple(children) + tuple(parents)) == 14,
            "global defect quotient ranks")

    reflected_children = tuple(reflect_block_vector(row, cell_count)
                               for row in children)
    reflected_parents = tuple(reflect_block_vector(row, cell_count)
                              for row in parents)
    require(rank(tuple(children) + reflected_children) == 14,
            "defect space reflection")
    require(rank(tuple(parents) + reflected_parents) == 6,
            "parent lines reflection")
    plus = tuple(add(row, reflected) for row, reflected
                 in zip(children, reflected_children))
    minus = tuple(subtract(row, reflected) for row, reflected
                  in zip(children, reflected_children))
    quotient_split = (
        rank(tuple(parents) + plus) - rank(parents),
        rank(tuple(parents) + minus) - rank(parents),
    )
    require(quotient_split == (4, 4),
            ("defect quotient reflection split", quotient_split))

    return (
        TENT_SECTION,
        drops,
        line_laplacian,
        NONROOT_DIGITS,
        tuple(stalk_records),
        (rank(children), rank(parents),
         rank(tuple(children) + tuple(parents)) - rank(parents)),
        quotient_split,
        digest_json(tuple(children)),
        digest_json(tuple(parents)),
    )


def diagonal_coefficient_record(cells, root_even, root_odd):
    kernels, source_ratio_record = BUNDLE.source_ratio_kernels(cells)
    actual = kernels[0]
    support = kernels[1]
    require(BUNDLE.POINTS == POINTS, "bundle point order")
    require(BUNDLE.PRIME == PRIME, "bundle field")

    pair_records = []
    unequal = []
    for left, right in ((0, 1), (2, 3), (4, 5)):
        comparable = equal = missing = 0
        for r0 in range(P):
            for r1 in range(P):
                first = actual[r0][r1][left]
                second = actual[r0][r1][right]
                if first is None or second is None:
                    missing += 1
                    continue
                comparable += 1
                if first == second:
                    equal += 1
                else:
                    unequal.append((left, right, r0, r1, first, second))
        pair_records.append((left, right, comparable, equal,
                             comparable - equal, missing))
    require(tuple(record[2:] for record in pair_records)
            == ((159, 159, 0, 10), (149, 143, 6, 20),
                (159, 159, 0, 10)),
            ("root-pair equality", pair_records))
    require(tuple((row[2], row[3]) for row in unequal)
            == ((0, 11), (0, 12), (6, 5), (6, 7), (12, 0), (12, 1)),
            ("central inequalities", unequal))

    sym = tuple(root_even)
    middle_current = root_odd[1]
    middle_space = sym + (middle_current,)
    nonmultiples = (1, 2, 4, 5, 7, 8, 10, 11)
    complete_spaces = []
    for r0 in nonmultiples + (6,):
        rows = tuple(tuple(actual[r0][r1][point] for point in range(6))
                     for r1 in range(P))
        require(all(value is not None for row in rows for value in row),
                ("unexpected missing coefficient", r0))
        target = sym if r0 in nonmultiples else middle_space
        require(rank(rows) == len(target) and rank(rows + target) == len(target),
                ("coefficient module", r0, rank(rows), rank(rows + target)))
        complete_spaces.append((r0, rank(rows), rowspace_digest(rows)))

    # At r0=3 and 9 the only missing source ratios are on one of the two
    # middle arcs.  Every outer pair is completely comparable and equal to
    # the output diagonal coefficient by b1baa781a.  Thus the completed
    # coefficient rows lie in the four-space ``middle_space``.  The clean-room
    # two-current conditional rank is four; injectivity of the six pointed
    # response rows forces equality with this four-space.
    conditional = TWO.EXPECTED_CONDITIONAL_RANKS[0][0]
    require(conditional == (4, 3, 3, 4, 3, 3, 4, 3, 3, 4, 3, 3, 4),
            ("conditional parent", conditional))
    for r0, missing_point in ((3, 3), (9, 2)):
        require(all(
            actual[r0][r1][point] is not None
            for r1 in range(P) for point in (0, 1, 4, 5)
        ), ("outer coefficient missing", r0))
        require(all(actual[r0][r1][0] == actual[r0][r1][1]
                    and actual[r0][r1][4] == actual[r0][r1][5]
                    for r1 in range(P)),
                ("outer pair mismatch", r0))
        require(tuple(point for point in range(6)
                      if any(actual[r0][r1][point] is None
                             for r1 in range(P))) == (missing_point,),
                ("middle missing point", r0))
        require(conditional[r0] == 4, ("conditional completion rank", r0))

    interior_modules = tuple(
        (r0, "root-even" if r0 % 3 else "root-even+middle-current",
         conditional[r0])
        for r0 in range(1, 12)
    )
    require(tuple(item[2] for item in interior_modules)
            == conditional[1:12], "interior module ranks")

    require(all(
        support[r0][r1][point] is not None
        for r0 in range(P) for r1 in range(P) for point in range(6)
    ), "support-normalized source kernel")

    return (
        tuple(pair_records),
        tuple(unequal),
        tuple(complete_spaces),
        rowspace_digest(sym),
        rowspace_digest(middle_space),
        interior_modules,
        source_ratio_record,
        BUNDLE.EXPECTED_SEMANTIC_SHA256,
        TWO.EXPECTED_SEMANTIC_SHA256,
    ), source_ratio_record


def main() -> None:
    require(PRIME == BUNDLE.PRIME == TWO.PRIME, "field mismatch")
    require(POINTS == ((0, 0), (1, 0), (1, 6),
                       (3, 6), (3, 12), (2, 12)), "point order")

    tensor, parent, profile_digest = reconstruct_pointed_tensor()
    relation_rows = tuple(
        tuple(sum(tensor[point][difference][relation]
                  for difference in range(P)) % PRIME
              for relation in range(P))
        for point in range(6)
    )
    (arc_record, edge_basis, root_even,
     root_odd, _state_tail) = p7_arc_record(relation_rows, parent)

    _profiles, _boundaries, cells, source_record = BUNDLE.T2.two_digit_context()
    require(len(cells) == 33, ("source cell count", len(cells)))
    coefficient_record, source_ratio_record = diagonal_coefficient_record(
        cells, root_even, root_odd
    )
    tent_record = tent_source_record(cells, source_ratio_record)

    # Static classification of the diagonal K bundle.  Its exact diagonal and
    # fixed-r0 partition-of-unity properties are pinned by the parent semantic;
    # this script does not rerun the five-coordinate endpoint computation.
    classification = (
        "six_scalar_cylinder_sections",
        "fixed_r0_sum_over_r1_equals_one",
        "diagonal_endomorphisms_of_trivial_arc_bundle",
        "not_probability_without_positivity",
        "not_local_system_noninvertible_and_no_transport",
        "not_H1_no_clock_no_composition_and_H1_P7_zero",
    )

    record = (
        PRIME,
        POINTS,
        STATE_PATH,
        STATE_NAMES,
        ARC_REVERSAL,
        PATH_REFLECTION,
        arc_record,
        coefficient_record,
        tent_record,
        classification,
        profile_digest,
        source_record[-2:],
        digest_json(tensor),
        digest_json(parent),
    )
    record_sha = digest_json(record[:-4])
    semantic = digest_json(record)
    if EXPECTED_RECORD_SHA256 != "TO_BE_PINNED":
        require(record_sha == EXPECTED_RECORD_SHA256,
                ("record drift", record_sha, EXPECTED_RECORD_SHA256))
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("semantic drift", semantic, EXPECTED_SEMANTIC_SHA256))

    print("== r=5 P7 / bidirected-P4 pointed-carrier audit ==")
    print(f"field=F_{PRIME};pointed_order={POINTS}")
    print(f"P7_vertices={P7_VERTICES}")
    print(f"directed_arcs_(tail,head)={arc_record[1]}")
    print(f"partial_tournament_pair_counts_(both,single,missing)={arc_record[3]}")
    print(f"adjacency_in_state_path_order_{STATE_PATH}={arc_record[2]}")
    print("edge_response_order=(rank,tail_star_rank,root_even_rank,relation_sha,parent_sha,root_even_sha)")
    print(f"edge_response={(arc_record[4], arc_record[5], arc_record[6], arc_record[7], arc_record[8], arc_record[9])}")
    print(f"representations_(arc_reversal,path_reflection,tail_star,root_even)={arc_record[10]}")
    print(f"P7_boundary_(rank,H1_dim)={arc_record[11]}")
    print(f"ordinary_P4_current_boundary_(rank,middle_current_boundary)={arc_record[12]}")
    print("root_pair_record_order=(left_arc,right_arc,comparable,equal,unequal,missing)")
    print(f"root_pair_source_to_output_records={coefficient_record[0]}")
    print(f"central_root_unequal_coefficients={coefficient_record[1]}")
    print(f"complete_coefficient_spaces_(r0,rank,rowspace_sha)={coefficient_record[2]}")
    print(f"coefficient_module_shas_(root_even,plus_middle_current)={(coefficient_record[3], coefficient_record[4])}")
    print(f"interior_r0_module_classification={coefficient_record[5]}")
    print(f"tent_section_(values,drops,line_laplacian,nonroot_digits)={tent_record[:4]}")
    print(f"defect_stalk_records_(raw_rank,quotient_rank,all_rank,digest)={tent_record[4]}")
    print(f"defect_global_(raw,parent,quotient)={tent_record[5]};reflection_+-={tent_record[6]}")
    print(f"classification={classification}")
    print(f"record_sha256={record_sha}")
    print(f"semantic_sha256={semantic}")
    print("verdict=pointed relation rows are the injective P7-edge/bidirected-P4 arc module;interior rank3/4 is arc-even versus arc-even plus the middle current")
    print("tent_verdict=exception support is a reflection-equivariant edge section with an 8D decorated-edge quotient defect;not a divisor without a valuation or zero-locus theorem")
    print("typing=arcs have owner-state tails and root-paired reversals;K4 six-edge cardinality is not an incidence-preserving identification")
    print("scope=derived finite exact static response sidecar;no chronology,local system,H1 cocycle,physical current,row exclusion,or LRC14")


if __name__ == "__main__":
    main()
