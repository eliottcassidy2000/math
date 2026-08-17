#!/usr/bin/env python3
"""Exact common-child/middle-dual refinement of the THM-3534 candidate.

The frozen 78-state bank has a rank-68 parent rowspace R and thirteen child
rowspaces C_t.  This postprocessor compares two a priori unrelated rank-two
objects:

    U = q(intersection_t C_t) in Q=W/R,
    ker(Q -> Q_A),  Q_A=W/(R+A R),

where A is pure reversal of the six pointed arcs.  It also evaluates the
canonical quotient-dual pairing against the two middle block annihilator
lines.  All arithmetic is exact in the pinned split prime field; no endpoint
integration is rerun.

The output is static representation data.  In particular, the compressed
arc operator found below is not an action of arc reversal, because A does not
preserve the common-child realization.  Nothing here constructs chronology,
a physical current, an H1 class, a D5 flux, or an LRC(14) argument.
"""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = (
    ROOT / "05-knowledge/results/"
    "lrc_r5_third_digit_78_state_exact_bank_20260816.json"
)
CERTIFICATE_SHA256 = (
    "472925c638de7ac90b1a7880184766f2e8acec0f8291e2c46cf747b37cd46712"
)
EXPECTED_PARENT_RAW_SHA256 = (
    "a227bc2f385d8a2eaecb27f317fa5ed66623c70938d8a97aba620298a8a7b61b"
)
EXPECTED_PARENT_SEMANTIC_SHA256 = (
    "3d1527fb4ce4931680e50d7135b9d1129c1816e3a9158645523e2728ddc71ec2"
)
EXPECTED_SEMANTIC_SHA256 = (
    "e734d40131404ac2ca0ea1a0fbd4154d998b7bb2e18200351fcbc3debc560ce1"
)

P = 13
ARCS = 6
SIZE = P * ARCS
ARC_REVERSAL = (1, 0, 3, 2, 5, 4)
CHAMBER_REVERSAL = (5, 4, 3, 2, 1, 0)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def digest(value: object) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":")
    ).encode("ascii")
    return sha256(encoded).hexdigest()


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def rref(rows, prime: int, columns: int | None = None):
    if columns is None:
        columns = len(rows[0]) if rows else 0
    work = [[value % prime for value in row] for row in rows]
    require(all(len(row) >= columns for row in work), "ragged matrix")
    pivot_row = 0
    pivots = []
    for column in range(columns):
        pivot = next(
            (row for row in range(pivot_row, len(work))
             if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        unit = pow(work[pivot_row][column], -1, prime)
        work[pivot_row] = [value * unit % prime
                           for value in work[pivot_row]]
        for row in range(len(work)):
            if row == pivot_row or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [
                (left - factor * right) % prime
                for left, right in zip(work[row], work[pivot_row])
            ]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == len(work):
            break
    return (
        tuple(tuple(row) for row in work[:pivot_row]),
        tuple(pivots),
    )


def rank(rows, prime: int, columns: int | None = None) -> int:
    return len(rref(rows, prime, columns)[0])


def row_basis(rows, prime: int, columns: int = SIZE):
    return rref(rows, prime, columns)[0]


def nullspace(rows, prime: int, columns: int | None = None):
    if columns is None:
        columns = len(rows[0]) if rows else 0
    basis, pivots = rref(rows, prime, columns)
    free = tuple(column for column in range(columns) if column not in pivots)
    answer = []
    for free_column in free:
        vector = [0] * columns
        vector[free_column] = 1
        for row, pivot in enumerate(pivots):
            vector[pivot] = -basis[row][free_column] % prime
        answer.append(tuple(vector))
    return tuple(answer)


def transpose(matrix):
    return tuple(
        tuple(matrix[row][column] for row in range(len(matrix)))
        for column in range(len(matrix[0]))
    ) if matrix else ()


def matmul(left, right, prime: int):
    if not left:
        return ()
    return tuple(
        tuple(
            sum(left[row][middle] * right[middle][column]
                for middle in range(len(right))) % prime
            for column in range(len(right[0]))
        )
        for row in range(len(left))
    )


def inverse(matrix, prime: int):
    size = len(matrix)
    require(all(len(row) == size for row in matrix), "inverse square")
    augmented = tuple(
        tuple(row) + tuple(int(i == j) for j in range(size))
        for i, row in enumerate(matrix)
    )
    reduced, pivots = rref(augmented, prime, size)
    require(pivots == tuple(range(size)), ("inverse singular", pivots))
    return tuple(tuple(row[size:]) for row in reduced)


def intersection_basis(left, right, prime: int, columns: int = SIZE):
    left = row_basis(left, prime, columns)
    right = row_basis(right, prime, columns)
    left_perp = nullspace(left, prime, columns)
    right_perp = nullspace(right, prime, columns)
    return row_basis(
        nullspace(left_perp + right_perp, prime, columns),
        prime,
        columns,
    )


def pair_rows(rows, dual_vectors, prime: int):
    return tuple(
        tuple(
            sum(row[column] * vector[column] for column in range(SIZE))
            % prime
            for vector in dual_vectors
        )
        for row in rows
    )


def build_rows(k2, k3, prime: int):
    parent = []
    children = [[] for _ in range(P)]
    for arc in range(ARCS):
        for r1 in range(P):
            parent_row = [0] * SIZE
            child_rows = [[0] * SIZE for _ in range(P)]
            for r0 in range(P):
                column = arc * P + r0
                value = k2[arc][r0][r1] % prime
                parent_row[column] = value
                for r2 in range(P):
                    child_rows[r2][column] = (
                        value * k3[arc][r0][r1][r2] % prime
                    )
            parent.append(tuple(parent_row))
            for r2 in range(P):
                children[r2].append(tuple(child_rows[r2]))
    return tuple(parent), tuple(tuple(child) for child in children)


def coordinate_permutation(point_map, digit_map):
    permutation = [None] * SIZE
    for point in range(ARCS):
        for digit in range(P):
            old = point * P + digit
            new = point_map[point] * P + digit_map[digit]
            permutation[old] = new
    require(sorted(permutation) == list(range(SIZE)), "permutation")
    return tuple(permutation)


def permute_row(row, permutation):
    answer = [0] * len(row)
    for old, new in enumerate(permutation):
        answer[new] = row[old]
    return tuple(answer)


def transform_rows(rows, permutation):
    return tuple(permute_row(row, permutation) for row in rows)


def embedded_block_vector(arc: int, vector):
    answer = [0] * SIZE
    answer[arc * P:(arc + 1) * P] = vector
    return tuple(answer)


def deterministic_quotient_lifts(rows, quotient_dual, dimension: int,
                                  prime: int):
    selected = []
    quotient_rows = []
    for row, quotient_row in zip(rows, pair_rows(rows, quotient_dual, prime)):
        if rank(tuple(quotient_rows) + (quotient_row,),
                prime, len(quotient_dual)) > len(quotient_rows):
            selected.append(row)
            quotient_rows.append(quotient_row)
        if len(selected) == dimension:
            break
    require(len(selected) == dimension, "quotient lifts")
    return tuple(selected), tuple(quotient_rows)


def load_bank():
    require(lf_sha256(CERTIFICATE) == CERTIFICATE_SHA256, "certificate hash")
    bank = json.loads(CERTIFICATE.read_text(encoding="ascii"))
    require(bank["schema"] == "r5-third-digit-78-diagonal-bank-v1",
            "certificate schema")
    require(bank["parent_raw_sha256"] == EXPECTED_PARENT_RAW_SHA256,
            "parent raw hash")
    require(bank["parent_semantic_sha256"]
            == EXPECTED_PARENT_SEMANTIC_SHA256,
            "parent semantic hash")
    return bank


def analyse(bank):
    prime = bank["prime"]
    k2 = bank["k2"]
    k3 = bank["k3"]
    parent_rows, children = build_rows(k2, k3, prime)
    parent = row_basis(parent_rows, prime)
    quotient_dual = nullspace(parent, prime, SIZE)
    require((len(parent), len(quotient_dual)) == (68, 10),
            "parent quotient dimensions")

    common_child = row_basis(children[0], prime)
    for child in children[1:]:
        common_child = intersection_basis(common_child, child, prime)
    common_kernel = intersection_basis(common_child, parent, prime)
    common_image = row_basis(
        pair_rows(common_child, quotient_dual, prime), prime, 10
    )
    require((len(common_child), len(common_kernel), len(common_image))
            == (14, 12, 2), "common-child subquotient")

    arc_permutation = coordinate_permutation(
        ARC_REVERSAL, tuple(range(P))
    )
    reversed_parent = transform_rows(parent, arc_permutation)
    arc_stable_parent = row_basis(parent + reversed_parent, prime)
    arc_kernel = row_basis(
        pair_rows(reversed_parent, quotient_dual, prime), prime, 10
    )
    require((len(arc_stable_parent), len(arc_kernel)) == (70, 2),
            "arc-stable kernel")
    require(rank(common_image + arc_kernel, prime, 10) == 2,
            "common image equals reversal kernel")

    arc_stable_dual = nullspace(arc_stable_parent, prime, SIZE)
    common_annihilator_coefficients = nullspace(common_image, prime, 10)
    common_annihilator = tuple(
        tuple(
            sum(coefficients[index] * quotient_dual[index][column]
                for index in range(10)) % prime
            for column in range(SIZE)
        )
        for coefficients in common_annihilator_coefficients
    )
    require(len(common_annihilator) == len(arc_stable_dual) == 8,
            "annihilator dimensions")
    require(rank(common_annihilator + arc_stable_dual, prime, SIZE) == 8,
            "common annihilator equals QA dual")

    blocks = tuple(
        tuple(
            tuple(k2[arc][r0][r1] % prime for r1 in range(P))
            for r0 in range(P)
        )
        for arc in range(ARCS)
    )
    block_duals = []
    block_dual_dimensions = []
    for arc, block in enumerate(blocks):
        local = nullspace(transpose(block), prime, P)
        block_dual_dimensions.append(len(local))
        block_duals.append(tuple(
            embedded_block_vector(arc, vector) for vector in local
        ))
    require(tuple(block_dual_dimensions) == (2, 2, 1, 1, 2, 2),
            block_dual_dimensions)
    require(rank(tuple(vector for block in block_duals for vector in block),
                 prime, SIZE) == 10, "block dual decomposition")

    ell_plus_local = block_duals[2][0][2 * P:3 * P]
    ell_minus_local = tuple(reversed(ell_plus_local))
    ell_plus = embedded_block_vector(2, ell_plus_local)
    ell_minus = embedded_block_vector(3, ell_minus_local)
    middle_dual = (ell_plus, ell_minus)
    require(rank(block_duals[3] + (ell_minus,), prime, SIZE) == 1,
            "chamber-linked middle minus line")

    block_pairing_ranks = tuple(
        rank(pair_rows(common_child, block_duals[arc], prime),
             prime, len(block_duals[arc]))
        for arc in range(ARCS)
    )
    require(block_pairing_ranks == (0, 0, 1, 1, 0, 0),
            block_pairing_ranks)

    common_lifts, _common_coordinates = deterministic_quotient_lifts(
        common_child, quotient_dual, 2, prime
    )
    middle_pairing = pair_rows(common_lifts, middle_dual, prime)
    middle_pairing_determinant = (
        middle_pairing[0][0] * middle_pairing[1][1]
        - middle_pairing[0][1] * middle_pairing[1][0]
    ) % prime
    require(rank(middle_pairing, prime, 2) == 2
            and middle_pairing_determinant,
            "perfect middle pairing")
    require(middle_pairing[0][1] == middle_pairing[1][0] == 0
            and middle_pairing[0][0] == middle_pairing[1][1],
            "chamber-linked scalar pairing")
    pairing_scalar = middle_pairing[0][0]
    dual_lifts = matmul(inverse(middle_pairing, prime), common_lifts, prime)
    require(pair_rows(dual_lifts, middle_dual, prime)
            == ((1, 0), (0, 1)), "dual normalized lifts")

    chamber_permutation = coordinate_permutation(
        CHAMBER_REVERSAL, tuple(P - 1 - digit for digit in range(P))
    )
    chamber_common = transform_rows(common_child, chamber_permutation)
    chamber_kernel = transform_rows(common_kernel, chamber_permutation)
    require(rank(common_child + chamber_common, prime, SIZE) == 14
            and rank(common_kernel + chamber_kernel, prime, SIZE) == 12,
            "chamber preserves common subquotient")
    require(transform_rows((ell_plus,), chamber_permutation)[0] == ell_minus
            and transform_rows((ell_minus,), chamber_permutation)[0] == ell_plus,
            "chamber swaps middle dual")
    chamber_matrix = pair_rows(
        transform_rows(dual_lifts, chamber_permutation), middle_dual, prime
    )
    swap = ((0, 1), (1, 0))
    require(chamber_matrix == swap, chamber_matrix)

    reversed_common = transform_rows(common_child, arc_permutation)
    reversed_kernel = transform_rows(common_kernel, arc_permutation)
    common_arc_intersection = intersection_basis(
        common_child, reversed_common, prime
    )
    require(rank(common_kernel + reversed_kernel, prime, SIZE) == 12,
            "arc preserves common kernel")
    require(rank(common_child + reversed_common, prime, SIZE) == 16
            and rank(common_arc_intersection + common_kernel,
                     prime, SIZE) == 12,
            "arc-transverse common lifts")
    require(rank(
        common_image
        + row_basis(pair_rows(reversed_common, quotient_dual, prime),
                    prime, 10),
        prime,
        10,
    ) == 2, "arc projects back to common image")
    compressed_arc = pair_rows(
        transform_rows(dual_lifts, arc_permutation), middle_dual, prime
    )
    require(compressed_arc[0][0] == compressed_arc[1][1] == 0
            and compressed_arc[0][1] == compressed_arc[1][0],
            "compressed arc is scalar swap")
    arc_scalar = compressed_arc[0][1]
    arc_scalar_square = arc_scalar * arc_scalar % prime
    require(arc_scalar not in (0, 1, prime - 1)
            and compressed_arc
            == tuple(tuple(arc_scalar * swap[row][column] % prime
                           for column in range(2)) for row in range(2)),
            "non-involutive compressed arc")
    require(matmul(compressed_arc, compressed_arc, prime)
            == ((arc_scalar_square, 0), (0, arc_scalar_square)),
            "compressed arc square")

    digit_permutation = coordinate_permutation(
        tuple(range(ARCS)), tuple(P - 1 - digit for digit in range(P))
    )
    digit_common = transform_rows(common_child, digit_permutation)
    digit_kernel = transform_rows(common_kernel, digit_permutation)
    digit_image = row_basis(
        pair_rows(digit_common, quotient_dual, prime), prime, 10
    )
    digit_gate = (
        rank(parent + digit_kernel, prime, SIZE),
        len(digit_image),
        rank(common_image + digit_image, prime, 10),
    )
    require(digit_gate == (76, 10, 10), digit_gate)

    record = (
        CERTIFICATE_SHA256,
        prime,
        (len(parent), len(quotient_dual)),
        (len(common_child), len(common_kernel), len(common_image)),
        (len(arc_stable_parent), len(arc_kernel)),
        digest(common_image),
        digest(arc_kernel),
        len(common_annihilator),
        digest(row_basis(common_annihilator, prime)),
        digest(row_basis(arc_stable_dual, prime)),
        tuple(block_dual_dimensions),
        block_pairing_ranks,
        middle_pairing,
        middle_pairing_determinant,
        pairing_scalar,
        chamber_matrix,
        (len(reversed_common), rank(common_child + reversed_common,
                                    prime, SIZE),
         len(common_arc_intersection)),
        compressed_arc,
        arc_scalar,
        arc_scalar_square,
        digit_gate,
    )
    semantic = digest(record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                (semantic, EXPECTED_SEMANTIC_SHA256))
    return record, semantic


def main() -> None:
    bank = load_bank()
    record, semantic = analyse(bank)
    (
        certificate_sha, prime, parent_q, common_dims, arc_dims,
        common_digest, arc_digest, annihilator_dimension,
        common_annihilator_digest, qa_dual_digest,
        block_dimensions, block_pairing_ranks,
        middle_pairing, middle_pairing_determinant, pairing_scalar,
        chamber_matrix, arc_lift_record, compressed_arc, arc_scalar,
        arc_scalar_square, digit_gate,
    ) = record
    print("== r5 common-child / middle-dual exact refinement ==")
    print(f"frozen_bank=(certificate_sha256={certificate_sha},prime={prime})")
    print(f"parent_quotient=(parent_rank,Q_dimension)={parent_q}")
    print(f"common_subquotient=(intersection_dimension,kernel_in_parent,image_in_Q)={common_dims}")
    print(f"arc_stabilization=(parent_plus_reversal_rank,kernel_dimension)={arc_dims}")
    print(f"core_identity=(common_image_digest={common_digest},arc_kernel_digest={arc_digest},equal=True)")
    print(f"dual_identity=(dimension={annihilator_dimension},common_annihilator_digest={common_annihilator_digest},QA_dual_digest={qa_dual_digest},equal=True)")
    print(f"block_duals=(dimensions={block_dimensions},common_pairing_ranks={block_pairing_ranks})")
    print(f"middle_pairing=(matrix={middle_pairing},det={middle_pairing_determinant},scalar={pairing_scalar},perfect=True,outer_blocks_zero=True)")
    print(f"chamber_subquotient=(matrix={chamber_matrix},genuine_involution=True)")
    print(f"arc_realization=(common_dimension,union_with_reversal,intersection)={arc_lift_record}")
    print(f"compressed_arc=(matrix={compressed_arc},scalar_times_chamber={arc_scalar},square_scalar={arc_scalar_square},involution=False)")
    print(f"digit_reflection_gate=(parent_plus_image_of_common_kernel_rank,image_rank_in_Q,union_with_U_rank)={digit_gate}")
    print(f"semantic_sha256={semantic}")
    print("connection=U=q(intersection_t C_t)=ker(Q10->QA8) pairs perfectly with L=Q10_middle_dual, hence U is canonically V_rel_dual after the THM-3534 contraction L->V_rel")
    print("preserved=all-child common quotient class,middle orientation defect,perfect evaluation pairing,coupled chamber involution")
    print("lost=common literal lift under arc reversal,digit reflection,endpoint D6,chronology,closure edge,current semantics,JC target predicate")
    print("sidecar=a marked self-duality is required to turn V_rel_dual into V_rel;D6 and lawful clock/chain data remain required for any current")
    print("tournament_xor=the intrinsic chamber relation is both-way;the projected arc correspondence is lambda*chamber with lambda not plus_or_minus_one,so XOR erases the exact obstruction")
    print("verdict=the frozen rank-two common child core is exactly the orientation-forgetting kernel and only a dual static response carrier;forcing arc reversal kills it,while compressing reversal produces nontrivial scalar holonomy rather than an action")
    print("scope=FINITE-EXACT static quotient-duality result;no physical current,no H1 bridge,no D5 flux,no row exclusion,no LRC14")
    print("PASS")


if __name__ == "__main__":
    main()
