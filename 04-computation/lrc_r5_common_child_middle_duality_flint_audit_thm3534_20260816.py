#!/usr/bin/env python3
"""Independent FLINT audit of the new THM-3534 finite maps.

This companion does not import the primary custom-RREF implementation.  It
reconstructs the frozen 78-state parent/child bank, computes intersections
through annihilator sums with python-flint, and enumerates chamber-stable
windows through a cached seven-orbit lattice.  It audits the common quotient
core, its middle perfect pairing, the compressed reversal correspondence, the
endpoint-retaining response representation, and the complete 127-window
character census.

The objects are static finite-field representations.  The calculation does
not construct chronology, a physical current, a bridge to a Jacobian flux, or
an LRC(14) closure.
"""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path

from flint import fmpz_mod_ctx, fmpz_mod_mat


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
EXPECTED_COMMON_DIGEST = (
    "584eb2a4022ffed65be254f73eda96facaaaa87b2cc1b80ff8f7714dd74bd65a"
)
EXPECTED_QA_DUAL_DIGEST = (
    "b5143e9ea23f3f981c1b6b334799e7924df0a4a544b69d6dcd18defc01d0f238"
)
EXPECTED_SEMANTIC_SHA256 = (
    "036823a7ba481528a3d02f9b36cecc08c2a76b796d69b437ac644ea62251c53e"
)

DIGITS = 13
ARCS = 6
SIZE = DIGITS * ARCS
ARC_REVERSAL = (1, 0, 3, 2, 5, 4)
CHAMBER_REVERSAL = (5, 4, 3, 2, 1, 0)
REFLECTION_ORBITS = (
    (0, 12), (1, 11), (2, 10), (3, 9), (4, 8), (5, 7), (6,),
)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def digest(value: object) -> str:
    data = json.dumps(
        value, sort_keys=True, separators=(",", ":")
    ).encode("ascii")
    return sha256(data).hexdigest()


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def ints(matrix) -> tuple[tuple[int, ...], ...]:
    return tuple(tuple(int(value) for value in row) for row in matrix.tolist())


def mat(rows, context):
    return fmpz_mod_mat([list(row) for row in rows], context)


def identity(size: int, context):
    return mat(
        tuple(tuple(int(row == column) for column in range(size))
              for row in range(size)),
        context,
    )


def stacked(*matrices, context):
    rows = []
    columns = None
    for matrix in matrices:
        if columns is None:
            columns = matrix.ncols()
        else:
            require(matrix.ncols() == columns, "stack column mismatch")
        rows.extend(ints(matrix))
    require(columns is not None, "empty stack input")
    if not rows:
        return fmpz_mod_mat(0, columns, context)
    return mat(tuple(rows), context)


def row_basis(matrix, context):
    reduced, matrix_rank = matrix.rref()
    if matrix_rank == 0:
        return fmpz_mod_mat(0, matrix.ncols(), context)
    return mat(ints(reduced)[:matrix_rank], context)


def kernel_rows(matrix, prime: int, context):
    """Return a row matrix whose rows form the right kernel of matrix."""
    reduced, matrix_rank = matrix.rref()
    reduced_rows = ints(reduced)
    columns = matrix.ncols()
    pivots = []
    for row in reduced_rows[:matrix_rank]:
        pivots.append(next(column for column, value in enumerate(row) if value))
    free = tuple(column for column in range(columns) if column not in pivots)
    if not free:
        return fmpz_mod_mat(0, columns, context)
    vectors = []
    for free_column in free:
        vector = [0] * columns
        vector[free_column] = 1
        for row_index, pivot in enumerate(pivots):
            vector[pivot] = -reduced_rows[row_index][free_column] % prime
        vectors.append(tuple(vector))
    answer = mat(tuple(vectors), context)
    require((matrix * answer.transpose()).rank() == 0, "kernel verification")
    # Keep the free-variable basis in increasing free-column order.  This is
    # the frozen quotient-dual chart used by the primary certificate; an extra
    # row reduction would preserve every subspace and rank while changing the
    # displayed pairing scalar.
    return answer


def intersection(left, right, prime: int, context):
    """Row(left) intersect Row(right), via (L^perp + R^perp)^perp."""
    left = row_basis(left, context)
    right = row_basis(right, context)
    left_annihilator = kernel_rows(left, prime, context)
    right_annihilator = kernel_rows(right, prime, context)
    annihilator_sum = row_basis(
        stacked(left_annihilator, right_annihilator, context=context),
        context,
    )
    return kernel_rows(annihilator_sum, prime, context)


def equal_rowspaces(left, right, context) -> bool:
    return (
        left.rank() == right.rank()
        and stacked(left, right, context=context).rank() == left.rank()
    )


def permutation_matrix(point_map, digit_map, context):
    rows = [[0] * SIZE for _ in range(SIZE)]
    for point in range(ARCS):
        for digit in range(DIGITS):
            old = point * DIGITS + digit
            new = point_map[point] * DIGITS + digit_map[digit]
            rows[old][new] = 1
    answer = mat(tuple(tuple(row) for row in rows), context)
    require(answer.rank() == SIZE, "coordinate permutation")
    return answer


def embedded_block_rows(arc: int, local_rows, context):
    rows = []
    for local in ints(local_rows):
        row = [0] * SIZE
        row[arc * DIGITS:(arc + 1) * DIGITS] = local
        rows.append(tuple(row))
    return mat(tuple(rows), context)


def build_rows(k2, k3, prime: int, context):
    parent_rows = []
    child_rows = [[] for _ in range(DIGITS)]
    for arc in range(ARCS):
        for r1 in range(DIGITS):
            parent = [0] * SIZE
            children = [[0] * SIZE for _ in range(DIGITS)]
            for r0 in range(DIGITS):
                column = arc * DIGITS + r0
                value = k2[arc][r0][r1] % prime
                parent[column] = value
                for r2 in range(DIGITS):
                    children[r2][column] = (
                        value * k3[arc][r0][r1][r2] % prime
                    )
            parent_rows.append(tuple(parent))
            for r2 in range(DIGITS):
                child_rows[r2].append(tuple(children[r2]))
    return (
        mat(tuple(parent_rows), context),
        tuple(mat(tuple(rows), context) for rows in child_rows),
    )


def select_quotient_lifts(rows, quotient_dual, dimension: int, context):
    """Greedy canonical lifts, with all rank decisions delegated to FLINT."""
    selected = []
    coordinates = []
    for row in ints(row_basis(rows, context)):
        row_matrix = mat((row,), context)
        coordinate = ints(row_matrix * quotient_dual.transpose())[0]
        trial = mat(tuple(coordinates + [coordinate]), context)
        if trial.rank() > len(coordinates):
            selected.append(row)
            coordinates.append(coordinate)
        if len(selected) == dimension:
            break
    require(len(selected) == dimension, "quotient lifts")
    return mat(tuple(selected), context), mat(tuple(coordinates), context)


def pivot_columns(matrix):
    reduced, matrix_rank = matrix.rref()
    return tuple(
        next(column for column, value in enumerate(row) if value)
        for row in ints(reduced)[:matrix_rank]
    )


def coordinates_in_basis(transformed, basis, context):
    """Coordinates of transformed rows in an independent row basis."""
    pivots = pivot_columns(basis)
    require(len(pivots) == basis.nrows(), "basis pivot count")
    basis_square = mat(
        tuple(tuple(int(basis[row, column]) for column in pivots)
              for row in range(basis.nrows())),
        context,
    )
    target_square = mat(
        tuple(tuple(int(transformed[row, column]) for column in pivots)
              for row in range(transformed.nrows())),
        context,
    )
    answer = target_square * basis_square.inv()
    require(equal_rowspaces(transformed, answer * basis, context),
            "coordinate reconstruction")
    return answer


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


def main() -> None:
    bank = load_bank()
    prime = bank["prime"]
    context = fmpz_mod_ctx(prime)
    parent_raw, child_raw = build_rows(bank["k2"], bank["k3"], prime, context)
    parent = row_basis(parent_raw, context)
    children = tuple(row_basis(child, context) for child in child_raw)
    quotient_dual = kernel_rows(parent, prime, context)
    require((parent.rank(), quotient_dual.rank()) == (68, 10),
            "parent quotient dimensions")

    common = children[0]
    for child in children[1:]:
        common = intersection(common, child, prime, context)
    common_kernel = intersection(common, parent, prime, context)
    common_image = row_basis(common * quotient_dual.transpose(), context)
    require((common.rank(), common_kernel.rank(), common_image.rank())
            == (14, 12, 2), "common-child subquotient")

    arc_permutation = permutation_matrix(
        ARC_REVERSAL, tuple(range(DIGITS)), context
    )
    reversed_parent = parent * arc_permutation
    arc_stabilized_parent = row_basis(
        stacked(parent, reversed_parent, context=context), context
    )
    reversal_kernel = row_basis(
        reversed_parent * quotient_dual.transpose(), context
    )
    require((arc_stabilized_parent.rank(), reversal_kernel.rank()) == (70, 2),
            "arc-stabilized dimensions")
    require(equal_rowspaces(common_image, reversal_kernel, context),
            "common image equals reversal kernel")
    common_digest = digest(ints(row_basis(common_image, context)))
    reversal_digest = digest(ints(row_basis(reversal_kernel, context)))
    require(common_digest == reversal_digest == EXPECTED_COMMON_DIGEST,
            (common_digest, reversal_digest))

    common_annihilator_coefficients = kernel_rows(
        common_image, prime, context
    )
    common_annihilator = row_basis(
        common_annihilator_coefficients * quotient_dual, context
    )
    qa_dual = kernel_rows(arc_stabilized_parent, prime, context)
    require(common_annihilator.rank() == qa_dual.rank() == 8,
            "common annihilator dimensions")
    require(equal_rowspaces(common_annihilator, qa_dual, context),
            "common annihilator equals QA dual")
    common_annihilator_digest = digest(ints(common_annihilator))
    qa_dual_digest = digest(ints(row_basis(qa_dual, context)))
    require(common_annihilator_digest == qa_dual_digest
            == EXPECTED_QA_DUAL_DIGEST,
            (common_annihilator_digest, qa_dual_digest))

    k2 = bank["k2"]
    blocks = tuple(
        mat(tuple(tuple(k2[arc][r0][r1] % prime
                        for r1 in range(DIGITS))
                  for r0 in range(DIGITS)), context)
        for arc in range(ARCS)
    )
    block_duals = tuple(
        embedded_block_rows(
            arc, kernel_rows(block.transpose(), prime, context), context
        )
        for arc, block in enumerate(blocks)
    )
    block_dimensions = tuple(block.rank() for block in block_duals)
    require(block_dimensions == (2, 2, 1, 1, 2, 2), block_dimensions)
    block_pairing_ranks = tuple(
        (common * block_dual.transpose()).rank()
        for block_dual in block_duals
    )
    require(block_pairing_ranks == (0, 0, 1, 1, 0, 0),
            block_pairing_ranks)

    ell_plus_local = ints(block_duals[2])[0][2 * DIGITS:3 * DIGITS]
    ell_minus_local = tuple(reversed(ell_plus_local))
    ell_plus = embedded_block_rows(
        2, mat((ell_plus_local,), context), context
    )
    ell_minus = embedded_block_rows(
        3, mat((ell_minus_local,), context), context
    )
    middle_dual = stacked(ell_plus, ell_minus, context=context)
    require(equal_rowspaces(block_duals[3], ell_minus, context),
            "chamber-linked middle minus")

    common_lifts, _ = select_quotient_lifts(
        common, quotient_dual, 2, context
    )
    middle_pairing = common_lifts * middle_dual.transpose()
    expected_pairing_scalar = 636675481197456361648540
    expected_pairing_determinant = 149750845022728455688979
    require(ints(middle_pairing) == (
        (expected_pairing_scalar, 0),
        (0, expected_pairing_scalar),
    ), ints(middle_pairing))
    require(int(middle_pairing.det()) == expected_pairing_determinant,
            int(middle_pairing.det()))
    dual_lifts = middle_pairing.inv() * common_lifts
    require(ints(dual_lifts * middle_dual.transpose()) == ((1, 0), (0, 1)),
            "normalized middle-dual lifts")

    difference = blocks[2] - blocks[3]
    ell_plus_row = mat((ell_plus_local,), context)
    ell_minus_row = mat((ell_minus_local,), context)
    response_source = mat(
        (ints(difference)[3], ints(difference)[9]), context
    )
    response_image = stacked(
        ell_plus_row * difference,
        ell_minus_row * difference,
        context=context,
    )
    endpoint_six = mat((ints(difference)[6],), context)
    response_extension = stacked(
        response_image, endpoint_six, context=context
    )
    require((response_source.rank(), response_image.rank(),
             response_extension.rank()) == (2, 2, 3),
            "response extension dimensions")
    require(stacked(response_source, response_image,
                    context=context).rank() == 3,
            "source-image intersection")
    require(stacked(response_source, endpoint_six,
                    context=context).rank() == 3,
            "source relative injection")
    require(stacked(response_image, endpoint_six,
                    context=context).rank() == 3,
            "image relative injection")

    digit_reflection_13 = permutation_matrix(
        tuple(range(ARCS)),
        tuple(DIGITS - 1 - digit for digit in range(DIGITS)),
        context,
    )
    local_reflection = mat(
        tuple(tuple(int(new == DIGITS - 1 - old)
                    for new in range(DIGITS))
              for old in range(DIGITS)),
        context,
    )
    response_chamber_images = mat(
        tuple(tuple(-int(value) % prime for value in row)
              for row in ints(response_extension * local_reflection)),
        context,
    )
    extension_chamber = coordinates_in_basis(
        response_chamber_images, response_extension, context
    )
    expected_extension_chamber = mat(
        ((0, 1, 0), (1, 0, 0), (0, 0, 1)), context
    )
    require(ints(extension_chamber) == ints(expected_extension_chamber),
            ints(extension_chamber))
    extension_boundary_rank = (
        extension_chamber - identity(3, context)
    ).rank()
    relative_boundary_rank = mat(
        ((-1 % prime, 1), (1, -1 % prime)), context
    ).rank()
    extension_h0_h1 = (
        3 - extension_boundary_rank, 3 - extension_boundary_rank
    )
    relative_h0_h1 = (
        2 - relative_boundary_rank, 2 - relative_boundary_rank
    )
    require((extension_boundary_rank, extension_h0_h1,
             relative_boundary_rank, relative_h0_h1)
            == (1, (2, 2), 1, (1, 1)),
            "twisted cohomology ranks")

    chamber_permutation = permutation_matrix(
        CHAMBER_REVERSAL,
        tuple(DIGITS - 1 - digit for digit in range(DIGITS)),
        context,
    )
    chamber_common = common * chamber_permutation
    chamber_kernel = common_kernel * chamber_permutation
    require(equal_rowspaces(common, chamber_common, context),
            "chamber preserves common realization")
    require(equal_rowspaces(common_kernel, chamber_kernel, context),
            "chamber preserves common kernel")
    require(equal_rowspaces(ell_plus * chamber_permutation,
                            ell_minus, context),
            "chamber swaps middle dual lines")
    require(equal_rowspaces(ell_minus * chamber_permutation,
                            ell_plus, context),
            "chamber swaps middle dual lines back")
    chamber_matrix = (
        dual_lifts * chamber_permutation * middle_dual.transpose()
    )
    swap = mat(((0, 1), (1, 0)), context)
    require(ints(chamber_matrix) == ints(swap), ints(chamber_matrix))

    reversed_common = common * arc_permutation
    reversed_kernel = common_kernel * arc_permutation
    common_arc_intersection = intersection(
        common, reversed_common, prime, context
    )
    require(equal_rowspaces(common_kernel, reversed_kernel, context),
            "arc preserves common kernel")
    arc_realization = (
        common.rank(),
        stacked(common, reversed_common, context=context).rank(),
        common_arc_intersection.rank(),
    )
    require(arc_realization == (14, 16, 12), arc_realization)
    require(equal_rowspaces(common_arc_intersection,
                            common_kernel, context),
            "arc intersection equals parent kernel")
    reversed_common_image = row_basis(
        reversed_common * quotient_dual.transpose(), context
    )
    require(equal_rowspaces(common_image, reversed_common_image, context),
            "arc returns to common quotient image")
    compressed_arc = (
        dual_lifts * arc_permutation * middle_dual.transpose()
    )
    arc_scalar = 145859431184888028092125
    arc_scalar_square = 17748677861075734903229
    expected_compressed_arc = mat(
        ((0, arc_scalar), (arc_scalar, 0)), context
    )
    require(ints(compressed_arc) == ints(expected_compressed_arc),
            ints(compressed_arc))
    require(ints(compressed_arc * compressed_arc)
            == ((arc_scalar_square, 0), (0, arc_scalar_square)),
            ints(compressed_arc * compressed_arc))
    require(arc_scalar_square != 1, "compressed arc is not an involution")

    digit_common = common * digit_reflection_13
    digit_kernel = common_kernel * digit_reflection_13
    digit_image = row_basis(digit_common * quotient_dual.transpose(), context)
    digit_gate = (
        stacked(parent, digit_kernel, context=context).rank(),
        digit_image.rank(),
        stacked(common_image, digit_image, context=context).rank(),
    )
    require(digit_gate == (76, 10, 10), digit_gate)

    quotient_lifts, quotient_coordinates = select_quotient_lifts(
        children[0], quotient_dual, 10, context
    )
    normalized_quotient_lifts = quotient_coordinates.inv() * quotient_lifts
    require(ints(normalized_quotient_lifts * quotient_dual.transpose())
            == ints(identity(10, context)),
            "normalized quotient lifts")
    global_chamber = (
        normalized_quotient_lifts
        * chamber_permutation
        * quotient_dual.transpose()
    )
    require(ints(global_chamber * global_chamber)
            == ints(identity(10, context)),
            "global chamber involution")
    chamber_plus = kernel_rows(
        (global_chamber - identity(10, context)).transpose(),
        prime,
        context,
    )
    chamber_minus = kernel_rows(
        (global_chamber + identity(10, context)).transpose(),
        prime,
        context,
    )
    require((chamber_plus.rank(), chamber_minus.rank()) == (5, 5),
            "global chamber characters")

    orbit_spaces = []
    for orbit in REFLECTION_ORBITS:
        space = children[orbit[0]]
        for digit in orbit[1:]:
            space = intersection(space, children[digit], prime, context)
        orbit_spaces.append(space)
    window_cache = {}
    stable_counts = {}
    exceptional_windows = []
    for mask in range(1, 1 << len(REFLECTION_ORBITS)):
        least_bit = mask & -mask
        orbit_index = least_bit.bit_length() - 1
        previous = mask ^ least_bit
        if previous:
            window = intersection(
                window_cache[previous], orbit_spaces[orbit_index],
                prime, context,
            )
        else:
            window = orbit_spaces[orbit_index]
        window_cache[mask] = window
        window_image = row_basis(window * quotient_dual.transpose(), context)
        plus_dimension = intersection(
            window_image, chamber_plus, prime, context
        ).rank()
        minus_dimension = intersection(
            window_image, chamber_minus, prime, context
        ).rank()
        character = (
            window_image.rank(), plus_dimension, minus_dimension
        )
        stable_counts[character] = stable_counts.get(character, 0) + 1
        require(plus_dimension + minus_dimension == window_image.rank()
                and plus_dimension == minus_dimension,
                ("unbalanced stable window", mask, character))
        if window_image.rank() > 2:
            digits = tuple(
                digit
                for index, orbit in enumerate(REFLECTION_ORBITS)
                if (mask >> index) & 1
                for digit in orbit
            )
            exceptional_windows.append(
                (window_image.rank(), digits, window.rank())
            )

    stable_histogram = tuple(sorted(stable_counts.items()))
    expected_histogram = (
        ((2, 1, 1), 119), ((6, 3, 3), 4), ((10, 5, 5), 4),
    )
    expected_exceptional = (
        (6, (0, 12), 20),
        (6, (2, 10), 18),
        (6, (0, 12, 2, 10), 18),
        (10, (3, 9), 34),
        (10, (5, 7), 30),
        (10, (6,), 50),
        (6, (3, 9, 6), 22),
        (10, (5, 7, 6), 22),
    )
    require(stable_histogram == expected_histogram, stable_histogram)
    require(tuple(exceptional_windows) == expected_exceptional,
            tuple(exceptional_windows))

    record = (
        CERTIFICATE_SHA256,
        prime,
        (parent.rank(), quotient_dual.rank()),
        (common.rank(), common_kernel.rank(), common_image.rank()),
        (arc_stabilized_parent.rank(), reversal_kernel.rank()),
        common_digest,
        reversal_digest,
        common_annihilator_digest,
        qa_dual_digest,
        block_dimensions,
        block_pairing_ranks,
        ints(middle_pairing),
        int(middle_pairing.det()),
        (relative_boundary_rank, relative_h0_h1),
        (extension_boundary_rank, extension_h0_h1),
        ints(extension_chamber),
        ints(chamber_matrix),
        arc_realization,
        ints(compressed_arc),
        arc_scalar,
        arc_scalar_square,
        digit_gate,
        stable_histogram,
        tuple(exceptional_windows),
    )
    semantic = digest(record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                (semantic, EXPECTED_SEMANTIC_SHA256))

    print("== independent FLINT audit: THM-3534 new finite maps ==")
    print(f"engine=python-flint.fmpz_mod_mat;certificate_sha256={CERTIFICATE_SHA256};prime={prime}")
    print(f"parent_quotient=(rank,Qdim)={(parent.rank(), quotient_dual.rank())}")
    print(f"common_core=(common,kernel,image)={(common.rank(), common_kernel.rank(), common_image.rank())};reversal_stabilization={(arc_stabilized_parent.rank(), reversal_kernel.rank())}")
    print(f"core_digests=(U={common_digest},reversal_kernel={reversal_digest},equal=True)")
    print(f"dual_digests=(U_perp={common_annihilator_digest},QA_dual={qa_dual_digest},equal=True)")
    print(f"block_duals=(dimensions={block_dimensions},pairing_ranks={block_pairing_ranks})")
    print(f"middle_pairing=(matrix={ints(middle_pairing)},det={int(middle_pairing.det())},perfect=True)")
    print(f"relative_response=(boundary_rank,H0,H1)={(relative_boundary_rank,) + relative_h0_h1};endpoint_extension={(extension_boundary_rank,) + extension_h0_h1};extension_chamber={ints(extension_chamber)}")
    print(f"chamber_on_U={ints(chamber_matrix)};arc_realization={arc_realization}")
    print(f"compressed_arc=(matrix={ints(compressed_arc)},lambda={arc_scalar},lambda_square={arc_scalar_square},involution=False)")
    print(f"digit_reflection_gate={digit_gate}")
    print(f"stable_windows=(histogram={stable_histogram},exceptional={tuple(exceptional_windows)})")
    print(f"semantic_sha256={semantic}")
    print("verdict=ACCEPT every new finite-rank/map/census claim by an independent algebra engine and a distinct window-lattice algorithm")
    print("typing=U=ker(Q10->QA8) pairs perfectly with the marked middle dual;the endpoint D6 class is not supplied by any rank-3 child window")
    print("scope=FINITE-EXACT static representation audit only;no bridge,no physical current,no D5 flux,no row exclusion,no LRC14")
    print("PASS")


if __name__ == "__main__":
    main()
