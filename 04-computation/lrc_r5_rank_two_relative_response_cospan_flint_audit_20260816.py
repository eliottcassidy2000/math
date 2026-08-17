#!/usr/bin/env python3
"""Independent FLINT audit of the THM-3534 finite matrix claims.

This intentionally does not import the primary custom-RREF companion.  It
loads the same immutable compact bank and delegates modular row reduction,
kernel extraction, inversion, and multiplication to python-flint's arbitrary
modulus matrix implementation.
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
EXPECTED_SEMANTIC_SHA256 = (
    "f1dfa230027e7fd8f48c7eb7bac2034d0188ee61a4ce6f24d2579860cf532020"
)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def digest(value: object) -> str:
    data = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    return sha256(data).hexdigest()


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def ints(matrix) -> tuple[tuple[int, ...], ...]:
    return tuple(tuple(int(value) for value in row) for row in matrix.tolist())


def mat(rows, context):
    return fmpz_mod_mat([list(row) for row in rows], context)


def stacked(*matrices, context):
    rows = []
    for matrix in matrices:
        rows.extend(ints(matrix))
    return mat(rows, context)


def kernel_basis(matrix, prime: int, context):
    reduced, matrix_rank = matrix.rref()
    rows = ints(reduced)
    columns = matrix.ncols()
    pivots = []
    for row in rows[:matrix_rank]:
        pivot = next(column for column, value in enumerate(row) if value)
        pivots.append(pivot)
    free = [column for column in range(columns) if column not in pivots]
    result = []
    for free_column in free:
        vector = [0] * columns
        vector[free_column] = 1
        for row_index, pivot in enumerate(pivots):
            vector[pivot] = -rows[row_index][free_column] % prime
        result.append(tuple(vector))
    for vector in result:
        product = matrix * mat(tuple((value,) for value in vector), context)
        require(not any(int(value) for value in product.entries()), "kernel check")
    return tuple(result)


def pivot_columns(matrix):
    reduced, matrix_rank = matrix.rref()
    rows = ints(reduced)
    return tuple(next(column for column, value in enumerate(row) if value)
                 for row in rows[:matrix_rank])


def main() -> None:
    require(lf_sha256(CERTIFICATE) == CERTIFICATE_SHA256, "certificate hash")
    bank = json.loads(CERTIFICATE.read_text(encoding="ascii"))
    prime = bank["prime"]
    context = fmpz_mod_ctx(prime)
    k2 = bank["k2"]
    blocks = tuple(
        mat(tuple(tuple(k2[arc][r0][r1] for r1 in range(13))
                  for r0 in range(13)), context)
        for arc in range(6)
    )
    plus, minus = blocks[2], blocks[3]
    difference = plus - minus

    block_ranks = tuple(block.transpose().rank() for block in blocks)
    ones_column = mat(tuple((1,) for _ in range(13)), context)
    ones_row = mat(((1,) * 13,), context)
    row_sum_controls = tuple(
        tuple(int(value) for value in (block * ones_column).entries())
        for block in blocks
    )
    require(row_sum_controls == ((1,) * 13,) * 6, row_sum_controls)
    constant_mode_quotient_ranks = tuple(
        stacked(block.transpose(), ones_row, context=context).rank()
        - block.transpose().rank()
        for block in blocks
    )
    require(constant_mode_quotient_ranks == (0,) * 6,
            constant_mode_quotient_ranks)
    difference_seams = tuple(
        int(value) for value in (difference * ones_column).entries()
    )
    require(difference_seams == (0,) * 13, difference_seams)
    support = tuple(
        (row, column)
        for row in range(13) for column in range(13)
        if int(difference[row, column])
    )
    source = mat((ints(difference)[3], ints(difference)[9]), context)
    endpoint = mat(
        (ints(difference)[0], ints(difference)[6], ints(difference)[12]),
        context,
    )
    endpoint_mid = mat((ints(difference)[6],), context)
    ranks = (
        difference.rank(), source.rank(), endpoint.rank(),
        stacked(source, endpoint, context=context).rank(),
    )
    require(block_ranks == (11, 11, 12, 12, 11, 11), block_ranks)
    require(ranks == (5, 2, 3, 5), ranks)

    left_plus = kernel_basis(plus.transpose(), prime, context)
    left_minus_raw = kernel_basis(minus.transpose(), prime, context)
    require(len(left_plus) == len(left_minus_raw) == 1, "dual line dimensions")
    ell_plus = left_plus[0]
    ell_minus = tuple(reversed(ell_plus))
    require(mat((ell_minus,), context).rank() == 1, "minus nonzero")
    require((minus.transpose()
             * mat(tuple((value,) for value in ell_minus), context)).rank() == 0,
            "linked minus kernel")

    image_plus = mat((ell_plus,), context) * difference
    image_minus = mat((ell_minus,), context) * difference
    image = stacked(image_plus, image_minus, context=context)
    source_image = stacked(source, image, context=context)
    source_image_endpoint = stacked(source_image, endpoint, context=context)
    intersection_dimensions = (
        source.rank() + image.rank() - source_image.rank(),
        source_image.rank() + endpoint.rank() - source_image_endpoint.rank(),
    )
    require((image.rank(), source_image.rank(), source_image_endpoint.rank())
            == (2, 3, 5), "response plane ranks")
    require(intersection_dimensions == (1, 1), intersection_dimensions)
    require(stacked(source_image, endpoint_mid, context=context).rank() == 3,
            "D6 lies in source+image")
    require(stacked(source, endpoint_mid, context=context).rank() == 3
            and stacked(image, endpoint_mid, context=context).rank() == 3,
            "relative maps injective")

    marked_basis = stacked(image, endpoint_mid, context=context)
    pivots = pivot_columns(marked_basis)
    require(len(pivots) == 3, pivots)
    square = mat(tuple(tuple(int(marked_basis[row, column])
                             for column in pivots)
                       for row in range(3)), context)
    source_rows = ints(source)
    coordinate_rows = []
    for source_row in source_rows:
        target = mat((tuple(source_row[column] for column in pivots),), context)
        coordinate_rows.append(tuple(int(value) for value in
                                     (target * square.inv()).entries()))
    coordinates = tuple(coordinate_rows)
    a, b, c = coordinates[0]
    require(coordinates[1] == (b, a, c), coordinates)
    determinant = (a * a - b * b) % prime
    require(determinant == 58994564309246164003190 and c,
            (determinant, c))

    swap_minus_identity = mat(((-1 % prime, 1), (1, -1 % prime)), context)
    twisted_rank = swap_minus_identity.rank()
    require(twisted_rank == 1, twisted_rank)

    middle_union_rank = stacked(
        plus.transpose(), minus.transpose(), context=context
    ).rank()
    require(middle_union_rank == 13, middle_union_rank)

    record = (
        CERTIFICATE_SHA256,
        prime,
        block_ranks,
        constant_mode_quotient_ranks,
        difference_seams,
        len(support),
        ranks,
        image.rank(),
        source_image.rank(),
        intersection_dimensions,
        coordinates,
        determinant,
        twisted_rank,
        middle_union_rank,
    )
    semantic = digest(record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                (semantic, EXPECTED_SEMANTIC_SHA256))

    print("== independent FLINT audit: r5 rank-two relative response cospan ==")
    print(f"engine=python-flint.fmpz_mod_mat;certificate_sha256={CERTIFICATE_SHA256};prime={prime}")
    print(f"parent_block_ranks={block_ranks};constant_mode_quotient_ranks={constant_mode_quotient_ranks};D_seams={difference_seams}")
    print(f"middle_support={len(support)};D_S_E_direct_sum={ranks}")
    print(f"dual_contraction=(rank={image.rank()},source_union_rank={source_image.rank()},intersection_dimensions={intersection_dimensions})")
    print(f"minimal_endpoint_sidecar=(r0=6,source_plus_line_rank={stacked(source_image, endpoint_mid, context=context).rank()})")
    print(f"marked_coordinates={coordinates};det={determinant}")
    print(f"twisted_boundary_rank={twisted_rank};twisted_H1_dimension={2-twisted_rank};arc_stable_middle_rank={middle_union_rank};arc_stable_middle_dual=0")
    print(f"semantic_sha256={semantic}")
    print("verdict=ACCEPT exact digit-cycle/rank/intersection/transition/twisted-H1 claims by independent algebra engine")
    print("scope=finite exact good-reduction matrix audit only;no physical current,no D5 flux,no LRC14")
    print("PASS")


if __name__ == "__main__":
    main()
