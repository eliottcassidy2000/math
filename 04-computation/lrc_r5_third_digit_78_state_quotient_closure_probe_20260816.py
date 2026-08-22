#!/usr/bin/env python3
"""Extract the static ten-dimensional quotient in the r2 pointed-P4 bank.

The expensive endpoint computation is not part of the normal reproduction
path.  A compact exact bank, extracted once from the hash-pinned third-digit
probe, stores only the six diagonal K2 and K3 weights.  This script rebuilds
the audited 78-column parent/child rows from that bank and studies

    Q = F^78 / Row(parent),       dim Q = 10.

All statements are finite-field, static rowspace statements.  In particular,
an ambient coordinate completion is not called a chronological transfer.

Use ``--refresh-certificate`` only to regenerate the compact bank from the
audited endpoint pipeline.  Ordinary normal/-O verification reads the bank
and performs no endpoint integration.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import importlib
from itertools import combinations
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
PARENT_PATH = ROOT / "04-computation/lrc_r5_third_current_digit_pointed_root_difference_diagonal_bundle_probe_20260816.py"
CERTIFICATE_PATH = ROOT / "05-knowledge/results/lrc_r5_third_digit_78_state_exact_bank_20260816.json"
PARENT_RAW_SHA256 = "a227bc2f385d8a2eaecb27f317fa5ed66623c70938d8a97aba620298a8a7b61b"
PARENT_SEMANTIC_SHA256 = "3d1527fb4ce4931680e50d7135b9d1129c1816e3a9158645523e2728ddc71ec2"
PARENT_K2_DIGEST = "185b2fb843d37a6e6f73be48375e40e93029ff507392df16f0058892184a1db2"
PARENT_K3_MATRIX_DIGEST = "b5b0b3bceb926ab8176c0fb5c9fd57c4ad5fa87379b8e29731d7e7094f552263"
PARENT_K3_DIAGONAL_DIGEST = "1a3a8c73c62a9a7293ed9b80337df7211b45db31362b8bf1a5d223a6f584bec6"
EXPECTED_CERTIFICATE_SHA256 = "472925c638de7ac90b1a7880184766f2e8acec0f8291e2c46cf747b37cd46712"
EXPECTED_SEMANTIC_SHA256 = "2cfe4a03f00b431ce1542cb23d96650a34819050d8f430559b531253ed7cf874"

P = 13
ARCS = 6
SIZE = P * ARCS
ARC_REVERSAL = (1, 0, 3, 2, 5, 4)
CHAMBER_REVERSAL = (5, 4, 3, 2, 1, 0)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def raw_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_json(value: object) -> str:
    return sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")).hexdigest()


def load_parent_module():
    require(raw_sha256(PARENT_PATH) == PARENT_RAW_SHA256,
            ("third-digit parent source drift", raw_sha256(PARENT_PATH)))
    parent_directory = str(PARENT_PATH.parent)
    if parent_directory not in sys.path:
        sys.path.insert(0, parent_directory)
    # The real stem is intentional: ProcessPoolExecutor's Windows workers
    # must be able to import the module named on operator_worker.__module__.
    module = importlib.import_module(PARENT_PATH.stem)
    require(module.EXPECTED_SEMANTIC_SHA256 == PARENT_SEMANTIC_SHA256,
            "third-digit parent semantic drift")
    return module


def full_diagonal(diagonal, prime):
    return tuple(tuple(diagonal[row] if row == column else 0
                       for column in range(ARCS)) for row in range(ARCS))


def bank_digest_record(k2, k3, prime):
    k2_matrix_digests = tuple(
        digest_json(full_diagonal(tuple(k2[point][r0][r1]
                                        for point in range(ARCS)), prime))
        for r0 in range(P) for r1 in range(P)
    )
    k3_matrix_digests = tuple(
        digest_json(full_diagonal(tuple(k3[point][r0][r1][r2]
                                        for point in range(ARCS)), prime))
        for r0 in range(P) for r1 in range(P) for r2 in range(P)
    )
    k3_diagonal = tuple(
        tuple(tuple(tuple(k3[point][r0][r1][r2]
                          for point in range(ARCS))
                    for r2 in range(P))
              for r1 in range(P))
        for r0 in range(P)
    )
    return (
        digest_json(k2_matrix_digests),
        digest_json(k3_matrix_digests),
        digest_json(k3_diagonal),
    )


def refresh_certificate() -> None:
    """One-time extraction; deliberately excluded from normal reproduction."""

    parent_module = load_parent_module()
    _profiles, boundaries, cells, _record = parent_module.three_digit_context()
    operator, _operator_record = parent_module.build_sparse_operator(boundaries)
    parent, indexed, sources = parent_module.assemble_parent_tensor(operator, cells)
    require(parent_module.digest_json(parent) == parent_module.PARENT_TENSOR_SHA256,
            "parent tensor extraction drift")
    k2_full, _k2_record = parent_module.rebuild_k2(parent)
    k3_full, k3_record = parent_module.transition_bank(parent, indexed, sources, "actual")

    k2 = tuple(tuple(tuple(k2_full[r0][r1][point][point]
                              for r1 in range(P)) for r0 in range(P))
                     for point in range(ARCS))
    k3 = tuple(tuple(tuple(tuple(k3_full[r0][r1][r2][point][point]
                                    for r2 in range(P)) for r1 in range(P))
                           for r0 in range(P)) for point in range(ARCS))
    digests = bank_digest_record(k2, k3, parent_module.PRIME)
    require(digests == (PARENT_K2_DIGEST, PARENT_K3_MATRIX_DIGEST,
                        PARENT_K3_DIAGONAL_DIGEST),
            ("extracted bank digest", digests))
    require((k3_record[17], k3_record[18]) == digests[1:],
            ("transition record digest", k3_record[17:19], digests[1:]))
    certificate = {
        "schema": "r5-third-digit-78-diagonal-bank-v1",
        "parent_raw_sha256": PARENT_RAW_SHA256,
        "parent_semantic_sha256": PARENT_SEMANTIC_SHA256,
        "prime": parent_module.PRIME,
        "zeta13": parent_module.C.context()["zeta"],
        "points": parent_module.POINTS,
        "k2": k2,
        "k3": k3,
        "bank_digests": digests,
    }
    CERTIFICATE_PATH.write_text(
        json.dumps(certificate, sort_keys=True, separators=(",", ":")) + "\n",
        encoding="ascii", newline="\n",
    )
    print(f"wrote={CERTIFICATE_PATH}")
    print(f"certificate_lf_sha256={lf_sha256(CERTIFICATE_PATH)}")


def load_certificate():
    require(CERTIFICATE_PATH.is_file(), ("missing exact bank", CERTIFICATE_PATH))
    require(raw_sha256(PARENT_PATH) == PARENT_RAW_SHA256,
            ("third-digit parent source drift", raw_sha256(PARENT_PATH)))
    cert_sha = lf_sha256(CERTIFICATE_PATH)
    if EXPECTED_CERTIFICATE_SHA256 != "TO_BE_PINNED":
        require(cert_sha == EXPECTED_CERTIFICATE_SHA256,
                ("certificate drift", cert_sha, EXPECTED_CERTIFICATE_SHA256))
    certificate = json.loads(CERTIFICATE_PATH.read_text(encoding="ascii"))
    require(certificate["schema"] == "r5-third-digit-78-diagonal-bank-v1",
            "certificate schema")
    require(certificate["parent_raw_sha256"] == PARENT_RAW_SHA256,
            "certificate parent source")
    require(certificate["parent_semantic_sha256"] == PARENT_SEMANTIC_SHA256,
            "certificate parent semantic")
    prime = certificate["prime"]
    k2 = certificate["k2"]
    k3 = certificate["k3"]
    digests = bank_digest_record(k2, k3, prime)
    require(tuple(certificate["bank_digests"]) == digests,
            ("certificate self digest", certificate["bank_digests"], digests))
    require(digests == (PARENT_K2_DIGEST, PARENT_K3_MATRIX_DIGEST,
                        PARENT_K3_DIAGONAL_DIGEST),
            ("audited bank digest", digests))
    return certificate, cert_sha


def rref(rows, prime, columns=None):
    if columns is None:
        columns = len(rows[0]) if rows else 0
    work = [[value % prime for value in row] for row in rows]
    require(all(len(row) == columns for row in work), "ragged matrix")
    rank = 0
    pivots = []
    for column in range(columns):
        pivot = next((row for row in range(rank, len(work)) if work[row][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = pow(work[rank][column], -1, prime)
        work[rank] = [value * inverse % prime for value in work[rank]]
        for row in range(len(work)):
            if row == rank or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [(value - factor * pivot_value) % prime
                         for value, pivot_value in zip(work[row], work[rank])]
        pivots.append(column)
        rank += 1
        if rank == len(work):
            break
    return tuple(tuple(row) for row in work[:rank]), tuple(pivots)


def rank(rows, prime, columns=None):
    return len(rref(rows, prime, columns)[0])


def nullspace(rows, prime, columns=None):
    basis, pivots = rref(rows, prime, columns)
    if columns is None:
        columns = len(rows[0]) if rows else 0
    free = tuple(column for column in range(columns) if column not in pivots)
    answer = []
    for free_column in free:
        vector = [0] * columns
        vector[free_column] = 1
        for row, pivot in enumerate(pivots):
            vector[pivot] = -basis[row][free_column] % prime
        answer.append(tuple(vector))
    return tuple(answer)


def matmul(left, right, prime):
    if not left:
        return ()
    return tuple(tuple(sum(left[row][middle] * right[middle][column]
                           for middle in range(len(right))) % prime
                       for column in range(len(right[0])))
                 for row in range(len(left)))


def transpose(matrix):
    return tuple(tuple(matrix[row][column] for row in range(len(matrix)))
                 for column in range(len(matrix[0]))) if matrix else ()


def inverse(matrix, prime):
    size = len(matrix)
    work = [[value % prime for value in row]
            + [int(row_index == column) for column in range(size)]
            for row_index, row in enumerate(matrix)]
    require(all(len(row) == 2 * size for row in work), "inverse shape")
    for column in range(size):
        pivot = next((row for row in range(column, size) if work[row][column]), None)
        require(pivot is not None, ("singular inverse", column))
        work[column], work[pivot] = work[pivot], work[column]
        scale = pow(work[column][column], -1, prime)
        work[column] = [value * scale % prime for value in work[column]]
        for row in range(size):
            if row == column or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [(value - factor * pivot_value) % prime
                         for value, pivot_value in zip(work[row], work[column])]
    return tuple(tuple(row[size:]) for row in work)


def identity(size):
    return tuple(tuple(int(row == column) for column in range(size))
                 for row in range(size))


def matrix_difference(left, right):
    for row in range(len(left)):
        for column in range(len(left[row])):
            if left[row][column] != right[row][column]:
                return row, column, left[row][column], right[row][column]
    return None


def pair_rows(rows, annihilator, prime):
    return tuple(tuple(sum(row[column] * vector[column]
                           for column in range(SIZE)) % prime
                       for vector in annihilator) for row in rows)


def row_basis(rows, prime, columns=SIZE):
    return rref(rows, prime, columns)[0]


def intersection_basis(left, right, prime, columns=SIZE):
    left = rref(left, prime, columns)[0]
    right = rref(right, prime, columns)[0]
    left_perp = nullspace(left, prime, columns)
    right_perp = nullspace(right, prime, columns)
    return rref(nullspace(left_perp + right_perp, prime, columns),
                prime, columns)[0]


def canonical_section(parent_basis, child_rows, annihilator, prime):
    selected = []
    quotient_rows = []
    current_rank = 0
    child_quotient = pair_rows(child_rows, annihilator, prime)
    for index, qrow in enumerate(child_quotient):
        new_rank = rank(tuple(quotient_rows) + (qrow,), prime, len(annihilator))
        if new_rank > current_rank:
            selected.append(child_rows[index])
            quotient_rows.append(qrow)
            current_rank = new_rank
        if current_rank == len(annihilator):
            break
    require(current_rank == len(annihilator), ("child quotient rank", current_rank))
    selected = tuple(selected)
    qmatrix = pair_rows(selected, annihilator, prime)
    normalized = matmul(inverse(qmatrix, prime), selected, prime)
    require(pair_rows(normalized, annihilator, prime) == identity(len(annihilator)),
            "normalized quotient section")
    require(rank(parent_basis + normalized, prime, SIZE) == SIZE,
            "parent/section completion")
    return selected, normalized


def coordinate_permutation(point_map, digit_map):
    permutation = [None] * SIZE
    for point in range(ARCS):
        for digit in range(P):
            old = point * P + digit
            new = point_map[point] * P + digit_map[digit]
            permutation[old] = new
    require(sorted(permutation) == list(range(SIZE)), "coordinate permutation")
    return tuple(permutation)


def permute_row(row, permutation):
    answer = [0] * len(row)
    for old, new in enumerate(permutation):
        answer[new] = row[old]
    return tuple(answer)


def transform_rows(rows, permutation):
    return tuple(permute_row(row, permutation) for row in rows)


def invariance_record(parent_basis, permutation, annihilator, prime):
    transformed = transform_rows(parent_basis, permutation)
    union_rank = rank(parent_basis + transformed, prime, SIZE)
    first = None
    if union_rank != len(parent_basis):
        for index, row in enumerate(transformed):
            qrow = pair_rows((row,), annihilator, prime)[0]
            if any(qrow):
                first = (index, next((column, value) for column, value in enumerate(qrow)
                                     if value), digest_json(qrow))
                break
    return union_rank == len(parent_basis), union_rank, first


def induced_action(section, annihilator, permutation, parent_basis, prime):
    invariant, union_rank, first = invariance_record(
        parent_basis, permutation, annihilator, prime)
    if not invariant:
        return None, (invariant, union_rank, first)
    matrix = pair_rows(transform_rows(section, permutation), annihilator, prime)
    square = matmul(matrix, matrix, prime)
    require(square == identity(len(annihilator)),
            ("involution on quotient", matrix_difference(square, identity(len(annihilator)))))
    plus = tuple(tuple((matrix[row][column] - int(row == column)) % prime
                       for column in range(len(matrix))) for row in range(len(matrix)))
    minus = tuple(tuple((matrix[row][column] + int(row == column)) % prime
                        for column in range(len(matrix))) for row in range(len(matrix)))
    eigens = (len(matrix) - rank(plus, prime, len(matrix)),
              len(matrix) - rank(minus, prime, len(matrix)))
    return matrix, (invariant, union_rank, first, eigens, digest_json(matrix))


def joint_eigenspaces(left, right, prime):
    require(matmul(left, right, prime) == matmul(right, left, prime),
            "quotient involutions commute")
    size = len(left)
    result = []
    for left_sign in (1, -1):
        for right_sign in (1, -1):
            left_equations = transpose(tuple(
                tuple((left[row][column] - left_sign * int(row == column)) % prime
                      for column in range(size)) for row in range(size)))
            right_equations = transpose(tuple(
                tuple((right[row][column] - right_sign * int(row == column)) % prime
                      for column in range(size)) for row in range(size)))
            result.append(((left_sign, right_sign),
                           size - rank(left_equations + right_equations, prime, size)))
    require(sum(value for _signs, value in result) == size,
            ("joint eigenspace dimensions", result))
    return tuple(result)


def histogram(values):
    counts = {}
    for value in values:
        counts[value] = counts.get(value, 0) + 1
    return tuple(sorted(counts.items()))


def build_rows(k2, k3, prime):
    parent = []
    children = [[] for _r2 in range(P)]
    for point in range(ARCS):
        for r1 in range(P):
            parent_row = [0] * SIZE
            child_rows = [[0] * SIZE for _r2 in range(P)]
            for r0 in range(P):
                value = k2[point][r0][r1] % prime
                column = point * P + r0
                parent_row[column] = value
                for r2 in range(P):
                    child_rows[r2][column] = value * k3[point][r0][r1][r2] % prime
            parent.append(tuple(parent_row))
            for r2 in range(P):
                children[r2].append(tuple(child_rows[r2]))
    return tuple(parent), tuple(tuple(rows) for rows in children)


def block_rank_record(parent_rows, prime):
    return tuple(rank(tuple(parent_rows[point * P:(point + 1) * P]), prime, SIZE)
                 for point in range(ARCS))


def fourier_vectors(mode, zeta, prime):
    pure = []
    for point in range(ARCS):
        row = [0] * SIZE
        for digit in range(P):
            row[point * P + digit] = pow(zeta, mode * digit, prime)
        pure.append(tuple(row))
    return tuple(pure)


def fourier_record(annihilator, zeta, prime):
    mode_rows = []
    records = []
    cumulative = []
    for mode in range(P):
        pure = fourier_vectors(mode, zeta, prime)
        quotient = pair_rows(pure, annihilator, prime)
        even = []
        odd = []
        for left, right in ((0, 1), (2, 3), (4, 5)):
            even.append(tuple((pure[left][column] + pure[right][column]) % prime
                              for column in range(SIZE)))
            odd.append(tuple((pure[left][column] - pure[right][column]) % prime
                             for column in range(SIZE)))
        q_even = pair_rows(tuple(even), annihilator, prime)
        q_odd = pair_rows(tuple(odd), annihilator, prime)
        cumulative.extend(quotient)
        mode_rows.append(quotient)
        records.append((mode, rank(quotient, prime, len(annihilator)),
                        rank(q_even, prime, len(annihilator)),
                        rank(q_odd, prime, len(annihilator)),
                        rank(tuple(cumulative), prime, len(annihilator))))
    pair_records = tuple((mode, rank(mode_rows[mode] + mode_rows[-mode % P],
                                    prime, len(annihilator)))
                         for mode in range((P + 1) // 2))
    return tuple(records), pair_records


def chamber_fourier_record(annihilator, chamber_permutation, zeta, prime):
    """Pair m and -m into exact chamber-reflection eigenvectors."""

    records = []
    for mode in range(1, (P + 1) // 2):
        positive = fourier_vectors(mode, zeta, prime)
        negative = fourier_vectors(-mode % P, zeta, prime)
        phase = pow(zeta, -mode, prime)
        signs = []
        for sign in (1, -1):
            vectors = []
            for point in range(ARCS):
                reflected = CHAMBER_REVERSAL[point]
                vector = tuple(
                    (positive[point][column]
                     + sign * phase * negative[reflected][column]) % prime
                    for column in range(SIZE)
                )
                transformed = permute_row(vector, chamber_permutation)
                expected = tuple(sign * value % prime for value in vector)
                require(transformed == expected,
                        ("chamber Fourier eigenvector", mode, sign, point))
                vectors.append(vector)
            signs.append(rank(pair_rows(tuple(vectors), annihilator, prime),
                              prime, len(annihilator)))
        records.append((mode, signs[0], signs[1], sum(signs)))
    return tuple(records)


def minimal_child_spanning_set(children, prime):
    for size in range(1, P + 1):
        witnesses = []
        for subset in combinations(range(P), size):
            rows = tuple(row for digit in subset for row in children[digit])
            if rank(rows, prime, SIZE) == SIZE:
                witnesses.append(subset)
        if witnesses:
            return size, len(witnesses), tuple(witnesses), digest_json(tuple(witnesses))
    return None


def first_containment_failure(children, prime):
    contained = 0
    first = None
    pair_union_ranks = []
    for source in range(P):
        source_rank = rank(children[source], prime, SIZE)
        for target in range(P):
            target_rank = rank(children[target], prime, SIZE)
            union_rank = rank(children[source] + children[target], prime, SIZE)
            if union_rank == source_rank:
                contained += 1
            elif first is None:
                first = (source, target, source_rank, target_rank, union_rank)
            if source < target:
                pair_union_ranks.append(union_rank)
    return contained, P * P, first, histogram(tuple(pair_union_ranks))


def analyse(certificate, cert_sha):
    prime = certificate["prime"]
    zeta = certificate["zeta13"]
    require(zeta != 1 and pow(zeta, P, prime) == 1,
            ("primitive thirteenth root", zeta, prime))
    k2 = certificate["k2"]
    k3 = certificate["k3"]
    parent_rows, children = build_rows(k2, k3, prime)
    parent_basis = row_basis(parent_rows, prime)
    require(len(parent_basis) == 68, ("parent rank", len(parent_basis)))
    annihilator = nullspace(parent_basis, prime, SIZE)
    require(len(annihilator) == 10, ("quotient dimension", len(annihilator)))
    block_ranks = block_rank_record(parent_rows, prime)
    require(block_ranks == (11, 11, 12, 12, 11, 11),
            ("parent block ranks", block_ranks))

    child_ranks = tuple(rank(child, prime, SIZE) for child in children)
    union_ranks = tuple(rank(parent_basis + child, prime, SIZE) for child in children)
    quotient_ranks = tuple(rank(pair_rows(child, annihilator, prime), prime, 10)
                           for child in children)
    require(union_ranks == (78,) * P and quotient_ranks == (10,) * P,
            ("fixed-r2 quotient", union_ranks, quotient_ranks))
    block_codimensions = tuple(P - value for value in block_ranks)
    child_block_quotients = []
    for child in children:
        profile = []
        for point in range(ARCS):
            parent_block = tuple(parent_rows[point * P:(point + 1) * P])
            child_block = tuple(child[point * P:(point + 1) * P])
            profile.append(rank(parent_block + child_block, prime, SIZE)
                           - block_ranks[point])
        child_block_quotients.append(tuple(profile))
    require(tuple(child_block_quotients) == (block_codimensions,) * P,
            ("fixed-r2 block quotient", child_block_quotients))

    sections = []
    normalized = []
    for child in children:
        selected, normalized_section = canonical_section(
            parent_basis, child, annihilator, prime)
        sections.append(selected)
        normalized.append(normalized_section)
    sections = tuple(sections)
    normalized = tuple(normalized)
    section_equal_pairs = 0
    section_union_ranks = []
    first_section_difference = None
    for left in range(P):
        for right in range(left + 1, P):
            union_rank = rank(sections[left] + sections[right], prime, SIZE)
            section_union_ranks.append(union_rank)
            if union_rank == 10:
                section_equal_pairs += 1
            elif first_section_difference is None:
                first_section_difference = (left, right, union_rank)

    gauge_ranks = []
    gauge_rows = []
    first_gauge = None
    for r2 in range(P):
        differences = tuple(tuple((normalized[r2][row][column]
                                   - normalized[0][row][column]) % prime
                                  for column in range(SIZE)) for row in range(10))
        require(all(not any(qrow) for qrow in pair_rows(differences, annihilator, prime)),
                ("section gauge leaves quotient", r2))
        gauge_rank = rank(differences, prime, SIZE)
        gauge_ranks.append(gauge_rank)
        gauge_rows.extend(differences)
        if gauge_rank and first_gauge is None:
            first_gauge = (r2, gauge_rank, digest_json(differences))

    common_child = row_basis(children[0], prime)
    common_section = row_basis(sections[0], prime)
    for r2 in range(1, P):
        common_child = intersection_basis(common_child, children[r2], prime)
        common_section = intersection_basis(common_section, sections[r2], prime)
    common_child_quotient_rank = rank(pair_rows(common_child, annihilator, prime), prime, 10)
    common_section_quotient_rank = rank(pair_rows(common_section, annihilator, prime), prime, 10)

    arc_perm = coordinate_permutation(ARC_REVERSAL, tuple(range(P)))
    digit_perm = coordinate_permutation(tuple(range(ARCS)), tuple(P - 1 - digit for digit in range(P)))
    chamber_perm = coordinate_permutation(CHAMBER_REVERSAL,
                                           tuple(P - 1 - digit for digit in range(P)))
    shift_perm = coordinate_permutation(tuple(range(ARCS)),
                                        tuple((digit + 1) % P for digit in range(P)))
    arc_matrix, arc_record = induced_action(
        normalized[0], annihilator, arc_perm, parent_basis, prime)
    digit_matrix, digit_record = induced_action(
        normalized[0], annihilator, digit_perm, parent_basis, prime)
    chamber_matrix, chamber_record = induced_action(
        normalized[0], annihilator, chamber_perm, parent_basis, prime)
    shift_record = invariance_record(parent_basis, shift_perm, annihilator, prime)
    joint = None
    if arc_matrix is not None and chamber_matrix is not None:
        joint = joint_eigenspaces(arc_matrix, chamber_matrix, prime)

    arc_pair_obstruction = []
    for pair in ((0, 1), (2, 3), (4, 5)):
        pair_parent = tuple(
            parent_rows[point * P + row]
            for point in pair for row in range(P)
        )
        pair_rank = rank(pair_parent, prime, SIZE)
        pair_union = rank(pair_parent + transform_rows(pair_parent, arc_perm),
                          prime, SIZE)
        arc_pair_obstruction.append((pair, pair_rank, pair_union,
                                     pair_union - pair_rank))

    # Arc reversal misses Q by exactly the two middle-root directions.  The
    # smallest repaired quotient kills that defect and has dimension eight.
    arc_stable_parent = row_basis(
        parent_basis + transform_rows(parent_basis, arc_perm), prime)
    arc_stable_annihilator = nullspace(arc_stable_parent, prime, SIZE)
    arc_stable_child_quotient_ranks = tuple(
        rank(pair_rows(child, arc_stable_annihilator, prime),
             prime, len(arc_stable_annihilator)) for child in children
    )
    require(len(arc_stable_parent) == 70
            and len(arc_stable_annihilator) == 8
            and arc_stable_child_quotient_ranks == (8,) * P,
            ("arc-stable quotient", len(arc_stable_parent),
             len(arc_stable_annihilator), arc_stable_child_quotient_ranks))
    _arc_selected, arc_stable_section = canonical_section(
        arc_stable_parent, children[0], arc_stable_annihilator, prime)
    arc_stable_matrix, arc_stable_record = induced_action(
        arc_stable_section, arc_stable_annihilator, arc_perm,
        arc_stable_parent, prime)
    chamber_stable_matrix, chamber_stable_record = induced_action(
        arc_stable_section, arc_stable_annihilator, chamber_perm,
        arc_stable_parent, prime)
    arc_stable_joint = joint_eigenspaces(
        arc_stable_matrix, chamber_stable_matrix, prime)

    chamber_child_matches = []
    first_chamber_child_failure = None
    arc_child_matches = []
    first_arc_child_failure = None
    for r2 in range(P):
        reflected = P - 1 - r2
        chamber_union = rank(transform_rows(children[r2], chamber_perm)
                             + children[reflected], prime, SIZE)
        chamber_match = (chamber_union == child_ranks[r2] == child_ranks[reflected])
        chamber_child_matches.append(chamber_match)
        if not chamber_match and first_chamber_child_failure is None:
            first_chamber_child_failure = (r2, reflected, child_ranks[r2],
                                           child_ranks[reflected], chamber_union)
        arc_union = rank(transform_rows(children[r2], arc_perm)
                         + children[r2], prime, SIZE)
        arc_match = arc_union == child_ranks[r2]
        arc_child_matches.append(arc_match)
        if not arc_match and first_arc_child_failure is None:
            first_arc_child_failure = (r2, child_ranks[r2], arc_union)

    fourier_modes, fourier_pairs = fourier_record(annihilator, zeta, prime)
    chamber_fourier = chamber_fourier_record(
        annihilator, chamber_perm, zeta, prime)
    arc_stable_fourier, arc_stable_fourier_pairs = fourier_record(
        arc_stable_annihilator, zeta, prime)
    containment = first_containment_failure(children, prime)
    cumulative_child_ranks = []
    accumulated = ()
    for child in children:
        accumulated = row_basis(accumulated + child, prime)
        cumulative_child_ranks.append(len(accumulated))
    minimal_spanning_set = minimal_child_spanning_set(children, prime)

    # A fixed parent+section basis gives unique static coordinates for every
    # child, but the child snapshots themselves are rank-deficient and do not
    # generally carry one another.
    common_basis = parent_basis + normalized[0]
    require(rank(common_basis, prime, SIZE) == SIZE, "common ambient basis")
    common_inverse = inverse(common_basis, prime)
    coordinate_maps = tuple(matmul(child, common_inverse, prime) for child in children)
    require(all(matmul(coordinates, common_basis, prime) == child
                for coordinates, child in zip(coordinate_maps, children)),
            "common static coordinate maps")
    child_plus_fixed_section_ranks = tuple(
        rank(child + normalized[0], prime, SIZE) for child in children
    )

    common_child_quotient = row_basis(
        pair_rows(common_child, annihilator, prime), prime, 10)
    common_child_chamber = tuple(
        tuple(sum(row[middle] * chamber_matrix[middle][column]
                  for middle in range(10)) % prime
              for column in range(10))
        for row in common_child_quotient
    )
    common_child_chamber_invariant = (
        rank(common_child_quotient + common_child_chamber, prime, 10)
        == len(common_child_quotient)
    )
    chamber_plus = nullspace(transpose(tuple(
        tuple((chamber_matrix[row][column] - int(row == column)) % prime
              for column in range(10)) for row in range(10))), prime, 10)
    chamber_minus = nullspace(transpose(tuple(
        tuple((chamber_matrix[row][column] + int(row == column)) % prime
              for column in range(10)) for row in range(10))), prime, 10)
    common_child_chamber_record = (
        common_child_chamber_invariant,
        len(intersection_basis(common_child_quotient, chamber_plus, prime, 10)),
        len(intersection_basis(common_child_quotient, chamber_minus, prime, 10)),
    )

    # Flat and support-normalized controls coincide in the parent extraction:
    # each active r2 fibre has all 13 digits.  Their child rows are P/13, so
    # they add no quotient direction.
    inverse13 = pow(P, -1, prime)
    hostile_child = tuple(tuple(value * inverse13 % prime for value in row)
                          for row in parent_rows)
    hostile = (
        rank(hostile_child, prime, SIZE),
        rank(parent_basis + hostile_child, prime, SIZE),
        rank(pair_rows(hostile_child, annihilator, prime), prime, 10),
        "support_equals_flat_due_full_13_support",
    )
    require(hostile[:3] == (68, 68, 0), ("flat/support hostile", hostile))

    record = (
        PARENT_RAW_SHA256, PARENT_SEMANTIC_SHA256, cert_sha,
        prime, zeta, tuple(tuple(point) for point in certificate["points"]),
        block_ranks, len(parent_basis), len(annihilator),
        child_ranks, union_ranks, quotient_ranks,
        block_codimensions, tuple(child_block_quotients),
        section_equal_pairs, histogram(tuple(section_union_ranks)),
        first_section_difference, tuple(gauge_ranks),
        rank(tuple(gauge_rows), prime, SIZE), first_gauge,
        (len(common_child), common_child_quotient_rank,
         len(common_section), common_section_quotient_rank),
        arc_record, tuple(arc_pair_obstruction), digit_record,
        chamber_record, shift_record, joint,
        (len(arc_stable_parent), len(arc_stable_annihilator),
         arc_stable_record, chamber_stable_record, arc_stable_joint,
         arc_stable_child_quotient_ranks),
        (sum(chamber_child_matches), first_chamber_child_failure,
         sum(arc_child_matches), first_arc_child_failure),
        fourier_modes, fourier_pairs, chamber_fourier,
        arc_stable_fourier, arc_stable_fourier_pairs,
        common_child_chamber_record, containment,
        tuple(cumulative_child_ranks), minimal_spanning_set,
        child_plus_fixed_section_ranks,
        tuple((rank(matrix, prime, SIZE), digest_json(matrix))
              for matrix in coordinate_maps),
        hostile,
        digest_json(parent_basis), digest_json(annihilator),
        digest_json(tuple(digest_json(child) for child in children)),
        digest_json(tuple(digest_json(section) for section in normalized)),
    )
    semantic = digest_json(record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("semantic drift", semantic, EXPECTED_SEMANTIC_SHA256))
    return record, semantic


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--refresh-certificate", action="store_true")
    args = parser.parse_args()
    if args.refresh_certificate:
        refresh_certificate()
        return

    certificate, cert_sha = load_certificate()
    record, semantic = analyse(certificate, cert_sha)
    (
        _parent_raw, _parent_semantic, _cert_sha, prime, zeta, points,
        block_ranks, parent_rank, quotient_dimension,
        child_ranks, union_ranks, quotient_ranks,
        block_codimensions, child_block_quotients,
        section_equal_pairs, section_union_hist, first_section_difference,
        gauge_ranks, gauge_span_rank, first_gauge,
        common_intersections,
        arc_record, arc_pair_obstruction, digit_record, chamber_record,
        shift_record, joint, arc_stable_record,
        child_symmetry, fourier_modes, fourier_pairs, chamber_fourier,
        arc_stable_fourier, arc_stable_fourier_pairs,
        common_child_chamber, containment,
        cumulative_child_ranks, minimal_spanning_set,
        child_plus_fixed_section_ranks, coordinate_maps, hostile,
        parent_basis_digest, annihilator_digest, child_digest, section_digest,
    ) = record
    print("== r=5 third-digit 78-state quotient closure extraction ==")
    print(f"parent=(raw_sha256={PARENT_RAW_SHA256},semantic_sha256={PARENT_SEMANTIC_SHA256})")
    print(f"certificate=(path={CERTIFICATE_PATH.name},sha256={cert_sha},bank_digests={(PARENT_K2_DIGEST, PARENT_K3_MATRIX_DIGEST, PARENT_K3_DIAGONAL_DIGEST)})")
    print(f"field=(prime={prime},zeta13={zeta});points={points}")
    print(f"canonical_quotient=(ambient=78,parent_rank={parent_rank},dimension={quotient_dimension},block_ranks={block_ranks},block_codimensions={block_codimensions},parent_basis_digest={parent_basis_digest},annihilator_digest={annihilator_digest})")
    print(f"fixed_r2=(child_ranks={child_ranks},parent_child_union_ranks={union_ranks},quotient_image_ranks={quotient_ranks},per_arc_quotient_profile={child_block_quotients[0]},same_for_all_r2={len(set(child_block_quotients)) == 1},child_digest={child_digest})")
    print(f"quotient_sections=(equal_pairs={section_equal_pairs}/78,union_rank_hist={section_union_hist},first_lift_difference={first_section_difference},normalized_gauge_ranks={gauge_ranks},gauge_span_rank={gauge_span_rank},first_gauge={first_gauge},section_digest={section_digest})")
    print(f"common_intersections=(all_child_dimension,quotient_rank,all_section_dimension,quotient_rank)={common_intersections}")
    print("symmetry_record=(invariant,parent_plus_transform_rank,first_failure[,quotient_(plus,minus),matrix_digest])")
    print(f"arc_reversal={arc_record}")
    print(f"arc_reversal_localization=(pair,parent_rank,stabilized_rank,defect)={arc_pair_obstruction}")
    print(f"pure_digit_reflection={digit_record}")
    print(f"chamber_digit_reflection={chamber_record}")
    print(f"digit_shift={shift_record}")
    print(f"joint_arc_chamber_characters={joint}")
    print("arc_stable_quotient=(stabilized_parent_rank,dimension,arc_record,chamber_record,joint_characters,child_image_ranks)")
    print(f"arc_stable_quotient={arc_stable_record}")
    print(f"child_symmetry=(chamber_matches,first_failure,arc_matches,first_failure)={child_symmetry}")
    print("fourier_mode_record=(mode,quotient_rank,arc_even_rank,arc_odd_rank,cumulative_rank)")
    print(f"digit_fourier_modes={fourier_modes}")
    print(f"digit_fourier_reflection_pairs=(mode,rank_of_m_and_minus_m)={fourier_pairs}")
    print(f"chamber_fourier_pair_eigenspaces=(mode,plus_rank,minus_rank,total)={chamber_fourier}")
    print(f"arc_stable_fourier_modes={arc_stable_fourier}")
    print(f"arc_stable_fourier_pairs={arc_stable_fourier_pairs}")
    print(f"common_child_quotient_chamber=(invariant,plus_rank,minus_rank)={common_child_chamber}")
    print(f"child_to_child=(containments,total,first_failure,pair_union_rank_hist)={containment}")
    print(f"cumulative_child_union_ranks={cumulative_child_ranks}")
    print(f"minimal_child_family_spanning_78=(size,count,witnesses,digest)={minimal_spanning_set}")
    print(f"child_plus_fixed_quotient_section_ranks={child_plus_fixed_section_ranks}")
    print(f"common_static_coordinate_maps=(rank,digest)={coordinate_maps}")
    print(f"flat_support_hostile=(child_rank,parent_union_rank,quotient_rank,note)={hostile}")
    print(f"semantic_sha256={semantic}")
    print("interpretation=Q=F^78/Row(parent) is the same abstract ten-space for every r2, but canonical child lifts differ by parent-row gauges; a parent-plus-section basis is only a static ambient coordinatization")
    print("transfer_gate=child-to-child containment is required for a child-only linear transfer; failure is not repaired by naming the common quotient")
    print("scope=FINITE-EXACT static rowspace quotient on the audited pointed-P4 bank;no chronology,no complete address,no current,no H1,no ancestry,no LRC14")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
