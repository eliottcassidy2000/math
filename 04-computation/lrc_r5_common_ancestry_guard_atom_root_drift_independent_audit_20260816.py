#!/usr/bin/env python3
"""Independent audit of the r=5 common-ancestry guard-atom package.

The submitted bucket implementation is hash-pinned but never imported.  This
verifier rebuilds the THM-3514 39-atom partition on the two THM-2471/2594
linked-node inputs before P_(13^5), uses an independent interval sweep and
product integrator, and constructs the root/address/common-offset tensors
directly from atom/root labels.

Owner/drift Fourier support is decided exactly in Z[zeta_13] by reducing a
cyclic coefficient vector modulo Phi_13; five split-prime evaluations are a
second certificate.  Scope is support and spectral capacity only.
"""

from __future__ import annotations

from bisect import bisect_right
from fractions import Fraction
from hashlib import sha256
from itertools import permutations
import importlib.util
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
P = 13
STAGE_PATH = ROOT / "04-computation/lrc14_stage2_theta_contraction_opus_20260728.py"
CANDIDATE_PATH = ROOT / "04-computation/lrc_r5_common_ancestry_guard_atom_root_drift_probe_20260816.py"
CANDIDATE_OUTPUT = ROOT / "05-knowledge/results/lrc_r5_common_ancestry_guard_atom_root_drift_probe_20260816.out"

STAGE_SHA256 = "09c43af0a0a56c7a0833bbfd13ed6a96bc5a7a3718aa1bc6b77a144bde101a06"
CANDIDATE_SHA256 = "83f1fa49ac4d02e21a1d76fed169d101715a6620342714ed05b9172ae967a730"
CANDIDATE_OUTPUT_SHA256 = "7727571722f69a9c59af0183c50c8bae0b7944d6acb7bbb1b3c0e9bd02d54565"
EXPECTED_PAIR_GAUGE_SHA256 = "e4d0d4fa674e1f54496e613f7a3e1f057af033fa8992322a5f414ea176e1c3d4"
EXPECTED_SEMANTIC_SHA256 = "ad81e945207703956f5d6ec300430d562c9f98a9ec8788011119a4857ab34e01"

CHAMBERS = ("left", "middle", "right")
ACTIVE = ("left", "right")
CHAMBER_INDEX = {chamber: index for index, chamber in enumerate(CHAMBERS)}
ACTIVE_PAIRS = (("left", "left"), ("left", "right"), ("right", "left"), ("right", "right"))
WALSH_SIGNS = (
    (1, 1, 1, 1),
    (1, 1, -1, -1),
    (1, -1, 1, -1),
    (1, -1, -1, 1),
)
SPLIT_PRIMES = (547, 911, 1093, 2003, 2549)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_json(value: object) -> str:
    body = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    return sha256(body).hexdigest()


def load_stage():
    require(lf_sha256(STAGE_PATH) == STAGE_SHA256, "THM-2594 source drift")
    name = "thm2594_for_common_ancestry_independent_audit"
    spec = importlib.util.spec_from_file_location(name, STAGE_PATH)
    require(spec is not None and spec.loader is not None, "stage loader")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


M = load_stage()
GRID = M.T_DEN
require(GRID % 91 == 0, "guard partition not resolved")
ATOM_UNIT = GRID // 91
ATOMS = tuple((sheet, chamber) for sheet in range(P) for chamber in CHAMBERS)
ATOM_INDEX = {atom: index for index, atom in enumerate(ATOMS)}


def atom_intervals() -> tuple[tuple[int, int], ...]:
    bounds = ((0, 1), (1, 6), (6, 7))
    intervals = tuple(
        ((7 * sheet + low) * ATOM_UNIT, (7 * sheet + high) * ATOM_UNIT)
        for sheet in range(P)
        for low, high in bounds
    )
    require(intervals[0][0] == 0 and intervals[-1][1] == GRID, "atom endpoints")
    require(all(left[1] == right[0] for left, right in zip(intervals, intervals[1:])), "atom gaps")
    return intervals


ATOM_INTERVALS = atom_intervals()


def profile_pieces(starts: list[int], values: list[int]) -> list[tuple[int, int, int]]:
    pieces = []
    for index, left in enumerate(starts):
        right = starts[index + 1] if index + 1 < len(starts) else GRID
        if values[index]:
            pieces.append((left, right, int(values[index])))
    return pieces


def profile_value(starts: list[int], values: list[int], point: int) -> int:
    return int(values[bisect_right(starts, point) - 1])


def profile_mass(starts: list[int], values: list[int]) -> int:
    return sum(
        int(values[index])
        * ((starts[index + 1] if index + 1 < len(starts) else GRID) - left)
        for index, left in enumerate(starts)
    )


def intersect_weighted_with_intervals(
    pieces: list[tuple[int, int, int]], intervals: list[tuple[int, int]]
) -> list[tuple[int, int, int]]:
    """Independent two-pointer intersection of two sorted disjoint families."""
    answer = []
    i = j = 0
    while i < len(pieces) and j < len(intervals):
        left, right, weight = pieces[i]
        cut_left, cut_right = intervals[j]
        overlap_left = max(left, cut_left)
        overlap_right = min(right, cut_right)
        if overlap_left < overlap_right and weight:
            answer.append((overlap_left, overlap_right, weight))
        if right <= cut_right:
            i += 1
        if cut_right <= right:
            j += 1
    return answer


def build_f_pieces(
    e_intervals: list[tuple[int, int]], q_intervals: list[tuple[int, int]]
) -> list[tuple[int, int, int]]:
    p2_starts, p2_values = M.weighted_fold(
        [(left, right, 1) for left, right in e_intervals], M.RPKT, GRID
    )
    return intersect_weighted_with_intervals(
        profile_pieces(p2_starts, p2_values), q_intervals
    )


def split_by_atoms(
    pieces: list[tuple[int, int, int]],
) -> tuple[list[tuple[int, int, int]], ...]:
    """Split before P_(13^5), using a forward atom sweep rather than bisecting each piece."""
    groups: list[list[tuple[int, int, int]]] = [[] for _ in ATOMS]
    atom = 0
    input_mass = 0
    output_mass = 0
    previous_right = 0
    for left, right, weight in pieces:
        require(0 <= left < right <= GRID and weight >= 0, (left, right, weight))
        require(left >= previous_right, ("piece overlap", previous_right, left))
        previous_right = right
        input_mass += weight * (right - left)
        while ATOM_INTERVALS[atom][1] <= left:
            atom += 1
        cursor = left
        local_atom = atom
        while cursor < right:
            stop = min(right, ATOM_INTERVALS[local_atom][1])
            if weight:
                groups[local_atom].append((cursor, stop, weight))
                output_mass += weight * (stop - cursor)
            cursor = stop
            if cursor < right:
                local_atom += 1
    require(input_mass == output_mass, ("atom partition mass", input_mass, output_mass))
    return tuple(groups)


def verify_profile_sum(
    parts: tuple[tuple[list[int], list[int]], ...],
    whole: tuple[list[int], list[int]],
    label: str,
) -> None:
    points = {0}
    points.update(whole[0])
    for starts, _values in parts:
        points.update(starts)
    for point in sorted(points):
        if point >= GRID:
            continue
        actual = sum(profile_value(starts, values, point) for starts, values in parts)
        expected = profile_value(whole[0], whole[1], point)
        require(actual == expected, (label, point, actual, expected))


def product_integral(
    left: tuple[list[int], list[int]], right: tuple[list[int], list[int]]
) -> int:
    """Integral numerator of two full-grid step profiles, independently merged."""
    left_starts, left_values = left
    right_starts, right_values = right
    require(left_starts[0] == right_starts[0] == 0, "profile origin")
    i = j = 0
    cursor = 0
    total = 0
    while cursor < GRID:
        left_stop = left_starts[i + 1] if i + 1 < len(left_starts) else GRID
        right_stop = right_starts[j + 1] if j + 1 < len(right_starts) else GRID
        stop = min(left_stop, right_stop)
        total += int(left_values[i]) * int(right_values[j]) * (stop - cursor)
        cursor = stop
        if cursor == left_stop and i + 1 < len(left_starts):
            i += 1
        if cursor == right_stop and j + 1 < len(right_starts):
            j += 1
    return total


def rank_mod(matrix: tuple[tuple[int, ...], ...], prime: int) -> int:
    rows = [list(value % prime for value in row) for row in matrix]
    rank = 0
    columns = len(rows[0]) if rows else 0
    for column in range(columns):
        pivot = next((row for row in range(rank, len(rows)) if rows[row][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][column], -1, prime)
        rows[rank] = [value * inverse % prime for value in rows[rank]]
        for row in range(len(rows)):
            if row == rank:
                continue
            factor = rows[row][column]
            if factor:
                rows[row] = [
                    (value - factor * pivot_value) % prime
                    for value, pivot_value in zip(rows[row], rows[rank])
                ]
        rank += 1
        if rank == len(rows):
            break
    return rank


def rank_q(matrix: tuple[tuple[int, ...], ...]) -> int:
    rows = [list(Fraction(value) for value in row) for row in matrix]
    rank = 0
    columns = len(rows[0]) if rows else 0
    for column in range(columns):
        pivot = next((row for row in range(rank, len(rows)) if rows[row][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        pivot_value = rows[rank][column]
        rows[rank] = [value / pivot_value for value in rows[rank]]
        for row in range(rank + 1, len(rows)):
            factor = rows[row][column]
            if factor:
                rows[row] = [
                    value - factor * pivot_entry
                    for value, pivot_entry in zip(rows[row], rows[rank])
                ]
        rank += 1
        if rank == len(rows):
            break
    return rank


def pivot_columns(matrix: tuple[tuple[int, ...], ...], prime: int) -> tuple[int, ...]:
    rows = [list(value % prime for value in row) for row in matrix]
    pivots = []
    pivot_row = 0
    for column in range(len(rows[0])):
        pivot = next((row for row in range(pivot_row, len(rows)) if rows[row][column]), None)
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        inverse = pow(rows[pivot_row][column], -1, prime)
        rows[pivot_row] = [value * inverse % prime for value in rows[pivot_row]]
        for row in range(len(rows)):
            if row == pivot_row:
                continue
            factor = rows[row][column]
            if factor:
                rows[row] = [
                    (value - factor * pivot_entry) % prime
                    for value, pivot_entry in zip(rows[row], rows[pivot_row])
                ]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == len(rows):
            break
    return tuple(pivots)


def primitive_root(prime: int) -> int:
    order = prime - 1
    factors = []
    value = order
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            factors.append(divisor)
            while value % divisor == 0:
                value //= divisor
        divisor += 1
    if value > 1:
        factors.append(value)
    for candidate in range(2, prime):
        if all(pow(candidate, order // factor, prime) != 1 for factor in factors):
            return candidate
    raise RuntimeError(("primitive root", prime))


def zeta13_mod(prime: int) -> int:
    require((prime - 1) % P == 0, prime)
    root = pow(primitive_root(prime), (prime - 1) // P, prime)
    require(pow(root, P, prime) == 1 and root != 1, (prime, root))
    return root


def cyclic_add(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(a + b for a, b in zip(left, right))


def cyclic_scale(poly: tuple[int, ...], scalar: int) -> tuple[int, ...]:
    return tuple(scalar * value for value in poly)


def cyclic_mul(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    answer = [0] * P
    for i, a in enumerate(left):
        if not a:
            continue
        for j, b in enumerate(right):
            if b:
                answer[(i + j) % P] += a * b
    return tuple(answer)


def cyclotomic_zero(poly: tuple[int, ...]) -> bool:
    """Zero at zeta_13 iff a cyclic representative is a multiple of Phi_13."""
    require(len(poly) == P, len(poly))
    return len(set(poly)) == 1


def eval_poly(poly: tuple[int, ...], root: int, prime: int) -> int:
    return sum(value % prime * pow(root, exponent, prime) for exponent, value in enumerate(poly)) % prime


def permutation_sign(permutation: tuple[int, ...]) -> int:
    inversions = sum(
        permutation[i] > permutation[j]
        for i in range(len(permutation))
        for j in range(i + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def determinant_poly(matrix: tuple[tuple[tuple[int, ...], ...], ...]) -> tuple[int, ...]:
    require(len(matrix) == 4 and all(len(row) == 4 for row in matrix), "4x4 polynomial matrix")
    answer = (0,) * P
    one = (1,) + (0,) * (P - 1)
    for permutation in permutations(range(4)):
        term = one
        for row in range(4):
            term = cyclic_mul(term, matrix[row][permutation[row]])
        answer = cyclic_add(answer, cyclic_scale(term, permutation_sign(permutation)))
    return answer


def determinant_mod(matrix: tuple[tuple[int, ...], ...], prime: int) -> int:
    rows = [list(value % prime for value in row) for row in matrix]
    answer = 1
    for column in range(len(rows)):
        pivot = next((row for row in range(column, len(rows)) if rows[row][column]), None)
        if pivot is None:
            return 0
        if pivot != column:
            rows[column], rows[pivot] = rows[pivot], rows[column]
            answer = -answer
        value = rows[column][column]
        answer = answer * value % prime
        inverse = pow(value, -1, prime)
        for row in range(column + 1, len(rows)):
            factor = rows[row][column] * inverse % prime
            if factor:
                rows[row] = [
                    (entry - factor * pivot_entry) % prime
                    for entry, pivot_entry in zip(rows[row], rows[column])
                ]
    return answer % prime


def main() -> None:
    candidate_hashes = (lf_sha256(CANDIDATE_PATH), lf_sha256(CANDIDATE_OUTPUT))
    require(
        candidate_hashes == (CANDIDATE_SHA256, CANDIDATE_OUTPUT_SHA256),
        ("candidate drift", candidate_hashes),
    )

    e_intervals = M.build_set(M.PAT_E, M.ZELL)
    q_intervals = M.build_set(M.PAT_QB, M.ZELL)
    f_pieces = build_f_pieces(e_intervals, q_intervals)
    e_pieces = [(left, right, 1) for left, right in e_intervals]
    require((len(f_pieces), len(e_pieces)) == (120234, 57072), "input piece census")

    f_groups = split_by_atoms(f_pieces)
    e_groups = split_by_atoms(e_pieces)
    require(
        (sum(bool(group) for group in f_groups), sum(bool(group) for group in e_groups))
        == (21, 20),
        "nonempty atom census",
    )

    # The atom indicator is multiplied into f and E here, before this fold.
    u_profiles = tuple(M.weighted_fold(list(group), M.DCOLL, GRID) for group in f_groups)
    v_profiles = tuple(M.weighted_fold(list(group), M.DCOLL, GRID) for group in e_groups)
    u_whole = M.weighted_fold(f_pieces, M.DCOLL, GRID)
    v_whole = M.weighted_fold(e_pieces, M.DCOLL, GRID)
    verify_profile_sum(u_profiles, u_whole, "arrival pre-fold atom partition")
    verify_profile_sum(v_profiles, v_whole, "source pre-fold atom partition")
    require(
        sum(profile_mass(*profile) for profile in u_profiles) == profile_mass(*u_whole)
        and sum(profile_mass(*profile) for profile in v_profiles) == profile_mass(*v_whole),
        "folded partition mass",
    )

    u_windows = tuple(
        tuple(M.extract_window(starts, values, root, P, GRID) for root in range(P))
        for starts, values in u_profiles
    )
    v_windows = tuple(
        tuple(M.extract_window(starts, values, root, P, GRID) for root in range(P))
        for starts, values in v_profiles
    )
    u_whole_windows = tuple(
        M.extract_window(u_whole[0], u_whole[1], root, P, GRID) for root in range(P)
    )
    v_whole_windows = tuple(
        M.extract_window(v_whole[0], v_whole[1], root, P, GRID) for root in range(P)
    )

    # bucket[C,D,d,s] and common_offset[C,D,d,c].
    bucket = [[[[0 for _s in range(P)] for _d in range(P)] for _right in CHAMBERS]
              for _left in CHAMBERS]
    common_offset = [[[[0 for _c in range(P)] for _d in range(P)] for _right in CHAMBERS]
                     for _left in CHAMBERS]
    pair_common_offset = [[[0 for _c in range(P)] for _right in ATOMS]
                          for _left in ATOMS]
    root_marginal = [0] * P
    atom_pair_active = [[False] * len(ATOMS) for _ in ATOMS]
    diagonal_atom_mass = 0
    for left_index, (left_sheet, left_chamber) in enumerate(ATOMS):
        ci = CHAMBER_INDEX[left_chamber]
        for right_index, (right_sheet, right_chamber) in enumerate(ATOMS):
            cj = CHAMBER_INDEX[right_chamber]
            drift = (right_sheet - left_sheet) % P
            for current_root in range(P):
                for source_root in range(P):
                    mass = product_integral(
                        u_windows[left_index][current_root],
                        v_windows[right_index][source_root],
                    )
                    if not mass:
                        continue
                    atom_pair_active[left_index][right_index] = True
                    root_drift = (current_root - source_root) % P
                    root_marginal[root_drift] += mass
                    bucket[ci][cj][drift][root_drift] += mass
                    if root_drift == 0:
                        diagonal_atom_mass += mass
                    if root_drift == (-drift) % P:
                        left_offset = (left_sheet - current_root) % P
                        right_offset = (right_sheet - source_root) % P
                        require(left_offset == right_offset, "gauge offset")
                        common_offset[ci][cj][drift][left_offset] += mass
                        pair_common_offset[left_index][right_index][left_offset] += mass

    atom_pair_support = sum(value for row in atom_pair_active for value in row)
    require(atom_pair_support == 420, atom_pair_support)
    require(diagonal_atom_mass == root_marginal[0] == 0, "s=0 atom law")
    require(all(value > 0 for value in root_marginal[1:]), "root colour support")

    independent_root_marginal = tuple(
        sum(
            product_integral(
                u_whole_windows[current_root],
                v_whole_windows[(current_root - root_drift) % P],
            )
            for current_root in range(P)
        )
        for root_drift in range(P)
    )
    require(tuple(root_marginal) == independent_root_marginal, "whole/atom root marginal")
    require(
        all(
            sum(
                bucket[left][right][drift][root_drift]
                for left in range(3)
                for right in range(3)
                for drift in range(P)
            ) == root_marginal[root_drift]
            for root_drift in range(P)
        ),
        "bucket/root marginal",
    )

    denominator = M.RPKT * M.DCOLL * M.DCOLL * GRID
    total_mass = sum(root_marginal)
    require(Fraction(total_mass, denominator) == 169 * M.I5, "169 I5 anchor")
    require(all(denominator % prime for prime in SPLIT_PRIMES), "split-prime denominator")

    # Exhaust the algebraic iff independently of support.
    gauge_truth_table = tuple(
        (
            left_sheet,
            right_sheet,
            current_root,
            source_root,
            (left_sheet - current_root) % P == (right_sheet - source_root) % P,
            (current_root - source_root) % P == (left_sheet - right_sheet) % P,
        )
        for left_sheet in range(P)
        for right_sheet in range(P)
        for current_root in range(P)
        for source_root in range(P)
    )
    require(all(left == right for *_labels, left, right in gauge_truth_table), "gauge iff")

    gauge_rows = []
    opposite_rows = []
    for left_chamber in CHAMBERS:
        ci = CHAMBER_INDEX[left_chamber]
        for right_chamber in CHAMBERS:
            cj = CHAMBER_INDEX[right_chamber]
            gauge = tuple(bucket[ci][cj][drift][(-drift) % P] for drift in range(P))
            opposite = tuple(bucket[ci][cj][drift][drift] for drift in range(P))
            require(
                all(sum(common_offset[ci][cj][drift]) == gauge[drift] for drift in range(P)),
                (left_chamber, right_chamber, "offset marginal"),
            )
            gauge_rows.append((left_chamber, right_chamber, gauge))
            opposite_rows.append((left_chamber, right_chamber, opposite))

    full_support = sum(value != 0 for _left, _right, row in gauge_rows for value in row)
    opposite_support = sum(value != 0 for _left, _right, row in opposite_rows for value in row)
    require((full_support, opposite_support) == (72, 72), (full_support, opposite_support))
    active_rows = tuple(
        next(row for left0, right0, row in gauge_rows if (left0, right0) == (left, right))
        for left, right in ACTIVE_PAIRS
    )
    active_support = sum(value != 0 for row in active_rows for value in row)
    require(active_support == 48 and all(row[0] == 0 for row in active_rows), "48/52 support")
    require(
        all(
            (value != 0) == (drift != 0)
            for row in active_rows
            for drift, value in enumerate(row)
        ),
        "active support mask",
    )

    k4_rank_q = rank_q(active_rows)
    k4_ranks = tuple(rank_mod(active_rows, prime) for prime in SPLIT_PRIMES)
    require(k4_rank_q == 4 and all(rank == 4 for rank in k4_ranks), (k4_rank_q, k4_ranks))
    walsh_rows = tuple(
        tuple(
            sum(sign * active_rows[row][drift] for row, sign in enumerate(signs))
            for drift in range(P)
        )
        for signs in WALSH_SIGNS
    )
    raw_walsh_support = tuple(sum(value != 0 for value in row) for row in walsh_rows)
    require(raw_walsh_support == (12, 12, 12, 12), raw_walsh_support)

    corner_zero_modes = tuple(sum(row) for row in active_rows)
    zero_modes = tuple(sum(row) for row in walsh_rows)
    require(corner_zero_modes[0] == corner_zero_modes[3], corner_zero_modes)
    require(corner_zero_modes[1] == corner_zero_modes[2], corner_zero_modes)
    require(tuple(value != 0 for value in zero_modes) == (True, False, False, True), zero_modes)

    # Retain c before Walsh contraction.  The exact two-character transform
    # is a cyclic polynomial in one primitive root.
    active_offsets = tuple(
        tuple(
            tuple(common_offset[CHAMBER_INDEX[left]][CHAMBER_INDEX[right]][drift][offset]
                  for offset in range(P))
            for drift in range(P)
        )
        for left, right in ACTIVE_PAIRS
    )
    walsh_offsets = tuple(
        tuple(
            tuple(
                sum(sign * active_offsets[row][drift][offset]
                    for row, sign in enumerate(signs))
                for offset in range(P)
            )
            for drift in range(P)
        )
        for signs in WALSH_SIGNS
    )
    require(
        tuple(tuple(sum(offsets) for offsets in row) for row in walsh_offsets) == walsh_rows,
        "offset marginal/Walsh",
    )

    coefficient_polynomials = tuple(
        tuple(
            tuple(
                tuple(
                    sum(
                        walsh_offsets[walsh][drift][offset]
                        for drift in range(P)
                        for offset in range(P)
                        if (-owner_frequency * offset - drift_frequency * drift) % P
                        == exponent
                    )
                    for exponent in range(P)
                )
                for drift_frequency in range(P)
            )
            for walsh in range(4)
        )
        for owner_frequency in range(P)
    )
    exact_support = tuple(
        tuple(
            sum(not cyclotomic_zero(coefficient_polynomials[k][walsh][frequency])
                for frequency in range(P))
            for walsh in range(4)
        )
        for k in range(P)
    )
    require(
        exact_support[0] == (13, 12, 12, 13)
        and all(counts == (13, 13, 13, 13) for counts in exact_support[1:]),
        exact_support,
    )

    prime_roots = tuple((prime, zeta13_mod(prime)) for prime in SPLIT_PRIMES)
    certificate_masks = []
    prime_support_counts = []
    for owner_frequency in range(P):
        owner_masks = []
        owner_prime_counts = []
        for walsh in range(4):
            row_masks = []
            per_prime = [0] * len(prime_roots)
            for drift_frequency in range(P):
                poly = coefficient_polynomials[owner_frequency][walsh][drift_frequency]
                mask = 0
                for prime_index, (prime, root) in enumerate(prime_roots):
                    value = eval_poly(poly, root, prime)
                    if value:
                        mask |= 1 << prime_index
                        per_prime[prime_index] += 1
                row_masks.append(mask)
                require((mask != 0) == (not cyclotomic_zero(poly)),
                        (owner_frequency, walsh, drift_frequency, mask, poly))
            owner_masks.append(tuple(row_masks))
            owner_prime_counts.append(tuple(per_prime))
        certificate_masks.append(tuple(owner_masks))
        prime_support_counts.append(tuple(owner_prime_counts))

    require(
        all(
            certificate_masks[k][walsh][frequency] != 0
            for k in range(1, P)
            for walsh in range(4)
            for frequency in range(P)
        ),
        "five-prime owner support coverage",
    )
    require(
        certificate_masks[0][1][0] == certificate_masks[0][2][0] == 0,
        "k0 exact zero reductions",
    )

    # Exact K4 rank at every owner frequency.  Split fields choose a minor;
    # its determinant is then recomputed as a cyclic integer polynomial.
    rank_certificates = []
    for owner_frequency in range(P):
        entry_polys = tuple(
            tuple(
                tuple(
                    sum(
                        active_offsets[row][drift][offset]
                        for offset in range(P)
                        if (-owner_frequency * offset) % P == exponent
                    )
                    for exponent in range(P)
                )
                for drift in range(P)
            )
            for row in range(4)
        )
        certificate = None
        for prime, root in prime_roots:
            evaluated = tuple(
                tuple(eval_poly(entry_polys[row][drift], root, prime) for drift in range(P))
                for row in range(4)
            )
            pivots = pivot_columns(evaluated, prime)
            if len(pivots) < 4:
                continue
            columns = pivots[:4]
            poly_matrix = tuple(
                tuple(entry_polys[row][column] for column in columns)
                for row in range(4)
            )
            determinant = determinant_poly(poly_matrix)
            evaluated_minor = tuple(
                tuple(evaluated[row][column] for column in columns)
                for row in range(4)
            )
            require(
                eval_poly(determinant, root, prime) == determinant_mod(evaluated_minor, prime),
                (owner_frequency, prime, columns, "determinant evaluation"),
            )
            require(not cyclotomic_zero(determinant), (owner_frequency, prime, columns))
            certificate = (prime, columns, digest_json(determinant))
            break
        require(certificate is not None, ("owner K4 rank", owner_frequency))
        rank_certificates.append((owner_frequency, certificate))

    defect_mass = tuple(
        sum(
            bucket[left][right][drift][root_drift]
            for left in range(3)
            for right in range(3)
            for drift in range(P)
            for root_drift in range(P)
            if (drift + root_drift) % P == defect
        )
        for defect in range(P)
    )
    require(all(defect_mass), "gauge defect support")
    pair_gauge_digest = digest_json(pair_common_offset)
    require(
        pair_gauge_digest == EXPECTED_PAIR_GAUGE_SHA256,
        (pair_gauge_digest, EXPECTED_PAIR_GAUGE_SHA256),
    )

    semantic_record = (
        STAGE_SHA256,
        candidate_hashes,
        (GRID, ATOM_UNIT, tuple(ATOMS), tuple(ATOM_INTERVALS)),
        (len(f_pieces), len(e_pieces)),
        tuple(tuple(group) for group in f_groups),
        tuple(tuple(group) for group in e_groups),
        tuple((tuple(profile[0]), tuple(profile[1])) for profile in u_profiles),
        tuple((tuple(profile[0]), tuple(profile[1])) for profile in v_profiles),
        atom_pair_support,
        tuple(root_marginal),
        tuple(
            tuple(tuple(tuple(row) for row in drifts) for drifts in right_rows)
            for right_rows in bucket
        ),
        tuple(
            tuple(tuple(tuple(row) for row in drifts) for drifts in right_rows)
            for right_rows in common_offset
        ),
        pair_gauge_digest,
        (total_mass, denominator),
        tuple(gauge_rows),
        tuple(opposite_rows),
        active_rows,
        (full_support, active_support, opposite_support),
        (k4_rank_q, k4_ranks),
        walsh_rows,
        raw_walsh_support,
        corner_zero_modes,
        zero_modes,
        exact_support,
        prime_roots,
        tuple(certificate_masks),
        tuple(prime_support_counts),
        tuple(rank_certificates),
        defect_mass,
        "support/spectral capacity only; endpoint weights and current absent",
    )
    semantic_hash = sha256(repr(semantic_record).encode("utf-8")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(
            semantic_hash == EXPECTED_SEMANTIC_SHA256,
            (semantic_hash, EXPECTED_SEMANTIC_SHA256),
        )

    mask_histogram = tuple(
        sorted(
            (mask, sum(
                certificate_masks[k][walsh][frequency] == mask
                for k in range(P)
                for walsh in range(4)
                for frequency in range(P)
            ))
            for mask in range(1 << len(SPLIT_PRIMES))
            if any(
                certificate_masks[k][walsh][frequency] == mask
                for k in range(P)
                for walsh in range(4)
                for frequency in range(P)
            )
        )
    )
    rank_prime_histogram = tuple(
        (prime, sum(certificate[0] == prime for _k, certificate in rank_certificates))
        for prime in SPLIT_PRIMES
    )

    print("R5 COMMON-ANCESTRY GUARD-ATOM/ROOT-DRIFT -- INDEPENDENT AUDIT")
    print("status=ACCEPT_SCOPED_FINITE_EXACT; support/spectral capacity only; LRC(14) OPEN")
    print(f"dependencies=(thm2594={STAGE_SHA256},candidate={candidate_hashes[0]},candidate_output={candidate_hashes[1]})")
    print(f"pre_fold_partition=(atoms=39,grid={GRID},input_pieces={(len(f_pieces),len(e_pieces))},nonempty_groups={(sum(bool(group) for group in f_groups),sum(bool(group) for group in e_groups))},pointwise_recombination=PASS)")
    print(f"linked_atom_root_tensor=(atom_pair_support={atom_pair_support}/1521,root_support={sum(value != 0 for value in root_marginal)}/13,s0={root_marginal[0]})")
    print(f"service_anchor=(value={Fraction(total_mass, denominator)},equals_169I5={Fraction(total_mass, denominator) == 169*M.I5})")
    print(f"gauge_iff=(truth_table={len(gauge_truth_table)}/{P**4},equation=a-u=b-q_iff_s=-d,offset_marginals=PASS)")
    print(f"gauge_support=(full={full_support}/117,active={active_support}/52,opposite_orientation={opposite_support}/117,active_mask=F13_star)")
    print(f"K4_rank=(rational={k4_rank_q},split_primes={tuple(zip(SPLIT_PRIMES,k4_ranks))})")
    print(f"k0_Walsh=(raw_support={raw_walsh_support},zero_mode_nonzero={tuple(int(value != 0) for value in zero_modes)},exact_drift_spectral_support={exact_support[0]})")
    print(f"owner_nonzero_modes_exact_support={tuple(sorted(set(exact_support[1:])))};all_12=(13,13,13,13)=PASS")
    print(f"split_prime_coordinate_certificate_masks={mask_histogram}")
    print(f"owner_K4_exact_rank_certificates=(count={len(rank_certificates)},prime_histogram={rank_prime_histogram},digest={digest_json(tuple(rank_certificates))})")
    print(f"pair_level_common_gauge_offset_sha256={pair_gauge_digest}")
    print(f"gauge_defect_mass_support={sum(value != 0 for value in defect_mass)}/13")
    print(f"semantic_sha256={semantic_hash}")
    print("typing=the same rational 39-atom partition is evaluated separately on linked X and Y before P_13^5; no one-point identification or endpoint diagonal")
    print("scope=no equality with frozen U_full endpoint weights, no q_H-q5 bridge value, physical current, grouped coefficient, row exclusion, or LRC(14)")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
