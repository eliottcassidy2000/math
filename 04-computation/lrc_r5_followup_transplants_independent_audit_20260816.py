#!/usr/bin/env python3
"""Clean-room audit of the two THM-2594 -> U_full follow-up probes.

The submitted fixed-absolute-root and root-difference probes are hash-pinned
but never imported.  This verifier reruns THM-2594's proved common-base table,
reconstructs the THM-3514/3518 H-q5 Walsh target directly from the independently
audited unguarded endpoint atoms, and uses separately written linear algebra.

Scope: finite exact arithmetic in the certified split field.  No common U_full
ancestry relation, physical current, absolute H1 class, bispectrum, row
exclusion, or LRC(14) consequence is asserted.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from contextlib import redirect_stdout
from fractions import Fraction
from itertools import product
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
N = 13

PRIMARY_PATH = ROOT / "04-computation/lrc14_stage2_theta_contraction_opus_20260728.py"
ATOM_PATH = ROOT / (
    "04-computation/"
    "lrc_ufull_owner_boundary_k4xf13_endpoint_factorization_"
    "independent_audit_20260816.py"
)
CANDIDATE_A_PATH = ROOT / "04-computation/lrc_r5_absolute_root_fourth_channel_probe_20260816.py"
CANDIDATE_A_OUTPUT = ROOT / "05-knowledge/results/lrc_r5_absolute_root_fourth_channel_probe_20260816.out"
CANDIDATE_B_PATH = ROOT / "04-computation/lrc_r5_root_difference_ufull_drift_alignment_probe_20260816.py"
CANDIDATE_B_OUTPUT = ROOT / "05-knowledge/results/lrc_r5_root_difference_ufull_drift_alignment_probe_20260816.out"

PRIMARY_SHA256 = "09c43af0a0a56c7a0833bbfd13ed6a96bc5a7a3718aa1bc6b77a144bde101a06"
ATOM_SHA256 = "f89be10c65bb77270199f9399b155d5a2c82c0da121b3e8589fe3c1f7e9824fc"
CANDIDATE_A_SHA256 = "113793a41d39e8ac1ad3d745b9c9238cba92721805a379ee4378a000bacca33c"
CANDIDATE_A_OUTPUT_SHA256 = "f367cc991df9edae8b27a312d16f2189f62899e8a961f22cbf24a339ba9b10b0"
CANDIDATE_B_SHA256 = "0448a111eb9f37fe8ee79c1c94d4be82e42911c15b72800f61d3316fa229fc86"
CANDIDATE_B_OUTPUT_SHA256 = "5d634711f53f07280b53f4d494800ca3e1dd396f93789c6a40706c1b3e2a8d69"
EXPECTED_A_SLAVED_DIGEST = "666200f5e79e48c02e7beaeec90ebfedccb708bd1e383d1a01997d0abf5073f7"
EXPECTED_A_ABSOLUTE_DIGEST = "9353c6e645b570c659a5d42491a1dd7c79825d4bd338a476aabbd0afe10e69f2"
EXPECTED_B_SHIFT_DIGEST = "ad63d008da3dfa8525db7ef4e724cf6e3a971b3f6efa3d6da9791fbd8003ed91"
EXPECTED_SEMANTIC_SHA256 = "69e4dd6d71290663a560a5809c9fc509a0649891fd28b403a69d92b822046fc5"

CHAMBERS = ("left", "middle", "right")
ATOMS = tuple((sheet, chamber) for sheet in range(N) for chamber in CHAMBERS)
BUCKETS = tuple(
    (left, right, drift)
    for left in CHAMBERS
    for right in CHAMBERS
    for drift in range(N)
)
BUCKET_INDEX = {bucket: index for index, bucket in enumerate(BUCKETS)}
CORNERS = (
    ("left", "left"),
    ("left", "right"),
    ("right", "left"),
    ("right", "right"),
)
Q_H = (1, 0, 1)
Q_Q5 = (1, 0, 0)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def json_digest(value: object) -> str:
    body = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    return hashlib.sha256(body).hexdigest()


def load_module(path: Path, name: str, expected_hash: str):
    actual = lf_sha256(path)
    require(actual == expected_hash, (name, actual, expected_hash))
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, (name, "loader"))
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


ATOM = load_module(ATOM_PATH, "thm3514_atoms_for_r5_followup_audit", ATOM_SHA256)


def atom_worker(alpha: int) -> tuple[object, ...]:
    return ATOM.worker(alpha)


def source_worker() -> tuple[tuple[tuple[tuple[tuple[int, ...], ...], ...], ...], int]:
    primary = load_module(PRIMARY_PATH, "thm2594_primary_for_r5_followup_audit", PRIMARY_SHA256)
    sink = io.StringIO()
    with redirect_stdout(sink):
        state = primary.main()
    joint = state[12]
    denominator = int(state[17])
    frozen = tuple(
        tuple(
            tuple(tuple(int(value) for value in theta_row) for theta_row in ell_bank)
            for ell_bank in q_bank
        )
        for q_bank in joint
    )
    require(len(frozen) == N and all(len(q_bank) == N for q_bank in frozen), "joint shape")
    require(
        all(len(ell_bank) == 7 and all(len(theta_row) == 3 for theta_row in ell_bank)
            for q_bank in frozen for ell_bank in q_bank),
        "joint inner shape",
    )
    return frozen, denominator


def rref(
    matrix: tuple[tuple[int, ...], ...], prime: int
) -> tuple[tuple[tuple[int, ...], ...], tuple[int, ...]]:
    rows = [list(value % prime for value in row) for row in matrix]
    if not rows:
        return (), ()
    columns = len(rows[0])
    require(all(len(row) == columns for row in rows), "ragged matrix")
    pivots: list[int] = []
    pivot_row = 0
    for column in range(columns):
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
                    (value - factor * pivot_value) % prime
                    for value, pivot_value in zip(rows[row], rows[pivot_row])
                ]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == len(rows):
            break
    return tuple(tuple(row) for row in rows), tuple(pivots)


def rank_mod(matrix: tuple[tuple[int, ...], ...], prime: int) -> int:
    return len(rref(matrix, prime)[1])


def nullspace(
    matrix: tuple[tuple[int, ...], ...], prime: int
) -> tuple[tuple[int, ...], ...]:
    require(bool(matrix), "nullspace needs equations")
    reduced, pivots = rref(matrix, prime)
    columns = len(matrix[0])
    free = tuple(column for column in range(columns) if column not in pivots)
    basis = []
    for free_column in free:
        vector = [0] * columns
        vector[free_column] = 1
        for row, pivot in enumerate(pivots):
            vector[pivot] = (-reduced[row][free_column]) % prime
        basis.append(tuple(vector))
    return tuple(basis)


def rank_q(matrix: tuple[tuple[int, ...], ...]) -> int:
    rows = [list(Fraction(value) for value in row) for row in matrix]
    if not rows:
        return 0
    columns = len(rows[0])
    pivot_row = 0
    for column in range(columns):
        pivot = next((row for row in range(pivot_row, len(rows)) if rows[row][column]), None)
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        pivot_value = rows[pivot_row][column]
        rows[pivot_row] = [value / pivot_value for value in rows[pivot_row]]
        for row in range(pivot_row + 1, len(rows)):
            factor = rows[row][column]
            if factor:
                rows[row] = [
                    value - factor * pivot_value
                    for value, pivot_value in zip(rows[row], rows[pivot_row])
                ]
        pivot_row += 1
        if pivot_row == len(rows):
            break
    return pivot_row


def dft(values: tuple[int, ...], root: int, prime: int) -> tuple[int, ...]:
    size = len(values)
    return tuple(
        sum(value * pow(root, (-frequency * index) % size, prime)
            for index, value in enumerate(values)) % prime
        for frequency in range(size)
    )


def septimal_profiles(
    table: tuple[tuple[int, ...], ...], zeta7: int, prime: int
) -> tuple[tuple[int, ...], ...]:
    require(len(table) == 7 and all(len(row) == N for row in table), "7x13 table")
    return tuple(
        tuple(
            sum(table[ell][coordinate] * pow(zeta7, (-mode * ell) % 7, prime)
                for ell in range(7)) % prime
            for coordinate in range(N)
        )
        for mode in range(7)
    )


def spectra(
    rows: tuple[tuple[int, ...], ...], zeta13: int, prime: int
) -> tuple[tuple[int, ...], ...]:
    return tuple(dft(row, zeta13, prime) for row in rows)


def normalize_table(
    table: tuple[tuple[int, ...], ...], denominator: int, prime: int
) -> tuple[tuple[int, ...], ...]:
    inverse = pow(denominator % prime, -1, prime)
    return tuple(tuple(value % prime * inverse % prime for value in row) for row in table)


def materialize_target(chunks: tuple[tuple[object, ...], ...]) -> tuple[object, ...]:
    """Rebuild q_H and q_q5 by direct tau slices and an independent kernel route."""
    require(tuple(chunk[0] for chunk in chunks) == tuple(range(N)), "atom worker order")
    require(tuple(ATOM.ATOMS) == ATOMS and tuple(ATOM.BUCKETS) == BUCKETS, "atom order")
    tables = tuple(row for chunk in chunks for row in chunk[1])
    require(len(tables) == N * N, len(tables))
    (_word, _t_den, nn, prime, root, zeta13, *_rest) = ATOM.context()

    direct_raw = {Q_H: [0] * len(BUCKETS), Q_Q5: [0] * len(BUCKETS)}
    kernel_raw = {Q_H: [0] * len(BUCKETS), Q_Q5: [0] * len(BUCKETS)}
    wrong_phase_raw = [0] * len(BUCKETS)
    kernels = {
        (left, right, drift, frequency): sum(
            int(ATOM.safe(left, sheet))
            * int(ATOM.safe(right, sheet + drift))
            * pow(zeta13, (-frequency * sheet) % N, prime)
            for sheet in range(N)
        ) % prime
        for left in CHAMBERS
        for right in CHAMBERS
        for drift in range(N)
        for frequency in (0, 1)
    }

    table_index = 0
    for alpha in range(N):
        for beta in range(N):
            stored_beta, ax_values, by_values = tables[table_index]
            table_index += 1
            require(stored_beta == beta, (alpha, beta, stored_beta))
            left_atoms = tuple((index, value) for index, value in enumerate(ax_values) if value)
            right_atoms = tuple((index, value) for index, value in enumerate(by_values) if value)

            for tau in range(N):
                slice_values = [0] * len(BUCKETS)
                for left_index, left_value in left_atoms:
                    left_sheet, left_chamber = ATOMS[left_index]
                    if not ATOM.safe(left_chamber, left_sheet + tau):
                        continue
                    for right_index, right_value in right_atoms:
                        right_sheet, right_chamber = ATOMS[right_index]
                        if not ATOM.safe(right_chamber, right_sheet + tau):
                            continue
                        bucket = BUCKET_INDEX[
                            (left_chamber, right_chamber, (right_sheet - left_sheet) % N)
                        ]
                        slice_values[bucket] = (
                            slice_values[bucket] + left_value * right_value
                        ) % prime
                for q in (Q_H, Q_Q5):
                    qa, qb, qt = q
                    phase = pow(
                        zeta13,
                        (beta - alpha * qa - beta * qb - tau * qt) % N,
                        prime,
                    )
                    for bucket, value in enumerate(slice_values):
                        direct_raw[q][bucket] = (
                            direct_raw[q][bucket] + phase * value
                        ) % prime

            for left_index, left_value in left_atoms:
                left_sheet, left_chamber = ATOMS[left_index]
                for right_index, right_value in right_atoms:
                    right_sheet, right_chamber = ATOMS[right_index]
                    drift = (right_sheet - left_sheet) % N
                    bucket = BUCKET_INDEX[(left_chamber, right_chamber, drift)]
                    atom_product = left_value * right_value % prime
                    for q in (Q_H, Q_Q5):
                        qa, qb, qt = q
                        base_phase = pow(
                            zeta13, (beta - alpha * qa - beta * qb) % N, prime
                        )
                        translated_phase = pow(zeta13, (qt * left_sheet) % N, prime)
                        term = (
                            base_phase
                            * translated_phase
                            * kernels[(left_chamber, right_chamber, drift, qt)]
                            * atom_product
                        ) % prime
                        kernel_raw[q][bucket] = (kernel_raw[q][bucket] + term) % prime
                        if q == Q_H:
                            wrong_term = (
                                base_phase
                                * pow(zeta13, (-qt * left_sheet) % N, prime)
                                * kernels[(left_chamber, right_chamber, drift, qt)]
                                * atom_product
                            ) % prime
                            wrong_phase_raw[bucket] = (
                                wrong_phase_raw[bucket] + wrong_term
                            ) % prime

    normalizer = pow(N**3, -1, prime)
    direct = {
        q: tuple(value * normalizer % prime for value in direct_raw[q])
        for q in (Q_H, Q_Q5)
    }
    kernel = {
        q: tuple(value * normalizer % prime for value in kernel_raw[q])
        for q in (Q_H, Q_Q5)
    }
    wrong_phase = tuple(value * normalizer % prime for value in wrong_phase_raw)
    require(direct == kernel, "direct/kernel target disagreement")
    wrong_phase_mismatches = sum(
        left != right for left, right in zip(direct[Q_H], wrong_phase)
    )
    require(wrong_phase_mismatches > 0, "opposite target sheet phase accidentally agrees")

    bridge = tuple((left - right) % prime for left, right in zip(direct[Q_H], direct[Q_Q5]))
    corners = tuple(
        tuple(bridge[BUCKET_INDEX[(left, right, drift)]] for drift in range(N))
        for left, right in CORNERS
    )
    walsh = (
        tuple(sum(corners[row][drift] for row in range(4)) % prime for drift in range(N)),
        tuple((corners[0][d] + corners[1][d] - corners[2][d] - corners[3][d]) % prime for d in range(N)),
        tuple((corners[0][d] - corners[1][d] + corners[2][d] - corners[3][d]) % prime for d in range(N)),
        tuple((corners[0][d] - corners[1][d] - corners[2][d] + corners[3][d]) % prime for d in range(N)),
    )
    target_spectra = spectra(walsh, zeta13, prime)
    require(rank_mod(walsh, prime) == 4, "target rank")
    require(all(all(row) for row in target_spectra), "target spectral support")
    require(sum(walsh[0]) % prime == 389266878372286537904, "bridge regression")
    return nn, prime, root, zeta13, direct, kernels, wrong_phase_mismatches, walsh, target_spectra


def wedge_equations(
    source: tuple[tuple[int, ...], ...],
    target: tuple[tuple[int, ...], ...],
    prime: int,
) -> tuple[tuple[int, ...], ...]:
    channels = len(source)
    equations = []
    for frequency in range(N):
        x = tuple(row[frequency] for row in source)
        y = tuple(row[frequency] for row in target)
        require(any(y), ("zero target", frequency))
        for left in range(4):
            for right in range(left + 1, 4):
                row = [0] * (4 * channels)
                for channel in range(channels):
                    row[right * channels + channel] = y[left] * x[channel] % prime
                    row[left * channels + channel] = (
                        -y[right] * x[channel]
                    ) % prime
                equations.append(tuple(row))
    return tuple(equations)


def anchor_equations(
    source: tuple[tuple[int, ...], ...],
    target: tuple[tuple[int, ...], ...],
    prime: int,
) -> tuple[tuple[int, ...], ...]:
    channels = len(source)
    equations = []
    for frequency in range(N):
        require(target[0][frequency], ("zero target anchor", frequency))
        for output in range(1, 4):
            row = [0] * (4 * channels)
            for channel in range(channels):
                row[output * channels + channel] = (
                    target[0][frequency] * source[channel][frequency]
                ) % prime
                row[channel] = (
                    -target[output][frequency] * source[channel][frequency]
                ) % prime
            equations.append(tuple(row))
    return tuple(equations)


def augmented_equations(
    source: tuple[tuple[int, ...], ...],
    target: tuple[tuple[int, ...], ...],
    prime: int,
) -> tuple[tuple[int, ...], ...]:
    channels = len(source)
    equations = []
    for frequency in range(N):
        for output in range(4):
            row = [0] * (4 * channels + N)
            for channel in range(channels):
                row[output * channels + channel] = source[channel][frequency]
            row[4 * channels + frequency] = (-target[output][frequency]) % prime
            equations.append(tuple(row))
    return tuple(equations)


def projective_signature(
    source: tuple[tuple[int, ...], ...],
    target: tuple[tuple[int, ...], ...],
    prime: int,
) -> tuple[int, int, int, int, int]:
    source_rank = rank_mod(source, prime)
    equation_rank = rank_mod(wedge_equations(source, target, prime), prime)
    nullity = 4 * len(source) - equation_rank
    annihilator = 4 * (len(source) - source_rank)
    require(nullity >= annihilator, (nullity, annihilator))
    return source_rank, equation_rank, nullity, annihilator, nullity - annihilator


def full_projective_audit(
    source: tuple[tuple[int, ...], ...],
    target: tuple[tuple[int, ...], ...],
    prime: int,
) -> tuple[object, ...]:
    signature = projective_signature(source, target, prime)
    wedge_rank = signature[1]
    anchor_rank = rank_mod(anchor_equations(source, target, prime), prime)
    require(anchor_rank == wedge_rank, ("anchor/wedge", anchor_rank, wedge_rank))
    augmented = augmented_equations(source, target, prime)
    augmented_rank = rank_mod(augmented, prime)
    basis = nullspace(augmented, prime)
    channels = len(source)
    lambda_projection = tuple(vector[4 * channels :] for vector in basis)
    lambda_rank = rank_mod(lambda_projection, prime) if lambda_projection else 0
    annihilates = all(
        sum(vector[output * channels + channel] * source[channel][frequency]
            for channel in range(channels)) % prime == 0
        for vector in basis
        for output in range(4)
        for frequency in range(N)
    )
    require(annihilates, "augmented basis contains non-annihilator")
    require(lambda_rank == 0, ("lambda projection", lambda_rank))
    require(len(basis) == signature[3], ("augmented nullity", len(basis), signature))
    return signature, anchor_rank, augmented_rank, len(basis), lambda_rank, annihilates


def dilate(
    rows: tuple[tuple[int, ...], ...], multiplier: int
) -> tuple[tuple[int, ...], ...]:
    require(1 <= multiplier < N, multiplier)
    return tuple(
        tuple(row[multiplier * frequency % N] for frequency in range(N))
        for row in rows
    )


def translate(
    rows: tuple[tuple[int, ...], ...], shift: int, zeta13: int, prime: int
) -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple(value * pow(zeta13, (frequency * shift) % N, prime) % prime
              for frequency, value in enumerate(row))
        for row in rows
    )


def gauge_signature(row: tuple[int, ...], zeta13: int, prime: int) -> tuple[int, ...]:
    require(all(row), "gauge signature needs units")
    normalized = tuple(value * pow(row[0], -1, prime) % prime for value in row)
    return min(
        tuple(
            normalized[frequency] * pow(zeta13, (frequency * shift) % N, prime) % prime
            for frequency in range(N)
        )
        for shift in range(N)
    )


def allocation_census(
    source: tuple[tuple[int, ...], ...],
    target: tuple[tuple[int, ...], ...],
    zeta13: int,
    prime: int,
) -> tuple[int, int, int, tuple[int, ...]]:
    exact = 0
    common_shift_gauge = 0
    classes_seen = set()
    count = 0
    for allocation in product(range(7), repeat=4):
        count += 1
        ratios = tuple(
            tuple(
                target[channel][frequency]
                * pow(source[allocation[channel]][frequency], -1, prime)
                % prime
                for frequency in range(N)
            )
            for channel in range(4)
        )
        exact += len(set(ratios)) == 1
        signatures = tuple(gauge_signature(row, zeta13, prime) for row in ratios)
        classes = len(set(signatures))
        classes_seen.add(classes)
        common_shift_gauge += classes == 1
    return count, exact, common_shift_gauge, tuple(sorted(classes_seen))


def primitive_root_small(prime: int) -> int:
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


def source_split_check(
    table: tuple[tuple[int, ...], ...], denominator: int, prime: int
) -> tuple[int, int]:
    require((prime - 1) % 91 == 0 and denominator % prime, (prime, denominator % prime))
    generator = primitive_root_small(prime)
    zeta91 = pow(generator, (prime - 1) // 91, prime)
    zeta7 = pow(zeta91, 13, prime)
    zeta13 = pow(zeta91, 7, prime)
    normalized = normalize_table(table, denominator, prime)
    profiles = septimal_profiles(normalized, zeta7, prime)
    transformed = spectra(profiles, zeta13, prime)
    return sum(value != 0 for row in transformed for value in row), rank_mod(table, prime)


def owner_bank(
    joint: tuple[tuple[tuple[tuple[int, ...], ...], ...], ...],
    exponent_sign: int,
    zeta13: int,
    prime: int,
) -> tuple[tuple[tuple[int, ...], ...], tuple[tuple[int, ...], ...]]:
    k0 = [[0] * N for _ in range(7)]
    k1 = [[0] * N for _ in range(7)]
    for u in range(N):
        phase = pow(zeta13, (exponent_sign * u) % N, prime)
        for shift in range(N):
            q = (u - shift) % N
            for ell in range(7):
                value = sum(joint[u][q][ell]) % prime
                k0[ell][shift] = (k0[ell][shift] + value) % prime
                k1[ell][shift] = (k1[ell][shift] + phase * value) % prime
    return tuple(tuple(row) for row in k0), tuple(tuple(row) for row in k1)


def main() -> None:
    pinned = (
        lf_sha256(CANDIDATE_A_PATH),
        lf_sha256(CANDIDATE_A_OUTPUT),
        lf_sha256(CANDIDATE_B_PATH),
        lf_sha256(CANDIDATE_B_OUTPUT),
    )
    require(
        pinned == (
            CANDIDATE_A_SHA256,
            CANDIDATE_A_OUTPUT_SHA256,
            CANDIDATE_B_SHA256,
            CANDIDATE_B_OUTPUT_SHA256,
        ),
        ("candidate drift", pinned),
    )

    with ProcessPoolExecutor(max_workers=5) as pool:
        source_future = pool.submit(source_worker)
        atom_futures = tuple(pool.submit(atom_worker, alpha) for alpha in range(N))
        joint, denominator = source_future.result()
        chunks = tuple(future.result() for future in atom_futures)

    (
        nn,
        prime,
        root,
        zeta13,
        target_buckets,
        target_kernels,
        wrong_target_phase_mismatches,
        walsh,
        target_spectra,
    ) = materialize_target(chunks)
    require(nn % 7 == 0 and denominator % prime, (nn, denominator % prime))
    zeta7 = pow(root, nn // 7, prime)
    require(pow(zeta7, 7, prime) == 1 and zeta7 != 1, "zeta7")
    require(pow(zeta13, N, prime) == 1 and zeta13 != 1, "zeta13")

    # Package A: relative theta and fixed absolute-root marginals.
    slaved_raw = [[0] * N for _ in range(7)]
    absolute_raw = [[0] * N for _ in range(7)]
    for u in range(N):
        for q in range(N):
            for ell in range(7):
                for theta in range(3):
                    value = joint[u][q][ell][theta]
                    slaved_raw[ell][theta] += value
                    absolute_raw[ell][(theta + 2 * u) % N] += value
    slaved_raw_t = tuple(tuple(row) for row in slaved_raw)
    absolute_raw_t = tuple(tuple(row) for row in absolute_raw)
    slaved_residue_t = tuple(tuple(value % prime for value in row) for row in slaved_raw_t)
    absolute_residue_t = tuple(tuple(value % prime for value in row) for row in absolute_raw_t)
    require(json_digest(slaved_residue_t) == EXPECTED_A_SLAVED_DIGEST, "A slaved digest")
    require(json_digest(absolute_residue_t) == EXPECTED_A_ABSOLUTE_DIGEST, "A absolute digest")
    require(
        (sum(value != 0 for row in slaved_raw_t for value in row),
         sum(value != 0 for row in absolute_raw_t for value in row)) == (12, 18),
        "A physical support",
    )
    a_q_ranks = (
        rank_q(slaved_raw_t),
        rank_q(absolute_raw_t),
        rank_q(slaved_raw_t + absolute_raw_t),
    )
    require(a_q_ranks == (3, 3, 4), a_q_ranks)

    slaved = normalize_table(slaved_raw_t, denominator, prime)
    absolute = normalize_table(absolute_raw_t, denominator, prime)
    slaved_profiles = septimal_profiles(slaved, zeta7, prime)
    absolute_profiles = septimal_profiles(absolute, zeta7, prime)
    slaved_spectra = spectra(slaved_profiles, zeta13, prime)
    absolute_spectra = spectra(absolute_profiles, zeta13, prime)
    slaved_residue_spectra = spectra(
        septimal_profiles(slaved_residue_t, zeta7, prime), zeta13, prime
    )
    absolute_residue_spectra = spectra(
        septimal_profiles(absolute_residue_t, zeta7, prime), zeta13, prime
    )
    a_support = (
        sum(value != 0 for row in slaved_spectra for value in row),
        sum(value != 0 for row in absolute_spectra for value in row),
    )
    require(a_support == (91, 91), a_support)
    require(
        (rank_mod(slaved_profiles, prime), rank_mod(absolute_profiles, prime),
         rank_mod(slaved_profiles + absolute_profiles, prime)) == (3, 3, 4),
        "A split ranks",
    )

    a_core = {
        "slaved7": full_projective_audit(slaved_spectra, target_spectra, prime),
        "absolute7": full_projective_audit(absolute_spectra, target_spectra, prime),
        "union14": full_projective_audit(
            slaved_spectra + absolute_spectra, target_spectra, prime
        ),
    }
    require(
        tuple(value[0] for value in a_core.values())
        == ((3, 12, 16, 16, 0), (3, 12, 16, 16, 0), (4, 16, 40, 40, 0)),
        a_core,
    )
    a_denominator_invariance = (
        projective_signature(
            slaved_residue_spectra + absolute_residue_spectra,
            target_spectra,
            prime,
        ),
        a_core["union14"][0],
    )
    require(a_denominator_invariance[0] == a_denominator_invariance[1], a_denominator_invariance)
    a_dilations = tuple(
        projective_signature(
            dilate(slaved_spectra + absolute_spectra, multiplier),
            target_spectra,
            prime,
        )
        for multiplier in range(1, N)
    )
    require(a_dilations == ((4, 16, 40, 40, 0),) * 12, a_dilations)
    a_translations = tuple(
        projective_signature(
            translate(slaved_spectra + absolute_spectra, shift, zeta13, prime),
            target_spectra,
            prime,
        )
        for shift in range(N)
    )
    require(a_translations == ((4, 16, 40, 40, 0),) * N, a_translations)
    a_relative_affine = tuple(
        (
            multiplier,
            relative_shift,
            projective_signature(
                dilate(slaved_spectra, multiplier)
                + translate(
                    dilate(absolute_spectra, multiplier),
                    relative_shift,
                    zeta13,
                    prime,
                ),
                target_spectra,
                prime,
            ),
        )
        for multiplier in range(1, N)
        for relative_shift in range(N)
    )
    require(len(a_relative_affine) == 156, len(a_relative_affine))
    require(
        all(record[2][4] == 0 for record in a_relative_affine),
        tuple(record for record in a_relative_affine if record[2][4] != 0),
    )

    folded_pairs = ((1, 6), (2, 5), (3, 4))
    named_signatures = []
    for absolute_mode in range(7):
        for choices in product((0, 1), repeat=3):
            selected = tuple(folded_pairs[index][choices[index]] for index in range(3))
            rows = tuple(slaved_spectra[index] for index in selected) + (
                absolute_spectra[absolute_mode],
            )
            named_signatures.append(projective_signature(rows, target_spectra, prime))
    require(len(named_signatures) == 56, len(named_signatures))
    require(tuple(named_signatures) == ((4, 16, 0, 0, 0),) * 56, set(named_signatures))

    binary_signatures = []
    zero_binary_masks = []
    for mask in range(1, 1 << 7):
        sidecar = tuple(
            sum(absolute_spectra[mode][frequency]
                for mode in range(7) if mask & (1 << mode)) % prime
            for frequency in range(N)
        )
        if not any(sidecar):
            zero_binary_masks.append(mask)
        rows = tuple(slaved_spectra[mode] for mode in (1, 2, 3)) + (sidecar,)
        binary_signatures.append(projective_signature(rows, target_spectra, prime))
    require(len(binary_signatures) == 127, len(binary_signatures))
    require(zero_binary_masks == [127], zero_binary_masks)
    require(binary_signatures.count((4, 16, 0, 0, 0)) == 126, "A binary rank4")
    require(binary_signatures.count((3, 12, 4, 4, 0)) == 1, "A binary rank3")

    # Package B: the typed root difference, theta refinement, and owner phases.
    shift_raw = [[0] * N for _ in range(7)]
    theta_raw = [[[0] * N for _ in range(3)] for _ in range(7)]
    for u in range(N):
        for shift in range(N):
            q = (u - shift) % N
            for ell in range(7):
                for theta in range(3):
                    value = joint[u][q][ell][theta]
                    shift_raw[ell][shift] += value
                    theta_raw[ell][theta][shift] += value
    shift_raw_t = tuple(tuple(row) for row in shift_raw)
    shift_residue_t = tuple(tuple(value % prime for value in row) for row in shift_raw_t)
    require(json_digest(shift_residue_t) == EXPECTED_B_SHIFT_DIGEST, "B shift digest")
    require(sum(value != 0 for row in shift_raw_t for value in row) == 72, "B physical support")
    require(all(row[0] == 0 for row in shift_raw_t), "B same-root column")
    require(rank_q(shift_raw_t) == 6, "B rational rank")
    theta_rows_q = tuple(
        tuple(theta_raw[ell][theta][shift] for shift in range(N))
        for theta in range(3)
        for ell in range(7)
    )
    require(rank_q(theta_rows_q) == 8, "B theta rational rank")

    shift = normalize_table(shift_raw_t, denominator, prime)
    source_profiles = septimal_profiles(shift, zeta7, prime)
    source_spectra = spectra(source_profiles, zeta13, prime)
    source_residue_spectra = spectra(
        septimal_profiles(shift_residue_t, zeta7, prime), zeta13, prime
    )
    require(all(all(row) for row in source_spectra), "B source spectral zero")
    require(rank_mod(source_profiles, prime) == 6, "B source split rank")
    theta_rows = tuple(
        tuple(
            sum(
                theta_raw[ell][theta][shift] % prime
                * pow(denominator, -1, prime)
                * pow(zeta7, (-mode * ell) % 7, prime)
                for ell in range(7)
            ) % prime
            for shift in range(N)
        )
        for theta in range(3)
        for mode in range(7)
    )
    theta_spectra = spectra(theta_rows, zeta13, prime)
    require(rank_mod(theta_rows, prime) == 8, "B theta split rank")

    b_core = {
        "root_difference7": full_projective_audit(source_spectra, target_spectra, prime),
        "theta_resolved21": full_projective_audit(theta_spectra, target_spectra, prime),
    }
    require(
        tuple(value[0] for value in b_core.values())
        == ((6, 24, 4, 4, 0), (8, 32, 52, 52, 0)),
        b_core,
    )
    b_denominator_invariance = (
        projective_signature(source_residue_spectra, target_spectra, prime),
        b_core["root_difference7"][0],
    )
    require(b_denominator_invariance[0] == b_denominator_invariance[1], b_denominator_invariance)
    b_dilations = tuple(
        projective_signature(dilate(source_spectra, multiplier), target_spectra, prime)
        for multiplier in range(1, N)
    )
    b_theta_dilations = tuple(
        projective_signature(dilate(theta_spectra, multiplier), target_spectra, prime)
        for multiplier in range(1, N)
    )
    require(b_dilations == ((6, 24, 4, 4, 0),) * 12, b_dilations)
    require(b_theta_dilations == ((8, 32, 52, 52, 0),) * 12, b_theta_dilations)

    owner_audits = {}
    for label, sign in (("minus", -1), ("plus", 1)):
        owner0_raw, owner1_raw = owner_bank(joint, sign, zeta13, prime)
        require(owner0_raw == shift_residue_t, label)
        owner0 = normalize_table(owner0_raw, denominator, prime)
        owner1 = normalize_table(owner1_raw, denominator, prime)
        profile0 = septimal_profiles(owner0, zeta7, prime)
        profile1 = septimal_profiles(owner1, zeta7, prime)
        difference = tuple(
            tuple((profile1[mode][shift] - profile0[mode][shift]) % prime for shift in range(N))
            for mode in range(7)
        )
        spectrum0 = spectra(profile0, zeta13, prime)
        spectrum1 = spectra(profile1, zeta13, prime)
        difference_spectrum = spectra(difference, zeta13, prime)
        union = spectrum0 + spectrum1
        ranks = (
            rank_mod(profile0, prime),
            rank_mod(profile1, prime),
            rank_mod(difference, prime),
            rank_mod(profile0 + profile1, prime),
        )
        systems = (
            full_projective_audit(spectrum0, target_spectra, prime),
            full_projective_audit(spectrum1, target_spectra, prime),
            full_projective_audit(difference_spectrum, target_spectra, prime),
            full_projective_audit(union, target_spectra, prime),
        )
        dilations = tuple(
            projective_signature(dilate(union, multiplier), target_spectra, prime)
            for multiplier in range(1, N)
        )
        owner_audits[label] = (ranks, systems, dilations)
    require(owner_audits["minus"][0] == (6, 6, 6, 9), owner_audits["minus"][0])
    require(owner_audits["plus"][0] == (6, 6, 6, 9), owner_audits["plus"][0])
    expected_owner_signatures = (
        (6, 24, 4, 4, 0),
        (6, 24, 4, 4, 0),
        (6, 24, 4, 4, 0),
        (9, 36, 20, 20, 0),
    )
    for label in owner_audits:
        require(
            tuple(value[0] for value in owner_audits[label][1]) == expected_owner_signatures,
            (label, owner_audits[label][1]),
        )
        require(
            owner_audits[label][2] == ((9, 36, 20, 20, 0),) * 12,
            (label, owner_audits[label][2]),
        )

    owner_all_modes = []
    owner_coupled_affine = []
    for exponent in range(1, N):
        owner0_raw, ownerk_raw = owner_bank(joint, exponent, zeta13, prime)
        owner0 = normalize_table(owner0_raw, denominator, prime)
        ownerk = normalize_table(ownerk_raw, denominator, prime)
        profile0 = septimal_profiles(owner0, zeta7, prime)
        profilek = septimal_profiles(ownerk, zeta7, prime)
        spectrum0 = spectra(profile0, zeta13, prime)
        spectrumk = spectra(profilek, zeta13, prime)
        union = spectrum0 + spectrumk
        ranks = (rank_mod(profile0, prime), rank_mod(profilek, prime), rank_mod(profile0 + profilek, prime))
        signature = projective_signature(union, target_spectra, prime)
        owner_all_modes.append((exponent, ranks, signature))
        for drift_multiplier in range(1, N):
            owner_coupled_affine.append(
                (
                    exponent,
                    drift_multiplier,
                    projective_signature(
                        dilate(union, drift_multiplier), target_spectra, prime
                    ),
                )
            )
    require(len(owner_all_modes) == 12 and len(owner_coupled_affine) == 144, "owner affine census")
    require(
        {record[1] for record in owner_all_modes} == {(6, 6, 9)},
        {record[1] for record in owner_all_modes},
    )
    require(
        {record[2] for record in owner_all_modes}
        == {record[2] for record in owner_coupled_affine}
        == {(9, 36, 20, 20, 0)},
        ({record[2] for record in owner_all_modes}, {record[2] for record in owner_coupled_affine}),
    )

    folded_signatures = []
    folded_selections = []
    for choices in product((0, 1), repeat=3):
        selected = (0,) + tuple(folded_pairs[index][choices[index]] for index in range(3))
        rows = tuple(source_spectra[mode] for mode in selected)
        folded_selections.append(selected)
        folded_signatures.append(projective_signature(rows, target_spectra, prime))
    require(len(set(folded_selections)) == 8, folded_selections)
    require(tuple(folded_signatures) == ((4, 16, 0, 0, 0),) * 8, folded_signatures)
    allocations = allocation_census(source_spectra, target_spectra, zeta13, prime)
    require(allocations == (2401, 0, 0, (4,)), allocations)

    small_primes = (547, 1093, 2003)
    source_split_checks = tuple(
        (
            small_prime,
            source_split_check(slaved_raw_t, denominator, small_prime),
            source_split_check(absolute_raw_t, denominator, small_prime),
            source_split_check(shift_raw_t, denominator, small_prime),
        )
        for small_prime in small_primes
    )
    require(
        all(record[1:] == ((91, 3), (91, 3), (91, 6)) for record in source_split_checks),
        source_split_checks,
    )

    # DFT sign reversals are relative dilations by -1.  The explicit checks
    # below also cover all translations, which are common nonzero frequency
    # scalars and hence absorb into the projective multipliers.
    dft_sign_checks = (
        projective_signature(dilate(source_spectra, N - 1), target_spectra, prime),
        projective_signature(source_spectra, dilate(target_spectra, N - 1), prime),
        projective_signature(
            tuple(source_spectra[(-mode) % 7] for mode in range(7)),
            target_spectra,
            prime,
        ),
    )
    require(dft_sign_checks == ((6, 24, 4, 4, 0),) * 3, dft_sign_checks)

    semantic = (
        (PRIMARY_SHA256, ATOM_SHA256),
        pinned,
        (prime, root, nn, zeta7, zeta13, denominator),
        target_buckets,
        target_kernels,
        wrong_target_phase_mismatches,
        walsh,
        target_spectra,
        slaved_raw_t,
        absolute_raw_t,
        slaved_residue_t,
        absolute_residue_t,
        a_q_ranks,
        a_support,
        a_core,
        a_denominator_invariance,
        a_dilations,
        a_translations,
        a_relative_affine,
        tuple(named_signatures),
        tuple(binary_signatures),
        tuple(zero_binary_masks),
        shift_raw_t,
        shift_residue_t,
        tuple(tuple(tuple(row) for row in ell_bank) for ell_bank in theta_raw),
        source_profiles,
        source_spectra,
        theta_rows,
        theta_spectra,
        b_core,
        b_denominator_invariance,
        b_dilations,
        b_theta_dilations,
        owner_audits,
        tuple(owner_all_modes),
        tuple(owner_coupled_affine),
        tuple(folded_selections),
        tuple(folded_signatures),
        allocations,
        source_split_checks,
        dft_sign_checks,
        "fixed mixing plus one common circulant only; broader typed relations remain open",
    )
    semantic_hash = hashlib.sha256(repr(semantic).encode("utf-8")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(
            semantic_hash == EXPECTED_SEMANTIC_SHA256,
            (semantic_hash, EXPECTED_SEMANTIC_SHA256),
        )

    print("R5 FOLLOW-UP TRANSPLANTS -- INDEPENDENT AUDIT")
    print("status=ACCEPT_SCOPED_FINITE_EXACT; split-field/common-circulant only; LRC(14) OPEN")
    print(f"dependencies=(thm2594={PRIMARY_SHA256},thm3514_atom_audit={ATOM_SHA256})")
    print(f"candidate_hashes=(A={pinned[0]},A_out={pinned[1]},B={pinned[2]},B_out={pinned[3]})")
    print(f"field=(prime={prime},zeta7={zeta7},zeta13={zeta13},denominator_unit={denominator % prime != 0})")
    print(f"target_reconstruction=(direct_equals_kernel=True,wrong_left_sheet_phase_mismatches={wrong_target_phase_mismatches},walsh_rank=4,spectral_support=52/52)")
    print(f"A_support=(physical_slaved=12/91,physical_absolute=18/91,spectral={a_support[0]}/91+{a_support[1]}/91)")
    print(f"A_exact_ranks=(slaved={a_q_ranks[0]},absolute={a_q_ranks[1]},union={a_q_ranks[2]},intersection={a_q_ranks[0]+a_q_ranks[1]-a_q_ranks[2]})")
    print(f"A_core_projective={a_core}")
    print(f"A_affine_gauge=(common_translations={len(a_translations)}/13,dilations={len(a_dilations)}/12,relative_phase_x_dilation={len(a_relative_affine)}/156,signature_classes={tuple(sorted({record[2] for record in a_relative_affine}))},all_zero_excess=True)")
    print(f"A_sidecars=(named_cases={len(named_signatures)},named_rank4={named_signatures.count((4,16,0,0,0))},binary_cases={len(binary_signatures)},binary_rank4={binary_signatures.count((4,16,0,0,0))},zero_mask={tuple(zero_binary_masks)},all_zero_excess=True)")
    print(f"B_support=(physical=72/91,same_root_zero=True,spectral=91/91)")
    print(f"B_exact_ranks=(root_difference=6,theta_resolved=8)")
    print(f"B_core_projective={b_core}")
    print(f"B_affine_gauge=(root_dilations={len(b_dilations)}/12,theta_dilations={len(b_theta_dilations)}/12,owner_minus_dilations={len(owner_audits['minus'][2])}/12,owner_plus_dilations={len(owner_audits['plus'][2])}/12,all_zero_excess=True)")
    print(f"B_owner_phase_minus=(ranks={owner_audits['minus'][0]},systems={owner_audits['minus'][1]})")
    print(f"B_owner_phase_plus=(ranks={owner_audits['plus'][0]},systems={owner_audits['plus'][1]})")
    print(f"B_owner_full_affine=(owner_modes={len(owner_all_modes)}/12,coupled_owner_mode_x_drift_dilation={len(owner_coupled_affine)}/144,rank_classes={tuple(sorted({record[1] for record in owner_all_modes}))},signature_classes={tuple(sorted({record[2] for record in owner_coupled_affine}))})")
    print(f"B_folded=(cases={len(folded_selections)},signatures={tuple(folded_signatures)})")
    print(f"B_allocations=(all_7pow4={allocations[0]},common_exact={allocations[1]},common_amplitude_shift={allocations[2]},kernel_class_counts={allocations[3]})")
    print(f"projective_formulation=(tuple_fields=source_rank,wedge_rank,nullity,annihilator_dimension,excess;anchor_equals_wedge=True;augmented_lambda_projection_zero=True)")
    print(f"source_split_prime_crosschecks={source_split_checks}")
    print("sign_and_gauge_ledger=(target_phase=+a forced by direct tau route; source_owner_+u_and_-u_both_audited; drift_DFT_reversal_is_dilation_12; translations_are_absorbed_frequency_phases)")
    print(f"denominator_ledger=(one common DENC scalar; invertible at certified and three hostile split primes; A_projective_invariance={a_denominator_invariance[0] == a_denominator_invariance[1]};B_projective_invariance={b_denominator_invariance[0] == b_denominator_invariance[1]})")
    print(f"semantic_sha256={semantic_hash}")
    print("scope=fixed frequency-independent channel mixing followed by one common C13-circulant/projective multiplier; arbitrary right operators, nonlinear support relations, and pre-marginal ancestry couplings are not excluded")
    print("nonconsequence=no U_full common ancestry, Boolean coupling, physical current, absolute H1, bispectrum, scalar exclusion, or LRC(14)")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
