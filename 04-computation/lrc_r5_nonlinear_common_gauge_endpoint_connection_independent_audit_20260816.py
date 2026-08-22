#!/usr/bin/env python3
"""Independent hostile audit of the nonlinear R5 common-gauge package.

This program does not import either submitted common-gauge probe.  It rebuilds
the source guard atoms before the P_(13^5) transfer from the hash-pinned
THM-2594 interval engine and rebuilds the target pair function from the
independently audited THM-3514 endpoint-atom engine.  Candidate arrays and
candidate digests enter only in a later, external comparison.

Scope: exact finite common-address pullbacks and split-field nonvanishing.
Nothing here is a one-integral temporal current, an absolute H1 class, a
scalar-row exclusion, or an LRC(14) theorem.
"""

from __future__ import annotations

from bisect import bisect_right
from concurrent.futures import ProcessPoolExecutor
from contextlib import redirect_stdout
from itertools import combinations
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
PRIMARY_PATH = ROOT / "04-computation/lrc14_stage2_theta_contraction_opus_20260728.py"
ATOM_PATH = ROOT / (
    "04-computation/"
    "lrc_ufull_owner_boundary_k4xf13_endpoint_factorization_"
    "independent_audit_20260816.py"
)
PRIMARY_SHA256 = "09c43af0a0a56c7a0833bbfd13ed6a96bc5a7a3718aa1bc6b77a144bde101a06"
ATOM_SHA256 = "f89be10c65bb77270199f9399b155d5a2c82c0da121b3e8589fe3c1f7e9824fc"
EXPECTED_ENDPOINT_FULL = 389266878372286537904
EXPECTED_ENDPOINT_DNEQ = 64768431058832930873
FINAL_COMPARISON_SOURCE_PAIR_SHA256 = "e4d0d4fa674e1f54496e613f7a3e1f057af033fa8992322a5f414ea176e1c3d4"
FINAL_COMPARISON_TARGET_PAIR_SHA256 = "c2d5911b287510335edc6aefa6d3b865c982568f678bb89ee9b82ee211962df1"
EXPECTED_SEMANTIC_SHA256: str | None = "1cf457b221107326b2cd7a5bb56a0cb1626934a2e5d8df57dcb04f04eab56d02"

P = 13
CHAMBERS = (("left", 0, 1), ("middle", 1, 6), ("right", 6, 7))
CHAMBER_NAMES = tuple(name for name, _lo, _hi in CHAMBERS)
ACTIVE = ("left", "right")
CORNERS = (("left", "left"), ("left", "right"), ("right", "left"), ("right", "right"))
ATOMS = tuple((sheet, chamber) for sheet in range(P) for chamber in CHAMBER_NAMES)
ATOM_INDEX = {atom: index for index, atom in enumerate(ATOMS)}
Q_H = (1, 0, 1)
Q_Q5 = (1, 0, 0)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()


def load_module(path: Path, name: str, expected_hash: str):
    actual = lf_sha256(path)
    require(actual == expected_hash, (name, actual, expected_hash))
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, (name, "loader"))
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def split_weighted_pieces(
    pieces: list[tuple[int, int, int]],
    atom_intervals: tuple[tuple[int, str, int, int], ...],
) -> tuple[tuple[tuple[int, int, int], ...], ...]:
    """Intersect disjoint weighted pieces with the contiguous 39-atom partition."""
    groups: list[list[tuple[int, int, int]]] = [[] for _ in atom_intervals]
    atom_index = 0
    mass_in = 0
    mass_out = 0
    for left, right, weight in pieces:
        require(left < right and weight >= 0, (left, right, weight))
        mass_in += weight * (right - left)
        while atom_intervals[atom_index][3] <= left:
            atom_index += 1
        scan = atom_index
        while scan < len(atom_intervals) and atom_intervals[scan][2] < right:
            _sheet, _chamber, atom_left, atom_right = atom_intervals[scan]
            lo = max(left, atom_left)
            hi = min(right, atom_right)
            if lo < hi:
                groups[scan].append((lo, hi, weight))
                mass_out += weight * (hi - lo)
            if atom_right >= right:
                break
            scan += 1
    require(mass_in == mass_out, ("weighted atom partition", mass_in, mass_out))
    return tuple(tuple(group) for group in groups)


def profile_value(starts: list[int], values: list[int], point: int) -> int:
    return values[bisect_right(starts, point) - 1]


def profiles_sum_equal(
    profiles: tuple[tuple[list[int], list[int]], ...],
    target: tuple[list[int], list[int]],
    grid: int,
) -> bool:
    breakpoints = {0, grid}
    breakpoints.update(target[0])
    for starts, _values in profiles:
        breakpoints.update(starts)
    for point in sorted(breakpoints)[:-1]:
        if sum(profile_value(starts, values, point) for starts, values in profiles) != profile_value(
            target[0], target[1], point
        ):
            return False
    return True


def source_worker() -> tuple[object, ...]:
    """Rebuild source atoms before transfer and retain the common gauge c."""
    M = load_module(PRIMARY_PATH, "thm2594_source_for_nonlinear_audit", PRIMARY_SHA256)
    capture = io.StringIO()
    with redirect_stdout(capture):
        state = M.main()
    E, _len_e, QB = state[0], state[1], state[2]
    pair_mass = state[11]
    denominator = int(state[17])

    n_starts, n_values = M.weighted_fold([(a, b, 1) for a, b in E], M.RPKT, M.T_DEN)
    f_pieces: list[tuple[int, int, int]] = []
    for q_left, q_right in QB:
        index = bisect_right(n_starts, q_left) - 1
        while True:
            p_left = n_starts[index]
            p_right = n_starts[index + 1] if index + 1 < len(n_starts) else M.T_DEN
            left = max(q_left, p_left)
            right = min(q_right, p_right)
            if left < right and n_values[index]:
                f_pieces.append((left, right, n_values[index]))
            if p_right >= q_right:
                break
            index += 1

    unit = M.T_DEN // 91
    atom_intervals = tuple(
        (sheet, chamber, (7 * sheet + low) * unit, (7 * sheet + high) * unit)
        for sheet in range(P)
        for chamber, low, high in CHAMBERS
    )
    require(atom_intervals[0][2] == 0 and atom_intervals[-1][3] == M.T_DEN, "atom endpoints")
    require(
        all(left[3] == right[2] for left, right in zip(atom_intervals, atom_intervals[1:])),
        "atom partition gaps",
    )

    f_groups = split_weighted_pieces(f_pieces, atom_intervals)
    e_groups = split_weighted_pieces([(a, b, 1) for a, b in E], atom_intervals)
    u_profiles = tuple(M.weighted_fold(list(group), M.DCOLL, M.T_DEN) for group in f_groups)
    v_profiles = tuple(M.weighted_fold(list(group), M.DCOLL, M.T_DEN) for group in e_groups)
    u_windows = tuple(
        tuple(M.extract_window(starts, values, root, P, M.T_DEN) for root in range(P))
        for starts, values in u_profiles
    )
    v_windows = tuple(
        tuple(M.extract_window(starts, values, root, P, M.T_DEN) for root in range(P))
        for starts, values in v_profiles
    )

    full_u = M.weighted_fold(f_pieces, M.DCOLL, M.T_DEN)
    full_v = M.weighted_fold([(a, b, 1) for a, b in E], M.DCOLL, M.T_DEN)
    restoration = []
    for root in range(P):
        target_u = M.extract_window(full_u[0], full_u[1], root, P, M.T_DEN)
        target_v = M.extract_window(full_v[0], full_v[1], root, P, M.T_DEN)
        restoration.append(
            (
                profiles_sum_equal(tuple(u_windows[index][root] for index in range(len(ATOMS))), target_u, M.T_DEN),
                profiles_sum_equal(tuple(v_windows[index][root] for index in range(len(ATOMS))), target_v, M.T_DEN),
            )
        )
    require(all(left and right for left, right in restoration), restoration)

    # Middle source atoms need not vanish: the 26-atom boundary is a target
    # owner-support restriction, not a theorem about the transferred source.
    middle_indices = tuple(index for index, (_sheet, chamber) in enumerate(ATOMS) if chamber == "middle")
    middle_nonzero_profiles = (
        sum(any(value != 0 for value in u_profiles[index][1]) for index in middle_indices),
        sum(any(value != 0 for value in v_profiles[index][1]) for index in middle_indices),
    )

    chamber_index = {name: index for index, name in enumerate(CHAMBER_NAMES)}
    gauge_offset = [[[[0 for _c in range(P)] for _d in range(P)] for _right in CHAMBERS]
                    for _left in CHAMBERS]
    pair_gauge_offset = [[[0 for _c in range(P)] for _right in ATOMS] for _left in ATOMS]
    pair_level_nonzero = 0
    for left_index, (left_sheet, left_chamber) in enumerate(ATOMS):
        ci = chamber_index[left_chamber]
        for right_index, (right_sheet, right_chamber) in enumerate(ATOMS):
            cj = chamber_index[right_chamber]
            drift = (right_sheet - left_sheet) % P
            pair_nonzero = False
            for current_root in range(P):
                source_root = (current_root + drift) % P
                root_shift = (current_root - source_root) % P
                require((root_shift + drift) % P == 0, (current_root, source_root, drift))
                common_offset = (left_sheet - current_root) % P
                require(
                    common_offset == (right_sheet - source_root) % P,
                    (left_sheet, right_sheet, current_root, source_root),
                )
                left_profile = u_windows[left_index][current_root]
                right_profile = v_windows[right_index][source_root]
                mass = M.product_cum(
                    left_profile[0], left_profile[1], right_profile[0], right_profile[1], M.T_DEN
                )[3]
                gauge_offset[ci][cj][drift][common_offset] += mass
                pair_gauge_offset[left_index][right_index][common_offset] += mass
                pair_nonzero = pair_nonzero or mass != 0
            pair_level_nonzero += int(pair_nonzero)

    full_gauge_support = sum(
        sum(gauge_offset[left][right][drift]) != 0
        for left in range(3)
        for right in range(3)
        for drift in range(P)
    )
    active_gauge_support = sum(
        sum(gauge_offset[chamber_index[left]][chamber_index[right]][drift]) != 0
        for left, right in CORNERS
        for drift in range(P)
    )
    require((full_gauge_support, active_gauge_support) == (72, 48),
            (full_gauge_support, active_gauge_support))
    pair_gauge_offset_tuple = tuple(
        tuple(tuple(offsets) for offsets in row) for row in pair_gauge_offset
    )
    require(
        digest(pair_gauge_offset_tuple) == FINAL_COMPARISON_SOURCE_PAIR_SHA256,
        (digest(pair_gauge_offset_tuple), FINAL_COMPARISON_SOURCE_PAIR_SHA256),
    )

    # Independent-gauge hostile: every drift occurs if the two offsets differ.
    independent_drifts = {
        ((q + c_right) - (u + c_left)) % P
        for u in range(P)
        for q in range(P)
        for c_left in range(P)
        for c_right in range(P)
    }
    require(independent_drifts == set(range(P)), independent_drifts)

    # Summing all atom pairs, not imposing a gauge diagonal, recovers THM-2594's P[u,q].
    all_pair_controls = []
    for u in range(P):
        for q in range(P):
            left_sum = M.extract_window(full_u[0], full_u[1], u, P, M.T_DEN)
            right_sum = M.extract_window(full_v[0], full_v[1], q, P, M.T_DEN)
            all_pair_controls.append(
                M.product_cum(left_sum[0], left_sum[1], right_sum[0], right_sum[1], M.T_DEN)[3]
                == pair_mass[u][q]
            )
    require(all(all_pair_controls), "pre-transfer atom restoration did not recover pair table")

    return (
        tuple(
            tuple(tuple(tuple(offsets) for offsets in drifts) for drifts in right_rows)
            for right_rows in gauge_offset
        ),
        pair_gauge_offset_tuple,
        denominator,
        tuple(restoration),
        len(f_pieces),
        tuple(len(group) for group in f_groups),
        tuple(len(group) for group in e_groups),
        middle_nonzero_profiles,
        pair_level_nonzero,
        (full_gauge_support, active_gauge_support),
    )


def endpoint_worker(alpha: int) -> tuple[object, ...]:
    A = load_module(ATOM_PATH, f"thm3514_endpoint_atom_worker_{alpha}", ATOM_SHA256)
    return A.worker(alpha)


def endpoint_pair_function(chunks: tuple[tuple[object, ...], ...]) -> tuple[object, ...]:
    A = load_module(ATOM_PATH, "thm3514_endpoint_for_nonlinear_audit", ATOM_SHA256)
    require(tuple(A.ATOMS) == ATOMS, "endpoint atom order")
    (_word, _t_den, nn, prime, root, zeta, *_rest) = A.context()
    tables = tuple(row for chunk in chunks for row in chunk[1])
    require(tuple(chunk[0] for chunk in chunks) == tuple(range(P)), "worker order")
    require(len(tables) == P * P, len(tables))
    raw_h = [[0] * len(ATOMS) for _ in ATOMS]
    raw_q5 = [[0] * len(ATOMS) for _ in ATOMS]
    table_index = 0
    for alpha in range(P):
        for beta in range(P):
            stored_beta, ax_values, by_values = tables[table_index]
            table_index += 1
            require(stored_beta == beta, (alpha, beta, stored_beta))
            for tau in range(P):
                left_active = tuple(
                    index
                    for index, value in enumerate(ax_values)
                    if value and A.safe(ATOMS[index][1], ATOMS[index][0] + tau)
                )
                right_active = tuple(
                    index
                    for index, value in enumerate(by_values)
                    if value and A.safe(ATOMS[index][1], ATOMS[index][0] + tau)
                )
                phase_h = pow(zeta, (beta - alpha - tau) % P, prime)
                phase_q5 = pow(zeta, (beta - alpha) % P, prime)
                for left_index in left_active:
                    left_value = ax_values[left_index]
                    for right_index in right_active:
                        product_value = left_value * by_values[right_index] % prime
                        raw_h[left_index][right_index] = (
                            raw_h[left_index][right_index] + phase_h * product_value
                        ) % prime
                        raw_q5[left_index][right_index] = (
                            raw_q5[left_index][right_index] + phase_q5 * product_value
                        ) % prime
    normalizer = pow(P**3, -1, prime)
    bridge = tuple(
        tuple((raw_h[i][j] - raw_q5[i][j]) * normalizer % prime for j in range(len(ATOMS)))
        for i in range(len(ATOMS))
    )
    return nn, prime, root, zeta, bridge


def dft(values: tuple[int, ...], zeta: int, prime: int) -> tuple[int, ...]:
    require(len(values) == P, len(values))
    return tuple(
        sum(value * pow(zeta, (-frequency * index) % P, prime) for index, value in enumerate(values))
        % prime
        for frequency in range(P)
    )


def walsh(rows: tuple[tuple[int, ...], ...], prime: int) -> tuple[tuple[int, ...], ...]:
    require(len(rows) == 4 and all(len(row) == P for row in rows), "K4 rows")
    return (
        tuple(sum(rows[index][d] for index in range(4)) % prime for d in range(P)),
        tuple((rows[0][d] + rows[1][d] - rows[2][d] - rows[3][d]) % prime for d in range(P)),
        tuple((rows[0][d] - rows[1][d] + rows[2][d] - rows[3][d]) % prime for d in range(P)),
        tuple((rows[0][d] - rows[1][d] - rows[2][d] + rows[3][d]) % prime for d in range(P)),
    )


def spectrum(rows: tuple[tuple[int, ...], ...], zeta: int, prime: int) -> tuple[tuple[int, ...], ...]:
    return tuple(dft(row, zeta, prime) for row in walsh(rows, prime))


def rref_rank(matrix: tuple[tuple[int, ...], ...], prime: int) -> int:
    rows = [list(value % prime for value in row) for row in matrix]
    if not rows:
        return 0
    columns = len(rows[0])
    require(all(len(row) == columns for row in rows), "ragged matrix")
    rank = 0
    for column in range(columns):
        pivot = next((row for row in range(rank, len(rows)) if rows[row][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][column], -1, prime)
        rows[rank] = [value * inverse % prime for value in rows[rank]]
        for row in range(len(rows)):
            if row == rank or rows[row][column] == 0:
                continue
            factor = rows[row][column]
            rows[row] = [(left - factor * right) % prime for left, right in zip(rows[row], rows[rank])]
        rank += 1
        if rank == len(rows):
            break
    return rank


def projective_system(
    source: tuple[tuple[int, ...], ...],
    target: tuple[tuple[int, ...], ...],
    prime: int,
) -> tuple[int, int]:
    """Unknowns are a 4x4 mixer followed by thirteen projective multipliers."""
    require(len(source) == len(target) == 4, "four spectral channels")
    equations = []
    for frequency in range(P):
        for output in range(4):
            equation = [0] * 29
            for column in range(4):
                equation[4 * output + column] = source[column][frequency]
            equation[16 + frequency] = -target[output][frequency] % prime
            equations.append(tuple(equation))
    rank = rref_rank(tuple(equations), prime)
    return rank, 29 - rank


def main() -> None:
    with ProcessPoolExecutor(max_workers=5) as pool:
        source_future = pool.submit(source_worker)
        endpoint_futures = tuple(pool.submit(endpoint_worker, alpha) for alpha in range(P))
        gauge_offset, pair_gauge_offset, denominator, restoration, f_piece_count, f_group_counts, e_group_counts, middle_nonzero_profiles, pair_level_nonzero, gauge_supports = (
            source_future.result()
        )
        endpoint_chunks = tuple(future.result() for future in endpoint_futures)

    nn, prime, root, zeta, endpoint_pair = endpoint_pair_function(endpoint_chunks)
    require(denominator % prime != 0, (denominator, prime))
    inverse_denominator = pow(denominator, -1, prime)
    require(pow(zeta, P, prime) == 1 and zeta != 1, "primitive root")

    chamber_index = {name: index for index, name in enumerate(CHAMBER_NAMES)}
    source_common = tuple(
        tuple(
            tuple(
                gauge_offset[chamber_index[left]][chamber_index[right]][drift][c]
                for drift in range(P)
            )
            for left, right in CORNERS
        )
        for c in range(P)
    )

    source_owner_rows = []
    source_owner_spectra = []
    for owner_character in range(P):
        rows = tuple(
            tuple(
                sum(
                    source_common[c][corner][drift]
                    * pow(zeta, (-owner_character * c) % P, prime)
                    for c in range(P)
                )
                % prime * inverse_denominator % prime
                for drift in range(P)
            )
            for corner in range(4)
        )
        source_owner_rows.append(rows)
        source_owner_spectra.append(spectrum(rows, zeta, prime))

    source_support = tuple(
        sum(source_owner_rows[0][corner][drift] != 0 for corner in range(4) for drift in range(P))
        for _once in range(1)
    )[0]
    require(source_support == 48, source_support)
    require(all(source_owner_rows[0][corner][0] == 0 for corner in range(4)), "same-root source bucket")
    source_k0_counts = tuple(sum(value != 0 for value in row) for row in source_owner_spectra[0])
    require(source_k0_counts == (13, 12, 12, 13), source_k0_counts)
    owner_full_counts = tuple(
        tuple(sum(value != 0 for value in row) for row in source_owner_spectra[k])
        for k in range(1, P)
    )
    require(all(counts == (13, 13, 13, 13) for counts in owner_full_counts), owner_full_counts)

    endpoint_full_rows = [[0] * P for _ in range(4)]
    for left_index, (left_sheet, left_chamber) in enumerate(ATOMS):
        if left_chamber not in ACTIVE:
            continue
        for right_index, (right_sheet, right_chamber) in enumerate(ATOMS):
            if right_chamber not in ACTIVE:
                continue
            corner = CORNERS.index((left_chamber, right_chamber))
            drift = (right_sheet - left_sheet) % P
            endpoint_full_rows[corner][drift] = (
                endpoint_full_rows[corner][drift] + endpoint_pair[left_index][right_index]
            ) % prime
    endpoint_full_rows_tuple = tuple(tuple(row) for row in endpoint_full_rows)
    endpoint_restricted_rows = tuple(
        tuple(0 if drift == 0 else value for drift, value in enumerate(row))
        for row in endpoint_full_rows_tuple
    )
    endpoint_full_total = sum(sum(row) for row in endpoint_full_rows_tuple) % prime
    endpoint_restricted_total = sum(sum(row) for row in endpoint_restricted_rows) % prime
    require(endpoint_full_total == EXPECTED_ENDPOINT_FULL, endpoint_full_total)
    require(endpoint_restricted_total == EXPECTED_ENDPOINT_DNEQ, endpoint_restricted_total)
    endpoint_full_spectrum = spectrum(endpoint_full_rows_tuple, zeta, prime)
    endpoint_restricted_spectrum = spectrum(endpoint_restricted_rows, zeta, prime)
    require(all(all(row) for row in endpoint_full_spectrum), "full endpoint spectrum")
    require(all(all(row) for row in endpoint_restricted_spectrum), "d!=0 endpoint spectrum")
    require(digest(endpoint_pair) == FINAL_COMPARISON_TARGET_PAIR_SHA256,
            (digest(endpoint_pair), FINAL_COMPARISON_TARGET_PAIR_SHA256))

    pulled_by_c = [0] * P
    pulled_rows_by_c = [[[0] * P for _ in range(4)] for _ in range(P)]
    for left_index, (left_sheet, left_chamber) in enumerate(ATOMS):
        for right_index, (right_sheet, right_chamber) in enumerate(ATOMS):
            endpoint_value = endpoint_pair[left_index][right_index]
            drift = (right_sheet - left_sheet) % P
            corner = CORNERS.index((left_chamber, right_chamber)) if (
                left_chamber in ACTIVE and right_chamber in ACTIVE
            ) else None
            for c in range(P):
                value = (
                    pair_gauge_offset[left_index][right_index][c]
                    % prime * inverse_denominator % prime * endpoint_value % prime
                )
                pulled_by_c[c] = (pulled_by_c[c] + value) % prime
                if corner is not None:
                    pulled_rows_by_c[c][corner][drift] = (
                        pulled_rows_by_c[c][corner][drift] + value
                    ) % prime
    pulled_owner_characters = dft(tuple(pulled_by_c), zeta, prime)
    require(all(pulled_owner_characters), pulled_owner_characters)

    pulled_owner_spectra = []
    for owner_character in range(P):
        pulled_rows = tuple(
            tuple(
                sum(
                    pulled_rows_by_c[c][corner][drift]
                    * pow(zeta, (-owner_character * c) % P, prime)
                    for c in range(P)
                )
                % prime
                for drift in range(P)
            )
            for corner in range(4)
        )
        pulled_owner_spectra.append(spectrum(pulled_rows, zeta, prime))

    systems = {
        target_name: tuple(
            projective_system(source_owner_spectra[owner_character], target, prime)
            for owner_character in range(P)
        )
        for target_name, target in (
            ("endpoint_full", endpoint_full_spectrum),
            ("endpoint_dneq", endpoint_restricted_spectrum),
        )
    }
    require(
        all(result == (29, 0) for census in systems.values() for result in census),
        systems,
    )

    proof = (
        "39_guard_atoms_partition_f_and_e_before_P_13^5_transfer",
        "atom_sums_restore_every_THM2594_root_profile_and_pair_mass",
        "one_common_torsor_offset_a=u+c_b=q+c_iff_root_shift_equals_minus_drift",
        "same_root_disjointness_removes_exactly_four_owner_active_buckets",
        "endpoint_pair_function_rebuilt_from_THM3514_atom_workers_before_drift_grouping",
        "finite_scalar_pair_pullback_retains_owner_offset_but_is_not_temporal_current",
        "29_variable_systems_retain_4x4_mixer_and_all_13_frequency_multipliers",
    )
    boundary = (
        "FINITE_EXACT_SPLIT_FIELD_ADDRESS_PULLBACK_ONLY",
        "not_one_integral_over_a_common_temporal_stalk",
        "no_absolute_H1_bispectrum_scalar_row_exclusion_or_LRC14_claim",
    )
    semantic_surface = (
        PRIMARY_SHA256,
        ATOM_SHA256,
        (prime, root, zeta, nn),
        denominator,
        digest(source_common),
        digest(pair_gauge_offset),
        pair_level_nonzero,
        gauge_supports,
        source_support,
        source_k0_counts,
        owner_full_counts,
        digest(endpoint_pair),
        digest(endpoint_full_rows_tuple),
        digest(endpoint_restricted_rows),
        (endpoint_full_total, endpoint_restricted_total),
        tuple(tuple(sum(value != 0 for value in row) for row in bank) for bank in (endpoint_full_spectrum, endpoint_restricted_spectrum)),
        pulled_owner_characters,
        systems,
        proof,
        boundary,
    )
    semantic_sha256 = hashlib.sha256(repr(semantic_surface).encode("utf-8")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_sha256 == EXPECTED_SEMANTIC_SHA256, (semantic_sha256, EXPECTED_SEMANTIC_SHA256))

    source = Path(__file__).resolve()
    print("R5 NONLINEAR COMMON-GAUGE ENDPOINT CONNECTION -- INDEPENDENT HOSTILE AUDIT")
    print("status=ACCEPT_SCOPED_FINITE_EXACT_PACKAGE;LRC14_OPEN")
    print(f"dependencies=(thm2594={PRIMARY_SHA256},thm3514_atom_engine={ATOM_SHA256})")
    print(f"field=(prime={prime},root={root},zeta13={zeta},order={nn})")
    print(f"source_pretransfer=(f_pieces={f_piece_count},f_atom_piece_counts={f_group_counts},e_atom_piece_counts={e_group_counts},restoration={restoration},middle_nonzero_profiles={middle_nonzero_profiles})")
    print("common_gauge=(a=u+c,b=q+c,therefore_s=u-q=-d;independent_offsets_destroy_relation)")
    print(f"source_common_gauge=(full_support={gauge_supports[0]}/117,owner_active_support={gauge_supports[1]}/52,pair_level_nonzero={pair_level_nonzero}/1521,pair_digest={digest(pair_gauge_offset)})")
    print(f"source_owner_active=(support={source_support}/52,k0_spectral_counts={source_k0_counts},k_nonzero_counts={tuple(sorted(set(owner_full_counts)))},digest={digest(source_common)})")
    print(f"endpoint=(full_bridge={endpoint_full_total},dneq_bridge={endpoint_restricted_total},full_spectral_counts={tuple(sum(value != 0 for value in row) for row in endpoint_full_spectrum)},dneq_spectral_counts={tuple(sum(value != 0 for value in row) for row in endpoint_restricted_spectrum)},pair_digest={digest(endpoint_pair)})")
    print(f"pair_pullback=(owner_character_values={pulled_owner_characters},nonzero={sum(value != 0 for value in pulled_owner_characters)}/13)")
    print(f"projective_systems=(variables=29,census={dict((str(key),(len(value),tuple(sorted(set(value))))) for key,value in systems.items())})")
    print(f"proof={proof}")
    print(f"boundary={boundary}")
    print(f"semantic_sha256={semantic_sha256}")
    print(f"script_sha256_lf={lf_sha256(source)}")
    print("reproducibility=no_candidate_imports;no_randomness;no_elapsed_fields;normal_and_O_transcripts_must_match")
    print("commands=python -B 04-computation/lrc_r5_nonlinear_common_gauge_endpoint_connection_independent_audit_20260816.py;python -B -O 04-computation/lrc_r5_nonlinear_common_gauge_endpoint_connection_independent_audit_20260816.py")
    print("PASS")


if __name__ == "__main__":
    main()
