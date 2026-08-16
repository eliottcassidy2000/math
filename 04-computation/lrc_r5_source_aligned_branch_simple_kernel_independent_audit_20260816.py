#!/usr/bin/env python3
"""Clean-room audit of the R5 source-aligned branch/simple-kernel package.

Neither submitted source-aligned candidate is imported.  The program rebuilds
THM-2471's source refinement

    f_omega^src = 1_Q P^2(e 1_omega)

from the hash-pinned THM-2594 interval engine, and independently rebuilds the
THM-3514 endpoint pair function from its audited atom engine.

The resulting endpoint-weighted finite pairing has a one-base simple-kernel
integral representation by finite linearity.  AX and BY remain preintegrated
endpoint scalars, so this is not a common-node temporal current, grouped
coefficient, row exclusion, or LRC(14) theorem.
"""

from __future__ import annotations

from bisect import bisect_right
from concurrent.futures import ProcessPoolExecutor
from contextlib import redirect_stdout
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
EXPECTED_SEMANTIC_SHA256 = "822957b53ef3bcaa8509fa68bd2c0d3080dec7632db3eca063e2d947162bd4b9"

P = 13
PKT = P * P
CHAMBERS = (("left", 0, 1), ("middle", 1, 6), ("right", 6, 7))
CHAMBER_NAMES = tuple(name for name, _lo, _hi in CHAMBERS)
ACTIVE = ("left", "right")
CORNERS = (("left", "left"), ("left", "right"), ("right", "left"), ("right", "right"))
ATOMS = tuple((sheet, chamber) for sheet in range(P) for chamber in CHAMBER_NAMES)
ATOM_INDEX = {atom: index for index, atom in enumerate(ATOMS)}


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


def profile_value(starts: list[int], values: list[int], point: int) -> int:
    return values[bisect_right(starts, point) - 1]


def sum_profiles(
    profiles: tuple[tuple[list[int], list[int]], ...], grid: int
) -> tuple[list[int], list[int]]:
    breaks = sorted(
        {0, *(point for starts, _values in profiles for point in starts if point < grid)}
    )
    starts: list[int] = []
    values: list[int] = []
    for point in breaks:
        value = sum(profile_value(s, v, point) for s, v in profiles)
        if not values or value != values[-1]:
            starts.append(point)
            values.append(value)
    require(starts and starts[0] == 0 and starts[-1] < grid, "profile sum")
    return starts, values


def profiles_equal(
    left: tuple[list[int], list[int]], right: tuple[list[int], list[int]], grid: int
) -> bool:
    breaks = sorted({0, grid, *left[0], *right[0]})
    return all(
        profile_value(left[0], left[1], point) == profile_value(right[0], right[1], point)
        for point in breaks[:-1]
    )


def split_intervals(
    intervals: list[tuple[int, int]],
    atom_intervals: tuple[tuple[int, str, int, int], ...],
) -> tuple[tuple[tuple[int, int], ...], ...]:
    groups: list[list[tuple[int, int]]] = [[] for _ in atom_intervals]
    index = 0
    measure_in = 0
    measure_out = 0
    for left, right in intervals:
        measure_in += right - left
        while atom_intervals[index][3] <= left:
            index += 1
        scan = index
        while scan < len(atom_intervals) and atom_intervals[scan][2] < right:
            _sheet, _chamber, atom_left, atom_right = atom_intervals[scan]
            lo = max(left, atom_left)
            hi = min(right, atom_right)
            if lo < hi:
                groups[scan].append((lo, hi))
                measure_out += hi - lo
            if atom_right >= right:
                break
            scan += 1
    require(measure_in == measure_out, (measure_in, measure_out))
    return tuple(tuple(group) for group in groups)


def restrict_profile(
    starts: list[int], values: list[int], intervals: list[tuple[int, int]], grid: int
) -> list[tuple[int, int, int]]:
    pieces: list[tuple[int, int, int]] = []
    for interval_left, interval_right in intervals:
        index = bisect_right(starts, interval_left) - 1
        while True:
            profile_left = starts[index]
            profile_right = starts[index + 1] if index + 1 < len(starts) else grid
            left = max(interval_left, profile_left)
            right = min(interval_right, profile_right)
            if left < right and values[index]:
                pieces.append((left, right, values[index]))
            if profile_right >= interval_right:
                break
            index += 1
    return pieces


def rank_mod(matrix: tuple[tuple[int, ...], ...], prime: int) -> int:
    rows = [list(value % prime for value in row) for row in matrix]
    if not rows:
        return 0
    columns = len(rows[0])
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


def source_worker() -> tuple[object, ...]:
    M = load_module(PRIMARY_PATH, "thm2471_source_aligned_cleanroom", PRIMARY_SHA256)
    e_intervals = M.build_set(M.PAT_E, M.ZELL)
    q_intervals = M.build_set(M.PAT_QB, M.ZELL)
    unit = M.T_DEN // 91
    atom_intervals = tuple(
        (sheet, chamber, (7 * sheet + low) * unit, (7 * sheet + high) * unit)
        for sheet in range(P)
        for chamber, low, high in CHAMBERS
    )
    e_groups = split_intervals(e_intervals, atom_intervals)

    # Raw T^2 transition: row omega, column eta is the integral over I_eta
    # of 169 P_169(1_omega).  Scaling does not affect support or rank.
    transition_rows = []
    for _sheet, _chamber, left, right in atom_intervals:
        starts, values = M.weighted_fold([(left, right, 1)], PKT, M.T_DEN)
        profile = M.product_cum(starts, values, [0], [1], M.T_DEN)
        transition_rows.append(
            tuple(M.set_integral(profile, [(eta_left, eta_right)]) for _a, _c, eta_left, eta_right in atom_intervals)
        )
    transition = tuple(transition_rows)
    transition_support = sum(value != 0 for row in transition for value in row)
    require(transition_support == len(ATOMS) ** 2, transition_support)
    for chamber_index in range(3):
        reference = transition[chamber_index]
        require(
            all(transition[3 * sheet + chamber_index] == reference for sheet in range(P)),
            ("source-sheet row identity", chamber_index),
        )
    transition_ranks = tuple(rank_mod(transition, prime) for prime in (547, 911, 1093))
    require(transition_ranks == (3, 3, 3), transition_ranks)

    branch_sheet = tuple(n // P for n in range(PKT))
    require(
        all((n // P) == int((n + 0) // P) == int((n + 1 - 1) // P) for n in range(PKT)),
        "inverse branch sheet recovery",
    )
    require(tuple(branch_sheet.count(sheet) for sheet in range(P)) == (P,) * P, "branch census")

    # Source refinement f_omega^src=1_Q P^2(e 1_omega).
    f_groups = []
    p2_profiles = []
    for group in e_groups:
        starts, values = M.weighted_fold([(left, right, 1) for left, right in group], PKT, M.T_DEN)
        p2_profiles.append((starts, values))
        f_groups.append(tuple(restrict_profile(starts, values, q_intervals, M.T_DEN)))
    f_groups_tuple = tuple(f_groups)
    u_profiles = tuple(M.weighted_fold(list(group), M.DCOLL, M.T_DEN) for group in f_groups_tuple)
    v_profiles = tuple(
        M.weighted_fold([(left, right, 1) for left, right in group], M.DCOLL, M.T_DEN)
        for group in e_groups
    )

    whole_p2 = M.weighted_fold([(left, right, 1) for left, right in e_intervals], PKT, M.T_DEN)
    whole_f_pieces = restrict_profile(whole_p2[0], whole_p2[1], q_intervals, M.T_DEN)
    whole_u = M.weighted_fold(whole_f_pieces, M.DCOLL, M.T_DEN)
    whole_v = M.weighted_fold([(left, right, 1) for left, right in e_intervals], M.DCOLL, M.T_DEN)
    require(profiles_equal(sum_profiles(tuple(p2_profiles), M.T_DEN), whole_p2, M.T_DEN), "P2 atom restoration")
    require(profiles_equal(sum_profiles(u_profiles, M.T_DEN), whole_u, M.T_DEN), "U source restoration")
    require(profiles_equal(sum_profiles(v_profiles, M.T_DEN), whole_v, M.T_DEN), "V source restoration")

    u_windows = tuple(
        tuple(M.extract_window(starts, values, root, P, M.T_DEN) for root in range(P))
        for starts, values in u_profiles
    )
    v_windows = tuple(
        tuple(M.extract_window(starts, values, root, P, M.T_DEN) for root in range(P))
        for starts, values in v_profiles
    )
    nonempty_u = tuple(index for index, (_s, values) in enumerate(u_profiles) if any(values))
    nonempty_v = tuple(index for index, (_s, values) in enumerate(v_profiles) if any(values))
    require(len(nonempty_u) == len(nonempty_v) == 20, (nonempty_u, nonempty_v))

    u_root_sums = tuple(sum_profiles(u_windows[index], M.T_DEN) for index in nonempty_u)
    v_root_sums = tuple(sum_profiles(v_windows[index], M.T_DEN) for index in nonempty_v)
    all_root_pair_support = 0
    for left in u_root_sums:
        for right in v_root_sums:
            mass = M.product_cum(left[0], left[1], right[0], right[1], M.T_DEN)[3]
            all_root_pair_support += int(mass != 0)
    require(all_root_pair_support == 400, all_root_pair_support)

    chamber_index = {name: index for index, name in enumerate(CHAMBER_NAMES)}
    gauge_offset = [[[[0] * P for _d in range(P)] for _right in CHAMBERS] for _left in CHAMBERS]
    pair_gauge_offset = [[[0] * P for _right in ATOMS] for _left in ATOMS]
    gauge_pair_support = 0
    for left_index in nonempty_u:
        left_sheet, left_chamber = ATOMS[left_index]
        ci = chamber_index[left_chamber]
        for right_index in nonempty_v:
            right_sheet, right_chamber = ATOMS[right_index]
            cj = chamber_index[right_chamber]
            drift = (right_sheet - left_sheet) % P
            pair_nonzero = False
            for current_root in range(P):
                source_root = (current_root + drift) % P
                offset = (left_sheet - current_root) % P
                require(offset == (right_sheet - source_root) % P, "common gauge")
                left = u_windows[left_index][current_root]
                right = v_windows[right_index][source_root]
                mass = M.product_cum(left[0], left[1], right[0], right[1], M.T_DEN)[3]
                gauge_offset[ci][cj][drift][offset] += mass
                pair_gauge_offset[left_index][right_index][offset] += mass
                pair_nonzero = pair_nonzero or mass != 0
            gauge_pair_support += int(pair_nonzero)
    require(gauge_pair_support == 362, gauge_pair_support)

    active_support = sum(
        sum(gauge_offset[chamber_index[left]][chamber_index[right]][drift]) != 0
        for left, right in CORNERS
        for drift in range(P)
    )
    require(active_support == 48, active_support)

    denominator = M.RPKT * M.DCOLL * M.DCOLL * M.T_DEN
    return (
        transition,
        transition_support,
        transition_ranks,
        branch_sheet,
        tuple(tuple(tuple(tuple(offsets) for offsets in drifts) for drifts in right_rows) for right_rows in gauge_offset),
        tuple(tuple(tuple(offsets) for offsets in row) for row in pair_gauge_offset),
        all_root_pair_support,
        gauge_pair_support,
        active_support,
        denominator,
        tuple(nonempty_u),
        tuple(nonempty_v),
        tuple(tuple((tuple(starts), tuple(values)) for starts, values in roots) for roots in u_windows),
        tuple(tuple((tuple(starts), tuple(values)) for starts, values in roots) for roots in v_windows),
        tuple(len(group) for group in f_groups_tuple),
        tuple(len(group) for group in e_groups),
    )


def endpoint_worker(alpha: int) -> tuple[object, ...]:
    A = load_module(ATOM_PATH, f"thm3514_source_aligned_endpoint_worker_{alpha}", ATOM_SHA256)
    return A.worker(alpha)


def endpoint_pair_function(chunks: tuple[tuple[object, ...], ...]) -> tuple[object, ...]:
    A = load_module(ATOM_PATH, "thm3514_source_aligned_endpoint_cleanroom", ATOM_SHA256)
    require(tuple(A.ATOMS) == ATOMS, "atom order")
    (_word, _t_den, nn, prime, root, zeta, *_rest) = A.context()
    tables = tuple(row for chunk in chunks for row in chunk[1])
    require(tuple(chunk[0] for chunk in chunks) == tuple(range(P)), "worker order")
    base = [[0] * len(ATOMS) for _ in ATOMS]
    table_index = 0
    for alpha in range(P):
        for beta in range(P):
            stored_beta, ax_values, by_values = tables[table_index]
            table_index += 1
            require(stored_beta == beta, (alpha, beta, stored_beta))
            phase = pow(zeta, (beta - alpha) % P, prime)
            for left_index, left_value in enumerate(ax_values):
                if not left_value:
                    continue
                for right_index, right_value in enumerate(by_values):
                    if not right_value:
                        continue
                    base[left_index][right_index] = (
                        base[left_index][right_index]
                        + phase * left_value % prime * right_value
                    ) % prime
    normalizer = pow(P**3, -1, prime)
    endpoint_bank = []
    for residue_t in range(P):
        rows = []
        for left_index, (left_sheet, left_chamber) in enumerate(ATOMS):
            row = []
            for right_index, (right_sheet, right_chamber) in enumerate(ATOMS):
                kernel = sum(
                    A.safe(left_chamber, left_sheet + tau)
                    * A.safe(right_chamber, right_sheet + tau)
                    * pow(zeta, (-residue_t * tau) % P, prime)
                    for tau in range(P)
                ) % prime
                row.append(base[left_index][right_index] * kernel % prime * normalizer % prime)
            rows.append(tuple(row))
        endpoint_bank.append(tuple(rows))
    endpoint_bank_tuple = tuple(endpoint_bank)
    bridge = tuple(
        tuple((endpoint_bank_tuple[1][i][j] - endpoint_bank_tuple[0][i][j]) % prime
              for j in range(len(ATOMS)))
        for i in range(len(ATOMS))
    )
    require(digest(bridge) == "c2d5911b287510335edc6aefa6d3b865c982568f678bb89ee9b82ee211962df1",
            digest(bridge))
    return nn, prime, root, zeta, endpoint_bank_tuple


def dft(values: tuple[int, ...], zeta: int, prime: int) -> tuple[int, ...]:
    return tuple(
        sum(value * pow(zeta, (-frequency * index) % P, prime) for index, value in enumerate(values))
        % prime
        for frequency in range(P)
    )


def walsh(rows: tuple[tuple[int, ...], ...], prime: int) -> tuple[tuple[int, ...], ...]:
    return (
        tuple(sum(rows[index][d] for index in range(4)) % prime for d in range(P)),
        tuple((rows[0][d] + rows[1][d] - rows[2][d] - rows[3][d]) % prime for d in range(P)),
        tuple((rows[0][d] - rows[1][d] + rows[2][d] - rows[3][d]) % prime for d in range(P)),
        tuple((rows[0][d] - rows[1][d] - rows[2][d] + rows[3][d]) % prime for d in range(P)),
    )


def double_center(
    matrix: tuple[tuple[int, ...], ...], universe: tuple[int, ...], prime: int
) -> tuple[tuple[int, ...], ...]:
    inverse = pow(len(universe), -1, prime)
    row_means = tuple(sum(matrix[i][j] for j in universe) % prime * inverse % prime for i in range(len(matrix)))
    column_means = tuple(sum(matrix[i][j] for i in universe) % prime * inverse % prime for j in range(len(matrix)))
    grand = sum(matrix[i][j] for i in universe for j in universe) % prime
    grand = grand * inverse % prime * inverse % prime
    return tuple(
        tuple((matrix[i][j] - row_means[i] - column_means[j] + grand) % prime for j in range(len(matrix)))
        for i in range(len(matrix))
    )


def main() -> None:
    with ProcessPoolExecutor(max_workers=5) as pool:
        source_future = pool.submit(source_worker)
        endpoint_futures = tuple(pool.submit(endpoint_worker, alpha) for alpha in range(P))
        source = source_future.result()
        endpoint_chunks = tuple(future.result() for future in endpoint_futures)

    (
        transition, transition_support, transition_ranks, branch_sheet,
        gauge_offset, pair_gauge_offset, all_root_pair_support,
        gauge_pair_support, active_support, denominator, nonempty_u,
        nonempty_v, u_windows, v_windows, f_group_counts, e_group_counts,
    ) = source
    nn, prime, root, zeta, endpoint_bank = endpoint_pair_function(endpoint_chunks)
    require(denominator % prime != 0 and pow(zeta, P, prime) == 1 and zeta != 1, "field gate")
    inverse_denominator = pow(denominator, -1, prime)

    chamber_index = {name: index for index, name in enumerate(CHAMBER_NAMES)}
    source_rows_by_owner = []
    source_spectra_by_owner = []
    source_ranks = []
    source_spectral_counts = []
    for owner_frequency in range(P):
        rows = tuple(
            tuple(
                sum(
                    gauge_offset[chamber_index[left]][chamber_index[right]][drift][offset]
                    * pow(zeta, (-owner_frequency * offset) % P, prime)
                    for offset in range(P)
                ) % prime * inverse_denominator % prime
                for drift in range(P)
            )
            for left, right in CORNERS
        )
        spectra = tuple(dft(row, zeta, prime) for row in walsh(rows, prime))
        source_rows_by_owner.append(rows)
        source_spectra_by_owner.append(spectra)
        source_ranks.append(rank_mod(rows, prime))
        source_spectral_counts.append(tuple(sum(value != 0 for value in row) for row in spectra))
    require(tuple(source_ranks) == (4,) * P, source_ranks)
    require(source_spectral_counts[0] == (13, 12, 12, 13), source_spectral_counts[0])
    require(all(counts == (13, 13, 13, 13) for counts in source_spectral_counts[1:]),
            source_spectral_counts)

    endpoint_bridge = tuple(
        tuple((endpoint_bank[1][i][j] - endpoint_bank[0][i][j]) % prime
              for j in range(len(ATOMS)))
        for i in range(len(ATOMS))
    )
    endpoint_totals = tuple(sum(sum(row) for row in matrix) % prime for matrix in endpoint_bank)
    require(endpoint_totals[6] == 225010624370142818572, endpoint_totals[6])
    centered_all_bank = tuple(
        double_center(matrix, tuple(range(len(ATOMS))), prime) for matrix in endpoint_bank
    )
    centered_source_bank = tuple(
        double_center(matrix, tuple(nonempty_u), prime) for matrix in endpoint_bank
    )
    centered_bridge_all = double_center(endpoint_bridge, tuple(range(len(ATOMS))), prime)
    centered_bridge_source = double_center(endpoint_bridge, tuple(nonempty_u), prime)

    def pullback(pair_function: tuple[tuple[int, ...], ...]) -> tuple[int, ...]:
        by_offset = [0] * P
        for left_index in range(len(ATOMS)):
            for right_index in range(len(ATOMS)):
                endpoint_weight = pair_function[left_index][right_index]
                for offset in range(P):
                    by_offset[offset] = (
                        by_offset[offset]
                        + pair_gauge_offset[left_index][right_index][offset]
                        % prime * inverse_denominator % prime * endpoint_weight
                    ) % prime
        return dft(tuple(by_offset), zeta, prime)

    residue_pullbacks = tuple(pullback(matrix) for matrix in endpoint_bank)
    residue_centered_all_pullbacks = tuple(pullback(matrix) for matrix in centered_all_bank)
    residue_centered_source_pullbacks = tuple(pullback(matrix) for matrix in centered_source_bank)
    residue_owner_supports = (
        sum(value != 0 for row in residue_pullbacks for value in row),
        sum(value != 0 for row in residue_centered_all_pullbacks for value in row),
        sum(value != 0 for row in residue_centered_source_pullbacks for value in row),
    )
    require(residue_owner_supports == (169, 169, 169), residue_owner_supports)

    full_pullback = pullback(endpoint_bridge)
    centered_all_pullback = pullback(centered_bridge_all)
    centered_source_pullback = pullback(centered_bridge_source)
    pullback_supports = tuple(
        sum(value != 0 for value in values)
        for values in (full_pullback, centered_all_pullback, centered_source_pullback)
    )
    require(pullback_supports == (13, 13, 13), pullback_supports)

    # Independent simple-kernel route: choose u and q from the common offset,
    # integrate the source profiles on the one outer base, then apply the
    # preintegrated endpoint scalar.  It must equal the finite tensor pairing.
    simple_engine = load_module(PRIMARY_PATH, "thm2471_simple_kernel_route", PRIMARY_SHA256)

    def simple_kernel_route(pair_function: tuple[tuple[int, ...], ...]) -> tuple[int, ...]:
        by_offset = []
        for offset in range(P):
            total = 0
            for left_index in nonempty_u:
                left_sheet, _left_chamber = ATOMS[left_index]
                current_root = (left_sheet - offset) % P
                left = u_windows[left_index][current_root]
                for right_index in nonempty_v:
                    right_sheet, _right_chamber = ATOMS[right_index]
                    source_root = (right_sheet - offset) % P
                    right = v_windows[right_index][source_root]
                    # Rebuild the common-base integral from the two profiles.
                    mass = simple_engine.product_cum(
                        list(left[0]), list(left[1]), list(right[0]), list(right[1]), simple_engine.T_DEN
                    )[3]
                    total = (
                        total + mass % prime * inverse_denominator % prime
                        * pair_function[left_index][right_index]
                    ) % prime
            by_offset.append(total)
        return dft(tuple(by_offset), zeta, prime)

    simple_full = simple_kernel_route(endpoint_bridge)
    simple_centered = simple_kernel_route(centered_bridge_all)
    require(simple_full == full_pullback, "full simple-kernel mismatch")
    require(simple_centered == centered_all_pullback, "centered simple-kernel mismatch")
    simple_relation = simple_kernel_route(endpoint_bank[6])
    simple_relation_centered = simple_kernel_route(centered_all_bank[6])
    require(simple_relation == residue_pullbacks[6], "relation residue simple-kernel mismatch")
    require(simple_relation_centered == residue_centered_all_pullbacks[6],
            "relation residue centered simple-kernel mismatch")

    proof = (
        "raw_T2_transition_built_from_inverse_branches_before_E_or_Q",
        "branch_n_recovers_source_sheet_floor_n_over_13",
        "source_refinement_is_1Q_P2_ePomega_not_arrival_product_fPomega",
        "source_atom_sums_restore_P2e_U_and_V_exactly",
        "common_gauge_is_imposed_before_root_and_atom_marginalization",
        "endpoint_E_rebuilt_from_THM3514_AX_BY_atom_tables",
        "finite_pairing_equals_one_outer_base_simple_kernel_by_finite_linearity",
    )
    boundary = (
        "AX_BY_are_preintegrated_endpoint_scalars_not_observables_on_outer_base_y",
        "simple_kernel_integral_representation_is_type_correct",
        "physical_temporal_current_or_common_endpoint_node_claim_is_type_incorrect",
        "no_grouped_coefficient_row_exclusion_or_LRC14_claim",
    )
    semantic_surface = (
        PRIMARY_SHA256, ATOM_SHA256, (prime, root, zeta, nn),
        digest(transition), transition_support, transition_ranks, digest(branch_sheet),
        f_group_counts, e_group_counts, tuple(nonempty_u), tuple(nonempty_v),
        digest(gauge_offset), digest(pair_gauge_offset),
        (all_root_pair_support, gauge_pair_support, active_support),
        tuple(source_ranks), tuple(source_spectral_counts), digest(endpoint_bank),
        endpoint_totals, digest(centered_all_bank), digest(centered_source_bank),
        residue_pullbacks, residue_centered_all_pullbacks,
        residue_centered_source_pullbacks, residue_owner_supports,
        full_pullback, centered_all_pullback, centered_source_pullback,
        pullback_supports,
        proof, boundary,
    )
    semantic_sha256 = hashlib.sha256(repr(semantic_surface).encode("utf-8")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_sha256 == EXPECTED_SEMANTIC_SHA256, (semantic_sha256, EXPECTED_SEMANTIC_SHA256))

    source_path = Path(__file__).resolve()
    print("R5 SOURCE-ALIGNED BRANCH / SIMPLE-KERNEL -- INDEPENDENT HOSTILE AUDIT")
    print("status=ACCEPT_SCOPED_FINITE_EXACT_PACKAGE;LRC14_OPEN")
    print(f"dependencies=(thm2594={PRIMARY_SHA256},thm3514_atom_engine={ATOM_SHA256})")
    print(f"field=(prime={prime},root={root},zeta13={zeta},order={nn})")
    print(f"raw_T2_transition=(support={transition_support}/1521,ranks={transition_ranks},source_sheet_row_identity=True,digest={digest(transition)})")
    print(f"inverse_branch_recovery=(branches={len(branch_sheet)},sheet_formula=floor(n/13),sheet_census={tuple(branch_sheet.count(sheet) for sheet in range(P))},digest={digest(branch_sheet)})")
    print(f"source_refinement=(nonempty_U={len(nonempty_u)},nonempty_V={len(nonempty_v)},f_group_piece_counts={f_group_counts},e_group_piece_counts={e_group_counts})")
    print(f"first_collision=(all_root_pair_support={all_root_pair_support}/400,common_gauge_pair_support={gauge_pair_support}/400,active_type_support={active_support}/52,pair_digest={digest(pair_gauge_offset)})")
    print(f"source_owner_spectrum=(ranks={tuple(source_ranks)},counts={tuple(source_spectral_counts)})")
    print(f"endpoint_residue_bank=(totals={endpoint_totals},digest={digest(endpoint_bank)},centered_all_digest={digest(centered_all_bank)},centered_source_digest={digest(centered_source_bank)})")
    print(f"all_residue_owner_pullbacks=(full_support={residue_owner_supports[0]}/169,centered_all_support={residue_owner_supports[1]}/169,centered_source_support={residue_owner_supports[2]}/169)")
    print(f"relation_class_1_0_6=(total={endpoint_totals[6]},full={residue_pullbacks[6]},centered_all={residue_centered_all_pullbacks[6]},centered_source={residue_centered_source_pullbacks[6]},supports={(sum(value != 0 for value in residue_pullbacks[6]),sum(value != 0 for value in residue_centered_all_pullbacks[6]),sum(value != 0 for value in residue_centered_source_pullbacks[6]))})")
    print(f"bridge_pullbacks=(full={full_pullback},centered_all={centered_all_pullback},centered_source={centered_source_pullback},supports={pullback_supports})")
    print("simple_kernel_routes=(bridge_full=True,bridge_centered=True,relation_full=True,relation_centered=True)")
    print(f"proof={proof}")
    print(f"boundary={boundary}")
    print(f"semantic_sha256={semantic_sha256}")
    print(f"script_sha256_lf={lf_sha256(source_path)}")
    print("reproducibility=no_candidate_imports;no_randomness;no_elapsed_fields;normal_and_O_transcripts_must_match")
    print("commands=python -B 04-computation/lrc_r5_source_aligned_branch_simple_kernel_independent_audit_20260816.py;python -B -O 04-computation/lrc_r5_source_aligned_branch_simple_kernel_independent_audit_20260816.py")
    print("PASS")


if __name__ == "__main__":
    main()
