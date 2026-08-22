#!/usr/bin/env python3
"""Clean-room audit of the source-aligned seven-cell residue spectrum.

The submitted seven-cell probe is neither imported nor read during the
reconstruction.  This program splits the THM-2594 source before transfer,
retains its seven owner cells before every marginal, independently rebuilds
the complete THM-3514 (1,0,t) endpoint bank, and audits the two distinct ANOVA
operations: endpoint-pair ANOVA before contraction and cell/owner ANOVA after
contraction.

AX and BY are preintegrated endpoint scalars.  The result is a finite
THM-2334 residue pushforward, not C(a;X,m), a physical current, a THM-2512
bridge, a row exclusion, or LRC(14).
"""

from __future__ import annotations

from bisect import bisect_right
from concurrent.futures import ProcessPoolExecutor
import hashlib
import importlib.util
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
EXPECTED_SEMANTIC_SHA256 = "a0534b8d995956c126ce204117a9488b222ad25dd8b483b269e37740ffe13ccb"

P = 13
Q = 7
PKT = P * P
CHAMBERS = (("left", 0, 1), ("middle", 1, 6), ("right", 6, 7))
CHAMBER_NAMES = tuple(name for name, _low, _high in CHAMBERS)
ATOMS = tuple((sheet, chamber) for sheet in range(P) for chamber in CHAMBER_NAMES)
ALL_ADDRESSES = tuple(range(len(ATOMS)))


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
        profile_value(left[0], left[1], point)
        == profile_value(right[0], right[1], point)
        for point in breaks[:-1]
    )


def split_intervals(
    intervals: list[tuple[int, int]],
    atom_intervals: tuple[tuple[int, str, int, int], ...],
) -> tuple[tuple[tuple[int, int], ...], ...]:
    groups: list[list[tuple[int, int]]] = [[] for _ in atom_intervals]
    atom_index = 0
    input_measure = 0
    output_measure = 0
    for interval_left, interval_right in intervals:
        input_measure += interval_right - interval_left
        while atom_intervals[atom_index][3] <= interval_left:
            atom_index += 1
        scan = atom_index
        while scan < len(atom_intervals) and atom_intervals[scan][2] < interval_right:
            _sheet, _chamber, atom_left, atom_right = atom_intervals[scan]
            left = max(interval_left, atom_left)
            right = min(interval_right, atom_right)
            if left < right:
                groups[scan].append((left, right))
                output_measure += right - left
            if atom_right >= interval_right:
                break
            scan += 1
    require(input_measure == output_measure, (input_measure, output_measure))
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


def source_worker() -> tuple[object, ...]:
    M = load_module(PRIMARY_PATH, "thm2594_seven_cell_cleanroom", PRIMARY_SHA256)

    cells = tuple(tuple(intervals) for intervals in M.build_cells())
    cell_measures = tuple(M.measure(list(cell)) for cell in cells)
    require(len(cells) == Q and cell_measures == (M.T_DEN // Q,) * Q, cell_measures)
    require(
        M.merge_touch([interval for cell in cells for interval in cell]) == [(0, M.T_DEN)],
        "seven cells do not tile the base",
    )
    for left in range(Q):
        for right in range(left + 1, Q):
            require(M.intersect_lists(list(cells[left]), list(cells[right])) == [],
                    ("cell overlap", left, right))

    e_intervals = M.build_set(M.PAT_E, M.ZELL)
    q_intervals = M.build_set(M.PAT_QB, M.ZELL)
    unit = M.T_DEN // 91
    atom_intervals = tuple(
        (sheet, chamber, (7 * sheet + low) * unit, (7 * sheet + high) * unit)
        for sheet in range(P)
        for chamber, low, high in CHAMBERS
    )
    e_groups = split_intervals(e_intervals, atom_intervals)

    p2_profiles = []
    f_groups = []
    for group in e_groups:
        profile = M.weighted_fold([(left, right, 1) for left, right in group], PKT, M.T_DEN)
        p2_profiles.append(profile)
        f_groups.append(tuple(restrict_profile(profile[0], profile[1], q_intervals, M.T_DEN)))
    u_profiles = tuple(M.weighted_fold(list(group), M.DCOLL, M.T_DEN) for group in f_groups)
    v_profiles = tuple(
        M.weighted_fold([(left, right, 1) for left, right in group], M.DCOLL, M.T_DEN)
        for group in e_groups
    )

    whole_p2 = M.weighted_fold([(left, right, 1) for left, right in e_intervals], PKT, M.T_DEN)
    whole_f = restrict_profile(whole_p2[0], whole_p2[1], q_intervals, M.T_DEN)
    whole_u = M.weighted_fold(whole_f, M.DCOLL, M.T_DEN)
    whole_v = M.weighted_fold([(left, right, 1) for left, right in e_intervals], M.DCOLL, M.T_DEN)
    require(profiles_equal(sum_profiles(tuple(p2_profiles), M.T_DEN), whole_p2, M.T_DEN),
            "P2 source restoration")
    require(profiles_equal(sum_profiles(u_profiles, M.T_DEN), whole_u, M.T_DEN),
            "U source restoration")
    require(profiles_equal(sum_profiles(v_profiles, M.T_DEN), whole_v, M.T_DEN),
            "V source restoration")

    u_windows = tuple(
        tuple(M.extract_window(starts, values, root, P, M.T_DEN) for root in range(P))
        for starts, values in u_profiles
    )
    v_windows = tuple(
        tuple(M.extract_window(starts, values, root, P, M.T_DEN) for root in range(P))
        for starts, values in v_profiles
    )
    nonempty_u = tuple(index for index, (_starts, values) in enumerate(u_profiles) if any(values))
    nonempty_v = tuple(index for index, (_starts, values) in enumerate(v_profiles) if any(values))
    require(len(nonempty_u) == len(nonempty_v) == 20, (nonempty_u, nonempty_v))

    # Literal coordinate order is [left atom][right atom][common offset][cell].
    ledger = [[[[0] * Q for _offset in range(P)] for _right in ATOMS] for _left in ATOMS]
    marginal = [[[0] * P for _right in ATOMS] for _left in ATOMS]
    for left_index in nonempty_u:
        left_sheet, _left_chamber = ATOMS[left_index]
        for right_index in nonempty_v:
            right_sheet, _right_chamber = ATOMS[right_index]
            for offset in range(P):
                current_root = (left_sheet - offset) % P
                source_root = (right_sheet - offset) % P
                left = u_windows[left_index][current_root]
                right = v_windows[right_index][source_root]
                product = M.product_cum(left[0], left[1], right[0], right[1], M.T_DEN)
                total = product[3]
                by_cell = tuple(M.set_integral(product, list(cell)) for cell in cells)
                require(sum(by_cell) == total, ("cell partition on pair", left_index, right_index, offset))
                marginal[left_index][right_index][offset] = total
                ledger[left_index][right_index][offset] = list(by_cell)

    ledger_tuple = tuple(
        tuple(tuple(tuple(cells_for_offset) for cells_for_offset in offsets) for offsets in row)
        for row in ledger
    )
    marginal_tuple = tuple(tuple(tuple(offsets) for offsets in row) for row in marginal)
    require(
        all(
            sum(ledger_tuple[i][j][offset]) == marginal_tuple[i][j][offset]
            for i in range(len(ATOMS)) for j in range(len(ATOMS)) for offset in range(P)
        ),
        "cell marginal recovery",
    )
    atom_pair_support = sum(
        any(ledger_tuple[i][j][offset][cell] != 0 for offset in range(P) for cell in range(Q))
        for i in range(len(ATOMS)) for j in range(len(ATOMS))
    )
    entry_support = sum(
        value != 0
        for row in ledger_tuple for offsets in row for cells_for_offset in offsets
        for value in cells_for_offset
    )
    require(atom_pair_support == 362, atom_pair_support)
    require(entry_support == 5150, entry_support)
    require(digest(marginal_tuple) ==
            "96337f74f2599044870e90878b8bcfb2ce4a04dc970c88b28f23a4f237b3ea53",
            digest(marginal_tuple))

    denominator = M.RPKT * M.DCOLL * M.DCOLL * M.T_DEN
    return (
        cells, cell_measures, ledger_tuple, marginal_tuple, atom_pair_support,
        entry_support, denominator, tuple(nonempty_u), tuple(nonempty_v),
    )


def endpoint_worker(alpha: int) -> tuple[object, ...]:
    A = load_module(ATOM_PATH, f"thm3514_seven_cell_endpoint_worker_{alpha}", ATOM_SHA256)
    return A.worker(alpha)


def endpoint_bank(chunks: tuple[tuple[object, ...], ...]) -> tuple[object, ...]:
    A = load_module(ATOM_PATH, "thm3514_seven_cell_endpoint_cleanroom", ATOM_SHA256)
    require(tuple(A.ATOMS) == ATOMS, "atom order")
    (_word, _t_den, nn, prime, root, zeta13, *_rest) = A.context()
    tables = tuple(row for chunk in chunks for row in chunk[1])
    require(tuple(chunk[0] for chunk in chunks) == tuple(range(P)), "worker order")

    base = [[0] * len(ATOMS) for _ in ATOMS]
    table_index = 0
    for alpha in range(P):
        for beta in range(P):
            stored_beta, ax_values, by_values = tables[table_index]
            table_index += 1
            require(stored_beta == beta, (alpha, beta, stored_beta))
            phase = pow(zeta13, (beta - alpha) % P, prime)
            for left_index, left_value in enumerate(ax_values):
                if not left_value:
                    continue
                for right_index, right_value in enumerate(by_values):
                    if right_value:
                        base[left_index][right_index] = (
                            base[left_index][right_index]
                            + phase * left_value % prime * right_value
                        ) % prime

    normalizer = pow(P**3, -1, prime)
    bank = []
    for residue in range(P):
        matrix = []
        for left_index, (left_sheet, left_chamber) in enumerate(ATOMS):
            row = []
            for right_index, (right_sheet, right_chamber) in enumerate(ATOMS):
                guard_kernel = sum(
                    A.safe(left_chamber, left_sheet + tau)
                    * A.safe(right_chamber, right_sheet + tau)
                    * pow(zeta13, (-residue * tau) % P, prime)
                    for tau in range(P)
                ) % prime
                row.append(base[left_index][right_index] * guard_kernel % prime * normalizer % prime)
            matrix.append(tuple(row))
        bank.append(tuple(matrix))
    bank_tuple = tuple(bank)
    require(digest(bank_tuple) ==
            "c28119c8b54f47e5b7a46f1508fbba604b0e3997eaadb05b03ad28edd9aed468",
            digest(bank_tuple))
    return nn, prime, root, zeta13, bank_tuple


def matrix_add(matrices: tuple[tuple[tuple[int, ...], ...], ...], prime: int):
    return tuple(
        tuple(sum(matrix[row][column] for matrix in matrices) % prime
              for column in range(len(matrices[0][0])))
        for row in range(len(matrices[0]))
    )


def pair_anova(matrix: tuple[tuple[int, ...], ...], prime: int):
    size = len(matrix)
    inverse = pow(size, -1, prime)
    row_means = tuple(sum(row) % prime * inverse % prime for row in matrix)
    column_means = tuple(
        sum(matrix[row][column] for row in range(size)) % prime * inverse % prime
        for column in range(size)
    )
    grand = sum(sum(row) for row in matrix) % prime * inverse % prime * inverse % prime
    flat = tuple(tuple(grand for _column in range(size)) for _row in range(size))
    left_only = tuple(
        tuple((row_means[row] - grand) % prime for _column in range(size))
        for row in range(size)
    )
    right_only = tuple(
        tuple((column_means[column] - grand) % prime for column in range(size))
        for _row in range(size)
    )
    interaction = tuple(
        tuple((matrix[row][column] - row_means[row] - column_means[column] + grand) % prime
              for column in range(size))
        for row in range(size)
    )
    components = (flat, left_only, right_only, interaction)
    require(matrix_add(components, prime) == matrix, "endpoint ANOVA reconstruction")
    require(all(sum(row) % prime == 0 for row in interaction), "interaction row sums")
    require(all(sum(interaction[row][column] for row in range(size)) % prime == 0
                for column in range(size)), "interaction column sums")
    return components


def contract(ledger, pair_function, inverse_denominator: int, prime: int):
    output = [[0] * P for _cell in range(Q)]
    for left in range(len(ATOMS)):
        for right in range(len(ATOMS)):
            weight = pair_function[left][right]
            if not weight:
                continue
            for offset in range(P):
                for cell in range(Q):
                    output[cell][offset] = (
                        output[cell][offset]
                        + ledger[left][right][offset][cell] % prime
                        * inverse_denominator % prime * weight
                    ) % prime
    return tuple(tuple(row) for row in output)


def output_add(outputs: tuple[tuple[tuple[int, ...], ...], ...], prime: int):
    return tuple(
        tuple(sum(output[cell][owner] for output in outputs) % prime for owner in range(P))
        for cell in range(Q)
    )


def dft2(output, zeta7: int, zeta13: int, prime: int):
    return tuple(
        tuple(
            sum(
                output[cell][owner]
                * pow(zeta7, (-cell_frequency * cell) % Q, prime)
                * pow(zeta13, (-owner_frequency * owner) % P, prime)
                for cell in range(Q) for owner in range(P)
            ) % prime
            for owner_frequency in range(P)
        )
        for cell_frequency in range(Q)
    )


def nonzero_coordinates(spectrum) -> frozenset[tuple[int, int]]:
    return frozenset(
        (cell_frequency, owner_frequency)
        for cell_frequency in range(Q) for owner_frequency in range(P)
        if spectrum[cell_frequency][owner_frequency] != 0
    )


def spectrum_shape(spectrum) -> tuple[int, int, int, int, int]:
    dc = int(spectrum[0][0] != 0)
    cell_axis = sum(spectrum[frequency][0] != 0 for frequency in range(1, Q))
    owner_axis = sum(spectrum[0][frequency] != 0 for frequency in range(1, P))
    mixed = sum(
        spectrum[cell_frequency][owner_frequency] != 0
        for cell_frequency in range(1, Q) for owner_frequency in range(1, P)
    )
    return dc + cell_axis + owner_axis + mixed, dc, cell_axis, owner_axis, mixed


def output_anova(output, prime: int):
    inverse_cells = pow(Q, -1, prime)
    inverse_owners = pow(P, -1, prime)
    row_means = tuple(sum(row) % prime * inverse_owners % prime for row in output)
    column_means = tuple(
        sum(output[cell][owner] for cell in range(Q)) % prime * inverse_cells % prime
        for owner in range(P)
    )
    grand = sum(sum(row) for row in output) % prime * inverse_cells % prime * inverse_owners % prime
    flat = tuple(tuple(grand for _owner in range(P)) for _cell in range(Q))
    cell_only = tuple(
        tuple((row_means[cell] - grand) % prime for _owner in range(P))
        for cell in range(Q)
    )
    owner_only = tuple(
        tuple((column_means[owner] - grand) % prime for owner in range(P))
        for _cell in range(Q)
    )
    mixed = tuple(
        tuple((output[cell][owner] - row_means[cell] - column_means[owner] + grand) % prime
              for owner in range(P))
        for cell in range(Q)
    )
    components = (flat, cell_only, owner_only, mixed)
    require(output_add(components, prime) == output, "output ANOVA reconstruction")
    return components


def main() -> None:
    with ProcessPoolExecutor(max_workers=5) as pool:
        source_future = pool.submit(source_worker)
        endpoint_futures = tuple(pool.submit(endpoint_worker, alpha) for alpha in range(P))
        source = source_future.result()
        endpoint_chunks = tuple(future.result() for future in endpoint_futures)

    (
        cells, cell_measures, ledger, marginal, atom_pair_support,
        entry_support, denominator, nonempty_u, nonempty_v,
    ) = source
    nn, prime, root, zeta13, bank = endpoint_bank(endpoint_chunks)
    zeta7 = pow(root, nn // Q, prime)
    require(pow(zeta7, Q, prime) == 1 and zeta7 != 1, "order-seven root")
    require(pow(zeta13, P, prime) == 1 and zeta13 != 1, "order-thirteen root")
    require(denominator % prime != 0, "normalization denominator")
    inverse_denominator = pow(denominator, -1, prime)

    support_shapes = []
    full_outputs = []
    full_spectra = []
    endpoint_interaction_outputs = []
    endpoint_component_supports = []
    endpoint_partition_unions = []
    for residue, matrix in enumerate(bank):
        components = pair_anova(matrix, prime)
        component_outputs = tuple(
            contract(ledger, component, inverse_denominator, prime) for component in components
        )
        full_output = contract(ledger, matrix, inverse_denominator, prime)
        require(output_add(component_outputs, prime) == full_output,
                ("contracted endpoint ANOVA", residue))

        full_spectrum = dft2(full_output, zeta7, zeta13, prime)
        component_spectra = tuple(
            dft2(output, zeta7, zeta13, prime) for output in component_outputs
        )
        component_coordinates = tuple(nonzero_coordinates(spectrum) for spectrum in component_spectra)
        union = frozenset().union(*component_coordinates)
        require(len(union) == Q * P and nonzero_coordinates(full_spectrum) == union,
                ("endpoint ANOVA full contraction", residue))
        full_coordinates = nonzero_coordinates(full_spectrum)
        shape = spectrum_shape(full_spectrum)
        require(shape == (91, 1, 6, 12, 72), (residue, shape))

        # Cell marginal returns the previous 13-owner pullback exactly.
        marginal_output = tuple(
            sum(
                marginal[left][right][owner] % prime * inverse_denominator % prime
                * matrix[left][right]
                for left in range(len(ATOMS)) for right in range(len(ATOMS))
            ) % prime
            for owner in range(P)
        )
        require(tuple(sum(full_output[cell][owner] for cell in range(Q)) % prime
                      for owner in range(P)) == marginal_output,
                ("cell erasure marginal", residue))

        output_mixed = output_anova(full_output, prime)[3]
        require(all(sum(output_mixed[cell][owner] for cell in range(Q)) % prime == 0
                    for owner in range(P)), ("cell erasure mixed", residue))
        require(all(sum(output_mixed[cell]) % prime == 0 for cell in range(Q)),
                ("owner erasure mixed", residue))

        support_shapes.append(shape)
        full_outputs.append(full_output)
        full_spectra.append(full_spectrum)
        endpoint_interaction_outputs.append(component_outputs[3])
        endpoint_component_supports.append(tuple(len(x) for x in component_coordinates))
        endpoint_partition_unions.append(len(union))

    fixed_output_mixed = output_anova(full_outputs[6], prime)[3]
    fixed_output_components = output_anova(fixed_output_mixed, prime)
    fixed_output_component_spectra = tuple(
        dft2(component, zeta7, zeta13, prime) for component in fixed_output_components
    )
    fixed_output_shape = (
        len(nonzero_coordinates(dft2(fixed_output_mixed, zeta7, zeta13, prime))),
    ) + tuple(len(nonzero_coordinates(spectrum)) for spectrum in fixed_output_component_spectra)
    require(fixed_output_shape == (72, 0, 0, 0, 72), fixed_output_shape)

    endpoint_totals = tuple(sum(sum(row) for row in matrix) % prime for matrix in bank)
    require(endpoint_totals[6] == 225010624370142818572, endpoint_totals[6])
    require(tuple(support_shapes) == ((91, 1, 6, 12, 72),) * P, support_shapes)
    require(tuple(endpoint_partition_unions) == (91,) * P, endpoint_partition_unions)
    require(endpoint_component_supports[6][3] == 91,
            ("fixed residue endpoint-pair interaction contraction", endpoint_component_supports[6]))

    endpoint_interaction_spectrum = dft2(
        endpoint_interaction_outputs[6], zeta7, zeta13, prime
    )
    endpoint_interaction_shape = spectrum_shape(endpoint_interaction_spectrum)
    require(endpoint_interaction_shape == (91, 1, 6, 12, 72),
            endpoint_interaction_shape)

    relation_output = full_outputs[6]
    inverse_cells = pow(Q, -1, prime)
    cell_erased = tuple(
        tuple(
            sum(relation_output[cell][owner] for cell in range(Q))
            % prime * inverse_cells % prime
            for owner in range(P)
        )
        for _cell in range(Q)
    )
    inverse_owners = pow(P, -1, prime)
    owner_erased = tuple(
        tuple(sum(relation_output[cell]) % prime * inverse_owners % prime for _owner in range(P))
        for cell in range(Q)
    )
    erasure_shapes = (
        spectrum_shape(dft2(cell_erased, zeta7, zeta13, prime)),
        spectrum_shape(dft2(owner_erased, zeta7, zeta13, prime)),
    )
    require(erasure_shapes == ((13, 1, 0, 12, 0), (7, 1, 6, 0, 0)), erasure_shapes)

    candidate_order_ledger = tuple(
        tuple(
            tuple(
                tuple(ledger[left][right][offset][cell] for offset in range(P))
                for cell in range(Q)
            )
            for right in range(len(ATOMS))
        )
        for left in range(len(ATOMS))
    )
    candidate_order_ledger_digest = digest(candidate_order_ledger)
    coordinate_bank_digest = digest(tuple(full_outputs))
    spectrum_bank_digest = digest(tuple(full_spectra))
    require(candidate_order_ledger_digest ==
            "39d7a0b4e5b2d8b85631d682ed1967091e44dc41e17b33a77e7184d3dc93e0cf",
            candidate_order_ledger_digest)
    require(coordinate_bank_digest ==
            "989dafc220a6d09aeacfce4af0e9a4fe13eedacc79fa66032ea39bc107fd8efb",
            coordinate_bank_digest)
    require(spectrum_bank_digest ==
            "5f173227c5e203309f61bdfd9d47cc64a3b49ae8f14abd0f7bfc469eda278533",
            spectrum_bank_digest)
    mixed_witness = full_spectra[6][1][1]
    require(mixed_witness == 218019411785559321795, mixed_witness)

    proof = (
        "seven_half_open_cells_partition_the_outer_base_before_marginalization",
        "cell_sum_recovers_every_39x39x13_source_entry",
        "complete_1_0_t_endpoint_bank_rebuilt_from_independent_THM3514_atoms",
        "full_output_has_DC_cell_owner_and_mixed_Fourier_sectors",
        "endpoint_pair_interaction_is_contracted_before_output_ANOVA",
        "either_cell_or_owner_erasure_annihilates_the_mixed_sector",
    )
    boundary = (
        "AX_BY_are_preintegrated_endpoint_scalars",
        "fixed_1_0_6_is_a_THM2334_residue_pushforward_not_exact_address_C_a_X_m",
        "not_a_physical_current_or_THM2512_bridge",
        "no_row_exclusion_or_LRC14_claim",
    )
    semantic_surface = (
        PRIMARY_SHA256, ATOM_SHA256, (prime, root, zeta7, zeta13, nn),
        tuple(cell_measures), digest(cells), digest(ledger), digest(marginal),
        atom_pair_support, entry_support, len(ATOMS) * len(ATOMS) * P * Q,
        tuple(nonempty_u), tuple(nonempty_v), digest(bank), endpoint_totals,
        tuple(support_shapes), tuple(endpoint_component_supports),
        tuple(endpoint_partition_unions), digest(tuple(full_outputs)),
        digest(tuple(full_spectra)), digest(tuple(endpoint_interaction_outputs)),
        fixed_output_shape, endpoint_interaction_shape, erasure_shapes,
        candidate_order_ledger_digest, coordinate_bank_digest,
        spectrum_bank_digest, mixed_witness, proof, boundary,
    )
    semantic_sha256 = hashlib.sha256(repr(semantic_surface).encode("utf-8")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
                (semantic_sha256, EXPECTED_SEMANTIC_SHA256))

    source_path = Path(__file__).resolve()
    print("R5 SOURCE-ALIGNED 7x13 RESIDUE SPECTRUM -- INDEPENDENT HOSTILE AUDIT")
    print("status=ACCEPT_SCOPED_FINITE_EXACT_PACKAGE;LRC14_OPEN")
    print(f"dependencies=(thm2594={PRIMARY_SHA256},thm3514_atom_engine={ATOM_SHA256})")
    print(f"field=(prime={prime},root={root},zeta7={zeta7},zeta13={zeta13},order={nn})")
    print(f"cell_partition=(count={len(cells)},measures={cell_measures},disjoint=True,tiles_base=True,digest={digest(cells)})")
    print(f"source_cell_ledger=(atom_pair_support={atom_pair_support}/1521,entry_support={entry_support}/138411,offset_cell_digest={digest(ledger)},cell_offset_digest={candidate_order_ledger_digest})")
    print(f"cell_marginal=(recovers_all_39x39x13_entries=True,digest={digest(marginal)})")
    print(f"endpoint_residue_bank=(digest={digest(bank)},totals={endpoint_totals})")
    print(f"residue_support_shapes={tuple(support_shapes)}")
    print(f"endpoint_pair_ANOVA=(component_supports={tuple(endpoint_component_supports)},partition_unions={tuple(endpoint_partition_unions)}/91)")
    print(f"fixed_relation_1_0_6=(endpoint_total={endpoint_totals[6]},mixed_witness_h1_k1={mixed_witness},output_ANOVA_shape={fixed_output_shape},endpoint_pair_ANOVA_shape={endpoint_interaction_shape})")
    print(f"erasures=(cell_shape={erasure_shapes[0]},owner_shape={erasure_shapes[1]},mixed_killed=True,cell_sum_recovers_prior_owner_pullback=True)")
    print(f"candidate_comparison=(coordinate_bank_digest={coordinate_bank_digest},spectrum_bank_digest={spectrum_bank_digest},all_exact=True)")
    print(f"proof={proof}")
    print(f"boundary={boundary}")
    print(f"semantic_sha256={semantic_sha256}")
    print(f"script_sha256_lf={lf_sha256(source_path)}")
    print("reproducibility=no_candidate_imports;no_randomness;no_elapsed_fields;normal_and_O_transcripts_must_match")
    print("commands=python -B 04-computation/lrc_r5_source_aligned_relation_residue_7x13_independent_audit_20260816.py;python -B -O 04-computation/lrc_r5_source_aligned_relation_residue_7x13_independent_audit_20260816.py")
    print("PASS")


if __name__ == "__main__":
    main()
