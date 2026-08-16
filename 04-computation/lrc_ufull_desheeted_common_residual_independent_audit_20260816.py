#!/usr/bin/env python3
"""Clean-room audit of the U_full desheeted common-residual package.

The submitted common-residual candidate is not imported.  Starting from the
hash-pinned THM-3514 endpoint interval engine and the conventions of the old
point-diagonal API, this program derives a residual coordinate, reconstructs
the old same-sheet point bank, and independently forms the ordered all-sheet
and cross-sheet residual couplings.

The endpoint OWNER condition descends exactly to the zeroth seven-cell.  Thus
the resulting 7x13 arrays are rank-one delta-cell lifts of one-dimensional
residue profiles, not genuine cell/residue mixing.  They remain alternate
finite diagonal couplings with no collision, horizon, source, chronology,
exact-address C(a;X,m), or LRC(14) meaning.
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
ENDPOINT_PATH = ROOT / (
    "04-computation/"
    "lrc_ufull_owner_boundary_k4xf13_endpoint_factorization_"
    "independent_audit_20260816.py"
)
POINT_PATH = ROOT / "04-computation/lrc_endpoint_ufull_atom_bridge_kernel_probe_20260816.py"
ENDPOINT_SHA256 = "f89be10c65bb77270199f9399b155d5a2c82c0da121b3e8589fe3c1f7e9824fc"
POINT_SHA256 = "d1182de4d777bab20a8d423cf942151ac3149014b67d9c34883cbce37a7b0a9f"
EXPECTED_SEMANTIC_SHA256 = "9d2070f27bcac8cf576a75bc542222be1e31679bc12e989c2f9bc276e8dd872c"

P = 13
Q = 7
EXPECTED_POINT_GAMMA_SHA256 = (
    "771545a5cb1f0f03459b8d351de668ad950ece5fcb985fa61d599d643de3303f"
)
EXPECTED_POINT_VALUES = (
    633668780131603861,
    405160484437854840264,
    167726070588785644466,
)
CARTESIAN_BRIDGE = 389266878372286537904
CONTROL_PAIRS = frozenset(((0, 0), (0, 1), (1, 0), (6, 6), (12, 0), (12, 12)))


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()


def digest_integers(values: tuple[int, ...]) -> str:
    return hashlib.sha256(
        ",".join(str(value) for value in values).encode("ascii")
    ).hexdigest()


def load_endpoint(name: str):
    require(lf_sha256(ENDPOINT_PATH) == ENDPOINT_SHA256, "endpoint engine hash drift")
    require(lf_sha256(POINT_PATH) == POINT_SHA256, "point API hash drift")
    spec = importlib.util.spec_from_file_location(name, ENDPOINT_PATH)
    require(spec is not None and spec.loader is not None, (name, "loader"))
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def merge_touch(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    intervals.sort()
    answer = [intervals[0]]
    for left, right in intervals[1:]:
        old_left, old_right = answer[-1]
        require(left >= old_right, ("overlapping intervals", answer[-1], (left, right)))
        if left == old_right:
            answer[-1] = (old_left, right)
        else:
            answer.append((left, right))
    return answer


def partition_intervals(
    intervals: list[tuple[int, int]],
    atoms: tuple[tuple[int, str, int, int], ...],
) -> tuple[tuple[tuple[int, int], ...], ...]:
    groups: list[list[tuple[int, int]]] = [[] for _atom in atoms]
    atom_index = 0
    input_measure = 0
    output_measure = 0
    for interval_left, interval_right in intervals:
        input_measure += interval_right - interval_left
        while atoms[atom_index][3] <= interval_left:
            atom_index += 1
        scan = atom_index
        while scan < len(atoms) and atoms[scan][2] < interval_right:
            _sheet, _chamber, atom_left, atom_right = atoms[scan]
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


def guard_safe(sheet: int, chamber: str, tau: int) -> int:
    midpoint_numerator = {"left": 1, "middle": 7, "right": 13}[chamber]
    doubled_s = (14 * (sheet + tau) + midpoint_numerator) % 182
    distance = min(doubled_s, 182 - doubled_s)
    return int(distance > 26)


def seven_cells(denominator: int) -> tuple[tuple[tuple[int, int], ...], ...]:
    require(denominator % 14 == 0, denominator)
    half = denominator // 14
    width = denominator // Q
    cells = []
    for cell in range(Q):
        left = (cell * width - half) % denominator
        right = left + width
        if right <= denominator:
            cells.append(((left, right),))
        else:
            cells.append(((0, right - denominator), (left, denominator)))
    flattened = sorted((left, right, cell) for cell, pieces in enumerate(cells) for left, right in pieces)
    require(flattened[0][0] == 0 and flattened[-1][1] == denominator, "cell endpoints")
    require(all(old[1] == new[0] for old, new in zip(flattened, flattened[1:])),
            "cell partition")
    return tuple(cells)


def literal_guard_control(A, alpha: int, beta: int, groups, atoms, t_den: int) -> tuple[int, ...]:
    records = []
    for tau in range(P):
        selected = [
            interval
            for index, group in enumerate(groups)
            if guard_safe(atoms[index][0], atoms[index][1], tau)
            for interval in group
        ]
        selected_merged = merge_touch(selected)
        ell = (tau, -alpha % P, -beta % P, 0, 0, 0, 0, alpha, beta)
        literal = A.M.fast_build_set(A.context()[0], t_den, A.M.PATTERN_E, ell)
        literal_merged = merge_touch(list(literal))
        require(selected_merged == literal_merged,
                ("literal guard restoration", alpha, beta, tau))
        records.append(sum(right - left for left, right in literal_merged))
    return tuple(records)


def residual_worker(alpha: int) -> tuple[object, ...]:
    A = load_endpoint(f"thm3514_desheeted_residual_worker_{alpha}")
    (
        word, t_den, nn, prime, root, zeta13, q_intervals, _q_starts,
        _embeddings, _tabs, atom_intervals,
    ) = A.context()
    require(tuple(A.ATOMS) == tuple((sheet, chamber) for sheet in range(P)
                                    for chamber in A.CHAMBER_NAMES), "atom order")
    sheet_width = t_den // P
    residual_denominator = P * t_den
    cells = seven_cells(residual_denominator)
    flat_cells = sorted(
        (left, right, cell) for cell, pieces in enumerate(cells) for left, right in pieces
    )
    cell_starts = [left for left, _right, _cell in flat_cells]

    # u=169(x-a*T/13), so y=u/(13T).  Q(169t) is Q(u mod T).
    q_lift = tuple(
        (left + branch * t_den, right + branch * t_den)
        for branch in range(P) for left, right in q_intervals
    )
    q_starts = tuple(left for left, _right in q_lift)
    require(q_lift[0][0] >= 0 and q_lift[-1][1] <= residual_denominator,
            "Q residual lift")

    target_frequency = word[A.M.TARGET_B]
    require(target_frequency == 742586 == P**2 * 4394, target_frequency)
    require(target_frequency // P == P * 4394, "frequency descent")
    require((target_frequency * t_den) % nn == 0, "residual phase periodicity")

    no_guard = dict(A.M.PATTERN_E)
    require(no_guard.pop(A.M.GUARD) == "guard_safe", "guard deletion")
    safe_masks = tuple(
        sum(
            (guard_safe(sheet, chamber, tau) << index)
            for index, (sheet, chamber, _left, _right) in enumerate(atom_intervals)
        )
        for tau in range(P)
    )
    require(
        all(
            guard_safe(sheet, chamber, tau) == A.safe(chamber, sheet + tau)
            for tau in range(P)
            for sheet, chamber, _left, _right in atom_intervals
        ),
        "derived guard table",
    )

    phase_cache: dict[int, int] = {}

    def phase(point: int) -> int:
        key = point % t_den
        value = phase_cache.get(key)
        if value is None:
            value = pow(root, (target_frequency * key) % nn, prime)
            phase_cache[key] = value
        return value

    rows = []
    control_records = []
    owner_cell_zero_exact = True
    for beta in range(P):
        ell = (0, -alpha % P, -beta % P, 0, 0, 0, 0, alpha, beta)
        unguarded = A.M.fast_build_set(word, t_den, no_guard, ell)
        groups = partition_intervals(unguarded, atom_intervals)
        require(all(not groups[index] for index, (_sheet, chamber, _left, _right)
                    in enumerate(atom_intervals) if chamber == "middle"),
                ("middle owner leak", alpha, beta))

        if (alpha, beta) in CONTROL_PAIRS:
            control_records.append((beta, literal_guard_control(
                A, alpha, beta, groups, atom_intervals, t_den
            )))

        events: dict[int, int] = {0: 0, residual_denominator: 0}
        for _left, right, _cell in flat_cells:
            events.setdefault(right, 0)
        endpoint_identity_checks = 0
        for index, group in enumerate(groups):
            sheet, _chamber, atom_left, atom_right = atom_intervals[index]
            require(atom_left // sheet_width == sheet or atom_left == sheet * sheet_width,
                    ("atom sheet", index))
            bit = 1 << index
            sheet_start = sheet * sheet_width
            for left, right in group:
                require(sheet_start <= left < right <= sheet_start + sheet_width,
                        ("sheet interval", index, left, right))
                residual_left = P**2 * (left - sheet_start)
                residual_right = P**2 * (right - sheet_start)
                require(0 <= residual_left < residual_right <= residual_denominator,
                        ("residual interval", index))
                # 169*x = 13*a*T + u, hence Q(169t)=Q(u mod T).
                require((P**2 * left - residual_left) % t_den == 0,
                        ("Q descent", index, left))
                require((P**2 * right - residual_right) % t_den == 0,
                        ("Q descent", index, right))
                # target_frequency is a multiple of 169, so the sheet phase vanishes.
                require(
                    pow(root, (target_frequency * P**2 * left) % nn, prime)
                    == phase(residual_left),
                    ("frequency descent", index, left),
                )
                events[residual_left] = events.get(residual_left, 0) ^ bit
                events[residual_right] = events.get(residual_right, 0) ^ bit
                endpoint_identity_checks += 2

        points = sorted(events)
        active = 0
        same = [[0] * Q for _tau in range(P)]
        full = [[0] * Q for _tau in range(P)]
        count_cache: dict[int, tuple[int, ...]] = {}
        for position, left in enumerate(points[:-1]):
            active ^= events[left]
            right = points[position + 1]
            if not active or left == right:
                continue
            cell_position = bisect_right(cell_starts, left) - 1
            cell_left, cell_right, cell = flat_cells[cell_position]
            require(cell_left <= left < right <= cell_right,
                    ("cell split", left, right, flat_cells[cell_position]))
            # OWNER has speed 13 and remains inside the unguarded endpoint
            # pattern.  Under t=(y+a)/13 it becomes
            #
            #     ||13t|| = ||y+a|| = ||y|| < 1/14,
            #
            # exactly the defining interval of cell_0.  This is an interval-
            # geometry statement over Q, before finite-field integration.
            require(cell == 0,
                    ("endpoint OWNER escaped cell zero", alpha, beta, left, right, cell))
            counts = count_cache.get(active)
            if counts is None:
                counts = tuple((active & safe_mask).bit_count() for safe_mask in safe_masks)
                count_cache[active] = counts

            q_index = bisect_right(q_starts, left) - 1
            if q_index < 0:
                q_index = 0
            while q_index < len(q_lift) and q_lift[q_index][1] <= left:
                q_index += 1
            while q_index < len(q_lift) and q_lift[q_index][0] < right:
                q_left, q_right = q_lift[q_index]
                overlap_left = max(left, q_left)
                overlap_right = min(right, q_right)
                if overlap_left < overlap_right:
                    jump = (phase(overlap_left) - phase(overlap_right)) % prime
                    for tau, count in enumerate(counts):
                        same[tau][cell] = (same[tau][cell] + count * jump) % prime
                        full[tau][cell] = (full[tau][cell] + count * count * jump) % prime
                q_index += 1

        character_phase = pow(zeta13, beta, prime)
        for tau in range(P):
            same_row = tuple(character_phase * value % prime for value in same[tau])
            full_row = tuple(character_phase * value % prime for value in full[tau])
            cross_row = tuple((all_value - diagonal) % prime
                              for all_value, diagonal in zip(full_row, same_row))
            rows.append((beta, tau, same_row, cross_row, full_row))
        require(endpoint_identity_checks > 0, ("empty endpoint checks", alpha, beta))

    return alpha, tuple(rows), tuple(control_records), owner_cell_zero_exact


def inverse_tables(gamma_cells, zeta13: int, prime: int):
    normalizer = pow(P**3, -1, prime)
    tables = []
    for kind in range(3):
        table = []
        for cell in range(Q):
            row = []
            for residue in range(P):
                total = 0
                index = 0
                for alpha in range(P):
                    for beta in range(P):
                        for tau in range(P):
                            exponent = -(alpha + tau * residue) % P
                            total = (
                                total + gamma_cells[index][kind][cell]
                                * pow(zeta13, exponent, prime)
                            ) % prime
                            index += 1
                row.append(total * normalizer % prime)
            table.append(tuple(row))
        tables.append(tuple(table))
    return tuple(tables)


def dft2(table, zeta7: int, zeta13: int, prime: int):
    return tuple(
        tuple(
            sum(
                table[cell][residue]
                * pow(zeta7, (-cell_frequency * cell) % Q, prime)
                * pow(zeta13, (-residue_frequency * residue) % P, prime)
                for cell in range(Q) for residue in range(P)
            ) % prime
            for residue_frequency in range(P)
        )
        for cell_frequency in range(Q)
    )


def support_shape(spectrum) -> tuple[int, int, int, int, int]:
    dc = int(spectrum[0][0] != 0)
    cell_axis = sum(spectrum[frequency][0] != 0 for frequency in range(1, Q))
    residue_axis = sum(spectrum[0][frequency] != 0 for frequency in range(1, P))
    mixed = sum(
        spectrum[cell_frequency][residue_frequency] != 0
        for cell_frequency in range(1, Q) for residue_frequency in range(1, P)
    )
    return dc + cell_axis + residue_axis + mixed, dc, cell_axis, residue_axis, mixed


def output_interaction(table, prime: int):
    inverse_cells = pow(Q, -1, prime)
    inverse_residues = pow(P, -1, prime)
    row_means = tuple(sum(row) % prime * inverse_residues % prime for row in table)
    column_means = tuple(
        sum(table[cell][residue] for cell in range(Q)) % prime * inverse_cells % prime
        for residue in range(P)
    )
    grand = sum(sum(row) for row in table) % prime
    grand = grand * inverse_cells % prime * inverse_residues % prime
    answer = tuple(
        tuple((table[cell][residue] - row_means[cell] - column_means[residue] + grand) % prime
              for residue in range(P))
        for cell in range(Q)
    )
    require(all(sum(row) % prime == 0 for row in answer), "interaction row sums")
    require(all(sum(answer[cell][residue] for cell in range(Q)) % prime == 0
                for residue in range(P)), "interaction column sums")
    return answer


def rank_mod(matrix, prime: int) -> int:
    rows = [[value % prime for value in row] for row in matrix]
    columns = len(rows[0]) if rows else 0
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
            rows[row] = [
                (value - factor * pivot_value) % prime
                for value, pivot_value in zip(rows[row], rows[rank])
            ]
        rank += 1
        if rank == len(rows):
            break
    return rank


def main() -> None:
    with ProcessPoolExecutor(max_workers=4) as pool:
        chunks = tuple(pool.map(residual_worker, range(P)))
    require(tuple(chunk[0] for chunk in chunks) == tuple(range(P)), "worker order")
    require(all(chunk[3] for chunk in chunks), "owner support escaped cell zero")

    A = load_endpoint("thm3514_desheeted_residual_main")
    word, t_den, nn, prime, root, zeta13, *_rest = A.context()
    require(word[A.M.OWNER] == P, ("OWNER speed", word[A.M.OWNER]))
    zeta7 = pow(root, nn // Q, prime)
    require(pow(zeta7, Q, prime) == 1 and zeta7 != 1, "order-seven root")
    require(pow(zeta13, P, prime) == 1 and zeta13 != 1, "order-thirteen root")

    rows = tuple(row for chunk in chunks for row in chunk[1])
    require(len(rows) == P**3, len(rows))
    gamma_cells = tuple((row[2], row[3], row[4]) for row in rows)
    require(
        all(value == 0 for row in gamma_cells for sector in row for value in sector[1:]),
        "character bank escaped exact owner cell",
    )
    same_gamma = tuple(sum(row[0]) % prime for row in gamma_cells)
    cross_gamma = tuple(sum(row[1]) % prime for row in gamma_cells)
    full_gamma = tuple(sum(row[2]) % prime for row in gamma_cells)
    require(all((same + cross) % prime == full
                for same, cross, full in zip(same_gamma, cross_gamma, full_gamma)),
            "same plus cross gamma")
    point_hash = digest_integers(same_gamma)
    require(point_hash == EXPECTED_POINT_GAMMA_SHA256, point_hash)
    cross_gamma_hash = digest_integers(cross_gamma)
    full_gamma_hash = digest_integers(full_gamma)
    require(cross_gamma_hash ==
            "f02e10902bd8cc5ec9b588f059c7ebb72eac14b5f9ee46f6055f2f2ee0a530c6",
            cross_gamma_hash)
    require(full_gamma_hash ==
            "c16e6c4e82b73aced92f5d9ab61779d535263deed17a9caa90fce0ca6c947dcf",
            full_gamma_hash)

    same_table, cross_table, full_table = inverse_tables(gamma_cells, zeta13, prime)
    tables = (same_table, cross_table, full_table)
    require(
        all((same_table[cell][residue] + cross_table[cell][residue]) % prime
            == full_table[cell][residue]
            for cell in range(Q) for residue in range(P)),
        "same plus cross tables",
    )
    require(
        all(value == 0 for table in tables for row in table[1:] for value in row),
        "inverse table escaped exact owner cell",
    )
    spectra = tuple(dft2(table, zeta7, zeta13, prime) for table in tables)
    table_hashes = tuple(digest(table) for table in tables)
    spectrum_hashes = tuple(digest(spectrum) for spectrum in spectra)
    require(table_hashes == (
        "b80423dad8b155bd02e79d98c981ec3ada15934d91cb72ebacc385ded6b7ed6d",
        "50d44c5db5cea5842ceada123d06c54f0221561d08744c83eaf3aebd55f2420c",
        "64adfc3850b8a8924546a75e051116aa10e3a62661ca409175aad174a9ebb926",
    ), table_hashes)
    require(spectrum_hashes == (
        "17e08d5c0a5d107bad1d9828e754ed197d2737ee4dc7c1fcbfdeccb410992f72",
        "14c5490ee52e4283848290feee3b0995a44791fee30a9aaef6ad2a945b9442e8",
        "398a2c2ccb57714aa9d735d3bca4f24c2212396a8be9e2bd07e8e4ec57b6c05a",
    ), spectrum_hashes)
    shapes = tuple(support_shape(spectrum) for spectrum in spectra)
    require(shapes == ((91, 1, 6, 12, 72),) * 3, shapes)
    interactions = tuple(output_interaction(table, prime) for table in tables)
    interaction_shapes = tuple(
        support_shape(dft2(interaction, zeta7, zeta13, prime))
        for interaction in interactions
    )
    require(interaction_shapes == ((72, 0, 0, 0, 72),) * 3, interaction_shapes)

    coordinate_ranks = tuple(rank_mod(table, prime) for table in tables)
    interaction_ranks = tuple(rank_mod(interaction, prime) for interaction in interactions)
    require(coordinate_ranks == (1, 1, 1), coordinate_ranks)
    require(interaction_ranks == (1, 1, 1), interaction_ranks)
    residue_profiles = tuple(table[0] for table in tables)
    residue_spectra = tuple(
        tuple(
            sum(
                profile[residue] * pow(zeta13, (-frequency * residue) % P, prime)
                for residue in range(P)
            ) % prime
            for frequency in range(P)
        )
        for profile in residue_profiles
    )
    require(all(all(value != 0 for value in spectrum) for spectrum in residue_spectra),
            residue_spectra)
    require(
        all(
            spectra[sector][cell_frequency][residue_frequency]
            == residue_spectra[sector][residue_frequency]
            for sector in range(3)
            for cell_frequency in range(Q)
            for residue_frequency in range(P)
        ),
        "delta-cell Fourier factorization",
    )
    inverse_cells = pow(Q, -1, prime)
    inverse_residues = pow(P, -1, prime)
    for sector, (profile, interaction) in enumerate(zip(residue_profiles, interactions)):
        residue_mean = sum(profile) % prime * inverse_residues % prime
        require(
            all(
                interaction[cell][residue]
                == (
                    (int(cell == 0) - inverse_cells)
                    * (profile[residue] - residue_mean)
                ) % prime
                for cell in range(Q)
                for residue in range(P)
            ),
            ("ANOVA outer product", sector),
        )
    residue_hashes = tuple(digest(profile) for profile in residue_profiles)
    residue_spectrum_hashes = tuple(digest(spectrum) for spectrum in residue_spectra)

    fixed_modes = tuple(
        tuple(
            sum(table[cell][6] * pow(zeta7, (-frequency * cell) % Q, prime)
                for cell in range(Q)) % prime
            for frequency in range(Q)
        )
        for table in tables
    )
    require(all(all(value != 0 for value in modes) for modes in fixed_modes), fixed_modes)
    exact_zero_cell_rows = tuple(
        tuple(cell for cell in range(1, Q) if all(table[cell][residue] == 0
                                                 for residue in range(P)))
        for table in tables
    )
    require(exact_zero_cell_rows == ((1, 2, 3, 4, 5, 6),) * 3,
            exact_zero_cell_rows)

    same_values = (
        sum(same_table[cell][1] for cell in range(Q)) % prime,
        sum(same_table[cell][0] for cell in range(Q)) % prime,
    )
    same_bridge = (same_values[0] - same_values[1]) % prime
    require((*same_values, same_bridge) == EXPECTED_POINT_VALUES,
            (*same_values, same_bridge))
    full_values = (
        sum(full_table[cell][1] for cell in range(Q)) % prime,
        sum(full_table[cell][0] for cell in range(Q)) % prime,
    )
    full_bridge = (full_values[0] - full_values[1]) % prime
    cross_bridge = (full_bridge - same_bridge) % prime
    require(full_values == (95783417057771114126, 124341028951542154618), full_values)
    require(full_bridge == 543695274352737840377, full_bridge)
    require(cross_bridge == 375969203763952195911, cross_bridge)
    require(full_bridge != CARTESIAN_BRIDGE, (full_bridge, CARTESIAN_BRIDGE))

    controls = tuple((chunk[0], chunk[2]) for chunk in chunks if chunk[2])
    proof = (
        "91t_equals_7a_plus_r_and_y_equals_13t_minus_a",
        "Q_169t_equals_Q_13y_on_every_transformed_endpoint",
        "742586_equals_13_squared_times_4394_kills_the_sheet_phase",
        "literal_guarded_controls_equal_unguarded_atom_restoration",
        "same_sheet_bank_recovers_the_old_point_diagonal_hash",
        "ordered_same_plus_cross_equals_full_before_every_transform",
        "OWNER_norm_13t_equals_norm_y_forces_exact_cell_0_support",
        "inverse_tables_are_delta_cell0_times_one_dimensional_residue_profiles",
        "ANOVA_is_a_rank_one_outer_product_not_genuine_cell_residue_mixing",
    )
    boundary = (
        "alternate_common_residual_diagonal_coupling",
        "full_bridge_differs_from_the_Cartesian_endpoint_bridge",
        "formal_7x13_Fourier_support_is_separable_and_not_a_mixing_certificate",
        "no_collision_horizon_source_chronology_exact_address_or_LRC14_claim",
    )
    semantic_surface = (
        ENDPOINT_SHA256, POINT_SHA256, (word, t_den, nn, prime, root, zeta7, zeta13),
        digest(controls), point_hash, cross_gamma_hash, full_gamma_hash,
        table_hashes, spectrum_hashes, shapes, interaction_shapes, coordinate_ranks,
        interaction_ranks, residue_hashes, residue_spectrum_hashes, fixed_modes,
        exact_zero_cell_rows, (*same_values, same_bridge), (*full_values, full_bridge),
        cross_bridge, proof, boundary,
    )
    semantic_sha256 = hashlib.sha256(repr(semantic_surface).encode("utf-8")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
                (semantic_sha256, EXPECTED_SEMANTIC_SHA256))

    source_path = Path(__file__).resolve()
    print("U_FULL DESHEETED COMMON-RESIDUAL -- INDEPENDENT HOSTILE AUDIT")
    print("status=ACCEPT_SCOPED_FINITE_EXACT_ALTERNATE_DIAGONAL;LRC14_OPEN")
    print(f"dependencies=(thm3514_engine={ENDPOINT_SHA256},old_point_api={POINT_SHA256})")
    print(f"field=(prime={prime},root={root},zeta7={zeta7},zeta13={zeta13},order={nn})")
    print(f"coordinate_descent=(91t=7a+r,y=13t-a,Q169t=Q13y,742586={P**2}*4394,residual_frequency={word[A.M.TARGET_B]//P})")
    print(f"guard_restoration=(atom_law_exhaustive=True,literal_controls={len(CONTROL_PAIRS)*P},digest={digest(controls)})")
    print(f"same_sheet_point_gamma=(hash={point_hash},values={same_values},bridge={same_bridge})")
    print(f"gamma_decomposition=(same_plus_cross_equals_full=True,same={point_hash},cross={cross_gamma_hash},full={full_gamma_hash})")
    print(f"table_spectra=(order=(same,cross,full),shapes={shapes},table_digests={table_hashes},spectrum_digests={spectrum_hashes})")
    print(f"cell_geometry=(OWNER_implies_exact_cell_0,character_occupancy_by_cell={(P**3, 0, 0, 0, 0, 0, 0)},exact_zero_rows={exact_zero_cell_rows})")
    print(f"rank_one_factorization=(table_ranks={coordinate_ranks},ANOVA_ranks={interaction_ranks},table=delta_cell0_times_R,ANOVA=(delta_cell0-1/7)_times_(R-mean_R))")
    print(f"residue_profiles=(profile_digests={residue_hashes},spectrum_digests={residue_spectrum_hashes},all_13_modes_nonzero=True)")
    print(f"formal_spectral_shapes_not_mixing=(raw={shapes},ANOVA={interaction_shapes})")
    print(f"fixed_class_1_0_6_delta_profiles_with_repeated_F7_modes={fixed_modes}")
    print(f"bridge_boundary=(full_values={full_values},full_bridge={full_bridge},cross_bridge={cross_bridge},Cartesian_bridge={CARTESIAN_BRIDGE},different=True)")
    print(f"proof={proof}")
    print(f"boundary={boundary}")
    print(f"semantic_sha256={semantic_sha256}")
    print(f"script_sha256_lf={lf_sha256(source_path)}")
    print("reproducibility=no_candidate_imports;no_randomness;no_elapsed_fields;normal_and_O_transcripts_must_match")
    print("commands=python -B 04-computation/lrc_ufull_desheeted_common_residual_independent_audit_20260816.py;python -B -O 04-computation/lrc_ufull_desheeted_common_residual_independent_audit_20260816.py")
    print("PASS")


if __name__ == "__main__":
    main()
