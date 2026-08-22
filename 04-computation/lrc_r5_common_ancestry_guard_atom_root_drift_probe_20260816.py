#!/usr/bin/env python3
"""Exact common-ancestry guard-atom/root-drift probe at the r=5 window.

THM-2471 equation (31) permits an arbitrary rational Boolean partition on
the linked arrival and source nodes X_(u,a), Y_(q,e').  Here that partition
is the actual 39-atom U_full guard partition I_(sheet,chamber) from THM-3514.
This is a lawful Boolean ancestry refinement; it never identifies X and Y as
one circle point and never replaces their common outer base by an endpoint
diagonal.

The output retains two independently typed drifts:

    collision-root drift s = u-q,
    guard-address drift  d = right_sheet-left_sheet.

The common torsor gauge is left_sheet-u = right_sheet-q, equivalently
s=-d.  We test this predicate before every marginalization, audit its full
3x3x13 bucket support, and then inspect the owner-active 4x13 K4 carrier.
No endpoint weight identity, physical current, row exclusion, or LRC(14)
conclusion is asserted.
"""

from __future__ import annotations

from bisect import bisect_right
from fractions import Fraction
from hashlib import sha256
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
P = 13
STAGE_PATH = ROOT / "04-computation/lrc14_stage2_theta_contraction_opus_20260728.py"
STAGE_SHA256 = "09c43af0a0a56c7a0833bbfd13ed6a96bc5a7a3718aa1bc6b77a144bde101a06"
EXPECTED_SEMANTIC_SHA256 = "3d8c88fb7b9762f41ef35c00d980b99fc435c8352baf5dddb9fe412d1baeace0"

CHAMBERS = ("left", "middle", "right")
ACTIVE = ("left", "right")
CHAMBER_INDEX = {name: index for index, name in enumerate(CHAMBERS)}
WALSH_SIGNS = (
    (1, 1, 1, 1),
    (1, 1, -1, -1),
    (1, -1, 1, -1),
    (1, -1, -1, 1),
)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_json(value: object) -> str:
    body = json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    return sha256(body).hexdigest()


def load_stage():
    require(lf_sha256(STAGE_PATH) == STAGE_SHA256, "THM-2594 stage source drift")
    spec = importlib.util.spec_from_file_location("thm2594_guard_atom_probe", STAGE_PATH)
    require(spec is not None and spec.loader is not None, "stage module loader")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


M = load_stage()
GRID = M.T_DEN
require(GRID % 91 == 0, "stage grid does not resolve guard atoms")
ATOM_UNIT = GRID // 91
ATOMS = tuple((sheet, chamber) for sheet in range(P) for chamber in CHAMBERS)
ATOM_INDEX = {atom: index for index, atom in enumerate(ATOMS)}


def atom_ranges() -> tuple[tuple[int, int], ...]:
    ranges = []
    chamber_bounds = ((0, 1), (1, 6), (6, 7))
    for sheet in range(P):
        for low, high in chamber_bounds:
            ranges.append(((7 * sheet + low) * ATOM_UNIT, (7 * sheet + high) * ATOM_UNIT))
    require(ranges[0][0] == 0 and ranges[-1][1] == GRID, "atom boundary endpoints")
    require(all(left[1] == right[0] for left, right in zip(ranges, ranges[1:])),
            "guard atoms do not partition the circle")
    return tuple(ranges)


RANGES = atom_ranges()
RANGE_STARTS = tuple(left for left, _right in RANGES)


def partition_weighted_pieces(
    pieces: list[tuple[int, int, int]],
) -> tuple[list[tuple[int, int, int]], ...]:
    """Intersect a disjoint weighted profile with the 39 guard atoms."""
    groups: list[list[tuple[int, int, int]]] = [[] for _atom in ATOMS]
    input_mass = 0
    output_mass = 0
    for left, right, weight in pieces:
        require(0 <= left < right <= GRID and weight >= 0, (left, right, weight))
        input_mass += weight * (right - left)
        atom_index = bisect_right(RANGE_STARTS, left) - 1
        require(atom_index >= 0, (left, atom_index))
        cursor = left
        while cursor < right:
            atom_right = RANGES[atom_index][1]
            stop = min(right, atom_right)
            if weight:
                groups[atom_index].append((cursor, stop, weight))
                output_mass += weight * (stop - cursor)
            cursor = stop
            if cursor < right:
                atom_index += 1
    require(input_mass == output_mass, ("partition mass", input_mass, output_mass))
    return tuple(groups)


def profile_pieces(starts: list[int], values: list[int]) -> list[tuple[int, int, int]]:
    answer = []
    for index, left in enumerate(starts):
        right = starts[index + 1] if index + 1 < len(starts) else GRID
        if values[index]:
            answer.append((left, right, values[index]))
    return answer


def profile_value(starts: list[int], values: list[int], point: int) -> int:
    return values[bisect_right(starts, point) - 1]


def require_profile_partition(
    parts: tuple[tuple[list[int], list[int]], ...],
    whole: tuple[list[int], list[int]],
    label: str,
) -> None:
    breakpoints = {0}
    breakpoints.update(whole[0])
    for starts, _values in parts:
        breakpoints.update(starts)
    for point in sorted(breakpoints):
        if point >= GRID:
            continue
        total = sum(profile_value(starts, values, point) for starts, values in parts)
        expected = profile_value(whole[0], whole[1], point)
        require(total == expected, (label, point, total, expected))


def build_f_pieces(
    e_intervals: list[tuple[int, int]],
    q_intervals: list[tuple[int, int]],
) -> list[tuple[int, int, int]]:
    n_starts, n_values = M.weighted_fold(
        [(left, right, 1) for left, right in e_intervals], M.RPKT, GRID
    )
    pieces: list[tuple[int, int, int]] = []
    for q_left, q_right in q_intervals:
        index = bisect_right(n_starts, q_left) - 1
        while True:
            p_left = n_starts[index]
            p_right = n_starts[index + 1] if index + 1 < len(n_starts) else GRID
            left = max(q_left, p_left)
            right = min(q_right, p_right)
            if left < right and n_values[index]:
                pieces.append((left, right, n_values[index]))
            if p_right >= q_right:
                break
            index += 1
    return pieces


def rank_mod(rows: list[list[int]], prime: int) -> int:
    matrix = [[value % prime for value in row] for row in rows]
    rank = 0
    columns = len(matrix[0]) if matrix else 0
    for column in range(columns):
        pivot = next((row for row in range(rank, len(matrix)) if matrix[row][column]), None)
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse = pow(matrix[rank][column], -1, prime)
        matrix[rank] = [value * inverse % prime for value in matrix[rank]]
        for row in range(len(matrix)):
            if row == rank or not matrix[row][column]:
                continue
            factor = matrix[row][column]
            matrix[row] = [
                (value - factor * pivot_value) % prime
                for value, pivot_value in zip(matrix[row], matrix[rank])
            ]
        rank += 1
        if rank == len(matrix):
            break
    return rank


def root_of_order_13(prime: int) -> int:
    for value in range(2, prime):
        if pow(value, P, prime) == 1 and value != 1:
            return value
    raise RuntimeError(("no order-13 root", prime))


def fourier_support(rows: list[list[int]], prime: int) -> tuple[int, ...]:
    root = root_of_order_13(prime)
    return tuple(
        sum(
            value * pow(root, -frequency * drift % P, prime)
            for drift, value in enumerate(row)
        ) % prime
        for row in rows
        for frequency in range(P)
    )


def main() -> None:
    e_intervals = M.build_set(M.PAT_E, M.ZELL)
    q_intervals = M.build_set(M.PAT_QB, M.ZELL)
    f_pieces = build_f_pieces(e_intervals, q_intervals)

    e_groups = partition_weighted_pieces(
        [(left, right, 1) for left, right in e_intervals]
    )
    f_groups = partition_weighted_pieces(f_pieces)
    nonempty_e_groups = sum(bool(group) for group in e_groups)
    nonempty_f_groups = sum(bool(group) for group in f_groups)
    require(nonempty_e_groups > 0 and nonempty_f_groups > 0,
            "both ancestry legs must retain at least one guard atom")

    u_profiles = tuple(M.weighted_fold(list(group), M.DCOLL, GRID) for group in f_groups)
    v_profiles = tuple(M.weighted_fold(list(group), M.DCOLL, GRID) for group in e_groups)
    u_whole = M.weighted_fold(f_pieces, M.DCOLL, GRID)
    v_whole = M.weighted_fold(
        [(left, right, 1) for left, right in e_intervals], M.DCOLL, GRID
    )
    require_profile_partition(u_profiles, u_whole, "arrival atom partition")
    require_profile_partition(v_profiles, v_whole, "source atom partition")

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

    # bucket[C][D][address drift d][root drift s]
    bucket = [[[[0 for _s in range(P)] for _d in range(P)] for _right in CHAMBERS]
              for _left in CHAMBERS]
    # On the common-gauge slice s=-d, retain the shared absolute torsor
    # offset c=left_sheet-current_root=right_sheet-source_root.
    gauge_offset = [[[[0 for _offset in range(P)] for _d in range(P)]
                     for _right in CHAMBERS] for _left in CHAMBERS]
    pair_gauge_offset = [[[0 for _offset in range(P)] for _right in ATOMS]
                         for _left in ATOMS]
    atom_pair_nonzero = 0
    for left_index, (left_sheet, left_chamber) in enumerate(ATOMS):
        ci = CHAMBER_INDEX[left_chamber]
        for right_index, (right_sheet, right_chamber) in enumerate(ATOMS):
            cj = CHAMBER_INDEX[right_chamber]
            drift = (right_sheet - left_sheet) % P
            pair_nonzero = False
            for root_drift in range(P):
                mass = 0
                for current_root in range(P):
                    source_root = (current_root - root_drift) % P
                    root_mass = M.product_cum(
                        u_windows[left_index][current_root][0],
                        u_windows[left_index][current_root][1],
                        v_windows[right_index][source_root][0],
                        v_windows[right_index][source_root][1],
                        GRID,
                    )[3]
                    mass += root_mass
                    if root_drift == (-drift) % P:
                        common_offset = (left_sheet - current_root) % P
                        require(
                            common_offset == (right_sheet - source_root) % P,
                            "common-gauge offset mismatch",
                        )
                        gauge_offset[ci][cj][drift][common_offset] += root_mass
                        pair_gauge_offset[left_index][right_index][common_offset] += root_mass
                bucket[ci][cj][drift][root_drift] += mass
                pair_nonzero = pair_nonzero or mass != 0
            atom_pair_nonzero += int(pair_nonzero)

    root_marginal = []
    for root_drift in range(P):
        mass = 0
        for current_root in range(P):
            source_root = (current_root - root_drift) % P
            mass += M.product_cum(
                u_whole_windows[current_root][0],
                u_whole_windows[current_root][1],
                v_whole_windows[source_root][0],
                v_whole_windows[source_root][1],
                GRID,
            )[3]
        root_marginal.append(mass)
        reconstructed = sum(
            bucket[left][right][drift][root_drift]
            for left in range(3)
            for right in range(3)
            for drift in range(P)
        )
        require(reconstructed == mass, ("root marginal", root_drift, reconstructed, mass))

    require(root_marginal[0] == 0, "first-collision diagonal did not vanish")
    require(all(value > 0 for value in root_marginal[1:]), "a nonzero root colour vanished")

    denominator = M.RPKT * M.DCOLL * M.DCOLL * GRID
    total_mass = sum(root_marginal)
    require(Fraction(total_mass, denominator) == 169 * M.I5,
            ("THM-2471 service anchor", Fraction(total_mass, denominator), 169 * M.I5))

    gauge_rows = []
    opposite_rows = []
    active_gauge_rows = []
    for left_chamber in CHAMBERS:
        for right_chamber in CHAMBERS:
            ci = CHAMBER_INDEX[left_chamber]
            cj = CHAMBER_INDEX[right_chamber]
            gauge = tuple(bucket[ci][cj][drift][(-drift) % P] for drift in range(P))
            opposite = tuple(bucket[ci][cj][drift][drift] for drift in range(P))
            gauge_rows.append((left_chamber, right_chamber, gauge))
            opposite_rows.append((left_chamber, right_chamber, opposite))
            require(gauge[0] == opposite[0] == 0, (left_chamber, right_chamber, "d=0"))
            if left_chamber in ACTIVE and right_chamber in ACTIVE:
                active_gauge_rows.append(list(gauge))
            require(
                all(sum(gauge_offset[ci][cj][drift]) == gauge[drift] for drift in range(P)),
                (left_chamber, right_chamber, "offset marginal"),
            )

    require(len(active_gauge_rows) == 4, "K4 corner row count")
    active_gauge_support = sum(value != 0 for row in active_gauge_rows for value in row)
    full_gauge_support = sum(value != 0 for _left, _right, row in gauge_rows for value in row)
    opposite_support = sum(value != 0 for _left, _right, row in opposite_rows for value in row)

    primes = (547, 911)
    ranks = tuple(rank_mod(active_gauge_rows, prime) for prime in primes)
    walsh_rows = [
        [
            sum(sign * active_gauge_rows[row][drift] for row, sign in enumerate(signs))
            for drift in range(P)
        ]
        for signs in WALSH_SIGNS
    ]
    raw_walsh_support = tuple(sum(value != 0 for value in row) for row in walsh_rows)
    spectra = tuple(fourier_support(walsh_rows, prime) for prime in primes)
    fourier_counts = tuple(
        tuple(sum(value != 0 for value in spectrum[row * P:(row + 1) * P]) for row in range(4))
        for spectrum in spectra
    )
    # Exact cyclotomic certificate: for a rational row P(x) of degree <=12,
    # P(zeta_13)=0 at a primitive root would force P=c*Phi_13.  Every row
    # has coefficient d=0 equal to zero and at least one other nonzero
    # coefficient, so it is not a constant-coefficient multiple of Phi_13.
    # Hence all twelve primitive modes are nonzero.  The zero mode is the
    # ordinary integer row sum and is checked directly.
    zero_mode_values = tuple(sum(row) for row in walsh_rows)
    require(
        zero_mode_values[0] != 0
        and zero_mode_values[1] == 0
        and zero_mode_values[2] == 0
        and zero_mode_values[3] != 0,
        ("Walsh zero-mode profile", zero_mode_values),
    )
    exact_fourier_counts = (13, 12, 12, 13)
    corner_zero_modes = tuple(sum(row) for row in active_gauge_rows)
    require(corner_zero_modes[0] == corner_zero_modes[3],
            ("diagonal chamber balance", corner_zero_modes))
    require(corner_zero_modes[1] == corner_zero_modes[2],
            ("off-diagonal chamber balance", corner_zero_modes))

    # Restore the common absolute offset before its character contraction.
    # This is the exact owner phase suggested by THM-3518.  Five good split
    # primes certify nonvanishing whenever at least one reduction survives.
    split_primes = (547, 911, 1093, 2003, 2549)
    owner_spectrum_records = []
    combined_owner_support = [[[False] * P for _walsh in range(4)] for _k in range(P)]
    for prime in split_primes:
        root = root_of_order_13(prime)
        prime_rows_by_k = []
        for owner_frequency in range(P):
            owner_rows = []
            active_pairs = ((0, 0), (0, 2), (2, 0), (2, 2))
            for left, right in active_pairs:
                owner_rows.append([
                    sum(
                        gauge_offset[left][right][drift][offset]
                        * pow(root, -owner_frequency * offset % P, prime)
                        for offset in range(P)
                    ) % prime
                    for drift in range(P)
                ])
            owner_walsh = [
                [
                    sum(sign * owner_rows[row][drift] for row, sign in enumerate(signs)) % prime
                    for drift in range(P)
                ]
                for signs in WALSH_SIGNS
            ]
            owner_spectra = [
                [
                    sum(
                        owner_walsh[row][drift] * pow(root, -frequency * drift % P, prime)
                        for drift in range(P)
                    ) % prime
                    for frequency in range(P)
                ]
                for row in range(4)
            ]
            for row in range(4):
                for frequency in range(P):
                    combined_owner_support[owner_frequency][row][frequency] |= (
                        owner_spectra[row][frequency] != 0
                    )
            prime_rows_by_k.append((
                owner_frequency,
                rank_mod(owner_rows, prime),
                tuple(sum(value != 0 for value in row) for row in owner_walsh),
                tuple(sum(value != 0 for value in row) for row in owner_spectra),
            ))
        owner_spectrum_records.append((prime, tuple(prime_rows_by_k)))

    combined_owner_counts = tuple(
        tuple(sum(combined_owner_support[k][row]) for row in range(4))
        for k in range(P)
    )
    # k=0 is decided exactly above; replace finite-reduction losses by the
    # Phi_13 certificate.  The positive owner modes are independently
    # certified coordinatewise by at least one good reduction.
    combined_owner_counts = (exact_fourier_counts,) + combined_owner_counts[1:]
    require(
        combined_owner_counts[0] == (13, 12, 12, 13)
        and all(counts == (13, 13, 13, 13) for counts in combined_owner_counts[1:]),
        ("owner-phase spectral closure", combined_owner_counts),
    )
    owner_rank_certified = tuple(
        max(prime_rows[k][1] for _prime, prime_rows in owner_spectrum_records)
        for k in range(P)
    )
    require(owner_rank_certified == (4,) * P,
            ("owner-phase K4 ranks", owner_rank_certified))

    defect_mass = []
    for defect in range(P):
        defect_mass.append(sum(
            bucket[left][right][drift][root_drift]
            for left in range(3)
            for right in range(3)
            for drift in range(P)
            for root_drift in range(P)
            if (drift + root_drift) % P == defect
        ))
    defect_support = sum(value != 0 for value in defect_mass)
    pair_gauge_digest = digest_json(pair_gauge_offset)

    record = (
        STAGE_SHA256,
        tuple(ATOMS),
        (nonempty_f_groups, nonempty_e_groups),
        tuple(len(group) for group in f_groups),
        tuple(len(group) for group in e_groups),
        tuple(len(profile[0]) for profile in u_profiles),
        tuple(len(profile[0]) for profile in v_profiles),
        atom_pair_nonzero,
        tuple(root_marginal),
        tuple(gauge_rows),
        tuple(opposite_rows),
        tuple(tuple(row) for row in active_gauge_rows),
        (full_gauge_support, active_gauge_support, opposite_support),
        ranks,
        tuple(tuple(row) for row in walsh_rows),
        raw_walsh_support,
        spectra,
        fourier_counts,
        zero_mode_values,
        exact_fourier_counts,
        corner_zero_modes,
        tuple(
            tuple(
                tuple(tuple(offsets) for offsets in drifts)
                for drifts in right_rows
            )
            for right_rows in gauge_offset
        ),
        pair_gauge_digest,
        tuple(owner_spectrum_records),
        combined_owner_counts,
        owner_rank_certified,
        tuple(defect_mass),
        total_mass,
        denominator,
    )
    semantic = digest_json(record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256, (semantic, EXPECTED_SEMANTIC_SHA256))

    print("== r=5 common-ancestry U_full guard-atom/root-drift probe ==")
    print(f"dependency={STAGE_PATH.name}:{STAGE_SHA256}")
    print("partition=39 actual guard atoms on each linked node X_(u,a),Y_(q,e'); one THM-2471 Boolean fibre")
    print(f"input_piece_counts=(f={len(f_pieces)},E={len(e_intervals)}); atom_groups_nonempty={(nonempty_f_groups,nonempty_e_groups)}/(39,39)")
    print(f"folded_profile_sizes=(arrival={tuple(len(profile[0]) for profile in u_profiles)},source={tuple(len(profile[0]) for profile in v_profiles)})")
    print(f"atom_pair_support={atom_pair_nonzero}/1521")
    print(f"root_colour_support={sum(value != 0 for value in root_marginal)}/13; s=0 mass={root_marginal[0]}")
    print(f"service_anchor=sum_M={Fraction(total_mass, denominator)}=169*I5: PASS")
    print("common_gauge=left_sheet-u=right_sheet-q iff s=-d; predicate applied before marginalization")
    print(f"common_gauge_support={full_gauge_support}/117; masks={tuple((left,right,tuple(int(value != 0) for value in row)) for left,right,row in gauge_rows)}")
    print(f"owner_active_common_gauge_support={active_gauge_support}/52; d=0 vanishes by source diagonal law")
    print(f"opposite_sign_s_equals_d_support={opposite_support}/117 (hostile orientation comparison)")
    print(f"active_K4_row_ranks_mod_{primes}={ranks}")
    print(f"active_Walsh_raw_support={raw_walsh_support}/13")
    print(f"active_Walsh_F13_Fourier_support_mod_{primes}={fourier_counts}")
    print(f"active_Walsh_zero_modes_nonzero={tuple(int(value != 0) for value in zero_mode_values)}")
    print(f"active_Walsh_F13_Fourier_support_exact={exact_fourier_counts} (Phi13 coefficient certificate)")
    print(f"active_corner_zero_modes={corner_zero_modes}; LL=RR and LR=RL exactly")
    print(f"common_offset_owner_phase_split_primes={split_primes}")
    print(f"pair_level_common_gauge_offset_sha256={pair_gauge_digest}")
    print(f"owner_frequency_x_Walsh_drift_Fourier_support_certified={combined_owner_counts}")
    print(f"owner_frequency_K4_ranks_certified={owner_rank_certified}")
    print(f"gauge_defect_mass_support={defect_support}/13")
    print(f"semantic_sha256={semantic}")
    print("scope=lawful guard-atom ancestry support and spectral capacity only; no endpoint weight identity, current, row exclusion, or LRC(14)")
    print("all exact checks passed")
    return {
        "gauge_offset": gauge_offset,
        "pair_gauge_offset": pair_gauge_offset,
        "active_gauge_rows": active_gauge_rows,
        "walsh_rows": walsh_rows,
        "root_marginal": root_marginal,
        "denominator": denominator,
        "semantic_sha256": semantic,
    }


if __name__ == "__main__":
    main()
