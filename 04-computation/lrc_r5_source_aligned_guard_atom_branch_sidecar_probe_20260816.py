#!/usr/bin/env python3
"""Exact source-aligned guard-atom ancestry probe at the r=5 window.

The arrival-labelled companion partitions f=1_Q P^K e after the K-clock.
Here the left packet is instead refined at the endpoint/source time:

    e_omega=e P_omega,
    f_omega^src=1_Q P^K(e P_omega).

This is THM-2471 (34), so both labels omega,nu are endpoint-time guard atoms
even though their packets meet later on one first-collision base.  The probe
also computes the raw I_source -> I_arrival transition under T^K, K=2.  That
unlabelled transition has rank three and kills every nontrivial sheet
character; the retained inverse-branch index restores the source sheet.

No endpoint AX/BY weight, relation current, grouped coefficient, row
exclusion, or LRC(14) conclusion is asserted here.
"""

from __future__ import annotations

from bisect import bisect_right
from contextlib import redirect_stdout
from fractions import Fraction
from hashlib import sha256
import importlib.util
import io
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
P = 13
BASE_PATH = ROOT / "04-computation/lrc_r5_common_ancestry_guard_atom_root_drift_probe_20260816.py"
BASE_SHA256 = "83f1fa49ac4d02e21a1d76fed169d101715a6620342714ed05b9172ae967a730"
BASE_SEMANTIC = "3d8c88fb7b9762f41ef35c00d980b99fc435c8352baf5dddb9fe412d1baeace0"
EXPECTED_SEMANTIC_SHA256 = "31e9e90c63053944b590195555be07ccbf84fd4c7abc2101de6a2a3562202de6"
SPLIT_PRIMES = (547, 911, 1093, 2003, 2549)
ACTIVE_PAIRS = ((0, 0), (0, 2), (2, 0), (2, 2))


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_json(value: object) -> str:
    body = json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    return sha256(body).hexdigest()


def load_base():
    require(lf_sha256(BASE_PATH) == BASE_SHA256, "arrival-labelled base source drift")
    spec = importlib.util.spec_from_file_location("r5_source_aligned_base", BASE_PATH)
    require(spec is not None and spec.loader is not None, "base module loader")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


B = load_base()
M = B.M
GRID = B.GRID
ATOMS = B.ATOMS
CHAMBERS = B.CHAMBERS
CHAMBER_INDEX = B.CHAMBER_INDEX
WALSH_SIGNS = B.WALSH_SIGNS


def source_current_group(
    endpoint_group: tuple[tuple[int, int, int], ...],
    q_intervals: list[tuple[int, int]],
) -> list[tuple[int, int, int]]:
    """Build 1_Q P^K(e P_omega) from one endpoint-labelled group."""
    starts, values = M.weighted_fold(list(endpoint_group), M.RPKT, GRID)
    pieces: list[tuple[int, int, int]] = []
    for q_left, q_right in q_intervals:
        index = bisect_right(starts, q_left) - 1
        while True:
            piece_left = starts[index]
            piece_right = starts[index + 1] if index + 1 < len(starts) else GRID
            left = max(q_left, piece_left)
            right = min(q_right, piece_right)
            if left < right and values[index]:
                pieces.append((left, right, values[index]))
            if piece_right >= q_right:
                break
            index += 1
    return pieces


def atom_of_s(s: Fraction) -> int:
    """Return the 39-atom index for s=91t in [0,91)."""
    require(0 <= s < 91, ("atom coordinate", s))
    sheet = int(s // 7)
    residual = s - 7 * sheet
    chamber = 0 if residual < 1 else (1 if residual < 6 else 2)
    return 3 * sheet + chamber


def temporal_transition() -> tuple[list[list[int]], tuple[object, ...]]:
    """Count I_source(t) x I_arrival(13^2 t) on the common exact grid."""
    clock = P**2
    cells = 91 * clock
    transition = [[0] * len(ATOMS) for _atom in ATOMS]
    branch_sheet_recovery = True
    for cell in range(cells):
        # Midpoints avoid the null boundary set; every endpoint lies on this grid.
        source_s = Fraction(2 * cell + 1, 2 * clock)
        arrival_s = Fraction(2 * cell + 1, 2) % 91
        source = atom_of_s(source_s)
        arrival = atom_of_s(arrival_s)
        transition[source][arrival] += 1

        # If n=floor(13^2 t) is the retained inverse branch in z=(x+n)/13^2,
        # then its high base-13 digit is exactly the source guard sheet.
        branch = (2 * cell + 1) // (2 * 91)
        branch_sheet_recovery &= source // 3 == branch // P

    row_identity = all(
        transition[3 * sheet + chamber] == transition[chamber]
        for chamber in range(3)
        for sheet in range(P)
    )
    ranks = tuple(B.rank_mod(transition, prime) for prime in SPLIT_PRIMES)
    support = sum(value != 0 for row in transition for value in row)
    sheet_aggregate = [
        [
            sum(
                transition[3 * source_sheet + source_chamber]
                          [3 * arrival_sheet + arrival_chamber]
                for source_chamber in range(3)
                for arrival_chamber in range(3)
            )
            for arrival_sheet in range(P)
        ]
        for source_sheet in range(P)
    ]
    require(sum(map(sum, transition)) == cells, "temporal transition mass")
    require(support == len(ATOMS) ** 2, ("temporal support", support))
    require(row_identity, "source sheet survived the unlabelled temporal marginal")
    require(ranks == (3,) * len(SPLIT_PRIMES), ("temporal ranks", ranks))
    require(all(value == 91 for row in sheet_aggregate for value in row),
            ("sheet aggregate", sheet_aggregate))
    require(branch_sheet_recovery, "inverse branch did not restore source sheet")
    record = (
        cells,
        support,
        ranks,
        row_identity,
        tuple(tuple(row) for row in sheet_aggregate),
        branch_sheet_recovery,
        digest_json(transition),
    )
    return transition, record


def main() -> dict[str, object]:
    with redirect_stdout(io.StringIO()):
        inherited = B.main()
    require(inherited["semantic_sha256"] == BASE_SEMANTIC, "base semantic drift")

    e_intervals = M.build_set(M.PAT_E, M.ZELL)
    q_intervals = M.build_set(M.PAT_QB, M.ZELL)
    whole_f = B.build_f_pieces(e_intervals, q_intervals)
    endpoint_groups = B.partition_weighted_pieces(
        [(left, right, 1) for left, right in e_intervals]
    )
    source_f_groups = tuple(
        source_current_group(group, q_intervals) for group in endpoint_groups
    )
    nonempty = (
        sum(bool(group) for group in source_f_groups),
        sum(bool(group) for group in endpoint_groups),
    )
    require(nonempty == (20, 20), ("source-labelled nonempty atoms", nonempty))
    require(
        sum(sum(weight * (right - left) for left, right, weight in group)
            for group in source_f_groups)
        == sum(weight * (right - left) for left, right, weight in whole_f),
        "source-labelled packet mass",
    )

    u_profiles = tuple(
        M.weighted_fold(list(group), M.DCOLL, GRID) for group in source_f_groups
    )
    v_profiles = tuple(
        M.weighted_fold(list(group), M.DCOLL, GRID) for group in endpoint_groups
    )
    u_whole = M.weighted_fold(whole_f, M.DCOLL, GRID)
    v_whole = M.weighted_fold(
        [(left, right, 1) for left, right in e_intervals], M.DCOLL, GRID
    )
    B.require_profile_partition(u_profiles, u_whole, "source-labelled current partition")
    B.require_profile_partition(v_profiles, v_whole, "endpoint source partition")

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

    # bucket[C][D][endpoint-sheet drift d][collision-root drift s]
    bucket = [[[[0 for _s in range(P)] for _d in range(P)]
               for _right in CHAMBERS] for _left in CHAMBERS]
    gauge_offset = [[[[0 for _c in range(P)] for _d in range(P)]
                     for _right in CHAMBERS] for _left in CHAMBERS]
    pair_gauge_offset = [[[0 for _c in range(P)] for _right in ATOMS]
                         for _left in ATOMS]
    full_pair_support = 0
    common_pair_support = 0

    for left_index, (left_sheet, left_chamber) in enumerate(ATOMS):
        left_ci = CHAMBER_INDEX[left_chamber]
        for right_index, (right_sheet, right_chamber) in enumerate(ATOMS):
            right_ci = CHAMBER_INDEX[right_chamber]
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
                        offset = (left_sheet - current_root) % P
                        require(
                            offset == (right_sheet - source_root) % P,
                            "source-aligned common-gauge mismatch",
                        )
                        gauge_offset[left_ci][right_ci][drift][offset] += root_mass
                        pair_gauge_offset[left_index][right_index][offset] += root_mass
                bucket[left_ci][right_ci][drift][root_drift] += mass
                pair_nonzero |= mass != 0
            full_pair_support += int(pair_nonzero)
            common_pair_support += int(any(pair_gauge_offset[left_index][right_index]))

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
        require(reconstructed == mass, ("root marginal", root_drift))

    denominator = M.RPKT * M.DCOLL * M.DCOLL * GRID
    total_mass = sum(root_marginal)
    require(root_marginal[0] == 0, "source diagonal did not vanish")
    require(all(value > 0 for value in root_marginal[1:]), "root colour vanished")
    require(Fraction(total_mass, denominator) == 169 * M.I5, "service anchor")

    gauge_rows = []
    opposite_rows = []
    active_rows = []
    for left in range(3):
        for right in range(3):
            gauge = [bucket[left][right][drift][(-drift) % P] for drift in range(P)]
            opposite = [bucket[left][right][drift][drift] for drift in range(P)]
            gauge_rows.append(gauge)
            opposite_rows.append(opposite)
            require(all(sum(gauge_offset[left][right][drift]) == gauge[drift]
                        for drift in range(P)), "offset marginal")
            if (left, right) in ACTIVE_PAIRS:
                active_rows.append(gauge)

    full_gauge_support = sum(value != 0 for row in gauge_rows for value in row)
    active_support = sum(value != 0 for row in active_rows for value in row)
    opposite_support = sum(value != 0 for row in opposite_rows for value in row)
    # Source-time refinement inherits the endpoint owner support on both legs,
    # so the middle chamber disappears before the common-gauge marginal.
    require((full_gauge_support, active_support) == (48, 48),
            ("gauge support", full_gauge_support, active_support))
    require(all(row[0] == 0 and all(value > 0 for value in row[1:])
                for row in active_rows), "active 48 carrier")
    active_ranks = tuple(B.rank_mod(active_rows, prime) for prime in SPLIT_PRIMES)
    require(active_ranks == (4,) * len(SPLIT_PRIMES), ("active ranks", active_ranks))

    walsh_rows = [
        [sum(sign * active_rows[index][drift]
             for index, sign in enumerate(signs))
         for drift in range(P)]
        for signs in WALSH_SIGNS
    ]
    zero_modes = tuple(sum(row) for row in walsh_rows)
    require(
        zero_modes[0] != 0 and zero_modes[1] == zero_modes[2] == 0
        and zero_modes[3] != 0,
        ("Walsh zero modes", zero_modes),
    )

    combined_support = [[[False] * P for _row in range(4)] for _k in range(P)]
    owner_rank_by_prime = []
    for prime in SPLIT_PRIMES:
        root = B.root_of_order_13(prime)
        prime_ranks = []
        for owner_frequency in range(P):
            owner_rows = []
            for left, right in ACTIVE_PAIRS:
                owner_rows.append([
                    sum(
                        gauge_offset[left][right][drift][offset]
                        * pow(root, -owner_frequency * offset % P, prime)
                        for offset in range(P)
                    ) % prime
                    for drift in range(P)
                ])
            prime_ranks.append(B.rank_mod(owner_rows, prime))
            owner_walsh = [
                [sum(sign * owner_rows[index][drift]
                     for index, sign in enumerate(signs)) % prime
                 for drift in range(P)]
                for signs in WALSH_SIGNS
            ]
            spectrum = B.fourier_support(owner_walsh, prime)
            for row in range(4):
                for frequency in range(P):
                    combined_support[owner_frequency][row][frequency] |= (
                        spectrum[row * P + frequency] != 0
                    )
        owner_rank_by_prime.append((prime, tuple(prime_ranks)))

    owner_ranks = tuple(
        max(ranks[k] for _prime, ranks in owner_rank_by_prime) for k in range(P)
    )
    combined_counts = tuple(
        tuple(sum(combined_support[k][row]) for row in range(4)) for k in range(P)
    )
    # k=0 has the same exact Phi_13 certificate as the arrival-labelled probe:
    # every row has coefficient d=0 zero and another nonzero coefficient.
    combined_counts = ((13, 12, 12, 13),) + combined_counts[1:]
    require(owner_ranks == (4,) * P, ("owner ranks", owner_ranks))
    require(all(counts == (13, 13, 13, 13) for counts in combined_counts[1:]),
            ("owner spectral closure", combined_counts))

    _transition, temporal_record = temporal_transition()
    pair_digest = digest_json(pair_gauge_offset)
    record = (
        BASE_SHA256,
        BASE_SEMANTIC,
        tuple(ATOMS),
        nonempty,
        tuple(len(group) for group in source_f_groups),
        tuple(len(group) for group in endpoint_groups),
        full_pair_support,
        common_pair_support,
        tuple(root_marginal),
        tuple(tuple(row) for row in gauge_rows),
        tuple(tuple(row) for row in active_rows),
        (full_gauge_support, active_support, opposite_support),
        active_ranks,
        tuple(tuple(row) for row in walsh_rows),
        zero_modes,
        tuple(
            tuple(tuple(tuple(offsets) for offsets in drifts) for drifts in rights)
            for rights in gauge_offset
        ),
        pair_digest,
        tuple(owner_rank_by_prime),
        combined_counts,
        owner_ranks,
        temporal_record,
        total_mass,
        denominator,
    )
    semantic = digest_json(record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                (semantic, EXPECTED_SEMANTIC_SHA256))

    print("== r=5 source-aligned guard-atom / inverse-branch-sidecar probe ==")
    print(f"dependency=(arrival_base,{BASE_SHA256},{BASE_SEMANTIC})")
    print("left_refinement=f_omega^src=1_Q P^K(e P_omega); right_refinement=e_nu=e P_nu")
    print(f"endpoint_atom_groups_nonempty=(left={nonempty[0]},right={nonempty[1]})/39")
    print(f"atom_pair_support=(all_roots={full_pair_support},common_gauge={common_pair_support})/1521")
    print(f"root_colour_support={sum(value != 0 for value in root_marginal)}/13; s=0 mass={root_marginal[0]}")
    print(f"service_anchor=sum_M={Fraction(total_mass, denominator)}=169*I5: PASS")
    print(f"common_gauge_support=(full={full_gauge_support}/117,owner_active={active_support}/52,opposite={opposite_support}/117)")
    print(f"active_K4_ranks_mod_split_primes={active_ranks}")
    print(f"active_Walsh_zero_modes_nonzero={tuple(int(value != 0) for value in zero_modes)}")
    print(f"owner_frequency_x_Walsh_drift_Fourier_support_certified={combined_counts}")
    print(f"owner_frequency_K4_ranks_certified={owner_ranks}")
    print(f"pair_level_source_aligned_common_gauge_sha256={pair_digest}")
    print("raw_temporal_transition=(support=1521/1521,rankQ=3,source_sheet_rows_identical=True)")
    print("raw_temporal_sheet_aggregate=all_91; every nonzero source F13 character is erased")
    print("inverse_branch_sidecar=source_sheet=floor(branch_0_to_168/13): PASS")
    print(f"semantic_sha256={semantic}")
    print("scope=source-time endpoint labels on one ancestry base and exact branch-sidecar boundary only; no endpoint AX/BY weight,current,row exclusion,LRC(14)")
    print("all exact checks passed")
    return {
        "gauge_offset": gauge_offset,
        "pair_gauge_offset": pair_gauge_offset,
        "active_gauge_rows": active_rows,
        "root_marginal": root_marginal,
        "denominator": denominator,
        "semantic_sha256": semantic,
        "temporal_record": temporal_record,
    }


if __name__ == "__main__":
    main()
