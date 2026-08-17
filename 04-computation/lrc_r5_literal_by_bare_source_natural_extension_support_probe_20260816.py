#!/usr/bin/env python3
"""Put the literal BY support on THM-2471's bare-source stalk node.

The preceding one-leg Fubini probe integrated the U_full BY endpoint variable
independently of the r=5 ancestry variable.  THM-2471 already contains the
lawful pullback for the *bare source* endpoint:

    v = Y_(q,e)(y) = ((y+q)/13 + e)/13^5
                    = (y + q + 13e)/13^6.

This companion inserts the literal U_full BY endpoint *support* at that node
before the 13^5 transfer, the first-collision product, the seven-cell
integration, and the common-offset marginal.  It retains the endpoint atom,
the U_full (alpha,beta) shift, every guard shift tau, and (ell,c).

The probe is deliberately support-only: it does not insert the BY endpoint
character, AX remains preintegrated, and the U_full speed word is not the
canonical r=5 speed word.  Hence this is a lawful measurable one-leg cospan,
not a THM-2334 current, THM-2449 table, THM-2512 bridge, row exclusion, or
LRC(14) conclusion.
"""

from __future__ import annotations

from contextlib import redirect_stdout
from hashlib import sha256
import importlib.util
import io
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
P = 13
Q = 7
SELECTED_ADDRESS = (1, 0)

SPECTRUM_PATH = (
    ROOT
    / "04-computation/lrc_r5_source_aligned_relation_residue_7x13_spectrum_probe_20260816.py"
)
SPECTRUM_SHA256 = "5f3fbf08bef6f9a61e684f0f7616e80e1dbbda4f6bb2ed4ca3788d3b8b53d65a"
SPECTRUM_SEMANTIC = "cd55336bb1dfe5f37f020c242c4bca5b7c6be339ec57e95d69e10bbe68d9dbaa"
SOURCE_PATH = (
    ROOT
    / "04-computation/lrc_r5_source_aligned_guard_atom_branch_sidecar_probe_20260816.py"
)
SOURCE_SHA256 = "22c5c748392817ccc36889a007c65bd5f44b26c10638df6f6aac48e917547f41"
SOURCE_SEMANTIC = "31e9e90c63053944b590195555be07ccbf84fd4c7abc2101de6a2a3562202de6"
TARGET_PATH = (
    ROOT
    / "04-computation/lrc_ufull_owner_boundary_k4xf13_endpoint_factorization_independent_audit_20260816.py"
)
TARGET_SHA256 = "f89be10c65bb77270199f9399b155d5a2c82c0da121b3e8589fe3c1f7e9824fc"
TARGET_SEMANTIC = "d52c9f0a56c14a83e1e6b175c7b725314c99f09d44509bc8582847a5857f7da6"

EXPECTED_SEMANTIC_SHA256 = "cce1850a9b287f4b612c6724bfc63441febc78b2223b19d6d72e6989b2590500"


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_json(value: object) -> str:
    body = json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    return sha256(body).hexdigest()


def load_module(path: Path, name: str, expected_hash: str):
    require(lf_sha256(path) == expected_hash, (name, "source hash drift"))
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, (name, "module loader"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def scale_plain(
    intervals: tuple[tuple[int, int], ...] | list[tuple[int, int]], factor: int
) -> tuple[tuple[int, int], ...]:
    return tuple((left * factor, right * factor) for left, right in intervals)


def scale_weighted(
    pieces: tuple[tuple[int, int, int], ...] | list[tuple[int, int, int]],
    factor: int,
) -> tuple[tuple[int, int, int], ...]:
    return tuple((left * factor, right * factor, weight) for left, right, weight in pieces)


def intersect_plain(
    left: tuple[tuple[int, int], ...] | list[tuple[int, int]],
    right: tuple[tuple[int, int], ...] | list[tuple[int, int]],
) -> tuple[tuple[int, int], ...]:
    """Intersect sorted disjoint half-open interval unions."""
    output: list[tuple[int, int]] = []
    i = 0
    j = 0
    while i < len(left) and j < len(right):
        start = max(left[i][0], right[j][0])
        stop = min(left[i][1], right[j][1])
        if start < stop:
            if output and output[-1][1] == start:
                output[-1] = (output[-1][0], stop)
            else:
                output.append((start, stop))
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return tuple(output)


def interval_mass(intervals: tuple[tuple[int, int], ...]) -> int:
    return sum(right - left for left, right in intervals)


def source_geometry(source, common_grid: int):
    base = source.B
    stage = source.M
    source_grid = source.GRID
    source_scale = common_grid // source_grid
    e_intervals = stage.build_set(stage.PAT_E, stage.ZELL)
    q_intervals = stage.build_set(stage.PAT_QB, stage.ZELL)
    endpoint_groups_raw = base.partition_weighted_pieces(
        [(left, right, 1) for left, right in e_intervals]
    )
    source_groups_raw = tuple(
        source.source_current_group(group, q_intervals)
        for group in endpoint_groups_raw
    )
    endpoint_groups = tuple(
        tuple((left * source_scale, right * source_scale)
              for left, right, weight in group if weight)
        for group in endpoint_groups_raw
    )
    source_groups = tuple(
        scale_weighted(group, source_scale) for group in source_groups_raw
    )
    cells = tuple(
        scale_plain(tuple(cell), source_scale) for cell in stage.build_cells()
    )
    return endpoint_groups, source_groups, cells, source_scale


def target_groups(target, alpha: int, beta: int, target_scale: int):
    word, t_den, _nn, *_rest, atom_intervals = target.context()
    no_guard = dict(target.M.PATTERN_E)
    require(no_guard.pop(target.M.GUARD) == "guard_safe", "removed non-guard")
    shift = (0, -alpha % P, -beta % P, 0, 0, 0, 0, alpha, beta)
    intervals = target.M.fast_build_set(word, t_den, no_guard, shift)
    groups = target.partition_two_pointer(intervals, atom_intervals)
    return tuple(scale_plain(group, target_scale) for group in groups), len(intervals)


def make_windows(stage, groups, multiplier: int, grid: int):
    profiles = tuple(
        stage.weighted_fold(
            [(left, right, 1) for left, right in group], multiplier, grid
        )
        for group in groups
    )
    return tuple(
        tuple(stage.extract_window(starts, values, root, P, grid)
              for root in range(P))
        for starts, values in profiles
    )


def make_left_windows(stage, source_groups, grid: int):
    profiles = tuple(
        stage.weighted_fold(list(group), stage.DCOLL, grid)
        for group in source_groups
    )
    return tuple(
        tuple(stage.extract_window(starts, values, root, P, grid)
              for root in range(P))
        for starts, values in profiles
    )


def build_tensor(
    stage,
    atoms,
    left_windows,
    right_windows,
    cells,
    grid: int,
    right_cell_masks=None,
):
    tensor = [
        [[[0 for _offset in range(P)] for _ell in range(Q)] for _right in atoms]
        for _left in atoms
    ]
    for left_index, (left_sheet, _left_chamber) in enumerate(atoms):
        for right_index, (right_sheet, _right_chamber) in enumerate(atoms):
            drift = (right_sheet - left_sheet) % P
            root_drift = (-drift) % P
            masks = cells if right_cell_masks is None else right_cell_masks[right_index]
            for current_root in range(P):
                source_root = (current_root - root_drift) % P
                profile = stage.product_cum(
                    left_windows[left_index][current_root][0],
                    left_windows[left_index][current_root][1],
                    right_windows[right_index][source_root][0],
                    right_windows[right_index][source_root][1],
                    grid,
                )
                if not profile[3]:
                    continue
                masses = tuple(stage.set_integral(profile, mask) for mask in masks)
                require(sum(masses) <= profile[3], "cell masks exceed profile")
                if right_cell_masks is None:
                    require(sum(masses) == profile[3], "cell partition failed")
                offset = (left_sheet - current_root) % P
                require(
                    offset == (right_sheet - source_root) % P,
                    "common torsor offset mismatch",
                )
                for ell, mass in enumerate(masses):
                    tensor[left_index][right_index][ell][offset] += mass
    return tensor


def tensor_scaled_equal(left, right, scale: int) -> bool:
    return all(
        left[i][j][ell][offset] == scale * right[i][j][ell][offset]
        for i in range(len(left))
        for j in range(len(left[i]))
        for ell in range(Q)
        for offset in range(P)
    )


def tensor_stats(tensor) -> tuple[int, int, int]:
    values = [
        tensor[i][j][ell][offset]
        for i in range(len(tensor))
        for j in range(len(tensor[i]))
        for ell in range(Q)
        for offset in range(P)
    ]
    return sum(value != 0 for value in values), sum(values), len(values)


def guarded_coordinate_table(target, tensor, tau: int, drop_right: bool = False):
    table = [[0 for _offset in range(P)] for _ell in range(Q)]
    for left_index, (left_sheet, left_chamber) in enumerate(target.ATOMS):
        left_safe = target.safe(left_chamber, left_sheet + tau)
        if not left_safe:
            continue
        for right_index, (right_sheet, right_chamber) in enumerate(target.ATOMS):
            right_safe = 1 if drop_right else target.safe(
                right_chamber, right_sheet + tau
            )
            if not right_safe:
                continue
            for ell in range(Q):
                for offset in range(P):
                    table[ell][offset] += tensor[left_index][right_index][ell][offset]
    return tuple(tuple(row) for row in table)


def mismatch_count(left, right) -> int:
    return sum(
        a != b
        for left_matrix, right_matrix in zip(left, right)
        for left_row, right_row in zip(left_matrix, right_matrix)
        for a, b in zip(left_row, right_row)
    )


def tau_transform(bank, zeta: int, prime: int):
    return tuple(
        tuple(
            tuple(
                sum(
                    bank[tau][ell][offset]
                    * pow(zeta, (-t * tau) % P, prime)
                    for tau in range(P)
                )
                % prime
                for offset in range(P)
            )
            for ell in range(Q)
        )
        for t in range(P)
    )


def main() -> None:
    spectrum = load_module(SPECTRUM_PATH, "by_section_spectrum", SPECTRUM_SHA256)
    require(
        spectrum.EXPECTED_SEMANTIC_SHA256 == SPECTRUM_SEMANTIC,
        "spectrum semantic drift",
    )
    source = load_module(SOURCE_PATH, "by_section_source", SOURCE_SHA256)
    target = load_module(TARGET_PATH, "by_section_target", TARGET_SHA256)
    require(source.EXPECTED_SEMANTIC_SHA256 == SOURCE_SEMANTIC,
            "source semantic drift")
    require(target.EXPECTED_SEMANTIC_SHA256 == TARGET_SEMANTIC,
            "target semantic drift")
    require(tuple(source.ATOMS) == tuple(target.ATOMS), "atom chart mismatch")

    with redirect_stdout(io.StringIO()):
        source_data = source.main()
    require(source_data["semantic_sha256"] == SOURCE_SEMANTIC,
            "source replay semantic drift")
    old_tensor, _marginal, _cell_record, _pair_support, _entry_support = (
        spectrum.build_cell_tensor(source, source_data)
    )

    word, target_grid, nn, prime, root, zeta, *_rest = target.context()
    source_grid = source.GRID
    common_grid = math.lcm(source_grid, target_grid)
    source_scale = common_grid // source_grid
    target_scale = common_grid // target_grid
    require((source_scale, target_scale) == (2501183, 1540),
            "common-grid scales")
    require(nn % (P * Q) == 0, "split embedding lacks 91st roots")
    eta = pow(root, nn // Q, prime)
    require(pow(eta, Q, prime) == 1 and eta != 1, "order-seven root")
    require(pow(zeta, P, prime) == 1 and zeta != 1, "order-thirteen root")

    endpoint_groups, source_groups, cells, recovered_scale = source_geometry(
        source, common_grid
    )
    require(recovered_scale == source_scale, "source scaling drift")
    require(sum(bool(group) for group in endpoint_groups) == 20,
            "source endpoint group census")

    # All 169 target shifts: exact endpoint-level natural-section support census.
    address_rows = []
    selected_target_groups = None
    total_target_intervals = 0
    marginal_false_positives = 0
    all_intersection_support = 0
    for alpha in range(P):
        for beta in range(P):
            groups, interval_count = target_groups(target, alpha, beta, target_scale)
            total_target_intervals += interval_count
            intersections = tuple(
                intersect_plain(endpoint_groups[index], groups[index])
                for index in range(len(target.ATOMS))
            )
            source_masses = tuple(interval_mass(group) for group in endpoint_groups)
            target_masses = tuple(interval_mass(group) for group in groups)
            diagonal_masses = tuple(interval_mass(group) for group in intersections)
            false_positives = sum(
                source_mass > 0 and target_mass > 0 and diagonal_mass == 0
                for source_mass, target_mass, diagonal_mass in zip(
                    source_masses, target_masses, diagonal_masses
                )
            )
            marginal_false_positives += false_positives
            all_intersection_support += sum(value > 0 for value in diagonal_masses)
            right_guard_rows = tuple(
                (
                    sum(
                        diagonal_masses[index]
                        for index, (sheet, chamber) in enumerate(target.ATOMS)
                        if target.safe(chamber, sheet + tau)
                    ),
                    sum(
                        diagonal_masses[index] > 0
                        for index, (sheet, chamber) in enumerate(target.ATOMS)
                        if target.safe(chamber, sheet + tau)
                    ),
                )
                for tau in range(P)
            )
            address_rows.append(
                (
                    alpha,
                    beta,
                    len(groups),
                    sum(value > 0 for value in target_masses),
                    sum(value > 0 for value in diagonal_masses),
                    sum(diagonal_masses),
                    false_positives,
                    right_guard_rows,
                )
            )
            if (alpha, beta) == SELECTED_ADDRESS:
                selected_target_groups = groups
    require(selected_target_groups is not None, "selected endpoint address absent")
    require(all_intersection_support > 0, "all natural-section intersections empty")

    selected_diagonal_groups = tuple(
        intersect_plain(endpoint_groups[index], selected_target_groups[index])
        for index in range(len(target.ATOMS))
    )
    require(any(selected_diagonal_groups), "selected diagonal address empty")

    stage = source.M
    left_windows = make_left_windows(stage, source_groups, common_grid)
    identity_right_windows = make_windows(
        stage, endpoint_groups, stage.DCOLL, common_grid
    )
    diagonal_right_windows = make_windows(
        stage, selected_diagonal_groups, stage.DCOLL, common_grid
    )
    identity_tensor = build_tensor(
        stage,
        source.ATOMS,
        left_windows,
        identity_right_windows,
        cells,
        common_grid,
    )
    require(tensor_scaled_equal(identity_tensor, old_tensor, source_scale),
            "natural-extension identity recovery")
    diagonal_tensor = build_tensor(
        stage,
        source.ATOMS,
        left_windows,
        diagonal_right_windows,
        cells,
        common_grid,
    )

    # Hostile: impose the U_full endpoint set directly at the outer base y.
    # The source bare endpoint remains in its old atom; this is precisely the
    # ill-typed v=y replacement that the natural-extension map avoids.
    fiat_masks = tuple(
        tuple(
            intersect_plain(cells[ell], selected_target_groups[right_index])
            for ell in range(Q)
        )
        for right_index in range(len(target.ATOMS))
    )
    fiat_tensor = build_tensor(
        stage,
        source.ATOMS,
        left_windows,
        identity_right_windows,
        cells,
        common_grid,
        right_cell_masks=fiat_masks,
    )

    identity_stats = tensor_stats(identity_tensor)
    diagonal_stats = tensor_stats(diagonal_tensor)
    fiat_stats = tensor_stats(fiat_tensor)
    require(diagonal_stats[0] > 0 and diagonal_stats[1] > 0,
            "selected diagonal support vanished")
    require(diagonal_tensor != fiat_tensor, "v=Y and v=y tensors coincided")

    natural_guard_bank = tuple(
        guarded_coordinate_table(target, diagonal_tensor, tau) for tau in range(P)
    )
    fiat_guard_bank = tuple(
        guarded_coordinate_table(target, fiat_tensor, tau) for tau in range(P)
    )
    dropped_right_bank = tuple(
        guarded_coordinate_table(target, diagonal_tensor, tau, drop_right=True)
        for tau in range(P)
    )
    fiat_mismatches = mismatch_count(natural_guard_bank, fiat_guard_bank)
    guard_mismatches = mismatch_count(natural_guard_bank, dropped_right_bank)
    require(fiat_mismatches > 0, "fiat-base hostile invisible")
    require(guard_mismatches > 0, "right-guard hostile invisible")

    natural_tau_bank = tau_transform(natural_guard_bank, zeta, prime)
    fiat_tau_bank = tau_transform(fiat_guard_bank, zeta, prime)
    dropped_tau_bank = tau_transform(dropped_right_bank, zeta, prime)
    natural_shapes = tuple(
        spectrum.support_shape(spectrum.fourier_2d(table, eta, zeta, prime))
        for table in natural_tau_bank
    )
    fiat_shapes = tuple(
        spectrum.support_shape(spectrum.fourier_2d(table, eta, zeta, prime))
        for table in fiat_tau_bank
    )
    dropped_shapes = tuple(
        spectrum.support_shape(spectrum.fourier_2d(table, eta, zeta, prime))
        for table in dropped_tau_bank
    )

    source_word = tuple(stage.W)
    target_word = tuple(word)
    differing_roles = tuple(
        index for index, (left, right) in enumerate(zip(source_word, target_word))
        if left != right
    )
    require(differing_roles == (1, 3, 5), "speed-word boundary changed")
    require(stage.DCOLL == P**5, "r=5 collision depth")
    require(
        all(
            (q + P * sheet) % (P**6) == q + P * sheet
            for q in range(P)
            for sheet in (0, 1, P**5 - 1)
        ),
        "depth-six branch chart",
    )

    address_digest = digest_json(address_rows)
    identity_digest = digest_json(identity_tensor)
    diagonal_digest = digest_json(diagonal_tensor)
    fiat_digest = digest_json(fiat_tensor)
    natural_guard_digest = digest_json(natural_guard_bank)
    fiat_guard_digest = digest_json(fiat_guard_bank)
    dropped_guard_digest = digest_json(dropped_right_bank)
    record = (
        SPECTRUM_SHA256,
        SPECTRUM_SEMANTIC,
        SOURCE_SHA256,
        SOURCE_SEMANTIC,
        TARGET_SHA256,
        TARGET_SEMANTIC,
        source_grid,
        target_grid,
        common_grid,
        source_scale,
        target_scale,
        total_target_intervals,
        address_digest,
        all_intersection_support,
        marginal_false_positives,
        SELECTED_ADDRESS,
        tuple(interval_mass(group) for group in selected_diagonal_groups),
        identity_stats,
        diagonal_stats,
        fiat_stats,
        identity_digest,
        diagonal_digest,
        fiat_digest,
        natural_guard_digest,
        fiat_guard_digest,
        dropped_guard_digest,
        fiat_mismatches,
        guard_mismatches,
        digest_json(natural_tau_bank),
        digest_json(fiat_tau_bank),
        digest_json(dropped_tau_bank),
        natural_shapes,
        fiat_shapes,
        dropped_shapes,
        source_word,
        target_word,
        differing_roles,
    )
    semantic = digest_json(record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("semantic", semantic, EXPECTED_SEMANTIC_SHA256))

    print("== r=5 literal-BY bare-source natural-extension support probe ==")
    print(
        "status=FINITE-EXACT LAWFUL COORDINATE-LEVEL ONE-LEG SUPPORT COSPAN; "
        "not a current or LRC(14) closure"
    )
    print(
        f"dependencies=((spectrum,{SPECTRUM_SHA256},{SPECTRUM_SEMANTIC}),"
        f"(source,{SOURCE_SHA256},{SOURCE_SEMANTIC}),"
        f"(target,{TARGET_SHA256},{TARGET_SEMANTIC}))"
    )
    print(
        "natural_extension=v=Y_(q,e')(y)=(y+q+13e')/13^6; "
        "branch_address=q+13e' in Z/13^6; v is not y"
    )
    print(
        f"grids=(source={source_grid},target={target_grid},lcm={common_grid},"
        f"scales={source_scale},{target_scale})"
    )
    print(
        f"all_address_endpoint_census=(alpha_beta=169,target_intervals="
        f"{total_target_intervals},diagonal_atom_support={all_intersection_support}/6591,"
        f"positive_marginal_but_empty_diagonal={marginal_false_positives},"
        f"sha256={address_digest})"
    )
    print(
        f"selected_address={SELECTED_ADDRESS}; identity_recovery="
        f"old_tensor*{source_scale}: PASS"
    )
    print(
        f"identity_tensor=(support={identity_stats[0]}/{identity_stats[2]},"
        f"mass={identity_stats[1]},sha256={identity_digest})"
    )
    print(
        f"natural_diagonal_tensor=(support={diagonal_stats[0]}/{diagonal_stats[2]},"
        f"mass={diagonal_stats[1]},sha256={diagonal_digest})"
    )
    print(
        f"hostile_v_equals_y=(support={fiat_stats[0]}/{fiat_stats[2]},"
        f"mass={fiat_stats[1]},guarded_coordinate_mismatches={fiat_mismatches}/1183,"
        f"sha256={fiat_digest})"
    )
    print(
        f"guarded_natural_bank=(tau_x_ell_x_c=1183,sha256={natural_guard_digest},"
        f"tau_fourier_sha256={digest_json(natural_tau_bank)})"
    )
    print(
        f"hostile_right_guard_deleted=(coordinate_mismatches={guard_mismatches}/1183,"
        f"sha256={dropped_guard_digest})"
    )
    print(f"natural_tau_spectrum_shapes={natural_shapes}")
    print(f"fiat_tau_spectrum_shapes={fiat_shapes}")
    print(f"right_guard_dropped_tau_spectrum_shapes={dropped_shapes}")
    print(
        f"tuple_boundary=(source_word={source_word},U_full_word={target_word},"
        f"differing_roles={differing_roles})"
    )
    print(
        "preserved=bare-source endpoint node Y,source and U_full endpoint atoms,"
        "alpha,beta,all tau guards,ell,c,depth-six branch sum"
    )
    print(
        "lost_or_missing=BY character,AX common-stalk lift,exact grouped address,"
        "U_full/source same-tuple theorem,chronology,physical current"
    )
    print(
        "nonconsequence=no THM-2449/2512 input,current,row exclusion,or LRC(14) closure"
    )
    print(f"semantic_sha256={semantic}")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
