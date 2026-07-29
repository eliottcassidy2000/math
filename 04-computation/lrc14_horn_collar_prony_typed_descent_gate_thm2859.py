#!/usr/bin/env python3
"""Exact typed-descent gate for the THM-2859 Z8 edge and THM-2863 probes.

The calculation deliberately stops before forming Prony coefficients unless
one fully typed q3/q11 current transports across the physical step 2 -> 68
edge.  It distinguishes the strict same-labelled-cell forest fibre from the
weaker construction that attaches a shifted q-current to the surviving base
root.

No executable Python assertion statement is used.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from hashlib import sha256
from math import gcd
from pathlib import Path
import ast
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_digest(payload):
    return sha256(payload.replace(b"\r\n", b"\n")).hexdigest()


LOCAL_PINS = {
    COMP / "lrc14_horn_collar_endpoint_carry_thm2859.py":
        "6e062f3cc57c80fcff372c272bc138e280205bb953e484f1cc267340774260f0",
    COMP / "lrc14_horn_collar_endpoint_orbit_action_thm2859.py":
        "a4c145892ec3caf1f199a13bac07472e726f5dada0b523e86274ad8c14ce2846",
    COMP / "lrc14_q3_q11_transverse_endpoint_horn_thm2847.py":
        "258659c5093d98eea84056bdd3599b32d2a244bcd37dfa5f22dc5b25946ffe72",
    COMP / "lrc14_literal_fixed_sheet_allocation_thm2806.py":
        "311d0d85500f0c65945ebe5913f09d34a16293119c942b42eeaa854fbf85f71e",
}
for path, expected in LOCAL_PINS.items():
    require(
        lf_digest(path.read_bytes()) == expected,
        f"pinned local dependency changed: {path.name}",
    )


PRONY_COMMIT = "ab30502ddcd54fa9b4d597dd101184f2d3916b52"
GIT_PINS = {
    "04-computation/lrc14_endpoint_prony_carry_character3_thm2863.py":
        "3d09e2f06cd17c38a34b00fa133f67963287173e5a69092e4e2d2b190f4459b9",
    "05-knowledge/results/lrc14_endpoint_prony_carry_character3_thm2863.out":
        "2cf1902b19fd6ae9c6780521da0e431c8c7e1887235ea94ef4a8a5ca7cf9c7e4",
}
for relative, expected in GIT_PINS.items():
    payload = subprocess.check_output(
        ["git", "show", f"{PRONY_COMMIT}:{relative}"],
        cwd=ROOT,
    )
    require(
        lf_digest(payload) == expected,
        f"pinned THM-2863 blob changed: {relative}",
    )


import lrc14_horn_collar_endpoint_carry_thm2859 as edge
import lrc14_horn_collar_endpoint_orbit_action_thm2859 as complete_forest
import lrc14_q3_q11_transverse_endpoint_horn_thm2847 as horn


P = 13
QS = (3, 11)
STEPS = (2, 68)
FRAMES = ("source", "target")
MULTIPLIERS = (38, 64, 90, 116)
Z8 = (0, 8)
FACTOR_NAMES = ("E3", "clock", "q1", "q2", "c2", "c3")


def factor_signature(interval, module, e3, clocks, clock, s, t):
    allocation = horn.allocation
    return (
        allocation.contained(interval, e3),
        allocation.contained(interval, clocks[clock]),
        horn.safe_comb_contains(
            interval, module, module.W[1], 182,
            -14 * s - 13, -14 * s + 13,
        ),
        horn.safe_comb_contains(
            interval, module, module.W[2], 182,
            -14 * t - 13, -14 * t + 13,
        ),
        horn.safe_comb_contains(
            interval, module, module.C2, 182,
            14 * s - 13, 14 * s + 13,
        ),
        horn.safe_comb_contains(
            interval, module, module.C3, 182,
            14 * t - 13, 14 * t + 13,
        ),
    )


def circular_shift(interval, amount, period):
    return horn.circular_shift_interval(interval, amount, period)


def cyclic_gap_necklace(values):
    ordered = tuple(sorted(values))
    require(ordered, "empty axis has no cyclic gap necklace")
    gaps = tuple(
        (ordered[(index + 1) % len(ordered)] - ordered[index]) % P
        for index in range(len(ordered))
    )
    return min(
        gaps[index:] + gaps[:index]
        for index in range(len(gaps))
    )


def ancestry_arguments():
    ancestry = edge.independent.ancestry
    e_intervals = tuple(
        ancestry.base.build_set(
            ancestry.base.PAT_E, ancestry.base.ZELL
        )
    )
    q_intervals = tuple(
        ancestry.base.build_set(
            ancestry.host.PAT_QB, ancestry.base.ZELL
        )
    )
    q_depth, q_starts = ancestry.scaled_intervals(
        q_intervals, ancestry.DEPTH
    )
    e_depth_pack, e_depth_pack_starts = ancestry.scaled_intervals(
        e_intervals, ancestry.DEPTH * ancestry.PACK
    )
    e_depth, e_depth_starts = ancestry.scaled_intervals(
        e_intervals, ancestry.DEPTH
    )
    return (
        q_depth,
        q_starts,
        e_depth_pack,
        e_depth_pack_starts,
        e_depth,
        e_depth_starts,
    )


def main():
    require(
        MULTIPLIERS == tuple(12 + 26 * m for m in range(1, 5))
        and all(gcd(value, 91) == 1 for value in MULTIPLIERS),
        "THM-2863 multiplier progression changed",
    )

    paths, _horn_rows, _horn_common_atoms = edge.common42_path_geometry()
    edge_result = edge.step68_physical_audit(paths)
    candidates = tuple(
        row for row in paths
        if (
            row["cell"] in edge.HORN_20
            and row["distinguished"]
            and len(row["vertices"]) > 68
        )
    )
    expected_cells = {
        (1, s, t)
        for s in (0, 3, 12)
        for t in (5, 6, 9, 10)
    }
    require(
        len(candidates) == 12
        and {row["cell"] for row in candidates} == expected_cells,
        "physical T8 edge label bank changed",
    )
    require(
        edge_result["literal_ancestry_set_differences"] == (0, 0)
        and edge_result["literal_ancestry_counts"] == (966606, 28534)
        and edge_result["literal_ancestry_digest"]
        == "15c804c7cea9f61feab3b641eccdc035d937142b446d1cc14e059210eb1534fd",
        "base-root literal ancestry certificate changed",
    )

    allocation = horn.allocation
    (
        _module,
        full,
        details,
        e3,
        clocks,
        _q_pairs,
    ) = allocation.build_geometry()
    period = full.T
    unit = period // P
    half_step = edge.H
    shift = allocation.physical.SHIFT
    atom = allocation.ATOM_INTERVAL
    translated_atom = (
        atom[0] + 66 * half_step,
        atom[1] + 66 * half_step,
    )
    require(
        period == 742586 * half_step
        and unit == 57122 * half_step
        and candidates[0]["vertices"][2][:2] == atom
        and candidates[0]["vertices"][68][:2] == translated_atom,
        "physical scale or step-2/68 intervals changed",
    )
    require(
        {
            row["vertices"][2][:2] for row in candidates
        } == {atom}
        and {
            row["vertices"][68][:2] for row in candidates
        } == {translated_atom},
        "twelve labels ceased to collapse to one physical edge",
    )

    # Strict interpretation: the q-pulled source interval itself must be a
    # vertex of the complete 685-component collar forest in the same
    # labelled cell.  Also test the stronger cell-forgetting universe.
    (
        complete_paths,
        _complete_cells,
        _complete_module,
        _complete_e3,
        _complete_clocks,
        complete_kinds,
    ) = complete_forest.build_paths()
    cell_path_intervals = defaultdict(set)
    all_path_intervals = set()
    for row in complete_paths:
        cell_path_intervals[row["cell"]].update(
            piece[:2] for piece in row["vertices"]
        )
        all_path_intervals.update(
            piece[:2] for piece in row["vertices"]
        )
    require(
        len(complete_paths) == 685
        and sum(len(row["vertices"]) for row in complete_paths) == 63_895
        and complete_kinds
        == Counter({
            "cofiber_rooted": 587,
            "common_only": 98,
        }),
        "complete collar forest changed",
    )
    strict_state_rows = []
    strict_any_cell_rows = []
    strict_positive_controls = []
    strict_witnesses = {}
    for record in candidates:
        for step in STEPS:
            base = record["vertices"][step][:2]
            q0_interval = circular_shift(base, 0, period)
            strict_positive_controls.append(
                q0_interval in cell_path_intervals[record["cell"]]
            )
            for q in QS:
                pulled = circular_shift(base, q * unit, period)
                survives = pulled in cell_path_intervals[record["cell"]]
                strict_state_rows.append(
                    (record["cell"], step, q, survives)
                )
                strict_any_cell_rows.append(pulled in all_path_intervals)
                strict_witnesses.setdefault(
                    (step, q),
                    (record["cell"], base, pulled, survives),
                )
    require(
        all(strict_positive_controls)
        and not any(row[-1] for row in strict_state_rows)
        and not any(strict_any_cell_rows)
        and len(strict_state_rows) == 48,
        "strict complete-forest fibre boundary changed",
    )

    q_displacements_h = {
        q: (
            q * (unit // half_step)
            if q * (unit // half_step) <= period // (2 * half_step)
            else q * (unit // half_step) - period // half_step
        )
        for q in QS
    }
    require(
        q_displacements_h == {3: 171366, 11: -114244},
        "centered q displacements changed",
    )

    # Relaxed interpretation: keep the physical base root and attach the
    # q-shifted source/target current.  This tests whether the failure is
    # merely root naming or persists after all factors and carriers survive.
    relaxed_state_rows = []
    expected_weight = allocation.ATOM_WEIGHT
    for record in candidates:
        clock, s, t = record["cell"]
        source_base, target_base = details[clock]
        for step in STEPS:
            base = record["vertices"][step][:2]
            native_target = circular_shift(
                (base[0] + shift, base[1] + shift), 0, period
            )
            for q in QS:
                pulled_source = circular_shift(base, q * unit, period)
                pulled_target = circular_shift(
                    native_target, q * unit, period
                )
                source_q = allocation.physical.overlap.shift_weighted(
                    source_base, q * unit
                )
                target_q = allocation.physical.overlap.shift_weighted(
                    target_base, q * unit
                )
                signatures = tuple(
                    factor_signature(
                        interval, full, e3, clocks, clock, s, t
                    )
                    for interval in (
                        base,
                        native_target,
                        pulled_source,
                        pulled_target,
                    )
                )
                source_hit = horn.weighted_hit(pulled_source, source_q)
                target_hit = horn.weighted_hit(pulled_target, target_q)
                survives = (
                    all(all(signature) for signature in signatures)
                    and source_hit == ((*pulled_source, expected_weight),)
                    and target_hit == ((*pulled_target, expected_weight),)
                )
                relaxed_state_rows.append(
                    (
                        record["cell"],
                        step,
                        q,
                        signatures,
                        source_hit,
                        target_hit,
                        survives,
                    )
                )
    require(
        len(relaxed_state_rows) == 48
        and all(row[-1] for row in relaxed_state_rows),
        "relaxed q-current factor/carrier gate changed",
    )

    # Literal ancestry labels are constant across the T8 displacement in
    # every q-shifted frame.  Their U component is nevertheless empty on the
    # q3/q11 pulled sheets, another warning that this is not the base root.
    ancestry = edge.independent.ancestry
    ancestry_args = ancestry_arguments()
    ancestry_bank = {}
    for q in QS:
        for frame in FRAMES:
            extra = 0 if frame == "source" else shift
            for step, base in zip(STEPS, (atom, translated_atom)):
                interval = circular_shift(
                    (base[0] + extra, base[1] + extra),
                    q * unit,
                    period,
                )
                coordinate = sum(interval) // 2
                labels = ancestry.contributor_sets(
                    coordinate, *ancestry_args
                )
                ancestry_bank[(q, frame, step)] = labels
    ancestry_summary = {}
    for q in QS:
        rows = tuple(
            ancestry_bank[(q, frame, step)]
            for frame in FRAMES
            for step in STEPS
        )
        require(
            rows[0] == rows[1] == rows[2] == rows[3],
            f"q{q} ancestry labels changed across frame or T8 edge",
        )
        ancestry_summary[q] = (
            tuple(map(len, rows[0])),
            ancestry.path_digest(*rows[0]),
        )
    require(
        ancestry_summary
        == {
            3: (
                (0, 28535),
                "bdf95bb641f8d0ab0d8ebe764fecda782b9a7fe4995031c7f23e2f02868de931",
            ),
            11: (
                (0, 28534),
                "89a69bdfdedb80eeabcf8eb3401012714a33cc2cb9690221f835c923d4cf9dd2",
            ),
        },
        "q-shifted literal ancestry bank changed",
    )

    # Full endpoint masks are the decisive obstruction even under the
    # relaxed base-root convention.  The q0 positive control transports by
    # the unique Z8 translation.  Neither q3 nor q11 mask is in the same
    # translation orbit at the two ends, in either endpoint frame.
    mask_bank = {}
    interval_bank = {}
    for q in (0,) + QS:
        for frame in FRAMES:
            extra = 0 if frame == "source" else shift
            for step, base in zip(STEPS, (atom, translated_atom)):
                interval = circular_shift(
                    (base[0] + extra, base[1] + extra),
                    q * unit,
                    period,
                )
                interval_bank[(q, frame, step)] = interval
                mask_bank[(q, frame, step)] = edge.endpoint_mask(interval)

    q0_translations = tuple(
        edge.translation_deltas(
            mask_bank[(0, frame, 2)],
            mask_bank[(0, frame, 68)],
        )
        for frame in FRAMES
    )
    require(
        q0_translations == ((Z8,), (Z8,)),
        "q0 Z8 endpoint positive control changed",
    )

    endpoint_rows = {}
    for q in QS:
        start_source = mask_bank[(q, "source", 2)]
        end_source = mask_bank[(q, "source", 68)]
        require(
            start_source == mask_bank[(q, "target", 2)]
            and end_source == mask_bank[(q, "target", 68)],
            f"q{q} source/target endpoint frames diverged",
        )
        translations = tuple(
            edge.translation_deltas(
                mask_bank[(q, frame, 2)],
                mask_bank[(q, frame, 68)],
            )
            for frame in FRAMES
        )
        z8_hamming = tuple(
            sum(
                left != right
                for left, right in zip(
                    edge.translate_mask(
                        mask_bank[(q, frame, 2)], Z8
                    ),
                    mask_bank[(q, frame, 68)],
                )
            )
            for frame in FRAMES
        )
        translated = edge.translate_mask(start_source, Z8)
        lost = tuple(
            address
            for address, left, right in zip(
                edge.KEYS, translated, end_source
            )
            if left and not right
        )
        gained = tuple(
            address
            for address, left, right in zip(
                edge.KEYS, translated, end_source
            )
            if not left and right
        )
        start_description = edge.cartesian_description(start_source)
        end_description = edge.cartesian_description(end_source)
        start_first, start_second, _start_mass, start_cartesian = (
            start_description
        )
        end_first, end_second, _end_mass, end_cartesian = end_description
        require(
            start_cartesian
            and end_cartesian
            and cyclic_gap_necklace(start_second)
            != cyclic_gap_necklace(end_second),
            f"q{q} second-axis gap obstruction vanished",
        )
        endpoint_rows[q] = {
            "translations": translations,
            "z8_hamming": z8_hamming,
            "lost_gained": (
                len(lost), len(gained), lost[0], gained[0]
            ),
            "start_description": start_description,
            "end_description": end_description,
            "first_axis_necklaces": (
                cyclic_gap_necklace(start_first),
                cyclic_gap_necklace(end_first),
            ),
            "second_axis_necklaces": (
                cyclic_gap_necklace(start_second),
                cyclic_gap_necklace(end_second),
            ),
        }
    require(
        endpoint_rows[3]["translations"] == ((), ())
        and endpoint_rows[11]["translations"] == ((), ())
        and endpoint_rows[3]["z8_hamming"] == (80, 80)
        and endpoint_rows[11]["z8_hamming"] == (72, 72)
        and endpoint_rows[3]["lost_gained"][:2] == (40, 40)
        and endpoint_rows[11]["lost_gained"][:2] == (36, 36),
        "q3/q11 endpoint translation-orbit obstruction changed",
    )

    expected_descriptions = {
        3: (
            (
                tuple(range(10)),
                (0, 3, 4, 5, 8, 9, 10, 11, 12),
                90,
                True,
            ),
            (
                tuple(range(10)),
                (0, 5, 6, 7, 8, 9, 10, 11, 12),
                90,
                True,
            ),
        ),
        11: (
            (
                (0, 1, 2, 3, 4, 5, 8, 9, 12),
                (0, 1, 2, 3, 4, 5, 8, 11, 12),
                81,
                True,
            ),
            (
                (0, 1, 2, 3, 4, 5, 8, 9, 12),
                (0, 1, 2, 5, 6, 7, 8, 11, 12),
                81,
                True,
            ),
        ),
    }
    require(
        all(
            (
                endpoint_rows[q]["start_description"],
                endpoint_rows[q]["end_description"],
            ) == expected_descriptions[q]
            for q in QS
        ),
        "displayed endpoint rectangles changed",
    )

    # The two named THM-2863 columns do not see the lost rectangle geometry:
    # their occupancy is unchanged.  This is support data only; no endpoint
    # sum or frequency-cleared Prony coefficient is formed after gate failure.
    origins = (
        allocation.RIGHT_ORIGIN,
        allocation.add(
            allocation.RIGHT_ORIGIN, allocation.TARGET_STEP
        ),
    )
    origin_occupancy = {
        q: tuple(
            tuple(
                mask_bank[(q, "target", step)][edge.KEY_INDEX[origin]]
                for origin in origins
            )
            for step in STEPS
        )
        for q in QS
    }
    require(
        origin_occupancy
        == {
            3: ((1, 0), (1, 0)),
            11: ((1, 1), (1, 1)),
        },
        "named-origin support shadow changed",
    )

    strict_edge_rows = sum(
        all(
            row[-1]
            for row in strict_state_rows
            if row[0] == record["cell"] and row[2] == q
        )
        for record in candidates
        for q in QS
    )
    relaxed_pre_endpoint_edge_rows = len(candidates) * len(QS)
    relaxed_typed_edge_rows = sum(
        all(endpoint_rows[q]["translations"])
        for _record in candidates
        for q in QS
    )
    require(
        strict_edge_rows == 0
        and relaxed_pre_endpoint_edge_rows == 24
        and relaxed_typed_edge_rows == 0,
        "typed edge census changed",
    )
    probe_gate = tuple(
        (
            multiplier,
            strict_edge_rows,
            relaxed_pre_endpoint_edge_rows,
            relaxed_typed_edge_rows,
        )
        for multiplier in MULTIPLIERS
    )

    source_tree = ast.parse(Path(__file__).read_text())
    require(
        not any(
            isinstance(node, ast.Assert)
            for node in ast.walk(source_tree)
        ),
        "executable assertion statement found",
    )

    print("THM-2859 Z8 / PRONY TYPED DESCENT GATE")
    print(
        "pinned_local="
        + repr(tuple((path.name, digest) for path, digest in LOCAL_PINS.items()))
    )
    print(
        f"pinned_prony_commit={PRONY_COMMIT};"
        f"pinned_prony_blobs={tuple(GIT_PINS.items())}"
    )
    print(
        f"edge_labels={len(candidates)};physical_edge_rank=1;"
        f"source={atom};target={translated_atom};delta_h=66;"
        "endpoint_value=Z8"
    )
    print(
        f"q_centered_displacements_in_h={q_displacements_h};"
        f"strict_complete_forest_same_cell_state_rows=0/"
        f"{len(strict_state_rows)};"
        f"strict_complete_forest_any_cell_state_rows=0/"
        f"{len(strict_any_cell_rows)};"
        f"strict_same_sheet_edge_rows={strict_edge_rows}/24"
    )
    print(f"strict_minimal_witnesses={strict_witnesses}")
    print(
        f"relaxed_factor_carrier_state_rows={len(relaxed_state_rows)}/"
        f"{len(relaxed_state_rows)};"
        f"relaxed_pre_endpoint_edge_rows="
        f"{relaxed_pre_endpoint_edge_rows}/24;"
        f"factor_order={FACTOR_NAMES};carrier_weight={expected_weight}"
    )
    print(
        "base_root_ancestry="
        f"counts={edge_result['literal_ancestry_counts']};"
        f"digest={edge_result['literal_ancestry_digest']};"
        "symmetric_differences=(0,0)"
    )
    print(f"q_shifted_ancestry={ancestry_summary};U_empty_on_q3_q11=1")
    print(f"q0_endpoint_translation_positive_control={q0_translations}")
    for q in QS:
        row = endpoint_rows[q]
        print(
            f"q{q}_endpoint_rectangles="
            f"{row['start_description']}->{row['end_description']};"
            f"all_translation_deltas={row['translations']};"
            f"Z8_hamming={row['z8_hamming']};"
            f"Z8_lost_gained_first={row['lost_gained']};"
            f"second_axis_gap_necklaces={row['second_axis_necklaces']}"
        )
    print(f"named_prony_origin_occupancy_before_after={origin_occupancy}")
    print(
        "probe_gate_rows=(Y,strict_typed,relaxed_pre_endpoint,"
        f"relaxed_fully_typed)={probe_gate}"
    )
    print(
        "FIRST_STRICT_LOSS=physical same-labelled-cell forest vertex;"
        " q3/q11 pulled intervals are not forest vertices and have empty "
        "U ancestry"
    )
    print(
        "FIRST_RELAXED_LOSS=full endpoint-address transport;"
        " q3/q11 rectangles lie in different translation orbits in both frames"
    )
    print(
        "PRONY_SPLIT=SKIPPED_BY_TYPE_GATE;"
        " no scalar endpoint sums or Pi3(delta8-delta0) comparison formed"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
