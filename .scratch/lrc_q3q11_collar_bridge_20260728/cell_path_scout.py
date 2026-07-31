#!/usr/bin/env python3
"""Exact cell/path scout for the THM-2847-to-THM-2825 bridge.

This is scratch evidence only.  It asks how the 42 q3/q11 cells and the
20-cell E3-only horn sit inside the 587 rooted half-step paths.
"""

from collections import Counter, defaultdict
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))

import lrc14_nearest_half_step_common_right_collar_thm2825 as collar


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


COMMON_42 = {
    (1, s, t)
    for s in (0, 3, 8, 9, 10, 11, 12)
    for t in (5, 6, 7, 8, 9, 10)
}
HORN_20 = {
    (1, s, t)
    for s in (0, 3, 8, 9, 12)
    for t in (5, 6, 9, 10)
}
SELECTED_I = (142004992589460, 142005019034340)


def shifted_left(left, steps):
    result = left + steps * collar.H
    require(
        0 <= result
        and result + collar.copies.LENGTH <= collar.copies.T,
        "path met the ambient circle seam",
    )
    return result


def shifted_interval(interval, shift):
    left = interval[0] + shift
    right = interval[1] + shift
    require(
        0 <= left < right <= collar.copies.T,
        "tested factor interval met the ambient seam",
    )
    return left, right


def safe_comb_contains(interval, module, speed, denominator, lo, hi):
    require(
        module.T % (denominator * speed) == 0,
        "comb grid ceased to resolve",
    )
    unit = module.T // (denominator * speed)
    step = denominator * unit
    length = (hi - lo) * unit
    base = (lo % denominator) * unit
    left, right = interval
    k0 = (left - base) // step
    for k in range(k0 - 1, k0 + 3):
        forbidden_left = base + k * step
        forbidden_right = forbidden_left + length
        if max(left, forbidden_left) < min(right, forbidden_right):
            return False
    return True


def factor_signature(interval, full_module, e3, clocks, clock, s, t):
    return (
        collar.copies.contained(interval, e3),
        collar.copies.contained(interval, clocks[clock]),
        safe_comb_contains(
            interval, full_module, full_module.W[1], 182,
            -14 * s - 13, -14 * s + 13,
        ),
        safe_comb_contains(
            interval, full_module, full_module.W[2], 182,
            -14 * t - 13, -14 * t + 13,
        ),
        safe_comb_contains(
            interval, full_module, full_module.C2, 182,
            14 * s - 13, 14 * s + 13,
        ),
        safe_comb_contains(
            interval, full_module, full_module.C3, 182,
            14 * t - 13, 14 * t + 13,
        ),
    )


def main():
    (
        _module,
        _rails,
        _present,
        details,
        full_module,
        e3,
        clocks,
        q_pairs,
        _delayed,
        _source_weight,
        _target_weight,
        _rail_common,
    ) = collar.copies.physical_setup()

    semantic_caches = tuple({} for _clock in range(7))
    records = []
    cell_records = defaultdict(list)
    cell_common_counts = {}
    nonempty_cells = set()

    for clock in range(7):
        for s in collar.COMMON_S:
            for t in collar.COMMON_T:
                _source, _target, common, right = collar.cell_objects(
                    details, full_module, e3, clocks, clock, s, t
                )
                if not right:
                    continue
                cell = (clock, s, t)
                nonempty_cells.add(cell)
                cell_common_counts[cell] = len(common)
                common_by_left = {piece[0]: piece for piece in common}
                require(
                    len(common_by_left) == len(common),
                    "common left endpoints collided",
                )
                for root_index, root in enumerate(right):
                    path = []
                    step = 1
                    while shifted_left(root[0], step) in common_by_left:
                        path.append(
                            common_by_left[shifted_left(root[0], step)]
                        )
                        step += 1
                    require(len(path) >= 2, "root lost collar ladder")
                    values = tuple(
                        collar.semantic_value(
                            piece,
                            q_pairs[clock],
                            semantic_caches[clock],
                        )
                        != (0, 0)
                        for piece in (root, path[0], path[1])
                    )
                    require(
                        values[0] == values[2]
                        and values[0] != values[1],
                        "semantic collar grading changed",
                    )
                    factor_rows = None
                    if cell in COMMON_42:
                        factor_rows = tuple(
                            factor_signature(
                                piece[:2],
                                full_module,
                                e3,
                                clocks,
                                clock,
                                s,
                                t,
                            )
                            for piece in (root, path[0], path[1])
                        )
                        shifted_factor_rows = tuple(
                            factor_signature(
                                shifted_interval(
                                    piece[:2], collar.copies.SHIFT
                                ),
                                full_module,
                                e3,
                                clocks,
                                clock,
                                s,
                                t,
                            )
                            for piece in (root, path[0], path[1])
                        )
                    record = {
                        "cell": cell,
                        "root_index": root_index,
                        "root_left": root[0],
                        "path_length": len(path),
                        "semantic": values,
                        "root_mod_13": (root[0] // collar.H) % 13,
                        "factor_rows": factor_rows,
                        "shifted_factor_rows": shifted_factor_rows
                        if cell in COMMON_42 else None,
                    }
                    records.append(record)
                    cell_records[cell].append(record)

    require(len(records) == 587, "global rooted-path count changed")
    require(
        COMMON_42 <= nonempty_cells and HORN_20 <= COMMON_42,
        "q3/q11 bank left the nonempty collar bank",
    )

    def scope_records(cells):
        return [record for record in records if record["cell"] in cells]

    bank = scope_records(COMMON_42)
    horn = scope_records(HORN_20)
    extra = scope_records(COMMON_42 - HORN_20)
    outside = scope_records(nonempty_cells - COMMON_42)

    def summarize(rows):
        def holes(mask):
            return tuple(
                collar.copies.FACTOR_NAMES[index]
                for index, value in enumerate(mask)
                if not value
            )

        summary = {
            "roots": len(rows),
            "root_multiplicity": Counter(
                len(cell_records[cell])
                for cell in {row["cell"] for row in rows}
            ),
            "path_lengths": Counter(row["path_length"] for row in rows),
            "semantic": Counter(row["semantic"] for row in rows),
            "root_mod_13": Counter(row["root_mod_13"] for row in rows),
            "survive_14h": Counter(
                row["path_length"] >= 14 for row in rows
            ),
        }
        typed_rows = [row for row in rows if row["factor_rows"] is not None]
        if typed_rows:
            summary.update({
                "root_holes": Counter(
                    holes(row["factor_rows"][0]) for row in typed_rows
                ),
                "m1_all_six": Counter(
                    all(row["factor_rows"][1]) for row in typed_rows
                ),
                "m2_all_six": Counter(
                    all(row["factor_rows"][2]) for row in typed_rows
                ),
                "root_shifted_all_six": Counter(
                    all(row["shifted_factor_rows"][0])
                    for row in typed_rows
                ),
            })
        return summary

    cell_types = Counter(
        (
            "horn" if cell in HORN_20 else "extra",
            len(rows),
            tuple(sorted(record["path_length"] for record in rows)),
        )
        for cell, rows in cell_records.items()
        if cell in COMMON_42
    )
    root_index_profiles = Counter(
        (
            "horn" if cell in HORN_20 else "extra",
            tuple(record["root_index"] for record in rows),
        )
        for cell, rows in cell_records.items()
        if cell in COMMON_42
    )
    distinguished_42 = []
    for cell in sorted(COMMON_42):
        hits = [
            row
            for row in cell_records[cell]
            if (
                row["root_left"] + 2 * collar.H,
                row["root_left"] + 2 * collar.H + collar.copies.LENGTH,
            )
            == SELECTED_I
        ]
        require(len(hits) == 1, f"selected I has {len(hits)} roots in {cell}")
        distinguished_42.append(hits[0])
    distinguished = [
        row for row in distinguished_42 if row["cell"] in HORN_20
    ]
    distinguished_physical_roots = {
        (
            row["root_left"],
            row["root_left"] + collar.copies.LENGTH,
        )
        for row in distinguished
    }
    require(
        len(distinguished) == 20
        and len(distinguished_physical_roots) == 1,
        "distinguished horn roots stopped collapsing to one physical copy",
    )
    selected_root = next(iter(distinguished_physical_roots))
    selected_m1 = tuple(endpoint + collar.H for endpoint in selected_root)
    selected_m2 = tuple(
        endpoint + 2 * collar.H for endpoint in selected_root
    )
    require(selected_m2 == SELECTED_I, "selected M2 ceased to be I")
    unit = collar.copies.T // 13
    allocation_atoms = {}
    for q in (0, 3, 7, 11):
        left = (SELECTED_I[0] + q * unit) % collar.copies.T
        allocation_atoms[q] = (
            left,
            left + collar.copies.LENGTH,
        )
    horn_piece_intervals = {
        (record["root_left"] + step * collar.H,
         record["root_left"] + step * collar.H + collar.copies.LENGTH)
        for cell in HORN_20
        for record in cell_records[cell]
        for step in range(record["path_length"] + 1)
    }
    allocation_membership = {
        q: interval in horn_piece_intervals
        for q, interval in allocation_atoms.items()
    }
    require(
        allocation_membership == {0: True, 3: False, 7: False, 11: False},
        "allocation-atom/collar membership boundary changed",
    )

    print("Q3/Q11 COLLAR CELL/PATH SCOUT")
    print(
        f"global_cells={len(nonempty_cells)};global_roots={len(records)};"
        f"common42_cells={len(COMMON_42)};horn20_cells={len(HORN_20)}"
    )
    for name, rows in (
        ("common42", bank),
        ("horn20", horn),
        ("extra22", extra),
        ("outside", outside),
    ):
        summary = summarize(rows)
        print(
            f"{name}=roots:{summary['roots']};"
            f"root_multiplicity:{tuple(sorted(summary['root_multiplicity'].items()))};"
            f"semantic:{tuple(sorted(summary['semantic'].items()))};"
            f"survive14h:{tuple(sorted(summary['survive_14h'].items()))}"
        )
        print(
            f"{name}_path_lengths="
            f"{tuple(sorted(summary['path_lengths'].items()))}"
        )
        print(
            f"{name}_root_mod13="
            f"{tuple(sorted(summary['root_mod_13'].items()))}"
        )
        if "root_holes" in summary:
            print(
                f"{name}_root_holes="
                f"{tuple(sorted(summary['root_holes'].items()))};"
                f"M1_all_six={tuple(sorted(summary['m1_all_six'].items()))};"
                f"M2_all_six={tuple(sorted(summary['m2_all_six'].items()))};"
                "R_shifted_all_six="
                f"{tuple(sorted(summary['root_shifted_all_six'].items()))}"
            )
    print(f"cell_types={tuple(sorted(cell_types.items()))}")
    print(
        f"root_index_profiles={tuple(sorted(root_index_profiles.items()))}"
    )
    print(
        f"horn20_common_atoms={sum(cell_common_counts[cell] for cell in HORN_20)};"
        f"common42_selected_I_roots={len(distinguished_42)};"
        f"distinguished_labelled_roots={len(distinguished)};"
        f"distinguished_physical_roots={len(distinguished_physical_roots)};"
        f"R={selected_root};M1={selected_m1};M2={selected_m2}"
    )
    print(
        f"distinguished_root_indices="
        f"{tuple(sorted(Counter(row['root_index'] for row in distinguished).items()))};"
        f"path_lengths="
        f"{tuple(sorted(Counter(row['path_length'] for row in distinguished).items()))};"
        f"root_holes="
        f"{tuple(sorted(Counter(tuple(collar.copies.FACTOR_NAMES[index] for index, value in enumerate(row['factor_rows'][0]) if not value) for row in distinguished).items()))}"
    )
    print(
        "distinguished_cell_lengths="
        f"{tuple((row['cell'], row['path_length']) for row in distinguished)}"
    )
    print(
        f"pulled_allocation_atom_membership={allocation_membership};"
        "q0_is_selected_M2=1"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
