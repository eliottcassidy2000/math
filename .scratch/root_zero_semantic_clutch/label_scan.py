#!/usr/bin/env python3
"""Scan common lawful full-target labels for a uniform semantic clutch unit."""

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / ".scratch/root_zero_semantic_clutch"))

import clutch_probe as probe


COMMON_S = (0, 1, 2, 3, 8, 9, 10, 11, 12)
COMMON_T = (3, 4, 5, 6, 7, 8, 9, 10, 11)
STATE = None


def initialize():
    global STATE
    module, _prefixes, _, _, rails, present, _starts = (
        probe.clutch.relative.lift.m.core.build_carrier_data()
    )
    pair_prefixes = probe.clutch.relative.private.build_pair_prefixes(module)
    fork = probe.semantic.deepest_fork(module)
    semantic_prefixes = probe.build_semantic_prefixes(module, fork)
    source_e3 = probe.semantic.exclusive_source(module, 3)
    q_source = probe.clutch.relative.semantic.frac(
        probe.clutch.relative.Z
        + probe.clutch.Fraction(7 * 6715, probe.clutch.R)
    )
    q_target = q_source + probe.clutch.Fraction(7, probe.clutch.R)
    STATE = (module, rails, present, pair_prefixes, semantic_prefixes,
             source_e3, q_source, q_target)


def sections_and_witness(label):
    s, t = label
    (module, _rails, _present, _pairs, _prefixes,
     source_e3, q_source, q_target) = STATE
    sections = probe.semantic_sections(module, source_e3, s, t)
    witness = (
        probe.semantic.strict_member(q_source, sections[1], module.T)
        and probe.semantic.strict_member(q_target, sections[1], module.T)
    )
    return sections, witness


def scan_rail_zero(label):
    module, rails, present, pairs, prefixes, *_rest = STATE
    sections, witness = sections_and_witness(label)
    source, target, _masses, _rail_pairs, _details = probe.restricted_vectors(
        module, pairs, rails, present, prefixes, sections, 0,
        equal_weight_only=True,
    )
    source_unit, source_profile = probe.primitive(
        source, 12, probe.BROADCAST_CONTENT
    )
    target_unit, target_profile = probe.primitive(
        target, 1, probe.BROADCAST_CONTENT
    )
    return (
        label, witness, source, target, source_unit, target_unit,
        source_profile[2] if source_profile else None,
        target_profile[2] if target_profile else None,
    )


def scan_all_rails(label):
    module, rails, present, pairs, prefixes, *_rest = STATE
    sections, witness = sections_and_witness(label)
    rows = []
    for rail_index in range(14):
        source, target, _masses, _rail_pairs, _details = (
            probe.restricted_vectors(
                module, pairs, rails, present, prefixes, sections,
                rail_index, equal_weight_only=True,
            )
        )
        source_unit, source_profile = probe.primitive(
            source, 12, probe.BROADCAST_CONTENT
        )
        target_unit, target_profile = probe.primitive(
            target, 1, probe.BROADCAST_CONTENT
        )
        rows.append((
            rail_index,
            tuple(index for index, value in enumerate(source) if value),
            source == target,
            source_unit,
            target_unit,
            source_profile[2] if source_profile else None,
            target_profile[2] if target_profile else None,
        ))
    return label, witness, tuple(rows)


def main():
    labels = tuple((s, t) for s in COMMON_S for t in COMMON_T)
    with ProcessPoolExecutor(max_workers=4, initializer=initialize) as pool:
        rail_zero = tuple(pool.map(scan_rail_zero, labels, chunksize=3))
        repairs = tuple(
            row[0] for row in rail_zero
            if row[1] and row[2] == row[3] and row[4] and row[5]
        )
        print(f"RAIL0_PHASE common_label_universe={len(labels)} "
              f"all_whole_witness={all(row[1] for row in rail_zero)} "
              f"repairs={len(repairs)} labels={repairs}", flush=True)
        all_rows = tuple(pool.map(scan_all_rails, repairs, chunksize=1))

    uniform = tuple(
        (label, rows)
        for label, witness, rows in all_rows
        if witness and all(
            row[2] and row[3] and row[4]
            for row in rows
        )
    )
    print(f"common_label_universe={len(labels)} all_whole_witness="
          f"{all(row[1] for row in rail_zero)}")
    print(f"rail0_unit_repairs={len(repairs)} labels={repairs}")
    print(f"uniform_all14_labels={len(uniform)} "
          f"labels={tuple(row[0] for row in uniform)}")
    if uniform:
        print(f"first_uniform_rows={uniform[0]}")
    else:
        print("first_uniform_rows=NONE")
    print("rail0_repair_rows="
          f"{tuple((row[0], tuple(i for i,v in enumerate(row[2]) if v), row[6], row[7]) for row in rail_zero if row[0] in repairs)}")


if __name__ == "__main__":
    main()
