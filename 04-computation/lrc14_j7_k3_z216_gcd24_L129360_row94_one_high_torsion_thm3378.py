#!/usr/bin/env python3
"""Scratch exact audit of the first post-THM-3361 z1=216 wall row.

This file deliberately lives outside the repository.  It imports the frozen
current screen and located-torsion engines, reconstructs the live intrinsic-
cost queue, and checks only atlas row 94.  Its conclusion is confined to the
necessary projected k=3,z1=216 screen.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.util
from collections import Counter
from fractions import Fraction
from math import gcd, lcm
from pathlib import Path


LEVEL = 216
ROW_INDEX = 94
EXPECTED_ROW = ((1, 3, 8, 10, 11, 14), 129360, 12741, True)
EXPECTED_COMPONENTS = 30
EXPECTED_FAMILY = (
    52778880,
    13,
    24,
    129360,
    (94, 161, 174, 215, 237, 263, 319, 354, 369, 399, 430, 443, 472),
)
EXPECTED_REMAINING_INDICES = EXPECTED_FAMILY[4][1:]
EXPECTED_SCREEN_COUNTS = (571, 292, 195, 84)
EXPECTED_RESIDUAL_SHA256 = (
    "6f1137bbe6f0541217d4c3854fd94fcc9bc5eaf106ba3510e7bf06edfae86f9f"
)
EXPECTED_TWO_GAP = Fraction(268369544, 94245446295)
EXPECTED_ZERO_HIGH = 78
EXPECTED_ONE_HIGH = 84
EXPECTED_LOW_PAIRS = 5
EXPECTED_QUALIFYING = ((2, 80), (3, 4))
EXPECTED_EFFECTIVE = ((2, 80), (3, 4))
EXPECTED_TERMINAL_SHA256 = (
    "5cd6e447deb65f698338345d2a54d45662f1add466db6ff2faae78f3970766fd"
)
EXPECTED_ONE_HIGH_SHA256 = (
    "39cd315557f1ee59791c7e5a819472575b6857d872c75b5ff8958ed7e2df5049"
)
EXPECTED_TORSION_SHA256 = (
    "f7734fd72f4bd2a837ec9dc44610a91608676735a10cfbc1eb1f0ae65730ffa4"
)


DEPENDENCIES = {
    "04-computation/lrc14_j7_k3_z216_sixteen_family_prefix_translated_two_high_closure_scout_20260812.py":
        "cfb020bfc6636a52f1eaf55f82a925e70c11c90da7f87f36b0bd77ece1ec6a62",
    "05-knowledge/results/lrc14_j7_k3_z216_sixteen_family_prefix_translated_two_high_closure_scout_20260812.out":
        "a88646fbd28d807a0cc9671c509c4424056a539b49d04a2076ba17de57ef5ee4",
    "04-computation/lrc14_j7_k3_z216_fourth_intrinsic_ruler_cost_prefix_affine_multicover_closure_scout_20260803.py":
        "b515c70174d58ad859a08c29949cdee36a4e04122451a66237732570ab5ee213",
    "05-knowledge/results/lrc14_j7_k3_z216_fourth_intrinsic_ruler_cost_prefix_affine_multicover_closure_scout_20260803.out":
        "d1845611f27d427a1d38afe349ed07bc964590fd39a3e88cd45e6ea34a86bc38",
    "04-computation/lrc14_j7_k3_z216_cost_prefix_independent_audit_20260803.py":
        "d3581bc937597270558a77c0a02398a4cfd368306f5a302be01240cdb4cef7c3",
    "04-computation/lrc14_j7_k3_z216_natural_wall_family_screen_descent_thm3281.py":
        "430dee7ba03e0d5c9ae0df72ac512500de4f7056cb4663d1c8468bfb93a49bfe",
    "04-computation/lrc14_j7_k3_z297_ray_status_torsion_closure_thm2941.py":
        "d062c7ac8ebf6a433c8fb1543293e941c85625e2eb40b82fcf05fc2404539b0a",
    "04-computation/lrc14_j7_k3_z378_ray_status_closure_thm2941.py":
        "2ef5e0639354c38b13e17e41f91acb4143c7f60973295b0e2dd0f57eb8f38db2",
    "04-computation/lrc14_j7_k3_projected_scalar_body_atlas_thm2941.py":
        "2af6d96882f336a409a8657070ed76a75c09a53b3789101b83103b051e864ded",
    "05-knowledge/results/lrc14_j7_k3_projected_scalar_body_atlas_thm2941.out":
        "cee82237ce1f51729813b9c916edd3353204c18172abe1d71278dee2c5562eda",
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_sha(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(payload).hexdigest()


def digest(value):
    return hashlib.sha256(repr(value).encode("ascii")).hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def direct_order_two_or_three(cell_for_residue, denominator):
    """Independent exact-shift search, not the inherited density routine."""
    residues = tuple(sorted(cell_for_residue))
    for order in (2, 3):
        if denominator % order:
            continue
        shift = denominator // order
        for first_residue in residues:
            second_residue = (first_residue + shift) % denominator
            if second_residue in cell_for_residue:
                return (
                    order,
                    cell_for_residue[first_residue],
                    cell_for_residue[second_residue],
                    first_residue,
                    second_residue,
                    shift,
                )
    return (None, None, None, None, None, None)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
    )
    args = parser.parse_args()
    root = args.root.resolve()

    for relative, expected in DEPENDENCIES.items():
        path = root / relative
        require(lf_sha(path) == expected, (relative, lf_sha(path), expected))

    source = root / next(iter(DEPENDENCIES))
    syntax = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(syntax)),
            "inherited wrapper contains assert")
    require(not any(isinstance(node, ast.Constant) and isinstance(node.value, float)
                    for node in ast.walk(syntax)),
            "inherited wrapper contains float literal")

    cumulative = load("row94_cumulative", source)
    rows, components, _inherited, ranked, prefix, prefix_indices, prefix_hash = (
        cumulative.reconstruct_queue()
    )
    require((len(ranked), len(prefix), len(prefix_indices)) == (29, 16, 236),
            "queue census")
    require(prefix_hash == cumulative.EXPECTED_PREFIX_SHA256, prefix_hash)
    require(ranked[16][:5] == (50450400, 3, 72, 720720, (191, 228, 332)),
            "THM-3361 boundary")
    family = ranked[17]
    require(family[:5] == EXPECTED_FAMILY, family[:5])
    require(rows[ROW_INDEX] == EXPECTED_ROW, rows[ROW_INDEX])
    require(components[ROW_INDEX] == EXPECTED_COMPONENTS, components[ROW_INDEX])

    index, screen, direct, legacy = cumulative.screen_worker(
        (ROW_INDEX, rows[ROW_INDEX])
    )
    require(index == ROW_INDEX, index)
    screen = tuple(screen)
    residual = tuple(screen[13])
    require(
        screen[:6]
        == (
            LEVEL,
            EXPECTED_ROW[0],
            EXPECTED_ROW[1],
            EXPECTED_ROW[2],
            5390,
            EXPECTED_ROW[3],
        ),
        screen[:6],
    )
    require(tuple(screen[position] for position in (9, 10, 11, 12))
            == EXPECTED_SCREEN_COUNTS, screen[9:13])
    require(screen[9] == screen[10] + screen[11] + screen[12], "screen partition")
    require(screen[16] == screen[11], "exact status replay")
    require((direct, legacy) == (0, 195), (direct, legacy))
    require(direct + legacy == screen[11], "Farkas partition")
    require(len(residual) == 84 and digest(residual) == EXPECTED_RESIDUAL_SHA256,
            "residual packet")
    first_d = screen[4]
    require(first_d == 5390, first_d)
    require(all(ds.count(first_d) == 1 and lcm(*ds) == EXPECTED_ROW[1]
                for ds in residual), "passport typing")

    audit = load("row94_audit", cumulative.AUDIT_SOURCE)
    natural = audit.load("row94_natural", audit.NATURAL_SOURCE)
    engine = audit.status_engine(natural, "row94_terminal")
    engine.FIRST = LEVEL
    engine.ray.FIRST = LEVEL
    stream = engine.ray.Stream(EXPECTED_ROW[0])
    require((stream.L, stream.high_floor, stream.first_d)
            == EXPECTED_ROW[1:3] + (first_d,), "stream row")
    require(LEVEL < stream.high_floor, "at-least-one-high gate")

    needed = {
        denominator
        for divisors in residual
        for denominator in engine.suffix_slots(divisors, stream.first_d)
    }
    low, high, sign_census, recurrence_checks = engine.build_literal_tables(
        stream, needed
    )
    two_gap, two_gap_witness = engine.duplicate_two_high_gap(
        stream, residual, low, high
    )
    zero_high = engine.zero_high_scalar_passes(stream, residual, low)
    one_high = engine.one_high_cases(stream, residual, low, high)
    require(two_gap == EXPECTED_TWO_GAP > 0, two_gap)
    require(len(zero_high) == EXPECTED_ZERO_HIGH, len(zero_high))
    require(len(one_high) == EXPECTED_ONE_HIGH, len(one_high))
    require(digest(one_high) == EXPECTED_ONE_HIGH_SHA256, digest(one_high))

    passport_counts = Counter(case[0] for case in one_high)
    require(tuple(sorted(passport_counts)) == tuple(sorted(residual)),
            "one-high passport coverage")
    require(set(passport_counts.values()) == {1}, passport_counts)
    low_pairs = {
        tuple(sorted(record[1] for record in low_records))
        for _ds, _high_d, low_records, _excess in one_high
    }
    require(len(low_pairs) == EXPECTED_LOW_PAIRS, low_pairs)

    cell_cache = {}
    inherited_cache = {}
    inherited_qualifying = Counter()
    inherited_effective = Counter()
    direct_orders = Counter()
    direct_packet = []
    for divisors, high_d, low_records, excess in one_high:
        labels = tuple(sorted(record[1] for record in low_records))
        if labels not in cell_cache:
            cell_cache[labels] = engine.fixed_safe_cells(stream, labels)
        cells = cell_cache[labels]
        key = (labels, high_d)
        if key not in inherited_cache:
            inherited_cache[key] = engine.torsion_pigeonhole(cells, high_d)
        witness = inherited_cache[key]
        require(witness[0] is not None, (divisors, key, witness))
        inherited_qualifying[witness[0]] += 1
        inherited_effective[witness[1]] += 1

        cell_for_residue = {}
        for cell in cells:
            cell_for_residue.setdefault(cell % high_d, cell)
        direct_witness = direct_order_two_or_three(cell_for_residue, high_d)
        order, first_cell, second_cell, first_residue, second_residue, shift = (
            direct_witness
        )
        require(order in (2, 3), (divisors, key, direct_witness))
        require(high_d // gcd(high_d, shift) == order, direct_witness)
        require((second_residue - first_residue) % high_d == shift,
                direct_witness)
        require(7 >= order, order)
        require(all(
            engine.cell_clean(cell, label, stream.L)
            for cell in (first_cell, second_cell)
            for label in (stream.first, *labels)
        ), (divisors, key, "direct cells not fixed-safe"))
        direct_orders[order] += 1
        direct_packet.append(
            (divisors, high_d, labels, excess, len(cells),
             len(cell_for_residue), direct_witness)
        )

    inherited_packet = tuple(sorted(inherited_cache.items()))
    require(len(inherited_cache) == len(one_high), "unique torsion keys")
    require(tuple(sorted(inherited_qualifying.items())) == EXPECTED_QUALIFYING,
            inherited_qualifying)
    require(tuple(sorted(inherited_effective.items())) == EXPECTED_EFFECTIVE,
            inherited_effective)
    require(tuple(sorted(direct_orders.items())) == EXPECTED_EFFECTIVE,
            direct_orders)
    require(digest(inherited_packet) == EXPECTED_TORSION_SHA256,
            digest(inherited_packet))

    terminal = cumulative.terminal_worker((ROW_INDEX, EXPECTED_ROW[0], residual))
    require(digest(terminal) == EXPECTED_TERMINAL_SHA256, digest(terminal))
    require(terminal[19:22] == (0, 0, 0), terminal[19:22])

    direct_packet = tuple(direct_packet)
    screen_record = (ROW_INDEX, tuple(screen[:19]), direct, legacy)
    semantic_packet = (
        tuple(DEPENDENCIES.items()),
        family[:5],
        rows[ROW_INDEX],
        components[ROW_INDEX],
        screen_record,
        residual,
        tuple(sorted(needed)),
        sign_census,
        recurrence_checks,
        two_gap,
        two_gap_witness,
        zero_high,
        one_high,
        passport_counts,
        inherited_packet,
        direct_packet,
        terminal,
    )

    print("LRC14 projected-k3 z216 first-live-family row94 phase terminal")
    print("status=FINITE-EXACT;scope=necessary_projected_k3_z1_216_row_only")
    print("dependency=" + ";".join(
        f"{Path(path).name}:{value}" for path, value in DEPENDENCIES.items()
    ))
    print(
        "queue=post_THM3361_families:12;walls:110;"
        f"first_family:gcd24/L129360;rows:13;cost:52778880;"
        f"indices:{','.join(map(str, EXPECTED_FAMILY[4]))}"
    )
    print(
        "ROW;index=94;E=1,3,8,10,11,14;L=129360;high=12741;"
        "components=30;component_cost=3880800;"
        f"states={screen[9]};crude={screen[10]};status={screen[11]};"
        f"residual={screen[12]};direct={direct};legacy={legacy};"
        f"residual_sha256={digest(residual)}"
    )
    print(
        f"EXACTLY_ONE_HIGH;inherited_gate={LEVEL}<{stream.high_floor};"
        f"zero_high_scalar_hostiles={len(zero_high)};"
        f"duplicate_permitting_two_high_gap={two_gap};"
        f"one_high_cases={len(one_high)};passports={len(passport_counts)};"
        f"passport_multiplicity_hist={tuple(sorted(Counter(passport_counts.values()).items()))};"
        f"low_pairs={len(low_pairs)};one_sha256={digest(one_high)}"
    )
    print(
        "LOCATED_PHASE;"
        f"inherited_least_r={tuple(sorted(inherited_qualifying.items()))};"
        f"inherited_effective={tuple(sorted(inherited_effective.items()))};"
        f"independent_direct_orders={tuple(sorted(direct_orders.items()))};"
        f"torsion_sha256={digest(inherited_packet)};"
        f"direct_sha256={digest(direct_packet)}"
    )
    print(
        "mechanism=effective_order_2_or_3_gives_nonzero_phase_gap_at_least_1/3;"
        "strict_open_danger_diameter=1/7;therefore_one_of_two_fixed_safe_cells_"
        "survives_every_translated_high_phase;projected_safe_mass=1>36/91"
    )
    print(
        "consequence=row94_screen_empty;projected_ledger:372914->372913;"
        "wall_rows:110->109;family_not_excluded;"
        f"remaining_family_indices:{','.join(map(str, EXPECTED_REMAINING_INDICES))};"
        "remaining_family_cost=48898080"
    )
    print("scope=no_physical_entry_or_arbitrary_k_or_rung_or_LRC14_conclusion")
    print(f"terminal_sha256={digest(terminal)}")
    print(f"semantic_sha256={digest(semantic_packet)}")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
