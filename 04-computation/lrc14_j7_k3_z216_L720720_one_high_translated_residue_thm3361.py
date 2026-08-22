#!/usr/bin/env python3
"""Exact companion for the reserved THM-3361 L720720 family.

The universe is exactly atlas rows 191, 228, and 332 in the projected
``k=3, z_1=216`` screen.  The certificate reconstructs the inherited exact
ray/common-status partition, proves that every possible completion has
exactly one high drift, and closes all scalar-admissible one-high cases by a
translation-uniform safe-residue cardinality gate.  A separately implemented
located-torsion search confirms every terminal with effective order two or
three.

This is a necessary projected-coordinate certificate only.  It proves no
physical entry, arbitrary-k statement, rung, or LRC(14) conclusion.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.util
import multiprocessing as mp
from collections import Counter
from fractions import Fraction
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INHERITED_SOURCE = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z216_sixteen_family_prefix_translated_two_high_"
    "closure_scout_20260812.py"
)
INHERITED_OUTPUT = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z216_sixteen_family_prefix_translated_two_high_"
    "closure_scout_20260812.out"
)
FOURTH_SOURCE = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z216_fourth_intrinsic_ruler_cost_prefix_"
    "affine_multicover_closure_scout_20260803.py"
)
FOURTH_OUTPUT = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z216_fourth_intrinsic_ruler_cost_prefix_"
    "affine_multicover_closure_scout_20260803.out"
)
AUDIT_SOURCE = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z216_cost_prefix_independent_audit_20260803.py"
)
NATURAL_SOURCE = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z216_natural_wall_family_screen_descent_thm3281.py"
)
TORSION_SOURCE = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z297_ray_status_torsion_closure_thm2941.py"
)

DEPENDENCIES = (
    (
        INHERITED_SOURCE,
        "cfb020bfc6636a52f1eaf55f82a925e70c11c90da7f87f36b0bd77ece1ec6a62",
        None,
    ),
    (
        INHERITED_OUTPUT,
        "a88646fbd28d807a0cc9671c509c4424056a539b49d04a2076ba17de57ef5ee4",
        "semantic_sha256="
        "3f332da31cc80395396aae575ae8aaffe48e82513a7a3fcf56d72b046ff0b63d",
    ),
    (
        FOURTH_SOURCE,
        "b515c70174d58ad859a08c29949cdee36a4e04122451a66237732570ab5ee213",
        None,
    ),
    (
        FOURTH_OUTPUT,
        "d1845611f27d427a1d38afe349ed07bc964590fd39a3e88cd45e6ea34a86bc38",
        "semantic_sha256="
        "c201276eb84f71806a7eb683e42d723b930f7cd0f79862e8d6b1e0ad07d37dd9",
    ),
    (
        AUDIT_SOURCE,
        "d3581bc937597270558a77c0a02398a4cfd368306f5a302be01240cdb4cef7c3",
        None,
    ),
    (
        NATURAL_SOURCE,
        "430dee7ba03e0d5c9ae0df72ac512500de4f7056cb4663d1c8468bfb93a49bfe",
        None,
    ),
    (
        TORSION_SOURCE,
        "d062c7ac8ebf6a433c8fb1543293e941c85625e2eb40b82fcf05fc2404539b0a",
        None,
    ),
)

LEVEL = 216
RULER = 720720
HIGH_FLOOR = 70981
DIVISOR_GCD = 72
EXPECTED_ATLAS = (480, 447, 33)
EXPECTED_SCREEN_TOTALS = (16900, 8947, 7234, 719)
EXPECTED_FARKAS_TOTALS = (0, 7234)
EXPECTED_ONE_HIGH = (820, 719, 54)
EXPECTED_LEDGER = (372917, 372914)
EXPECTED_WALL = (113, 110)
EXPECTED_FAMILIES = (13, 12)

EXPECTED_ROWS = {
    191: {
        "atlas": ((1, 5, 8, 9, 11, 13), RULER, HIGH_FLOOR, True),
        "components": 22,
        "screen": (8288, 4002, 3899, 387, 0, 3899),
        "first_d": 10010,
        "residual_sha256": (
            "ef31239f5cbe2f903ca96174d109fa24ea857843d7b50db8d92c40904975710a"
        ),
        "screen_sha256": (
            "10cf9f55584ffa83d0b9766f083168fdd10ad1ba4c1cca389a84d0dc8dd8acad"
        ),
        "sign_census": ((-1, 31704), (1, 31617)),
        "recurrence_checks": 644905,
        "two_gap": Fraction(179568616148, 133555340093175),
        "zero_high": 369,
        "zero_sha256": (
            "f351ea0876ccebe59064e6b4737e78fd874e4bc4414675fb9fa19e68173edcc0"
        ),
        "one_high": 462,
        "one_passports": 387,
        "low_pairs": 26,
        "one_sha256": (
            "143bc7c679e1f72dbe865a094de52b1c898e99ccfb968d120d88851d7c797066"
        ),
        "qualifying": ((2, 409), (3, 51), (6, 2)),
        "effective": ((2, 409), (3, 53)),
        "torsion_sha256": (
            "7a67dbfd35b39ab489f9e96f4deac01cfea2831228e67df5c6a4fab6db95972c"
        ),
        "cell_range": (124168, 130496),
        "minimum_translated": (
            1,
            (2, 9009, 10010, 13104),
            2,
            (275, 320),
            2,
        ),
        "minimum_torsion": (
            1,
            (2, 9009, 10010, 13104),
            2,
            (275, 320),
            2,
            2,
            2,
            1,
        ),
        "classification_sha256": (
            "076480013a03db4cdc8de12d9ceef558e3014effd223a3207f625f6d8be2e89d"
        ),
    },
    228: {
        "atlas": ((1, 8, 9, 10, 11, 13), RULER, HIGH_FLOOR, True),
        "components": 26,
        "screen": (8135, 4706, 3109, 320, 0, 3109),
        "first_d": 10010,
        "residual_sha256": (
            "3b45b650b46b5703ce619f7de0a70266a9d65abe797ed48b59b92d937c1310a6"
        ),
        "screen_sha256": (
            "eb7d0e1f46f6f534dd0fd909763436be54f1280f4b5e01a9c707dc655de5705e"
        ),
        "sign_census": ((-1, 11589), (1, 11506)),
        "recurrence_checks": 235223,
        "two_gap": Fraction(445315966867, 211281696826200),
        "zero_high": 301,
        "zero_sha256": (
            "be77b8c8ec1165c4a110aa8f60b6145a163ee3165132579151e359a9199edbb5"
        ),
        "one_high": 346,
        "one_passports": 320,
        "low_pairs": 26,
        "one_sha256": (
            "375d1c6ce3d98da8284e53feb29485880232d4ee12c141c88fede5f0d7ef6535"
        ),
        "qualifying": ((2, 331), (3, 15)),
        "effective": ((2, 331), (3, 15)),
        "torsion_sha256": (
            "7232302a05cb2a3a38eafc97817e10f34dfd075022a1078a78ae4f0d928d751d"
        ),
        "cell_range": (130334, 139040),
        "minimum_translated": (
            1,
            (2, 7280, 10010, 90090),
            2,
            (297, 328),
            2,
        ),
        "minimum_torsion": (
            1,
            (2, 7280, 10010, 90090),
            2,
            (297, 328),
            2,
            2,
            2,
            1,
        ),
        "classification_sha256": (
            "984c82d44a81bf879c15a3023befdfc80486b73ba9af9d38264ba59eac0cf250"
        ),
    },
    332: {
        "atlas": ((2, 5, 8, 9, 11, 13), RULER, HIGH_FLOOR, True),
        "components": 22,
        "screen": (477, 239, 226, 12, 0, 226),
        "first_d": 10010,
        "residual_sha256": (
            "c15ebeae8c13dae20349a540ce9d08f26c7447411b588afada331cdccae5de82"
        ),
        "screen_sha256": (
            "c38e5ef4abf9408f4be8b841ead9296dd9da5974c60ee9c1cbdd1ba6e4ef36e0"
        ),
        "sign_census": ((-1, 3464), (1, 3413)),
        "recurrence_checks": 70016,
        "two_gap": Fraction(8767934779397, 3096873978798600),
        "zero_high": 12,
        "zero_sha256": (
            "b775adabb060b8a2967f20f9360bab00798be71c550b408f203db5d9fe2dfbc8"
        ),
        "one_high": 12,
        "one_passports": 12,
        "low_pairs": 2,
        "one_sha256": (
            "a225efc2453c408d7cb581bfc96201fc533f58a89bdf25a3558223387b857ca3"
        ),
        "qualifying": ((2, 12),),
        "effective": ((2, 12),),
        "torsion_sha256": (
            "62bbbfccc3ab4c0445e93416e23b0183c3e7a8e7a1721650b0e8eaf348fd6a0b"
        ),
        "cell_range": (136112, 138346),
        "minimum_translated": (
            452,
            (528, 10010, 13104, 80080),
            528,
            (243, 275),
            528,
        ),
        "minimum_torsion": (
            264,
            (528, 10010, 13104, 80080),
            528,
            (243, 275),
            2,
            2,
            528,
            264,
        ),
        "classification_sha256": (
            "25e6cfdf8dc8d2716a73bfc00b79d42a2bfe75369bd4b0398320b88e638e0985"
        ),
    },
}

# Filled after the first complete exact replay, then enforced in every run.
EXPECTED_NEW_ROW_SHA256 = {
    191: "155fc937d5217ab68b70b7d25c799d071e16f201a5a43d4fd8f45be6d1cfc290",
    228: "4ffee124df2e2a6cf0979808bf970da3b1041605b81631b1d7278dfb8c19a7e7",
    332: "fbfebc866707217b779a97a41527e645a964c9eb3b1cfa1278f86b5e0cd75d65",
}
EXPECTED_SEMANTIC_SHA256 = (
    "927e52562c64c35833428ea735829ba236b7dcaa5f9e4d2db24f750e1a2db77c"
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_sha(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n")
    require(b"\r" not in payload, (path, "bare CR"))
    return hashlib.sha256(payload).hexdigest()


def digest(value):
    return hashlib.sha256(repr(value).encode("ascii")).hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def direct_cell_clean(cell, label, ruler):
    residue = (label * cell) % ruler
    return 14 * residue >= ruler and 14 * (residue + label) <= 13 * ruler


def direct_order_two_or_three(cell_for_residue, denominator):
    """Find a concrete safe-cell pair of exact order two or three.

    This search does not call the inherited density-pigeonhole routine.  It
    tests the exact order shifts directly in the independently rebuilt safe
    residue dictionary.
    """
    residues = tuple(sorted(cell_for_residue))
    for order in (2, 3):
        if denominator % order:
            continue
        shift = denominator // order
        for first_residue in residues:
            second_residue = (first_residue + shift) % denominator
            if second_residue not in cell_for_residue:
                continue
            effective = denominator // gcd(denominator, shift)
            require(effective == order, (denominator, shift, effective, order))
            return (
                order,
                cell_for_residue[first_residue],
                cell_for_residue[second_residue],
                first_residue,
                second_residue,
                shift,
            )
    return (None, None, None, None, None, None)


def translation_equality_hostile():
    """An aligned-safe support covered by one translated length-d/7 band."""
    denominator = 14
    support = (2, 3)
    capacity = (denominator + 6) // 7
    aligned_danger = tuple(
        residue
        for residue in support
        if 14 * min(residue, denominator - residue) < denominator
    )
    # The open interval (3/2, 7/2) has length 2=d/7 and contains 2,3.
    translated_danger = tuple(
        residue for residue in support if 3 < 2 * residue < 7
    )
    require(not aligned_danger, aligned_danger)
    require(translated_danger == support, translated_danger)
    require(len(support) == capacity, (support, capacity))
    return (
        denominator,
        support,
        capacity,
        aligned_danger,
        (3, 7, 2),
        translated_danger,
    )


def analyze_row(index):
    expected = EXPECTED_ROWS[index]
    inherited = load(f"thm3361_inherited_{index}", INHERITED_SOURCE)
    require(inherited.LEVEL == LEVEL, (index, inherited.LEVEL))
    require(inherited.AUDIT_SOURCE == AUDIT_SOURCE, (index, "audit path"))
    require(inherited.NATURAL_SOURCE == NATURAL_SOURCE, (index, "natural path"))
    require(inherited.TORSION_SOURCE == TORSION_SOURCE, (index, "torsion path"))

    natural = load(f"thm3361_atlas_{index}", NATURAL_SOURCE)
    rows, components = natural.atlas_rows()
    atlas_census = (
        len(rows),
        sum(row[3] for row in rows),
        sum(not row[3] for row in rows),
    )
    require(atlas_census == EXPECTED_ATLAS, (index, atlas_census))
    require(rows[index] == expected["atlas"], (index, rows[index]))
    require(components[index] == expected["components"], (index, components[index]))
    body, ruler, high_floor, wall = rows[index]
    require(
        (ruler, high_floor, wall, gcd(LEVEL, ruler))
        == (RULER, HIGH_FLOOR, True, DIVISOR_GCD),
        (index, rows[index]),
    )

    returned_index, screen, direct, legacy = inherited.screen_worker(
        (index, rows[index])
    )
    screen = tuple(screen)
    require(returned_index == index, (index, returned_index))
    first_d = screen[4]
    residual = tuple(screen[13])
    screen_counts = (
        screen[9],
        screen[10],
        screen[11],
        screen[12],
        direct,
        legacy,
    )
    require(screen_counts == expected["screen"], (index, screen_counts))
    require(
        screen[:6] == (LEVEL, body, ruler, high_floor, first_d, wall),
        (index, screen[:6]),
    )
    require(first_d == expected["first_d"], (index, first_d))
    require(screen[12] == len(residual), (index, "residual length"))
    require(direct + legacy == screen[11], (index, "Farkas partition"))
    require(screen[16] == screen[11], (index, "status replay"))
    require(
        all(ds.count(first_d) == 1 and lcm(*ds) == RULER for ds in residual),
        (index, "passport typing"),
    )
    screen_record = (index, tuple(screen[:19]), direct, legacy)
    require(digest(residual) == expected["residual_sha256"], (index, "residual hash"))
    require(digest(screen_record) == expected["screen_sha256"], (index, "screen hash"))

    suffixes = tuple(tuple(d for d in ds if d != first_d) for ds in residual)
    denominator_frequency = tuple(
        sorted(Counter(d for ds in suffixes for d in ds).items())
    )
    repeat_shapes = tuple(
        sorted(
            Counter(
                tuple(sorted(Counter(ds).values())) for ds in suffixes
            ).items()
        )
    )
    modulus_frequency = tuple(sorted(Counter(lcm(*ds) for ds in residual).items()))

    audit = inherited.load(f"thm3361_audit_{index}", inherited.AUDIT_SOURCE)
    terminal_natural = audit.load(
        f"thm3361_terminal_natural_{index}", audit.NATURAL_SOURCE
    )
    engine = audit.status_engine(terminal_natural, f"thm3361_terminal_{index}")
    engine.FIRST = LEVEL
    engine.ray.FIRST = LEVEL
    stream = engine.ray.Stream(body)
    require(
        (stream.L, stream.high_floor, stream.first_d)
        == (ruler, high_floor, first_d),
        (index, "stream atlas"),
    )
    needed = {
        d
        for ds in residual
        for d in engine.suffix_slots(ds, stream.first_d)
    }
    low, high, sign_census, recurrence_checks = engine.build_literal_tables(
        stream, needed
    )
    require(sign_census == expected["sign_census"], (index, sign_census))
    require(
        recurrence_checks == expected["recurrence_checks"],
        (index, recurrence_checks),
    )

    two_gap, two_gap_witness = engine.duplicate_two_high_gap(
        stream, residual, low, high
    )
    require(two_gap == expected["two_gap"] > 0, (index, two_gap))
    zero_high = engine.zero_high_scalar_passes(stream, residual, low)
    require(len(zero_high) == expected["zero_high"], (index, len(zero_high)))
    require(digest(zero_high) == expected["zero_sha256"], (index, "zero hash"))

    one_high = engine.one_high_cases(stream, residual, low, high)
    one_high_passports = tuple(sorted({case[0] for case in one_high}))
    low_pair_count = len(
        {
            tuple(sorted(record[1] for record in low_records))
            for _ds, _high_d, low_records, _excess in one_high
        }
    )
    require(len(one_high) == expected["one_high"], (index, len(one_high)))
    require(
        len(one_high_passports) == expected["one_passports"],
        (index, len(one_high_passports)),
    )
    require(one_high_passports == tuple(sorted(residual)), (index, "passport coverage"))
    require(low_pair_count == expected["low_pairs"], (index, low_pair_count))
    require(digest(one_high) == expected["one_sha256"], (index, "one-high hash"))
    high_denominator_histogram = tuple(
        sorted(Counter(case[1] for case in one_high).items())
    )

    cell_cache = {}
    clean_audit_cache = set()
    inherited_witness_cache = {}
    translated_packet = []
    direct_packet = []
    qualifying = Counter()
    effective = Counter()
    direct_orders = Counter()
    minimum_torsion = None
    minimum_translated = None
    for ds, high_d, low_records, excess in one_high:
        labels = tuple(sorted(record[1] for record in low_records))
        if labels not in cell_cache:
            cell_cache[labels] = engine.fixed_safe_cells(stream, labels)
        cells = cell_cache[labels]
        if labels not in clean_audit_cache:
            require(
                all(
                    direct_cell_clean(cell, label, ruler)
                    for cell in cells
                    for label in (stream.first, *labels)
                ),
                (index, labels, "direct clean reconstruction"),
            )
            clean_audit_cache.add(labels)

        cell_for_residue = {}
        for cell in cells:
            cell_for_residue.setdefault(cell % high_d, cell)
        support = tuple(sorted(cell_for_residue))
        danger_capacity = (high_d + 6) // 7
        translated_margin = len(support) - danger_capacity
        require(translated_margin > 0, (index, ds, high_d, labels, translated_margin))
        translated_record = (
            translated_margin,
            ds,
            high_d,
            labels,
            len(support),
        )
        if minimum_translated is None or translated_record < minimum_translated:
            minimum_translated = translated_record
        translated_packet.append(
            (
                ds,
                high_d,
                labels,
                excess,
                len(cells),
                len(support),
                danger_capacity,
                translated_margin,
            )
        )

        key = (labels, high_d)
        if key not in inherited_witness_cache:
            inherited_witness_cache[key] = engine.torsion_pigeonhole(cells, high_d)
        inherited_witness = inherited_witness_cache[key]
        require(
            inherited_witness[0] is not None and inherited_witness[1] in (2, 3),
            (index, key, inherited_witness),
        )
        require(
            (inherited_witness[2], inherited_witness[10])
            == (len(support), len(cells)),
            (index, key, "inherited support"),
        )
        qualifying[inherited_witness[0]] += 1
        effective[inherited_witness[1]] += 1
        torsion_record = (
            inherited_witness[2] - inherited_witness[3],
            ds,
            high_d,
            labels,
            inherited_witness[0],
            inherited_witness[1],
            inherited_witness[2],
            inherited_witness[3],
        )
        if minimum_torsion is None or torsion_record < minimum_torsion:
            minimum_torsion = torsion_record

        direct_witness = direct_order_two_or_three(cell_for_residue, high_d)
        require(direct_witness[0] in (2, 3), (index, key, direct_witness))
        order, first_cell, second_cell, first_residue, second_residue, shift = (
            direct_witness
        )
        require(
            (first_cell % high_d, second_cell % high_d)
            == (first_residue, second_residue),
            (index, key, "direct residues"),
        )
        require(
            (second_residue - first_residue) % high_d == shift,
            (index, key, "direct shift"),
        )
        require(
            all(
                direct_cell_clean(cell, label, ruler)
                for cell in (first_cell, second_cell)
                for label in (stream.first, *labels)
            ),
            (index, key, "located cells not clean"),
        )
        require(7 >= order, (index, key, "phase separation below one seventh"))
        direct_orders[order] += 1
        direct_packet.append(
            (
                ds,
                high_d,
                labels,
                len(cells),
                len(support),
                direct_witness,
            )
        )

    inherited_torsion_packet = tuple(sorted(inherited_witness_cache.items()))
    translated_packet = tuple(translated_packet)
    direct_packet = tuple(direct_packet)
    require(len(inherited_witness_cache) == len(one_high), (index, "torsion keys"))
    require(tuple(sorted(qualifying.items())) == expected["qualifying"], (index, qualifying))
    require(tuple(sorted(effective.items())) == expected["effective"], (index, effective))
    require(
        digest(inherited_torsion_packet) == expected["torsion_sha256"],
        (index, "inherited torsion hash"),
    )
    require(
        (
            min(map(len, cell_cache.values()), default=0),
            max(map(len, cell_cache.values()), default=0),
        )
        == expected["cell_range"],
        (index, "cell range"),
    )
    require(minimum_translated == expected["minimum_translated"], (index, minimum_translated))
    require(minimum_torsion == expected["minimum_torsion"], (index, minimum_torsion))

    classification_hash = digest(
        (
            screen_record,
            residual,
            denominator_frequency,
            repeat_shapes,
            modulus_frequency,
            sign_census,
            recurrence_checks,
            two_gap,
            two_gap_witness,
            zero_high,
            one_high,
            one_high_passports,
            high_denominator_histogram,
            inherited_torsion_packet,
            (),
            (),
            minimum_torsion,
            minimum_translated,
            (),
        )
    )
    require(
        classification_hash == expected["classification_sha256"],
        (index, classification_hash),
    )

    translated_hash = digest(translated_packet)
    direct_hash = digest(direct_packet)
    row_packet = (
        index,
        body,
        components[index],
        screen_counts,
        digest(screen_record),
        digest(residual),
        two_gap,
        digest(two_gap_witness),
        len(zero_high),
        digest(zero_high),
        len(one_high),
        len(one_high_passports),
        low_pair_count,
        digest(one_high),
        minimum_translated,
        translated_hash,
        tuple(sorted(qualifying.items())),
        tuple(sorted(effective.items())),
        tuple(sorted(direct_orders.items())),
        minimum_torsion,
        digest(inherited_torsion_packet),
        direct_hash,
        classification_hash,
    )
    row_hash = digest(row_packet)
    expected_new = EXPECTED_NEW_ROW_SHA256[index]
    if expected_new is not None:
        require(row_hash == expected_new, (index, row_hash))
    return {
        "index": index,
        "body": body,
        "components": components[index],
        "screen": screen_counts,
        "screen_sha256": digest(screen_record),
        "residual_sha256": digest(residual),
        "two_gap": two_gap,
        "two_gap_witness_sha256": digest(two_gap_witness),
        "zero_high": len(zero_high),
        "zero_sha256": digest(zero_high),
        "one_high": len(one_high),
        "one_passports": len(one_high_passports),
        "low_pairs": low_pair_count,
        "one_sha256": digest(one_high),
        "torsion_keys": len(inherited_witness_cache),
        "qualifying": tuple(sorted(qualifying.items())),
        "effective": tuple(sorted(effective.items())),
        "direct_orders": tuple(sorted(direct_orders.items())),
        "minimum_translated": minimum_translated,
        "translated_sha256": translated_hash,
        "minimum_torsion": minimum_torsion,
        "inherited_torsion_sha256": digest(inherited_torsion_packet),
        "direct_torsion_sha256": direct_hash,
        "cell_range": expected["cell_range"],
        "classification_sha256": classification_hash,
        "row_sha256": row_hash,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=3)
    args = parser.parse_args()
    require(1 <= args.processes <= len(EXPECTED_ROWS), args.processes)

    for path, expected_hash, needle in DEPENDENCIES:
        require(lf_sha(path) == expected_hash, (path, "dependency changed"))
        if needle is not None:
            require(needle in path.read_text(encoding="utf-8"), (path, needle))

    syntax = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    assertion_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(syntax))
    float_nodes = sum(
        isinstance(node, ast.Constant) and isinstance(node.value, float)
        for node in ast.walk(syntax)
    )
    require((assertion_nodes, float_nodes) == (0, 0), "truth-gate syntax")
    require(HIGH_FLOOR == 13 * RULER // 132 + 1, "high-floor formula")
    require(LEVEL < HIGH_FLOOR, "inherited projected high gate inactive")
    hostile = translation_equality_hostile()

    indices = tuple(EXPECTED_ROWS)
    if args.processes == 1:
        results = tuple(analyze_row(index) for index in indices)
    else:
        with mp.get_context("spawn").Pool(args.processes) as pool:
            results = tuple(pool.imap(analyze_row, indices, chunksize=1))
    results = tuple(sorted(results, key=lambda row: row["index"]))
    require(tuple(row["index"] for row in results) == indices, "row order")

    screen_totals = tuple(
        sum(row["screen"][position] for row in results) for position in range(4)
    )
    farkas_totals = (
        sum(row["screen"][4] for row in results),
        sum(row["screen"][5] for row in results),
    )
    one_high_totals = (
        sum(row["one_high"] for row in results),
        sum(row["one_passports"] for row in results),
        sum(row["low_pairs"] for row in results),
    )
    require(screen_totals == EXPECTED_SCREEN_TOTALS, screen_totals)
    require(farkas_totals == EXPECTED_FARKAS_TOTALS, farkas_totals)
    require(one_high_totals == EXPECTED_ONE_HIGH, one_high_totals)
    require(
        screen_totals[0] == sum(screen_totals[1:]),
        (screen_totals, "screen partition"),
    )

    zero_high_total = sum(row["zero_high"] for row in results)
    require(zero_high_total == 682 > 0, zero_high_total)
    require(all(row["two_gap"] > 0 for row in results), "two-high gap")
    qualifying = Counter()
    effective = Counter()
    direct_orders = Counter()
    for row in results:
        qualifying.update(dict(row["qualifying"]))
        effective.update(dict(row["effective"]))
        direct_orders.update(dict(row["direct_orders"]))
    require(
        tuple(sorted(qualifying.items())) == ((2, 752), (3, 66), (6, 2)),
        qualifying,
    )
    require(
        tuple(sorted(effective.items())) == ((2, 752), (3, 68)),
        effective,
    )
    require(
        sum(direct_orders.values()) == EXPECTED_ONE_HIGH[0]
        and set(direct_orders).issubset({2, 3}),
        direct_orders,
    )
    translated_minimum = min(row["minimum_translated"] for row in results)
    torsion_minimum = min(row["minimum_torsion"] for row in results)
    require(translated_minimum[0] == 1, translated_minimum)
    require(torsion_minimum[0] == 1, torsion_minimum)

    semantic_rows = tuple(
        (
            row["index"],
            row["body"],
            row["components"],
            row["screen"],
            row["screen_sha256"],
            row["residual_sha256"],
            row["two_gap"],
            row["two_gap_witness_sha256"],
            row["zero_high"],
            row["zero_sha256"],
            row["one_high"],
            row["one_passports"],
            row["low_pairs"],
            row["one_sha256"],
            row["torsion_keys"],
            row["qualifying"],
            row["effective"],
            row["direct_orders"],
            row["minimum_translated"],
            row["translated_sha256"],
            row["minimum_torsion"],
            row["inherited_torsion_sha256"],
            row["direct_torsion_sha256"],
            row["cell_range"],
            row["classification_sha256"],
            row["row_sha256"],
        )
        for row in results
    )
    semantic_packet = (
        tuple((path.name, expected_hash) for path, expected_hash, _needle in DEPENDENCIES),
        LEVEL,
        RULER,
        HIGH_FLOOR,
        DIVISOR_GCD,
        tuple(
            (index, EXPECTED_ROWS[index]["atlas"], EXPECTED_ROWS[index]["components"])
            for index in indices
        ),
        semantic_rows,
        screen_totals,
        farkas_totals,
        one_high_totals,
        zero_high_total,
        tuple(sorted(qualifying.items())),
        tuple(sorted(effective.items())),
        tuple(sorted(direct_orders.items())),
        translated_minimum,
        torsion_minimum,
        hostile,
        EXPECTED_LEDGER,
        EXPECTED_WALL,
        EXPECTED_FAMILIES,
    )
    semantic_hash = digest(semantic_packet)
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, semantic_hash)

    lines = [
        "LRC14 projected-k3 z216 gcd72/L720720 one-high translated-residue certificate",
        "status=FINITE-EXACT;theorem=THM-3361;"
        "scope=necessary_projected_k3_z1_216_family_only",
        "dependency="
        + ";".join(
            f"{path.name}:{expected_hash}"
            for path, expected_hash, _needle in DEPENDENCIES
        ),
        f"source_ast=assert_nodes:{assertion_nodes};float_literals:{float_nodes}",
        f"universe=atlas:{EXPECTED_ATLAS};family:gcd{DIVISOR_GCD}/L{RULER};"
        f"indices:{','.join(map(str, indices))};components:"
        f"{','.join(str(row['components']) for row in results)};"
        f"intrinsic_cost:{sum(row['components'] * RULER for row in results)};"
        f"high_floor:{HIGH_FLOOR};projected_first:{LEVEL}",
    ]
    for row in results:
        lines.append(
            f"ROW;index={row['index']};E={','.join(map(str, row['body']))};"
            f"components={row['components']};states={row['screen'][0]};"
            f"crude={row['screen'][1]};status={row['screen'][2]};"
            f"residual={row['screen'][3]};direct={row['screen'][4]};"
            f"legacy={row['screen'][5]};zero_high_hostiles={row['zero_high']};"
            f"one_high_cases={row['one_high']};one_high_passports={row['one_passports']};"
            f"low_pairs={row['low_pairs']};two_high_gap={row['two_gap']};"
            f"translated_min={row['minimum_translated']};"
            f"torsion_min={row['minimum_torsion']};"
            f"effective_orders={row['effective']};direct_orders={row['direct_orders']};"
            f"screen_sha256={row['screen_sha256']};"
            f"residual_sha256={row['residual_sha256']};"
            f"one_sha256={row['one_sha256']};"
            f"translated_sha256={row['translated_sha256']};"
            f"inherited_torsion_sha256={row['inherited_torsion_sha256']};"
            f"direct_torsion_sha256={row['direct_torsion_sha256']};"
            f"classification_sha256={row['classification_sha256']};"
            f"row_sha256={row['row_sha256']}"
        )
    lines.extend(
        [
            f"SCREEN_TOTAL;states={screen_totals[0]};crude={screen_totals[1]};"
            f"status={screen_totals[2]};residual={screen_totals[3]};"
            f"identity={screen_totals[0]}={screen_totals[1]}+"
            f"{screen_totals[2]}+{screen_totals[3]};"
            f"direct={farkas_totals[0]};legacy={farkas_totals[1]}",
            f"EXACTLY_ONE_HIGH;inherited_gate={LEVEL}<{HIGH_FLOOR};"
            f"zero_high_scalar_hostiles={zero_high_total};"
            f"duplicate_permitting_two_high_gaps=all_strict_positive;"
            f"one_high_cases={one_high_totals[0]};"
            f"passports={one_high_totals[1]};low_pairs={one_high_totals[2]}",
            f"TRANSLATED_CARDINALITY;cases={one_high_totals[0]};failures=0;"
            f"gate=|C_mod_d|>ceil(d/7);minimum={translated_minimum};"
            "unit_action=permutation;danger_band=arbitrary_translation",
            f"LOCATED_TORSION_CROSSCHECK;cases={one_high_totals[0]};failures=0;"
            f"inherited_qualifying={tuple(sorted(qualifying.items()))};"
            f"inherited_effective={tuple(sorted(effective.items()))};"
            f"independent_direct_orders={tuple(sorted(direct_orders.items()))};"
            f"minimum={torsion_minimum}",
            f"HOSTILE_TRANSLATION;d={hostile[0]};support={hostile[1]};"
            f"ceil_d_over_7={hostile[2]};aligned_danger={hostile[3]};"
            f"translated_open_interval_doubled={hostile[4]};"
            f"translated_danger={hostile[5]};"
            "verdict=equality_or_aligned_only_is_insufficient",
            f"consequence=family_excluded;ledger:{EXPECTED_LEDGER[0]}->"
            f"{EXPECTED_LEDGER[1]};wall:{EXPECTED_WALL[0]}->{EXPECTED_WALL[1]};"
            f"families:{EXPECTED_FAMILIES[0]}->{EXPECTED_FAMILIES[1]};"
            "projected_k3_cap=216_unchanged",
            "scope=no_physical_entry_or_arbitrary_k_or_rung_or_LRC14_conclusion",
            f"semantic_sha256={semantic_hash}",
            "all_exact_controls=PASS",
        ]
    )
    print("\n".join(lines))


if __name__ == "__main__":
    main()
