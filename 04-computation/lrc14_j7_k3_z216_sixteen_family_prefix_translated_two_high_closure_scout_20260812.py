#!/usr/bin/env python3
"""Exact cumulative post-THM-3320 z216 prefix closure.

Reconstruct the deterministic complete-family queue after THM-3320, take its
first sixteen families (through gcd72/L35280), and replay the exact ray and
common-status screens.  Residual denominator passports are closed by the
inherited located-torsion terminal.  Five L35280 bodies have a nonpositive
duplicate-permitting two-high scalar gap; for them this file enumerates the
finite enlarged two-high scalar bank and uses translation-uniform exact
unit-orbit counts on fixed-safe cells.  Three denominator-two equality cases
use a separate projected-measure lemma.  No search horizon is used.

This is a projected-k3, z1=216 necessary-screen result only.  It is not a
physical-entry, arbitrary-k, rung, or LRC(14) theorem.
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

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
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

FOURTH_SOURCE_SHA256 = (
    "b515c70174d58ad859a08c29949cdee36a4e04122451a66237732570ab5ee213"
)
FOURTH_OUTPUT_SHA256 = (
    "d1845611f27d427a1d38afe349ed07bc964590fd39a3e88cd45e6ea34a86bc38"
)
FOURTH_SEMANTIC_SHA256 = (
    "c201276eb84f71806a7eb683e42d723b930f7cd0f79862e8d6b1e0ad07d37dd9"
)
AUDIT_SOURCE_SHA256 = (
    "d3581bc937597270558a77c0a02398a4cfd368306f5a302be01240cdb4cef7c3"
)
NATURAL_SOURCE_SHA256 = (
    "430dee7ba03e0d5c9ae0df72ac512500de4f7056cb4663d1c8468bfb93a49bfe"
)
TORSION_SOURCE_SHA256 = (
    "d062c7ac8ebf6a433c8fb1543293e941c85625e2eb40b82fcf05fc2404539b0a"
)

LEVEL = 216
PREFIX_FAMILIES = 16
EXPECTED_ATLAS = (480, 447, 33)
EXPECTED_INHERITED_PREFIX_ROWS = (5, 18, 4, 4)
EXPECTED_QUEUE_STATE = (349, 29)
EXPECTED_PREFIX_STATE = (236, 16)
EXPECTED_SCREEN_TOTALS = (227986, 109375, 118093, 518)
EXPECTED_FARKAS_TOTALS = (2, 118091)
EXPECTED_LEDGER = (373153, 372917)
EXPECTED_WALL = (349, 113)
EXPECTED_FAMILY_COUNT = (29, 13)
EXPECTED_LAST_BOUNDARY = (42124320, 50450400)
EXPECTED_PREFIX_SHA256 = (
    "7c665728bf21bf641f70d5e683fc615b1a24d3b3cdb0595462cdd17d0daa7721"
)
EXPECTED_SCREEN_SHA256 = (
    "d8b0e7e18794bea30c038193db1949585e07e38f3011f311b71570e26f635963"
)
EXPECTED_TERMINAL_SHA256 = (
    "10165a2c7a68387bff2be57501b7a2f115126a6855a89a5d9e4841bd071eea95"
)
EXPECTED_SEMANTIC_SHA256 = (
    "3f332da31cc80395396aae575ae8aaffe48e82513a7a3fcf56d72b046ff0b63d"
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


def screen_worker(item):
    index, row = item
    natural = load(f"prefix_screen_{index}", NATURAL_SOURCE)
    returned_index, screen, direct, legacy = natural.screen_worker(
        (index, (LEVEL, *row))
    )
    require(returned_index == index, (index, returned_index))
    return index, tuple(screen), direct, legacy


def common_torsion_capacity(cells, denominators):
    modulus = lcm(*denominators)
    support = len({cell % modulus for cell in cells})
    capacities = tuple(
        ((denominator + 6) // 7) * (modulus // denominator)
        for denominator in denominators
    )
    return modulus, support, capacities, support - sum(capacities)


def translated_safe_minimum(cells_tuple, denominator):
    """Minimum safe-cell count over every unit and translated danger band.

    At a fixed local coordinate, a high label of denominator ``d`` makes the
    cell residues lie in an arbitrary translated open cyclic interval of
    length ``d/7``.  Integer residues contained in such an interval have
    cyclic span at most ``floor((d-1)/7)``.  The block calculation below
    maximizes the exact weighted count in every corresponding cyclic window.
    """
    cells = np.asarray(cells_tuple, dtype=np.int64) % denominator
    units = np.asarray(
        tuple(
            unit
            for unit in range(1, denominator)
            if gcd(unit, denominator) == 1
        ),
        dtype=np.int64,
    )
    require(len(units), denominator)
    width = (denominator - 1) // 7 + 1
    worst_danger = -1
    worst_unit = None
    for start in range(0, len(units), 32):
        block = units[start : start + 32]
        residues = (cells[:, None] * block[None, :]) % denominator
        offsets = denominator * np.arange(len(block), dtype=np.int64)[:, None]
        encoded = residues.T + offsets
        counts = np.bincount(
            encoded.ravel(), minlength=len(block) * denominator
        ).reshape(len(block), denominator)
        extended = np.concatenate((counts, counts[:, : width - 1]), axis=1)
        cumulative = np.concatenate(
            (
                np.zeros((len(block), 1), dtype=np.int64),
                np.cumsum(extended, axis=1),
            ),
            axis=1,
        )
        danger_counts = np.max(
            cumulative[:, width : width + denominator]
            - cumulative[:, :denominator],
            axis=1,
        )
        offset = int(np.argmax(danger_counts))
        candidate = int(danger_counts[offset])
        candidate_unit = int(block[offset])
        if candidate > worst_danger or (
            candidate == worst_danger
            and (worst_unit is None or candidate_unit < worst_unit)
        ):
            worst_danger = candidate
            worst_unit = candidate_unit
    require(worst_unit is not None and worst_danger >= 0, denominator)
    return len(cells_tuple) - worst_danger, worst_unit, len(units), worst_danger


def two_high_scalar_cases(engine, stream, residual, low, high):
    """Finite superset of every scalar-admissible packet with >=2 highs."""
    required = stream.lower - stream.first_delta
    cases = set()
    three_high_masks = 0
    for divisors in residual:
        slots = engine.suffix_slots(divisors, stream.first_d)
        for mask in range(8):
            if mask.bit_count() < 2:
                continue
            high_ds = tuple(
                denominator
                for position, denominator in enumerate(slots)
                if (mask >> position) & 1
            )
            low_ds = tuple(
                denominator
                for position, denominator in enumerate(slots)
                if not ((mask >> position) & 1)
            )
            high_upper = sum((high[denominator][0] for denominator in high_ds), Fraction())
            if not low_ds:
                if high_upper >= required:
                    cases.add((divisors, high_ds, None, high_upper - required))
                    three_high_masks += 1
                continue
            low_d = low_ds[0]
            threshold = required - high_upper
            for value, label in low[low_d]:
                if value < threshold:
                    break
                cases.add((divisors, high_ds, label, value - threshold))
    return tuple(sorted(cases)), three_high_masks


def terminal_worker(item):
    index, body, residual = item
    audit = load(f"prefix_audit_{index}", AUDIT_SOURCE)
    natural = audit.load(f"prefix_natural_{index}", audit.NATURAL_SOURCE)
    engine = audit.status_engine(natural, f"prefix_terminal_{index}")
    engine.FIRST = LEVEL
    engine.ray.FIRST = LEVEL
    stream = engine.ray.Stream(body)
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

    cell_cache = {}
    torsion_cache = {}
    qualifying = Counter()
    effective = Counter()
    one_failures = []
    for divisors, high_d, low_records, excess in one_high:
        low_labels = tuple(sorted(record[1] for record in low_records))
        cell_cache.setdefault(
            low_labels, engine.fixed_safe_cells(stream, low_labels)
        )
        key = (low_labels, high_d)
        torsion_cache.setdefault(
            key, engine.torsion_pigeonhole(cell_cache[low_labels], high_d)
        )
        witness = torsion_cache[key]
        if witness[0] is None:
            one_failures.append(
                (divisors, high_d, low_records, excess, witness)
            )
        else:
            qualifying[witness[0]] += 1
            effective[witness[1]] += 1
    require(not one_failures, (index, "one-high torsion failures", one_failures))

    two_cases = ()
    common_closed = ()
    point_closed = ()
    if two_gap <= 0:
        two_cases, three_high_masks = two_high_scalar_cases(
            engine, stream, residual, low, high
        )
        require(not three_high_masks, (index, "three-high scalar survivor"))
        common_records = []
        point_records = []
        point_cache = {}
        for divisors, high_ds, low_label, excess in two_cases:
            require(low_label is not None and len(high_ds) == 2, (index, high_ds))
            fixed = (low_label,)
            cell_cache.setdefault(fixed, engine.fixed_safe_cells(stream, fixed))
            cells = cell_cache[fixed]
            common = common_torsion_capacity(cells, high_ds)
            if common[-1] > 0:
                common_records.append(
                    (divisors, high_ds, low_label, excess, len(cells), common)
                )
                continue
            minima = []
            for denominator in high_ds:
                key = (fixed, denominator)
                point_cache.setdefault(
                    key, translated_safe_minimum(cells, denominator)
                )
                minima.append(point_cache[key])
            minima = tuple(minima)
            lower = minima[0][0] + minima[1][0] - len(cells)
            if lower > 0:
                certificate = ("translation-uniform-full-projection", lower)
            else:
                require(
                    high_ds[0] == high_ds[1]
                    and high_ds[0] == 2
                    and len(cells) > 0,
                    (index, divisors, high_ds, minima, lower),
                )
                # Write z_i=L/2+h_i L.  On any one fixed-safe cell, the
                # inverse image in the local coordinate y of a high danger
                # comb has measure at most
                #
                #   (h_i+1)/(7(h_i+1/2)) <= 2/7.
                #
                # Indeed h_i full periods contribute h_i/7 and the remaining
                # half-period contributes at most 1/7.  Two such combs leave
                # projected measure at least 1-4/7=3/7=39/91>36/91.
                certificate = (
                    "denominator-two-single-cell-measure",
                    len(cells),
                    Fraction(3, 7),
                    Fraction(3, 91),
                )
            point_records.append(
                (
                    divisors,
                    high_ds,
                    low_label,
                    excess,
                    len(cells),
                    common,
                    minima,
                    certificate,
                )
            )
        common_closed = tuple(common_records)
        point_closed = tuple(point_records)
        require(
            len(common_closed) + len(point_closed) == len(two_cases),
            (index, "two-high partition"),
        )

    torsion_packet = tuple(sorted(torsion_cache.items()))
    return (
        index,
        body,
        stream.L,
        stream.high_floor,
        stream.first_d,
        len(residual),
        tuple(sorted(needed)),
        two_gap,
        two_gap_witness,
        len(zero_high),
        len(one_high),
        len({tuple(sorted(record[1] for record in low_records))
             for _ds, _high_d, low_records, _excess in one_high}),
        sign_census,
        recurrence_checks,
        len(torsion_cache),
        tuple(sorted(qualifying.items())),
        tuple(sorted(effective.items())),
        digest(one_high),
        digest(torsion_packet),
        len(two_cases),
        len(common_closed),
        len(point_closed),
        digest(two_cases),
        digest(common_closed),
        digest(point_closed),
        min(
            (
                record[-1][1]
                for record in point_closed
                if record[-1][0] == "translation-uniform-full-projection"
            ),
            default=None,
        ),
        tuple(
            sorted(
                record[-1][1]
                for record in point_closed
                if record[-1][0] == "denominator-two-single-cell-measure"
            )
        ),
    )


def reconstruct_queue():
    fourth = load("cumulative_fourth", FOURTH_SOURCE)
    third = load("cumulative_third", fourth.THIRD_SOURCE)
    audit = load("cumulative_audit", third.AUDIT_SOURCE)
    natural = audit.load("cumulative_natural", audit.NATURAL_SOURCE)
    order = audit.load("cumulative_order", audit.ORDER_SOURCE)
    rows, components = natural.atlas_rows()
    order_rows, order_components = order.atlas_rows()
    require((rows, components) == (order_rows, order_components), "atlas parsers")
    census = (
        len(rows),
        sum(row[3] for row in rows),
        sum(not row[3] for row in rows),
    )
    require(census == EXPECTED_ATLAS, census)

    closed = {
        index
        for index, row in enumerate(rows)
        if gcd(LEVEL, row[1]) in (8, 18)
        or not row[3]
        or (gcd(LEVEL, row[1]), row[1]) in audit.NATURAL_FAMILY_KEYS
    }
    inherited = []
    for prefix_number in range(1, 5):
        live = tuple(
            index
            for index, row in enumerate(rows)
            if row[3] and index not in closed
        )
        ranked = audit.ranked_families(rows, components, live)
        selected = audit.through_next_nonsingleton(ranked)
        indices = tuple(index for family in selected for index in family[4])
        inherited.append((prefix_number, len(live), len(ranked), selected, indices))
        closed.update(indices)
    require(
        tuple(len(record[4]) for record in inherited)
        == EXPECTED_INHERITED_PREFIX_ROWS,
        inherited,
    )
    require(
        inherited[-1][3] == fourth.EXPECTED_SELECTED_FAMILIES,
        "THM3320 selected families",
    )
    require(
        inherited[-1][4] == fourth.EXPECTED_SELECTED_INDICES,
        "THM3320 selected rows",
    )

    live = tuple(
        index
        for index, row in enumerate(rows)
        if row[3] and index not in closed
    )
    ranked = audit.ranked_families(rows, components, live)
    require((len(live), len(ranked)) == EXPECTED_QUEUE_STATE, (len(live), len(ranked)))
    prefix = tuple(ranked[:PREFIX_FAMILIES])
    prefix_indices = tuple(index for family in prefix for index in family[4])
    require(
        (len(prefix_indices), len(prefix)) == EXPECTED_PREFIX_STATE,
        (len(prefix_indices), len(prefix)),
    )
    require(
        (prefix[-1][0], ranked[PREFIX_FAMILIES][0]) == EXPECTED_LAST_BOUNDARY,
        (prefix[-1][0], ranked[PREFIX_FAMILIES][0]),
    )
    prefix_hash = digest(prefix)
    if EXPECTED_PREFIX_SHA256 is not None:
        require(prefix_hash == EXPECTED_PREFIX_SHA256, prefix_hash)
    return rows, components, inherited, ranked, prefix, prefix_indices, prefix_hash


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=4)
    args = parser.parse_args()
    require(args.processes >= 1, args.processes)

    dependencies = (
        (FOURTH_SOURCE, FOURTH_SOURCE_SHA256, None),
        (FOURTH_OUTPUT, FOURTH_OUTPUT_SHA256,
         f"semantic_sha256={FOURTH_SEMANTIC_SHA256}"),
        (AUDIT_SOURCE, AUDIT_SOURCE_SHA256, None),
        (NATURAL_SOURCE, NATURAL_SOURCE_SHA256, None),
        (TORSION_SOURCE, TORSION_SOURCE_SHA256, None),
    )
    for path, expected, needle in dependencies:
        require(lf_sha(path) == expected, (path, "dependency changed"))
        if needle is not None:
            require(needle in path.read_text(encoding="utf-8"), (path, needle))

    syntax = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    assertion_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(syntax))
    float_nodes = sum(
        isinstance(node, ast.Constant) and isinstance(node.value, float)
        for node in ast.walk(syntax)
    )
    require((assertion_nodes, float_nodes) == (0, 0), "truth-gate syntax")

    (
        rows,
        components,
        inherited,
        ranked,
        prefix,
        prefix_indices,
        prefix_hash,
    ) = reconstruct_queue()
    tasks = tuple((index, rows[index]) for index in prefix_indices)
    if args.processes == 1:
        screens = tuple(screen_worker(task) for task in tasks)
    else:
        with mp.get_context("spawn").Pool(min(args.processes, len(tasks))) as pool:
            screens = tuple(pool.imap(screen_worker, tasks, chunksize=1))
    require(tuple(row[0] for row in screens) == prefix_indices, "screen order")

    canonical_screens = []
    screen_by_index = {}
    for index, screen, direct, legacy in screens:
        source = rows[index]
        divisor_gcd = gcd(LEVEL, source[1])
        require(
            screen[:6]
            == (
                LEVEL,
                source[0],
                source[1],
                source[2],
                source[1] // divisor_gcd,
                source[3],
            ),
            (index, screen[:6]),
        )
        require(screen[16] == screen[11], (index, "status verification"))
        require(direct + legacy == screen[11], (index, "Farkas partition"))
        require(screen[12] == len(screen[13]), (index, "residual count"))
        record = (index, tuple(screen[:19]), direct, legacy)
        canonical_screens.append(record)
        screen_by_index[index] = record
    canonical_screens = tuple(canonical_screens)
    totals = tuple(
        sum(record[1][position] for record in canonical_screens)
        for position in (9, 10, 11, 12)
    )
    farkas = (
        sum(record[2] for record in canonical_screens),
        sum(record[3] for record in canonical_screens),
    )
    require(totals == EXPECTED_SCREEN_TOTALS, totals)
    require(farkas == EXPECTED_FARKAS_TOTALS, farkas)
    screen_hash = digest(canonical_screens)
    if EXPECTED_SCREEN_SHA256 is not None:
        require(screen_hash == EXPECTED_SCREEN_SHA256, screen_hash)

    terminal_tasks = tuple(
        (index, rows[index][0], record[1][13])
        for index, record in sorted(screen_by_index.items())
        if record[1][12]
    )
    if args.processes == 1:
        terminals = tuple(terminal_worker(task) for task in terminal_tasks)
    else:
        with mp.get_context("spawn").Pool(
            min(args.processes, len(terminal_tasks))
        ) as pool:
            terminals = tuple(pool.imap(terminal_worker, terminal_tasks, chunksize=1))
    require(
        tuple(record[0] for record in terminals)
        == tuple(task[0] for task in terminal_tasks),
        "terminal order",
    )
    require(sum(record[5] for record in terminals) == totals[3], "terminal passports")
    terminal_hash = digest(terminals)
    if EXPECTED_TERMINAL_SHA256 is not None:
        require(terminal_hash == EXPECTED_TERMINAL_SHA256, terminal_hash)

    excluded = len(prefix_indices)
    require(EXPECTED_LEDGER[0] - excluded == EXPECTED_LEDGER[1], "ledger")
    require(EXPECTED_WALL[0] - excluded == EXPECTED_WALL[1], "wall")
    require(
        EXPECTED_FAMILY_COUNT[0] - len(prefix) == EXPECTED_FAMILY_COUNT[1],
        "family count",
    )

    family_records = []
    family_lines = []
    for rank, family in enumerate(prefix, 1):
        cost, count, divisor_gcd, ruler, indices, _packet = family
        family_screens = tuple(screen_by_index[index] for index in indices)
        family_totals = tuple(
            sum(record[1][position] for record in family_screens)
            for position in (9, 10, 11, 12)
        )
        family_farkas = (
            sum(record[2] for record in family_screens),
            sum(record[3] for record in family_screens),
        )
        record = (
            rank,
            cost,
            count,
            divisor_gcd,
            ruler,
            indices,
            family_totals,
            family_farkas,
            digest(family_screens),
        )
        family_records.append(record)
        family_lines.append(
            f"FAMILY;rank={rank};gcd={divisor_gcd};L={ruler};rows={count};"
            f"cost={cost};indices={','.join(map(str, indices))};"
            f"states={family_totals[0]};crude={family_totals[1]};"
            f"status={family_totals[2]};residual={family_totals[3]};"
            f"direct={family_farkas[0]};legacy={family_farkas[1]};"
            f"screen_sha256={record[-1]}"
        )
    family_records = tuple(family_records)

    positive_gaps = sum(record[7] > 0 for record in terminals)
    nonpositive_gaps = sum(record[7] <= 0 for record in terminals)
    terminal_summary = (
        len(terminals),
        sum(record[5] for record in terminals),
        positive_gaps,
        nonpositive_gaps,
        sum(record[9] for record in terminals),
        sum(record[10] for record in terminals),
        sum(record[14] for record in terminals),
        sum(record[19] for record in terminals),
        sum(record[20] for record in terminals),
        sum(record[21] for record in terminals),
    )
    require(terminal_summary[:4] == (35, 518, 30, 5), terminal_summary)
    require(terminal_summary[7:] == (21, 3, 18), terminal_summary)

    semantic_packet = (
        tuple((path.name, expected) for path, expected, _needle in dependencies),
        tuple((record[0], record[1], record[2], record[4]) for record in inherited),
        prefix,
        prefix_hash,
        canonical_screens,
        screen_hash,
        family_records,
        terminals,
        terminal_hash,
        totals,
        farkas,
        terminal_summary,
        EXPECTED_LEDGER,
        EXPECTED_WALL,
        EXPECTED_FAMILY_COUNT,
        ranked[PREFIX_FAMILIES],
    )
    semantic_hash = digest(semantic_packet)
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, semantic_hash)

    lines = [
        "LRC14 projected-k3 z216 cumulative intrinsic-ruler prefix closure",
        "status=FINITE-EXACT;scope=necessary_projected_k3_z1_216_atlas_only",
        "dependency="
        + ";".join(f"{path.name}:{expected}" for path, expected, _needle in dependencies),
        f"runtime=python:{'.'.join(map(str, tuple(__import__('sys').version_info[:3])))};"
        f"numpy:{np.__version__};processes:{args.processes}",
        f"universe=atlas:{EXPECTED_ATLAS[0]};wall:{EXPECTED_ATLAS[1]};"
        f"order:{EXPECTED_ATLAS[2]};after_THM3320_wall:{len(ranked) and EXPECTED_QUEUE_STATE[0]};"
        f"families:{len(ranked)};prefix_families:{len(prefix)};prefix_rows:{len(prefix_indices)};"
        f"prefix_sha256:{prefix_hash}",
        f"queue_boundary=last_cost:{prefix[-1][0]};last:gcd{prefix[-1][2]}/L{prefix[-1][3]};"
        f"next_cost:{ranked[PREFIX_FAMILIES][0]};next:gcd{ranked[PREFIX_FAMILIES][2]}/"
        f"L{ranked[PREFIX_FAMILIES][3]};next_indices:{','.join(map(str, ranked[PREFIX_FAMILIES][4]))}",
        *family_lines,
        f"SCREEN_TOTAL;states={totals[0]};crude={totals[1]};status={totals[2]};"
        f"residual={totals[3]};direct={farkas[0]};legacy={farkas[1]};"
        f"screen_sha256={screen_hash}",
    ]
    for record in terminals:
        (
            index,
            body,
            ruler,
            high_floor,
            first_d,
            residual_count,
            _needed,
            gap,
            _gap_witness,
            zero_count,
            one_count,
            pair_count,
            _signs,
            recurrence_checks,
            witness_count,
            qualifying_histogram,
            effective_histogram,
            one_hash,
            witness_hash,
            two_count,
            common_count,
            point_count,
            two_hash,
            common_hash,
            point_hash,
            point_lower,
            sole_unit_sizes,
        ) = record
        lines.append(
            f"TERMINAL;index={index};E={','.join(map(str, body))};L={ruler};"
            f"high={high_floor};first_d={first_d};residual={residual_count};"
            f"two_gap={gap};zero_high_hostiles={zero_count};one_high={one_count};"
            f"low_pairs={pair_count};ray_checks={recurrence_checks};"
            f"torsion_keys={witness_count};qualifying={qualifying_histogram};"
            f"effective={effective_histogram};one_sha256={one_hash};"
            f"torsion_sha256={witness_hash};two_high_cases={two_count};"
            f"common_closed={common_count};translated_closed={point_count};"
            f"translated_min_lower={point_lower};d2_fixed_cell_counts={sole_unit_sizes};"
            f"two_sha256={two_hash};common_sha256={common_hash};translated_sha256={point_hash}"
        )
    lines.extend(
        [
            f"TERMINAL_TOTAL;bodies={terminal_summary[0]};passports={terminal_summary[1]};"
            f"positive_two_gap={terminal_summary[2]};nonpositive_two_gap={terminal_summary[3]};"
            f"zero_high_hostiles={terminal_summary[4]};one_high_cases={terminal_summary[5]};"
            f"torsion_keys={terminal_summary[6]};two_high_cases={terminal_summary[7]};"
            f"common_closed={terminal_summary[8]};translated_closed={terminal_summary[9]};"
            f"terminal_sha256={terminal_hash}",
            "translated_orbit_lemma=for fixed-safe C and denominator d,"
            "B_C(d)=max over units u and translated open cyclic d/7 bands J "
            "of #{c in C:uc mod d in J}, and m_C^tr(d)=|C|-B_C(d);"
            "m_C^tr(d1)+m_C^tr(d2)>|C| forces every local coordinate to have "
            "a common cell safe from both high combs, hence P_(E,Z)=T",
            "d2_measure_lemma=for z=L/2+hL one fixed-safe cell loses at most "
            "2/7 of local-coordinate mass; two d2 highs leave at least "
            "3/7=39/91>36/91, contradicting the aligned three-comb cap",
            f"consequence=all_prefix_rows_excluded;ledger:{EXPECTED_LEDGER[0]}->{EXPECTED_LEDGER[1]};"
            f"wall:{EXPECTED_WALL[0]}->{EXPECTED_WALL[1]};"
            f"families:{EXPECTED_FAMILY_COUNT[0]}->{EXPECTED_FAMILY_COUNT[1]};"
            "projected_k3_cap=216_unchanged",
            "scope=no_physical_entry_or_arbitrary_k_or_rung_or_LRC14_conclusion",
            f"semantic_sha256={semantic_hash}",
            "all_exact_controls=PASS",
        ]
    )
    print("\n".join(lines))


if __name__ == "__main__":
    main()
