#!/usr/bin/env python3
"""Exact projected-k=3 descent at z1=246,245,244.

This compositor derives its row universe from the pinned THM-2941 atlas.  It
inherits THM-2981's exact state/status engine but not its height list.  Rows
with first label below the integer high floor use the wall + two-high gap +
one-high translated-band gate.  Rows with first label already at the floor
use strict label order + the same two-high gap.  The affine max-gap test is
invoked only when the sharp cardinality test |S|>ceil(d/7) fails.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
import re
from collections import Counter, defaultdict
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE_SOURCE = ROOT / "04-computation/lrc14_j7_k3_z270_to_z247_cardinality_torsion_descent_thm2941.py"
BASE_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_k3_z270_to_z247_cardinality_torsion_descent_thm2941.out"
OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_k3_z246_to_z244_translated_band_descent_thm2941.out"

BASE_SOURCE_SHA256 = "f41fbdec062ad8c0d8b22ba1bab87d47902efd3591fd5fb6714aeb96f07de651"
BASE_OUTPUT_SHA256 = "45762f14e62ca5c9381f2e283f78408f80c8fe48e860648e53407b29322d9125"
BASE_SEMANTIC_SHA256 = "5b9789fbb40a557837ddcad3f3e2d6b5527640601d22d27636066920704367bd"
SEMANTIC_SHA256 = None

LEVELS = (246, 245, 244)
ROW_COUNTS = {246: 194, 245: 1, 244: 2}
# Per level: (first below high floor, first at or above high floor).
HIGH_FLOOR_ROW_COUNTS = {246: (156, 38), 245: (1, 0), 244: (0, 2)}
NEXT_HEIGHT = 243
NEXT_COUNT = 154
INHERITED_LEDGER = 375251
FINAL_LEDGER = 375054

# Filled only after the development scout derives the complete exact census.
LEVEL_TOTALS = None
TOTALS = None
TERMINAL_TOTALS = None
MIN_TWO_HIGH_GAP = None
ORDER_LEVEL_TOTALS = None


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    payload = path.read_bytes()
    require(b"\r" not in payload.replace(b"\r\n", b""), f"bare CR in {path}")
    return hashlib.sha256(payload.replace(b"\r\n", b"\n")).hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


require(sha(BASE_SOURCE) == BASE_SOURCE_SHA256, "THM-2981 source changed")
require(sha(BASE_OUTPUT) == BASE_OUTPUT_SHA256, "THM-2981 output changed")
require(
    f"semantic_sha256={BASE_SEMANTIC_SHA256}" in BASE_OUTPUT.read_text(),
    "THM-2981 semantic changed",
)
require(
    "projected k3 cap<=246;next exact frontier=246" in BASE_OUTPUT.read_text(),
    "THM-2981 conclusion changed",
)
base = load("lrc14_k3_z246_z244_base", BASE_SOURCE)
eng = base.eng
ATLAS = base.ATLAS


def atlas_tasks():
    pattern = re.compile(
        r"^row=E=([0-9,]+);.*;L=([0-9]+);high=([0-9]+);z1=([0-9]+);"
    )
    counts = Counter()
    tasks = []
    order = []
    for line in ATLAS.read_text().splitlines():
        match = pattern.match(line)
        if match is None:
            continue
        body = tuple(map(int, match.group(1).split(",")))
        L, high, first = map(int, match.group(2, 3, 4))
        counts[first] += 1
        if first in LEVELS:
            require(high == base.WALL.numerator * L // base.WALL.denominator + 1,
                    (first, body, "high floor"))
            tasks.append((first, body, L, high, first < high))
            order.append((first, body, L, high))
    derived = {z: counts[z] for z in LEVELS}
    require(derived == ROW_COUNTS, ("derived row universe", derived))
    split = {
        z: (sum(first == z and gate for first, _body, _L, _high, gate in tasks),
            sum(first == z and not gate for first, _body, _L, _high, gate in tasks))
        for z in LEVELS
    }
    require(split == HIGH_FLOOR_ROW_COUNTS, ("high-floor row split", split))
    for z in range(NEXT_HEIGHT + 1, max(LEVELS) + 1):
        require(counts[z] == ROW_COUNTS.get(z, 0), (z, counts[z]))
    require(counts[NEXT_HEIGHT] == NEXT_COUNT, (NEXT_HEIGHT, counts[NEXT_HEIGHT]))
    require(len(tasks) == sum(ROW_COUNTS.values()), len(tasks))
    return tuple(tasks), tuple(order), tuple((z, counts[z]) for z in range(NEXT_HEIGHT, max(LEVELS) + 1))


def evaluate(task):
    """Replay the exact screen while preserving both high-floor branches."""
    first, body, atlas_L, atlas_high, gate = task
    eng.FIRST = first
    eng.ray.FIRST = first
    stream = eng.ray.Stream(body)
    require((stream.L, stream.high_floor) == (atlas_L, atlas_high),
            (first, body, "atlas"))
    trials, states, checks, signs = eng.ray.ray_quotient_states(stream)
    crude, status, residual = eng.exact_common_status_screen(stream, states)
    require(set(states) == set(crude) | set(status) | set(residual),
            (first, body, "partition"))
    require(not (set(crude) & set(status)), (first, body, "overlap"))
    for ds, witness in crude.items():
        gap, q, M, target, capacity = witness
        D = lcm(*ds)
        require(M == D // q and target == dict(stream.target_data(D))[q],
                (first, body, ds, "crude data"))
        require(capacity == sum(eng.ray.local.fibre_cap(D, d, q) for d in ds),
                (first, body, ds, "capacity"))
        require(gap == target - capacity and gap > 0,
                (first, body, ds, "crude gap"))
    instance_digest = hashlib.sha256()
    representative = None
    arcs_cache = {}
    histogram_cache = {}
    verified = 0
    for ds, witness in sorted(status.items()):
        q, M, marginals, cap_set, histogram, certificate = witness
        D = lcm(*ds)
        require(M == D // q, (first, body, ds, "status cofactor"))
        rebuilt, capacities = eng.ray.local.hunter_status_data(D, ds, q)
        require(rebuilt == marginals and tuple(sorted(set(capacities))) == cap_set,
                (first, body, ds, "status data"))
        if D not in arcs_cache:
            arcs_cache[D] = eng.ray.fibre.projected_support_arcs(D, stream.ranges)
        key = (D, q)
        if key not in histogram_cache:
            histogram_cache[key] = eng.ray.fibre.residue_load_histogram(arcs_cache[D], q)
        require(histogram_cache[key] == histogram, (first, body, ds, "histogram"))
        eng.independent_farkas_check(q, marginals, capacities, histogram, certificate)
        instance = witness[:-1]
        if representative is None:
            representative = (first, body, ds, instance)
        instance_digest.update(f"{first}|{body}|{ds}|{instance}\n".encode())
        verified += 1
    stage_digest = hashlib.sha256()
    for ds in sorted(states):
        if ds in crude:
            stage, witness = "crude", crude[ds]
        elif ds in status:
            stage, witness = "status", status[ds][:-1]
        else:
            stage, witness = "residual", None
        stage_digest.update(repr((ds, base.frozen_state(states[ds]), stage, witness)).encode() + b"\n")
    for ds in residual:
        labels = states[ds]["labels"]
        # The quotient maximizer intentionally relaxes global label order.
        # Order belongs to the original packet, not to this stored upper-bound
        # representative.  In the wall branch the DP itself must nevertheless
        # have selected a high representative; in the order branch no such
        # assertion is sound or needed.
        require(labels[0] == first and len(labels) == len(set(labels)) == 4,
                (first, body, ds, "label identity"))
        require(tuple(sorted(stream.L // gcd(stream.L, z) for z in labels)) == ds,
                (first, body, ds, "orders"))
        high_count = sum(z >= stream.high_floor for z in labels[1:])
        if gate:
            require(first < stream.high_floor and high_count >= 1,
                    (first, body, ds, "wall quotient grammar", labels))
        else:
            require(first >= stream.high_floor,
                    (first, body, ds, "order branch", labels))
    sign_totals = {
        sign: sum(n for (_d, tag), n in signs.items() if tag == sign)
        for sign in (-1, 0, 1)
    }
    require(sign_totals[-1] == sign_totals[1], (first, body, "sign symmetry"))
    return (
        first, body, stream.L, stream.high_floor, stream.first_d, gate,
        trials, checks, tuple(sorted(sign_totals.items())), len(states),
        len(crude), len(status), len(residual), tuple(residual),
        tuple(sorted(Counter(w[1] for w in status.values()).items())),
        instance_digest.hexdigest(), verified, representative,
        stage_digest.hexdigest(),
    )


def cyclic_gaps(residues, d):
    require(residues, (d, "empty residue set"))
    rows = tuple(sorted(residues))
    if len(rows) == 1:
        return (d,)
    return tuple(rows[i + 1] - rows[i] for i in range(len(rows) - 1)) + (d + rows[0] - rows[-1],)


def translated_certificate(cells, d):
    residues = tuple(sorted({cell % d for cell in cells}))
    kappa = (d + 6) // 7
    slack = len(residues) - kappa
    residue_sha = hashlib.sha256(repr(residues).encode()).hexdigest()
    if slack > 0:
        return (d, len(residues), kappa, slack, "cardinality", True, 0, 0,
                None, len(cells), residue_sha, None)
    if not residues:
        return (d, 0, kappa, slack, "max-gap", False, 0, 1,
                ("empty-fixed-safe-residue-set",), len(cells), residue_sha, None)
    threshold = d - kappa + 1
    blockers = []
    shape_digest = hashlib.sha256()
    unit_checks = 0
    units = (0,) if d == 1 else tuple(unit for unit in range(1, d) if gcd(unit, d) == 1)
    for unit in units:
        image = tuple(sorted({unit * residue % d for residue in residues}))
        gaps = cyclic_gaps(image, d)
        maximum = max(gaps)
        shape_digest.update(repr((unit, image, gaps, maximum)).encode() + b"\n")
        unit_checks += 1
        if maximum >= threshold:
            blockers.append((unit, maximum, gaps, image))
    return (d, len(residues), kappa, slack, "max-gap", not blockers,
            unit_checks, len(blockers), blockers[0] if blockers else None,
            len(cells), residue_sha, shape_digest.hexdigest())


def affine_small_group_exhaustive():
    """Independently compare max-gap with literal affine-block containment."""
    checks = 0
    for d in range(1, 13):
        kappa = (d + 6) // 7
        units = (0,) if d == 1 else tuple(u for u in range(1, d) if gcd(u, d) == 1)
        blocks = tuple({(a + j) % d for j in range(kappa)} for a in range(d))
        for mask in range(1 << d):
            residues = tuple(r for r in range(d) if (mask >> r) & 1)
            obstructed = any(
                {u * r % d for r in residues} <= block
                for u in units
                for block in blocks
            )
            certificate = translated_certificate(residues, d)
            require(certificate[5] == (not obstructed),
                    (d, residues, certificate, obstructed, "affine exhaustive"))
            checks += 1
    require(checks == 8190, checks)
    return checks


def terminal(task):
    first, body, residual, gate = task
    eng.FIRST = first
    eng.ray.FIRST = first
    stream = eng.ray.Stream(body)
    needed = {d for ds in residual for d in eng.suffix_slots(ds, stream.first_d)}
    low, high, low_signs, ray_checks = eng.build_literal_tables(stream, needed)
    gap, gap_witness = eng.duplicate_two_high_gap(stream, residual, low, high)
    require(gap > 0, (first, body, "two-high gap", gap))
    if not gate:
        require(first >= stream.high_floor, (first, body, "order branch"))
        require(not any(low.values()), (first, body, "unexpected strict-order lows"))
        return (first, body, stream.L, stream.high_floor, gate,
                stream.lower - stream.first_delta, len(residual), gap, gap_witness,
                "strict-order-two-high", 0, 0, 0, low_signs,
                ray_checks, 0, 0, 0, 0, None, None, 0, None, None)

    require(first < stream.high_floor, (first, body, "wall branch"))
    require(dict(low_signs).get(-1, 0) > 0, (first, body, "negative lows lost"))
    zero = eng.zero_high_scalar_passes(stream, residual, low)
    cases = eng.one_high_cases(stream, residual, low, high)
    low_pairs = {tuple(sorted(z for _d, z in rows)) for _ds, _hd, rows, _e in cases}
    clean_cache = {}
    certificate_cache = {}
    cardinality_cases = shape_cases = failed_cases = unit_checks = 0
    minimum_slack = None
    first_obstruction = None
    digest = hashlib.sha256()
    for ds, high_d, low_rows, excess in cases:
        labels = tuple(sorted(z for _d, z in low_rows))
        if labels not in clean_cache:
            clean_cache[labels] = eng.fixed_safe_cells(stream, labels)
        key = (labels, high_d)
        if key not in certificate_cache:
            certificate_cache[key] = translated_certificate(clean_cache[labels], high_d)
        certificate = certificate_cache[key]
        minimum_slack = certificate[3] if minimum_slack is None else min(minimum_slack, certificate[3])
        unit_checks += certificate[6]
        if certificate[4] == "cardinality":
            cardinality_cases += 1
        elif certificate[5]:
            shape_cases += 1
        else:
            failed_cases += 1
            if first_obstruction is None:
                first_obstruction = (ds, high_d, low_rows, excess, certificate)
        digest.update(repr((ds, high_d, low_rows, excess, certificate)).encode() + b"\n")
    return (first, body, stream.L, stream.high_floor, gate,
            stream.lower - stream.first_delta, len(residual), gap, gap_witness,
            "wall-one-high", len(zero), len(cases), len(low_pairs), low_signs,
            ray_checks, cardinality_cases, shape_cases, failed_cases, unit_checks,
            minimum_slack, first_obstruction, len(certificate_cache),
            digest.hexdigest(), None)


def ft(value):
    return "NONE" if value is None else f"{value.numerator}/{value.denominator}"


def render(records, terminals, row_order, atlas_counts, development):
    totals = tuple(sum(row[i] for row in records) for i in (9, 10, 11, 12))
    level_totals = {
        z: tuple(sum(row[i] for row in records if row[0] == z) for i in (9, 10, 11, 12))
        for z in LEVELS
    }
    order_level_totals = {
        z: tuple(sum(row[i] for row in records if row[0] == z and not row[5])
                 for i in (9, 10, 11, 12))
        for z in LEVELS
    }
    terminal_totals = (
        len(terminals),
        sum(row[6] for row in terminals if not row[4]),
        sum(row[10] for row in terminals if row[4]),
        sum(row[11] for row in terminals if row[4]),
        sum(row[15] for row in terminals),
        sum(row[16] for row in terminals),
        sum(row[17] for row in terminals),
        sum(row[18] for row in terminals),
    )
    minimum_gap = min((row[7] for row in terminals), default=None)
    affine_exhaustive_checks = affine_small_group_exhaustive()
    if not development:
        require(LEVEL_TOTALS is not None and level_totals == LEVEL_TOTALS, level_totals)
        require(TOTALS is not None and totals == TOTALS, totals)
        require(TERMINAL_TOTALS is not None and terminal_totals == TERMINAL_TOTALS, terminal_totals)
        require(MIN_TWO_HIGH_GAP is not None and minimum_gap == MIN_TWO_HIGH_GAP, minimum_gap)
        require(ORDER_LEVEL_TOTALS is not None and order_level_totals == ORDER_LEVEL_TOTALS,
                order_level_totals)
    require(all(row[5] or row[12] == 0 for row in records),
            "a first-at-or-above-floor row survived the exact screen")
    require(all(row[4] for row in terminals),
            "terminal translated-band/max-gap proof escaped the first-below-floor branch")
    require(INHERITED_LEDGER - len(records) == FINAL_LEDGER, "ledger")
    failed = sum(row[17] for row in terminals)
    semantic_payload = (
        LEVELS, ROW_COUNTS, records, terminals, row_order, atlas_counts,
        level_totals, totals, terminal_totals, minimum_gap,
        HIGH_FLOOR_ROW_COUNTS, order_level_totals, affine_exhaustive_checks,
        BASE_SOURCE_SHA256, BASE_OUTPUT_SHA256, BASE_SEMANTIC_SHA256,
        INHERITED_LEDGER, FINAL_LEDGER, NEXT_HEIGHT, NEXT_COUNT,
    )
    semantic = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if SEMANTIC_SHA256 is not None:
        require(semantic == SEMANTIC_SHA256, "semantic digest changed")
    lines = [
        "LRC14 projected k=3 z1=246 through z1=244 translated-band descent",
        f"base_source_sha256={sha(BASE_SOURCE)}",
        f"base_output_sha256={sha(BASE_OUTPUT)}",
        f"base_semantic_sha256={BASE_SEMANTIC_SHA256}",
        "dependency_hash_basis=SHA-256 after CRLF-to-LF normalization;bare CR rejected",
        f"derived_universe=rows:{len(records)};counts:{ROW_COUNTS};row_order_sha256:{hashlib.sha256(repr(row_order).encode()).hexdigest()};next_height:{NEXT_HEIGHT};next_rows:{NEXT_COUNT}",
        f"frontier_totals=states:{totals[0]};crude:{totals[1]};status:{totals[2]};residual:{totals[3]}",
        f"independent_exact_farkas_checks={totals[2]}/{totals[2]}:PASS;solver_basis_not_frozen",
        f"high_floor_dichotomy=first_below:{sum(row[5] for row in records)};first_at_or_above:{sum(not row[5] for row in records)};per_level:{HIGH_FLOOR_ROW_COUNTS};wall branch forces at least one later high;order branch forces every later label high",
        f"strict_order_screen=rows:{sum(not row[5] for row in records)};per_level_state_totals:{order_level_totals};residual:0;all forty first-at-or-above-floor rows close before the terminal translated-band/max-gap proof",
        f"terminal_totals=residual_bodies:{terminal_totals[0]};order_closed_states:{terminal_totals[1]};zero_high_hostile:{terminal_totals[2]};one_high_cases:{terminal_totals[3]};cardinality_cases:{terminal_totals[4]};maxgap_cases:{terminal_totals[5]};failed_cases:{terminal_totals[6]};maxgap_unit_checks:{terminal_totals[7]};minimum_two_high_gap:{ft(minimum_gap)}",
        "translated_gate=complete-cell residue count |S|>ceil(d/7) is tried first;only nonpositive slack invokes the exact affine maximum-gap criterion;centered beta and pair selection are unused",
        f"affine_maxgap_independent_control=all subsets and all units modulo d=1..12;literal affine consecutive-block containment versus cyclic maximum-gap;checks:{affine_exhaustive_checks}:PASS;includes empty,singleton,and wraparound boundaries",
    ]
    for z in LEVELS:
        lines.append(f"LEVEL;z1={z};rows={ROW_COUNTS[z]};states={level_totals[z][0]};crude={level_totals[z][1]};status={level_totals[z][2]};residual={level_totals[z][3]}")
    for row in records:
        first, body, L, high, d1, gate, trials, checks, signs, states, crude, status, residual, residual_ds, mhist, instance_digest, verified, representative, stage_digest = row
        require(verified == status, (first, body, "verified status"))
        lines.append(f"BODY;z1={first};E={body};L={L};high={high};d1={d1};wall_gate={gate};trials={trials};checks={checks};signs={dict(signs)};states={states};crude={crude};status={status};residual={residual};M={dict(mhist)};instance_sha256={instance_digest};verified={verified};representative={representative};stage_sha256={stage_digest};residual_sha256={hashlib.sha256(repr(residual_ds).encode()).hexdigest()}")
    for row in terminals:
        first, body, L, high, gate, required, residual, gap, witness, mode, zero, cases, pairs, low_signs, checks, card, shape, failures, unit_checks, slack, obstruction, certificates, digest, _ = row
        lines.append(f"TERMINAL;z1={first};E={body};L={L};high={high};wall_gate={gate};mode={mode};required={ft(required)};residual={residual};two_high_gap={ft(gap)};gap_witness={witness};zero_high={zero};cases={cases};low_pairs={pairs};low_signs={dict(low_signs)};ray_checks={checks};cardinality_cases={card};maxgap_cases={shape};failed_cases={failures};maxgap_unit_checks={unit_checks};minimum_cardinality_slack={slack};certificates={certificates};case_sha256={digest};first_obstruction={obstruction}")
    conclusion = (
        f"all projected k3 rows at heights246..244 are empty;projected k3 cap<={NEXT_HEIGHT};next exact frontier={NEXT_HEIGHT}"
        if failed == 0
        else f"OPEN;translated sidecars leave {failed} terminal cases;no cap improvement"
    )
    lines += [
        f"ledger_decrement={INHERITED_LEDGER}-{len(records)}={FINAL_LEDGER};decrement counts body rows,not quotient states",
        "nonconsequence=projected scalar-atlas descent only;does not close z1=243,k<=1,the rung,or LRC(14)",
        f"conclusion={conclusion}",
        f"semantic_sha256={semantic}",
        "all_exact_controls=PASS",
    ]
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=8)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--development", action="store_true")
    args = parser.parse_args()
    require(args.processes >= 1, "positive process count required")
    require(args.development or all(x is not None for x in (LEVEL_TOTALS, TOTALS, TERMINAL_TOTALS, MIN_TWO_HIGH_GAP, ORDER_LEVEL_TOTALS, SEMANTIC_SHA256)),
            "unfrozen scout; pass --development")
    tasks, row_order, atlas_counts = atlas_tasks()
    if args.processes == 1:
        records = tuple(evaluate(task) for task in tasks)
    else:
        # Long-period bodies dominate the exact status replay.  Dispatch them
        # first so the final workers do not sit idle behind one late large-L
        # row, then restore canonical atlas order before hashing/rendering.
        scheduled_tasks = tuple(sorted(tasks, key=lambda task: (-task[2], task[0], task[1])))
        with mp.get_context("spawn").Pool(min(args.processes, len(tasks))) as pool:
            scheduled_records = tuple(pool.map(evaluate, scheduled_tasks, chunksize=1))
        record_by_key = {(row[0], row[1]): row for row in scheduled_records}
        require(len(record_by_key) == len(tasks), "duplicate atlas task key")
        records = tuple(record_by_key[(task[0], task[1])] for task in tasks)
    terminal_tasks = tuple((row[0], row[1], row[13], row[5]) for row in records if row[12])
    if args.processes == 1:
        terminals = tuple(terminal(task) for task in terminal_tasks)
    else:
        scheduled_terminals = tuple(
            sorted(terminal_tasks, key=lambda task: (-eng.ray.Stream(task[1]).L, task[0], task[1]))
        )
        with mp.get_context("spawn").Pool(min(args.processes, len(terminal_tasks))) as pool:
            scheduled_terminal_rows = tuple(pool.map(terminal, scheduled_terminals, chunksize=1))
        terminal_by_key = {(row[0], row[1]): row for row in scheduled_terminal_rows}
        require(len(terminal_by_key) == len(terminal_tasks), "duplicate terminal task key")
        terminals = tuple(terminal_by_key[(task[0], task[1])] for task in terminal_tasks)
    payload = render(records, terminals, row_order, atlas_counts, args.development)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(payload)
    print(payload, end="")


if __name__ == "__main__":
    main()
