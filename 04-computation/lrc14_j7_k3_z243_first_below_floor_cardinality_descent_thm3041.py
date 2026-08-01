#!/usr/bin/env python3
"""Exact candidate closure of all THM-3033 first-below-floor z1=243 rows.

The promoted THM-3033 evaluator is imported unchanged for the exact
state/common-status screen.  This file changes only the atlas row selector.
Its terminal proof preserves the canonical two-high, one-high cardinality,
and one-high maximum-gap typing and reports any failure rather than assuming
the next frontier closes.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
import os
import pickle
import re
import sys
import tempfile
from collections import Counter
from fractions import Fraction as Q
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE_SOURCE = ROOT / "04-computation/lrc14_j7_k3_z246_to_z244_translated_band_descent_thm2941.py"
BASE_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_k3_z246_to_z244_translated_band_descent_thm2941.out"
OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_k3_z243_first_below_floor_cardinality_descent_thm3041.out"

BASE_SOURCE_SHA256 = "c063e1d94e1f30c4ccd63acfabb172988d23d96f680a961d9042b3f11862fa79"
BASE_OUTPUT_SHA256 = "d012f84d8e380019b4da8efbe9afd6695079c585f47e4aa2c474828ca62014e4"
BASE_SEMANTIC_SHA256 = "3f9945d333b13adef8aa8e7b162960bfde97f1e5a77945edcc735013d8ce4231"

FIRST = 243
EXPECTED_ALL_ROWS = 154
EXPECTED_HIGH_FLOOR_ROWS = 3
EXPECTED_ROWS = 151


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    payload = path.read_bytes()
    require(b"\r" not in payload.replace(b"\r\n", b""), f"bare CR in {path}")
    return hashlib.sha256(payload.replace(b"\r\n", b"\n")).hexdigest()


SCRIPT_SHA256 = sha(Path(__file__).resolve())


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


require(sha(BASE_SOURCE) == BASE_SOURCE_SHA256, "THM-3033 source changed")
require(sha(BASE_OUTPUT) == BASE_OUTPUT_SHA256, "THM-3033 output changed")
require(
    f"semantic_sha256={BASE_SEMANTIC_SHA256}" in BASE_OUTPUT.read_text(),
    "THM-3033 semantic digest changed",
)
thm = load("thm3033_z243_probe_base", BASE_SOURCE)
eng = thm.eng
ATLAS = thm.ATLAS


def atlas_tasks():
    pattern = re.compile(
        r"^row=E=([0-9,]+);.*;L=([0-9]+);high=([0-9]+);z1=([0-9]+);"
    )
    all_rows = []
    tasks = []
    for line in ATLAS.read_text().splitlines():
        match = pattern.match(line)
        if match is None:
            continue
        body = tuple(map(int, match.group(1).split(",")))
        L, high, first = map(int, match.group(2, 3, 4))
        if first != FIRST:
            continue
        require(high == thm.base.WALL.numerator * L // thm.base.WALL.denominator + 1,
                (body, L, high, "high floor"))
        task = (first, body, L, high, first < high)
        all_rows.append(task)
        if first < high:
            tasks.append(task)
    require(len(all_rows) == EXPECTED_ALL_ROWS, ("all z243 rows", len(all_rows)))
    require(len(all_rows) - len(tasks) == EXPECTED_HIGH_FLOOR_ROWS,
            ("high-floor rows", len(all_rows) - len(tasks)))
    require(len(tasks) == EXPECTED_ROWS, ("first-below-floor rows", len(tasks)))
    require(len({(task[0], task[1]) for task in tasks}) == len(tasks),
            "duplicate body key")
    return tuple(tasks), tuple(all_rows)


def worker_evaluate(task):
    return task, thm.evaluate(task)


def terminal_probe(task):
    first, body, residual = task
    eng.FIRST = first
    eng.ray.FIRST = first
    stream = eng.ray.Stream(body)
    require(first < stream.high_floor, (body, first, stream.high_floor))
    needed = {d for ds in residual for d in eng.suffix_slots(ds, stream.first_d)}
    low, high, low_signs, ray_checks = eng.build_literal_tables(stream, needed)
    gap, gap_witness = eng.duplicate_two_high_gap(stream, residual, low, high)
    zero = eng.zero_high_scalar_passes(stream, residual, low)
    cases = eng.one_high_cases(stream, residual, low, high)

    clean_cache = {}
    certificate_cache = {}
    cardinality_cases = 0
    maxgap_cases = 0
    failed_cases = []
    unit_checks = 0
    minimum_slack = None
    digest = hashlib.sha256()
    for ds, high_d, low_rows, excess in cases:
        labels = tuple(sorted(z for _d, z in low_rows))
        if labels not in clean_cache:
            clean_cache[labels] = eng.fixed_safe_cells(stream, labels)
        key = (labels, high_d)
        if key not in certificate_cache:
            certificate_cache[key] = thm.translated_certificate(
                clean_cache[labels], high_d
            )
        certificate = certificate_cache[key]
        minimum_slack = (
            certificate[3]
            if minimum_slack is None
            else min(minimum_slack, certificate[3])
        )
        unit_checks += certificate[6]
        if certificate[4] == "cardinality":
            cardinality_cases += 1
            require(certificate[5], (body, ds, high_d, "cardinality failed"))
        elif certificate[5]:
            maxgap_cases += 1
        else:
            failed_cases.append((ds, high_d, low_rows, excess, certificate))
        digest.update(
            repr((ds, high_d, low_rows, excess, certificate)).encode() + b"\n"
        )

    two_high_closed = gap > 0
    closed = two_high_closed and not failed_cases
    return (
        first,
        body,
        stream.L,
        stream.high_floor,
        stream.lower - stream.first_delta,
        len(residual),
        gap,
        gap_witness,
        len(zero),
        len(cases),
        len({tuple(sorted(z for _d, z in rows)) for _ds, _hd, rows, _e in cases}),
        low_signs,
        ray_checks,
        cardinality_cases,
        maxgap_cases,
        len(failed_cases),
        unit_checks,
        minimum_slack,
        len(certificate_cache),
        digest.hexdigest(),
        tuple(failed_cases),
        two_high_closed,
        closed,
    )


def worker_terminal(task):
    return task, terminal_probe(task)


def checkpoint_path(checkpoint_dir, label, task):
    digest = hashlib.sha256(repr((label, task)).encode()).hexdigest()
    return checkpoint_dir / f"{label}-{digest}.pickle"


def load_checkpoint(checkpoint_dir, label, task):
    path = checkpoint_path(checkpoint_dir, label, task)
    if not path.exists():
        return None
    payload = pickle.loads(path.read_bytes())
    expected = (
        label,
        task,
        SCRIPT_SHA256,
        BASE_SOURCE_SHA256,
        BASE_OUTPUT_SHA256,
        BASE_SEMANTIC_SHA256,
    )
    require(payload[:6] == expected, (path, "checkpoint provenance"))
    return payload[6]


def save_checkpoint(checkpoint_dir, label, task, row):
    path = checkpoint_path(checkpoint_dir, label, task)
    payload = (
        label,
        task,
        SCRIPT_SHA256,
        BASE_SOURCE_SHA256,
        BASE_OUTPUT_SHA256,
        BASE_SEMANTIC_SHA256,
        row,
    )
    temporary = path.with_name(f"{path.name}.tmp-{os.getpid()}")
    temporary.write_bytes(pickle.dumps(payload, protocol=5))
    os.replace(temporary, path)


def run_pool(function, tasks, processes, label, checkpoint_dir):
    if not tasks:
        return ()
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    pending = []
    for task in tasks:
        cached = load_checkpoint(checkpoint_dir, label, task)
        if cached is None:
            pending.append(task)
        else:
            rows.append(cached)
    print(
        f"{label}_checkpoint_recovered={len(rows)};pending={len(pending)}",
        file=sys.stderr,
        flush=True,
    )

    def retain(envelope):
        task, row = envelope
        save_checkpoint(checkpoint_dir, label, task, row)
        rows.append(row)
        completed = len(rows)
        if completed % 5 == 0 or completed == len(tasks):
            print(
                f"{label}_progress={completed}/{len(tasks)}",
                file=sys.stderr,
                flush=True,
            )

    if processes == 1:
        for task in pending:
            retain(function(task))
        return tuple(rows)
    with mp.get_context("spawn").Pool(min(processes, len(tasks))) as pool:
        for envelope in pool.imap_unordered(function, pending, chunksize=1):
            retain(envelope)
    return tuple(rows)


def ft(value):
    return "NONE" if value is None else f"{value.numerator}/{value.denominator}"


def render(tasks, all_rows, records, terminals):
    record_by_key = {(row[0], row[1]): row for row in records}
    require(len(record_by_key) == len(tasks), "lost or duplicate evaluated row")
    records = tuple(record_by_key[(task[0], task[1])] for task in tasks)
    terminal_by_key = {(row[0], row[1]): row for row in terminals}
    require(len(terminal_by_key) == len(terminals), "duplicate terminal row")
    terminal_tasks = tuple(row for row in records if row[12])
    terminals = tuple(terminal_by_key[(row[0], row[1])] for row in terminal_tasks)

    totals = tuple(sum(row[index] for row in records) for index in (9, 10, 11, 12))
    require(totals[0] == totals[1] + totals[2] + totals[3], totals)
    require(sum(row[16] for row in records) == totals[2], "Farkas verification count")
    screen_closed = sum(row[12] == 0 for row in records)
    terminal_closed = sum(row[22] for row in terminals)
    survivors = tuple(row for row in terminals if not row[22])
    closed_rows = screen_closed + terminal_closed
    require(closed_rows + len(survivors) == EXPECTED_ROWS, "row accounting")

    terminal_totals = (
        len(terminals),
        sum(row[21] for row in terminals),
        sum(not row[21] for row in terminals),
        sum(row[8] for row in terminals),
        sum(row[9] for row in terminals),
        sum(row[13] for row in terminals),
        sum(row[14] for row in terminals),
        sum(row[15] for row in terminals),
        sum(row[16] for row in terminals),
    )
    maxgap_control = thm.affine_small_group_exhaustive()
    boundary_controls = thm.affine_boundary_controls()
    row_order = tuple((task[1], task[2], task[3]) for task in tasks)
    semantic_payload = (
        FIRST,
        SCRIPT_SHA256,
        BASE_SOURCE_SHA256,
        BASE_OUTPUT_SHA256,
        BASE_SEMANTIC_SHA256,
        tasks,
        all_rows,
        records,
        terminals,
        totals,
        terminal_totals,
        survivors,
        maxgap_control,
        boundary_controls,
    )
    semantic = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()

    lines = [
        "LRC14 THM3041 first-below-floor z1=243 exact candidate",
        f"script_sha256={SCRIPT_SHA256}",
        f"base_source_sha256={sha(BASE_SOURCE)}",
        f"base_output_sha256={sha(BASE_OUTPUT)}",
        f"base_semantic_sha256={BASE_SEMANTIC_SHA256}",
        "dependency_hash_basis=SHA-256 after CRLF-to-LF normalization;bare CR rejected",
        f"universe=z1:{FIRST};all_rows:{len(all_rows)};high_floor_rows_excluded:{len(all_rows)-len(tasks)};first_below_floor_rows:{len(tasks)};unique_bodies:{len({task[1] for task in tasks})};L_min:{min(task[2] for task in tasks)};L_max:{max(task[2] for task in tasks)};row_order_sha256:{hashlib.sha256(repr(row_order).encode()).hexdigest()}",
        f"screen_totals=states:{totals[0]};crude:{totals[1]};status:{totals[2]};residual:{totals[3]};screen_closed_rows:{screen_closed};residual_bodies:{len(terminals)}",
        f"independent_exact_farkas_checks={totals[2]}/{totals[2]}:PASS;solver_basis_not_frozen",
        f"terminal_totals=residual_bodies:{terminal_totals[0]};two_high_closed:{terminal_totals[1]};two_high_open:{terminal_totals[2]};zero_high_hostile:{terminal_totals[3]};one_high_cases:{terminal_totals[4]};cardinality_cases:{terminal_totals[5]};maxgap_cases:{terminal_totals[6]};failed_one_high_cases:{terminal_totals[7]};maxgap_unit_checks:{terminal_totals[8]}",
        f"closure=closed_rows:{closed_rows};surviving_rows:{len(survivors)};all_151_close:{not survivors}",
        f"affine_maxgap_control=checks:{maxgap_control};boundary_controls:{boundary_controls}",
    ]
    for row in records:
        first, body, L, high, d1, gate, trials, checks, signs, states, crude, status, residual, residual_ds, mhist, instance_digest, verified, representative, stage_digest = row
        require(gate and d1 == L // gcd(L, first), (body, "wall task typing"))
        lines.append(
            f"BODY;E={body};L={L};high={high};d1={d1};trials={trials};checks={checks};signs={dict(signs)};states={states};crude={crude};status={status};residual={residual};M={dict(mhist)};instance_sha256={instance_digest};verified={verified};representative={representative};stage_sha256={stage_digest};residual_sha256={hashlib.sha256(repr(residual_ds).encode()).hexdigest()}"
        )
    for row in terminals:
        first, body, L, high, required, residual, gap, witness, zero, cases, pairs, low_signs, checks, card, shape, failures, unit_checks, slack, certificates, digest, failed_cases, two_high_closed, closed = row
        lines.append(
            f"TERMINAL;E={body};L={L};high={high};required={ft(required)};residual={residual};two_high_gap={ft(gap)};two_high_closed={two_high_closed};gap_witness={witness};zero_high={zero};one_high_cases={cases};low_pairs={pairs};low_signs={dict(low_signs)};ray_checks={checks};cardinality_cases={card};maxgap_cases={shape};failed_one_high_cases={failures};maxgap_unit_checks={unit_checks};minimum_cardinality_slack={slack};certificates={certificates};case_sha256={digest};closed={closed}"
        )
    for row in survivors:
        body = row[1]
        reason = []
        if not row[21]:
            reason.append(("two-high", row[6], row[7]))
        if row[20]:
            reason.append(("one-high", row[20]))
        lines.append(f"SURVIVOR;E={body};reasons={tuple(reason)}")
    lines += [
        "scope=projected necessary-sector candidate only;cap and ledger consequence inactive before independent audit;not LRC14",
        f"semantic_sha256={semantic}",
        "all_exact_controls=PASS",
    ]
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=8)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument(
        "--checkpoint-dir",
        type=Path,
        default=Path(tempfile.gettempdir())
        / "lrc14-thm3041-z243-cardinality-checkpoints",
    )
    args = parser.parse_args()
    require(args.processes >= 1, "positive process count required")
    tasks, all_rows = atlas_tasks()
    scheduled = tuple(sorted(tasks, key=lambda task: (-task[2], task[1])))
    records = run_pool(
        worker_evaluate,
        scheduled,
        args.processes,
        "screen",
        args.checkpoint_dir,
    )
    terminal_tasks = tuple(
        (row[0], row[1], row[13]) for row in records if row[12]
    )
    scheduled_terminals = tuple(
        sorted(terminal_tasks, key=lambda task: (-eng.ray.Stream(task[1]).L, task[1]))
    )
    terminals = run_pool(
        worker_terminal,
        scheduled_terminals,
        args.processes,
        "terminal",
        args.checkpoint_dir,
    )
    payload = render(tasks, all_rows, records, terminals)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
