#!/usr/bin/env python3
"""Exact combined scratch compositor for projected k=3 levels 242, 241, 240.

It evaluates every occupied z1=242 and z1=240 row with the promoted THM-3033
screen/terminal, independently certifies that z1=241 has no atlas row, and
binds checkpoints to an explicit algorithm/dependency fingerprint.
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
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE_SOURCE = ROOT / "04-computation/lrc14_j7_k3_z246_to_z244_translated_band_descent_thm2941.py"
BASE_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_k3_z246_to_z244_translated_band_descent_thm2941.out"
OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_k3_z242_gap241_z240_compositional_descent_thm3052.out"
BASE_SOURCE_SHA256 = "c063e1d94e1f30c4ccd63acfabb172988d23d96f680a961d9042b3f11862fa79"
BASE_OUTPUT_SHA256 = "d012f84d8e380019b4da8efbe9afd6695079c585f47e4aa2c474828ca62014e4"
BASE_SEMANTIC_SHA256 = "3f9945d333b13adef8aa8e7b162960bfde97f1e5a77945edcc735013d8ce4231"
ATLAS_SHA256 = "cee82237ce1f51729813b9c916edd3353204c18172abe1d71278dee2c5562eda"
LEVELS = (242, 240)
GAP_LEVEL = 241
NEIGHBOR_COUNTS = {
    243: (154, 151, 3),
    242: (20, 16, 4),
    241: (0, 0, 0),
    240: (52, 49, 3),
    239: (4, 2, 2),
}
EXPECTED_SCREEN = {242: (1430, 279, 1125, 26), 240: (8019, 4504, 3504, 11)}
EXPECTED_ORDER = {242: (21, 5, 16, 0), 240: (3, 2, 1, 0)}
EXPECTED_TERMINAL = {
    242: (4, 4, 0, 26, 26, 26, 0, 0, 0),
    240: (2, 2, 0, 10, 11, 11, 0, 0, 0),
}
EXPECTED_ROW_DIGEST = {
    242: "e68ff08a2b5c7777c4dab61dbab063c3869a50d63427f90e955a8d0e29f6c86f",
    240: "2ff7c5c9818d7d33519f5950379f473330858e134e908c5112dfde19785cd588",
    "combined": "c7f760f8571b2510ecd218a0172e2a5ec13b0b58d65c015446b61c19257c99e5",
}
EXPECTED_SEMANTIC_SHA256 = "8e2e6a1ed90033f3d85ec9cd03acf326ffb0e67453bbc206e8f2ece843c7cd68"
CACHE_SCHEMA = "lrc14-k3-z242-gap241-z240-combined-v1"


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


require(sha(BASE_SOURCE) == BASE_SOURCE_SHA256, "THM-3033 source changed")
require(sha(BASE_OUTPUT) == BASE_OUTPUT_SHA256, "THM-3033 output changed")
require(f"semantic_sha256={BASE_SEMANTIC_SHA256}" in BASE_OUTPUT.read_text(),
        "THM-3033 semantic changed")
thm = load("thm3033_z242_gap241_z240_combined_base", BASE_SOURCE)
eng = thm.eng
ATLAS = thm.ATLAS
require(sha(ATLAS) == ATLAS_SHA256, "THM-2941 projected body atlas changed")


def atlas_tasks():
    pattern = re.compile(
        r"^row=E=([0-9,]+);.*;L=([0-9]+);high=([0-9]+);z1=([0-9]+);"
    )
    atlas_lines = ATLAS.read_text().splitlines()
    row_lines = tuple(line for line in atlas_lines if line.startswith("row="))
    require(len(row_lines) == 6060, ("atlas row-line count", len(row_lines)))
    rows = {}
    tasks = []
    for line in row_lines:
        match = pattern.match(line)
        require(match is not None, ("unparsed atlas row", line))
        body = tuple(map(int, match.group(1).split(",")))
        L, high, first = map(int, match.group(2, 3, 4))
        require(
            high == thm.base.WALL.numerator * L // thm.base.WALL.denominator + 1,
            (body, L, high, "high floor"),
        )
        rows.setdefault(first, []).append((body, L, high, first < high))
        if first in LEVELS:
            tasks.append((first, body, L, high, first < high))
    derived = {
        z: (
            len(rows.get(z, ())),
            sum(row[3] for row in rows.get(z, ())),
            sum(not row[3] for row in rows.get(z, ())),
        )
        for z in NEIGHBOR_COUNTS
    }
    require(derived == NEIGHBOR_COUNTS, ("neighbor census", derived))
    literal = {
        z: sum(f";z1={z};" in line for line in row_lines)
        for z in NEIGHBOR_COUNTS
    }
    require(
        literal == {z: values[0] for z, values in NEIGHBOR_COUNTS.items()},
        ("literal neighbor census", literal),
    )
    require(GAP_LEVEL not in rows, "z1=241 unexpectedly occupied")
    require(len(tasks) == 72, ("combined occupied rows", len(tasks)))
    for level in LEVELS:
        layer = tuple(task for task in tasks if task[0] == level)
        expected = NEIGHBOR_COUNTS[level]
        require(len(layer) == expected[0], (level, "row count"))
        require(sum(task[4] for task in layer) == expected[1], (level, "wall split"))
        require(sum(not task[4] for task in layer) == expected[2], (level, "order split"))
        row_order = tuple((task[1], task[2], task[3], task[4]) for task in layer)
        require(
            hashlib.sha256(repr(row_order).encode()).hexdigest()
            == EXPECTED_ROW_DIGEST[level],
            (level, "row-order digest"),
        )
    combined_order = tuple((task[1], task[2], task[3], task[4]) for task in tasks)
    require(
        hashlib.sha256(repr(combined_order).encode()).hexdigest()
        == EXPECTED_ROW_DIGEST["combined"],
        "combined row-order digest",
    )
    return tuple(tasks), derived


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
        minimum_slack = certificate[3] if minimum_slack is None else min(
            minimum_slack, certificate[3]
        )
        unit_checks += certificate[6]
        if certificate[4] == "cardinality":
            cardinality_cases += 1
            require(certificate[5], (body, ds, high_d, "cardinality failed"))
        elif certificate[5]:
            maxgap_cases += 1
        else:
            failed_cases.append((ds, high_d, low_rows, excess, certificate))
        digest.update(repr((ds, high_d, low_rows, excess, certificate)).encode() + b"\n")
    return (
        first, body, stream.L, stream.high_floor,
        stream.lower - stream.first_delta, len(residual), gap, gap_witness,
        len(zero), len(cases),
        len({tuple(sorted(z for _d, z in rows)) for _ds, _hd, rows, _e in cases}),
        low_signs, ray_checks, cardinality_cases, maxgap_cases,
        len(failed_cases), unit_checks, minimum_slack, len(certificate_cache),
        digest.hexdigest(), tuple(failed_cases), gap > 0,
        gap > 0 and not failed_cases,
    )


def worker_terminal(task):
    return task, terminal_probe(task)


def checkpoint_fingerprint():
    return hashlib.sha256(
        repr(
            (
                CACHE_SCHEMA,
                BASE_SOURCE_SHA256,
                BASE_OUTPUT_SHA256,
                BASE_SEMANTIC_SHA256,
                ATLAS_SHA256,
                EXPECTED_SCREEN,
                EXPECTED_ORDER,
                EXPECTED_TERMINAL,
                EXPECTED_ROW_DIGEST,
            )
        ).encode()
    ).hexdigest()


CACHE_FINGERPRINT = checkpoint_fingerprint()


def checkpoint_path(directory, label, task):
    task_digest = hashlib.sha256(repr(task).encode()).hexdigest()
    return directory / f"{label}-{CACHE_FINGERPRINT[:16]}-{task_digest}.pickle"


def run_pool(function, tasks, processes, label, directory):
    directory.mkdir(parents=True, exist_ok=True)
    rows = []
    pending = []
    for task in tasks:
        path = checkpoint_path(directory, label, task)
        if path.exists():
            saved_fingerprint, saved_task, row = pickle.loads(path.read_bytes())
            require(saved_fingerprint == CACHE_FINGERPRINT, (path, "fingerprint mismatch"))
            require(saved_task == task, (path, "task mismatch"))
            rows.append(row)
        else:
            pending.append(task)
    print(f"{label}_recovered={len(rows)};pending={len(pending)}", file=sys.stderr, flush=True)

    def retain(envelope):
        task, row = envelope
        path = checkpoint_path(directory, label, task)
        temporary = path.with_name(f"{path.name}.tmp-{os.getpid()}")
        temporary.write_bytes(
            pickle.dumps((CACHE_FINGERPRINT, task, row), protocol=5)
        )
        os.replace(temporary, path)
        rows.append(row)
        print(
            f"{label}_progress={len(rows)}/{len(tasks)};z1={task[0]};body={task[1]}",
            file=sys.stderr,
            flush=True,
        )

    if processes == 1:
        for task in pending:
            retain(function(task))
    else:
        with mp.get_context("spawn").Pool(min(processes, len(pending))) as pool:
            for envelope in pool.imap_unordered(function, pending, chunksize=1):
                retain(envelope)
    return tuple(rows)


def four_totals(rows):
    return tuple(sum(row[index] for row in rows) for index in (9, 10, 11, 12))


def terminal_totals(rows):
    return (
        len(rows),
        sum(row[21] for row in rows),
        sum(not row[21] for row in rows),
        sum(row[8] for row in rows),
        sum(row[9] for row in rows),
        sum(row[13] for row in rows),
        sum(row[14] for row in rows),
        sum(row[15] for row in rows),
        sum(row[16] for row in rows),
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=4)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument(
        "--checkpoint-dir",
        type=Path,
        default=Path(tempfile.gettempdir()) / "lrc14-z242-gap241-z240-combined",
    )
    args = parser.parse_args()
    tasks, neighbor_census = atlas_tasks()
    scheduled = tuple(sorted(tasks, key=lambda task: (-task[2], -task[0], task[1])))
    records = run_pool(
        worker_evaluate, scheduled, args.processes, "screen", args.checkpoint_dir
    )
    by_key = {(row[0], row[1]): row for row in records}
    require(len(by_key) == len(tasks), "lost screen row")
    records = tuple(by_key[(task[0], task[1])] for task in tasks)
    terminal_tasks = tuple(
        (row[0], row[1], row[13]) for row in records if row[5] and row[12]
    )
    terminals = run_pool(
        worker_terminal, terminal_tasks, args.processes, "terminal", args.checkpoint_dir
    )
    t_by_key = {(row[0], row[1]): row for row in terminals}
    require(len(t_by_key) == len(terminal_tasks), "lost terminal row")
    terminals = tuple(t_by_key[(task[0], task[1])] for task in terminal_tasks)

    level_screen = {}
    level_order = {}
    level_terminal = {}
    for level in LEVELS:
        layer_records = tuple(row for row in records if row[0] == level)
        layer_terminals = tuple(row for row in terminals if row[0] == level)
        level_screen[level] = four_totals(layer_records)
        level_order[level] = four_totals(tuple(row for row in layer_records if not row[5]))
        level_terminal[level] = terminal_totals(layer_terminals)
        require(level_screen[level] == EXPECTED_SCREEN[level], (level, "screen", level_screen[level]))
        require(level_order[level] == EXPECTED_ORDER[level], (level, "order", level_order[level]))
        require(level_terminal[level] == EXPECTED_TERMINAL[level], (level, "terminal", level_terminal[level]))

    totals = four_totals(records)
    order = four_totals(tuple(row for row in records if not row[5]))
    all_terminal = terminal_totals(terminals)
    require(totals == (9449, 4783, 4629, 37), ("combined screen", totals))
    require(order == (24, 7, 17, 0), ("combined order", order))
    require(all_terminal == (6, 6, 0, 36, 37, 37, 0, 0, 0), ("combined terminal", all_terminal))
    survivors = tuple(row for row in terminals if not row[22])
    require(not survivors, ("terminal survivors", survivors))

    semantic_packet = (
        tasks,
        neighbor_census,
        records,
        terminals,
        level_screen,
        level_order,
        level_terminal,
        totals,
        order,
        all_terminal,
        survivors,
        CACHE_FINGERPRINT,
    )
    semantic = hashlib.sha256(repr(semantic_packet).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, ("semantic digest", semantic))

    lines = [
        "LRC14 projected k3 z1=242, gap241, z1=240 combined exact scout",
        f"dependency=THM3033_source:{BASE_SOURCE_SHA256};output:{BASE_OUTPUT_SHA256};semantic:{BASE_SEMANTIC_SHA256};atlas:{ATLAS_SHA256}",
        f"checkpoint_fingerprint={CACHE_FINGERPRINT}",
        "neighbor_census=" + ";".join(
            f"z{z}:rows{values[0]}/wall{values[1]}/order{values[2]}"
            for z, values in sorted(neighbor_census.items(), reverse=True)
        ),
        f"universe=occupied_rows:{len(tasks)};wall:{sum(t[4] for t in tasks)};order:{sum(not t[4] for t in tasks)};gap_z241_rows:0;row_order_sha256:{EXPECTED_ROW_DIGEST['combined']}",
        f"screen=states:{totals[0]};crude:{totals[1]};status:{totals[2]};residual:{totals[3]};residual_bodies:{sum(row[12]>0 for row in records)}",
        f"order_screen=states:{order[0]};crude:{order[1]};status:{order[2]};residual:{order[3]}",
        f"terminal=residual_bodies:{all_terminal[0]};two_high_closed:{all_terminal[1]};two_high_open:{all_terminal[2]};zero_high:{all_terminal[3]};one_high:{all_terminal[4]};cardinality:{all_terminal[5]};maxgap:{all_terminal[6]};failed:{all_terminal[7]};unit_checks:{all_terminal[8]}",
    ]
    for level in LEVELS:
        s = level_screen[level]
        o = level_order[level]
        t = level_terminal[level]
        lines.append(
            f"LEVEL;z1={level};rows={NEIGHBOR_COUNTS[level][0]};wall={NEIGHBOR_COUNTS[level][1]};order={NEIGHBOR_COUNTS[level][2]};"
            f"states={s[0]};crude={s[1]};status={s[2]};residual={s[3]};"
            f"order_states={o[0]};order_crude={o[1]};order_status={o[2]};order_residual={o[3]};"
            f"terminal_bodies={t[0]};two_high_closed={t[1]};zero_high={t[3]};one_high={t[4]};cardinality={t[5]};failed={t[7]};"
            f"row_order_sha256={EXPECTED_ROW_DIGEST[level]}"
        )
    lines.append("GAP;z1=241;rows=0;states=0;mechanism=pinned_atlas_absence")
    for row in records:
        lines.append(
            f"BODY;z1={row[0]};E={row[1]};L={row[2]};high={row[3]};wall={row[5]};"
            f"states={row[9]};crude={row[10]};status={row[11]};residual={row[12]};"
            f"stage_sha256={row[18]};residual_sha256={hashlib.sha256(repr(row[13]).encode()).hexdigest()}"
        )
    for row in terminals:
        lines.append(
            f"TERMINAL;z1={row[0]};E={row[1]};residual={row[5]};gap={row[6]};"
            f"two_high={row[21]};zero_high={row[8]};one_high={row[9]};"
            f"cardinality={row[13]};maxgap={row[14]};failed={row[15]};"
            f"slack={row[17]};closed={row[22]};case_sha256={row[19]}"
        )
    lines += [
        "first_failed_implication=crude_plus_status closes 9412 of 9449 states, not the 37 residuals",
        "repair=all six residual bodies have positive two-high gap and all 37 one-high cases pass strict complete-cell cardinality",
        "safety=no null-set inference;literal covers map to the relaxation and terminal cells retain actual safe residues",
        "conditional_only=if THM3041 is promoted then ledger 374900-72=374828 and projected cap z1<=239",
        "scope=scratch exact candidate only;no theorem id,status,navigation,ledger mutation,or LRC consequence",
        f"semantic_sha256={semantic}",
        "all_exact_controls=PASS",
    ]
    payload = "\n".join(lines) + "\n"
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()



