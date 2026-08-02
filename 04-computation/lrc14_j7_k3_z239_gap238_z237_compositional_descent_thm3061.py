#!/usr/bin/env python3
"""Exact projected-k3 compositor for levels 239, 238, and 237."""

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
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE_SOURCE = (
    ROOT
    / "04-computation/lrc14_j7_k3_z242_gap241_z240_compositional_descent_thm3052.py"
)
BASE_OUTPUT = (
    ROOT
    / "05-knowledge/results/lrc14_j7_k3_z242_gap241_z240_compositional_descent_thm3052.out"
)
OUTPUT = (
    ROOT
    / "05-knowledge/results/lrc14_j7_k3_z239_gap238_z237_compositional_descent_thm3061.out"
)
BASE_SOURCE_SHA256 = "41d40df237d34ca92e06bb55eed022534810866cce9759d225f1031b49f3f2fb"
BASE_OUTPUT_SHA256 = "61a7a18531128226a967755dfe60c0012ca993adae27542afc9d47c81d7feb7b"
BASE_SEMANTIC_SHA256 = "8e2e6a1ed90033f3d85ec9cd03acf326ffb0e67453bbc206e8f2ece843c7cd68"
ATLAS_SHA256 = "cee82237ce1f51729813b9c916edd3353204c18172abe1d71278dee2c5562eda"
LEVELS = (239, 237)
GAP_LEVEL = 238
NEIGHBOR_COUNTS = {
    240: (52, 49, 3),
    239: (4, 2, 2),
    238: (0, 0, 0),
    237: (44, 33, 11),
    236: (1, 1, 0),
}
EXPECTED_SCREEN = {
    239: (13, 1, 12, 0),
    237: (4005, 1822, 1813, 370),
}
EXPECTED_ORDER = {
    239: (3, 0, 3, 0),
    237: (79, 71, 8, 0),
}
EXPECTED_TERMINAL = {
    239: (0, 0, 0, 0, 0, 0, 0, 0, 0),
    237: (14, 14, 14, 341, 377, 377, 0, 0, 0),
}
EXPECTED_ROW_DIGEST = {
    239: "0fd8cd4faf90222c7a7f92bcfbb483bded78a0165e4f41dfe175de106cfcc015",
    237: "e7760ae55940b1746d7d124a0b830b598ec001ec6da11fb652f8fff41c56ed1a",
    "combined": "76af1f045822eb75d507e98ecaae9eaf0b65c67ae6a1b8ce8432c9be40d05787",
}
EXPECTED_SEMANTIC_SHA256 = "2b825c5f12a59048b497b73dba234b79301d78e2db65325c63a580870eaccc88"
CACHE_SCHEMA = "lrc14-k3-z239-gap238-z237-combined-v1"


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


require(sha(BASE_SOURCE) == BASE_SOURCE_SHA256, "THM-3052 source changed")
require(sha(BASE_OUTPUT) == BASE_OUTPUT_SHA256, "THM-3052 output changed")
require(
    f"semantic_sha256={BASE_SEMANTIC_SHA256}" in BASE_OUTPUT.read_text(),
    "THM-3052 semantic changed",
)
base = load("thm3052_z239_gap238_z237_base", BASE_SOURCE)
require(sha(base.ATLAS) == ATLAS_SHA256, "THM-2941 projected body atlas changed")


def atlas_tasks():
    pattern = re.compile(
        r"^row=E=([0-9,]+);.*;L=([0-9]+);high=([0-9]+);z1=([0-9]+);"
    )
    row_lines = tuple(
        line for line in base.ATLAS.read_text().splitlines() if line.startswith("row=")
    )
    require(len(row_lines) == 6060, ("atlas row-line count", len(row_lines)))
    rows = {}
    tasks = []
    for line in row_lines:
        match = pattern.match(line)
        require(match is not None, ("unparsed atlas row", line))
        body = tuple(map(int, match.group(1).split(",")))
        L, high, first = map(int, match.group(2, 3, 4))
        require(
            high
            == base.thm.base.WALL.numerator
            * L
            // base.thm.base.WALL.denominator
            + 1,
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
        z: sum(f";z1={z};" in line for line in row_lines) for z in NEIGHBOR_COUNTS
    }
    require(
        literal == {z: values[0] for z, values in NEIGHBOR_COUNTS.items()},
        ("literal neighbor census", literal),
    )
    require(GAP_LEVEL not in rows, "z1=238 unexpectedly occupied")
    require(len(tasks) == 48, ("combined occupied rows", len(tasks)))
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
    return task, base.thm.evaluate(task)


def worker_terminal(task):
    return task, base.terminal_probe(task)


def checkpoint_fingerprint():
    return hashlib.sha256(
        repr(
            (
                CACHE_SCHEMA,
                BASE_SOURCE_SHA256,
                BASE_OUTPUT_SHA256,
                BASE_SEMANTIC_SHA256,
                ATLAS_SHA256,
                NEIGHBOR_COUNTS,
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
    if not tasks:
        return ()
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
    print(
        f"{label}_recovered={len(rows)};pending={len(pending)}",
        file=sys.stderr,
        flush=True,
    )

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
    elif pending:
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
        sum(row[22] for row in rows),
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
        default=Path(tempfile.gettempdir()) / "lrc14-k3-z239-gap238-z237-combined",
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
    terminal_by_key = {(row[0], row[1]): row for row in terminals}
    require(len(terminal_by_key) == len(terminal_tasks), "lost terminal row")
    terminals = tuple(
        terminal_by_key[(task[0], task[1])] for task in terminal_tasks
    )

    level_screen = {}
    level_order = {}
    level_terminal = {}
    for level in LEVELS:
        layer_records = tuple(row for row in records if row[0] == level)
        layer_terminals = tuple(row for row in terminals if row[0] == level)
        level_screen[level] = four_totals(layer_records)
        level_order[level] = four_totals(
            tuple(row for row in layer_records if not row[5])
        )
        level_terminal[level] = terminal_totals(layer_terminals)
        require(
            level_screen[level] == EXPECTED_SCREEN[level],
            (level, "screen", level_screen[level]),
        )
        require(
            level_order[level] == EXPECTED_ORDER[level],
            (level, "order", level_order[level]),
        )
        require(
            level_terminal[level] == EXPECTED_TERMINAL[level],
            (level, "terminal", level_terminal[level]),
        )

    totals = four_totals(records)
    order = four_totals(tuple(row for row in records if not row[5]))
    all_terminal = terminal_totals(terminals)
    require(totals == (4018, 1823, 1825, 370), ("combined screen", totals))
    require(order == (82, 71, 11, 0), ("combined order", order))
    require(
        all_terminal == (14, 14, 14, 341, 377, 377, 0, 0, 0),
        ("combined terminal", all_terminal),
    )
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
        require(
            semantic == EXPECTED_SEMANTIC_SHA256,
            ("semantic digest", semantic),
        )

    lines = [
        "LRC14 projected k3 z1=239, gap238, z1=237 combined exact scout",
        f"dependency=THM3052_source:{BASE_SOURCE_SHA256};output:{BASE_OUTPUT_SHA256};semantic:{BASE_SEMANTIC_SHA256};atlas:{ATLAS_SHA256}",
        f"checkpoint_fingerprint={CACHE_FINGERPRINT}",
        "neighbor_census="
        + ";".join(
            f"z{z}:rows{values[0]}/wall{values[1]}/order{values[2]}"
            for z, values in sorted(neighbor_census.items(), reverse=True)
        ),
        f"universe=occupied_rows:{len(tasks)};wall:{sum(t[4] for t in tasks)};order:{sum(not t[4] for t in tasks)};gap_z238_rows:0;row_order_sha256:{EXPECTED_ROW_DIGEST['combined']}",
        f"screen=states:{totals[0]};crude:{totals[1]};status:{totals[2]};residual:{totals[3]};residual_bodies:{len(terminal_tasks)}",
        f"order_screen=states:{order[0]};crude:{order[1]};status:{order[2]};residual:{order[3]}",
        f"terminal=residual_bodies:{all_terminal[0]};two_high_closed:{all_terminal[1]};closed:{all_terminal[2]};zero_high:{all_terminal[3]};one_high:{all_terminal[4]};cardinality:{all_terminal[5]};maxgap:{all_terminal[6]};failed:{all_terminal[7]};unit_checks:{all_terminal[8]}",
    ]
    for level in LEVELS:
        s = level_screen[level]
        o = level_order[level]
        t = level_terminal[level]
        lines.append(
            f"LEVEL;z1={level};rows={NEIGHBOR_COUNTS[level][0]};wall={NEIGHBOR_COUNTS[level][1]};order={NEIGHBOR_COUNTS[level][2]};"
            f"states={s[0]};crude={s[1]};status={s[2]};residual={s[3]};"
            f"order_states={o[0]};order_crude={o[1]};order_status={o[2]};order_residual={o[3]};"
            f"terminal_bodies={t[0]};two_high_closed={t[1]};closed={t[2]};zero_high={t[3]};one_high={t[4]};cardinality={t[5]};failed={t[7]};"
            f"row_order_sha256={EXPECTED_ROW_DIGEST[level]}"
        )
    lines.append("GAP;z1=238;rows=0;states=0;mechanism=pinned_atlas_absence")
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
        "first_failed_implication=crude_plus_status closes 3648 of 4018 states, not the 370 residuals",
        "repair=all 14 residual bodies have positive two-high gap and all 377 one-high cases pass strict complete-cell cardinality",
        "safety=no null-set inference;literal covers map to the relaxation and terminal cells retain actual safe residues",
        "promotion_consequence=ledger 374828-48=374780;projected cap z1<=236",
        "scope=projected necessary sector only;this computation makes no navigation/ledger mutation or LRC consequence",
        f"semantic_sha256={semantic}",
        "all_exact_controls=PASS",
    ]
    payload = "\n".join(lines) + "\n"
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
