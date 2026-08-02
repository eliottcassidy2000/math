#!/usr/bin/env python3
"""Exact projected-k3 compositor for occupied levels 236 and 235.

This is a typed continuation of THM-3061 inside the same pinned lossless
projected necessary atlas.  It closes exactly the one occupied row at
``z_1=236`` and the eleven occupied rows at ``z_1=235``.  It does not inspect
or close the next occupied layer ``z_1=234`` and has no physical-cover or
LRC(14) consequence.
"""

from __future__ import annotations

import argparse
import hashlib
import re
import sys
import tempfile
from pathlib import Path


HERE = Path(__file__).resolve()
ROOT = next(parent for parent in HERE.parents if (parent / "00-navigation").is_dir())
BASE_SOURCE = (
    ROOT
    / "04-computation/lrc14_j7_k3_z239_gap238_z237_compositional_descent_thm3061.py"
)
BASE_OUTPUT = (
    ROOT
    / "05-knowledge/results/lrc14_j7_k3_z239_gap238_z237_compositional_descent_thm3061.out"
)
OUTPUT = (
    ROOT
    / "05-knowledge/results/lrc14_j7_k3_z236_z235_compositional_descent_thm3071.out"
)

BASE_SOURCE_SHA256 = "93a90c48ebed4bc782bca31f378ebb0d4f7ee19ef471a60220fcda3e8927e2fb"
BASE_OUTPUT_SHA256 = "b30961378985a543a86a82681f8556353243703aca9bacf09bb6b4f61648274c"
BASE_SEMANTIC_SHA256 = "2b825c5f12a59048b497b73dba234b79301d78e2db65325c63a580870eaccc88"
ATLAS_SHA256 = "cee82237ce1f51729813b9c916edd3353204c18172abe1d71278dee2c5562eda"

LEVELS = (236, 235)
NEIGHBOR_COUNTS = {
    239: (4, 2, 2),
    238: (0, 0, 0),
    237: (44, 33, 11),
    236: (1, 1, 0),
    235: (11, 5, 6),
    234: (381, 330, 51),
}
EXPECTED_SCREEN = {
    236: (2, 2, 0, 0),
    235: (577, 99, 338, 140),
}
EXPECTED_ORDER = {
    236: (0, 0, 0, 0),
    235: (27, 9, 18, 0),
}
EXPECTED_TERMINAL = {
    236: (0, 0, 0, 0, 0, 0, 0, 0, 0),
    235: (4, 4, 4, 138, 140, 140, 0, 0, 0),
}
EXPECTED_ROW_DIGEST = {
    236: "d37a864a1cb2c8a0116583330fb8a45fa71d02ede2b7190d5b393cd86301184e",
    235: "40eb51949079b203ba7552e272741a4dea6e09c38eb5fe5114e26a50e3ba4766",
    "combined": "560fd7ba77b5170e234e2cde486438583791b1ab7696a516ffe23204e980349e",
}
EXPECTED_SEMANTIC_SHA256 = "d4cfb4c7e497007f044041e9f7e56d3b109107d033421d592f7f3f74e0995f02"
CACHE_SCHEMA = "lrc14-k3-z236-z235-combined-v1"

PREDECESSOR_LEDGER = 374780
DECREMENT = 12
UPDATED_LEDGER = 374768
UPDATED_CAP = 234


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    payload = path.read_bytes()
    require(b"\r" not in payload.replace(b"\r\n", b""), f"bare CR in {path}")
    return hashlib.sha256(payload.replace(b"\r\n", b"\n")).hexdigest()


require(sha(BASE_SOURCE) == BASE_SOURCE_SHA256, "THM-3061 source changed")
require(sha(BASE_OUTPUT) == BASE_OUTPUT_SHA256, "THM-3061 output changed")
base_output_text = BASE_OUTPUT.read_text()
require(
    f"semantic_sha256={BASE_SEMANTIC_SHA256}" in base_output_text,
    "THM-3061 semantic changed",
)
require(
    "promotion_consequence=ledger 374828-48=374780;projected cap z1<=236"
    in base_output_text,
    "THM-3061 ledger consequence changed",
)
sys.path.insert(0, str(ROOT / "04-computation"))
import lrc14_j7_k3_z239_gap238_z237_compositional_descent_thm3061 as base

require(sha(base.base.ATLAS) == ATLAS_SHA256, "THM-2941 projected body atlas changed")


def atlas_tasks():
    pattern = re.compile(
        r"^row=E=([0-9,]+);.*;L=([0-9]+);high=([0-9]+);z1=([0-9]+);"
    )
    row_lines = tuple(
        line
        for line in base.base.ATLAS.read_text().splitlines()
        if line.startswith("row=")
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
            == base.base.thm.base.WALL.numerator
            * L
            // base.base.thm.base.WALL.denominator
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
    require(len(tasks) == DECREMENT, ("combined occupied rows", len(tasks)))
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
    require(NEIGHBOR_COUNTS[237][0] == 44, "THM-3061 predecessor row count changed")
    require(
        NEIGHBOR_COUNTS[236][0] + NEIGHBOR_COUNTS[235][0] == DECREMENT,
        "additive set difference changed",
    )
    require(NEIGHBOR_COUNTS[234][0] > 0, "next occupied layer changed")
    require(PREDECESSOR_LEDGER - DECREMENT == UPDATED_LEDGER, "ledger arithmetic")
    return tuple(tasks), derived


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
                PREDECESSOR_LEDGER,
                DECREMENT,
                UPDATED_LEDGER,
                UPDATED_CAP,
            )
        ).encode()
    ).hexdigest()


CACHE_FINGERPRINT = checkpoint_fingerprint()
# Reuse THM-3061's audited checkpoint machinery while separating this cache
# bank by the stronger fingerprint above.
base.CACHE_FINGERPRINT = CACHE_FINGERPRINT


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=4)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument(
        "--checkpoint-dir",
        type=Path,
        default=Path(tempfile.gettempdir()) / "lrc14-k3-z236-z235-combined",
    )
    args = parser.parse_args()

    tasks, neighbor_census = atlas_tasks()
    scheduled = tuple(sorted(tasks, key=lambda task: (-task[2], -task[0], task[1])))
    records = base.run_pool(
        base.worker_evaluate,
        scheduled,
        args.processes,
        "screen",
        args.checkpoint_dir,
    )
    by_key = {(row[0], row[1]): row for row in records}
    require(len(by_key) == len(tasks), "lost screen row")
    records = tuple(by_key[(task[0], task[1])] for task in tasks)

    terminal_tasks = tuple(
        (row[0], row[1], row[13]) for row in records if row[5] and row[12]
    )
    terminals = base.run_pool(
        base.worker_terminal,
        terminal_tasks,
        args.processes,
        "terminal",
        args.checkpoint_dir,
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
        level_screen[level] = base.four_totals(layer_records)
        level_order[level] = base.four_totals(
            tuple(row for row in layer_records if not row[5])
        )
        level_terminal[level] = base.terminal_totals(layer_terminals)
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

    totals = base.four_totals(records)
    order = base.four_totals(tuple(row for row in records if not row[5]))
    all_terminal = base.terminal_totals(terminals)
    require(totals == (579, 101, 338, 140), ("combined screen", totals))
    require(order == (27, 9, 18, 0), ("combined order", order))
    require(
        all_terminal == (4, 4, 4, 138, 140, 140, 0, 0, 0),
        ("combined terminal", all_terminal),
    )
    require(all(row[6] > 0 for row in terminals), "nonpositive two-high gap")
    require(min(row[17] for row in terminals) == 1, "minimum slack changed")
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
        (PREDECESSOR_LEDGER, DECREMENT, UPDATED_LEDGER, UPDATED_CAP),
    )
    semantic = hashlib.sha256(repr(semantic_packet).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, ("semantic digest", semantic))

    lines = [
        "LRC14 projected k3 z1=236 and z1=235 compositional descent",
        f"dependency=THM3061_source:{BASE_SOURCE_SHA256};output:{BASE_OUTPUT_SHA256};semantic:{BASE_SEMANTIC_SHA256};atlas:{ATLAS_SHA256}",
        f"checkpoint_fingerprint={CACHE_FINGERPRINT}",
        "neighbor_census="
        + ";".join(
            f"z{z}:rows{values[0]}/wall{values[1]}/order{values[2]}"
            for z, values in sorted(neighbor_census.items(), reverse=True)
        ),
        f"universe=occupied_rows:{len(tasks)};wall:{sum(t[4] for t in tasks)};order:{sum(not t[4] for t in tasks)};row_order_sha256:{EXPECTED_ROW_DIGEST['combined']}",
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
    lines.append(
        "BOUNDARY;z1=234;rows=381;wall=330;order=51;status=NEXT_OCCUPIED_LAYER_NOT_INCLUDED"
    )
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
        "first_failed_implication=crude_plus_status closes 439 of 579 states, not the 140 residuals",
        "repair=all 4 residual bodies have positive two-high gap and all 140 one-high cases pass strict complete-cell cardinality",
        "safety=no null-set inference;literal covers map to the inherited relaxation and terminal cells retain actual safe residues",
        "set_difference=THM3061 closed z239:4 and z237:44;current disjoint pinned layers z236:1 and z235:11",
        f"promotion_consequence=ledger {PREDECESSOR_LEDGER}-{DECREMENT}={UPDATED_LEDGER};projected cap z1<={UPDATED_CAP}",
        "scope=projected necessary sector only;z234 is not inspected or closed;no physical-cover or LRC14 consequence",
        f"source_sha256={sha(HERE)}",
        f"semantic_sha256={semantic}",
        "all_exact_controls=PASS",
    ]
    payload = "\n".join(lines) + "\n"
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
