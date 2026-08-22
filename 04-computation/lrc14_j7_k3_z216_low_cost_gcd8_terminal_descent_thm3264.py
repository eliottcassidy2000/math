#!/usr/bin/env python3
"""Exact referee for the low-cost z216 gcd-eight terminal descent."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import re
from collections import Counter
from fractions import Fraction
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_3139 = ROOT / "04-computation/lrc14_j7_k3_z225_terminal_z224_screen_double_layer_descent_thm3139.py"
OUTPUT_3139 = ROOT / "05-knowledge/results/lrc14_j7_k3_z225_terminal_z224_screen_double_layer_descent_thm3139.out"
SOURCE_3261 = ROOT / "04-computation/lrc14_j7_k3_z216_unique_gcd18_terminal_descent_thm3261.py"
OUTPUT_3261 = ROOT / "05-knowledge/results/lrc14_j7_k3_z216_unique_gcd18_terminal_descent_thm3261.out"
SOURCE_3139_SHA256 = "92e1c22088998d37db89020ac1ebbb7d1ee17e3539152ab205f8b0cd92532e36"
OUTPUT_3139_SHA256 = "4f0ae623d0134406f83f466c0a9f1353525a992acb8ee82c2e92f2f0537f32c5"
SEMANTIC_3139_SHA256 = "1eb7e8faefafe41fd6f6cbc108ba09295dbcbb7fe8016d7ec0cbdc258adf358a"
SOURCE_3261_SHA256 = "308f57b090ea9586c92905d9113b79063ef2105c3e7ed4c90c7406d15400d376"
OUTPUT_3261_SHA256 = "d7dae9cd7e8f0305824b30a2b8683abe6d7d8c9bd0509a0d84026322fd65344c"
SEMANTIC_3261_SHA256 = "b6f83cd9cc405b4e5006ecfbd8ddfb25c115d7a55c270abc55d8d3b2a97b9a18"
ATLAS_SHA256 = "cee82237ce1f51729813b9c916edd3353204c18172abe1d71278dee2c5562eda"

LEVEL = 216
COST_LIMIT = 2_000_000
EXPECTED_CENSUS = (480, 447, 33)
EXPECTED_ROW_SHA256 = "53db9e1d3df2cf2b0398847682d909da81705e43a53ae2553d102fd152337649"
EXPECTED_STRATA = ((8, 19), (18, 1), (24, 135), (36, 15), (72, 310))
EXPECTED_SELECTED = (8, 23, 39, 57, 64, 66, 115, 142, 150, 152, 197, 272, 277, 300, 304, 306, 338)
EXPECTED_EXCLUDED = (238, 370)
EXPECTED_COST = 12_050_080
EXPECTED_TOTALS = (1746, 832, 895, 19)
EXPECTED_FARKAS = (0, 895)
EXPECTED_SCREEN_SHA256 = "7ad0a1a8036fefa37a9035d2d2b319206b170dcbea2f7f86b7cef96cd7c1b3f1"
EXPECTED_RESIDUAL_SHA256 = "ae75d7b0f6947dbc51dfacaa8b34dea9d2c985b04014eb09ac923ae720a0452d"
EXPECTED_RESIDUAL = (
    (64, (1, 2, 8, 10, 11, 14), 8),
    (152, (1, 4, 8, 10, 13, 14), 9),
    (197, (1, 5, 8, 10, 13, 14), 2),
)
EXPECTED_GAPS = (
    (64, Fraction(945034906, 306786044565)),
    (152, Fraction(35479656293, 19991919811500)),
    (197, Fraction(368677303, 146893907160)),
)
EXPECTED_TERMINAL_AGGREGATE = (3, 19, 19, 28, 12, 28, 0, 0, 0, 5, 3, 3)
EXPECTED_TERMINAL_SEMANTIC_SHA256 = "9328d556544c230503b6860ea27958d71bfa3f8af7816d02d350da360e92acac"
LEDGER_BEFORE = 373283
LEDGER_AFTER = 373266


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_sha(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n")
    require(b"\r" not in payload, (path, "bare CR"))
    return hashlib.sha256(payload).hexdigest()


def check_dependency(path, digest, semantic=None):
    require(lf_sha(path) == digest, (path, "hash changed"))
    if semantic is not None:
        require(f"semantic_sha256={semantic}" in path.read_text(encoding="utf-8"),
                (path, "semantic changed"))


def load(name, path=SOURCE_3139):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def terminal_semantic(row):
    return (*row[:6], *row[7:])


def ftext(value):
    return f"{value.numerator}/{value.denominator}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    check_dependency(SOURCE_3139, SOURCE_3139_SHA256)
    check_dependency(OUTPUT_3139, OUTPUT_3139_SHA256, SEMANTIC_3139_SHA256)
    check_dependency(SOURCE_3261, SOURCE_3261_SHA256)
    check_dependency(OUTPUT_3261, OUTPUT_3261_SHA256, SEMANTIC_3261_SHA256)

    entry = load("thm3264_entry")
    source = entry.load("thm3264_source")
    driver = source.load("thm3264_driver")
    predecessor = driver.load("thm3264_predecessor")
    bridge = predecessor.load("thm3264_bridge")
    base = bridge.load("thm3264_base")
    require(base.thm.ATLAS_SHA256 == ATLAS_SHA256, base.thm.ATLAS_SHA256)
    rows = entry.atlas_rows(base, LEVEL)
    census = (len(rows), sum(row[3] for row in rows), sum(not row[3] for row in rows))
    require(census == EXPECTED_CENSUS, census)
    require(hashlib.sha256(repr(rows).encode()).hexdigest() == EXPECTED_ROW_SHA256,
            "row order")
    strata = Counter(gcd(LEVEL, row[1]) for row in rows)
    require(tuple(sorted(strata.items())) == EXPECTED_STRATA, strata)

    pattern = re.compile(
        r"^row=E=([0-9,]+);h=[^;]+;r=([0-9]+);L=([0-9]+);high=([0-9]+);z1=([0-9]+);"
    )
    components = []
    for line in base.thm.ATLAS.read_text(encoding="utf-8").splitlines():
        if not line.startswith("row="):
            continue
        match = pattern.match(line)
        require(match is not None, line)
        if int(match.group(5)) == LEVEL:
            components.append(int(match.group(2)))
    require(len(components) == len(rows), len(components))
    g8 = tuple(i for i, row in enumerate(rows) if gcd(LEVEL, row[1]) == 8)
    selected = tuple(i for i in g8 if rows[i][1] * components[i] <= COST_LIMIT)
    excluded = tuple(i for i in g8 if i not in selected)
    require(selected == EXPECTED_SELECTED and excluded == EXPECTED_EXCLUDED,
            (selected, excluded))
    cost = sum(rows[i][1] * components[i] for i in selected)
    require(cost == EXPECTED_COST, cost)

    canonical = []
    farkas = [0, 0]
    for index in selected:
        task = (LEVEL, *rows[index])
        returned, full = entry.screen_worker(task)
        require(returned == task, (index, returned, task))
        row = tuple(full[:19])
        require(row[:6] == (LEVEL, rows[index][0], rows[index][1], rows[index][2],
                            rows[index][1] // 8, rows[index][3]), (index, row[:6]))
        require(row[16] == row[11], (index, "status verification"))
        canonical.append(row)
        farkas[0] += full[19]
        farkas[1] += full[20]
    canonical = tuple(canonical)
    totals = tuple(sum(row[i] for row in canonical) for i in (9, 10, 11, 12))
    require(totals == EXPECTED_TOTALS and tuple(farkas) == EXPECTED_FARKAS,
            (totals, farkas))
    require(hashlib.sha256(repr(canonical).encode()).hexdigest() == EXPECTED_SCREEN_SHA256,
            "screen digest")
    residual = tuple((row[1], row[13]) for row in canonical if row[12])
    require(hashlib.sha256(repr(residual).encode()).hexdigest() == EXPECTED_RESIDUAL_SHA256,
            "residual digest")
    residual_rows = tuple(
        (selected[pos], row[1], len(row[13]))
        for pos, row in enumerate(canonical) if row[12]
    )
    require(residual_rows == EXPECTED_RESIDUAL, residual_rows)
    quotient = Counter()
    for row in canonical:
        for ds in row[13]:
            divisor = lcm(*ds)
            require(divisor % row[4] == 0 and row[2] % divisor == 0,
                    (row[1], ds, divisor))
            quotient[divisor // row[4]] += 1
    require(tuple(sorted(quotient.items())) == ((8, 19),), quotient)

    terminals = []
    for index, body, _masks in EXPECTED_RESIDUAL:
        pos = selected.index(index)
        terminal = base.thm.terminal_probe((LEVEL, body, canonical[pos][13]))
        terminals.append((index, terminal))
    gaps = tuple((index, row[5]) for index, row in terminals)
    require(gaps == EXPECTED_GAPS, gaps)
    aggregate = (
        len(terminals), sum(row[4] for _i, row in terminals),
        sum(row[7] for _i, row in terminals), sum(row[8] for _i, row in terminals),
        sum(row[9] for _i, row in terminals), sum(row[10] for _i, row in terminals),
        sum(row[11] for _i, row in terminals), sum(row[12] for _i, row in terminals),
        sum(row[14] for _i, row in terminals), min(row[16] for _i, row in terminals),
        sum(bool(row[19]) for _i, row in terminals),
        sum(bool(row[20]) for _i, row in terminals),
    )
    require(aggregate == EXPECTED_TERMINAL_AGGREGATE, aggregate)
    terminal_packet = tuple((index, terminal_semantic(row)) for index, row in terminals)
    terminal_sha = hashlib.sha256(repr(terminal_packet).encode()).hexdigest()
    require(terminal_sha == EXPECTED_TERMINAL_SEMANTIC_SHA256, terminal_sha)
    require(LEDGER_BEFORE - len(selected) == LEDGER_AFTER, "ledger arithmetic")

    semantic_packet = (
        "lrc14-k3-z216-low-cost-g8-terminal-v1",
        (SOURCE_3139_SHA256, OUTPUT_3139_SHA256, SEMANTIC_3139_SHA256,
         SOURCE_3261_SHA256, OUTPUT_3261_SHA256, SEMANTIC_3261_SHA256,
         ATLAS_SHA256),
        (LEVEL, COST_LIMIT, EXPECTED_CENSUS, EXPECTED_ROW_SHA256, EXPECTED_STRATA,
         EXPECTED_SELECTED, EXPECTED_EXCLUDED, EXPECTED_COST, EXPECTED_TOTALS,
         EXPECTED_FARKAS, EXPECTED_SCREEN_SHA256, EXPECTED_RESIDUAL_SHA256,
         EXPECTED_RESIDUAL, ((8, 19),)),
        (EXPECTED_GAPS, EXPECTED_TERMINAL_AGGREGATE,
         EXPECTED_TERMINAL_SEMANTIC_SHA256),
        (LEDGER_BEFORE, LEDGER_AFTER, 216, 462, 429, 33),
    )
    semantic = hashlib.sha256(repr(semantic_packet).encode()).hexdigest()
    lines = [
        "LRC14 THM3264 projected-k3 z216 low-cost gcd8 terminal descent",
        f"engine=THM3139_source:{SOURCE_3139_SHA256};output:{OUTPUT_3139_SHA256};semantic:{SEMANTIC_3139_SHA256}",
        f"dependency=THM3261_source:{SOURCE_3261_SHA256};output:{OUTPUT_3261_SHA256};semantic:{SEMANTIC_3261_SHA256}",
        f"atlas=sha256:{ATLAS_SHA256};rows:6060",
        f"layer=z1:216;rows:480;wall:447;order:33;row_sha256:{EXPECTED_ROW_SHA256};strata:{EXPECTED_STRATA}",
        f"selection=gcd:8;cost_limit:{COST_LIMIT};selected:{EXPECTED_SELECTED};excluded:{EXPECTED_EXCLUDED};cost:{cost}",
        f"screen=rows:{len(selected)};states:{totals[0]};crude:{totals[1]};status:{totals[2]};residual:{totals[3]};direct_farkas:{farkas[0]};legacy_farkas:{farkas[1]};screen_sha256:{EXPECTED_SCREEN_SHA256}",
        f"residual=bodies:{EXPECTED_RESIDUAL};masks:19;residual_sha256:{EXPECTED_RESIDUAL_SHA256};quotient:((8,19),)",
        "terminal=rows:3;masks:19;zero_high_hostiles:19;one_high_cases:28;low_label_sets:12;coarse:28;exact:0;maxgap:0;failures:0;minimum_slack:5;positive_two_high_gap:3;closed:3",
        f"gaps:{tuple((i, ftext(g)) for i,g in EXPECTED_GAPS)};semantic_sha256:{terminal_sha}",
        "direction=the_exact_screens_and_complete_one_high_terminal_banks_are_upper_relaxations;their_exclusion_closes_all_17_selected_rows",
        "modular_boundary=every_residual_is_top_q8_in_a_four_point_divisor_chain_but_no_iterated_C2_action_common_carrier_or_parity_tower_is_constructed",
        f"consequence=ledger:{LEDGER_BEFORE}-17={LEDGER_AFTER};projected_k3_cap:z1<=216;remaining_z216:462=429wall+33order",
        "scope=17_low_cost_gcd8_rows_only;the_two_high_cost_gcd8_rows_and_all_other_z216_strata_remain_open;no_physical_cover_k_le_1_rung_or_LRC14_claim",
        f"semantic_sha256={semantic}",
        "all_exact_controls=PASS",
    ]
    text = "\n".join(lines) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8", newline="\n")
    print(text, end="")


if __name__ == "__main__":
    main()
