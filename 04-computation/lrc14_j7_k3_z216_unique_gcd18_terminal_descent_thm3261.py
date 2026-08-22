#!/usr/bin/env python3
"""Canonical candidate: unique z216 gcd-18 screen and terminal closure."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
from collections import Counter
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_3139 = ROOT / "04-computation/lrc14_j7_k3_z225_terminal_z224_screen_double_layer_descent_thm3139.py"
OUTPUT_3139 = ROOT / "05-knowledge/results/lrc14_j7_k3_z225_terminal_z224_screen_double_layer_descent_thm3139.out"
SOURCE_3251 = ROOT / "04-computation/lrc14_j7_k3_z218_terminal_descent_composed_cap216_thm3251.py"
OUTPUT_3251 = ROOT / "05-knowledge/results/lrc14_j7_k3_z218_terminal_descent_composed_cap216_thm3251.out"
SOURCE_3139_SHA256 = "92e1c22088998d37db89020ac1ebbb7d1ee17e3539152ab205f8b0cd92532e36"
OUTPUT_3139_SHA256 = "4f0ae623d0134406f83f466c0a9f1353525a992acb8ee82c2e92f2f0537f32c5"
SEMANTIC_3139_SHA256 = "1eb7e8faefafe41fd6f6cbc108ba09295dbcbb7fe8016d7ec0cbdc258adf358a"
SOURCE_3251_SHA256 = "d92a6b825268d8fa7147ebcd5229bf8fa1d4cd2ffa3e9c1ec140024bd926832f"
OUTPUT_3251_SHA256 = "66f7f6bb81992659d5cedf01aa2f22b12872810888a022460eb6a5ce1c05d674"
SEMANTIC_3251_SHA256 = "fe1f0e0df1c652be79b1e9d1e964bc4db64f1930e25b8fb8a6f7d3c5f47ccedb"
ATLAS_SHA256 = "cee82237ce1f51729813b9c916edd3353204c18172abe1d71278dee2c5562eda"

LEVEL = 216
EXPECTED_CENSUS = (480, 447, 33)
EXPECTED_ROW_SHA256 = "53db9e1d3df2cf2b0398847682d909da81705e43a53ae2553d102fd152337649"
EXPECTED_STRATA = ((8, 19), (18, 1), (24, 135), (36, 15), (72, 310))
ATLAS_INDEX = 185
BODY = (1, 5, 7, 9, 11, 13)
RULER = 630630
HIGH = 62108
EXPECTED_TOTALS = (143, 0, 113, 30)
EXPECTED_FARKAS = (0, 113)
EXPECTED_SCREEN_SHA256 = "5e2721f49b5f1d9f7dadd156fa7277cccb3a0693b648c2375d714bf47765c4c5"
EXPECTED_ROW_SEMANTIC_SHA256 = "33ba863d05950923a9a10719f4c17ab4c4b6df77cf8ca3f355b61618b3250d25"
EXPECTED_ROW_FAILURE_SHA256 = "e45b0879387d0df2907a0d17e021659a1f4eafd63959d040cae3d8ca78d7a1e6"
EXPECTED_MASK_SHA256 = "2d753b2dc75f3ba79d3ae0b3f79999669807f6906562ff5072c71bb0ed2fd45e"
EXPECTED_RESIDUAL_SHA256 = "652e914b0e9e11ec997c8f168286e6684864770d42b04b94676c4299ddecb812"
EXPECTED_QUOTIENT = (18, 18, 30)
EXPECTED_TERMINAL_SUMMARY = (1, 30, 30, 30, 1, 30, 0, 0, 0, 8, 1, 1)
EXPECTED_TERMINAL_SEMANTIC_SHA256 = "607e7cd3166f078ace15835d71dd3f7f0e732624c957a9273c3f7a277b475e0f"
EXPECTED_CLOSURE_SHA256 = "88af81e3326b00bceb14c280d2557c22210ba9c33050074814baeaaec9e340af"
EXPECTED_CASE_SHA256 = "4dd6438ccc4109f81388a4eb8ee50015e76b6b9050e1521c00763167408b0148"
EXPECTED_GAP = (2739297996619, 899643669274350)
LEDGER_BEFORE = 373284
LEDGER_AFTER = 373283


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_sha(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n")
    require(b"\r" not in payload, (path, "bare CR"))
    return hashlib.sha256(payload).hexdigest()


def load(name, path=SOURCE_3139):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def check_dependency(path, digest, semantic=None):
    require(lf_sha(path) == digest, (path, "hash changed"))
    if semantic is not None:
        require(f"semantic_sha256={semantic}" in path.read_text(encoding="utf-8"),
                (path, "semantic changed"))


def terminal_semantic(row):
    return (*row[:6], *row[7:])


def closure_row(row):
    return (
        row[1], row[4], (row[5].numerator, row[5].denominator), row[7], row[8],
        row[9], row[10], row[11], row[12], row[14], row[16], row[17], row[19], row[20],
    )


def ftext(value):
    return f"{value.numerator}/{value.denominator}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    check_dependency(SOURCE_3139, SOURCE_3139_SHA256)
    check_dependency(OUTPUT_3139, OUTPUT_3139_SHA256, SEMANTIC_3139_SHA256)
    check_dependency(SOURCE_3251, SOURCE_3251_SHA256)
    check_dependency(OUTPUT_3251, OUTPUT_3251_SHA256, SEMANTIC_3251_SHA256)

    entry = load("thm3261_entry")
    source = entry.load("thm3261_source")
    driver = source.load("thm3261_driver")
    predecessor = driver.load("thm3261_predecessor")
    bridge = predecessor.load("thm3261_bridge")
    base = bridge.load("thm3261_base")
    require(base.thm.ATLAS_SHA256 == ATLAS_SHA256, base.thm.ATLAS_SHA256)
    rows = entry.atlas_rows(base, LEVEL)
    census = (len(rows), sum(row[3] for row in rows), sum(not row[3] for row in rows))
    require(census == EXPECTED_CENSUS, census)
    require(hashlib.sha256(repr(rows).encode()).hexdigest() == EXPECTED_ROW_SHA256, "row order")
    strata = Counter(gcd(LEVEL, row[1]) for row in rows)
    require(tuple(sorted(strata.items())) == EXPECTED_STRATA, strata)
    selected = tuple((i, row) for i, row in enumerate(rows) if gcd(LEVEL, row[1]) == 18)
    require(selected == ((ATLAS_INDEX, (BODY, RULER, HIGH, True)),), selected)

    task = (LEVEL, BODY, RULER, HIGH, True)
    returned_task, full_screen = entry.screen_worker(task)
    require(returned_task == task, (returned_task, task))
    screen = tuple(full_screen[:19])
    require(screen[:6] == (LEVEL, BODY, RULER, HIGH, RULER // 18, True), screen[:6])
    require(screen[16] == screen[11], "status verification")
    totals = tuple(screen[i] for i in (9, 10, 11, 12))
    farkas = (full_screen[19], full_screen[20])
    require(totals == EXPECTED_TOTALS and farkas == EXPECTED_FARKAS, (totals, farkas))
    require(hashlib.sha256(repr((screen,)).encode()).hexdigest() == EXPECTED_SCREEN_SHA256,
            "screen digest")
    require(screen[15] == EXPECTED_ROW_SEMANTIC_SHA256, screen[15])
    require(screen[18] == EXPECTED_ROW_FAILURE_SHA256, screen[18])
    require(hashlib.sha256(repr(screen[13]).encode()).hexdigest() == EXPECTED_MASK_SHA256,
            "mask digest")
    residual = ((BODY, screen[13]),)
    require(hashlib.sha256(repr(residual).encode()).hexdigest() == EXPECTED_RESIDUAL_SHA256,
            "residual digest")
    quotient = Counter()
    for ds in screen[13]:
        divisor = lcm(*ds)
        require(divisor % screen[4] == 0 and RULER % divisor == 0, ds)
        quotient[(gcd(LEVEL, RULER), divisor // screen[4])] += 1
    quotient_packet = tuple((g, q, count) for (g, q), count in sorted(quotient.items()))
    require(quotient_packet == (EXPECTED_QUOTIENT,), quotient_packet)

    terminal = base.thm.terminal_probe((LEVEL, BODY, screen[13]))
    summary = (
        1, terminal[4], terminal[7], terminal[8], terminal[9], terminal[10], terminal[11],
        terminal[12], terminal[14], terminal[16], int(terminal[19]), int(terminal[20]),
    )
    require(summary == EXPECTED_TERMINAL_SUMMARY, summary)
    require((terminal[5].numerator, terminal[5].denominator) == EXPECTED_GAP, terminal[5])
    terminal_semantic_sha = hashlib.sha256(
        repr((terminal_semantic(terminal),)).encode()
    ).hexdigest()
    require(terminal_semantic_sha == EXPECTED_TERMINAL_SEMANTIC_SHA256,
            terminal_semantic_sha)
    closure_sha = hashlib.sha256(
        repr(("z216-g18-terminal-v1", (closure_row(terminal),))).encode()
    ).hexdigest()
    require(closure_sha == EXPECTED_CLOSURE_SHA256, closure_sha)
    case_sha = hashlib.sha256(repr(((BODY, terminal[17]),)).encode()).hexdigest()
    require(case_sha == EXPECTED_CASE_SHA256, case_sha)
    require(LEDGER_BEFORE - 1 == LEDGER_AFTER, "ledger arithmetic")

    semantic_packet = (
        "lrc14-k3-z216-unique-g18-terminal-v1",
        (SOURCE_3139_SHA256, OUTPUT_3139_SHA256, SEMANTIC_3139_SHA256,
         SOURCE_3251_SHA256, OUTPUT_3251_SHA256, SEMANTIC_3251_SHA256, ATLAS_SHA256),
        (LEVEL, EXPECTED_CENSUS, EXPECTED_ROW_SHA256, EXPECTED_STRATA, ATLAS_INDEX,
         BODY, RULER, HIGH, EXPECTED_TOTALS, EXPECTED_FARKAS, EXPECTED_SCREEN_SHA256,
         EXPECTED_ROW_SEMANTIC_SHA256, EXPECTED_ROW_FAILURE_SHA256,
         EXPECTED_MASK_SHA256, EXPECTED_RESIDUAL_SHA256, EXPECTED_QUOTIENT),
        (EXPECTED_TERMINAL_SUMMARY, EXPECTED_TERMINAL_SEMANTIC_SHA256,
         EXPECTED_CLOSURE_SHA256, EXPECTED_CASE_SHA256, EXPECTED_GAP),
        (LEDGER_BEFORE, LEDGER_AFTER, 216, 479, 446, 33),
    )
    semantic = hashlib.sha256(repr(semantic_packet).encode()).hexdigest()
    lines = [
        "LRC14 THM3261 projected-k3 z216 unique gcd18 terminal descent",
        f"engine=THM3139_source:{SOURCE_3139_SHA256};output:{OUTPUT_3139_SHA256};semantic:{SEMANTIC_3139_SHA256}",
        f"dependency=THM3251_source:{SOURCE_3251_SHA256};output:{OUTPUT_3251_SHA256};semantic:{SEMANTIC_3251_SHA256}",
        f"atlas=sha256:{ATLAS_SHA256};rows:6060",
        f"layer=z1:216;rows:480;wall:447;order:33;row_sha256:{EXPECTED_ROW_SHA256};strata:{EXPECTED_STRATA}",
        f"selected=atlas:{ATLAS_INDEX};E:{BODY};L:{RULER};high:{HIGH};gcd:18;first_d:{RULER//18}",
        f"screen=states:{totals[0]};crude:{totals[1]};status:{totals[2]};residual:{totals[3]};direct_farkas:{farkas[0]};legacy_farkas:{farkas[1]};screen_sha256:{EXPECTED_SCREEN_SHA256};row_semantic:{screen[15]};row_failure:{screen[18]}",
        f"residual=masks:30;mask_sha256:{EXPECTED_MASK_SHA256};residual_sha256:{EXPECTED_RESIDUAL_SHA256};quotient:{quotient_packet}",
        f"terminal=rows:1;masks:{terminal[4]};positive_two_high_gap:{int(terminal[19])};closed:{int(terminal[20])};zero_high_hostiles:{terminal[7]};one_high_cases:{terminal[8]};low_label_sets:{terminal[9]};coarse:{terminal[10]};exact:{terminal[11]};maxgap:{terminal[12]};failures:{terminal[14]};minimum_slack:{terminal[16]};gap:{ftext(terminal[5])};semantic_sha256:{terminal_semantic_sha};closure_sha256:{closure_sha};case_sha256:{case_sha}",
        "direction=the_exact_screen_and_complete_one_high_terminal_bank_are_upper_relaxations;their_exclusion_closes_the_unique_gcd18_row",
        "modular_boundary=Div18_is_P2_times_P3_and_every_residual_is_top_q18_but_no_C2_or_C3_action_common_carrier_or_free_factor_is_constructed",
        f"consequence=ledger:{LEDGER_BEFORE}-1={LEDGER_AFTER};projected_k3_cap:z1<=216;remaining_z216:479=446wall+33order",
        "scope=one_projected_k3_necessary_atlas_row_only;no_physical_cover_classification_k_le_1_final_rung_or_LRC14_claim",
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
