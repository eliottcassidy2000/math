#!/usr/bin/env python3
"""Exact scratch probe for the two costly projected-k3 z216 gcd-8 wall rows.

This intentionally imports the promoted THM-3139 screening/terminal engine
without changing it.  It freezes the two rows omitted only by THM-3264's
work-invoice cutoff, runs the full upper-relaxation screen, then applies the
same complete terminal probe used by that theorem.  The output is concise and
deterministic so normal and optimized runs can be compared byte-for-byte.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
from collections import Counter
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ENTRY_PATH = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z225_terminal_z224_screen_double_layer_descent_thm3139.py"
)
ENTRY_OUTPUT_PATH = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z225_terminal_z224_screen_double_layer_descent_thm3139.out"
)
THM3264_PATH = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z216_low_cost_gcd8_terminal_descent_thm3264.py"
)
THM3264_OUTPUT_PATH = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z216_low_cost_gcd8_terminal_descent_thm3264.out"
)
ENTRY_SHA256 = "92e1c22088998d37db89020ac1ebbb7d1ee17e3539152ab205f8b0cd92532e36"
ENTRY_OUTPUT_SHA256 = "4f0ae623d0134406f83f466c0a9f1353525a992acb8ee82c2e92f2f0537f32c5"
ENTRY_SEMANTIC_SHA256 = "1eb7e8faefafe41fd6f6cbc108ba09295dbcbb7fe8016d7ec0cbdc258adf358a"
THM3264_SHA256 = "d9841a86850c6609e62d0522b50ea38722cd5850f7173374a3b143ef577e0e3e"
THM3264_OUTPUT_SHA256 = "2eebe5f7acb2d1c02fa126ee03166cd274ed209cc52d4b8a729a8fbc5f0a9782"
THM3264_SEMANTIC_SHA256 = "64604458a454f7da3468b07ec5697a6dd62ac4413625a92dbd8e8ffff15e1a7e"
LEVEL = 216
TARGETS = (238, 370)
EXPECTED_COMPONENT_COUNT = 34
EXPECTED_ROWS = {
    238: ((1, 8, 10, 11, 13, 14), 560560, 55207, True),
    370: ((2, 8, 10, 11, 13, 14), 560560, 55207, True),
}


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def lf_sha(path: Path) -> str:
    payload = path.read_bytes().replace(b"\r\n", b"\n")
    require(b"\r" not in payload, (path, "bare CR"))
    return hashlib.sha256(payload).hexdigest()


def check_dependency(path: Path, digest: str, semantic: str | None = None) -> None:
    require(lf_sha(path) == digest, (path, "hash changed"))
    if semantic is not None:
        require(
            f"semantic_sha256={semantic}" in path.read_text(encoding="utf-8"),
            (path, "semantic changed"),
        )


def engine(prefix: str):
    entry = load(f"{prefix}_entry", ENTRY_PATH)
    source = entry.load(f"{prefix}_source")
    driver = source.load(f"{prefix}_driver")
    predecessor = driver.load(f"{prefix}_predecessor")
    bridge = predecessor.load(f"{prefix}_bridge")
    base = bridge.load(f"{prefix}_base")
    return entry, base


def screen_worker(item):
    index, task = item
    entry = load(f"costly_g8_screen_entry_{index}", ENTRY_PATH)
    returned, full = entry.screen_worker(task)
    require(returned == task, (index, returned, task))
    return index, tuple(full)


def terminal_worker(item):
    index, body, bank = item
    _entry, base = engine(f"costly_g8_terminal_{index}")
    return index, tuple(base.thm.terminal_probe((LEVEL, body, bank)))


def frac_text(value) -> str:
    return f"{value.numerator}/{value.denominator}"


def terminal_semantic(row):
    # The inherited engine documents row[6] as a chosen duplicate-gap witness;
    # omit it so the digest is independent of witness choice.
    return (*row[:6], *row[7:])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=2)
    args = parser.parse_args()
    require(args.processes >= 1, args.processes)

    check_dependency(ENTRY_PATH, ENTRY_SHA256)
    check_dependency(ENTRY_OUTPUT_PATH, ENTRY_OUTPUT_SHA256, ENTRY_SEMANTIC_SHA256)
    check_dependency(THM3264_PATH, THM3264_SHA256)
    check_dependency(
        THM3264_OUTPUT_PATH, THM3264_OUTPUT_SHA256, THM3264_SEMANTIC_SHA256
    )

    entry, base = engine("costly_g8_main")
    rows = entry.atlas_rows(base, LEVEL)
    require(len(rows) == 480, len(rows))
    require(
        tuple(i for i, row in enumerate(rows) if gcd(LEVEL, row[1]) == 8)
        == (8, 23, 39, 57, 64, 66, 115, 142, 150, 152, 197, 238,
            272, 277, 300, 304, 306, 338, 370),
        "gcd-8 stratum changed",
    )
    for index in TARGETS:
        require(rows[index] == EXPECTED_ROWS[index], (index, rows[index]))
    require(560560 * EXPECTED_COMPONENT_COUNT == 19059040, "invoice changed")

    tasks = tuple((index, (LEVEL, *rows[index])) for index in TARGETS)
    if args.processes == 1:
        screened_pairs = tuple(screen_worker(task) for task in tasks)
    else:
        with mp.get_context("spawn").Pool(min(args.processes, len(tasks))) as pool:
            screened_pairs = tuple(pool.map(screen_worker, tasks))
    screened_map = dict(screened_pairs)
    screened = tuple(screened_map[index] for index in TARGETS)

    for index, row in zip(TARGETS, screened):
        expected = EXPECTED_ROWS[index]
        require(row[:6] == (LEVEL, expected[0], expected[1], expected[2], 70070, True),
                (index, row[:6]))
        require(row[16] == row[11], (index, "status audit"))
        require(row[12] == len(row[13]), (index, "residual count"))
        for mask in row[13]:
            divisor = lcm(*mask)
            require(divisor % row[4] == 0 and row[2] % divisor == 0,
                    (index, mask, divisor))

    terminal_tasks = tuple(
        (index, row[1], row[13]) for index, row in zip(TARGETS, screened)
    )
    if args.processes == 1:
        terminal_pairs = tuple(terminal_worker(task) for task in terminal_tasks)
    else:
        with mp.get_context("spawn").Pool(min(args.processes, len(tasks))) as pool:
            terminal_pairs = tuple(pool.map(terminal_worker, terminal_tasks))
    terminal_map = dict(terminal_pairs)
    terminals = tuple(terminal_map[index] for index in TARGETS)

    print("LRC14 costly gcd8 z216 exact scratch probe")
    print(f"entry={ENTRY_PATH.relative_to(ROOT).as_posix()}")
    print(
        "dependencies="
        f"THM3139:{ENTRY_SHA256}/{ENTRY_OUTPUT_SHA256}/{ENTRY_SEMANTIC_SHA256};"
        f"THM3264:{THM3264_SHA256}/{THM3264_OUTPUT_SHA256}/{THM3264_SEMANTIC_SHA256}"
    )
    print(f"universe=level:{LEVEL};atlas_rows:{len(rows)};targets:{TARGETS}")
    for index, screen, terminal in zip(TARGETS, screened, terminals):
        quotient = Counter(lcm(*mask) // screen[4] for mask in screen[13])
        mask_sha = hashlib.sha256(repr(screen[13]).encode()).hexdigest()
        canonical_sha = hashlib.sha256(repr(screen[:19]).encode()).hexdigest()
        semantic = terminal_semantic(terminal)
        terminal_sha = hashlib.sha256(repr(semantic).encode()).hexdigest()
        print(
            "row="
            f"{index};body:{screen[1]};L:{screen[2]};"
            f"components:{EXPECTED_COMPONENT_COUNT};high_floor:{screen[3]};"
            f"first_d:{screen[4]};wall:{screen[5]}"
        )
        print(
            "screen="
            f"states:{screen[9]};crude:{screen[10]};status:{screen[11]};"
            f"residual:{screen[12]};farkas_direct:{screen[19]};"
            f"farkas_inherited:{screen[20]};quotient:{tuple(sorted(quotient.items()))};"
            f"mask_sha256:{mask_sha};canonical_sha256:{canonical_sha}"
        )
        print(
            "terminal="
            f"masks:{terminal[4]};two_high_gap:{frac_text(terminal[5])};"
            f"zero_high_hostiles:{terminal[7]};one_high_cases:{terminal[8]};"
            f"low_label_sets:{terminal[9]};coarse:{terminal[10]};"
            f"exact:{terminal[11]};maxgap:{terminal[12]};failures:{terminal[14]};"
            f"minimum_slack:{terminal[16]};positive_two_high_gap:{bool(terminal[19])};"
            f"closed:{bool(terminal[20])};semantic_sha256:{terminal_sha}"
        )

    mask_intersection = set(screened[0][13]).intersection(screened[1][13])
    print(
        "comparison="
        f"same_screen_counts:{tuple(screened[0][9:13]) == tuple(screened[1][9:13])};"
        f"same_mask_bank:{screened[0][13] == screened[1][13]};"
        f"mask_intersection:{len(mask_intersection)};"
        f"same_terminal_semantic:{terminal_semantic(terminals[0]) == terminal_semantic(terminals[1])}"
    )
    print(
        "direction=the_screen_and_terminal_banks_are_upper_relaxations;"
        "closed:true_excludes_the_corresponding_projected_wall_row"
    )
    print(
        "hostiles=zero_high_relaxations_retained;duplicate_two_high_gap_checked;"
        "the_two_leading_labels_compared_without_assuming_symmetry"
    )
    print(
        "scope=the_two_THM3264_cost_cutoff_rows_only;no_claim_for_gcd24_36_72,"
        "physical_cover,k_le_1,rung,or_LRC14"
    )


if __name__ == "__main__":
    main()
