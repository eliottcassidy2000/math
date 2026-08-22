#!/usr/bin/env python3
"""Exact referee for the complete projected-k3 z216 order-row screen descent."""

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
SOURCE_3139 = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z225_terminal_z224_screen_double_layer_descent_thm3139.py"
)
OUTPUT_3139 = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z225_terminal_z224_screen_double_layer_descent_thm3139.out"
)
SOURCE_3264 = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z216_low_cost_gcd8_terminal_descent_thm3264.py"
)
OUTPUT_3264 = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z216_low_cost_gcd8_terminal_descent_thm3264.out"
)

SOURCE_3139_SHA256 = "92e1c22088998d37db89020ac1ebbb7d1ee17e3539152ab205f8b0cd92532e36"
OUTPUT_3139_SHA256 = "4f0ae623d0134406f83f466c0a9f1353525a992acb8ee82c2e92f2f0537f32c5"
SEMANTIC_3139_SHA256 = "1eb7e8faefafe41fd6f6cbc108ba09295dbcbb7fe8016d7ec0cbdc258adf358a"
SOURCE_3264_SHA256 = "d9841a86850c6609e62d0522b50ea38722cd5850f7173374a3b143ef577e0e3e"
OUTPUT_3264_SHA256 = "2eebe5f7acb2d1c02fa126ee03166cd274ed209cc52d4b8a729a8fbc5f0a9782"
SEMANTIC_3264_SHA256 = "64604458a454f7da3468b07ec5697a6dd62ac4413625a92dbd8e8ffff15e1a7e"
ATLAS_SHA256 = "cee82237ce1f51729813b9c916edd3353204c18172abe1d71278dee2c5562eda"

LEVEL = 216
EXPECTED_CENSUS = (480, 447, 33)
EXPECTED_ROW_SHA256 = "53db9e1d3df2cf2b0398847682d909da81705e43a53ae2553d102fd152337649"
EXPECTED_STRATA = ((8, 19), (18, 1), (24, 135), (36, 15), (72, 310))
EXPECTED_ORDER_INDICES = (
    0, 1, 2, 4, 12, 13, 16, 20, 22, 29, 48, 51, 70, 73, 74, 85, 86,
    99, 100, 114, 128, 130, 172, 241, 244, 253, 255, 267, 287, 289,
    291, 317, 397,
)
EXPECTED_FAMILIES = (
    (24, 1176, 1, 25_872),
    (24, 1680, 23, 836_640),
    (72, 504, 1, 9_072),
    (72, 1008, 8, 157_248),
)
EXPECTED_BY_GCD = ((24, 24, 862_512), (72, 9, 166_320))
EXPECTED_COST = 1_028_832
EXPECTED_ORDER_PACKET_SHA256 = "ce13ef797ac9962b3a5cee0fa41450ccb6a2545e5df03463c375b620080c6060"
EXPECTED_PRIOR_INDICES = (
    8, 23, 39, 57, 64, 66, 115, 142, 150, 152, 197, 272, 277, 300,
    304, 306, 338,
)
EXPECTED_TOTALS = (481, 463, 18, 0)
EXPECTED_FARKAS = (0, 18)
EXPECTED_SCREEN_SHA256 = "f35dcbe7658ea48df9cfe08c18602ec17242afaeacfcdc88f21d310bfd328ced"
EXPECTED_EMPTY_RESIDUAL_SHA256 = "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d"
LEDGER_BEFORE = 373_266
LEDGER_AFTER = 373_233


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
        text = path.read_text(encoding="utf-8")
        require(f"semantic_sha256={semantic}" in text, (path, "semantic changed"))


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def atlas_rows():
    entry = load("order_atlas_entry", SOURCE_3139)
    source = entry.load("order_atlas_source")
    driver = source.load("order_atlas_driver")
    predecessor = driver.load("order_atlas_predecessor")
    bridge = predecessor.load("order_atlas_bridge")
    base = bridge.load("order_atlas_base")
    require(base.thm.ATLAS_SHA256 == ATLAS_SHA256, base.thm.ATLAS_SHA256)
    pattern = re.compile(
        r"^row=E=([0-9,]+);h=[^;]+;r=([0-9]+);"
        r"L=([0-9]+);high=([0-9]+);z1=([0-9]+);"
    )
    rows = []
    components = []
    for line in base.thm.ATLAS.read_text(encoding="utf-8").splitlines():
        if not line.startswith("row="):
            continue
        match = pattern.match(line)
        require(match is not None, line)
        if int(match.group(5)) != LEVEL:
            continue
        body = tuple(map(int, match.group(1).split(",")))
        components.append(int(match.group(2)))
        ruler = int(match.group(3))
        high = int(match.group(4))
        rows.append((body, ruler, high, LEVEL < high))
    return tuple(rows), tuple(components)


def worker(indexed_task):
    index, task = indexed_task
    entry = load(f"order_screen_worker_{index}", SOURCE_3139)
    returned, full = entry.screen_worker(task)
    require(returned == task, (index, returned, task))
    return index, tuple(full[:19]), int(full[19]), int(full[20])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=6)
    args = parser.parse_args()
    require(args.processes >= 1, args.processes)

    check_dependency(SOURCE_3139, SOURCE_3139_SHA256)
    check_dependency(OUTPUT_3139, OUTPUT_3139_SHA256, SEMANTIC_3139_SHA256)
    check_dependency(SOURCE_3264, SOURCE_3264_SHA256)
    check_dependency(OUTPUT_3264, OUTPUT_3264_SHA256, SEMANTIC_3264_SHA256)

    rows, components = atlas_rows()
    census = (len(rows), sum(row[3] for row in rows), sum(not row[3] for row in rows))
    require(census == EXPECTED_CENSUS, census)
    require(hashlib.sha256(repr(rows).encode()).hexdigest() == EXPECTED_ROW_SHA256,
            "row order")
    strata = Counter(gcd(LEVEL, row[1]) for row in rows)
    require(tuple(sorted(strata.items())) == EXPECTED_STRATA, strata)

    order_indices = tuple(index for index, row in enumerate(rows) if not row[3])
    require(order_indices == EXPECTED_ORDER_INDICES, order_indices)
    order_packet = tuple(
        (index, rows[index][0], gcd(LEVEL, rows[index][1]), rows[index][1],
         components[index], rows[index][1] * components[index])
        for index in order_indices
    )
    require(hashlib.sha256(repr(order_packet).encode()).hexdigest()
            == EXPECTED_ORDER_PACKET_SHA256, "order packet")
    require(sum(item[-1] for item in order_packet) == EXPECTED_COST, "cost")

    families = defaultdict(lambda: [0, 0])
    by_gcd = defaultdict(lambda: [0, 0])
    for _index, _body, divisor_gcd, ruler, _components, cost in order_packet:
        families[(divisor_gcd, ruler)][0] += 1
        families[(divisor_gcd, ruler)][1] += cost
        by_gcd[divisor_gcd][0] += 1
        by_gcd[divisor_gcd][1] += cost
    family_packet = tuple(
        (divisor_gcd, ruler, *families[(divisor_gcd, ruler)])
        for divisor_gcd, ruler in sorted(families)
    )
    gcd_packet = tuple((divisor_gcd, *by_gcd[divisor_gcd]) for divisor_gcd in sorted(by_gcd))
    require(family_packet == EXPECTED_FAMILIES, family_packet)
    require(gcd_packet == EXPECTED_BY_GCD, gcd_packet)

    prior = load("order_prior_thm3264", SOURCE_3264)
    require(tuple(prior.EXPECTED_SELECTED) == EXPECTED_PRIOR_INDICES,
            prior.EXPECTED_SELECTED)
    prior_rederived = tuple(
        index for index, row in enumerate(rows)
        if gcd(LEVEL, row[1]) == 8
        and row[1] * components[index] <= prior.COST_LIMIT
    )
    require(prior_rederived == EXPECTED_PRIOR_INDICES, prior_rederived)
    require(all(rows[index][3] for index in prior_rederived), "prior row not wall")
    require(all(not rows[index][3] for index in order_indices), "order row is wall")
    require(set(prior_rederived).isdisjoint(order_indices), "prior/order overlap")

    tasks = tuple((index, (LEVEL, *rows[index])) for index in order_indices)
    if args.processes == 1:
        screened = tuple(worker(task) for task in tasks)
    else:
        with mp.get_context("spawn").Pool(min(args.processes, len(tasks))) as pool:
            screened = tuple(pool.map(worker, tasks, chunksize=1))
    require(tuple(index for index, _row, _f0, _f1 in screened) == order_indices,
            "screen order")

    canonical = []
    direct_farkas = 0
    legacy_farkas = 0
    for index, row, direct, legacy in screened:
        source_row = rows[index]
        divisor_gcd = gcd(LEVEL, source_row[1])
        require(row[:6] == (
            LEVEL, source_row[0], source_row[1], source_row[2],
            source_row[1] // divisor_gcd, source_row[3],
        ), (index, row[:6]))
        require(row[16] == row[11], (index, "status verification"))
        require(row[12] == 0 and row[13] == (),
                (index, "unexpected residual"))
        canonical.append(row)
        direct_farkas += direct
        legacy_farkas += legacy
    canonical = tuple(canonical)
    totals = tuple(sum(row[position] for row in canonical) for position in (9, 10, 11, 12))
    farkas = (direct_farkas, legacy_farkas)
    require(totals == EXPECTED_TOTALS, totals)
    require(farkas == EXPECTED_FARKAS, farkas)
    screen_sha = hashlib.sha256(repr(canonical).encode()).hexdigest()
    require(screen_sha == EXPECTED_SCREEN_SHA256, screen_sha)
    residual = tuple((row[1], row[13]) for row in canonical if row[12])
    residual_sha = hashlib.sha256(repr(residual).encode()).hexdigest()
    require(residual == () and residual_sha == EXPECTED_EMPTY_RESIDUAL_SHA256,
            (residual, residual_sha))

    require(LEDGER_BEFORE - len(order_indices) == LEDGER_AFTER, "ledger")
    require(462 - len(order_indices) == 429, "layer arithmetic")
    require(EXPECTED_CENSUS[1] - 18 == 429, "wall arithmetic")

    semantic_packet = (
        "lrc14-k3-z216-order-row-screen-descent-v1",
        (SOURCE_3139_SHA256, OUTPUT_3139_SHA256, SEMANTIC_3139_SHA256,
         SOURCE_3264_SHA256, OUTPUT_3264_SHA256, SEMANTIC_3264_SHA256,
         ATLAS_SHA256),
        (LEVEL, EXPECTED_CENSUS, EXPECTED_ROW_SHA256, EXPECTED_STRATA,
         EXPECTED_ORDER_INDICES, EXPECTED_FAMILIES, EXPECTED_BY_GCD,
         EXPECTED_COST, EXPECTED_ORDER_PACKET_SHA256),
        (EXPECTED_PRIOR_INDICES, 0, EXPECTED_TOTALS, EXPECTED_FARKAS,
         EXPECTED_SCREEN_SHA256, EXPECTED_EMPTY_RESIDUAL_SHA256),
        (LEDGER_BEFORE, LEDGER_AFTER, 216, 429, 429, 0),
    )
    semantic = hashlib.sha256(repr(semantic_packet).encode()).hexdigest()
    lines = [
        "LRC14 THM3270 projected-k3 z216 complete order-row screen descent",
        f"engine=THM3139_source:{SOURCE_3139_SHA256};output:{OUTPUT_3139_SHA256};semantic:{SEMANTIC_3139_SHA256}",
        f"dependency=THM3264_source:{SOURCE_3264_SHA256};output:{OUTPUT_3264_SHA256};semantic:{SEMANTIC_3264_SHA256}",
        f"atlas=sha256:{ATLAS_SHA256};rows:6060",
        f"layer=z1:216;rows:480;wall:447;order:33;row_sha256:{EXPECTED_ROW_SHA256};strata:{EXPECTED_STRATA}",
        f"selection=order_rows:{EXPECTED_ORDER_INDICES};families:{EXPECTED_FAMILIES};by_gcd:{EXPECTED_BY_GCD};cost:{EXPECTED_COST};packet_sha256:{EXPECTED_ORDER_PACKET_SHA256}",
        f"disjointness=THM3264_wall_indices:{EXPECTED_PRIOR_INDICES};intersection:0;prior_wall:17;selected_order:33",
        f"screen=rows:33;states:{totals[0]};crude:{totals[1]};status:{totals[2]};residual:{totals[3]};direct_farkas:{farkas[0]};legacy_farkas:{farkas[1]};screen_sha256:{screen_sha}",
        f"residual=bodies:0;masks:0;residual_sha256:{residual_sha}",
        "direction=the_complete_THM3139_screen_is_an_upper_relaxation_and_excludes_every_order_row_without_terminal_analysis",
        f"consequence=ledger:{LEDGER_BEFORE}-33={LEDGER_AFTER};projected_k3_cap:z1<=216;remaining_z216:429=429wall+0order",
        "scope=complete_z216_order_sublayer_only;no_terminal_probe_divisor_action_parity_tower_common_carrier_physical_cover_or_LRC14_claim",
        f"semantic_sha256={semantic}",
        "all_exact_controls=PASS",
    ]
    print("\n".join(lines))


if __name__ == "__main__":
    mp.freeze_support()
    main()
