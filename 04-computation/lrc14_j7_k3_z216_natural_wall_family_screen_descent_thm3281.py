#!/usr/bin/env python3
"""Exact companion for three intrinsic z216 wall-family closures."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
import re
from collections import Counter
from math import gcd
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
SOURCE_3270 = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z216_order_row_screen_descent_thm3270.py"
)
OUTPUT_3270 = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z216_order_row_screen_descent_thm3270.out"
)

SOURCE_3139_SHA256 = "92e1c22088998d37db89020ac1ebbb7d1ee17e3539152ab205f8b0cd92532e36"
OUTPUT_3139_SHA256 = "4f0ae623d0134406f83f466c0a9f1353525a992acb8ee82c2e92f2f0537f32c5"
SEMANTIC_3139_SHA256 = "1eb7e8faefafe41fd6f6cbc108ba09295dbcbb7fe8016d7ec0cbdc258adf358a"
SOURCE_3264_SHA256 = "d9841a86850c6609e62d0522b50ea38722cd5850f7173374a3b143ef577e0e3e"
OUTPUT_3264_SHA256 = "2eebe5f7acb2d1c02fa126ee03166cd274ed209cc52d4b8a729a8fbc5f0a9782"
SEMANTIC_3264_SHA256 = "64604458a454f7da3468b07ec5697a6dd62ac4413625a92dbd8e8ffff15e1a7e"
SOURCE_3270_SHA256 = "dbcc60e3d691483e023486cdb9b935381e5172cea2b828b51ed3da99560fe2ab"
OUTPUT_3270_SHA256 = "0fa283cc68ac831579e17b4e9d817cd5dd13a4c58363d51335e0c1911897d1f0"
SEMANTIC_3270_SHA256 = "a8eb670dbb06bedc9c9180dc23041a77cb5c5f717f1d2ba2b1de52f8322bc602"
ATLAS_SHA256 = "cee82237ce1f51729813b9c916edd3353204c18172abe1d71278dee2c5562eda"

LEVEL = 216
EXPECTED_CENSUS = (480, 447, 33)
EXPECTED_ROW_SHA256 = "53db9e1d3df2cf2b0398847682d909da81705e43a53ae2553d102fd152337649"
EXPECTED_STRATA = ((8, 19), (18, 1), (24, 135), (36, 15), (72, 310))
EXPECTED_FAMILIES = (
    (
        24,
        2352,
        (14, 24, 44, 45, 54, 76, 88, 124, 125, 132, 204, 242, 246,
         257, 284, 285, 294, 346, 375, 423),
        1_086_624,
        "f4101ec3457957d91aef1fefcb2771cb08f400d8b2252e8a680519f192e252e7",
        (710, 309, 401, 0),
        (0, 401),
        "44984dba0ca096cf43eafb6479f1977eb4c7ea9eefbb9c804c1844837008e459",
    ),
    (
        36,
        8820,
        (32, 41, 178, 323),
        829_080,
        "a03cf6fbcf98781424ccc44b0b92c55c9b22625cbe698ea1098997534459e777",
        (179, 145, 34, 0),
        (0, 34),
        "671d9823cd1ef155cd58d30cdc0a56b7ec4023b6b13d1d9258f1083ec71b2c93",
    ),
    (
        72,
        2520,
        (6, 9, 10, 15, 25, 31, 40, 55, 68, 69, 71, 77, 102, 103,
         117, 135, 177, 269, 270, 279, 295, 322, 402),
        1_144_080,
        "3a1c6a5d45f9fb5ef8c763ff0d32540dd4b03496a99c6f17f19793faa7b5e553",
        (410, 400, 10, 0),
        (0, 10),
        "ede315ae68c29463aff5f6e02dbcca37e0283c449b5936510bffa4209fba11a0",
    ),
)
EXPECTED_UNION_INDICES = (
    6, 9, 10, 14, 15, 24, 25, 31, 32, 40, 41, 44, 45, 54, 55, 68,
    69, 71, 76, 77, 88, 102, 103, 117, 124, 125, 132, 135, 177, 178,
    204, 242, 246, 257, 269, 270, 279, 284, 285, 294, 295, 322, 323,
    346, 375, 402, 423,
)
EXPECTED_UNION_COST = 3_059_784
EXPECTED_UNION_PACKET_SHA256 = "967d9ef712a12d081394837460abff934e0cc81f7a356e1f6da6bd677d077560"
EXPECTED_TOTALS = (1299, 854, 445, 0)
EXPECTED_FARKAS = (0, 445)
EXPECTED_SCREEN_SHA256 = "05ee24dbeeaec5e67e4a2ae93482453493f1d6edc226d80187c6b32e4acea8de"
EXPECTED_EMPTY_RESIDUAL_SHA256 = "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d"
EXPECTED_SEMANTIC_SHA256 = "2330ce22ee3b318da26e70ce2c6c01b020c8f2b78bd958d529640a599514cd32"
LEDGER_BEFORE = 373_233
LEDGER_AFTER = 373_186
WALL_BEFORE = 429
WALL_AFTER = 382


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
        require(
            f"semantic_sha256={semantic}" in path.read_text(encoding="utf-8"),
            (path, "semantic changed"),
        )


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def atlas_rows():
    entry = load("natural_wall_atlas_entry", SOURCE_3139)
    source = entry.load("natural_wall_atlas_source")
    driver = source.load("natural_wall_atlas_driver")
    predecessor = driver.load("natural_wall_atlas_predecessor")
    bridge = predecessor.load("natural_wall_atlas_bridge")
    base = bridge.load("natural_wall_atlas_base")
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
        rows.append(
            (
                tuple(map(int, match.group(1).split(","))),
                int(match.group(3)),
                int(match.group(4)),
                LEVEL < int(match.group(4)),
            )
        )
        components.append(int(match.group(2)))
    return tuple(rows), tuple(components)


def screen_worker(indexed_task):
    index, task = indexed_task
    entry = load(f"natural_wall_worker_{index}", SOURCE_3139)
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
    check_dependency(SOURCE_3270, SOURCE_3270_SHA256)
    check_dependency(OUTPUT_3270, OUTPUT_3270_SHA256, SEMANTIC_3270_SHA256)

    rows, components = atlas_rows()
    census = (len(rows), sum(row[3] for row in rows), sum(not row[3] for row in rows))
    require(census == EXPECTED_CENSUS, census)
    require(hashlib.sha256(repr(rows).encode()).hexdigest() == EXPECTED_ROW_SHA256,
            "row order")
    strata = Counter(gcd(LEVEL, row[1]) for row in rows)
    require(tuple(sorted(strata.items())) == EXPECTED_STRATA, strata)

    family_keys = {(family[0], family[1]) for family in EXPECTED_FAMILIES}
    union_indices = tuple(
        index for index, row in enumerate(rows)
        if (gcd(LEVEL, row[1]), row[1]) in family_keys
    )
    require(union_indices == EXPECTED_UNION_INDICES, union_indices)
    require(all(rows[index][3] for index in union_indices), "selected non-wall row")

    union_packet = tuple(
        (
            index,
            rows[index][0],
            gcd(LEVEL, rows[index][1]),
            rows[index][1],
            components[index],
            rows[index][1] * components[index],
            rows[index][3],
        )
        for index in union_indices
    )
    require(sum(item[5] for item in union_packet) == EXPECTED_UNION_COST,
            "union cost")
    require(
        hashlib.sha256(repr(union_packet).encode()).hexdigest()
        == EXPECTED_UNION_PACKET_SHA256,
        "union packet",
    )

    prior = load("natural_wall_prior_3264", SOURCE_3264)
    order = load("natural_wall_prior_3270", SOURCE_3270)
    prior_indices = tuple(prior.EXPECTED_SELECTED)
    order_indices = tuple(order.EXPECTED_ORDER_INDICES)
    require(set(union_indices).isdisjoint(prior_indices), "THM3264 overlap")
    require(set(union_indices).isdisjoint(order_indices), "THM3270 overlap")

    tasks = tuple((index, (LEVEL, *rows[index])) for index in union_indices)
    if args.processes == 1:
        screened = tuple(screen_worker(task) for task in tasks)
    else:
        with mp.get_context("spawn").Pool(min(args.processes, len(tasks))) as pool:
            screened = tuple(pool.map(screen_worker, tasks, chunksize=1))
    require(tuple(index for index, _row, _direct, _legacy in screened)
            == union_indices, "screen order")

    canonical = []
    direct = 0
    legacy = 0
    for index, row, direct_count, legacy_count in screened:
        source_row = rows[index]
        divisor_gcd = gcd(LEVEL, source_row[1])
        require(
            row[:6]
            == (
                LEVEL,
                source_row[0],
                source_row[1],
                source_row[2],
                source_row[1] // divisor_gcd,
                source_row[3],
            ),
            (index, row[:6]),
        )
        require(row[16] == row[11], (index, "status verification"))
        require(direct_count + legacy_count == row[11],
                (index, "Farkas count"))
        require(row[12] == 0 and row[13] == (), (index, "residual"))
        canonical.append(row)
        direct += direct_count
        legacy += legacy_count
    canonical = tuple(canonical)
    totals = tuple(sum(row[position] for row in canonical)
                   for position in (9, 10, 11, 12))
    farkas = (direct, legacy)
    screen_sha = hashlib.sha256(repr(canonical).encode()).hexdigest()
    residual = tuple((row[1], row[13]) for row in canonical if row[12])
    residual_sha = hashlib.sha256(repr(residual).encode()).hexdigest()
    require(totals == EXPECTED_TOTALS, totals)
    require(farkas == EXPECTED_FARKAS, farkas)
    require(screen_sha == EXPECTED_SCREEN_SHA256, screen_sha)
    require(residual_sha == EXPECTED_EMPTY_RESIDUAL_SHA256, residual_sha)

    family_lines = []
    family_semantics = []
    by_index = {index: row for index, row, _direct, _legacy in screened}
    count_by_index = {
        index: (direct_count, legacy_count)
        for index, _row, direct_count, legacy_count in screened
    }
    for (
        divisor_gcd,
        ruler,
        indices,
        expected_cost,
        packet_sha,
        expected_totals,
        expected_farkas,
        expected_screen_sha,
    ) in EXPECTED_FAMILIES:
        derived_indices = tuple(
            index for index in union_indices
            if (gcd(LEVEL, rows[index][1]), rows[index][1])
            == (divisor_gcd, ruler)
        )
        require(derived_indices == indices, (divisor_gcd, ruler, derived_indices))
        packet = tuple(
            (
                index,
                rows[index][0],
                divisor_gcd,
                ruler,
                components[index],
                ruler * components[index],
                rows[index][3],
            )
            for index in indices
        )
        cost = sum(item[5] for item in packet)
        require(cost == expected_cost, (divisor_gcd, ruler, cost))
        require(hashlib.sha256(repr(packet).encode()).hexdigest() == packet_sha,
                (divisor_gcd, ruler, "packet"))
        family_canonical = tuple(by_index[index] for index in indices)
        family_totals = tuple(sum(row[position] for row in family_canonical)
                              for position in (9, 10, 11, 12))
        family_farkas = (
            sum(count_by_index[index][0] for index in indices),
            sum(count_by_index[index][1] for index in indices),
        )
        family_screen_sha = hashlib.sha256(
            repr(family_canonical).encode()
        ).hexdigest()
        require(family_totals == expected_totals,
                (divisor_gcd, ruler, family_totals))
        require(family_farkas == expected_farkas,
                (divisor_gcd, ruler, family_farkas))
        require(family_screen_sha == expected_screen_sha,
                (divisor_gcd, ruler, family_screen_sha))
        family_lines.append(
            "family="
            f"gcd{divisor_gcd}/L{ruler};rows:{len(indices)};cost:{cost};"
            f"indices:{','.join(map(str, indices))};states:{family_totals[0]};"
            f"crude:{family_totals[1]};status:{family_totals[2]};"
            f"residual:{family_totals[3]};direct_farkas:{family_farkas[0]};"
            f"legacy_farkas:{family_farkas[1]};packet_sha256:{packet_sha};"
            f"screen_sha256:{family_screen_sha}"
        )
        family_semantics.append(
            (
                divisor_gcd,
                ruler,
                packet,
                family_canonical,
                family_totals,
                family_farkas,
                family_screen_sha,
            )
        )

    require(LEDGER_BEFORE - len(union_indices) == LEDGER_AFTER, "ledger")
    require(WALL_BEFORE - len(union_indices) == WALL_AFTER, "wall count")
    semantic_packet = (
        SOURCE_3139_SHA256,
        OUTPUT_3139_SHA256,
        SEMANTIC_3139_SHA256,
        SOURCE_3264_SHA256,
        OUTPUT_3264_SHA256,
        SEMANTIC_3264_SHA256,
        SOURCE_3270_SHA256,
        OUTPUT_3270_SHA256,
        SEMANTIC_3270_SHA256,
        ATLAS_SHA256,
        EXPECTED_ROW_SHA256,
        census,
        tuple(family_semantics),
        union_packet,
        totals,
        farkas,
        screen_sha,
        residual_sha,
        prior_indices,
        order_indices,
        LEDGER_BEFORE,
        LEDGER_AFTER,
        WALL_BEFORE,
        WALL_AFTER,
    )
    semantic = hashlib.sha256(repr(semantic_packet).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, semantic)

    lines = [
        "LRC14 THM3281 projected-k3 z216 three natural wall-family screen descent",
        (
            "dependency="
            f"THM3139:{SOURCE_3139_SHA256}/{OUTPUT_3139_SHA256}/{SEMANTIC_3139_SHA256};"
            f"THM3264:{SOURCE_3264_SHA256}/{OUTPUT_3264_SHA256}/{SEMANTIC_3264_SHA256};"
            f"THM3270:{SOURCE_3270_SHA256}/{OUTPUT_3270_SHA256}/{SEMANTIC_3270_SHA256}"
        ),
        (
            f"universe=z1:{LEVEL};rows:{census[0]};wall:{census[1]};"
            f"order:{census[2]};row_sha256:{EXPECTED_ROW_SHA256};"
            f"strata:{EXPECTED_STRATA}"
        ),
        *family_lines,
        (
            f"union=rows:{len(union_indices)};cost:{EXPECTED_UNION_COST};"
            f"indices:{','.join(map(str, union_indices))};"
            f"packet_sha256:{EXPECTED_UNION_PACKET_SHA256};states:{totals[0]};"
            f"crude:{totals[1]};status:{totals[2]};residual:{totals[3]};"
            f"direct_farkas:{farkas[0]};legacy_farkas:{farkas[1]};"
            f"screen_sha256:{screen_sha};residual_sha256:{residual_sha}"
        ),
        "disjoint=THM3264_indices:PASS;THM3270_indices:PASS;all_selected_wall:PASS",
        (
            f"consequence=ledger:{LEDGER_BEFORE}-{len(union_indices)}={LEDGER_AFTER};"
            f"z216_wall:{WALL_BEFORE}-{len(union_indices)}={WALL_AFTER};order:0"
        ),
        (
            "direction=the_THM3139_ray_quotient_is_an_upper_relaxation_of_every_"
            "actual_projected_packet;crude_upper_exclusions_and_legacy_full-table_"
            "exact_Farkas_contradictions_close_the_enlarged_state_set;no_converse_"
            "from_the_quotient_is_used"
        ),
        (
            "evidence_boundary=all_445_status_exclusions_are_verified_inside_the_"
            "inherited_worker;the_screen_digest_uses_only_basis-invariant_row[:19]_"
            "and_branch_counts_not_solver-selected_dual_representatives"
        ),
        (
            "scope=projected_k3_z216_necessary_wall_rows_only;no_terminal_probe;"
            "no_divisor_action;no_common_carrier_or_phase;no_physical_cover_or_LRC14_claim"
        ),
        f"semantic_sha256={semantic}",
        "all_exact_controls=PASS",
    ]
    payload = "\n".join(lines) + "\n"
    print(payload, end="")


if __name__ == "__main__":
    mp.freeze_support()
    main()
