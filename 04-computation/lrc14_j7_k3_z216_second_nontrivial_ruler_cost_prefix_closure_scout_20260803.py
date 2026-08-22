#!/usr/bin/env python3
"""Exact continuation through the second non-singleton z216 ruler family.

Starting after the first cost-prefix closure, this reconstructs and ranks the
375 remaining projected-k3 wall rows, screens the prefix through the next
non-singleton family, and classifies every exact-status obstruction without
persisting solver-selected Farkas representatives.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import itertools
import multiprocessing as mp
from collections import Counter, defaultdict
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FIRST_SOURCE = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z216_first_nontrivial_ruler_cost_prefix_closure_scout_20260803.py"
)
FIRST_OUTPUT = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z216_first_nontrivial_ruler_cost_prefix_closure_scout_20260803.out"
)
NATURAL_SOURCE = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z216_natural_wall_family_screen_descent_thm3281.py"
)
ORDER_SOURCE = ROOT / (
    "04-computation/lrc14_j7_k3_z216_order_row_screen_descent_thm3270.py"
)

FIRST_SOURCE_SHA256 = "e3718e71d2f46ea691e269aaf7e05c292b930c4d6a9b6d6d89a835bc38df75fe"
FIRST_OUTPUT_SHA256 = "6dba76e65300371399376aede7e79b9a6fa865b6ebb58475a14ef8e262328be1"
FIRST_SEMANTIC_SHA256 = "44f7e8cf786ed428297f3f389214cab7c0b28158a856be63806a3394e594158f"
NATURAL_SOURCE_SHA256 = "430dee7ba03e0d5c9ae0df72ac512500de4f7056cb4663d1c8468bfb93a49bfe"

LEVEL = 216
EXPECTED_CENSUS = (480, 447, 33, 375, 36)
EXPECTED_ROW_SHA256 = "53db9e1d3df2cf2b0398847682d909da81705e43a53ae2553d102fd152337649"
EXPECTED_LIVE_SHA256 = "d080f2f34cff4d1115ff492beffdc941ce9d34b3f98431a8f4356d48b5e310b2"
EXPECTED_RANK_SHA256 = "139717dbfdfeceaf25ff7b71465815f4eac19145d03cd8c97b36f6df9c246bc5"
EXPECTED_RANK_HEAD = (
    (
        609_840,
        1,
        72,
        27_720,
        (134,),
        ((134, (1, 4, 6, 9, 10, 11), 22, 609_840),),
    ),
    (
        655_200,
        1,
        72,
        32_760,
        (121,),
        ((121, (1, 4, 5, 9, 12, 13), 20, 655_200),),
    ),
    (
        1_361_808,
        16,
        72,
        3_528,
        (17, 27, 46, 56, 79, 89, 126, 138, 205, 248, 258, 286, 298,
         347, 376, 425),
        (
            (17, (1, 2, 4, 6, 9, 14), 22, 77_616),
            (27, (1, 2, 4, 9, 12, 14), 26, 91_728),
            (46, (1, 2, 6, 7, 9, 12), 20, 70_560),
            (56, (1, 2, 6, 9, 12, 14), 26, 91_728),
            (79, (1, 3, 4, 9, 12, 14), 24, 84_672),
            (89, (1, 3, 6, 9, 12, 14), 24, 84_672),
            (126, (1, 4, 6, 7, 9, 12), 20, 70_560),
            (138, (1, 4, 6, 9, 12, 14), 26, 91_728),
            (205, (1, 6, 7, 9, 12, 14), 22, 77_616),
            (248, (2, 3, 4, 9, 12, 14), 26, 91_728),
            (258, (2, 3, 6, 9, 12, 14), 26, 91_728),
            (286, (2, 4, 6, 7, 9, 14), 22, 77_616),
            (298, (2, 4, 6, 9, 12, 14), 28, 98_784),
            (347, (2, 6, 7, 9, 12, 14), 24, 84_672),
            (376, (3, 4, 6, 9, 12, 14), 26, 91_728),
            (425, (4, 6, 7, 9, 12, 14), 24, 84_672),
        ),
    ),
)
EXPECTED_NEXT_FAMILIES = (
    (1_729_728, 1, 72, 72_072, (157,)),
    (2_121_504, 3, 24, 25_872, (53, 293, 357)),
)
EXPECTED_SELECTED_INDICES = (
    134, 121, 17, 27, 46, 56, 79, 89, 126, 138, 205, 248, 258, 286,
    298, 347, 376, 425,
)
EXPECTED_SELECTED_COST = 2_626_848
EXPECTED_SELECTED_PACKET_SHA256 = (
    "ca9f2d6ed1b342e04ebfaf0230764b471143b64bb94b185b123068c121ba4559"
)
EXPECTED_SELECTED_TOTALS = (750, 479, 271, 0)
EXPECTED_SELECTED_FARKAS = (0, 271)
EXPECTED_SELECTED_SCREEN_SHA256 = (
    "1f626e370c995d088b279173eaefa04ea1e7bf4c893dfff135228a636d47e86d"
)
EXPECTED_EMPTY_RESIDUAL_SHA256 = (
    "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d"
)
EXPECTED_TAXONOMY_TOTALS = (248, 1, 11, 11, 264, 1, 11, 271)
EXPECTED_TAXONOMY_BY_ROW = (
    (134, (0, 0, 0, 1)),
    (121, (15, 0, 0, 1)),
    (17, (12, 0, 0, 1)),
    (27, (16, 0, 0, 2)),
    (46, (0, 0, 0, 0)),
    (56, (23, 0, 0, 1)),
    (79, (0, 0, 0, 0)),
    (89, (0, 0, 0, 0)),
    (126, (0, 0, 0, 0)),
    (138, (89, 1, 5, 3)),
    (205, (3, 0, 0, 0)),
    (248, (0, 0, 0, 0)),
    (258, (0, 0, 0, 0)),
    (286, (2, 0, 0, 0)),
    (298, (81, 0, 6, 2)),
    (347, (3, 0, 0, 0)),
    (376, (0, 0, 0, 0)),
    (425, (4, 0, 0, 0)),
)
EXPECTED_TAXONOMY_SHA256 = (
    "4cccaba60ecab8fae74a13c74b515d4c3cb92e4692d7f741de8622ac10b1834b"
)
EXPECTED_WEIGHTED_CORE_SHA256 = (
    "6e75e3c6cd9809ccc4b56c4829bee183c92f4a3328d1bcb1ca2a6a1f3d642ebc"
)
LEDGER_BEFORE = 373_179
LEDGER_AFTER = 373_161
WALL_BEFORE = 375
WALL_AFTER = 357
EXPECTED_SEMANTIC_SHA256 = (
    "a73b57558c5fd74736fa9551c72f94d20b2c053ef7b8f36c27222ae993f8dcd8"
)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


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


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def status_engine(natural, suffix: str):
    entry = load(f"second_ruler_status_entry_{suffix}", natural.SOURCE_3139)
    source = entry.load(f"second_ruler_status_source_{suffix}")
    driver = source.load(f"second_ruler_status_driver_{suffix}")
    predecessor = driver.load(f"second_ruler_status_predecessor_{suffix}")
    bridge = predecessor.load(f"second_ruler_status_bridge_{suffix}")
    base = bridge.load(f"second_ruler_status_base_{suffix}")
    return base.thm.screen_engine.eng


def screen_worker(item):
    index, task = item
    natural = load(f"second_ruler_natural_screen_{index}", NATURAL_SOURCE)
    return natural.screen_worker((index, task))


def coordinate_union_certificates(capacities, marginals, histogram):
    certificates = []
    for threshold, _count in histogram:
        if threshold <= 0:
            continue
        demand = sum(
            count for load_value, count in histogram if load_value >= threshold
        )
        for coordinate_mask in range(1, 16):
            if not all(
                (capacities[pattern] >= threshold)
                == bool(pattern & coordinate_mask)
                for pattern in range(16)
            ):
                continue
            coordinates = tuple(
                coordinate for coordinate in range(4)
                if coordinate_mask & (1 << coordinate)
            )
            supply = sum(marginals[coordinate] for coordinate in coordinates)
            if demand > supply:
                certificates.append(
                    (threshold, coordinate_mask, coordinates, demand, supply,
                     demand - supply)
                )
    return tuple(certificates)


def two_fan_certificates(capacities, marginals, histogram):
    """Find B union (A intersect (C union D)) tail obstructions."""
    certificates = []
    for threshold, _count in histogram:
        if threshold <= 0:
            continue
        demand = sum(
            count for load_value, count in histogram if load_value >= threshold
        )
        for a, b in itertools.permutations(range(4), 2):
            c, d = tuple(sorted(set(range(4)).difference((a, b))))
            if not all(
                (capacities[pattern] >= threshold)
                == (
                    bool(pattern & (1 << b))
                    or (
                        bool(pattern & (1 << a))
                        and bool(pattern & ((1 << c) | (1 << d)))
                    )
                )
                for pattern in range(16)
            ):
                continue
            supply = marginals[b] + min(
                marginals[a], marginals[c] + marginals[d]
            )
            if demand > supply:
                certificates.append(
                    (threshold, (a, b, c, d), demand, supply, demand - supply)
                )
    return tuple(certificates)


def zero_reduced_union_certificates(capacities, marginals, histogram):
    """Find coordinate unions after deleting patterns forced off by zero marginals."""
    forbidden_mask = sum(
        1 << coordinate for coordinate, marginal in enumerate(marginals)
        if marginal == 0
    )
    if not forbidden_mask:
        return ()
    certificates = []
    for threshold, _count in histogram:
        if threshold <= 0:
            continue
        demand = sum(
            count for load_value, count in histogram if load_value >= threshold
        )
        for coordinate_mask in range(1, 16):
            if coordinate_mask & forbidden_mask:
                continue
            if not all(
                (capacities[pattern] >= threshold)
                == bool(pattern & coordinate_mask)
                for pattern in range(16)
                if not (pattern & forbidden_mask)
            ):
                continue
            coordinates = tuple(
                coordinate for coordinate in range(4)
                if coordinate_mask & (1 << coordinate)
            )
            supply = sum(marginals[coordinate] for coordinate in coordinates)
            if demand > supply:
                certificates.append(
                    (
                        threshold, forbidden_mask, coordinate_mask, coordinates,
                        demand, supply, demand - supply,
                    )
                )
    return tuple(certificates)


def taxonomy_worker(item):
    index, task, expected_counts = item
    natural = load(f"second_ruler_natural_taxonomy_{index}", NATURAL_SOURCE)
    eng = status_engine(natural, str(index))
    eng.FIRST = LEVEL
    eng.ray.FIRST = LEVEL
    _first, body, ruler, high, wall = task
    stream = eng.ray.Stream(body)
    require((stream.L, stream.high_floor) == (ruler, high), (index, "atlas"))
    _trials, states, _checks, _signs = eng.ray.ray_quotient_states(stream)
    crude, status, residual = eng.exact_common_status_screen(stream, states)
    counts = (len(states), len(crude), len(status), len(residual))
    require(counts == expected_counts, (index, counts, expected_counts))
    require(wall and not residual, (index, "unexpected residual"))

    records = []
    weighted_core = []
    branch_counts = Counter()
    union_certificate_count = 0
    zero_reduced_certificate_count = 0
    fan_certificate_count = 0
    for divisors, witness in sorted(status.items()):
        q, cofactor, marginals, capacity_set, histogram, certificate = witness
        divisor_lcm = lcm(*divisors)
        require(divisor_lcm == q * cofactor, (index, divisors, "cofactor"))
        rebuilt_marginals, capacities = eng.ray.local.hunter_status_data(
            divisor_lcm, divisors, q
        )
        require(rebuilt_marginals == marginals, (index, divisors, "marginals"))
        require(tuple(sorted(set(capacities))) == capacity_set,
                (index, divisors, "capacities"))
        eng.independent_farkas_check(
            q, marginals, capacities, histogram, certificate
        )
        unions = coordinate_union_certificates(capacities, marginals, histogram)
        zero_reduced = () if unions else zero_reduced_union_certificates(
            capacities, marginals, histogram
        )
        fans = () if (unions or zero_reduced) else two_fan_certificates(
            capacities, marginals, histogram
        )
        if unions:
            branch = "coordinate_union"
            evidence = unions
            union_certificate_count += len(unions)
        elif zero_reduced:
            branch = "zero_reduced_union"
            evidence = zero_reduced
            zero_reduced_certificate_count += len(zero_reduced)
        elif fans:
            branch = "two_fan"
            evidence = fans
            fan_certificate_count += len(fans)
        else:
            branch = "weighted_core"
            evidence = ()
            weighted_core.append(
                (
                    index, divisors, q, cofactor, marginals, capacity_set,
                    histogram, tuple(capacities),
                )
            )
        branch_counts[branch] += 1
        records.append(
            (
                index, divisors, q, cofactor, marginals, capacity_set,
                histogram, branch, evidence,
            )
        )
    return (
        index,
        counts,
        (
            branch_counts["coordinate_union"],
            branch_counts["zero_reduced_union"],
            branch_counts["two_fan"],
            branch_counts["weighted_core"],
        ),
        union_certificate_count,
        zero_reduced_certificate_count,
        fan_certificate_count,
        tuple(records),
        tuple(weighted_core),
    )


def totals(rows):
    return tuple(sum(row[position] for row in rows) for position in (9, 10, 11, 12))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=3)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    require(args.processes >= 1, args.processes)

    check_dependency(FIRST_SOURCE, FIRST_SOURCE_SHA256)
    check_dependency(FIRST_OUTPUT, FIRST_OUTPUT_SHA256, FIRST_SEMANTIC_SHA256)
    check_dependency(NATURAL_SOURCE, NATURAL_SOURCE_SHA256)
    first = load("second_ruler_first_main", FIRST_SOURCE)
    natural = load("second_ruler_natural_main", NATURAL_SOURCE)
    order = load("second_ruler_order_main", ORDER_SOURCE)
    rows, components = natural.atlas_rows()
    require(hashlib.sha256(repr(rows).encode()).hexdigest() == EXPECTED_ROW_SHA256,
            "row order")

    prior_closed = (
        set(natural.EXPECTED_UNION_INDICES)
        | set(order.EXPECTED_ORDER_INDICES)
        | set(first.EXPECTED_SELECTED_INDICES)
        | {
            index for index, row in enumerate(rows)
            if gcd(LEVEL, row[1]) in (8, 18)
        }
    )
    require(len(prior_closed) == 105, len(prior_closed))
    live = tuple(
        index for index, row in enumerate(rows)
        if row[3] and index not in prior_closed
    )
    families = defaultdict(list)
    for index in live:
        families[(gcd(LEVEL, rows[index][1]), rows[index][1])].append(index)
    census = (
        len(rows), sum(row[3] for row in rows), sum(not row[3] for row in rows),
        len(live), len(families),
    )
    require(census == EXPECTED_CENSUS, census)
    live_packet = tuple((index, rows[index], components[index]) for index in live)
    require(
        hashlib.sha256(repr(live_packet).encode()).hexdigest()
        == EXPECTED_LIVE_SHA256,
        "live packet",
    )

    ranked = []
    for (divisor_gcd, ruler), raw_indices in families.items():
        indices = tuple(raw_indices)
        packet = tuple(
            (index, rows[index][0], components[index], ruler * components[index])
            for index in indices
        )
        ranked.append(
            (sum(item[3] for item in packet), len(indices), divisor_gcd,
             ruler, indices, packet)
        )
    ranked = tuple(sorted(ranked))
    require(ranked[:3] == EXPECTED_RANK_HEAD, ranked[:3])
    require(
        tuple((row[0], row[1], row[2], row[3], row[4]) for row in ranked[3:5])
        == EXPECTED_NEXT_FAMILIES,
        ranked[3:5],
    )
    require(
        hashlib.sha256(repr(ranked).encode()).hexdigest() == EXPECTED_RANK_SHA256,
        "family rank",
    )
    stopping_rank = next(
        rank for rank, family in enumerate(ranked) if family[1] >= 2
    )
    require(stopping_rank == 2, stopping_rank)
    selected_families = ranked[: stopping_rank + 1]
    selected_indices = tuple(
        index for family in selected_families for index in family[4]
    )
    require(selected_indices == EXPECTED_SELECTED_INDICES, selected_indices)
    family_rank = {
        index: rank + 1
        for rank, family in enumerate(selected_families)
        for index in family[4]
    }
    selected_packet = tuple(
        (
            index, rows[index][0], gcd(LEVEL, rows[index][1]), rows[index][1],
            components[index], rows[index][1] * components[index],
            rows[index][3], family_rank[index],
        )
        for index in selected_indices
    )
    require(sum(item[5] for item in selected_packet) == EXPECTED_SELECTED_COST,
            "selected cost")
    require(
        hashlib.sha256(repr(selected_packet).encode()).hexdigest()
        == EXPECTED_SELECTED_PACKET_SHA256,
        "selected packet",
    )

    screen_tasks = tuple(
        (index, (LEVEL, *rows[index])) for index in selected_indices
    )
    if args.processes == 1:
        screened = tuple(screen_worker(task) for task in screen_tasks)
    else:
        with mp.get_context("spawn").Pool(
            min(args.processes, len(screen_tasks))
        ) as pool:
            screened = tuple(pool.map(screen_worker, screen_tasks, chunksize=1))
    require(tuple(item[0] for item in screened) == selected_indices, "screen order")
    canonical = []
    direct = 0
    legacy = 0
    count_by_index = {}
    for index, row, direct_count, legacy_count in screened:
        source = rows[index]
        divisor_gcd = gcd(LEVEL, source[1])
        require(
            row[:6]
            == (
                LEVEL, source[0], source[1], source[2],
                source[1] // divisor_gcd, source[3],
            ),
            (index, row[:6]),
        )
        require(row[16] == row[11], (index, "status verification"))
        require(direct_count + legacy_count == row[11], (index, "Farkas count"))
        require(row[12] == 0 and row[13] == (), (index, "residual"))
        canonical.append((index, tuple(row[:19])))
        direct += direct_count
        legacy += legacy_count
        count_by_index[index] = tuple(row[position] for position in (9, 10, 11, 12))
    canonical = tuple(canonical)
    selected_rows = tuple(row for _index, row in canonical)
    selected_totals = totals(selected_rows)
    selected_farkas = (direct, legacy)
    screen_sha = hashlib.sha256(repr(canonical).encode()).hexdigest()
    residual = tuple((index, row[13]) for index, row in canonical if row[12])
    residual_sha = hashlib.sha256(repr(residual).encode()).hexdigest()
    require(selected_totals == EXPECTED_SELECTED_TOTALS, selected_totals)
    require(selected_farkas == EXPECTED_SELECTED_FARKAS, selected_farkas)
    require(screen_sha == EXPECTED_SELECTED_SCREEN_SHA256, screen_sha)
    require(not residual and residual_sha == EXPECTED_EMPTY_RESIDUAL_SHA256,
            residual_sha)

    taxonomy_tasks = tuple(
        (index, (LEVEL, *rows[index]), count_by_index[index])
        for index in selected_indices
    )
    if args.processes == 1:
        taxonomy_rows = tuple(taxonomy_worker(task) for task in taxonomy_tasks)
    else:
        with mp.get_context("spawn").Pool(
            min(args.processes, len(taxonomy_tasks))
        ) as pool:
            taxonomy_rows = tuple(
                pool.map(taxonomy_worker, taxonomy_tasks, chunksize=1)
            )
    require(tuple(row[0] for row in taxonomy_rows) == selected_indices,
            "taxonomy order")
    taxonomy_by_row = tuple((row[0], row[2]) for row in taxonomy_rows)
    taxonomy_totals = (
        sum(row[2][0] for row in taxonomy_rows),
        sum(row[2][1] for row in taxonomy_rows),
        sum(row[2][2] for row in taxonomy_rows),
        sum(row[2][3] for row in taxonomy_rows),
        sum(row[3] for row in taxonomy_rows),
        sum(row[4] for row in taxonomy_rows),
        sum(row[5] for row in taxonomy_rows),
        sum(sum(row[2]) for row in taxonomy_rows),
    )
    taxonomy_packet = tuple(
        record for row in taxonomy_rows for record in row[6]
    )
    weighted_core = tuple(
        record for row in taxonomy_rows for record in row[7]
    )
    taxonomy_sha = hashlib.sha256(repr(taxonomy_packet).encode()).hexdigest()
    weighted_sha = hashlib.sha256(repr(weighted_core).encode()).hexdigest()
    if EXPECTED_TAXONOMY_TOTALS is not None:
        require(taxonomy_totals == EXPECTED_TAXONOMY_TOTALS, taxonomy_totals)
    if EXPECTED_TAXONOMY_BY_ROW is not None:
        require(taxonomy_by_row == EXPECTED_TAXONOMY_BY_ROW, taxonomy_by_row)
    if EXPECTED_TAXONOMY_SHA256 is not None:
        require(taxonomy_sha == EXPECTED_TAXONOMY_SHA256, taxonomy_sha)
    if EXPECTED_WEIGHTED_CORE_SHA256 is not None:
        require(weighted_sha == EXPECTED_WEIGHTED_CORE_SHA256, weighted_sha)

    require(LEDGER_BEFORE - len(selected_indices) == LEDGER_AFTER, "ledger")
    require(WALL_BEFORE - len(selected_indices) == WALL_AFTER, "wall count")
    semantic_packet = (
        FIRST_SOURCE_SHA256, FIRST_OUTPUT_SHA256, FIRST_SEMANTIC_SHA256,
        NATURAL_SOURCE_SHA256, census, live_packet, ranked, selected_packet,
        canonical, selected_totals, selected_farkas, residual_sha,
        taxonomy_by_row, taxonomy_totals, taxonomy_packet, weighted_core,
        LEDGER_BEFORE, LEDGER_AFTER, WALL_BEFORE, WALL_AFTER,
    )
    semantic = hashlib.sha256(repr(semantic_packet).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, semantic)

    lines = [
        "LRC14 projected-k3 z216 second nontrivial ruler cost-prefix closure",
        (
            f"dependency=first_prefix:{FIRST_SOURCE_SHA256}/"
            f"{FIRST_OUTPUT_SHA256}/{FIRST_SEMANTIC_SHA256};"
            f"THM3281_source:{NATURAL_SOURCE_SHA256}"
        ),
        (
            f"universe=z1:{LEVEL};atlas_rows:{census[0]};wall:{census[1]};"
            f"order:{census[2]};prior_closed:105;live_wall:{census[3]};"
            f"live_families:{census[4]};row_sha256:{EXPECTED_ROW_SHA256};"
            f"live_sha256:{EXPECTED_LIVE_SHA256};rank_sha256:"
            f"{EXPECTED_RANK_SHA256}"
        ),
    ]
    for rank, family in enumerate(ranked[:5], start=1):
        lines.append(
            f"rank={rank};cost:{family[0]};rows:{family[1]};gcd:{family[2]};"
            f"L:{family[3]};indices:{','.join(map(str, family[4]))}"
        )
    lines.append(
        "selection=cost_ordered_complete_family_prefix_through_second_"
        f"non_singleton;families:3;rows:{len(selected_indices)};"
        f"indices:{','.join(map(str, selected_indices))};"
        f"cost:{EXPECTED_SELECTED_COST};packet_sha256:"
        f"{EXPECTED_SELECTED_PACKET_SHA256}"
    )
    for family in selected_families:
        family_counts = tuple(
            sum(count_by_index[index][position] for index in family[4])
            for position in range(4)
        )
        lines.append(
            f"family=gcd{family[2]}/L{family[3]};rows:{family[1]};"
            f"cost:{family[0]};states:{family_counts[0]};"
            f"crude:{family_counts[1]};status:{family_counts[2]};"
            f"residual:{family_counts[3]}"
        )
    lines.extend(
        [
            (
                f"selected_screen=states:{selected_totals[0]};"
                f"crude:{selected_totals[1]};status:{selected_totals[2]};"
                f"residual:{selected_totals[3]};direct_farkas:{direct};"
                f"legacy_farkas:{legacy};screen_sha256:{screen_sha};"
                f"residual_sha256:{residual_sha}"
            ),
            (
                f"status_taxonomy=coordinate_union:{taxonomy_totals[0]};"
                f"zero_reduced_union:{taxonomy_totals[1]};"
                f"two_fan:{taxonomy_totals[2]};weighted_core:{taxonomy_totals[3]};"
                f"union_certificates:{taxonomy_totals[4]};"
                f"zero_reduced_certificates:{taxonomy_totals[5]};"
                f"fan_certificates:{taxonomy_totals[6]};"
                f"exact_Farkas_checked:{taxonomy_totals[7]};"
                f"by_row:{taxonomy_by_row};taxonomy_sha256:{taxonomy_sha};"
                f"weighted_core_sha256:{weighted_sha}"
            ),
            (
                "taxonomy_mechanism=coordinate_union_uses_tail_demand_greater_"
                "than_the_sum_of_named_marginals;zero_reduced_union_first_"
                "deletes_patterns_forced_off_by_a_zero_marginal;two_fan_uses_B_union_A_"
                "intersection_C_union_D_with_upper_bound_mB_plus_min_mA_"
                "mC_plus_mD;weighted_core_retains_exact_multirow_Farkas_checks"
            ),
            (
                f"consequence=ledger:{LEDGER_BEFORE}-{len(selected_indices)}="
                f"{LEDGER_AFTER};z216_wall:{WALL_BEFORE}-"
                f"{len(selected_indices)}={WALL_AFTER};order:0;cap:216"
            ),
            (
                "direction=empty_exact_upper_relaxations_exclude_the_selected_"
                "projected_wall_rows;screen_feasibility_is_not_used_as_a_"
                "converse_and_family_cost_is_only_a_work_order_not_a_safety_"
                "invariant"
            ),
            (
                "hostiles=the_pinned_first_prefix_replays_a_35_status_positive_"
                "control_an_8_residual_negative_control_and_a_rejected_negative_"
                "Farkas_multiplier;this_continuation_preserves_all_filters_and_"
                "checks_every_one_of_271_duals_exactly"
            ),
            (
                "stopping_boundary=357_z216_wall_rows_in_33_families_remain;"
                "next_cost_prefix_is_gcd72_L72072_singleton_then_gcd24_L25872_"
                "three_row_family"
            ),
            (
                "scope=projected_k3_z216_necessary_wall_rows_only;no_terminal_"
                "probe_divisor_action_common_carrier_phase_physical_cover_"
                "arbitrary_k_le_1_rung_or_LRC14_claim"
            ),
            f"semantic_sha256={semantic}",
            "all_exact_controls=PASS",
        ]
    )
    payload = "\n".join(lines) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    mp.freeze_support()
    main()
