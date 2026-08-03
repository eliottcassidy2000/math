#!/usr/bin/env python3
"""Exact z216 cost-prefix probe through the first non-singleton ruler family.

The maintained projected-k3 atlas is only a necessary quotient.  This probe
inherits its exact upper-relaxation screen, reconstructs the current 380-row
wall remainder, ranks every complete ``(gcd(216,L),L)`` family by ``sum L*r``,
and screens the cost-ordered prefix ending at the first family with more than
one row.  Persisted hashes omit solver-selected Farkas representatives.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
from collections import Counter, defaultdict
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_3281 = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z216_natural_wall_family_screen_descent_thm3281.py"
)
OUTPUT_3281 = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z216_natural_wall_family_screen_descent_thm3281.out"
)
SOURCE_3270 = ROOT / (
    "04-computation/lrc14_j7_k3_z216_order_row_screen_descent_thm3270.py"
)
COSTLY_REPORT = ROOT / (
    "05-knowledge/results/lrc14-z216-costly-gcd8-closure-opus-20260803.md"
)
COSTLY_AUDIT_SOURCE = ROOT / (
    ".scratch/lrc14_z216_gcd8_independent_audit_gcd8_independent_audit.py"
)
COSTLY_AUDIT_OUTPUT = ROOT / (
    ".scratch/lrc14_z216_gcd8_independent_audit_gcd8_independent_audit.out"
)

SOURCE_3281_SHA256 = "430dee7ba03e0d5c9ae0df72ac512500de4f7056cb4663d1c8468bfb93a49bfe"
OUTPUT_3281_SHA256 = "b9c23b21bf9766efbbc14aa97e2bb4268ddb6abb09b54a2b7424b8744ff686b2"
SEMANTIC_3281_SHA256 = "2330ce22ee3b318da26e70ce2c6c01b020c8f2b78bd958d529640a599514cd32"
COSTLY_REPORT_SHA256 = "2dfdaebf4ccd4b051d6839450c9039db4b2112d5839470a65ccc067f7b9d8fd5"
COSTLY_AUDIT_SOURCE_SHA256 = (
    "e6d6c8428c6216e27d3e0581231884bbfa74dc9bf5c1311314d2c37b7ee6e35e"
)
COSTLY_AUDIT_OUTPUT_SHA256 = (
    "1432ffaca0025137da91ed392313c0a9da48d6d657c7b6faefcc157660c5e92a"
)

LEVEL = 216
EXPECTED_CENSUS = (480, 447, 33, 380, 39)
EXPECTED_ROW_SHA256 = "53db9e1d3df2cf2b0398847682d909da81705e43a53ae2553d102fd152337649"
EXPECTED_LIVE_SHA256 = "23378edbda9d5a7d26ea0869820a541620aaf9f8d278383fa6407e6d2239fc99"
EXPECTED_STRATA = ((8, 19), (18, 1), (24, 135), (36, 15), (72, 310))
EXPECTED_LIVE_BY_GCD = (
    (24, 12, 91, 319_446_960),
    (36, 2, 11, 156_008_160),
    (72, 25, 278, 1_609_301_232),
)
EXPECTED_RANK_SHA256 = "7b58989c9b1cb66c28228393710d7204f417ccee54cc3fea426cc8b682b2b2b7"
EXPECTED_RANK_HEAD = (
    (
        243_936,
        1,
        72,
        11_088,
        (210,),
        ((210, (1, 6, 8, 9, 11, 12), 22, 243_936),),
    ),
    (
        443_520,
        1,
        24,
        18_480,
        (50,),
        ((50, (1, 2, 6, 8, 10, 11), 24, 443_520),),
    ),
    (
        517_440,
        3,
        24,
        5_880,
        (18, 140, 299),
        (
            (18, (1, 2, 4, 6, 10, 14), 26, 152_880),
            (140, (1, 4, 6, 10, 12, 14), 30, 176_400),
            (299, (2, 4, 6, 10, 12, 14), 32, 188_160),
        ),
    ),
)
EXPECTED_SELECTED_INDICES = (210, 50, 18, 140, 299)
EXPECTED_SELECTED_COST = 1_204_896
EXPECTED_SELECTED_PACKET_SHA256 = (
    "4c0b70020eaea2db08cb6a6995fcdef0f0e991ae686c32855876e943a58ea18c"
)
EXPECTED_TRIANGLE = (
    (4, 6, 10, 14),
    (1, 2, 12),
    ((1, 2), (1, 12), (2, 12)),
)
EXPECTED_FAMILY_TOTALS = (
    ((72, 11_088), (8, 8, 0, 0), (0, 0)),
    ((24, 18_480), (22, 21, 1, 0), (0, 1)),
    ((24, 5_880), (12, 9, 3, 0), (0, 3)),
)
EXPECTED_SELECTED_TOTALS = (42, 38, 4, 0)
EXPECTED_SELECTED_FARKAS = (0, 4)
EXPECTED_ELEMENTARY_UNION_CERTIFICATES = (
    (50, ((3, 10, (1, 3), 896, 528, 368),)),
    (
        18,
        (
            (2, 11, (0, 1, 3), 486, 252, 234),
            (3, 11, (0, 1, 3), 272, 252, 20),
        ),
    ),
    (140, ((2, 7, (0, 1, 2), 608, 420, 188),)),
    (299, ((4, 15, (0, 1, 2, 3), 80, 60, 20),)),
)
EXPECTED_SELECTED_SCREEN_SHA256 = (
    "7a01d44843af06f0c8054926baf8007363bf309a316c10a5bc6b429135f044cc"
)
EXPECTED_EMPTY_RESIDUAL_SHA256 = (
    "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d"
)
CONTROL_INDICES = (14, 64)
EXPECTED_CONTROL_TOTALS = (
    (14, (59, 24, 35, 0), (0, 35), 0),
    (64, (115, 50, 57, 8), (0, 57), 8),
)
EXPECTED_CONTROL_SCREEN_SHA256 = (
    "041c3668888c19a2d0689e127205f7983062bb6ef586f3acb103e2dbaac00c09"
)
EXPECTED_RESIDUAL_CONTROL_MASK_SHA256 = (
    "098c78ce3def914ca1c0ad26458558898fe624040222f429e6163b9a8e9855c4"
)
LEDGER_BEFORE = 373_184
LEDGER_AFTER = 373_179
WALL_BEFORE = 380
WALL_AFTER = 375
EXPECTED_SEMANTIC_SHA256 = (
    "44f7e8cf786ed428297f3f389214cab7c0b28158a856be63806a3394e594158f"
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


def screen_worker(item):
    index, task = item
    natural = load(f"ruler_prefix_natural_{index}", SOURCE_3281)
    return natural.screen_worker((index, task))


def status_engine(natural):
    """Recover the pinned exact status data engine through the promoted chain."""
    entry = load("ruler_prefix_status_entry", natural.SOURCE_3139)
    source = entry.load("ruler_prefix_status_source")
    driver = source.load("ruler_prefix_status_driver")
    predecessor = driver.load("ruler_prefix_status_predecessor")
    bridge = predecessor.load("ruler_prefix_status_bridge")
    base = bridge.load("ruler_prefix_status_base")
    return base.thm.screen_engine.eng


def totals(rows):
    return tuple(sum(row[position] for row in rows) for position in (9, 10, 11, 12))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=3)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    require(args.processes >= 1, args.processes)

    check_dependency(SOURCE_3281, SOURCE_3281_SHA256)
    check_dependency(OUTPUT_3281, OUTPUT_3281_SHA256, SEMANTIC_3281_SHA256)
    check_dependency(COSTLY_REPORT, COSTLY_REPORT_SHA256)
    check_dependency(COSTLY_AUDIT_SOURCE, COSTLY_AUDIT_SOURCE_SHA256)
    check_dependency(COSTLY_AUDIT_OUTPUT, COSTLY_AUDIT_OUTPUT_SHA256)
    audit_text = COSTLY_AUDIT_OUTPUT.read_text(encoding="utf-8")
    require("all_exact_controls=PASS" in audit_text, "costly gcd8 audit status")
    require("negative_Farkas_alpha_rejected" in audit_text, "Farkas hostile missing")

    natural = load("ruler_prefix_natural_main", SOURCE_3281)
    order = load("ruler_prefix_order_main", SOURCE_3270)
    rows, components = natural.atlas_rows()
    require(hashlib.sha256(repr(rows).encode()).hexdigest() == EXPECTED_ROW_SHA256,
            "row order")
    strata = Counter(gcd(LEVEL, row[1]) for row in rows)
    require(tuple(sorted(strata.items())) == EXPECTED_STRATA, strata)

    gcd8 = {index for index, row in enumerate(rows) if gcd(LEVEL, row[1]) == 8}
    gcd18 = {index for index, row in enumerate(rows) if gcd(LEVEL, row[1]) == 18}
    order_rows = set(order.EXPECTED_ORDER_INDICES)
    natural_rows = set(natural.EXPECTED_UNION_INDICES)
    require(
        not (
            gcd8.intersection(gcd18)
            or gcd8.intersection(order_rows)
            or gcd8.intersection(natural_rows)
            or gcd18.intersection(order_rows)
            or gcd18.intersection(natural_rows)
            or order_rows.intersection(natural_rows)
        ),
        "prior closures overlap",
    )
    closed = gcd8 | gcd18 | order_rows | natural_rows
    require(len(closed) == 100, len(closed))
    live = tuple(
        index for index, row in enumerate(rows) if row[3] and index not in closed
    )
    census = (
        len(rows),
        sum(row[3] for row in rows),
        sum(not row[3] for row in rows),
        len(live),
        0,
    )

    families = defaultdict(list)
    for index in live:
        families[(gcd(LEVEL, rows[index][1]), rows[index][1])].append(index)
    census = (*census[:4], len(families))
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
            (
                index,
                rows[index][0],
                components[index],
                ruler * components[index],
            )
            for index in indices
        )
        ranked.append(
            (
                sum(item[3] for item in packet),
                len(indices),
                divisor_gcd,
                ruler,
                indices,
                packet,
            )
        )
    ranked = tuple(sorted(ranked))
    require(len(ranked) == 39, len(ranked))
    require(ranked[:3] == EXPECTED_RANK_HEAD, ranked[:3])
    require(
        hashlib.sha256(repr(ranked).encode()).hexdigest() == EXPECTED_RANK_SHA256,
        "family rank",
    )
    by_gcd = tuple(
        (
            divisor_gcd,
            sum(family[2] == divisor_gcd for family in ranked),
            sum(family[1] for family in ranked if family[2] == divisor_gcd),
            sum(family[0] for family in ranked if family[2] == divisor_gcd),
        )
        for divisor_gcd in (24, 36, 72)
    )
    require(by_gcd == EXPECTED_LIVE_BY_GCD, by_gcd)

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
            index,
            rows[index][0],
            gcd(LEVEL, rows[index][1]),
            rows[index][1],
            components[index],
            rows[index][1] * components[index],
            rows[index][3],
            family_rank[index],
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

    triangle_indices = EXPECTED_RANK_HEAD[2][4]
    common_core = set(rows[triangle_indices[0]][0])
    for index in triangle_indices[1:]:
        common_core.intersection_update(rows[index][0])
    outer = set().union(*(set(rows[index][0]) for index in triangle_indices))
    outer.difference_update(common_core)
    triangle = (
        tuple(sorted(common_core)),
        tuple(sorted(outer)),
        tuple(
            tuple(sorted(set(rows[index][0]).difference(common_core)))
            for index in triangle_indices
        ),
    )
    require(triangle == EXPECTED_TRIANGLE, triangle)

    all_indices = selected_indices + CONTROL_INDICES
    tasks = tuple((index, (LEVEL, *rows[index])) for index in all_indices)
    if args.processes == 1:
        screened = tuple(screen_worker(task) for task in tasks)
    else:
        with mp.get_context("spawn").Pool(min(args.processes, len(tasks))) as pool:
            screened = tuple(pool.map(screen_worker, tasks, chunksize=1))
    require(tuple(item[0] for item in screened) == all_indices, "screen order")

    canonical = {}
    farkas = {}
    for index, row, direct, legacy in screened:
        source = rows[index]
        divisor_gcd = gcd(LEVEL, source[1])
        require(
            row[:6]
            == (
                LEVEL,
                source[0],
                source[1],
                source[2],
                source[1] // divisor_gcd,
                source[3],
            ),
            (index, row[:6]),
        )
        require(row[16] == row[11], (index, "status verification"))
        require(direct + legacy == row[11], (index, "Farkas count"))
        canonical[index] = tuple(row)
        farkas[index] = (direct, legacy)

    selected_rows = tuple(canonical[index] for index in selected_indices)
    selected_totals = totals(selected_rows)
    selected_farkas = (
        sum(farkas[index][0] for index in selected_indices),
        sum(farkas[index][1] for index in selected_indices),
    )
    selected_screen = tuple(
        (index, tuple(canonical[index][:19])) for index in selected_indices
    )
    selected_screen_sha = hashlib.sha256(repr(selected_screen).encode()).hexdigest()
    selected_residual = tuple(
        (index, canonical[index][13])
        for index in selected_indices
        if canonical[index][12]
    )
    selected_residual_sha = hashlib.sha256(
        repr(selected_residual).encode()
    ).hexdigest()
    require(selected_totals == EXPECTED_SELECTED_TOTALS, selected_totals)
    require(selected_farkas == EXPECTED_SELECTED_FARKAS, selected_farkas)
    require(selected_screen_sha == EXPECTED_SELECTED_SCREEN_SHA256,
            selected_screen_sha)
    require(not selected_residual, selected_residual)
    require(selected_residual_sha == EXPECTED_EMPTY_RESIDUAL_SHA256,
            selected_residual_sha)

    family_totals = []
    for family in selected_families:
        indices = family[4]
        rows_here = tuple(canonical[index] for index in indices)
        family_totals.append(
            (
                (family[2], family[3]),
                totals(rows_here),
                (
                    sum(farkas[index][0] for index in indices),
                    sum(farkas[index][1] for index in indices),
                ),
            )
        )
    family_totals = tuple(family_totals)
    require(family_totals == EXPECTED_FAMILY_TOTALS, family_totals)

    # Each non-crude selected row has a single exact status state.  Rebuild a
    # solver-free explanation: at a load threshold, the capacity-good binary
    # patterns are exactly the union of a named coordinate subset, while the
    # histogram tail demand exceeds the sum of those coordinate marginals.
    # This is the elementary union bound underlying the corresponding Farkas
    # contradiction, and is independent of the optimizer's chosen dual basis.
    eng = status_engine(natural)
    eng.FIRST = LEVEL
    eng.ray.FIRST = LEVEL
    elementary_certificates = []
    for index in selected_indices:
        row = canonical[index]
        if not row[11]:
            require(row[17] is None, (index, "unexpected representative"))
            continue
        require(row[11] == 1 and row[17] is not None,
                (index, "non-unique status state"))
        first, body, divisors, instance = row[17]
        require((first, body) == (LEVEL, rows[index][0]), (index, row[17]))
        q, cofactor, marginals, capacity_set, histogram = instance
        divisor_lcm = lcm(*divisors)
        require(divisor_lcm == q * cofactor, (index, "status cofactor"))
        rebuilt_marginals, capacities = eng.ray.local.hunter_status_data(
            divisor_lcm, divisors, q
        )
        require(rebuilt_marginals == marginals, (index, "status marginals"))
        require(tuple(sorted(set(capacities))) == capacity_set,
                (index, "status capacities"))
        certificates = []
        for threshold, _count in histogram:
            if threshold <= 0:
                continue
            demand = sum(
                count for load_value, count in histogram
                if load_value >= threshold
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
                        (
                            threshold,
                            coordinate_mask,
                            coordinates,
                            demand,
                            supply,
                            demand - supply,
                        )
                    )
        require(certificates, (index, "no elementary union certificate"))
        elementary_certificates.append((index, tuple(certificates)))
    elementary_certificates = tuple(elementary_certificates)
    require(
        elementary_certificates == EXPECTED_ELEMENTARY_UNION_CERTIFICATES,
        elementary_certificates,
    )

    controls = tuple(
        (index, tuple(canonical[index][:19])) for index in CONTROL_INDICES
    )
    control_sha = hashlib.sha256(repr(controls).encode()).hexdigest()
    require(control_sha == EXPECTED_CONTROL_SCREEN_SHA256, control_sha)
    for index, expected, expected_farkas, residual in EXPECTED_CONTROL_TOTALS:
        require(totals((canonical[index],)) == expected, (index, "control totals"))
        require(farkas[index] == expected_farkas, (index, "control Farkas"))
        require(canonical[index][12] == residual, (index, "control residual"))
    residual_control = canonical[64][13]
    require(
        hashlib.sha256(repr(residual_control).encode()).hexdigest()
        == EXPECTED_RESIDUAL_CONTROL_MASK_SHA256,
        "residual control masks",
    )
    residual_quotient = Counter(
        lcm(*mask) // canonical[64][4] for mask in residual_control
    )
    require(tuple(sorted(residual_quotient.items())) == ((8, 8),),
            residual_quotient)

    require(LEDGER_BEFORE - len(selected_indices) == LEDGER_AFTER, "ledger")
    require(WALL_BEFORE - len(selected_indices) == WALL_AFTER, "wall count")
    semantic_packet = (
        SOURCE_3281_SHA256,
        OUTPUT_3281_SHA256,
        SEMANTIC_3281_SHA256,
        COSTLY_REPORT_SHA256,
        COSTLY_AUDIT_SOURCE_SHA256,
        COSTLY_AUDIT_OUTPUT_SHA256,
        EXPECTED_CENSUS,
        EXPECTED_ROW_SHA256,
        EXPECTED_STRATA,
        live_packet,
        ranked,
        by_gcd,
        selected_packet,
        triangle,
        selected_screen,
        selected_totals,
        selected_farkas,
        elementary_certificates,
        selected_residual_sha,
        controls,
        tuple(sorted(residual_quotient.items())),
        LEDGER_BEFORE,
        LEDGER_AFTER,
        WALL_BEFORE,
        WALL_AFTER,
    )
    semantic = hashlib.sha256(repr(semantic_packet).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, semantic)

    lines = [
        "LRC14 projected-k3 z216 first nontrivial ruler cost-prefix closure",
        (
            "dependencies="
            f"THM3281:{SOURCE_3281_SHA256}/{OUTPUT_3281_SHA256}/"
            f"{SEMANTIC_3281_SHA256};costly_g8_report:{COSTLY_REPORT_SHA256};"
            f"costly_g8_audit:{COSTLY_AUDIT_SOURCE_SHA256}/"
            f"{COSTLY_AUDIT_OUTPUT_SHA256}"
        ),
        (
            f"universe=z1:{LEVEL};atlas_rows:{census[0]};wall:{census[1]};"
            f"order:{census[2]};prior_closed:100;live_wall:{census[3]};"
            f"live_families:{census[4]};row_sha256:{EXPECTED_ROW_SHA256};"
            f"live_sha256:{EXPECTED_LIVE_SHA256};strata:{EXPECTED_STRATA}"
        ),
        (
            f"live_by_gcd={by_gcd};total_invoice:"
            f"{sum(item[3] for item in by_gcd)};rank_sha256:{EXPECTED_RANK_SHA256}"
        ),
    ]
    for rank, family in enumerate(ranked, start=1):
        lines.append(
            f"rank={rank};cost:{family[0]};rows:{family[1]};gcd:{family[2]};"
            f"L:{family[3]};indices:{','.join(map(str, family[4]))}"
        )
    lines.extend(
        [
            (
                "selection=cost_ordered_complete_family_prefix_through_first_"
                f"non_singleton;families:3;rows:{len(selected_indices)};"
                f"indices:{','.join(map(str, selected_indices))};"
                f"cost:{EXPECTED_SELECTED_COST};"
                f"packet_sha256:{EXPECTED_SELECTED_PACKET_SHA256}"
            ),
            (
                f"triangle=common_core:{triangle[0]};outer_triple:{triangle[1]};"
                f"two_subsets:{triangle[2]}"
            ),
        ]
    )
    for family, summary in zip(selected_families, family_totals):
        lines.append(
            f"family=gcd{family[2]}/L{family[3]};rows:{family[1]};"
            f"indices:{','.join(map(str, family[4]))};cost:{family[0]};"
            f"states:{summary[1][0]};crude:{summary[1][1]};"
            f"status:{summary[1][2]};residual:{summary[1][3]};"
            f"direct_farkas:{summary[2][0]};legacy_farkas:{summary[2][1]}"
        )
    lines.extend(
        [
            (
                f"selected_screen=states:{selected_totals[0]};"
                f"crude:{selected_totals[1]};status:{selected_totals[2]};"
                f"residual:{selected_totals[3]};"
                f"direct_farkas:{selected_farkas[0]};"
                f"legacy_farkas:{selected_farkas[1]};"
                f"screen_sha256:{selected_screen_sha};"
                f"residual_sha256:{selected_residual_sha}"
            ),
            (
                "elementary_union_certificates="
                f"{elementary_certificates};mechanism:capacity_good_patterns_"
                "equal_a_coordinate_union_but_histogram_tail_demand_exceeds_"
                "the_sum_of_its_marginals"
            ),
            (
                "controls=prior_gcd24_L2352_index14:"
                "59=24crude+35status+0residual;"
                "known_gcd8_index64:115=50crude+57status+8residual;"
                f"residual_quotient:{tuple(sorted(residual_quotient.items()))};"
                f"control_sha256:{control_sha};"
                f"residual_mask_sha256:{EXPECTED_RESIDUAL_CONTROL_MASK_SHA256};"
                "negative_Farkas_alpha_mutation:rejected_by_pinned_independent_audit"
            ),
            (
                f"consequence=ledger:{LEDGER_BEFORE}-{len(selected_indices)}="
                f"{LEDGER_AFTER};z216_wall:{WALL_BEFORE}-"
                f"{len(selected_indices)}={WALL_AFTER};order:0;cap:216"
            ),
            (
                "direction=the_THM3139_ray_screen_is_a_necessary_upper_"
                "relaxation;empty_enlarged_screens_exclude_the_five_selected_"
                "projected_wall_rows;no_converse_from_screen_feasibility_is_used"
            ),
            (
                "evidence_boundary=all_four_status_certificates_are_checked_over_"
                "exact_rationals_but_persisted_hashes_bind_only_canonical_"
                "row_prefixes_and_counts_not_solver_selected_duals_or_magnitudes"
            ),
            (
                "scope=cost_prefix_of_three_complete_intrinsic_ruler_families_in_"
                "the_projected_k3_z216_necessary_wall_remainder_only;the_other_"
                "375_wall_rows_arbitrary_k_le_1_rung_physical_cover_and_LRC14_"
                "remain_open"
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
