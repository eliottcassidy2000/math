#!/usr/bin/env python3
r"""Exact child top-four census on every open THM-2904 centre.

THM-2904 allocates a hypothetical five-cover of a marked carrier C to the
earliest maximum-coverage label x in its finite hostile-centre core.  On
that branch the remaining four labels must cover

    R = C minus D_x

and must avoid the inherited gate prefix, x itself, and every earlier
hostile centre.  This verifier reconstructs all 51,222 such marked children
and finds their exact global top four allowed singleton coverages.

The literal head begins at 15 and is initially scanned through 2,000.  If
the fourth retained value q_4 exceeds h_R/7, THM-735 gives

    c_R(w) < h_R/7 + gamma/w,       gamma = (99/70) r_R / 7.

Consequently

    N_0 = ceil(gamma/(q_4-h_R/7)) - 1

is a sufficient finite horizon.  The verifier scans through
max(2,000,N_0), recomputes q_4, and checks the final omitted-label majorant
including its equality boundary.  No arbitrary search cap is used.

The resulting child closures are composed branchwise with the already
proved THM-2904 fully closed rows, THM-2911 finite-H1 rows, THM-2907
pair-exception rows, and THM-2913 one-H3 toothpick rows.  Root projection
and exact set difference are taken both against the canonical 314-root
union through THM-2912 and the 351-root union through THM-2913.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.util
import multiprocessing as mp
import os
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PARENT_SOURCE = (
    ROOT
    / "04-computation/"
    "lrc14_j6_all_hard_ranked_h1_hunter_pivot_census_codex_20260729.py"
)
PARENT_LEDGER = (
    ROOT
    / "05-knowledge/results/lrc14_j6_ranked_h1_thm2911/"
    "thm2904_hostile_centre.ledger.out"
)
THM2911_VERIFY = (
    ROOT / "04-computation/lrc14_j6_ranked_h1_thm2911/verify.py"
)
THM2911_OUT = (
    ROOT / "05-knowledge/results/lrc14_j6_ranked_h1_thm2911/locked.out"
)
THM2907_OUT = (
    ROOT
    / "05-knowledge/results/"
    "lrc14_j6_paircap_exception_h4_link_child_census_codex_20260729.out"
)
THM2907_ENDPOINT_OUT = (
    ROOT
    / "05-knowledge/results/"
    "lrc14_j6_paircap_exception_endpoint_convention_audit_codex_20260729.out"
)
THM2912_SOURCE = (
    ROOT
    / "04-computation/"
    "lrc14_j6_one_h3_row_child_top4_census_codex_20260729.py"
)
THM2912_OUT = (
    ROOT
    / "05-knowledge/results/"
    "lrc14_j6_one_h3_row_child_top4_census_codex_20260729.out"
)
THM2913_SOURCE = (
    ROOT
    / "04-computation/"
    "lrc14_j6_one_h3_row_pair_hunter_toothpick_closure_codex_20260729.py"
)
THM2913_OUT = (
    ROOT
    / "05-knowledge/results/"
    "lrc14_j6_one_h3_row_pair_hunter_toothpick_closure_codex_20260729.out"
)

# These are refreshed only after a clean LF-policy replay.
DEPENDENCY_SHA256: dict[Path, str] = {
    PARENT_SOURCE: "644104b0de90654466e75c6531109736b0445aadb357eee2413e8787ac3a53fa",
    PARENT_LEDGER: "bec35518329b5d9e6ba2c9a8c87bfb20234a0c07dc1a5c5f2babec21888d452a",
    THM2911_VERIFY: "e0ac67539f7ff09376645a62beef0a9ac7d0766a2e749666f94d1fd4d6487b15",
    THM2911_OUT: "e5c58cc2eb325928c00839c2593450ea7cce8945b3835898ec83c6c5f42fac9b",
    THM2907_OUT: "1df929e106cd16c094886d59f3702ba9bafa395ee906527fe4592a1552e9b458",
    THM2907_ENDPOINT_OUT: "a93a4a724dac6c55806f3358c2f5ab25de8f0261c92906a0161414781a717d20",
    THM2912_SOURCE: "d2810560a7d002d7eeadecc6a50a7733c90585527295aa5e85e72775739b839b",
    THM2912_OUT: "454d87c8beeb81405b031cce4b40bdda0d385cfcd9c48e6fcf4eb810cfc00c5a",
    THM2913_SOURCE: "14e56e124197cd1bdae841efa195a1e7c282e7ea368a610e5f4d56509431858b",
    THM2913_OUT: "3604644a9691b13e7fa245249b68c9586ec2775996834f7761f32eb0b89f3e64",
}

FIRST_EXTERNAL = 15
INITIAL_HEAD_LIMIT = 2_000
OFFICIAL_ROOT_COUNT = 3_432
THM2911_MAX_CUTOFF = 15_000

# Locked after the first theorem-grade full replay.
EXPECTED_SUMMARY: dict[str, object] | None = {
    "parent_rows": 11_842,
    "pair_exception_rows": 52,
    "hostile_centres": 55_293,
    "THM2904_closed_pivots": 4_071,
    "open_pivots": 51_222,
    "parent_closed_rows": 279,
    "open_parent_rows": 11_563,
    "closed_children": 46_356,
    "failed_children": 4_866,
    "equality_children": 0,
    "sealed_children": 51_222,
    "extended_children": 1,
    "maximum_head_limit": 2_134,
    "maximum_final_required": 2_134,
    "T_rows": 8_112,
    "remaining_open_rows": 3_451,
    "H_rows": 3_090,
    "E_rows": 52,
    "P_rows": 279,
    "Q_rows": 38,
    "branch_atoms": (
        ("-", 3_375),
        ("E", 52),
        ("H", 24),
        ("HP", 215),
        ("HT", 2_851),
        ("P", 64),
        ("T", 5_261),
    ),
    "branch_atoms_with_Q": (
        ("-", 3_338),
        ("E", 52),
        ("H", 23),
        ("HP", 215),
        ("HQ", 1),
        ("HT", 2_851),
        ("P", 64),
        ("Q", 37),
        ("T", 5_261),
    ),
    "H_intersection_T": 2_851,
    "H_intersection_P": 215,
    "E_intersection_T": 0,
    "Q_intersection_H": 1,
    "Q_intersection_T": 0,
    "survivor_bodies": 3_411,
    "no_G5_bodies": 21,
    "P_roots": 10,
    "live_route_bodies": 3_401,
    "live_route_roots": 1_219,
    "live_route_failures": 2_182,
    "pure_roots": 1_229,
    "pure_overlap314": 277,
    "pure_additive314": 952,
    "pure_union314": 1_266,
    "pure_residual314": 2_166,
    "EHP_composed_roots": 1_248,
    "EHP_overlap314": 284,
    "EHP_additive314": 964,
    "EHP_union314": 1_278,
    "EHP_residual314": 2_154,
    "E_additive_gain": 6,
    "H_additive_gain": 6,
    "mixed_EH_gain": 0,
    "EHPQ_composed_roots": 1_285,
    "Q_root_gain": 37,
    "Q_extra_synergy": 0,
    "current_baseline": 351,
    "current_overlap": 321,
    "current_additive": 964,
    "current_union": 1_315,
    "current_residual": 2_117,
    "exception_children": 400,
    "exception_closed_children": 272,
    "exception_failed_children": 128,
    "exception_completed_rows": 0,
    "exception_fail_histogram": ((1, 5), (2, 26), (3, 14), (4, 6), (5, 1)),
    "one_H3_children": 807,
    "one_H3_closed_children": 765,
    "one_H3_failed_children": 42,
    "one_H3_slice_roots": 172,
    "one_H3_full_route_roots": 170,
    "one_H3_composed_roots": 172,
    "child_component_histogram": (
        (16, 3),
        (18, 38),
        (20, 322),
        (22, 1_210),
        (24, 2_633),
        (26, 4_209),
        (28, 5_580),
        (30, 6_309),
        (32, 6_225),
        (34, 5_665),
        (36, 4_886),
        (38, 3_953),
        (40, 2_971),
        (42, 2_251),
        (44, 1_543),
        (46, 1_049),
        (48, 776),
        (50, 510),
        (52, 358),
        (54, 249),
        (56, 179),
        (58, 108),
        (60, 66),
        (62, 42),
        (64, 26),
        (66, 21),
        (68, 13),
        (70, 6),
        (72, 10),
        (74, 2),
        (76, 2),
        (78, 2),
        (80, 2),
        (82, 1),
        (88, 2),
    ),
    "open_centres_per_parent": (
        (1, 737),
        (2, 1_428),
        (3, 1_822),
        (4, 2_033),
        (5, 2_183),
        (6, 1_629),
        (7, 969),
        (8, 446),
        (9, 206),
        (10, 82),
        (11, 20),
        (12, 7),
        (13, 1),
    ),
    "failed_children_per_failed_root": (
        (1, 939),
        (2, 572),
        (3, 333),
        (4, 151),
        (5, 73),
        (6, 58),
        (7, 24),
        (8, 14),
        (9, 6),
        (10, 5),
        (11, 3),
        (12, 3),
        (14, 1),
    ),
}
EXPECTED_ROOT_DIGESTS: dict[str, str] | None = {
    "live_route1219": "b25f51db595f01bd6c24db6f6d25ed9a230f6c7abf10841a1ef4aacdb793f371",
    "live_route_overlap267": "befd2efb71c9ecd4c73ef746dffd01a11aebe244692265080a1dac0e3ca2c40b",
    "pure_additive952": "b7e0783936dbaa4bc9b26e260528e8a46db236fbb32e8540657c5c35d03d11c7",
    "pure_union1266": "f6956e378cb953f831d87720541c7568a0a80509986588aedecd583800a6da23",
    "composed_roots1248": "86f997f438c3d70175375c0205cdb6f8fe80568391be49935e8298d360d7611e",
    "composed_overlap284": "b91ecc9601043c762e90989041757e39669dfc87326fb439b1cb3fa00a78a0ba",
    "composed_additive964": "4b2dcf6945aa80e8896f22115df5096f028dce5d9936b3148b0083c259657254",
    "composed_union1278": "ca4578e8638a0aaf978156ad07ada4b406fdd24fee1562acb63cf4fcab5ef496",
    "E_additive_gain6": "a9dc50a6955859191d0e77224490e2f853cb98fb094b33cabe1646d2e155f84d",
    "H_additive_gain6": "4f829186336283bf707859a774a2af5ae6e9e6514827f774242d79f88f5d6ab7",
    "Q_roots38": "60ff2ba4341a6a8c07947d3d053af7afb995f371ecd56f83bc4c08983be0526d",
    "final_roots1285": "e34db4c3fac3aea27031df52f5e2577d21c5fbf779967c7e24258dc44c73d481",
    "current_overlap321": "5aa24e243f98612be05e4e61b990e5aeb655c24ad69ebefbe9bf05f25c6a2dbd",
    "current_additive964": "4b2dcf6945aa80e8896f22115df5096f028dce5d9936b3148b0083c259657254",
    "current_union1315": "47c0b23646ae5744f5354d6475aa23283c275341e6665ba5532640cddea0c41f",
}
EXPECTED_OPEN_PIVOT_KEY_DIGEST: str | None = (
    "7da02bd0669380e5e9e2554b4f38118fefd99db6f5058c59b449d360816185c9"
)
EXPECTED_LEDGER_SEMANTIC_SHA256: str | None = (
    "274114167a4e173d242fd0d62b980593df5fffb6adfc9cf116ec171fe628b1ff"
)
EXPECTED_LEDGER_FILE_SHA256: str | None = (
    "798cd660ab60e2021b28074a1390af3f6b1367c99f2d0ab63a581513f7871071"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def ceiling(value: F) -> int:
    return (value.numerator + value.denominator - 1) // value.denominator


def tuple_ints(text: str) -> tuple[int, ...]:
    return () if not text else tuple(map(int, text.split(",")))


def unique_output_text(path: Path, prefix: str) -> str:
    matches = [
        line.removeprefix(prefix)
        for line in path.read_text().splitlines()
        if line.startswith(prefix)
    ]
    require(len(matches) == 1, f"{path.name}: expected one {prefix!r} line")
    return matches[0]


def literal_output_value(path: Path, prefix: str) -> object:
    return ast.literal_eval(unique_output_text(path, prefix))


def body_digest(values: set[tuple[int, ...]]) -> str:
    payload = "".join(
        ",".join(map(str, body)) + "\n"
        for body in sorted(values)
    )
    return hashlib.sha256(payload.encode()).hexdigest()


def key_digest(values: set[tuple[object, ...]], domain: bytes) -> str:
    digest = hashlib.sha256(domain + b"\n")
    for value in sorted(values):
        digest.update((repr(value) + "\n").encode())
    return digest.hexdigest()


def branch_key(row: dict[str, object]) -> tuple[object, ...]:
    return row["body"], row["rank"], row["apex"], row["prefix"]


def child_key(row: dict[str, object]) -> tuple[object, ...]:
    return (
        row["body"],
        row["rank"],
        row["apex"],
        row["prefix"],
        row["center"],
        row["earlier"],
    )


def load_parent():
    require(
        file_sha256(PARENT_SOURCE) == DEPENDENCY_SHA256[PARENT_SOURCE],
        "THM2904 source changed",
    )
    spec = importlib.util.spec_from_file_location("thm2915_parent", PARENT_SOURCE)
    require(spec is not None and spec.loader is not None, "cannot load THM2904")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


H = load_parent()


def check_dependencies() -> None:
    for path, expected in DEPENDENCY_SHA256.items():
        require(file_sha256(path) == expected, f"{path.name}: dependency hash changed")
    require(
        H.S2 == F(99, 70) and H.S2 * H.S2 > 2,
        "THM735 rational discrepancy majorant changed",
    )
    require(
        unique_output_text(THM2911_OUT, "mode=") == "LOCKED_FINITE_EXACT"
        and unique_output_text(THM2911_OUT, "all_exact_controls=") == "PASS"
        and "proved_union=181,residual=3251"
        in unique_output_text(THM2911_OUT, "current_composition="),
        "THM2911 locked composition changed",
    )
    h4_counts = literal_output_value(THM2907_OUT, "counts=")
    require(
        h4_counts[:10]
        == (18_290, 16_357, 1_933, 179, 529, 1_225, 18_280, 10, 51, 1)
        and h4_counts[24:26] == (52, 52)
        and unique_output_text(THM2907_OUT, "mode=") == "LOCKED"
        and unique_output_text(THM2907_OUT, "all_exact_controls=") == "PASS"
        and unique_output_text(THM2907_ENDPOINT_OUT, "mode=") == "LOCKED"
        and unique_output_text(THM2907_ENDPOINT_OUT, "all_exact_controls=")
        == "PASS",
        "THM2907 pair-exception discharge changed",
    )
    counts2912 = literal_output_value(THM2912_OUT, "counts=")
    require(
        counts2912[7:18]
        == (210, 807, 1_599_450, 20, 56, 807, 765, 0, 42, 172, 38)
        and counts2912[18:23] == (181, 39, 133, 314, 3_118)
        and unique_output_text(THM2912_OUT, "mode=") == "LOCKED"
        and unique_output_text(THM2912_OUT, "all_exact_controls=") == "PASS",
        "THM2912 canonical union changed",
    )
    counts2913 = literal_output_value(THM2913_OUT, "counts=")
    require(
        counts2913[:7] == (11_842, 11_511, 3_344, 210, 807, 765, 42)
        and counts2913[7:13] == (22, 20, 10, 10, 29, 29)
        and counts2913[13:20] == (38, 1, 37, 181, 314, 351, 3_081)
        and counts2913[26:29] == (0, 0, 0)
        and unique_output_text(THM2913_OUT, "mode=") == "LOCKED"
        and unique_output_text(THM2913_OUT, "all_exact_controls=") == "PASS",
        "THM2913 repaired toothpick closure changed",
    )


def parse_ranked_pairs(text: str, width: int) -> tuple[tuple[str, ...], ...]:
    if not text:
        return ()
    rows = tuple(tuple(item.split(":")) for item in text.split(","))
    require(all(len(row) == width for row in rows), "ranked field width changed")
    return rows


def parse_parent_ledger() -> list[dict[str, object]]:
    lines = PARENT_LEDGER.read_text().splitlines()
    require(
        lines[-2]
        == "ledger_sha256=ec878244b922ba5f48633614a86a1f9706c1fbdd0ebd6c61f020291cfd737bab"
        and lines[-1].endswith("not LRC14"),
        "THM2904 hostile-centre ledger footer changed",
    )
    rows: list[dict[str, object]] = []
    for line in lines:
        if not line.startswith("H1;"):
            continue
        fields = dict(part.split("=", 1) for part in line.split(";")[1:])
        core_rows = tuple(
            (F(value), int(label))
            for label, value in parse_ranked_pairs(fields["H1"], 2)
        )
        pivot_rows = tuple(
            (
                int(label),
                F(bound),
                F(margin),
                bool(int(closed)),
                bool(int(repaired)),
            )
            for label, bound, margin, closed, repaired
            in parse_ranked_pairs(fields["pivot"], 5)
        )
        require(
            tuple(label for _, label in core_rows)
            == tuple(label for label, *_ in pivot_rows),
            "THM2904 core/pivot order changed",
        )
        row: dict[str, object] = {
            "body": tuple_ints(fields["E"]),
            "rank": int(fields["rank"]),
            "apex": int(fields["a"]),
            "prefix": tuple_ints(fields["P"]),
            "stratum": fields["S"],
            "mass": F(fields["h"]),
            "pair_cap": F(fields["B2"]),
            "threshold": F(fields["lambda"]),
            "pair_exception": bool(int(fields["exception"])),
            "core_rows": core_rows,
            "pivot_rows": pivot_rows,
        }
        row["branch_closed"] = all(pivot[3] for pivot in pivot_rows)
        rows.append(row)
    require(
        len(rows) == 11_842
        and sum(len(row["core_rows"]) for row in rows) == 55_293
        and sum(
            not pivot[3]
            for row in rows
            for pivot in row["pivot_rows"]
        )
        == 51_222,
        "THM2904 parent/open-pivot universe changed",
    )
    return rows


def joined_parent_rows() -> list[dict[str, object]]:
    ledger_rows = parse_parent_ledger()
    source_rows = H.survivor_inputs()
    source_by_key = {branch_key(row): row for row in source_rows}
    require(
        len(source_by_key) == len(source_rows) == len(ledger_rows),
        "THM2904 source branch key collision",
    )
    seen: set[tuple[object, ...]] = set()
    for row in ledger_rows:
        key = branch_key(row)
        require(key in source_by_key and key not in seen, "source/ledger key join")
        seen.add(key)
        source = source_by_key[key]
        require(
            row["stratum"] == source["stratum"]
            and row["mass"] == source["mass"]
            and row["pair_cap"] == source["pair_cap"]
            and row["threshold"] == source["threshold"]
            and row["pair_exception"] == source["pair_exception"],
            "source/ledger mathematical fields changed",
        )
        row["components"] = source["components"]
        row["qs"] = source["qs"]
    require(seen == set(source_by_key), "source/ledger join is not surjective")
    ledger_rows.sort(key=lambda row: (row["body"], row["rank"], row["apex"]))
    return ledger_rows


def child_row(
    row: dict[str, object],
    index: int,
    carrier: list[tuple[F, F]],
) -> dict[str, object]:
    coverage, center = row["core_rows"][index]
    _, _, parent_margin, parent_closed, _ = row["pivot_rows"][index]
    require(not parent_closed, "closed THM2904 pivot entered child census")

    child = H.R.subtract_local(carrier, center)
    child_mass = sum((right - left for left, right in child), F(0))
    require(
        child_mass == row["mass"] - coverage and child_mass > 0,
        "literal child mass changed",
    )
    direct, direct_components, direct_mass = H.R.CORE.good_norm(
        tuple(sorted((*row["body"], row["apex"], center)))
    )
    require(
        child == direct
        and len(child) == direct_components
        and child_mass == direct_mass,
        "sequential child differs from direct full-family reconstruction",
    )

    earlier = tuple(label for _, label in row["core_rows"][:index])
    forbidden = frozenset((*row["prefix"], center, *earlier))
    require(
        len(forbidden) == len(row["prefix"]) + 1 + len(earlier)
        and center not in row["prefix"]
        and center not in earlier
        and not (set(earlier) & set(row["prefix"]))
        and min(forbidden) >= FIRST_EXTERNAL
        and max(forbidden) <= INITIAL_HEAD_LIMIT,
        "inherited prefix/centre/earlier-centre exclusion changed",
    )

    initial_labels = [
        label
        for label in range(FIRST_EXTERNAL, INITIAL_HEAD_LIMIT + 1)
        if label not in forbidden
    ]
    head = H.exact_coverages(child, initial_labels)
    require(len(head) >= 4, "initial child head has fewer than four labels")
    top_four_initial = tuple(sorted(head, reverse=True)[:4])
    require(
        len({label for _, label in top_four_initial}) == 4,
        "initial top four lost label distinctness",
    )

    gamma = H.S2 * len(child) / 7
    initial_excess = top_four_initial[-1][0] - child_mass / 7
    require(
        initial_excess > 0,
        "fourth head singleton does not give a finite THM735 horizon",
    )
    initial_required = max(
        INITIAL_HEAD_LIMIT,
        ceiling(gamma / initial_excess) - 1,
    )
    extra_labels = [
        label
        for label in range(INITIAL_HEAD_LIMIT + 1, initial_required + 1)
        if label not in forbidden
    ]
    if extra_labels:
        head.extend(H.exact_coverages(child, extra_labels))
    top_four = tuple(sorted(head, reverse=True)[:4])
    require(
        len({label for _, label in top_four}) == 4
        and all(label not in forbidden for _, label in top_four),
        "final top four violated distinctness or inherited exclusion",
    )
    final_excess = top_four[-1][0] - child_mass / 7
    require(final_excess > 0, "final fourth singleton lost finite excess")
    final_required = max(
        INITIAL_HEAD_LIMIT,
        ceiling(gamma / final_excess) - 1,
    )
    tail_bound = child_mass / 7 + gamma / (initial_required + 1)
    tail_gap = top_four[-1][0] - tail_bound
    require(
        final_required <= initial_required and tail_gap >= 0,
        "dynamic THM735 tail horizon failed",
    )
    top_four_sum = sum((value for value, _ in top_four), F(0))
    margin = child_mass - top_four_sum
    return {
        "body": row["body"],
        "rank": row["rank"],
        "apex": row["apex"],
        "prefix": row["prefix"],
        "pair_exception": row["pair_exception"],
        "parent_mass": row["mass"],
        "parent_components": row["components"],
        "center_index": index + 1,
        "center": center,
        "center_coverage": coverage,
        "earlier": earlier,
        "parent_margin": parent_margin,
        "child_mass": child_mass,
        "child_components": len(child),
        "initial_q4": top_four_initial[-1][0],
        "initial_required": initial_required,
        "head_limit": initial_required,
        "head_scanned": len(head),
        "final_required": final_required,
        "tail_bound": tail_bound,
        "tail_gap": tail_gap,
        "top_four": top_four,
        "margin": margin,
        "closed": margin > 0,
    }


def profile_parent(row: dict[str, object]) -> list[dict[str, object]]:
    carrier, components, mass = H.R.CORE.good_norm(
        tuple(sorted((*row["body"], row["apex"])))
    )
    require(
        components == row["components"]
        and mass == row["mass"]
        and mass > 0,
        "parent carrier reconstruction changed",
    )
    return [
        child_row(row, index, carrier)
        for index, pivot in enumerate(row["pivot_rows"])
        if not pivot[3]
    ]


def ledger_line(row: dict[str, object]) -> str:
    return (
        f"CHILD;E={','.join(map(str, row['body']))};rank={row['rank']};"
        f"a={row['apex']};P={','.join(map(str, row['prefix']))};"
        f"exception={int(row['pair_exception'])};"
        f"hC={ftext(row['parent_mass'])};rC={row['parent_components']};"
        f"i={row['center_index']};x={row['center']};"
        f"cx={ftext(row['center_coverage'])};"
        f"earlier={','.join(map(str, row['earlier']))};"
        f"pmargin={ftext(row['parent_margin'])};"
        f"hR={ftext(row['child_mass'])};rR={row['child_components']};"
        f"q4_2000={ftext(row['initial_q4'])};"
        f"N0={row['initial_required']};M={row['head_limit']};"
        f"scan={row['head_scanned']};N={row['final_required']};"
        f"sealed=1;tail={ftext(row['tail_bound'])};"
        f"tailgap={ftext(row['tail_gap'])};"
        + "top4="
        + ",".join(
            f"{label}:{ftext(value)}"
            for value, label in row["top_four"]
        )
        + f";margin={ftext(row['margin'])};closed={int(row['closed'])}"
    )


def proved_unions() -> tuple[
    set[tuple[int, ...]],
    set[tuple[int, ...]],
    set[tuple[int, ...]],
    set[tuple[int, ...]],
    set[tuple[int, ...]],
]:
    through_2905 = H.current_proved_union()
    require(len(through_2905) == 82, "THM2905 union changed")
    through_2904 = through_2905 | set(H.EXPECTED_ADDITIVE_ROOTS)
    require(len(through_2904) == 88, "THM2904 union changed")
    additive_2911 = {
        tuple_ints(line.removeprefix("ADDITIVE_ROOT="))
        for line in THM2911_OUT.read_text().splitlines()
        if line.startswith("ADDITIVE_ROOT=")
    }
    require(
        len(additive_2911) == 93
        and not (through_2904 & additive_2911),
        "THM2911 additive root bank changed",
    )
    through_2911 = through_2904 | additive_2911
    require(len(through_2911) == 181, "THM2911 proved union changed")
    roots_2912 = set(literal_output_value(THM2912_OUT, "closed_roots="))
    require(
        len(roots_2912) == 172
        and len(through_2911 & roots_2912) == 39,
        "THM2912 root bank/overlap changed",
    )
    through_2912 = through_2911 | roots_2912
    require(len(through_2912) == 314, "canonical union through THM2912 changed")
    roots_2913 = set(literal_output_value(THM2913_OUT, "closed_roots="))
    additive_2913 = set(literal_output_value(THM2913_OUT, "additive_roots="))
    require(
        len(roots_2913) == 38
        and len(additive_2913) == 37
        and roots_2913 - through_2912 == additive_2913,
        "THM2913 root bank/set difference changed",
    )
    through_2913 = through_2912 | roots_2913
    require(len(through_2913) == 351, "canonical union through THM2913 changed")
    return through_2911, roots_2912, through_2912, roots_2913, through_2913


def finite_h1_keys(
    rows: list[dict[str, object]],
) -> set[tuple[object, ...]]:
    keys: set[tuple[object, ...]] = set()
    for row in rows:
        rank_four = sum(row["qs"][:4], F(0))
        epsilon = 6 * row["mass"] / 7 - rank_four
        if epsilon <= 0:
            continue
        cutoff = ceiling(
            H.S2 * row["components"] / (7 * epsilon)
        ) - 1
        if cutoff <= THM2911_MAX_CUTOFF:
            keys.add(branch_key(row))
    require(len(keys) == 3_090, "THM2911 finite-H1 survivor slice changed")
    return keys


def root_projection(
    rows_by_body: dict[tuple[int, ...], set[tuple[object, ...]]],
    closed_keys: set[tuple[object, ...]],
) -> set[tuple[int, ...]]:
    return {
        body
        for body, keys in rows_by_body.items()
        if keys <= closed_keys
    }


def write_ledger(
    path: Path,
    children: list[dict[str, object]],
) -> tuple[str, str]:
    semantic = hashlib.sha256(
        b"LRC14/j6/all-open-centres/child-top4/THM2915/v1\n"
    )
    lines = ["LRC14 j6 all-open-centre exact child-top-four ledger"]
    for row in children:
        line = ledger_line(row)
        semantic.update((line + "\n").encode())
        lines.append(line)
    semantic_sha256 = semantic.hexdigest()
    lines.extend(
        (
            f"ledger_semantic_sha256={semantic_sha256}",
            "scope=all 51222 open THM2904 centres;exact dynamic-tail child top4;not LRC14",
        )
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write("\n".join(lines) + "\n")
    return semantic_sha256, file_sha256(path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--workers",
        type=int,
        default=min(os.cpu_count() or 1, 8),
    )
    parser.add_argument("--ledger", type=Path, required=True)
    args = parser.parse_args()
    require(args.workers >= 1, "worker count must be positive")
    check_dependencies()

    parents = joined_parent_rows()
    open_parents = [row for row in parents if not row["branch_closed"]]
    parent_closed = [row for row in parents if row["branch_closed"]]
    require(
        len(open_parents) == 11_563
        and len(parent_closed) == 279,
        "THM2904 open/closed parent partition changed",
    )

    if args.workers == 1:
        nested = list(map(profile_parent, open_parents))
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            nested = pool.map(profile_parent, open_parents, chunksize=1)
    children = [row for group in nested for row in group]
    children.sort(
        key=lambda row: (
            row["body"],
            row["rank"],
            row["apex"],
            row["prefix"],
            row["center_index"],
        )
    )
    require(len(children) == 51_222, "child universe changed")

    expected_child_keys: set[tuple[object, ...]] = set()
    for row in open_parents:
        for index, pivot in enumerate(row["pivot_rows"]):
            if pivot[3]:
                continue
            key = (
                row["body"],
                row["rank"],
                row["apex"],
                row["prefix"],
                pivot[0],
                tuple(label for _, label in row["core_rows"][:index]),
            )
            require(key not in expected_child_keys, "open-pivot identity collision")
            expected_child_keys.add(key)
    actual_child_keys = {child_key(row) for row in children}
    require(
        len(actual_child_keys) == len(children)
        and actual_child_keys == expected_child_keys,
        "51,222-child identity join failed",
    )
    short_keys = {
        (row["body"], row["rank"], row["apex"], row["center"])
        for row in children
    }
    require(
        len(short_keys) == len(children),
        "short child identity unexpectedly collided",
    )
    open_key_digest = key_digest(
        actual_child_keys,
        b"THM2915/open-pivot-keys/v1",
    )
    if EXPECTED_OPEN_PIVOT_KEY_DIGEST is not None:
        require(
            open_key_digest == EXPECTED_OPEN_PIVOT_KEY_DIGEST,
            "open-pivot identity digest changed",
        )

    children_by_parent: dict[
        tuple[object, ...], list[dict[str, object]]
    ] = defaultdict(list)
    children_by_body: dict[
        tuple[int, ...], list[dict[str, object]]
    ] = defaultdict(list)
    for row in children:
        children_by_parent[branch_key(row)].append(row)
        children_by_body[row["body"]].append(row)
    require(
        set(children_by_parent) == {branch_key(row) for row in open_parents},
        "child/open-parent identity join failed",
    )

    T = {
        key
        for key, group in children_by_parent.items()
        if all(row["closed"] for row in group)
    }
    P = {branch_key(row) for row in parent_closed}
    E = {branch_key(row) for row in parents if row["pair_exception"]}
    H1 = finite_h1_keys(parents)
    require(
        len(T) == 8_112
        and len(P) == 279
        and len(E) == 52
        and not (P & T)
        and not (E & (H1 | P))
        and not (E & T),
        "E/H/P/T branch routes changed",
    )

    atoms = Counter()
    for row in parents:
        key = branch_key(row)
        atom = "".join(
            label
            for label, route in (("E", E), ("H", H1), ("P", P), ("T", T))
            if key in route
        ) or "-"
        atoms[atom] += 1
    branch_atoms = tuple(sorted(atoms.items()))

    rows_by_body: dict[
        tuple[int, ...], set[tuple[object, ...]]
    ] = defaultdict(set)
    for row in parents:
        rows_by_body[row["body"]].add(branch_key(row))
    all_bodies = set(combinations(range(1, 15), 7))
    survivor_bodies = set(rows_by_body)
    no_g5_bodies = all_bodies - survivor_bodies
    P_roots = root_projection(rows_by_body, P)
    pure_roots = root_projection(rows_by_body, P | T)
    e_roots = root_projection(rows_by_body, E | P | T)
    h_roots = root_projection(rows_by_body, H1 | P | T)
    composed_roots = root_projection(rows_by_body, E | H1 | P | T)
    live_route_bodies = set(children_by_body)
    live_route_roots = {
        body
        for body, group in children_by_body.items()
        if all(row["closed"] for row in group)
    }
    require(
        len(all_bodies) == OFFICIAL_ROOT_COUNT
        and len(survivor_bodies) == 3_411
        and len(no_g5_bodies) == 21
        and len(P_roots) == 10
        and len(live_route_bodies) == 3_401
        and len(live_route_roots) == 1_219
        and pure_roots == P_roots | live_route_roots
        and len(pure_roots) == 1_229,
        "whole-root route quantifiers changed",
    )

    (
        through_2911,
        roots_2912,
        baseline314,
        roots_2913,
        baseline351,
    ) = proved_unions()
    require(
        no_g5_bodies <= baseline314 and P_roots <= baseline314,
        "outside-live-route roots escaped the canonical baseline",
    )
    pure_additive = pure_roots - baseline314
    pure_union314 = pure_roots | baseline314
    composed_additive = composed_roots - baseline314
    composed_union314 = composed_roots | baseline314
    e_additive_gain = (e_roots - pure_roots) - baseline314
    h_additive_gain = (h_roots - pure_roots) - baseline314
    mixed_gain = composed_roots - (e_roots | h_roots)
    require(
        len(pure_additive) == 952
        and len(pure_union314) == 1_266
        and len(composed_roots) == 1_248
        and len(composed_roots & baseline314) == 284
        and len(composed_additive) == 964
        and len(composed_union314) == 1_278
        and len(e_additive_gain) == 6
        and len(h_additive_gain) == 6
        and not mixed_gain
        and not (e_additive_gain & h_additive_gain),
        "strong branchwise root composition changed",
    )

    roots_2912_failed = set(literal_output_value(THM2912_OUT, "failed_roots="))
    one_h3_bodies = roots_2912 | roots_2912_failed
    one_h3_parent_keys: set[tuple[object, ...]] = set()
    for body in one_h3_bodies:
        candidates = [
            branch_key(row)
            for row in open_parents
            if row["body"] == body and not row["pair_exception"]
        ]
        require(
            len(candidates) == 1,
            f"{body}: one-H3 ordinary parent identity changed",
        )
        one_h3_parent_keys.add(candidates[0])
    one_h3_children = [
        row for row in children
        if branch_key(row) in one_h3_parent_keys
    ]
    one_h3_slice_closed = {
        body
        for body in one_h3_bodies
        if all(
            row["closed"]
            for row in one_h3_children
            if row["body"] == body
        )
    }
    one_h3_slice_only = roots_2912 - live_route_roots
    require(
        len(one_h3_parent_keys) == 210
        and len(one_h3_children) == 807
        and sum(row["closed"] for row in one_h3_children) == 765
        and sum(not row["closed"] for row in one_h3_children) == 42
        and one_h3_slice_closed == roots_2912
        and len(one_h3_slice_only) == 2
        and one_h3_slice_only
        == {
            (2, 3, 6, 10, 12, 13, 14),
            (2, 6, 7, 8, 9, 10, 11),
        }
        and roots_2912 <= composed_roots,
        "THM2912 one-H3 slice reproduction changed",
    )

    Q = {
        key
        for key in one_h3_parent_keys
        if key[0] in roots_2913
    }
    require(
        len(Q) == 38
        and len(Q & H1) == 1
        and not (Q & (E | P | T))
        and roots_2913 == roots_2912_failed,
        "THM2913 branch-key embedding changed",
    )
    atoms_with_q = Counter()
    for row in parents:
        key = branch_key(row)
        atom = "".join(
            label
            for label, route in (
                ("E", E),
                ("H", H1),
                ("P", P),
                ("Q", Q),
                ("T", T),
            )
            if key in route
        ) or "-"
        atoms_with_q[atom] += 1
    branch_atoms_with_q = tuple(sorted(atoms_with_q.items()))

    final_roots = root_projection(rows_by_body, E | H1 | P | Q | T)
    q_root_gain = final_roots - composed_roots
    current_additive = final_roots - baseline351
    current_overlap = final_roots & baseline351
    current_union = final_roots | baseline351
    require(
        len(final_roots) == 1_285
        and q_root_gain == roots_2913 - baseline314
        and len(q_root_gain) == 37
        and len(current_overlap) == 321
        and current_additive == composed_additive
        and len(current_additive) == 964
        and len(current_union) == 1_315,
        "current THM2913-aware branch composition changed",
    )

    exception_children = [
        row for row in children if branch_key(row) in E
    ]
    exception_fail_hist = Counter(
        sum(not child["closed"] for child in children_by_parent[key])
        for key in E
    )
    closed_children = [row for row in children if row["closed"]]
    failed_children = [row for row in children if not row["closed"]]
    closest_positive = min(
        closed_children,
        key=lambda row: (
            row["margin"],
            row["body"],
            row["rank"],
            row["apex"],
            row["center"],
        ),
    )
    closest_failure = max(
        failed_children,
        key=lambda row: (
            row["margin"],
            tuple(-value for value in row["body"]),
            -row["rank"],
            -row["apex"],
            -row["center"],
        ),
    )
    failed_children_per_failed_root = Counter(
        sum(not row["closed"] for row in group)
        for body, group in children_by_body.items()
        if body not in live_route_roots
    )

    summary: dict[str, object] = {
        "parent_rows": len(parents),
        "pair_exception_rows": len(E),
        "hostile_centres": sum(len(row["core_rows"]) for row in parents),
        "THM2904_closed_pivots": sum(
            pivot[3] for row in parents for pivot in row["pivot_rows"]
        ),
        "open_pivots": len(children),
        "parent_closed_rows": len(P),
        "open_parent_rows": len(open_parents),
        "closed_children": len(closed_children),
        "failed_children": len(failed_children),
        "equality_children": sum(row["margin"] == 0 for row in children),
        "sealed_children": sum(row["tail_gap"] >= 0 for row in children),
        "extended_children": sum(
            row["head_limit"] > INITIAL_HEAD_LIMIT for row in children
        ),
        "maximum_head_limit": max(row["head_limit"] for row in children),
        "maximum_final_required": max(row["final_required"] for row in children),
        "T_rows": len(T),
        "remaining_open_rows": len(open_parents) - len(T),
        "H_rows": len(H1),
        "E_rows": len(E),
        "P_rows": len(P),
        "Q_rows": len(Q),
        "branch_atoms": branch_atoms,
        "branch_atoms_with_Q": branch_atoms_with_q,
        "H_intersection_T": len(H1 & T),
        "H_intersection_P": len(H1 & P),
        "E_intersection_T": len(E & T),
        "Q_intersection_H": len(Q & H1),
        "Q_intersection_T": len(Q & T),
        "survivor_bodies": len(survivor_bodies),
        "no_G5_bodies": len(no_g5_bodies),
        "P_roots": len(P_roots),
        "live_route_bodies": len(live_route_bodies),
        "live_route_roots": len(live_route_roots),
        "live_route_failures": len(live_route_bodies - live_route_roots),
        "pure_roots": len(pure_roots),
        "pure_overlap314": len(pure_roots & baseline314),
        "pure_additive314": len(pure_additive),
        "pure_union314": len(pure_union314),
        "pure_residual314": OFFICIAL_ROOT_COUNT - len(pure_union314),
        "EHP_composed_roots": len(composed_roots),
        "EHP_overlap314": len(composed_roots & baseline314),
        "EHP_additive314": len(composed_additive),
        "EHP_union314": len(composed_union314),
        "EHP_residual314": OFFICIAL_ROOT_COUNT - len(composed_union314),
        "E_additive_gain": len(e_additive_gain),
        "H_additive_gain": len(h_additive_gain),
        "mixed_EH_gain": len(mixed_gain),
        "EHPQ_composed_roots": len(final_roots),
        "Q_root_gain": len(q_root_gain),
        "Q_extra_synergy": len(q_root_gain - (roots_2913 - baseline314)),
        "current_baseline": len(baseline351),
        "current_overlap": len(current_overlap),
        "current_additive": len(current_additive),
        "current_union": len(current_union),
        "current_residual": OFFICIAL_ROOT_COUNT - len(current_union),
        "exception_children": len(exception_children),
        "exception_closed_children": sum(
            row["closed"] for row in exception_children
        ),
        "exception_failed_children": sum(
            not row["closed"] for row in exception_children
        ),
        "exception_completed_rows": len(E & T),
        "exception_fail_histogram": tuple(sorted(exception_fail_hist.items())),
        "one_H3_children": len(one_h3_children),
        "one_H3_closed_children": sum(
            row["closed"] for row in one_h3_children
        ),
        "one_H3_failed_children": sum(
            not row["closed"] for row in one_h3_children
        ),
        "one_H3_slice_roots": len(one_h3_slice_closed),
        "one_H3_full_route_roots": len(live_route_roots & roots_2912),
        "one_H3_composed_roots": len(composed_roots & roots_2912),
        "child_component_histogram": tuple(
            sorted(Counter(row["child_components"] for row in children).items())
        ),
        "open_centres_per_parent": tuple(
            sorted(Counter(len(group) for group in children_by_parent.values()).items())
        ),
        "failed_children_per_failed_root": tuple(
            sorted(failed_children_per_failed_root.items())
        ),
    }
    if EXPECTED_SUMMARY is not None:
        require(summary == EXPECTED_SUMMARY, "locked THM2915 summary changed")

    root_sets = {
        "live_route1219": live_route_roots,
        "live_route_overlap267": live_route_roots & baseline314,
        "pure_additive952": pure_additive,
        "pure_union1266": pure_union314,
        "composed_roots1248": composed_roots,
        "composed_overlap284": composed_roots & baseline314,
        "composed_additive964": composed_additive,
        "composed_union1278": composed_union314,
        "E_additive_gain6": e_additive_gain,
        "H_additive_gain6": h_additive_gain,
        "Q_roots38": roots_2913,
        "final_roots1285": final_roots,
        "current_overlap321": current_overlap,
        "current_additive964": current_additive,
        "current_union1315": current_union,
    }
    root_digests = {
        name: body_digest(values) for name, values in root_sets.items()
    }
    if EXPECTED_ROOT_DIGESTS is not None:
        require(
            root_digests == EXPECTED_ROOT_DIGESTS,
            "locked THM2915 root sets changed",
        )

    ledger_semantic, ledger_file = write_ledger(args.ledger, children)
    if EXPECTED_LEDGER_SEMANTIC_SHA256 is not None:
        require(
            ledger_semantic == EXPECTED_LEDGER_SEMANTIC_SHA256,
            "child ledger semantic digest changed",
        )
    if EXPECTED_LEDGER_FILE_SHA256 is not None:
        require(
            ledger_file == EXPECTED_LEDGER_FILE_SHA256,
            "child ledger file hash changed",
        )

    print("LRC14 all-open-centre exact child-top-four closure")
    print(f"initial_head_limit={INITIAL_HEAD_LIMIT}")
    print(f"summary={summary}")
    print(f"root_digests={root_digests}")
    print(f"open_pivot_key_digest={open_key_digest}")
    print(
        "closest_positive="
        f"{closest_positive['body']};rank={closest_positive['rank']};"
        f"a={closest_positive['apex']};x={closest_positive['center']};"
        f"margin={ftext(closest_positive['margin'])}"
    )
    print(
        "closest_failure="
        f"{closest_failure['body']};rank={closest_failure['rank']};"
        f"a={closest_failure['apex']};x={closest_failure['center']};"
        f"margin={ftext(closest_failure['margin'])}"
    )
    print(
        f"minimum_tail_gap={ftext(min(row['tail_gap'] for row in children))}"
    )
    print(f"one_H3_slice_only={tuple(sorted(one_h3_slice_only))}")
    print(f"E_additive_gain_roots={tuple(sorted(e_additive_gain))}")
    print(f"H_additive_gain_roots={tuple(sorted(h_additive_gain))}")
    print(f"Q_root_gain_roots={tuple(sorted(q_root_gain))}")
    print(f"live_route_roots={tuple(sorted(live_route_roots))}")
    print(f"pure_additive_roots={tuple(sorted(pure_additive))}")
    print(f"composed_roots={tuple(sorted(composed_roots))}")
    print(f"composed_additive_roots={tuple(sorted(composed_additive))}")
    print(f"final_roots={tuple(sorted(final_roots))}")
    print(f"current_union_roots={tuple(sorted(current_union))}")
    print(f"ledger_semantic_sha256={ledger_semantic}")
    print(f"ledger_file_sha256={ledger_file}")
    print(
        "mode=DISCOVERY"
        if any(
            value is None
            for value in (
                EXPECTED_SUMMARY,
                EXPECTED_ROOT_DIGESTS,
                EXPECTED_LEDGER_SEMANTIC_SHA256,
                EXPECTED_LEDGER_FILE_SHA256,
            )
        )
        else "mode=LOCKED"
    )
    print(
        "scope=all 51222 open THM2904 centre children;"
        "dynamic exact top4 plus E/H/P/Q/T branch composition;not LRC14"
    )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
