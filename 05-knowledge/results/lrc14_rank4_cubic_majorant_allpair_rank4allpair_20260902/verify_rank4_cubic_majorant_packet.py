#!/usr/bin/env python3
"""Closed verifier for the all-residual-pair rank-four cubic packet.

The verifier does not reconstruct walls or repeat the optimizer.  It pins the
inherited universe and rank-three ledger, audits every row and exact integer
identity, recomputes both the economical and exact-all aggregate verdicts,
and derives the pairwise/uniform appender cutoffs using Python big integers.
Literal-wall and optimizer independence are supplied by the separate hostile
audit described in this packet's README.
"""

from __future__ import annotations

import csv
import hashlib
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
PACKET = Path(__file__).resolve().parent
UNIVERSE = (
    ROOT
    / "05-knowledge/results/"
    "lrc14_rank2_wall_graph_complete_pair_closure_thm4326/"
    "thm4231_remainder181194.csv"
)
RANK3 = (
    ROOT
    / "05-knowledge/results/lrc14_rank3_surplus_third_tail_thm4333/"
    "pair_rank3_2over27_screen_all.csv"
)
LEDGER = PACKET / "rank4_cubic_majorant_exactall.csv"
SUMMARY = PACKET / "rank4_cubic_majorant_exactall.out"

UNIVERSE_SHA256 = (
    "9dfbf0a8948bf23016ae40f40f9118020c9429ac60421681a9286fc4d34041a1"
)
RANK3_SHA256 = (
    "5e2d740e4f5d8be925be7d244661c8148fc045362e9165b1c9c1a2f5a35566e5"
)
LEDGER_SHA256 = (
    "60dab8a471065dee132e61a5695e6a827616e0e45e7c6a67cd11b426fd86623a"
)
POOL = (
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290,
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def labels(mask: int) -> str:
    return ";".join(str(POOL[i]) for i in range(30) if mask >> i & 1)


def strict_appender_cutoff(q: int, r: int, grid: int, lower14: int) -> int:
    # If m=lower14/(14D), interval discrepancy requires
    # s > (6 C/49)/((6/7)m-4/63)
    #   = 54 C D/(27 lower14-28D), C=1891+q+r.
    component_bound = 1891 + q + r
    numerator = 54 * component_bound * grid
    denominator = 27 * lower14 - 28 * grid
    require(denominator > 0, f"{q},{r}: nonpositive appender denominator")
    cutoff = numerator // denominator + 1
    m = Fraction(lower14, 14 * grid)
    lower_at_cutoff = (
        Fraction(6, 7) * m
        - Fraction(6 * component_bound, 49 * cutoff)
    )
    require(
        lower_at_cutoff > Fraction(4, 63),
        f"{q},{r}: appender cutoff is not strict",
    )
    if cutoff > 1:
        lower_before = (
            Fraction(6, 7) * m
            - Fraction(6 * component_bound, 49 * (cutoff - 1))
        )
        require(
            lower_before <= Fraction(4, 63),
            f"{q},{r}: appender cutoff is not least for certificate",
        )
    return cutoff


def main() -> None:
    require(sha256(UNIVERSE) == UNIVERSE_SHA256, "universe SHA256 changed")
    require(sha256(RANK3) == RANK3_SHA256, "rank-three ledger SHA256 changed")
    require(sha256(LEDGER) == LEDGER_SHA256, "rank-four ledger SHA256 changed")

    with UNIVERSE.open(newline="") as handle:
        universe = [tuple(map(int, row)) for row in csv.reader(handle)]
    require(len(universe) == 181_194, "universe cardinality changed")
    require(len(set(universe)) == len(universe), "universe has duplicates")
    require(all(0 < q < r for q, r in universe), "universe ordering changed")

    with RANK3.open(newline="") as handle:
        rank3_rows = list(csv.DictReader(handle))
    require(len(rank3_rows) == len(universe), "rank-three row count changed")

    with LEDGER.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    require(len(rows) == len(universe), "rank-four row count changed")

    degree_positive = 0
    exact_nodes = 0
    exact_prunes = 0
    economical_nodes = 0
    economical_prunes = 0
    exact_worst: tuple[Fraction, int, int, dict[str, str]] | None = None
    economical_worst: tuple[Fraction, int, int, dict[str, str]] | None = None
    exact_cutoff_worst: tuple[int, int, int] | None = None
    economical_cutoff_worst: tuple[int, int, int] | None = None
    max_nodes: tuple[int, int, int] | None = None
    controls: dict[tuple[int, int], dict[str, str]] = {}

    for index, (row, inherited, pair) in enumerate(
        zip(rows, rank3_rows, universe, strict=True)
    ):
        q, r = pair
        require(int(row["index"]) == index, f"row {index}: index changed")
        require((int(row["q"]), int(row["r"])) == pair, f"row {index}: pair")
        require(
            (int(inherited["q"]), int(inherited["r"])) == pair,
            f"row {index}: rank-three pair",
        )
        grid = int(row["grid"])
        require(grid == int(inherited["grid"]), f"{q},{r}: grid mismatch")
        masses = [int(row[f"rank{k}"]) for k in range(5)]
        require(all(value >= 0 for value in masses), f"{q},{r}: negative mass")
        require(
            sum(masses[:4]) == int(inherited["rank3_total"]),
            f"{q},{r}: THM-4333 rank-three prefix mismatch",
        )
        retained = int(row["retained"])
        require(sum(masses) == retained, f"{q},{r}: retained partition")

        coarse_lower14 = int(row["coarse_lower14"])
        coarse_ticks = int(row["coarse_ticks"])
        require(
            coarse_ticks == 81 * coarse_lower14 - 98 * grid,
            f"{q},{r}: coarse tick identity",
        )
        degree_positive += coarse_ticks > 0

        require(
            row["certificate_type"] == "CUBIC_EXACT",
            f"{q},{r}: exact-all certificate type changed",
        )
        lower14 = int(row["certificate_lower14"])
        ticks = int(row["certificate_ticks"])
        require(
            ticks == 81 * lower14 - 98 * grid,
            f"{q},{r}: cubic tick identity",
        )
        require(ticks > 0, f"{q},{r}: nonpositive cubic certificate")
        mask = int(row["certificate_body"], 16)
        require(mask.bit_count() == 8, f"{q},{r}: body rank changed")
        require(row["certificate_labels"] == labels(mask), f"{q},{r}: labels")
        direct_l4 = int(row["direct_l4"])
        require(
            14 * direct_l4 >= lower14,
            f"{q},{r}: direct optimizer-body control below certificate",
        )
        nodes = int(row["nodes"])
        prunes = int(row["prunes"])
        require(nodes > 0 and 0 <= prunes < nodes, f"{q},{r}: search totals")
        exact_nodes += nodes
        exact_prunes += prunes
        if coarse_ticks <= 0:
            economical_nodes += nodes
            economical_prunes += prunes

        exact_ratio = Fraction(lower14, 14 * grid)
        require(exact_ratio > Fraction(7, 81), f"{q},{r}: exact ratio target")
        exact_key = (exact_ratio, q, r, row)
        if exact_worst is None or exact_key[:3] < exact_worst[:3]:
            exact_worst = exact_key

        economical_lower14 = lower14 if coarse_ticks <= 0 else coarse_lower14
        economical_ratio = Fraction(economical_lower14, 14 * grid)
        require(
            economical_ratio > Fraction(7, 81),
            f"{q},{r}: economical ratio target",
        )
        economical_key = (economical_ratio, q, r, row)
        if economical_worst is None or economical_key[:3] < economical_worst[:3]:
            economical_worst = economical_key

        exact_cutoff = strict_appender_cutoff(q, r, grid, lower14)
        cutoff_key = (exact_cutoff, q, r)
        if exact_cutoff_worst is None or cutoff_key > exact_cutoff_worst:
            exact_cutoff_worst = cutoff_key
        economical_cutoff = strict_appender_cutoff(
            q, r, grid, economical_lower14
        )
        economical_cutoff_key = (economical_cutoff, q, r)
        if (
            economical_cutoff_worst is None
            or economical_cutoff_key > economical_cutoff_worst
        ):
            economical_cutoff_worst = economical_cutoff_key

        node_key = (nodes, q, r)
        if max_nodes is None or node_key > max_nodes:
            max_nodes = node_key
        if pair in {(50, 70), (509, 640), (50, 140)}:
            controls[pair] = row

    require(degree_positive == 29_595, "degree-positive split changed")
    require(exact_nodes == 344_646_646, "exact-all node total changed")
    require(exact_prunes == 171_098_718, "exact-all prune total changed")
    require(economical_nodes == 270_641_887, "economical node total changed")
    require(economical_prunes == 134_308_799, "economical prune total changed")
    require(max_nodes == (8_625, 78, 218), "max-node hostile changed")
    require(exact_worst is not None and economical_worst is not None, "no rows")
    require(exact_cutoff_worst is not None, "no exact cutoff")
    require(economical_cutoff_worst is not None, "no economical cutoff")

    exact_ratio, exact_q, exact_r, exact_row = exact_worst
    require((exact_q, exact_r) == (50, 70), "exact normalized hostile changed")
    require(
        int(exact_row["certificate_lower14"]) == 121_289_023_967_200,
        "exact normalized hostile lower bound changed",
    )
    require(
        exact_ratio == Fraction(5_227_975_171, 55_037_980_998),
        "exact normalized hostile ratio changed",
    )
    economical_ratio, economical_q, economical_r, _ = economical_worst
    require(
        (economical_q, economical_r) == (499, 636),
        "economical normalized hostile changed",
    )
    require(
        exact_cutoff_worst == (13_737, 50, 70),
        "cubic-only uniform appender cutoff changed",
    )
    require(
        economical_cutoff_worst == (38_442, 721, 746),
        "economical uniform appender cutoff changed",
    )

    # THM-4336 independently proves sharper exact L4 bounds on (50,70) and
    # (509,640), whose appender cutoffs are 6021 and 5295.  Replacing those
    # two cubic rows makes the maximum of all remaining pairwise cubic
    # cutoffs occur at (50,140).
    thm4336_cutoffs = {(50, 70): 6_021, (509, 640): 5_295}
    remaining_cutoffs = []
    for row in rows:
        q, r = int(row["q"]), int(row["r"])
        if (q, r) in thm4336_cutoffs:
            remaining_cutoffs.append((thm4336_cutoffs[(q, r)], q, r))
        elif (q, r) == (509, 640):
            remaining_cutoffs.append((5_295, q, r))
        else:
            remaining_cutoffs.append(
                (
                    strict_appender_cutoff(
                        q,
                        r,
                        int(row["grid"]),
                        int(row["certificate_lower14"]),
                    ),
                    q,
                    r,
                )
            )
    require(
        max(remaining_cutoffs) == (12_274, 50, 140),
        "THM-4336 hybrid uniform cutoff changed",
    )

    summary = SUMMARY.read_text(encoding="ascii")
    for needle in (
        "INPUT_PAIRS 181194 INPUT_NORMALIZED_FNV 95425eabee50378c",
        "SLICE_START 0 SLICE_COUNT 181194 EXACT_ALL 1",
        "DEGREE_POSITIVE 29595 DEGREE_CERTIFICATES 0 "
        "CUBIC_EXACT_POSITIVE 181194 CUBIC_EXACT_NONPOSITIVE 0",
        "CUBIC_NODES 344646646 CUBIC_PRUNES 171098718",
        "VERDICT PASS",
    ):
        require(needle in summary, f"summary missing {needle!r}")

    print("LRC14_RANK4_CUBIC_MAJORANT_PACKET_VERIFY_V1")
    print(
        f"UNIVERSE_ROWS {len(universe)} SHA256 {UNIVERSE_SHA256} "
        "ORDER_UNIQUE PASS"
    )
    print(f"RANK3_LEDGER_SHA256 {RANK3_SHA256} ALL_PREFIXES_MATCH PASS")
    print(f"RANK4_LEDGER_SHA256 {LEDGER_SHA256} ALL_IDENTITIES_MATCH PASS")
    print(
        "ECONOMICAL DEGREE_POSITIVE 29595 CUBIC_FALLBACK 151599 "
        f"NODES {economical_nodes} PRUNES {economical_prunes} ALL_PASS"
    )
    print(
        "EXACT_ALL CUBIC_POSITIVE 181194 NODES "
        f"{exact_nodes} PRUNES {exact_prunes} ALL_PASS"
    )
    print(
        f"EXACT_NORMALIZED_MIN PAIR {exact_q},{exact_r} RATIO {exact_ratio} "
        f"SURPLUS {exact_ratio-Fraction(7,81)}"
    )
    print(
        "ECONOMICAL_NORMALIZED_MIN PAIR "
        f"{economical_q},{economical_r} RATIO {economical_ratio}"
    )
    print(
        "APPENDER_CUTOFF CUBIC_ONLY 13737 PAIR 50,70 "
        "THM4336_HYBRID 12274 PAIR 50,140 "
        "ECONOMICAL 38442 PAIR 721,746"
    )
    print(
        "DIRECT_L4_SCOPE POSITIVE_CONTROL_AT_REPORTED_MAJORANT_MAXIMIZER_"
        "NOT_GLOBAL_L4_MINIMUM"
    )
    print("SCOPE FIXED_POOL_THM4231_REMAINDER_NO_ARBITRARY_ENTRY_NO_LRC14")
    print("VERDICT PASS")


if __name__ == "__main__":
    main()
