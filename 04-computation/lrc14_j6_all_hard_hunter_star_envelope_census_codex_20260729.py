#!/usr/bin/env python3
r"""Exact Hunter-star postprocessing of all 14,806 scalar-hard LRC14 rows.

THM-2897 proves that, on a literal carrier C with singleton ranks

    q_1 >= ... >= q_5

and an attained global pair-union cap B_2, every five-cover is bounded by

    G_5 = max_{0 <= a <= min(q_1,B_2)}
          a + sum_{r=2}^5 min(a,q_r,B_2-a).

The objective is piecewise linear. Its maximum is therefore attained at an
endpoint or at a=q_r, a=B_2-q_r, or a=B_2/2. This verifier joins the
hash-pinned THM-2899 hard ledger to the hash-pinned THM-2901 pair ledger,
evaluates all clipped breakpoints exactly, and recomposes whole roots against
the proved 76-root union through THM-2903. This is a scoped census, not LRC(14).
"""

from __future__ import annotations

import hashlib
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
HARD_LEDGER = (
    ROOT
    / "05-knowledge/results/lrc14_j6_all_root_ranked_suffix_scalar_hard_ledger_codex_20260729.out"
)
PAIR_LEDGER = (
    ROOT
    / "05-knowledge/results/lrc14_j6_all_hard_global_pair_cap_census_codex_20260729.ledger.out"
)
HARD_LEDGER_SHA256 = (
    "6be9a6c9218f3b42b2eea733c9050f5d35160664af0f19390337b3c5be57cb37"
)
PAIR_LEDGER_SHA256 = (
    "5dea0eaa45dd52fbf1bef7cfcc328899a4789bc277b6e1e8ac2f4bdf192b85e4"
)

THM2895_ROOTS = {
    (2, 8, 9, 10, 11, 13, 14),
    (1, 3, 9, 10, 11, 12, 14),
    (2, 5, 9, 11, 12, 13, 14),
    (2, 3, 4, 5, 6, 7, 8),
}
THM2898_ROOTS = {(1, 8, 10, 11, 12, 13, 14)}
THM2899_ROOTS = {
    (1, 2, 3, 4, 5, 6, 13),
    (1, 2, 3, 4, 6, 7, 14),
    (1, 2, 3, 4, 6, 11, 13),
    (1, 2, 3, 4, 6, 12, 13),
    (7, 8, 9, 11, 12, 13, 14),
}
THM2901_ROOTS = {
    (1, 2, 3, 4, 6, 11, 12),
    (1, 2, 3, 5, 6, 10, 13),
    (1, 2, 4, 5, 6, 12, 13),
    (1, 3, 4, 5, 6, 7, 14),
    (5, 7, 8, 10, 11, 13, 14),
}
THM2902_ROOTS = {
    (1, 2, 3, 4, 5, 10, 12),
    (1, 2, 3, 4, 5, 12, 13),
}

EXPECTED_COUNTS: tuple[object, ...] = (
    14_806,
    2_964,
    1_835,
    1_129,
    0,
    52,
    0,
    3_427,
    16,
    5,
    11,
    61,
    5,
    76,
    6,
    82,
    3_350,
    (("low", 3_053, 609), ("one", 7_853, 1_571), ("both", 3_900, 784)),
    ((2, 5), (3, 1)),
)
EXPECTED_STAR_ROOTS = (
    (1, 2, 3, 4, 5, 12, 13),
    (1, 2, 3, 4, 6, 7, 13),
    (1, 2, 3, 4, 6, 8, 13),
    (1, 2, 3, 4, 6, 11, 12),
    (1, 2, 3, 4, 8, 9, 11),
    (1, 2, 3, 5, 6, 10, 13),
    (1, 2, 4, 5, 6, 10, 13),
    (1, 2, 4, 5, 6, 12, 13),
    (1, 3, 4, 5, 6, 7, 14),
    (4, 5, 6, 8, 10, 11, 12),
    (5, 6, 7, 9, 10, 11, 14),
    (5, 6, 8, 9, 10, 11, 12),
    (5, 6, 8, 9, 10, 12, 13),
    (5, 6, 8, 10, 11, 12, 13),
    (5, 7, 8, 10, 11, 13, 14),
    (6, 7, 9, 11, 12, 13, 14),
)
EXPECTED_ADDITIVE_ROOTS = (
    (1, 2, 3, 4, 8, 9, 11),
    (1, 2, 4, 5, 6, 10, 13),
    (5, 6, 7, 9, 10, 11, 14),
    (5, 6, 8, 9, 10, 12, 13),
    (5, 6, 8, 10, 11, 12, 13),
    (6, 7, 9, 11, 12, 13, 14),
)
EXPECTED_ROW_DIGEST = (
    "c1f60bdbc3669202dde60d8f44c9db887807b179cb42e02137d475c0f282d066"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_fraction(text: str) -> F:
    numerator, denominator = text.split("/")
    return F(int(numerator), int(denominator))


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def fields(line: str) -> dict[str, str]:
    return {
        item.split("=", 1)[0]: item.split("=", 1)[1]
        for item in line.split(";")[1:]
        if "=" in item
    }


def row_key(
    row: dict[str, str], apex_name: str, prefix_name: str
) -> tuple[tuple[int, ...], int, int, tuple[int, ...]]:
    return (
        tuple(map(int, row["E"].split(","))),
        int(row["rank"]),
        int(row[apex_name]),
        tuple(map(int, row[prefix_name].split(","))),
    )


def parse_hard_rows(
) -> dict[tuple[tuple[int, ...], int, int, tuple[int, ...]], dict[str, object]]:
    rows: dict[
        tuple[tuple[int, ...], int, int, tuple[int, ...]], dict[str, object]
    ] = {}
    for line in HARD_LEDGER.read_text().splitlines():
        if not line.startswith("HARD;"):
            continue
        row = fields(line)
        key = row_key(row, "apex", "prefix")
        top5 = tuple(
            (int(item.split(":", 1)[0]), parse_fraction(item.split(":", 1)[1]))
            for item in row["top5"].split(",")
        )
        require(len(top5) == 5, "hard row lost a singleton rank")
        require(
            all(top5[index][1] >= top5[index + 1][1] for index in range(4)),
            "hard singleton ranks are not nonincreasing",
        )
        require(len({label for label, _ in top5}) == 5, "hard top five repeat a label")
        require(key not in rows, "duplicate hard row")
        rows[key] = {
            "stratum": row["S"],
            "body": key[0],
            "rank": key[1],
            "apex": key[2],
            "prefix": key[3],
            "gate_size": int(row["K"]),
            "mass": parse_fraction(row["m"]),
            "components": int(row["r"]),
            "top5": top5,
        }
    require(len(rows) == 14_806, "hard-ledger row count changed")
    return rows


def clipped(value: F, upper: F) -> F:
    return max(F(0), min(upper, value))


def hunter_star(qs: tuple[F, ...], pair_cap: F) -> tuple[F, F]:
    require(len(qs) == 5, "Hunter star expects five singleton ranks")
    require(
        all(qs[index] >= qs[index + 1] for index in range(4)),
        "singleton ranks are not nonincreasing",
    )
    require(pair_cap >= qs[0], "pair cap fell below a singleton cap")
    upper = min(qs[0], pair_cap)
    candidates = {F(0), upper, clipped(pair_cap / 2, upper)}
    for singleton in qs[1:]:
        candidates.add(clipped(singleton, upper))
        candidates.add(clipped(pair_cap - singleton, upper))

    def objective(center: F) -> F:
        return center + sum(
            (min(center, singleton, pair_cap - center) for singleton in qs[1:]),
            F(0),
        )

    evaluated = tuple(sorted(candidates))
    for left, right in zip(evaluated, evaluated[1:]):
        midpoint = (left + right) / 2
        require(
            2 * objective(midpoint) == objective(left) + objective(right),
            "candidate list missed a Hunter-star breakpoint",
        )
    winner = min(evaluated, key=lambda center: (-objective(center), center))
    return objective(winner), winner


def parse_pair_rows(
    hard_rows: dict[
        tuple[tuple[int, ...], int, int, tuple[int, ...]], dict[str, object]
    ]
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    seen: set[tuple[tuple[int, ...], int, int, tuple[int, ...]]] = set()
    for line in PAIR_LEDGER.read_text().splitlines():
        if not line.startswith("PAIR;"):
            continue
        row = fields(line)
        key = row_key(row, "a", "P")
        require(key in hard_rows, "pair row has no hard-row mate")
        require(key not in seen, "duplicate pair row")
        seen.add(key)
        hard = hard_rows[key]
        top5 = hard["top5"]
        qs = tuple(value for _, value in top5)
        mass = parse_fraction(row["h"])
        pair_cap = parse_fraction(row["B2"])
        pair_margin = parse_fraction(row["mB2"])
        direct_margin = parse_fraction(row["mdirect"])
        require(mass == hard["mass"], "joined carrier mass changed")
        require(row["S"] == hard["stratum"], "joined stratum changed")
        require(int(row["K"]) == hard["gate_size"], "joined gate size changed")
        require(int(row["r"]) == hard["components"], "joined component count changed")
        require(
            (
                parse_fraction(row["q1"]),
                parse_fraction(row["q2"]),
                parse_fraction(row["q3"]),
                parse_fraction(row["q5"]),
            )
            == (qs[0], qs[1], qs[2], qs[4]),
            "joined singleton ranks changed",
        )
        require(pair_margin == 4 * mass / 7 - pair_cap, "pair margin identity failed")
        require(
            direct_margin == mass - qs[4] - 2 * pair_cap,
            "direct partition margin identity failed",
        )
        envelope, center = hunter_star(qs, pair_cap)
        require(envelope <= sum(qs, F(0)), "star exceeded scalar top-five invoice")
        require(
            envelope <= qs[4] + 2 * pair_cap,
            "star exceeded old pair-partition invoice",
        )
        rows.append(
            {
                **hard,
                "qs": qs,
                "pair_cap": pair_cap,
                "pair_margin": pair_margin,
                "direct_margin": direct_margin,
                "envelope": envelope,
                "center": center,
                "margin": mass - envelope,
            }
        )
    require(len(rows) == 14_806 and seen == set(hard_rows), "ledger join changed")
    return sorted(rows, key=lambda row: (row["body"], row["rank"], row["apex"]))


def digest_line(row: dict[str, object]) -> str:
    return (
        f"E={','.join(map(str, row['body']))};S={row['stratum']};"
        f"rank={row['rank']};a={row['apex']};"
        f"P={','.join(map(str, row['prefix']))};"
        f"h={ftext(row['mass'])};B2={ftext(row['pair_cap'])};"
        f"q={','.join(ftext(value) for value in row['qs'])};"
        f"aStar={ftext(row['center'])};G5={ftext(row['envelope'])};"
        f"margin={ftext(row['margin'])};"
        f"old={int(row['direct_margin'] > 0)};"
        f"exception={int(row['pair_margin'] <= 0)}\n"
    )


def main() -> None:
    require(file_sha256(HARD_LEDGER) == HARD_LEDGER_SHA256, "hard ledger changed")
    require(file_sha256(PAIR_LEDGER) == PAIR_LEDGER_SHA256, "pair ledger changed")
    rows = parse_pair_rows(parse_hard_rows())
    by_body: dict[tuple[int, ...], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        by_body[row["body"]].append(row)

    star_roots = tuple(
        sorted(
            body
            for body, body_rows in by_body.items()
            if all(row["margin"] > 0 for row in body_rows)
        )
    )
    old_roots = tuple(
        sorted(
            body
            for body, body_rows in by_body.items()
            if all(row["direct_margin"] > 0 for row in body_rows)
        )
    )
    require(set(old_roots) == THM2901_ROOTS, "old direct-terminal roots changed")
    one_hard = {
        body
        for body, body_rows in by_body.items()
        if len(body_rows) == 1 and body_rows[0]["direct_margin"] <= 0
    }
    prior_fifteen = (
        THM2895_ROOTS | THM2898_ROOTS | THM2899_ROOTS | THM2901_ROOTS
    )
    require(len(prior_fifteen) == 15, "prior fifteen-root union changed")
    require(len(one_hard) == 61, "THM2903 one-hard bank changed")
    require(one_hard & THM2902_ROOTS == THM2902_ROOTS, "one-hard overlap changed")
    current_union = prior_fifteen | THM2902_ROOTS | one_hard
    require(len(current_union) == 76, "proved union through THM2903 changed")
    additive_roots = tuple(sorted(set(star_roots) - current_union))
    new_union = current_union | set(star_roots)

    counts = (
        len(rows),
        sum(row["margin"] > 0 for row in rows),
        sum(row["direct_margin"] > 0 for row in rows),
        sum(row["margin"] > 0 and row["direct_margin"] <= 0 for row in rows),
        sum(row["margin"] == 0 for row in rows),
        sum(row["pair_margin"] <= 0 for row in rows),
        sum(row["pair_margin"] <= 0 and row["margin"] > 0 for row in rows),
        len(by_body),
        len(star_roots),
        len(old_roots),
        len(set(star_roots) - set(old_roots)),
        len(one_hard),
        len(set(star_roots) & one_hard),
        len(current_union),
        len(additive_roots),
        len(new_union),
        3432 - len(new_union),
        tuple(
            (
                stratum,
                sum(row["stratum"] == stratum for row in rows),
                sum(row["stratum"] == stratum and row["margin"] > 0 for row in rows),
            )
            for stratum in ("low", "one", "both")
        ),
        tuple(sorted(Counter(len(by_body[body]) for body in additive_roots).items())),
    )
    require(counts == EXPECTED_COUNTS, "Hunter-star counts changed")
    require(star_roots == EXPECTED_STAR_ROOTS, "Hunter-star root list changed")
    require(
        additive_roots == EXPECTED_ADDITIVE_ROOTS,
        "additive Hunter-star root list changed",
    )

    digest = hashlib.sha256()
    digest.update(b"LRC14/j6/all-hard/Hunter-star-G5/v1\n")
    for row in rows:
        digest.update(digest_line(row).encode())
    row_digest = digest.hexdigest()
    require(row_digest == EXPECTED_ROW_DIGEST, "Hunter-star row digest changed")

    positive_rows = [row for row in rows if row["margin"] > 0]
    failed_rows = [row for row in rows if row["margin"] <= 0]
    exception_rows = [row for row in rows if row["pair_margin"] <= 0]
    closest_positive = min(
        positive_rows,
        key=lambda row: (row["margin"], row["body"], row["rank"], row["apex"]),
    )
    closest_failure = max(
        failed_rows,
        key=lambda row: (row["margin"], tuple(-x for x in row["body"])),
    )
    best_exception = max(
        exception_rows,
        key=lambda row: (row["margin"], tuple(-x for x in row["body"])),
    )
    print("LRC14 all-hard Hunter-star envelope census")
    print(f"counts={counts}")
    print(f"star_roots={star_roots}")
    print(f"additive_roots={additive_roots}")
    print(
        "closest_positive="
        f"{closest_positive['body']};rank={closest_positive['rank']};"
        f"a={closest_positive['apex']};margin={ftext(closest_positive['margin'])};"
        f"G5={ftext(closest_positive['envelope'])};"
        f"aStar={ftext(closest_positive['center'])}"
    )
    print(
        "closest_failure="
        f"{closest_failure['body']};rank={closest_failure['rank']};"
        f"a={closest_failure['apex']};margin={ftext(closest_failure['margin'])};"
        f"G5={ftext(closest_failure['envelope'])};"
        f"aStar={ftext(closest_failure['center'])}"
    )
    print(
        "best_pair_cap_exception="
        f"{best_exception['body']};rank={best_exception['rank']};"
        f"a={best_exception['apex']};margin={ftext(best_exception['margin'])};"
        f"G5={ftext(best_exception['envelope'])};"
        f"aStar={ftext(best_exception['center'])}"
    )
    print(f"row_digest={row_digest}")
    print("mode=LOCKED")
    print(
        "scope=14806 scalar-hard marked suffixes;exact Hunter-star "
        "postprocessing;not LRC14"
    )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
