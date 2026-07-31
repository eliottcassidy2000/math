#!/usr/bin/env python3
"""Join the ranked-H1 cutoff tails to the exact THM-2905 G5 rows."""

from __future__ import annotations

import hashlib
import importlib.util
from collections import Counter
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
G5_PATH = (
    ROOT
    / "04-computation/lrc14_j6_all_hard_hunter_star_envelope_census_thm2905.py"
)
EXPECTED_G5_SHA256 = (
    "2755c0dbfa2c2a59099e17fc94b1fb702a81b9d74ce6a1de66abc6022654fde8"
)
S2 = F(99, 70)
MAX_CUTOFF = 15_000


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ceiling(value: F) -> int:
    return -((-value.numerator) // value.denominator)


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def load_g5():
    require(G5_PATH.is_file(), "THM-2905 candidate verifier is absent")
    require(sha256(G5_PATH) == EXPECTED_G5_SHA256, "THM-2905 candidate changed")
    spec = importlib.util.spec_from_file_location("thm2905_g5", G5_PATH)
    require(spec is not None and spec.loader is not None, "cannot load G5 module")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def row_identity(row: dict[str, object]) -> tuple[object, ...]:
    return row["body"], row["K"], row["rank"], row["apex"], row["prefix"]


def row_line(row: dict[str, object]) -> str:
    return (
        f"E={row['body']};K={row['K']};rank={row['rank']};"
        f"a={row['apex']};P={row['prefix']};N={row['cutoff']};"
        f"h={ftext(row['h'])};"
        f"q={tuple(ftext(value) for value in row['qs'])};"
        f"B2={ftext(row['b2'])};G5={ftext(row['g5'])};"
        f"astar={ftext(row['maximizing_a'])};"
        f"margin={ftext(row['g5_margin'])};"
        f"direct={int(row['direct_closed'])};G={int(row['g5_closed'])}\n"
    )


def main() -> None:
    module = load_g5()
    all_rows = module.join_rows()
    eligible: list[dict[str, object]] = []
    tails: list[dict[str, object]] = []
    for row in all_rows:
        epsilon = 6 * row["h"] / 7 - sum(row["qs"][:4], F(0))
        if epsilon <= 0:
            continue
        cutoff = ceiling(S2 * row["components"] / (7 * epsilon)) - 1
        record = {**row, "epsilon": epsilon, "cutoff": cutoff}
        eligible.append(record)
        if cutoff > MAX_CUTOFF:
            tails.append(record)

    require(len(eligible) == 6_180, "rank-four/H1 eligible count changed")
    require(len(tails) == 181, "cutoff-tail count changed")
    g5_closed = [row for row in tails if row["g5_closed"]]
    survivors = [row for row in tails if not row["g5_closed"]]
    direct_closed = [row for row in tails if row["direct_closed"]]
    require(
        all(row["g5_closed"] for row in direct_closed),
        "G5 failed to dominate a direct tail closure",
    )
    gains = [row for row in g5_closed if not row["direct_closed"]]
    roots = {row["body"] for row in tails}
    closed_roots = {
        body
        for body in roots
        if all(
            row["g5_closed"]
            for row in tails
            if row["body"] == body
        )
    }

    canonical = hashlib.sha256(b"LRC14/j6/H1-cutoff-tail-G5/v1\n")
    for row in sorted(tails, key=row_identity):
        canonical.update(row_line(row).encode())
    survivor_hash = hashlib.sha256(b"LRC14/j6/H1-cutoff-tail-G5-survivors/v1\n")
    for row in sorted(survivors, key=row_identity):
        survivor_hash.update(row_line(row).encode())

    closest_closed = min(
        g5_closed, key=lambda row: (row["g5_margin"], row_identity(row))
    )
    closest_open = max(
        survivors, key=lambda row: (row["g5_margin"], row_identity(row))
    )
    largest_cutoff_open = max(
        survivors, key=lambda row: (row["cutoff"], row_identity(row))
    )
    print("LRC14 j6 rank-four/H1 181-cutoff-tail exact G5 partition")
    print(
        "counts="
        + repr(
            (
                len(eligible),
                len(tails),
                len(direct_closed),
                len(g5_closed),
                len(gains),
                len(survivors),
                len(roots),
                len(closed_roots),
            )
        )
    )
    print(
        "ranks="
        + repr(
            tuple(
                sorted(
                    (
                        rank,
                        sum(row["rank"] == rank for row in tails),
                        sum(
                            row["rank"] == rank and row["g5_closed"]
                            for row in tails
                        ),
                        sum(
                            row["rank"] == rank and not row["g5_closed"]
                            for row in tails
                        ),
                    )
                    for rank in Counter(row["rank"] for row in tails)
                )
            )
        )
    )
    print(f"closest_closed={row_line(closest_closed).strip()}")
    print(f"closest_open={row_line(closest_open).strip()}")
    print(f"largest_cutoff_open={row_line(largest_cutoff_open).strip()}")
    print(f"tail_ledger_sha256={canonical.hexdigest()}")
    print(f"survivor_ledger_sha256={survivor_hash.hexdigest()}")
    print("G5_SURVIVORS")
    for row in sorted(
        survivors,
        key=lambda item: (
            -item["g5_margin"],
            item["cutoff"],
            row_identity(item),
        ),
    ):
        print(row_line(row).strip())
    print("scope=181 rank-four/H1 cutoff tails;exact G5 join;not LRC14")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
