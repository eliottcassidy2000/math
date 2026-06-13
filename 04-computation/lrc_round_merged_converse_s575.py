#!/usr/bin/env python3
"""
lrc_round_merged_converse_s575.py

oraclebox1 connector, 2026-06-03 S575

Exact connector audit for the LRC round-tournament body under the repo's
primary converse/merged quotient.

Background:
  HYP-1998/S574 establish that open-time LRC runner sub-tournaments on m
  vertices are exactly ROUND tournaments, counted by A000016.  HYP-2086/HYP-2087
  say the hard LRC regime should be read on the Burnside fixed/dihedral side.

New connector step:
  Time reversal sends the half-turn runner tournament T(t) to T(t)^op, so the
  open round body should be quotiented by converse before being attached to the
  merged metagraph G_m/Z_2.  This script reuses the already validated S574
  round generator and exact individualization-refinement canonical labeling,
  then pairs every round class with its converse.

For n=14 LRC, m=n-1=13.  The exact result is:
  round classes = 316
  self-converse round classes = 64
  converse-merged round classes = 190
"""

from __future__ import annotations

import importlib.util
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
S574 = ROOT / "04-computation" / "lrc_round_count_m89_s574.py"


@dataclass(frozen=True)
class Row:
    m: int
    valid_dvectors: int
    round_classes: int
    a000016: int
    self_converse_round: int
    predicted_self_converse: int
    merged_converse_round: int


def load_s574():
    spec = importlib.util.spec_from_file_location("s574_rounds", S574)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {S574}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def opposite_adj(adj: list[list[int]]) -> list[list[int]]:
    m = len(adj)
    return [[0 if i == j else adj[j][i] for j in range(m)] for i in range(m)]


def class_representatives(rounds, m: int) -> dict[tuple[int, ...], list[list[int]]]:
    reps: dict[tuple[int, ...], list[list[int]]] = {}
    for d in rounds.valid_dvectors(m):
        adj = rounds.build_adj(d, m)
        reps.setdefault(rounds.canon(adj, m), adj)
    return reps


def row(rounds, m: int) -> Row:
    reps = class_representatives(rounds, m)
    self_converse = 0
    merged: set[tuple[int, ...]] = set()
    for canon, adj in reps.items():
        converse = rounds.canon(opposite_adj(adj), m)
        if canon == converse:
            self_converse += 1
        merged.add(min(canon, converse))

    return Row(
        m=m,
        valid_dvectors=len(rounds.valid_dvectors(m)),
        round_classes=len(reps),
        a000016=rounds.A000016(m),
        self_converse_round=self_converse,
        predicted_self_converse=2 ** ((m - 1) // 2),
        merged_converse_round=len(merged),
    )


def main() -> None:
    rounds = load_s574()
    rows = [row(rounds, m) for m in range(3, 14)]

    print("Round LRC body under converse/merged quotient")
    print("oraclebox1 connector, 2026-06-03 S575")
    print()
    print(
        f"{'m':>3} {'valid-d':>8} {'round':>7} {'A000016':>8} "
        f"{'SC-round':>9} {'2^floor':>9} {'merged':>7}"
    )
    for r in rows:
        print(
            f"{r.m:>3} {r.valid_dvectors:>8} {r.round_classes:>7} "
            f"{r.a000016:>8} {r.self_converse_round:>9} "
            f"{r.predicted_self_converse:>9} {r.merged_converse_round:>7}"
        )

    print()
    print("Checks:")
    print(f"  round == A000016 for all rows: {all(r.round_classes == r.a000016 for r in rows)}")
    print(
        "  SC-round == 2^floor((m-1)/2) through m=13: "
        f"{all(r.self_converse_round == r.predicted_self_converse for r in rows)}"
    )
    print(
        "  merged == (round + SC-round)/2 for all rows: "
        f"{all(r.merged_converse_round * 2 == r.round_classes + r.self_converse_round for r in rows)}"
    )
    print()
    target = rows[-1]
    print("n=14 handoff (m=13 runner sub-tournament):")
    print(
        f"  {target.round_classes} round/A000016 classes collapse to "
        f"{target.merged_converse_round} converse-merged classes;"
    )
    print(
        f"  {target.self_converse_round} are fixed by converse, and "
        f"{(target.round_classes - target.self_converse_round) // 2} "
        "are non-fixed converse pairs."
    )
    print()
    print("Connector interpretation:")
    print("  The open LRC necklace body should enter the merged metagraph through")
    print("  this 190-node quotient at n=14, not through all 316 round classes and")
    print("  not through the full A000568(13) ambient set.  Boundary/fixed classes")
    print("  are the remaining seam where HYP-2088's D/U/N blocker labels should be")
    print("  attached before binary time-word symmetry forgets ownership.")


if __name__ == "__main__":
    main()
