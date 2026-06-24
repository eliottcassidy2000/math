#!/usr/bin/env python3
"""S131: Farey products as complete-bipartite obstruction ledgers for LRC14.

The prompt's new dictionary is

    p/q  ->  K_{p,q},      p*q = |E(K_{p,q})|.

This script keeps the S130 guardrail: p*q is not the theorem denominator for
LRC14.  The theorem scale is still q in

    M(S)-1/14 = (14p-q)/(14q).

The product ledger becomes useful only after it is read as a graph:
K_{p,q} is nonplanar exactly when min(p,q) >= 3, because it then contains
K_{3,3}.  On the LRC14 unit-excess chain p/(14p-1), this means the first
nonplanar product ledger is p=3, namely 3/41, the known near-miss row
{1,...,11,13,36}.  That isolates a new proof split:

  tight floor       e=0
  coarse parent     e=1, p=1     K_{1,13}, planar
  two-block strip   e=1, p=2     K_{2,27}, planar
  K3,3 wall         e=1, p>=3    nonplanar product-minor certificate

The computation below audits that split on the S130 row bank and records a
Tournament Analysis whose vertices are proof ledgers rather than runners.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from math import gcd
from pathlib import Path
import sys


REPO = Path(__file__).resolve().parents[1]
THR = F(1, 14)


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


s130 = load_module(
    "s130_mutated_farey",
    REPO / "04-computation" / "lrc14_mutated_farey_tournament_codex_s130.py",
)
s127 = s130.s127
s124 = s130.s124


@dataclass(frozen=True)
class BipartiteLedger:
    p: int
    q: int

    @property
    def edges(self) -> int:
        return self.p * self.q

    @property
    def contains_k33(self) -> bool:
        return min(self.p, self.q) >= 3

    @property
    def planar(self) -> bool:
        return not self.contains_k33

    @property
    def rank(self) -> str:
        if self.contains_k33:
            return "K33-wall"
        if min(self.p, self.q) == 2:
            return "two-block"
        if min(self.p, self.q) == 1:
            return "star"
        return "empty"


def farey(order: int) -> list[F]:
    out: list[F] = []
    for q in range(1, order + 1):
        for p in range(q + 1):
            if gcd(p, q) == 1:
                out.append(F(p, q))
    return sorted(set(out))


def farey_level_table(max_order: int = 8) -> None:
    print("[Ordinary Farey levels: first complete-bipartite nonplanarity]")
    prev: set[F] = set()
    first_wall: tuple[int, F] | None = None
    for order in range(1, max_order + 1):
        nodes = set(farey(order))
        new = sorted(x for x in nodes - prev if x > 0)
        walls = [x for x in new if BipartiteLedger(x.numerator, x.denominator).contains_k33]
        if walls and first_wall is None:
            first_wall = (order, walls[0])
        summary = ", ".join(str(x) for x in new) if new else "-"
        wall_summary = ", ".join(str(x) for x in walls) if walls else "-"
        print(f"  F_{order:<2d} new={{ {summary} }}  K33-new={{ {wall_summary} }}")
        prev = nodes
    assert first_wall is not None
    order, frac = first_wall
    led = BipartiteLedger(frac.numerator, frac.denominator)
    print(
        f"  first wall: F_{order} has {frac} -> K_{{{led.p},{led.q}}} "
        f"with {led.edges} edges and a K_{{3,3}} subgraph"
    )
    print()


def lrc_child_chain(limit: int = 8) -> None:
    print("[LRC14 unit-excess chain as K_{p,q}]")
    print("  e=14p-q=1, so q=14p-1 and M=p/q is a Farey neighbor of 1/14.")
    print(f"  {'p':>2s} {'M':>7s} {'q':>4s} {'p*q':>5s} {'rank':>10s} {'proof role'}")
    for p in range(1, limit + 1):
        q = 14 * p - 1
        M = F(p, q)
        led = BipartiteLedger(p, q)
        if p == 1:
            role = "right parent / q-threshold loose"
        elif p == 2:
            role = "planar two-block strip"
        elif p == 3:
            role = "first K3,3 wall; S128 near-miss 12->36"
        else:
            role = "higher K3,3 wall"
        print(f"  {p:2d} {str(M):>7s} {q:4d} {led.edges:5d} {led.rank:>10s}  {role}")
    print()


def row_label(row) -> str:
    return row.label


def row_bank_bipartite_analysis() -> None:
    rows = s130.candidate_rows()
    Ms, _analysis = s130.row_bank_analysis(rows)
    buckets: Counter[tuple[str, str]] = Counter()
    examples: dict[tuple[str, str], list[str]] = {}
    unit_examples: list[tuple[str, F, BipartiteLedger]] = []

    for row, M in zip(rows, Ms, strict=True):
        e = s130.farey_excess(M)
        led = BipartiteLedger(M.numerator, M.denominator)
        if e == 0:
            gap_class = "tight-floor"
        elif e == 1:
            gap_class = "unit-excess"
        else:
            gap_class = "nonunit-excess"
        key = (gap_class, led.rank)
        buckets[key] += 1
        examples.setdefault(key, [])
        if len(examples[key]) < 4:
            examples[key].append(row.label)
        if e == 1:
            unit_examples.append((row.label, M, led))

    print("[S130 row-bank bipartite ledger: 749 AP/GW/petal/single-replacement rows]")
    print("  rows are bucketed by exact Farey excess and K_{p,q} product rank.")
    for key in sorted(buckets):
        gap_class, rank = key
        print(
            f"  {gap_class:15s} {rank:10s} count={buckets[key]:4d} "
            f"examples={examples[key]}"
        )
    print()

    print("[Unit-excess rows in the bank]")
    print(f"  {'row':24s} {'M':>7s} {'K_{p,q}':>10s} {'edges':>6s} {'rank':>10s}")
    for label, M, led in sorted(unit_examples, key=lambda x: (x[1].numerator, x[1].denominator, x[0])):
        print(
            f"  {label:24s} {str(M):>7s} K_{{{led.p},{led.q}}} "
            f"{led.edges:6d} {led.rank:>10s}"
        )
    print()


def selected_rows() -> None:
    ap = tuple(range(1, 14))
    rows = [
        ("AP", s130.Row("AP", ap, "known tight")),
        ("GW 12->24", s130.Row("GW 12->24", tuple(list(range(1, 12)) + [13, 24]), "known tight")),
        ("residue-liar 12->26", s130.Row("residue-liar 12->26", tuple(list(range(1, 12)) + [13, 26]), "loose")),
        ("near-miss 12->36", s130.Row("near-miss 12->36", tuple(list(range(1, 12)) + [13, 36]), "loose")),
        ("petal 10->20", s130.Row("petal 10->20", tuple(sorted((set(ap) - {10}) | {20})), "loose")),
        ("petal 13->26", s130.Row("petal 13->26", tuple(sorted((set(ap) - {13}) | {26})), "loose")),
    ]
    print("[Selected LRC14 rows]")
    print(f"  {'row':24s} {'M':>7s} {'e':>3s} {'q(S)':>4s} {'K_{p,q}':>10s} {'rank':>10s}")
    for _label, row in rows:
        M = s130.exact_M(row.speeds)
        e = s130.farey_excess(M)
        led = BipartiteLedger(M.numerator, M.denominator)
        print(
            f"  {row.label:24s} {str(M):>7s} {e:3d} {s124.q_threshold(row.speeds):4d} "
            f"K_{{{led.p},{led.q}}} {led.rank:>10s}"
        )
    print()


def ledger_tournament() -> None:
    names = [
        "q-binding",
        "sum-recursion",
        "Kpq-product",
        "K33-wall",
        "power-stress",
        "raw-iso",
    ]
    # Scores: theorem safety, additive locality, multiplicative side-channel,
    # obstruction signal, magnitude stress, and known false-positive resistance.
    scores = {
        "q-binding": (5, 4, 2, 2, 1, 5),
        "sum-recursion": (4, 5, 2, 1, 1, 4),
        "Kpq-product": (2, 2, 5, 3, 2, 2),
        "K33-wall": (1, 1, 4, 5, 1, 2),
        "power-stress": (0, 0, 2, 1, 5, 1),
        "raw-iso": (0, 1, 0, 0, 0, 0),
    }
    mask = 0
    bit = 0
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            if scores[names[i]] >= scores[names[j]]:
                mask |= 1 << bit
            bit += 1
    fp = s127.tournament_fingerprint(mask, len(names))
    print("[Tournament Analysis on proof ledgers]")
    print("  vertices: proof ledgers, not runners")
    print("  pairwise observable: lexicographic role score")
    print("    (theorem safety, additive locality, product side-channel, obstruction signal,")
    print("     magnitude stress, false-positive resistance)")
    print("  switch/gauge: larger role score wins; ties use listed Hamiltonian path")
    print("  preserves: which quotient can serve which part of the LRC14 proof")
    print("  destroys: exact runner geometry and continuous safe intervals")
    print(
        f"  fingerprint: score_hist={fp['score_hist']} c3={fp['c3']} "
        f"scc={fp['scc']} hp={fp['hp']}"
    )
    print(f"  Hamiltonian order: {' > '.join(names)}")
    print()


def proof_readout() -> None:
    print("[Proof readout]")
    print("  The product payload p*q now has a precise job: it is the edge count of")
    print("  K_{p,q}.  It must not replace q in the analytic gap, but it can mark")
    print("  the incidence complexity of a Farey escape.")
    print("  General Farey: F_4 is the first level where the product graph sees the")
    print("  Kuratowski K_{3,3} side, via 3/4 -> K_{3,4}.")
    print("  LRC14 unit-excess chain: the first nonplanar product graph is 3/41 ->")
    print("  K_{3,41}.  That is exactly the S128 near-miss 12->36.")
    print("  New local proof split:")
    print("    e=0: AP/GW floor candidates;")
    print("    e=1,p=1: coarse q-threshold parent 1/13;")
    print("    e=1,p=2: planar two-block strip 2/27, needs the petal/two-block kill;")
    print("    e=1,p>=3: K3,3 product-minor wall, a three-owner obstruction target.")
    print("  The next theorem should not claim nonplanarity is a contradiction by itself.")
    print("  It should prove that any remaining q=14 non-AP/GW atom either falls into")
    print("  the p=2 two-block strip handled by petal rigidity, or reaches p>=3 where")
    print("  the K3,3 ledger supplies the finite three-owner obstruction packet.")


def main() -> None:
    print("S131 LRC14 FAREY-BIPARTITE OBSTRUCTION LEDGER")
    print("=" * 78)
    print("[Assumption challenge]")
    print("  considered vertices: fractions, bipartition sides, incidence edges,")
    print("    K3,3 minors, LRC rows, residues, wall crossings, and proof ledgers.")
    print("  chosen quotient: M(S)=p/q -> K_{p,q}, while retaining excess e=14p-q.")
    print("  preserved predicate: product-side incidence depth of a Farey escape.")
    print("  destroyed information: exact safe-time intervals and row-specific geometry.")
    print()
    farey_level_table()
    lrc_child_chain()
    selected_rows()
    row_bank_bipartite_analysis()
    ledger_tournament()
    proof_readout()


if __name__ == "__main__":
    main()
