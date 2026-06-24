#!/usr/bin/env python3
"""HYP-2984/T1068: Farey mutations as LRC14 certificate schedulers.

The prompt asks for four Farey-sequence mutations:

  1. numerator multiplied by denominator,
  2. numerator and denominator added,
  3. denominator raised to the numerator,
  4. numerator raised to the denominator.

S130/HYP-2931 treated these as payloads attached to a fraction p/q.  This S164
pass adds the literal mutated-value reading and splices it into the current
HYP-2981/2982/2983 packet-certificate stack:

  product-value      (p*q)/q = p
  sum-value          (p+q)/q = 1 + p/q
  denominator-power  p/(q^p)
  numerator-power    (p^q)/q

The new observation is intentionally narrow.  The product-value collapse p is
not a theorem denominator, but after the unit-excess gate e=14p-q=1 it is a
certificate scheduler:

  p=1      q-parent / right-neighbor branch,
  p=2      C=27 petal or two-block branch,
  p>=3     K33 / state-lift branch.

The other mutations are guardrails: sum-value is an affine copy of M, while the
power mutations stress-test magnitude-blind quotients.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import log
from pathlib import Path
import sys


REPO = Path(__file__).resolve().parents[1]
THR = Fraction(1, 14)
AP = tuple(range(1, 14))


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


s130 = load_module(
    "s164_mutated_farey_s130",
    REPO / "04-computation" / "lrc14_mutated_farey_tournament_codex_s130.py",
)
s162 = load_module(
    "s164_fejer_scaffold_s162",
    REPO / "04-computation" / "lrc14_packet_fejer_interval_scaffold_codex_s162.py",
)


@dataclass(frozen=True)
class NamedRow:
    name: str
    source: str
    speeds: tuple[int, ...]


@dataclass(frozen=True)
class RowValue:
    row: NamedRow
    M: Fraction
    e: int
    p: int
    q: int

    @property
    def gap(self) -> Fraction:
        return self.M - THR

    @property
    def scheduler(self) -> str:
        if self.e == 0:
            return "boundary AP/GW"
        if self.e == 1 and self.p == 1:
            return "q-parent star"
        if self.e == 1 and self.p == 2:
            return "C27 petal/two-block"
        if self.e == 1 and self.p >= 3:
            return "K33/state-lift"
        if self.e > 1:
            return f"nonunit excess e={self.e}"
        return f"below floor e={self.e}"


def replace_one(drop: int, add: int) -> tuple[int, ...]:
    return tuple(sorted((set(AP) - {drop}) | {add}))


def replace_many(drops: tuple[int, ...], adds: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted((set(AP) - set(drops)) | set(adds)))


def exact_value(row: NamedRow) -> RowValue:
    M = s130.exact_M(row.speeds)
    return RowValue(row=row, M=M, e=14 * M.numerator - M.denominator, p=M.numerator, q=M.denominator)


def fmt_frac(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def format_big_log(row: RowValue, transform: str) -> str:
    p, q = row.p, row.q
    if transform == "product_value":
        return str(p)
    if transform == "sum_value":
        return fmt_frac(Fraction(p + q, q))
    if transform == "denominator_power":
        return f"log={log(p) - p * log(q):.3f}" if p > 0 else "-inf"
    if transform == "numerator_power":
        return f"log={q * log(p) - log(q):.3f}" if p > 0 else "-inf"
    raise ValueError(transform)


def transform_key(row: RowValue, transform: str) -> float | Fraction | int:
    p, q = row.p, row.q
    if transform == "product_value":
        return p
    if transform == "sum_value":
        return Fraction(p + q, q)
    if transform == "denominator_power":
        return log(p) - p * log(q) if p > 0 else float("-inf")
    if transform == "numerator_power":
        return q * log(p) - log(q) if p > 0 else float("-inf")
    raise ValueError(transform)


TRANSFORMS = (
    "product_value",
    "sum_value",
    "denominator_power",
    "numerator_power",
)


def s130_bank() -> list[RowValue]:
    rows = [NamedRow(r.label, r.family, r.speeds) for r in s130.candidate_rows()]
    return [exact_value(row) for row in rows]


def selected_rows() -> list[NamedRow]:
    rows = [
        NamedRow("AP", "boundary", AP),
        NamedRow("GW 12->24", "boundary", replace_one(12, 24)),
        NamedRow("residue liar 12->26", "q-witness loose", replace_one(12, 26)),
        NamedRow("near/K33 12->36", "K33 state-lift", replace_one(12, 36)),
        NamedRow("petal 10->20", "C27 petal", replace_one(10, 20)),
        NamedRow("petal 13->26", "C27 petal", replace_one(13, 26)),
        NamedRow("P10+GW", "two-block splice", replace_many((10, 12), (20, 24))),
        NamedRow("P10+K33", "two-block K33 splice", replace_many((10, 12), (20, 36))),
    ]
    for row in s162.default_rows():
        rows.append(NamedRow(row.name, row.source_family, row.speeds))

    seen: set[tuple[int, ...]] = set()
    out: list[NamedRow] = []
    for row in rows:
        if row.speeds not in seen:
            seen.add(row.speeds)
            out.append(row)
    return out


def rank_score(rows: list[RowValue], transform: str) -> tuple[int, int, int, Fraction]:
    agree = flip = tie = 0
    for i, j in combinations(range(len(rows)), 2):
        a, b = rows[i], rows[j]
        td = (a.gap > b.gap) - (a.gap < b.gap)
        ka = transform_key(a, transform)
        kb = transform_key(b, transform)
        pd = (ka > kb) - (ka < kb)
        if pd == 0:
            tie += 1
        elif pd == td:
            agree += 1
        else:
            flip += 1
    total = agree + flip + tie
    score = Fraction(2 * agree + tie, 2 * total) if total else Fraction(1)
    return agree, flip, tie, score


def low_frontier(rows: list[RowValue]) -> None:
    print("Low-frontier scheduler by literal product-value p=(p*q)/q")
    print("  filters: exact M in S130 AP/GW/petal/single-replacement bank")
    for bound_name, bound in (("M<=3/41", Fraction(3, 41)), ("M<=2/27", Fraction(2, 27))):
        front = sorted([r for r in rows if r.M <= bound], key=lambda r: (r.M, r.row.name))
        counts = Counter((r.p, r.scheduler) for r in front)
        print(f"  {bound_name}: rows={len(front)}")
        for key, count in sorted(counts.items()):
            print(f"    p={key[0]:<2d} {key[1]:24s} count={count}")
        for r in front[:12]:
            print(f"      {r.row.name:24s} M={fmt_frac(r.M):>7s} e={r.e:<2d} p={r.p:<2d} q={r.q:<4d} {r.scheduler}")
    print()


def transform_order_audit(rows: list[RowValue]) -> None:
    print("Mutation order audit on S130 row bank")
    print(f"  rows={len(rows)}; true key is M-1/14=e/(14q)")
    print("  product-value is literal (p*q)/q=p, not the old p*q payload.")
    for transform in TRANSFORMS:
        agree, flip, tie, score = rank_score(rows, transform)
        print(f"  {transform:19s} agree={agree:6d} flip={flip:6d} tie={tie:5d} score={fmt_frac(score)}")
    print("  readout: sum-value is theorem-safe because it is 1+M; product-value is not")
    print("  order-safe but it is a branch scheduler after e=1; power mutations are")
    print("  magnitude stress tests and should not be proof denominators.")
    print()


def selected_table() -> None:
    print("Selected packet rows under the four literal mutations")
    print(
        f"  {'row':34s} {'M':>8s} {'e':>3s} {'p':>3s} {'q':>5s} "
        f"{'product':>8s} {'sum':>9s} {'denpow':>12s} {'numpow':>12s} route"
    )
    for rv in sorted([exact_value(row) for row in selected_rows()], key=lambda r: (r.M, r.row.name)):
        print(
            f"  {rv.row.name:34s} {fmt_frac(rv.M):>8s} {rv.e:3d} {rv.p:3d} {rv.q:5d} "
            f"{format_big_log(rv, 'product_value'):>8s} "
            f"{format_big_log(rv, 'sum_value'):>9s} "
            f"{format_big_log(rv, 'denominator_power'):>12s} "
            f"{format_big_log(rv, 'numerator_power'):>12s} {rv.scheduler}"
        )
    print()


def unit_excess_chain(limit: int = 10) -> None:
    print("Unit-excess chain p/(14p-1)")
    print("  p is the product-value collapse and becomes a route index after e=1.")
    print(f"  {'p':>2s} {'M':>7s} {'q':>4s} {'sum':>4s} {'old p*q':>7s} {'literal route'}")
    for p in range(1, limit + 1):
        q = 14 * p - 1
        if p == 1:
            route = "q-parent star"
        elif p == 2:
            route = "C27 strip / petal"
        elif p == 3:
            route = "first K33 wall"
        else:
            route = "higher K33/state-lift"
        print(f"  {p:2d} {str(Fraction(p, q)):>7s} {q:4d} {p+q:4d} {p*q:7d} {route}")
    print()


def tournament_fingerprint(mask: int, n: int) -> dict[str, object]:
    def edge(i: int, j: int) -> bool:
        if i > j:
            return not edge(j, i)
        bit = 0
        for a, b in combinations(range(n), 2):
            if a == i and b == j:
                return bool((mask >> bit) & 1)
            bit += 1
        raise AssertionError

    outdeg = [sum(1 for j in range(n) if i != j and edge(i, j)) for i in range(n)]
    c3 = 0
    for a, b, c in combinations(range(n), 3):
        if edge(a, b) and edge(b, c) and edge(c, a):
            c3 += 1
        if edge(a, c) and edge(c, b) and edge(b, a):
            c3 += 1
    return {"score_hist": dict(sorted(Counter(outdeg).items())), "c3": c3}


def carrier_tournament() -> None:
    names = ["sum_value", "product_value", "denominator_power", "numerator_power"]
    # theorem_order, branch_scheduler, packet_side_channel, magnitude_stress,
    # false_quotient_detector.  Bigger is better; ties follow the listed order.
    scores = {
        "sum_value": (5, 1, 1, 0, 1),
        "product_value": (1, 5, 4, 1, 3),
        "denominator_power": (0, 2, 2, 4, 4),
        "numerator_power": (0, 3, 3, 5, 5),
    }
    mask = 0
    bit = 0
    for i, j in combinations(range(len(names)), 2):
        if scores[names[i]] >= scores[names[j]]:
            mask |= 1 << bit
        bit += 1
    fp = tournament_fingerprint(mask, len(names))
    print("Tournament Analysis")
    print("  vertices: the four literal Farey mutations")
    print("  pairwise observable: retained proof role vector")
    print("    (order safety, branch scheduling, packet side-channel, magnitude stress, false-quotient detection)")
    print("  switch/gauge: larger role vector wins; ties use listed Hamiltonian order")
    print(f"  fingerprint: score_hist={fp['score_hist']} directed_3cycles={fp['c3']}")
    print("  Hamiltonian path: sum_value > product_value > numerator_power > denominator_power")
    print("  challenged assumption: the numerator-times-denominator mutation is useless")
    print("    as a value, but after the e=1 quotient it is exactly the numerator p,")
    print("    which is the C27/K33 route switch.")
    print()


def proof_readout() -> None:
    print("Proof-route readout")
    print("  1. Keep exact M=p/q and e=14p-q first.")
    print("  2. If e=0, route to AP/GW boundary equality.")
    print("  3. If e=1, the literal product mutation collapses to p and schedules:")
    print("       p=1 -> q-parent/right-neighbor, p=2 -> C27 petal/two-block,")
    print("       p>=3 -> K33/state-lift/Fejer packet.")
    print("  4. Use sum-value only as an affine check that a quotient preserved M.")
    print("  5. Use the two power mutations as adversarial stress tests: a claimed")
    print("     tournament or scalar invariant surviving those distortions is probably")
    print("     too magnitude-blind and needs endpoint/Ramanujan/Fejer labels restored.")
    print("  This integrates HYP-2931/HYP-2940 with HYP-2981/HYP-2982/HYP-2983:")
    print("  Farey mutations are front-end schedulers; Fejer/Ramanujan/Kaczynski")
    print("  packets are the certificate back end.")


def main() -> None:
    print("HYP-2984/T1068 LRC14 FAREY MUTATION CERTIFICATE SCHEDULER")
    print("=" * 78)
    print("Assumption challenge")
    print("  considered vertices: runners, Farey fractions, mutated values, denominators,")
    print("    endpoint owners, Fejer packet fibers, Ramanujan phases, and proof obligations.")
    print("  chosen vertices: four literal Farey mutations used as proof-route carriers.")
    print("  preserved predicate: binding gap M-1/14 after exact e=14p-q is retained.")
    print("  destroyed information: wall geometry; it must be restored by packet labels.")
    print()
    unit_excess_chain()
    rows = s130_bank()
    low_frontier(rows)
    transform_order_audit(rows)
    selected_table()
    carrier_tournament()
    proof_readout()


if __name__ == "__main__":
    main()
