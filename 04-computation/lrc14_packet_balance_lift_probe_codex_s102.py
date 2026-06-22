#!/usr/bin/env python3
"""
Lift the HYP-2883 repeated-packet balance into actual LRC reciprocal tails.

The finite HYP-2632 packet graph has exact local balance

    loop(a) + sum_b edge(a,b) = 0

on residue vertices {0,2,3,4,5,6}.  This script asks what remains true after
choosing integer lifts of those residues and evaluating the actual support-six
reciprocal hyperplane sums.

The point is not to prove the tail bound.  The point is to expose the next
LRC proof object: a divergence defect measuring how far the finite signed
current fails to lift before low-height wall deletion and Abel summation.
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_support6_residue_cusp_codex_s12 as s12  # noqa: E402
import lrc14_two_large_dedekind_phase_codex_s23 as s23  # noqa: E402

AMBIENT_D = 9
HMAX = 12
CHECKPOINTS = (6, 8, 10, 12)
VERTICES = (0, 2, 3, 4, 5, 6)
CORE = (1, 8, 15, 22)
LOW_HEIGHT = 2


@dataclass(frozen=True)
class LiftGauge:
    name: str
    edge_layer: int
    loop_first_layer: int

    def lift(self, residue: int, layer: int) -> int:
        if residue == 0:
            return 7 * (layer + 1)
        return residue + 7 * layer

    def loop_support(self, residue: int) -> tuple[int, ...]:
        a = self.lift(residue, self.loop_first_layer)
        b = self.lift(residue, self.loop_first_layer + 1)
        return tuple(sorted(CORE + (a, b)))

    def edge_support(self, a: int, b: int) -> tuple[int, ...]:
        x = self.lift(a, self.edge_layer)
        y = self.lift(b, self.edge_layer)
        return tuple(sorted(CORE + (x, y)))


GAUGES = (
    LiftGauge("start_aligned", edge_layer=0, loop_first_layer=0),
    LiftGauge("raised_pair", edge_layer=1, loop_first_layer=1),
)


def chi7(x: int) -> int:
    x %= 7
    if x == 0:
        return 0
    return 1 if x in {1, 2, 4} else -1


def q_selector(a: int, b: int) -> int:
    return (a * b * (1 + 3 * ((a + b) % 7)) - 1) % 7


def loop_weight(residue: int) -> int:
    if residue == 0:
        return -4
    return (-43 - 7 * chi7(residue)) // 2


def edge_weight(a: int, b: int) -> int:
    if (a + b) % 7 == 2:
        return 0
    return 8 if chi7(q_selector(a, b)) == 1 else 1


def section(title: str) -> None:
    print("\n" + "=" * 92)
    print(title)
    print("=" * 92)


@lru_cache(None)
def lift_rows(support: tuple[int, ...]) -> tuple[tuple[int, int, complex, float, float], ...]:
    by_shell, _ = s12.six_support_shells(support, AMBIENT_D, HMAX)
    return tuple(s12.cumulative_rows(by_shell, CHECKPOINTS))


def lift_at(support: tuple[int, ...], h: int) -> complex:
    for row_h, _count, signed, _abs_k, _ratio in lift_rows(support):
        if row_h == h:
            return signed
    raise KeyError(h)


def low_height_relation_count(support: tuple[int, ...], height: int = LOW_HEIGHT) -> int:
    count = 0
    coeffs = range(-height, height + 1)

    def rec(i: int, dot: int, nonzero: bool) -> None:
        nonlocal count
        if i == len(support):
            if nonzero and dot == 0:
                count += 1
            return
        for c in coeffs:
            rec(i + 1, dot + c * support[i], nonzero or c != 0)

    rec(0, 0, False)
    # Count +/- relation as one wall direction.
    return count // 2


def finite_balance_report() -> None:
    section("FINITE PACKET BALANCE")
    print("All weights are in HYP-2632 units U before reciprocal lifting.")
    print(f"{'v':>2} {'loop':>6} {'incident':>8} {'balance':>8}")
    failures = []
    for v in VERTICES:
        incident = sum(edge_weight(min(v, u), max(v, u)) for u in VERTICES if u != v)
        balance = loop_weight(v) + incident
        if balance:
            failures.append((v, balance))
        print(f"{v:>2} {loop_weight(v):>6} {incident:>8} {balance:>8}")
    print(f"finite balance failures: {failures}")


def gauge_support_report() -> None:
    section("INTEGER LIFT GAUGES")
    for gauge in GAUGES:
        print(f"\n{gauge.name}")
        print(f"{'packet':>10} {'support':>28} {'low-height walls':>17}")
        for v in VERTICES:
            support = gauge.loop_support(v)
            print(f"{'loop '+str(v):>10} {str(support):>28} {low_height_relation_count(support):>17}")
        for i, a in enumerate(VERTICES):
            for b in VERTICES[i + 1 :]:
                if edge_weight(a, b) == 0:
                    continue
                support = gauge.edge_support(a, b)
                print(
                    f"{str((a,b)):>10} {str(support):>28} "
                    f"{low_height_relation_count(support):>17}"
                )


def divergence_report() -> None:
    section("LIFTED LOCAL DIVERGENCE")
    print(
        "For each residue vertex a, compute div_H(a) = lifted_loop_H(a) + "
        "sum_b lifted_edge_H(a,b).  Finite packets have div=0 exactly; nonzero "
        "values here are the reciprocal-lift defect that a proof must delete or "
        "sum by parts."
    )
    print()
    for gauge in GAUGES:
        print(f"Gauge: {gauge.name}")
        header = f"{'H':>3} {'max|div|':>13} {'L1 div':>13} {'sum div':>13}"
        for v in VERTICES:
            header += f" {str(v):>12}"
        print(header)
        for h in CHECKPOINTS:
            divs: dict[int, complex] = {}
            for v in VERTICES:
                div = lift_at(gauge.loop_support(v), h)
                for u in VERTICES:
                    if u == v:
                        continue
                    a, b = sorted((u, v))
                    if edge_weight(a, b) == 0:
                        continue
                    div += lift_at(gauge.edge_support(a, b), h)
                divs[v] = div
            max_div = max(abs(z.real) for z in divs.values())
            l1_div = sum(abs(z.real) for z in divs.values())
            sum_div = sum(z.real for z in divs.values())
            line = f"{h:>3} {max_div:>13.6g} {l1_div:>13.6g} {sum_div:>13.6g}"
            for v in VERTICES:
                line += f" {divs[v].real:>12.5g}"
            print(line)
        print()


def weighted_balance_report() -> None:
    section("GLOBAL INCIDENCE BALANCE")
    print(
        "The finite graph has sum(loop)+2*sum(edge)=0.  The lifted analogue "
        "tests whether incidence-counting, not scalar edge-counting, is the "
        "right analytic invariant."
    )
    print(f"{'gauge':>14} {'H':>3} {'loop sum':>13} {'2*edge sum':>13} {'total':>13}")
    for gauge in GAUGES:
        for h in CHECKPOINTS:
            loop_sum = sum(lift_at(gauge.loop_support(v), h).real for v in VERTICES)
            edge_sum = 0.0
            for i, a in enumerate(VERTICES):
                for b in VERTICES[i + 1 :]:
                    if edge_weight(a, b) == 0:
                        continue
                    edge_sum += lift_at(gauge.edge_support(a, b), h).real
            print(
                f"{gauge.name:>14} {h:>3} {loop_sum:>13.6g} "
                f"{2*edge_sum:>13.6g} {loop_sum + 2*edge_sum:>13.6g}"
            )


def tournament_analysis_note() -> None:
    section("TOURNAMENT ANALYSIS NOTE")
    print("Candidate vertices considered:")
    print(
        "  runners, raw speeds, residue packets, integer lifts, finite loops, "
        "finite edges, low-height wall motifs, additive-frequency shells, "
        "divergence defects, and proof obligations."
    )
    print("Pairwise observable:")
    print("  whether a quotient preserves local signed-current balance after reciprocal lift.")
    print("Switch/gauge:")
    print("  start-aligned lift versus raised-pair lift; both keep the mod-7 packet graph fixed.")
    print("Tie Hamiltonian path:")
    print(
        "  low_height_wall_deletion > lifted_divergence_defect > "
        "additive_frequency_Abel_sum > finite_packet_balance > raw_runner_labels"
    )
    print(
        "Challenged assumption: the exact finite balance automatically survives "
        "integer lifting.  The data below measure the defect that must be "
        "absorbed by wall deletion and summation by parts."
    )


def main() -> None:
    section("LRC14 PACKET BALANCE LIFT PROBE - CODEX S102")
    print(f"ambient d={AMBIENT_D}; Hmax={HMAX}; checkpoints={CHECKPOINTS}")
    print(f"finite packet unit U={s23.UNIT:.12g}")
    finite_balance_report()
    gauge_support_report()
    divergence_report()
    weighted_balance_report()
    tournament_analysis_note()
    section("S102 READING")
    print(
        "The finite signed-current identity is real, but integer lifts introduce "
        "a measurable divergence defect.  This is good news: it turns the vague "
        "residue-lift equidistribution gap into a concrete local-divergence "
        "lemma.  The LRC proof should delete low-height wall directions, then "
        "prove that the remaining lifted divergence is Abel-summable in the "
        "additive-frequency shells of HYP-2636."
    )


if __name__ == "__main__":
    main()
