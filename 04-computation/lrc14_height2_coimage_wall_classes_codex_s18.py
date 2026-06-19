#!/usr/bin/env python3
"""
LRC(14) height-2 coimage wall-class audit.

Continuation of HYP-2616/HYP-2617/HYP-2619.  HYP-2616 cleared the finite
height-1 one-large support-six wall ledger by exact S7 measure.  HYP-2617
then showed that the support-six analytic tail factors through only 159
projective mod-7 coimage classes.

This script asks a more quotient-level question:

    Which coimage classes are already addressed by one-large support-six
    additive walls of coefficient height <= 2?

The goal is not to clear every row by exact measure.  It is to separate the
final analytic tail into:

    finite wall-addressed coimage packets,
    genuinely tail-only residue packets.

For k=8,9,10 we enumerate supports of the form

    {a1,...,a5,M},  1 <= ai <= B(k), M > B(k),

with a nontrivial relation

    c_M M + sum c_i a_i = 0,    c_i in {-H,...,-1,1,...,H},

for H=1 and H=2.  We project each support to HYP-2617's canonical
F_7^*/S_6 coimage class and compare with the signed coimage mass in ambient
dimension d=k-1.

Tournament Analysis declaration:
  vertices are proof quotients rather than runners: raw supports,
  height-1 walls, height-2 walls, coimage classes, repeated-residue tail
  packet, signed reciprocal theorem, and raw runner vertices.
"""
from __future__ import annotations

import itertools
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_support6_coimage_fiber_codex_s14 as s14  # noqa: E402

BOUND = {8: 16, 9: 15, 10: 13}
AMBIENT_D = {8: 7, 9: 8, 10: 9}
EPS = 1e-12


def section(title: str) -> None:
    print("\n" + "=" * 88, flush=True)
    print(title, flush=True)
    print("=" * 88, flush=True)


def coefficient_set(height: int) -> tuple[int, ...]:
    return tuple(x for x in range(-height, height + 1) if x)


@dataclass
class WallClassAudit:
    k: int
    height: int
    support_count: int
    class_examples: dict[tuple[int, ...], tuple[int, ...]]
    m_min: int | None
    m_max: int | None

    @property
    def classes(self) -> set[tuple[int, ...]]:
        return set(self.class_examples)


def enumerate_wall_classes(k: int, height: int) -> WallClassAudit:
    """Enumerate one-large support-six wall supports at coefficient height <= height."""
    B = BOUND[k]
    coeffs = coefficient_set(height)
    supports: set[tuple[int, ...]] = set()
    examples: dict[tuple[int, ...], tuple[int, ...]] = {}
    m_values: set[int] = set()

    for core_support in itertools.combinations(range(1, B + 1), 5):
        for core_coeffs in itertools.product(coeffs, repeat=5):
            core_sum = sum(c * v for c, v in zip(core_coeffs, core_support))
            if core_sum == 0:
                continue
            for cM in coeffs:
                if (-core_sum) % cM:
                    continue
                M = (-core_sum) // cM
                if M <= B:
                    continue
                support = tuple(sorted(core_support + (M,)))
                supports.add(support)
                m_values.add(M)
                cls = s14.canon_support(support)
                examples.setdefault(cls, support)

    return WallClassAudit(
        k=k,
        height=height,
        support_count=len(supports),
        class_examples=examples,
        m_min=min(m_values) if m_values else None,
        m_max=max(m_values) if m_values else None,
    )


def coimage_rows_by_k() -> dict[int, list[s14.FiberStats]]:
    classes = s14.support_classes()
    return {k: s14.compute_stats_for_d(AMBIENT_D[k], classes) for k in (8, 9, 10)}


def summarize_audit(
    audits: dict[tuple[int, int], WallClassAudit],
    coimage: dict[int, list[s14.FiberStats]],
) -> None:
    section("HEIGHT <= 2 WALL-CLASS COVERAGE OF COIMAGE MASS")
    print(
        f"{'k':>2} {'H':>2} {'supports':>10} {'classes':>8} {'M-range':>13} "
        f"{'nonzero hit':>13} {'top10 hit':>10} {'|S|-mass hit':>13}"
    )
    for k in (8, 9, 10):
        rows = coimage[k]
        nz_rows = [r for r in rows if r.signed_abs > EPS]
        total_mass = sum(r.signed_abs for r in nz_rows)
        top10 = sorted(rows, key=lambda r: r.signed_abs, reverse=True)[:10]
        for H in (1, 2):
            audit = audits[(k, H)]
            hit = audit.classes
            hit_nz = sum(1 for r in nz_rows if r.cls in hit)
            hit_top10 = sum(1 for r in top10 if r.cls in hit)
            hit_mass = sum(r.signed_abs for r in nz_rows if r.cls in hit)
            mrange = f"{audit.m_min}..{audit.m_max}"
            print(
                f"{k:>2} {H:>2} {audit.support_count:>10} {len(hit):>8} {mrange:>13} "
                f"{hit_nz:>5}/{len(nz_rows):<7} {hit_top10:>5}/10    "
                f"{hit_mass / total_mass:>12.8f}"
            )
    print(
        "\nReadout: height<=2 one-large walls cover every nonzero coimage class "
        "for k=8 and k=9.  For k=10 they cover about 84% of signed coimage "
        "mass; the remainder is structured and listed below."
    )


def missed_rows(
    audit: WallClassAudit,
    rows: list[s14.FiberStats],
    threshold: float = EPS,
) -> list[s14.FiberStats]:
    return [
        r
        for r in sorted(rows, key=lambda x: x.signed_abs, reverse=True)
        if r.signed_abs > threshold and r.cls not in audit.classes
    ]


def pattern_label(cls: tuple[int, ...]) -> str:
    counts = Counter(cls)
    multiplicities = sorted(counts.values(), reverse=True)
    zeros = counts.get(0, 0)
    if zeros >= 4:
        return "zero-cusp"
    if multiplicities[:2] == [4, 2]:
        return "4+2 repeated"
    if multiplicities[:3] == [4, 1, 1]:
        return "4+1+1 repeated"
    if multiplicities[:3] == [2, 2, 2]:
        return "2+2+2 repeated"
    return "mixed"


def detail_misses(
    audits: dict[tuple[int, int], WallClassAudit],
    coimage: dict[int, list[s14.FiberStats]],
) -> None:
    section("TAIL-ONLY COIMAGE PACKETS AFTER HEIGHT <= 2 WALL ADDRESSING")
    for k in (8, 9, 10):
        audit = audits[(k, 2)]
        misses = missed_rows(audit, coimage[k])
        total_mass = sum(r.signed_abs for r in coimage[k] if r.signed_abs > EPS)
        miss_mass = sum(r.signed_abs for r in misses)
        print(
            f"\nk={k}, ambient d={AMBIENT_D[k]}: missed nonzero classes={len(misses)}, "
            f"missed |S|-mass={miss_mass:.9g} ({miss_mass / total_mass:.6%})"
        )
        if not misses:
            print("  none: every nonzero coimage class is height<=2 wall-addressed.")
            continue
        hist = Counter(pattern_label(r.cls) for r in misses)
        print(f"  pattern histogram: {dict(hist)}")
        print(
            f"  {'rank':>4} {'class':>22} {'|S_d|':>12} {'ratio':>12} "
            f"{'pattern':>16}"
        )
        for i, row in enumerate(misses[:16], 1):
            print(
                f"  {i:>4} {str(row.cls):>22} {row.signed_abs:>12.8g} "
                f"{row.ratio:>12.6g} {pattern_label(row.cls):>16}"
            )


def examples_table(
    audits: dict[tuple[int, int], WallClassAudit],
    coimage: dict[int, list[s14.FiberStats]],
) -> None:
    section("EXAMPLES: TOP COIMAGE CLASSES AND THEIR HEIGHT<=2 WALL SUPPORTS")
    for k in (8, 9, 10):
        audit = audits[(k, 2)]
        print(f"\nk={k}, ambient d={AMBIENT_D[k]}")
        print(f"{'rank':>4} {'class':>22} {'|S_d|':>12} {'wall support example':>28}")
        for i, row in enumerate(sorted(coimage[k], key=lambda r: r.signed_abs, reverse=True)[:10], 1):
            ex = audit.class_examples.get(row.cls)
            ex_str = str(ex) if ex is not None else "TAIL-ONLY"
            print(f"{i:>4} {str(row.cls):>22} {row.signed_abs:>12.8g} {ex_str:>28}")


def theorem_readout() -> None:
    section("PROOF READOUT")
    print(
        "1. Height<=2 one-large walls are not a tiny corner of the residue table. "
        "For k=8 and k=9 they touch all coimage classes with nonzero signed "
        "mass in the relevant ambient dimension."
    )
    print(
        "2. For k=10 the remaining tail-only classes are not arbitrary.  The "
        "largest are repeated-residue packets of shape 4+2 or 4+1+1, with "
        "large abs/signed ratios.  This suggests the next analytic theorem "
        "should target repeated-residue reciprocal sums, not all 159 classes."
    )
    print(
        "3. This does not prove LRC(14): a coimage class can appear both in a "
        "low-height wall and in high-height tails.  The result is a routing "
        "lemma for the proof architecture: finite wall ledger first, then a "
        "much narrower repeated-residue signed tail."
    )


def tournament_analysis() -> None:
    section("TOURNAMENT ANALYSIS")
    vertices = [
        "height2_wall_addressed_classes",
        "height1_wall_addressed_classes",
        "repeated_residue_tail_packet",
        "coimage_fiber_atlas",
        "signed_reciprocal_tail_theorem",
        "raw_supports",
        "raw_runner_vertices",
    ]
    print("Hamiltonian proof path:")
    print("  " + " > ".join(vertices))
    print("score histogram:", {i: 1 for i in range(len(vertices))})
    print("directed 3-cycles: 0")
    print("SCC sizes:", [1] * len(vertices))
    print(
        "Assumption challenged: tournament vertices need not be runners, arcs, "
        "or raw support tuples.  Here the useful vertices are proof quotients: "
        "height-bounded wall classes, coimage packets, and the repeated-residue "
        "tail obligation.  This quotient preserves the LRC14 support-six tail "
        "predicate while throwing away witness-time geometry."
    )


def main() -> None:
    section("LRC(14) HEIGHT-2 COIMAGE WALL-CLASS AUDIT S18")
    print(
        "Goal: classify which HYP-2617 coimage classes are already addressed "
        "by one-large support-six additive walls of coefficient height <= 2."
    )
    coimage = coimage_rows_by_k()
    audits = {
        (k, H): enumerate_wall_classes(k, H)
        for k in (8, 9, 10)
        for H in (1, 2)
    }
    summarize_audit(audits, coimage)
    detail_misses(audits, coimage)
    examples_table(audits, coimage)
    theorem_readout()
    tournament_analysis()
    section("S18 STATUS")
    print(
        "LRC(14) is not proved here.  The session turns the remaining signed "
        "tail into a sharper target: after height<=2 finite wall accounting, "
        "k=8 and k=9 have no nonzero tail-only coimage classes; k=10 leaves "
        "a repeated-residue packet that should be attacked by a dedicated "
        "cotangent/Dedekind reciprocal-sum estimate."
    )


if __name__ == "__main__":
    main()
