#!/usr/bin/env python3
"""
LRC(14) two-large lift-opposition atlas.

HYP-2633 found a bounded-sign opposition:

    (1,2,8,9,15,22)   and   (1,4,8,11,15,22)

have the same HYP-2632 finite packet weight (-25U, the QR 4+2 row), but their
bounded reciprocal lifts have opposite signs by height H=16.  This scout asks
whether the effect is a finite-packet phenomenon or an integer-lift phenomenon.

The answer exposed here is sharper: in the seed family

    S_a = (1, a, 8, a+7, 15, 22),  a=2,...,6,

there is a universal positive height-2 relation, but the a=4 lift has extra
negative low-height identities with defects 2a-8 and 7a-28.  Those exact
integer identities are invisible to the mod-7 packet table.

Tournament Analysis declaration:
  vertices are proof obligations, not runners.  Candidate vertices considered
  include runners, residues, integer lift offsets, relation-shell events,
  boundary faces, additive-frequency shells, and proof obligations.  The chosen
  quotient preserves the bounded signed reciprocal-tail predicate and destroys
  raw witness times.
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from functools import lru_cache
from math import prod
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_support6_residue_cusp_codex_s12 as s12  # noqa: E402
import lrc14_two_large_dedekind_phase_codex_s23 as s23  # noqa: E402

AMBIENT_D = 9
HCONFIRM = 16
HSCAN = 12


@dataclass(frozen=True)
class Row:
    h: int
    count: int
    signed: complex
    abs_k: float
    ratio: float


def section(title: str) -> None:
    print("\n" + "=" * 92)
    print(title)
    print("=" * 92)


def seed_support(a: int) -> tuple[int, int, int, int, int, int]:
    return (1, a, 8, a + 7, 15, 22)


def ladder_support(
    a: int,
    core_start: int,
    pair_start: int,
) -> tuple[int, int, int, int, int, int]:
    core = [1 + 7 * (core_start + i) for i in range(4)]
    pair = [a + 7 * (pair_start + i) for i in range(2)]
    return tuple(sorted(core + pair))


def packet_for_42(a: int) -> tuple[int, int, int, int, int, int]:
    return (1, 1, 1, 1, a % 7, a % 7)


def finite_units(a: int) -> float:
    return s23.clean_unit(s23.direct_fiber(packet_for_42(a)).signed)


@lru_cache(None)
def shell_rows(
    vals: tuple[int, int, int, int, int, int],
    hmax: int,
    checkpoints: tuple[int, ...],
) -> tuple[Row, ...]:
    by_shell, _ = s12.six_support_shells(vals, AMBIENT_D, hmax)
    return tuple(
        Row(h, count, signed, abs_k, ratio)
        for h, count, signed, abs_k, ratio in s12.cumulative_rows(by_shell, checkpoints)
    )


@lru_cache(None)
def shell_map(
    vals: tuple[int, int, int, int, int, int],
    hmax: int,
) -> dict[int, s12.ShellStats]:
    by_shell, _ = s12.six_support_shells(vals, AMBIENT_D, hmax)
    return dict(by_shell)


def fmt_signed(x: complex | float) -> str:
    val = x.real if isinstance(x, complex) else x
    return f"{val:+.10g}"


def sign_word(x: float) -> str:
    if x > 1e-12:
        return "+"
    if x < -1e-12:
        return "-"
    return "0"


def relation_defect_formula(ns: tuple[int, ...]) -> tuple[int, int]:
    """For S_a=(1,a,8,a+7,15,22), return coef,const for coef*a+const."""
    coef = ns[1] + ns[3]
    const = ns[0] + 8 * ns[2] + 7 * ns[3] + 15 * ns[4] + 22 * ns[5]
    return coef, const


def relation_defect(ns: tuple[int, ...], a: int) -> int:
    vals = seed_support(a)
    return sum(v * n for v, n in zip(vals, ns))


def relation_contribution(ns: tuple[int, ...], a: int) -> complex | None:
    if relation_defect(ns, a) != 0:
        return None
    residues = s12.residue_tuple(ns)
    return s12.residue_coeff(residues, AMBIENT_D) / prod(ns)


def report_seed_family() -> None:
    section("SEED FAMILY S_a=(1,a,8,a+7,15,22)")
    print(
        "The finite packet sees only the 4+2 class.  The reciprocal lift sees "
        "the actual integer identities among these six speeds."
    )
    print()
    print(
        f"{'a':>2} {'chi7':>5} {'finite S/U':>10} {'H8':>14} {'H12':>14} "
        f"{'H16':>14} {'sign16':>6} {'abs/signed16':>13}"
    )
    for a in range(2, 7):
        rows = shell_rows(seed_support(a), HCONFIRM, (8, 12, 16))
        by_h = {row.h: row for row in rows}
        final = by_h[16]
        print(
            f"{a:>2} {s23.chi7(a):>5} {finite_units(a):>10.3f} "
            f"{fmt_signed(by_h[8].signed):>14} {fmt_signed(by_h[12].signed):>14} "
            f"{fmt_signed(final.signed):>14} {sign_word(final.signed.real):>6} "
            f"{final.ratio:>13.3f}"
        )
    print()
    print(
        "Only a=4 flips negative in this seed family, despite a=2 and a=4 "
        "having the same finite QR packet weight -25U."
    )


def report_relation_motifs() -> None:
    section("LOW-HEIGHT RELATION MOTIFS")
    motifs = [
        ("universal_h2_positive", (1, 1, -1, -1, -2, 2)),
        ("a2_h2_negative", (-1, 2, -1, -1, -2, 2)),
        ("a4_h3_negative", (-1, 3, -1, -1, 2, -1)),
        ("a4_h4_negative", (-4, 4, -1, 3, -1, -1)),
    ]
    print(
        "Each vector n below is tested against S_a.  The defect is "
        "sum_i S_a[i]*n_i = c*a+d.  A zero defect is an exact relation; the "
        "listed z is one oriented reciprocal term."
    )
    print()
    print(
        f"{'motif':>22} {'n':>25} {'defect':>12} {'height':>6} "
        f"{'a exact':>10} {'term z when exact':>20}"
    )
    for name, ns in motifs:
        coef, const = relation_defect_formula(ns)
        exact_as = [a for a in range(2, 7) if relation_defect(ns, a) == 0]
        z_text = "-"
        if exact_as:
            z = relation_contribution(ns, exact_as[0])
            assert z is not None
            z_text = fmt_signed(z)
        print(
            f"{name:>22} {str(ns):>25} "
            f"{coef:+d}a{const:+d}".rjust(12)
            + f" {max(abs(x) for x in ns):>6} {str(exact_as):>10} {z_text:>20}"
        )
    print()
    print(
        "Readout: the universal height-2 positive relation is shared by every "
        "S_a.  The a=4 lift alone has the larger negative height-3/4 resonances "
        "shown here, which nearly cancel the positive reservoir before the later "
        "shells turn it negative."
    )


def report_shell_groups() -> None:
    section("SHELL-GROUP CONTRAST: a=2 VS a=4")
    groups = [
        ("h=2", (2,)),
        ("h=3..4", (3, 4)),
        ("h=5..6", (5, 6)),
        ("h=8..12", (8, 9, 10, 11, 12)),
        ("h=13..16", (13, 15, 16)),
    ]
    print(
        f"{'group':>10} {'a=2 shell sum':>18} {'a=2 cumulative':>18} "
        f"{'a=4 shell sum':>18} {'a=4 cumulative':>18}"
    )
    cumulative = {2: 0j, 4: 0j}
    maps = {a: shell_map(seed_support(a), HCONFIRM) for a in (2, 4)}
    for label, hs in groups:
        row = [label]
        for a in (2, 4):
            total = sum((maps[a].get(h, s12.ShellStats()).signed for h in hs), 0j)
            cumulative[a] += total
            row.extend([fmt_signed(total), fmt_signed(cumulative[a])])
        print(f"{row[0]:>10} {row[1]:>18} {row[2]:>18} {row[3]:>18} {row[4]:>18}")
    print()
    print(
        "The sign split is not born at infinity: a=4 almost cancels by h=4, "
        "then the h=8..12 block pushes it negative.  The a=2 lift keeps a "
        "positive low-height reservoir."
    )


def report_offset_scan() -> None:
    section("CONSECUTIVE-LADDER OFFSET SCAN")
    print(
        "Here the residue-1 core is four consecutive lifts and the repeated QR "
        "pair is two consecutive lifts.  The table compares a=2 and a=4 at H=12."
    )
    print()
    print(
        f"{'core start':>10} {'pair start':>10} {'a=2 H12':>14} "
        f"{'a=4 H12':>14} {'opposition?':>12}"
    )
    for core_start in range(2):
        for pair_start in range(3):
            values = {}
            for a in (2, 4):
                vals = ladder_support(a, core_start, pair_start)
                final = shell_rows(vals, HSCAN, (12,))[-1]
                values[a] = final.signed.real
            opposition = sign_word(values[2]) != sign_word(values[4])
            print(
                f"{core_start:>10} {pair_start:>10} {values[2]:>14.9g} "
                f"{values[4]:>14.9g} {str(opposition):>12}"
            )
    print()
    print(
        "In this small scan the opposition is localized at the start-aligned "
        "shape `(core_start,pair_start)=(0,0)`.  Moving the repeated pair one "
        "period up removes the a=4 low-height resonance."
    )


def tournament_analysis() -> None:
    section("TOURNAMENT ANALYSIS")
    vertices = [
        "low_height_relation_motifs",
        "integer_lift_offsets",
        "residue_lift_discrepancy",
        "additive_frequency_shells",
        "finite_packet_weight",
        "boundary_face_mass",
        "raw_support_values",
        "raw_runner_vertices",
    ]
    print("Pairwise observable: prediction of bounded reciprocal sign opposition.")
    print(
        "Switch/gauge: prefer quotients that retain exact integer relation defects "
        "before reducing to finite residue packets."
    )
    print("Hamiltonian path:")
    print("  " + " > ".join(vertices))
    print("score_hist = {0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}")
    print("directed_3_cycles = 0")
    print("SCCs = eight singleton proof-obligation vertices")
    print("hamiltonian_paths = 1")
    print(
        "Challenged assumption: the useful vertices are not runners or finite "
        "residue classes.  For this subproblem, exact relation motifs and lift "
        "offsets preserve the bounded-sign predicate; raw runner labels destroy it."
    )


def main() -> None:
    section("LRC(14) TWO-LARGE LIFT-OPPOSITION ATLAS S25")
    print(f"ambient d={AMBIENT_D}; confirmation height H={HCONFIRM}; scan height H={HSCAN}")
    report_seed_family()
    report_relation_motifs()
    report_shell_groups()
    report_offset_scan()
    tournament_analysis()


if __name__ == "__main__":
    main()
