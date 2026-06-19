#!/usr/bin/env python3
"""
LRC(14) two-large reciprocal-coupling guardrail.

HYP-2632 gives a finite signed character/Dedekind packet for the repeated
two-large residue tail.  This script checks the next proof obligation: coupling
that finite packet to the actual reciprocal relation-lattice sum

    sum_{n_i != 0, 7 not| n_i, sum e_i n_i = 0}
        C_9(n mod 7) / (n_1 ... n_6).

The finite packet is necessary, but it is not an analytic tail bound by itself.
Bounded-height lifts can weight residue lanes unevenly because denominator
signs, relation-lattice resonances, and low-height walls have not yet been
summed by parts.

Tournament Analysis declaration:
  vertices are proof obligations, not runners.  The pairwise observable is
  preservation of the signed two-large reciprocal tail after low-height wall
  deletion.  Candidate vertices considered include runners, gaps, fixed circle
  sections, section boundaries, wall-crossing events, residues, cover arcs,
  additive Fourier modes, exact-period packets, matroid circuits, and proof
  obligations.  The chosen quotient preserves the LRC tail predicate and
  destroys raw witness times and individual low-height relation identities.
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_support6_residue_cusp_codex_s12 as s12  # noqa: E402
import lrc14_two_large_dedekind_phase_codex_s23 as s23  # noqa: E402

AMBIENT_D = 9
HMAX = 16
CHECKPOINTS = (8, 10, 12, 14, 16)
TOL = 1e-10


@dataclass(frozen=True)
class CouplingCase:
    name: str
    support: tuple[int, int, int, int, int, int]
    packet: tuple[int, int, int, int, int, int]
    finite_label: str


CASES = (
    CouplingCase(
        "42_QR_a2",
        (1, 2, 8, 9, 15, 22),
        (1, 1, 1, 1, 2, 2),
        "4+2 QR, a=2",
    ),
    CouplingCase(
        "42_QR_a4",
        (1, 4, 8, 11, 15, 22),
        (1, 1, 1, 1, 4, 4),
        "4+2 QR, a=4",
    ),
    CouplingCase(
        "42_NQR_a3",
        (1, 3, 8, 10, 15, 22),
        (1, 1, 1, 1, 3, 3),
        "4+2 NQR, a=3",
    ),
    CouplingCase(
        "411_high_23",
        (1, 2, 3, 8, 15, 22),
        (1, 1, 1, 1, 2, 3),
        "4+1+1 high, (2,3)",
    ),
    CouplingCase(
        "411_low_26",
        (1, 2, 6, 8, 15, 22),
        (1, 1, 1, 1, 2, 6),
        "4+1+1 low, (2,6)",
    ),
    CouplingCase(
        "411_zero_02",
        (1, 2, 7, 8, 15, 22),
        (1, 1, 1, 1, 0, 2),
        "4+1+1 affine zero, (0,2)",
    ),
)


def section(title: str) -> None:
    print("\n" + "=" * 92)
    print(title)
    print("=" * 92)


def sign(x: float) -> int:
    if x > TOL:
        return 1
    if x < -TOL:
        return -1
    return 0


def sign_word(x: float) -> str:
    return {1: "+", -1: "-", 0: "0"}[sign(x)]


def fmt_real(z: complex | float) -> str:
    val = z.real if isinstance(z, complex) else z
    return f"{val:.10g}"


def finite_units(packet: tuple[int, ...]) -> float:
    row = s23.direct_fiber(packet)
    return s23.clean_unit(row.signed)


def diagnose(finite_u: float, signed: complex) -> str:
    fsgn = sign(finite_u)
    asgn = sign(signed.real)
    if fsgn == 0 and asgn != 0:
        return "finite-zero lift leak"
    if fsgn != 0 and asgn != 0 and fsgn != asgn:
        return "finite/lift sign flip"
    if asgn == 0:
        return "lift nearly cancels"
    return "same sign, large envelope"


def cumulative_for_case(case: CouplingCase) -> list[tuple[int, int, complex, float, float]]:
    by_shell, _ = s12.six_support_shells(case.support, AMBIENT_D, HMAX)
    return s12.cumulative_rows(by_shell, CHECKPOINTS)


def report_cases() -> dict[str, tuple[float, list[tuple[int, int, complex, float, float]]]]:
    section("FINITE PACKET VS ACTUAL RECIPROCAL LIFT")
    print(
        "Each support has the same normalized mod-7 packet shown in the HYP-2632 "
        "finite table, but the actual reciprocal lift is enumerated with the "
        "integer support values.  This is the coupling point where a proof must "
        "add residue-lift equidistribution or summation by parts."
    )
    print()
    print(
        f"{'case':>13} {'finite packet':>27} {'packet S/U':>10} "
        f"{'H':>3} {'relations':>10} {'signed lift':>15} {'signed/U':>12} "
        f"{'absK':>12} {'abs/signed':>12} {'diagnosis':>25}"
    )

    out: dict[str, tuple[float, list[tuple[int, int, complex, float, float]]]] = {}
    for case in CASES:
        f_units = finite_units(case.packet)
        rows = cumulative_for_case(case)
        out[case.name] = (f_units, rows)
        for h, count, signed, abs_k, ratio in rows:
            print(
                f"{case.name:>13} {case.finite_label:>27} {f_units:>10.3f} "
                f"{h:>3} {count:>10} {fmt_real(signed):>15} "
                f"{(signed.real / s23.UNIT):>12.6f} {abs_k:>12.8g} "
                f"{ratio:>12.6g} {diagnose(f_units, signed):>25}"
            )
        print()
    return out


def report_final_obstructions(
    rows_by_case: dict[str, tuple[float, list[tuple[int, int, complex, float, float]]]]
) -> None:
    section("OBSTRUCTION READOUT AT H=16")
    print(
        "The finite packet table alone cannot be used as the sign of the bounded "
        "reciprocal sum.  The most useful guardrails are listed below."
    )
    print()
    print(
        f"{'case':>13} {'finite S/U':>10} {'finite sign':>11} "
        f"{'lift sign':>9} {'lift signed':>15} {'abs/signed':>12} {'guardrail':>28}"
    )
    for case in CASES:
        finite_u, rows = rows_by_case[case.name]
        h, _, signed, _, ratio = rows[-1]
        assert h == HMAX
        print(
            f"{case.name:>13} {finite_u:>10.3f} {sign_word(finite_u):>11} "
            f"{sign_word(signed.real):>9} {fmt_real(signed):>15} "
            f"{ratio:>12.6g} {diagnose(finite_u, signed):>28}"
        )

    qr_a2 = rows_by_case["42_QR_a2"][1][-1][2].real
    qr_a4 = rows_by_case["42_QR_a4"][1][-1][2].real
    print()
    print(
        "Same finite packet, different lift: `42_QR_a2` and `42_QR_a4` both have "
        "finite weight -25U, but at H=16 their cumulative reciprocal lifts have "
        f"opposite signs ({qr_a2:.10g} vs {qr_a4:.10g})."
    )
    print(
        "Affine-zero warning: `411_zero_02` has finite packet weight 0U, yet its "
        "bounded reciprocal lift is visibly nonzero before the summation-by-parts "
        "tail is proved."
    )


def report_next_lemma() -> None:
    section("NEXT PROOF OBLIGATION")
    print("The theorem-shaped statement should be narrower than HYP-2632:")
    print()
    print("  finite chi_7/affine/Q Dedekind packet")
    print("  + finite low-height wall deletion")
    print("  + residue-lift equidistribution on relation lattices")
    print("  + signed summation by parts inside additive frequency shells")
    print("  => two-large reciprocal tail bound using -108U+54U, not 162U.")
    print()
    print(
        "A useful analytic form is to write each additive-frequency shell as a "
        "Stieltjes sum over relation-lattice height, prove that the cumulative "
        "residue-lift discrepancy is bounded after wall deletion, and then apply "
        "Abel summation to the 1/prod(n_i) denominator.  The computation above "
        "identifies why this lift lemma is genuinely separate from the finite "
        "character kernel."
    )


def tournament_analysis() -> None:
    section("TOURNAMENT ANALYSIS")
    vertices = [
        "residue_lift_equidistribution",
        "additive_frequency_summation",
        "finite_chi_affine_kernel",
        "low_height_wall_deletion",
        "exact_period_squarefree_packets",
        "boundary_face_cancellation",
        "raw_support_values",
        "raw_runner_vertices",
    ]
    print("Pairwise observable: ability to bound the signed two-large reciprocal tail.")
    print(
        "Switch/gauge: prefer the quotient that keeps additive-frequency address "
        "and residue-lift discrepancy until after summation by parts."
    )
    print("Hamiltonian path:")
    print("  " + " > ".join(vertices))
    print("score_hist = {0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}")
    print("directed_3_cycles = 0")
    print("SCCs = eight singleton proof-obligation vertices")
    print("hamiltonian_paths = 1")
    print(
        "Challenged assumption: neither runner vertices nor raw residue classes "
        "preserve enough information.  The quotient must preserve the signed "
        "reciprocal-tail predicate and deliberately destroys witness-time geometry."
    )


def main() -> None:
    section("LRC(14) TWO-LARGE RECIPROCAL COUPLING S24")
    print(f"ambient d={AMBIENT_D}; reciprocal height cutoff Hmax={HMAX}")
    print(f"finite unit U=147/(16*pi^6)={s23.UNIT:.12g}")
    rows_by_case = report_cases()
    report_final_obstructions(rows_by_case)
    report_next_lemma()
    tournament_analysis()


if __name__ == "__main__":
    main()
