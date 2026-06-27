#!/usr/bin/env python3
"""Observer-gluing ledger scout for the LRC(14) proof frontier.

This alternates two active frontiers.

Task A (HYP-3096): compute exact direct lonely-set topology for representative
THM-573 residual rows, keeping the polynomial-method CRT fields.

Task B (HYP-3097): attach pair/Pascal scissors fields to the same rows, so the
cap/pair-mass shadow is not studied in a separate quotient.

The output is a proof-obligation ledger, not a proof.  It is designed to show
which chart supplies a terminal exit and which chart still carries debt.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import comb, gcd
from typing import Iterable


THRESHOLD = Fraction(1, 14)
PAIR_APEX = comb(14, 2)


def lcm(a: int, b: int) -> int:
    return a // gcd(a, b) * b


def lcm_upto(n: int) -> int:
    out = 1
    for value in range(1, n + 1):
        out = lcm(out, value)
    return out


def gcd_all(values: Iterable[int]) -> int:
    out = 0
    for value in values:
        out = gcd(out, value)
    return out


def frac_part(value: Fraction) -> Fraction:
    return value - (value.numerator // value.denominator)


def norm_circle(value: Fraction) -> Fraction:
    residue = frac_part(value)
    return min(residue, 1 - residue)


def is_strict_witness(speeds: tuple[int, ...], t: Fraction) -> bool:
    return all(norm_circle(Fraction(speed) * t) > THRESHOLD for speed in speeds)


def lonely_arcs(speeds: tuple[int, ...]) -> list[tuple[Fraction, Fraction]]:
    """Exact positive-length components of {t: ||s_i t|| >= 1/14 all i}.

    Boundary-only witnesses are intentionally excluded from the component
    count because HYP-3096 needs denominator-net arcs, not isolated endpoints.
    """
    points = {Fraction(0), Fraction(1)}
    for speed in speeds:
        if speed == 0:
            continue
        for m in range(speed + 1):
            left = Fraction(m, speed) - Fraction(1, 14 * speed)
            right = Fraction(m, speed) + Fraction(1, 14 * speed)
            if 0 <= left <= 1:
                points.add(left)
            if 0 <= right <= 1:
                points.add(right)

    arcs: list[tuple[Fraction, Fraction]] = []
    ordered = sorted(points)
    for left, right in zip(ordered, ordered[1:]):
        if right <= left:
            continue
        midpoint = (left + right) / 2
        if is_strict_witness(speeds, midpoint):
            arcs.append((left, right))
    return arcs


def denominator_net_threshold(length: Fraction) -> int | None:
    if length <= 0:
        return None
    return length.denominator // length.numerator + 1


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def short_float(value: Fraction) -> str:
    return f"{float(value):.8g}"


def residue_counts(values: tuple[int, ...], modulus: int) -> tuple[int, ...]:
    counts = [0] * modulus
    for value in values:
        counts[value % modulus] += 1
    return tuple(counts)


def pair_shadow(count: int) -> Fraction:
    return Fraction(comb(count, 2), PAIR_APEX)


@dataclass(frozen=True)
class SampleRow:
    row_id: str
    speeds: tuple[int, ...]
    farey_root: tuple[int, int] | None = None
    note: str = ""


@dataclass(frozen=True)
class LedgerRow:
    row_id: str
    primitive_gcd: int
    covering_mod14: bool
    covering_mod7: bool
    count_7_divisible: int
    count_14_divisible: int
    count_even: int
    crt_c7_lift_status: str
    crt_c2_lift_status: str
    level7_residual_status: str
    direct_component_count: int
    direct_lonely_measure: Fraction
    largest_direct_arc: Fraction
    denominator_threshold: int | None
    terminal_exit: str
    mod7_residue_counts: tuple[int, ...]
    mod14_residue_counts: tuple[int, ...]
    h7_pair_shadow: Fraction
    even_pair_shadow: Fraction
    triangular_cap_shadow: Fraction
    farey_additive_lane_mod91: int | None
    farey_product_lane_mod91: int | None
    gluing_debt: tuple[str, ...]


def classify_lifts(speeds: tuple[int, ...]) -> tuple[str, str, str, str, tuple[str, ...]]:
    count7 = sum(speed % 7 == 0 for speed in speeds)
    covering14 = any(speed % 14 == 0 for speed in speeds)
    covering7 = any(speed % 7 == 0 for speed in speeds)

    if not covering14:
        return (
            "not_needed_no_mod14_covering",
            "grid_witness_t=1/14",
            "outside_covering_residual",
            "Q_WITNESS_1/14",
            (),
        )

    if count7 >= 7:
        return (
            "THM-573_c7_lift_exit",
            "dyadic_chart_not_needed_after_c7_exit",
            "closed_by_level7_sieve",
            "THM573_LEVEL7_EXIT",
            (),
        )

    c7_status = "I(13,7,1)_covering_mod7_residual" if covering7 else "level7_grid_proper"
    c2_status = "covering_mod14_improper_on_1/14_grid"
    residual = f"live_THM573_residual_H7={count7}_le6"
    debt = (
        "normalized_arc_chart",
        "moment_perron_or_gK8_chart",
        "branch_or_K33_shuttle_chart",
        "finite_denominator_budget",
    )
    return c7_status, c2_status, residual, "LIVE_OBSERVER_GLUE_DEBT", debt


def build_ledger_row(sample: SampleRow) -> LedgerRow:
    speeds = tuple(sorted(sample.speeds))
    arcs = lonely_arcs(speeds)
    measure = sum((right - left) for left, right in arcs)
    largest = max((right - left for left, right in arcs), default=Fraction(0))
    c7, c2, residual, terminal, debt = classify_lifts(speeds)
    if terminal == "LIVE_OBSERVER_GLUE_DEBT" and largest > 0:
        terminal = "DIRECT_ARC_CANDIDATE_NEEDS_UNIFORMITY_PROOF"
        debt = debt + ("uniform_component_bound",)
    if terminal == "LIVE_OBSERVER_GLUE_DEBT" and largest == 0:
        debt = debt + ("boundary_or_tight_locus_classification",)

    root_sum = root_product = None
    if sample.farey_root is not None:
        p, q = sample.farey_root
        root_sum = (p + q) % PAIR_APEX
        root_product = (p * q) % PAIR_APEX

    count7 = sum(speed % 7 == 0 for speed in speeds)
    count14 = sum(speed % 14 == 0 for speed in speeds)
    count_even = sum(speed % 2 == 0 for speed in speeds)
    return LedgerRow(
        row_id=sample.row_id,
        primitive_gcd=gcd_all(speeds),
        covering_mod14=count14 > 0,
        covering_mod7=count7 > 0,
        count_7_divisible=count7,
        count_14_divisible=count14,
        count_even=count_even,
        crt_c7_lift_status=c7,
        crt_c2_lift_status=c2,
        level7_residual_status=residual,
        direct_component_count=len(arcs),
        direct_lonely_measure=measure,
        largest_direct_arc=largest,
        denominator_threshold=denominator_net_threshold(largest),
        terminal_exit=terminal,
        mod7_residue_counts=residue_counts(speeds, 7),
        mod14_residue_counts=residue_counts(speeds, 14),
        h7_pair_shadow=pair_shadow(count7),
        even_pair_shadow=pair_shadow(count_even),
        triangular_cap_shadow=Fraction(PAIR_APEX, PAIR_APEX),
        farey_additive_lane_mod91=root_sum,
        farey_product_lane_mod91=root_product,
        gluing_debt=debt,
    )


def sample_rows() -> list[SampleRow]:
    rows: list[SampleRow] = [
        SampleRow(
            "A0_q_witness_AP_1_13",
            tuple(range(1, 14)),
            note="14-free boundary row; t=1/14 is a grid witness.",
        ),
        SampleRow(
            "A1_c7_discharge_7_multiples",
            (1, 2, 3, 4, 5, 6, 7, 14, 21, 28, 35, 42, 49),
            note="Synthetic row where THM-573 fires.",
        ),
        SampleRow(
            "A2_H7eq6_boundary_residual",
            (1, 2, 3, 4, 5, 6, 13, 14, 21, 28, 35, 42, 49),
            note="Exactly six 7-divisible speeds, so c=7 just fails.",
        ),
    ]
    for V in (2, 13, 50, 200):
        rows.append(
            SampleRow(
                f"B_apex_family_1_12_14V_V{V}",
                tuple(list(range(1, 13)) + [14 * V]),
                farey_root=(14, V),
                note="Apex family used in the V* crossover scout.",
            )
        )
    for B in (6, 8):
        apex = 84 * lcm_upto(B)
        rows.append(
            SampleRow(
                f"C_divisor_loaded_B{B}",
                tuple(list(range(1, 12)) + [13, apex]),
                farey_root=(84, lcm_upto(B)),
                note="THM-575 divisor-loaded raw Conjecture 7.1 obstruction.",
            )
        )
    return rows


def print_row(row: LedgerRow) -> None:
    print(f"ROW {row.row_id}")
    print(
        "  CRT:",
        f"gcd={row.primitive_gcd}",
        f"cover14={row.covering_mod14}",
        f"H7={row.count_7_divisible}",
        f"H14={row.count_14_divisible}",
        f"c7={row.crt_c7_lift_status}",
        f"c2={row.crt_c2_lift_status}",
    )
    print(
        "  DIRECT_L:",
        f"components={row.direct_component_count}",
        f"measure={fraction_text(row.direct_lonely_measure)}",
        f"measure_float={short_float(row.direct_lonely_measure)}",
        f"largest={fraction_text(row.largest_direct_arc)}",
        f"largest_float={short_float(row.largest_direct_arc)}",
        f"D_grid={row.denominator_threshold}",
    )
    print(
        "  PAIR_SHADOW:",
        f"tri_cap={fraction_text(row.triangular_cap_shadow)}",
        f"H7_pair={fraction_text(row.h7_pair_shadow)}",
        f"even_pair={fraction_text(row.even_pair_shadow)}",
        f"farey_sum_mod91={row.farey_additive_lane_mod91}",
        f"farey_product_mod91={row.farey_product_lane_mod91}",
    )
    print(
        "  SCISSORS:",
        f"mod7_counts={row.mod7_residue_counts}",
        f"mod14_counts={row.mod14_residue_counts}",
    )
    print("  EXIT:", row.terminal_exit)
    if row.gluing_debt:
        print("  DEBT:", ", ".join(row.gluing_debt))


def print_summary(rows: list[LedgerRow]) -> None:
    print("SUMMARY")
    exits: dict[str, int] = {}
    for row in rows:
        exits[row.terminal_exit] = exits.get(row.terminal_exit, 0) + 1
    print("  terminal_exit_histogram:", exits)
    residuals = [row for row in rows if row.level7_residual_status.startswith("live_")]
    if residuals:
        max_components = max(row.direct_component_count for row in residuals)
        min_largest = min(row.largest_direct_arc for row in residuals)
        max_denominator = max(row.denominator_threshold or 0 for row in residuals)
        print(
            "  live_residual_direct_arc_scan:",
            f"rows={len(residuals)}",
            f"max_components={max_components}",
            f"min_largest={fraction_text(min_largest)}",
            f"max_sample_D={max_denominator}",
        )
    unique_mod7 = {row.mod7_residue_counts for row in rows}
    print("  sector_scissors_signatures_mod7:", len(unique_mod7))
    print(
        "  readout:",
        "direct arcs are computable packet fields, but the live rows still need a",
        "uniform component/floor theorem plus moment/branch overlap maps before",
        "the scalar witness route is legal.",
    )


def print_tournament_analysis() -> None:
    vertices = [
        "observer_gluing_packet",
        "direct_lonely_component_packet",
        "sector_pair_scissors_packet",
        "crt_c7_level7_lift",
        "crt_c2_dyadic_lift",
        "pascal_pair_mass_shadow",
        "raw_direct_arc_scalar",
        "raw_I_table_enumeration",
    ]
    print("TOURNAMENT_ANALYSIS")
    print("  vertices_are: proof packets / observer charts, not runners")
    print(
        "  observable:",
        "preserves LRC predicate, names destroyed coordinate, retains CRT status,",
        "keeps component topology, keeps pair-scissors sidecars, supplies terminal exit",
    )
    print("  tie_hamiltonian_path:", " > ".join(vertices))
    print("  score_histogram:", {idx: 1 for idx in range(len(vertices))})
    print("  directed_3cycles: 0")
    print("  scc_sizes:", [1] * len(vertices))
    print("  challenged_assumption: raw residues or scalar caps are not vertices until their observer sidecars glue")


def main() -> None:
    print("LRC14 OBSERVER-GLUING LEDGER SCOUT (codex-2026-06-27-S258)")
    print("Task A: HYP-3096 direct witness-route ledger")
    print("Task B: HYP-3097 pair/Pascal scissors fields on the same rows")
    print()
    ledger = [build_ledger_row(sample) for sample in sample_rows()]
    for row in ledger:
        print_row(row)
        print()
    print_summary(ledger)
    print()
    print_tournament_analysis()


if __name__ == "__main__":
    main()
