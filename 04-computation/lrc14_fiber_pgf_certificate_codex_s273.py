#!/usr/bin/env python3
"""HYP-3140 scout: a 14-sheet fiber PGF certificate for the LRC14 Rprime floor.

This is an exact-rational generating-function reframe of the HYP-3136
multi-far factorization, specializing its remaining HYP-3129/HYP-3125
Rprime factor.  For S = R union 14Q, write u = 14t and count

    N_R(u) = #{a in {0,...,13}: (u+a)/14 is R-safe}.

Then Q-lonely is a function of u alone, and

    Rprime = meas(R-safe and Q-lonely)/(meas(R-safe)*meas(Q-lonely))
           = E[N_R(u) | Q-lonely] / E[N_R(u)].

Thus the SPEC covariance can be studied as a coefficient/conditional-mean
problem for the sheet-count PGF, rather than as a free-floating scalar.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations
from math import comb


F = Fraction


def merge(intervals: list[tuple[F, F]]) -> list[tuple[F, F]]:
    intervals = sorted((a, b) for a, b in intervals if a < b)
    out: list[tuple[F, F]] = []
    for a, b in intervals:
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out


def intersect(a_intervals: list[tuple[F, F]], b_intervals: list[tuple[F, F]]) -> list[tuple[F, F]]:
    a_intervals = merge(a_intervals)
    b_intervals = merge(b_intervals)
    out: list[tuple[F, F]] = []
    i = j = 0
    while i < len(a_intervals) and j < len(b_intervals):
        lo = max(a_intervals[i][0], b_intervals[j][0])
        hi = min(a_intervals[i][1], b_intervals[j][1])
        if lo < hi:
            out.append((lo, hi))
        if a_intervals[i][1] < b_intervals[j][1]:
            i += 1
        else:
            j += 1
    return out


def measure(intervals: list[tuple[F, F]]) -> F:
    return sum((b - a for a, b in intervals), F(0))


def safe_arcs(speed: int, threshold: F = F(1, 14)) -> list[tuple[F, F]]:
    """Return {t in [0,1): ||speed*t|| >= threshold}."""
    arcs = []
    for j in range(speed):
        lo = F(j, speed) + threshold / speed
        hi = F(j + 1, speed) - threshold / speed
        if lo < hi:
            arcs.append((lo, hi))
    return merge(arcs)


def safe_set(speeds: tuple[int, ...], threshold: F = F(1, 14)) -> list[tuple[F, F]]:
    arcs = [(F(0), F(1))]
    for speed in speeds:
        arcs = intersect(arcs, safe_arcs(speed, threshold))
        if not arcs:
            break
    return arcs


def sheet_safe_arcs(r_speeds: tuple[int, ...], sheet: int) -> list[tuple[F, F]]:
    """R-safe intervals in u for the sheet t=(u+sheet)/14."""
    t_arcs = [(F(sheet, 14), F(sheet + 1, 14))]
    for speed in r_speeds:
        t_arcs = intersect(t_arcs, safe_arcs(speed))
        if not t_arcs:
            break
    return merge([(14 * lo - sheet, 14 * hi - sheet) for lo, hi in t_arcs])


def contains(intervals: list[tuple[F, F]], point: F) -> bool:
    return any(lo <= point < hi for lo, hi in intervals)


def direct_covering_metrics(r_speeds: tuple[int, ...], q_speeds: tuple[int, ...]) -> tuple[F, F, F, F]:
    r_arcs = safe_set(r_speeds)
    q_lift_arcs = safe_set(tuple(14 * q for q in q_speeds))
    both = intersect(r_arcs, q_lift_arcs)
    r_meas = measure(r_arcs)
    q_meas = measure(q_lift_arcs)
    both_meas = measure(both)
    rprime = both_meas / (r_meas * q_meas)
    return r_meas, q_meas, both_meas, rprime


def fiber_pgf_metrics(r_speeds: tuple[int, ...], q_speeds: tuple[int, ...]) -> dict[str, object]:
    sheets = [sheet_safe_arcs(r_speeds, sheet) for sheet in range(14)]
    q_arcs = safe_set(q_speeds)

    breakpoints = {F(0), F(1)}
    for arcs in sheets + [q_arcs]:
        for lo, hi in arcs:
            breakpoints.add(lo)
            breakpoints.add(hi)
    points = sorted(breakpoints)

    pgf: dict[int, F] = defaultdict(F)
    q_pgf: dict[int, F] = defaultdict(F)
    total_sheet_mass = F(0)
    q_sheet_mass = F(0)
    q_mass = F(0)
    q_zero_mass = F(0)
    min_q_count = 15
    max_q_count = -1

    for lo, hi in zip(points, points[1:]):
        if lo >= hi:
            continue
        midpoint = (lo + hi) / 2
        width = hi - lo
        sheet_count = sum(1 for arcs in sheets if contains(arcs, midpoint))
        q_ok = contains(q_arcs, midpoint)
        pgf[sheet_count] += width
        total_sheet_mass += sheet_count * width
        if q_ok:
            q_pgf[sheet_count] += width
            q_mass += width
            q_sheet_mass += sheet_count * width
            min_q_count = min(min_q_count, sheet_count)
            max_q_count = max(max_q_count, sheet_count)
            if sheet_count == 0:
                q_zero_mass += width

    mean_sheets = total_sheet_mass
    q_mean_sheets = q_sheet_mass / q_mass
    rprime_fiber = q_mean_sheets / mean_sheets
    direct_r, direct_q, direct_both, direct_rprime = direct_covering_metrics(r_speeds, q_speeds)

    return {
        "pgf": dict(sorted(pgf.items())),
        "q_pgf": dict(sorted(q_pgf.items())),
        "mean_sheets": mean_sheets,
        "q_mean_sheets": q_mean_sheets,
        "q_mass": q_mass,
        "q_zero_mass": q_zero_mass,
        "min_q_count": min_q_count,
        "max_q_count": max_q_count,
        "rprime_fiber": rprime_fiber,
        "direct_r": direct_r,
        "direct_q": direct_q,
        "direct_both": direct_both,
        "direct_rprime": direct_rprime,
        "identity_ok": direct_r == mean_sheets / 14
        and direct_q == q_mass
        and direct_both == q_sheet_mass / 14
        and direct_rprime == rprime_fiber,
    }


def fmt_fraction(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def fmt_pgf(pgf: dict[int, F]) -> str:
    return " + ".join(f"{fmt_fraction(weight)}*y^{degree}" for degree, weight in pgf.items() if weight)


def print_case(label: str, r_speeds: tuple[int, ...], q_speeds: tuple[int, ...]) -> dict[str, object]:
    metrics = fiber_pgf_metrics(r_speeds, q_speeds)
    print("=" * 100)
    print(label)
    print(f"R={r_speeds}")
    print(f"Q={q_speeds}   far speeds={tuple(14 * q for q in q_speeds)}")
    print("-" * 100)
    print(f"F_R(y)        = {fmt_pgf(metrics['pgf'])}")
    print(f"F_R,Q(y)      = {fmt_pgf(metrics['q_pgf'])}")
    print(f"E[N_R]        = {fmt_fraction(metrics['mean_sheets'])} = {float(metrics['mean_sheets']):.8f}")
    print(f"E[N_R | Q]    = {fmt_fraction(metrics['q_mean_sheets'])} = {float(metrics['q_mean_sheets']):.8f}")
    print(f"Q mass        = {fmt_fraction(metrics['q_mass'])} = {float(metrics['q_mass']):.8f}")
    print(
        "Q sheet range = "
        f"{metrics['min_q_count']}..{metrics['max_q_count']}; "
        f"Q zero sheet mass={fmt_fraction(metrics['q_zero_mass'])}"
    )
    print(f"Rprime fiber  = {fmt_fraction(metrics['rprime_fiber'])} = {float(metrics['rprime_fiber']):.8f}")
    print(
        "direct check  = "
        f"measR {fmt_fraction(metrics['direct_r'])}, "
        f"measQ {fmt_fraction(metrics['direct_q'])}, "
        f"both {fmt_fraction(metrics['direct_both'])}, "
        f"Rprime {fmt_fraction(metrics['direct_rprime'])}; "
        f"identity_ok={metrics['identity_ok']}"
    )
    return metrics


def proof_carrier_tournament() -> None:
    carriers = [
        (
            "fiber_pgf_conditional_mean",
            {
                "preserves_predicate",
                "exact_identity",
                "coefficient_data",
                "finite_checkable",
                "quotient_legal",
                "spec_transfer",
                "sidecar_memory",
            },
        ),
        (
            "signed_SPEC_resonance_lattice",
            {
                "preserves_predicate",
                "exact_identity",
                "finite_checkable",
                "spec_transfer",
                "tail_bound",
            },
        ),
        (
            "edge_child_global_quotient",
            {
                "preserves_predicate",
                "quotient_legal",
                "sidecar_memory",
                "finite_checkable",
                "edge_payload",
            },
        ),
        (
            "lee_yang_root_pgf",
            {
                "coefficient_data",
                "zero_free_signal",
                "sidecar_memory",
                "finite_checkable",
            },
        ),
        (
            "macwilliams_weight_enumerator",
            {
                "coefficient_data",
                "transform_dual",
                "finite_checkable",
                "sidecar_memory",
            },
        ),
        (
            "q_pochhammer_modular_tail",
            {
                "coefficient_data",
                "tail_bound",
                "transform_dual",
            },
        ),
        (
            "moser_fibbinary_partial_cube",
            {
                "sidecar_memory",
                "quotient_legal",
                "finite_checkable",
            },
        ),
        (
            "raw_scalar_Rprime",
            {
                "preserves_predicate",
            },
        ),
    ]
    tie_order = {name: i for i, (name, _) in enumerate(carriers)}
    edges: dict[str, set[str]] = {name: set() for name, _ in carriers}
    scores = Counter()
    flips_against_raw = 0

    for i, (a_name, a_axes) in enumerate(carriers):
        for b_name, b_axes in carriers[i + 1 :]:
            a_score = len(a_axes - b_axes)
            b_score = len(b_axes - a_axes)
            if a_score == b_score:
                winner, loser = (a_name, b_name) if tie_order[a_name] < tie_order[b_name] else (b_name, a_name)
            elif a_score > b_score:
                winner, loser = a_name, b_name
            else:
                winner, loser = b_name, a_name
            edges[winner].add(loser)
            scores[winner] += 1
            if loser == "raw_scalar_Rprime":
                flips_against_raw += 1

    hist = Counter(scores[name] for name, _ in carriers)

    def has_edge(a: str, b: str) -> bool:
        return b in edges[a]

    names = [name for name, _ in carriers]
    tri_cycles = 0
    for a, b, c in combinations(names, 3):
        if has_edge(a, b) and has_edge(b, c) and has_edge(c, a):
            tri_cycles += 1
        if has_edge(a, c) and has_edge(c, b) and has_edge(b, a):
            tri_cycles += 1

    selected_path = sorted(names, key=lambda name: (-scores[name], tie_order[name]))

    print("=" * 100)
    print("TOURNAMENT ANALYSIS OVER GENERATING-FUNCTION PROOF CARRIERS")
    print("vertices=proof carriers, not runners, arcs, or raw speed rows")
    print(
        "pairwise observable=axis majority over predicate retention, exact identity, "
        "coefficient data, finite checkability, quotient legality, SPEC transfer, "
        "tail bounds, zero-free signal, and sidecar memory"
    )
    print("switch=give A->B when A wins more axes; ties use the listed guardrail order")
    print(f"score_hist={dict(sorted(hist.items()))}")
    print(f"directed_3cycles={tri_cycles}")
    print(f"edge_flips_against_raw_scalar={flips_against_raw}")
    print(f"selected_hamiltonian_path={' -> '.join(selected_path)}")


def main() -> None:
    print("HYP-3140 / codex-2026-06-27-S273")
    print("14-sheet fiber PGF certificate for the LRC14 Rprime floor")
    print("Refines the HYP-3136 integrated multi-far factorization at its Rprime factor")
    print()
    print("ASSUMPTION CHALLENGE")
    print(
        "Vertices need not be runners.  This scout uses u-fibers and sheet-count "
        "coefficients as vertices/observables; runners are remembered only through "
        "the R-safe sheet packet.  The preserved predicate is the LRC14 covering "
        "floor Rprime>0.  The destroyed information is the identity of individual "
        "safe sheets unless the PGF coefficient or global-consistency sidecar keeps it."
    )
    print()
    print("IDENTITY")
    print("N_R(u)=number of R-safe lifts (u+a)/14, a=0..13")
    print("Rprime = E[N_R(u) | Q-lonely] / E[N_R(u)]")
    print()

    cases = [
        ("canonical r=1: drop 12, Q={1}", tuple([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13]), (1,)),
        ("canonical r=2: drop 13, Q={1,2} (HYP-3129 worst)", tuple(range(1, 13)), (1, 2)),
        ("canonical r=3: R={1..11}, Q={1,2,3}", tuple(range(1, 12)), (1, 2, 3)),
        ("canonical r=4: R={1..10}, Q={1,2,3,4}", tuple(range(1, 11)), (1, 2, 3, 4)),
        ("canonical r=5: R={1..9}, Q={1..5}", tuple(range(1, 10)), (1, 2, 3, 4, 5)),
        ("canonical r=6: R={1..8}, Q={1..6}", tuple(range(1, 9)), (1, 2, 3, 4, 5, 6)),
        ("spread far: R={1..9}, Q={1,3,5}", tuple(range(1, 10)), (1, 3, 5)),
        ("gcd-loaded R: even R, Q={1,3,5}", (2, 4, 6, 8, 10, 12), (1, 3, 5)),
        ("gcd-loaded R: triples, Q={1,2,5}", (3, 6, 9, 12), (1, 2, 5)),
        ("coprime-ish sparse R, Q={1,2,3}", (5, 9, 11, 13), (1, 2, 3)),
    ]

    rows = []
    for label, r_speeds, q_speeds in cases:
        metrics = print_case(label, r_speeds, q_speeds)
        rows.append((label, metrics))

    print("=" * 100)
    print("SUMMARY")
    print(f"{'case':<48}{'Rprime':>10}{'E[N]':>12}{'E[N|Q]':>12}{'Q zero':>12}{'Q range':>10}")
    min_row = None
    for label, metrics in rows:
        rprime = metrics["rprime_fiber"]
        if min_row is None or rprime < min_row[1]["rprime_fiber"]:
            min_row = (label, metrics)
        q_range = f"{metrics['min_q_count']}..{metrics['max_q_count']}"
        print(
            f"{label[:48]:<48}"
            f"{float(rprime):>10.5f}"
            f"{float(metrics['mean_sheets']):>12.5f}"
            f"{float(metrics['q_mean_sheets']):>12.5f}"
            f"{float(metrics['q_zero_mass']):>12.5f}"
            f"{q_range:>10}"
        )
    assert min_row is not None
    print()
    print(
        "Worst targeted row: "
        f"{min_row[0]} with Rprime={float(min_row[1]['rprime_fiber']):.8f}. "
        "It is visible as a low-degree fiber-PGF defect: Q-safe mass still "
        "contains y^0 and y^1 coefficients, so min-sheet positivity is too crude; "
        "the right finite theorem is a conditional first-moment bound."
    )

    print()
    print("PROOF TARGET EMITTED")
    print(
        "For every legal residual row, prove coefficientwise or moment-form "
        "inequality E[N_R | Q-lonely] >= c * E[N_R] for the 14-sheet PGF. "
        "HYP-3129's SPEC certificate is exactly the Fourier transform of this "
        "finite PGF statement; HYP-3134 supplies the legality condition for "
        "forgetting individual sheets after the global-consistency class is named."
    )
    proof_carrier_tournament()


if __name__ == "__main__":
    main()
