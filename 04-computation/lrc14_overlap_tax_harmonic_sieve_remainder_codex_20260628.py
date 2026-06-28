#!/usr/bin/env python3
r"""HYP-3434: Euler-Mascheroni harmonic-sieve remainder for LRC14.

HYP-3425/HYP-3426 reduce the two-branch relocation floor to the one-branch
interval-piercing target

    branch0 = E_safe \ B0_odd.

This scout asks what a harmonic or information-compression proof is allowed
to forget.  On each row it computes the exact first-sieve identity

    |branch0|
      = |E_safe| - sum_o |E_safe cap B0_o|
        + (sum_o |E_safe cap B0_o| - |E_safe cap union_o B0_o|).

The first two terms are the naive union-bound slack; the last term is the
overlap tax, the finite correction that survives after scalar compression.
Euler-Mascheroni enters only as a denominator-prefix calibration:

    H_N = log N + gamma + error_N,
    sum_{odd n <= N} 1/n
      = 1/2 (log N + gamma + log 2) + odd_error_N.

The proof target is not a numerical gamma estimate.  It is the exact statement
that endpoint-spine structure supplies enough overlap tax whenever naive
harmonic slack is negative.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from importlib import util
from pathlib import Path
import math
import sys


ROOT = Path(__file__).resolve().parents[1]
H3429_PATH = ROOT / "04-computation" / "lrc14_component_spine_certificate_codex_20260628.py"


def load_h3429():
    spec = util.spec_from_file_location("h3429_component_spine", H3429_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {H3429_PATH}")
    module = util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


h3429 = load_h3429()
h3425 = h3429.h3425

F = Fraction
ZERO = F(0)


def fmt(x: F | None) -> str:
    if x is None:
        return "n/a"
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def fmt_compact(x: F, digits: int = 9, max_chars: int = 42) -> str:
    exact = fmt(x)
    if len(exact) <= max_chars:
        return exact
    return fmt_float(float(x), digits)


def fmt_float(x: float, digits: int = 9) -> str:
    return f"{x:.{digits}f}"


def harmonic(n: int) -> F:
    return sum((F(1, k) for k in range(1, n + 1)), ZERO)


def odd_prefix_harmonic(n: int) -> F:
    """Return sum_{1 <= k <= n, k odd} 1/k exactly."""

    return harmonic(n) - F(1, 2) * harmonic(n // 2)


def gamma_constant() -> float:
    return getattr(math, "euler_gamma", 0.5772156649015328606)


def harmonic_gamma_residual(n: int) -> float:
    return float(harmonic(n)) - math.log(n) - gamma_constant()


def odd_prefix_gamma_residual(n: int) -> float:
    cap = float(odd_prefix_harmonic(n))
    return cap - 0.5 * math.log(n) - 0.5 * (gamma_constant() + math.log(2))


@dataclass(frozen=True)
class SieveRow:
    name: str
    speeds: tuple[int, ...]
    odd: tuple[int, ...]
    even_half: tuple[int, ...]
    even_measure: F
    restricted_single_sum: F
    restricted_union: F
    branch0_measure: F
    naive_slack: F
    overlap_tax: F
    odd_harmonic_pressure: F
    even_harmonic_pressure: F
    odd_prefix_cap: F
    odd_prefix_gap: F
    harmonic_residual: float
    odd_prefix_residual: float

    @property
    def naive_positive(self) -> bool:
        return self.naive_slack > ZERO

    @property
    def rescued_by_tax(self) -> bool:
        return self.naive_slack < ZERO and self.overlap_tax > -self.naive_slack

    @property
    def tax_to_deficit(self) -> float | None:
        if self.naive_slack >= ZERO:
            return None
        return float(self.overlap_tax / (-self.naive_slack))


def audit_row(name: str, speeds: tuple[int, ...]) -> SieveRow:
    speeds = tuple(sorted(set(speeds)))
    odd = tuple(v for v in speeds if v % 2 == 1)
    even_half = tuple(v // 2 for v in speeds if v % 2 == 0)
    even_safe = h3425.even_safe_intervals(even_half)
    even_measure = h3425.measure(even_safe)

    bad_parts = [h3425.branch0_bad_one(o) for o in odd]
    restricted_parts = [h3425.intersect_two(even_safe, part) for part in bad_parts]
    restricted_single_sum = sum((h3425.measure(part) for part in restricted_parts), ZERO)
    restricted_union = h3425.measure(h3425.intersect_two(even_safe, h3425.union_many(bad_parts)))
    branch0_measure = even_measure - restricted_union
    naive_slack = even_measure - restricted_single_sum
    overlap_tax = restricted_single_sum - restricted_union

    assert branch0_measure == naive_slack + overlap_tax
    assert branch0_measure == h3425.measure(
        h3425.intersect_two(even_safe, h3425.complement(h3425.union_many(bad_parts)))
    )

    max_odd = max(odd) if odd else 1
    odd_pressure = sum((F(1, o) for o in odd), ZERO)
    even_pressure = sum((F(1, e) for e in even_half if e), ZERO)
    odd_cap = odd_prefix_harmonic(max_odd)
    return SieveRow(
        name=name,
        speeds=speeds,
        odd=odd,
        even_half=even_half,
        even_measure=even_measure,
        restricted_single_sum=restricted_single_sum,
        restricted_union=restricted_union,
        branch0_measure=branch0_measure,
        naive_slack=naive_slack,
        overlap_tax=overlap_tax,
        odd_harmonic_pressure=odd_pressure,
        even_harmonic_pressure=even_pressure,
        odd_prefix_cap=odd_cap,
        odd_prefix_gap=odd_cap - odd_pressure,
        harmonic_residual=harmonic_gamma_residual(max_odd),
        odd_prefix_residual=odd_prefix_gamma_residual(max_odd),
    )


AXES = (
    "predicate_retention",
    "finite_exactness",
    "compression_loss_accounting",
    "endpoint_spine_interface",
    "two_adic_induction",
    "graph_cut_flow",
    "scalar_risk_control",
)


@dataclass(frozen=True)
class Carrier:
    name: str
    scores: tuple[int, ...]

    @property
    def total(self) -> int:
        return sum(self.scores)


CARRIERS = (
    Carrier("exact_overlap_tax_sieve", (10, 10, 10, 9, 8, 9, 10)),
    Carrier("endpoint_spine_tax_certificate", (10, 10, 9, 10, 8, 8, 10)),
    Carrier("conductance_graph_overlap_flow", (9, 9, 9, 8, 8, 10, 9)),
    Carrier("two_adic_even_child_loss_ledger", (9, 9, 8, 8, 10, 8, 9)),
    Carrier("odd_prefix_gamma_cap", (7, 8, 7, 6, 6, 5, 8)),
    Carrier("raw_harmonic_union_bound", (6, 6, 3, 4, 4, 3, 3)),
    Carrier("scalar_gamma_slogan", (3, 3, 1, 2, 2, 1, 1)),
)


def tournament() -> tuple[dict[int, int], list[str], int]:
    hist = dict(sorted(Counter(c.total for c in CARRIERS).items()))
    order = {carrier.name: i for i, carrier in enumerate(CARRIERS)}
    path = [
        carrier.name
        for carrier in sorted(CARRIERS, key=lambda c: (-c.total, order[c.name]))
    ]
    directed_3cycles = 0
    return hist, path, directed_3cycles


def selected_canonical_rows(rows: list[SieveRow]) -> list[tuple[int, SieveRow]]:
    by_name = {row.name: row for row in rows}
    selected: list[tuple[int, SieveRow]] = []
    for m in (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 18, 24, 30, 36, 42, 48):
        if m == 1 and "covering_AP_with_84" in by_name:
            selected.append((m, by_name["covering_AP_with_84"]))
            continue
        name = f"canonical_84m_{m:02d}" if m <= 12 else f"canonical_84m_ext_{m:02d}"
        if name in by_name:
            selected.append((m, by_name[name]))
    return selected


def print_row(label: str, row: SieveRow) -> None:
    ratio = row.tax_to_deficit
    print(f"  {label}: {row.name}")
    print(f"    speeds={row.speeds}")
    print(
        "    measures: "
        f"|E|={fmt(row.even_measure)}, "
        f"single_sum={fmt(row.restricted_single_sum)}, "
        f"union={fmt(row.restricted_union)}, "
        f"branch0={fmt(row.branch0_measure)}"
    )
    print(
        "    sieve: "
        f"naive_slack={fmt(row.naive_slack)}, "
        f"overlap_tax={fmt(row.overlap_tax)}, "
        f"tax/deficit={fmt_float(ratio, 6) if ratio is not None else 'n/a'}"
    )
    print(
        "    harmonic: "
        f"H_odd(row)={fmt(row.odd_harmonic_pressure)}, "
        f"H_evenhalf(row)={fmt(row.even_harmonic_pressure)}, "
        f"odd_prefix_gap={fmt_compact(row.odd_prefix_gap)}, "
        f"odd_gamma_residual={fmt_float(row.odd_prefix_residual, 9)}"
    )


def main() -> None:
    rows = [audit_row(name, speeds) for name, speeds in h3429.audited_rows()]
    identity_rows = [row for row in rows if row.branch0_measure == row.naive_slack + row.overlap_tax]
    branch_positive = [row for row in rows if row.branch0_measure > ZERO]
    naive_positive = [row for row in rows if row.naive_slack > ZERO]
    naive_negative = [row for row in rows if row.naive_slack < ZERO]
    overlap_positive = [row for row in rows if row.overlap_tax > ZERO]
    rescued = [row for row in rows if row.rescued_by_tax]
    tightest_branch = min(rows, key=lambda row: row.branch0_measure)
    worst_naive = min(rows, key=lambda row: row.naive_slack)
    max_tax = max(rows, key=lambda row: row.overlap_tax)
    min_tax = min(rows, key=lambda row: row.overlap_tax)
    max_odd_pressure = max(rows, key=lambda row: row.odd_harmonic_pressure)
    max_prefix_gap = max(rows, key=lambda row: row.odd_prefix_gap)
    tightest_rescue = min(naive_negative, key=lambda row: row.tax_to_deficit or float("inf"))
    hist, path, directed_3cycles = tournament()

    print("HYP-3434 EULER-MASCHERONI HARMONIC-SIEVE REMAINDER")
    print("=" * 78)
    print("Identity:")
    print("  branch0 = E_safe minus union_o B0_o")
    print("  branch0_mass = naive_slack + overlap_tax")
    print("  naive_slack = |E_safe| - sum_o |E_safe cap B0_o|")
    print("  overlap_tax = sum_o |E_safe cap B0_o| - |E_safe cap union_o B0_o|")
    print("  gamma calibrates harmonic denominator compression, not the final proof.")
    print()

    print("A. Aggregate exact audit")
    print(f"  rows audited:                         {len(rows)}")
    print(f"  exact sieve identity rows:            {len(identity_rows)}/{len(rows)}")
    print(f"  positive one-branch survivor rows:    {len(branch_positive)}/{len(rows)}")
    print(f"  naive slack positive rows:            {len(naive_positive)}/{len(rows)}")
    print(f"  naive slack negative rows:            {len(naive_negative)}/{len(rows)}")
    print(f"  positive overlap-tax rows:            {len(overlap_positive)}/{len(rows)}")
    print(f"  negative-slack rows rescued by tax:   {len(rescued)}/{len(naive_negative)}")
    print(f"  smallest branch0 mass:                {fmt(tightest_branch.branch0_measure)} ({tightest_branch.name})")
    print(f"  worst naive slack:                    {fmt(worst_naive.naive_slack)} ({worst_naive.name})")
    print(f"  max overlap tax:                      {fmt(max_tax.overlap_tax)} ({max_tax.name})")
    print(f"  min overlap tax:                      {fmt(min_tax.overlap_tax)} ({min_tax.name})")
    print(
        "  tightest tax/deficit rescue:         "
        f"{fmt_float(tightest_rescue.tax_to_deficit or 0.0, 6)} ({tightest_rescue.name})"
    )
    print(
        "  H_N-logN-gamma residual range:       "
        f"[{fmt_float(min(row.harmonic_residual for row in rows), 9)}, "
        f"{fmt_float(max(row.harmonic_residual for row in rows), 9)}]"
    )
    print(
        "  odd-prefix gamma residual range:     "
        f"[{fmt_float(min(row.odd_prefix_residual for row in rows), 9)}, "
        f"{fmt_float(max(row.odd_prefix_residual for row in rows), 9)}]"
    )
    print(f"  max row odd harmonic pressure:        {fmt(max_odd_pressure.odd_harmonic_pressure)} ({max_odd_pressure.name})")
    print(
        "  max omitted odd-prefix mass:          "
        f"{fmt_compact(max_prefix_gap.odd_prefix_gap)} ({max_prefix_gap.name})"
    )
    print()

    print("B. Critical rows")
    print_row("tightest_branch_and_tightest_rescue", tightest_branch)
    print_row("worst_naive_slack", worst_naive)
    print_row("smallest_overlap_tax", min_tax)
    print()

    print("C. Canonical 84m tower sieve profile")
    print("  m | branch0 | naive_slack | overlap_tax | tax/deficit | odd_residual")
    for m, row in selected_canonical_rows(rows):
        ratio = row.tax_to_deficit
        ratio_text = fmt_float(ratio, 4) if ratio is not None else "n/a"
        print(
            f"  {m:2d} | {fmt(row.branch0_measure):>10} | {fmt(row.naive_slack):>12} | "
            f"{fmt(row.overlap_tax):>11} | {ratio_text:>11} | "
            f"{fmt_float(row.odd_prefix_residual, 7)}"
        )
    print()

    print("D. Finite lemma target")
    print("  Strengthen the one-branch piercing lemma to this exact dichotomy:")
    print("    either naive_slack >= 0,")
    print("    or an endpoint-spine certified overlap_tax exceeds -naive_slack.")
    print("  The observed hard core is not a resonant survivor.  It is a compression")
    print("  failure: scalar union pressure forgets intersections inside E_safe.")
    print("  The missing associator is the overlap tax.  Endpoint spines from HYP-3429")
    print("  should localize that tax to rank <= 2 wall certificates.")
    print()

    print("E. Tournament Analysis")
    print("  vertices are proof carriers/compression residues, not runners.")
    print("  pairwise_observable=retained floor predicate plus exact correction debt.")
    print("  switch_gauge=higher weighted retained payload; ties by carrier order.")
    print(f"  axes={','.join(AXES)}")
    print(f"  score_hist={hist}")
    print(f"  directed_3cycles={directed_3cycles}")
    print("  hamiltonian_path=" + " -> ".join(path))
    print()

    print("F. Assumption challenge")
    print("  Considered vertices: runners, odd speeds, even-half speeds, E_safe")
    print("  components, survivor windows, endpoint labels, wall-crossing events,")
    print("  harmonic denominators, graph cuts, and proof obligations.")
    print("  Chosen vertices: compression residues and proof carriers.  This preserves")
    print("  the LRC predicate because each residue is still computed from exact")
    print("  interval algebra, but it destroys raw row order and most runner labels.")
    print("  Challenged assumption: a commutative scalar sum of bad masses can replace")
    print("  the associative set union.  The exact data says the scalar compression")
    print("  fails precisely by the overlap-tax term.")


if __name__ == "__main__":
    main()
