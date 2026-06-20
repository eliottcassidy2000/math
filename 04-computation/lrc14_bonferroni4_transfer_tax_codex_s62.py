#!/usr/bin/env python3
"""HYP-2696 / T932: Bonferroni4 transfer tax for the LRC14 sector DP.

THM-556 says the six-sector level-4 Bonferroni upper bound is

    U4 = p0 + p5 + 5*p6.

THM-555/HYP-2691 says a one-speed insertion can only keep a missed-sector
state fixed or delete one missed sector.  Combining the two gives an exact
local recurrence:

    Delta U4(P,e) = mass(1 -> 0) - mass(5 -> 4) - 4*mass(6 -> 5).

Thus every positive movement of the final Bonferroni4 certificate is a genuine
one-missed closure, while high-missed states pay a signed transfer tax.  This
script verifies the identity on the shared HYP-2691 row bank and records which
insertion schedules keep the proof burden smallest.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, permutations
import sys

from lrc14_state_transfer_dp_codex_s58 import (
    OrderRun,
    TransferStep,
    orders_for,
    row_bank,
    run_order,
)
from lrc14_wide_branch_ridge_codex_s47 import CAP, Row, fmt, missed_distribution

sys.stdout.reconfigure(line_buffering=True)


def u4_from_dist(dist: tuple[Fraction, ...]) -> Fraction:
    return dist[0] + dist[5] + 5 * dist[6]


def level4_weight(missed_count: int) -> int:
    """Coefficient of a state with this missed count in U4."""

    if missed_count == 0:
        return 1
    if missed_count == 5:
        return 1
    if missed_count == 6:
        return 5
    return 0


@dataclass(frozen=True)
class StepTax:
    step: TransferStep
    delta_u4: Fraction
    close_one: Fraction
    hit_five: Fraction
    hit_six: Fraction
    transfer_tax: Fraction
    hit_by_count: tuple[Fraction, ...]
    fixed_by_count: tuple[Fraction, ...]

    @property
    def positive_delta(self) -> Fraction:
        return max(self.delta_u4, Fraction(0))


@dataclass(frozen=True)
class RunTaxAudit:
    run: OrderRun
    steps: tuple[StepTax, ...]

    @property
    def final_dist(self) -> tuple[Fraction, ...]:
        return self.run.final_profile.p_by_count

    @property
    def final_u4(self) -> Fraction:
        return u4_from_dist(self.final_dist)

    @property
    def final_p0(self) -> Fraction:
        return self.final_dist[0]

    @property
    def final_tail45(self) -> Fraction:
        return self.final_dist[5] + 5 * self.final_dist[6]

    @property
    def total_close_one(self) -> Fraction:
        return sum((step.close_one for step in self.steps), Fraction(0))

    @property
    def total_hit_five(self) -> Fraction:
        return sum((step.hit_five for step in self.steps), Fraction(0))

    @property
    def total_hit_six(self) -> Fraction:
        return sum((step.hit_six for step in self.steps), Fraction(0))

    @property
    def total_tax(self) -> Fraction:
        return sum((step.transfer_tax for step in self.steps), Fraction(0))

    @property
    def positive_u4_area(self) -> Fraction:
        return sum((step.positive_delta for step in self.steps), Fraction(0))

    @property
    def max_positive_delta(self) -> Fraction:
        return max((step.positive_delta for step in self.steps), default=Fraction(0))

    @property
    def first_closure_step(self) -> int:
        for idx, step in enumerate(self.steps, 1):
            if step.close_one:
                return idx
        return 0

    @property
    def last_positive_step(self) -> int:
        out = 0
        for idx, step in enumerate(self.steps, 1):
            if step.delta_u4 > 0:
                out = idx
        return out

    @property
    def tax_before_first_closure(self) -> Fraction:
        first = self.first_closure_step
        if not first:
            return Fraction(0)
        return sum((step.transfer_tax for step in self.steps[: first - 1]), Fraction(0))

    @property
    def max_prefix_u4_after_first(self) -> Fraction:
        prefix_u4 = Fraction(5)
        best = prefix_u4
        for step in self.steps:
            prefix_u4 += step.delta_u4
            best = max(best, prefix_u4)
        return best


def tax_step(step: TransferStep) -> StepTax:
    hit_by_count = [Fraction(0) for _ in range(7)]
    fixed_by_count = [Fraction(0) for _ in range(7)]
    direct_delta = Fraction(0)
    close_one = Fraction(0)
    hit_five = Fraction(0)
    hit_six = Fraction(0)

    for before, after, mass in step.transitions:
        b = len(before)
        a = len(after)
        direct_delta += Fraction(level4_weight(a) - level4_weight(b)) * mass
        if b == a:
            fixed_by_count[b] += mass
            continue
        assert b - a == 1
        hit_by_count[b] += mass
        if b == 1 and a == 0:
            close_one += mass
        elif b == 5 and a == 4:
            hit_five += mass
        elif b == 6 and a == 5:
            hit_six += mass

    formula_delta = close_one - hit_five - 4 * hit_six
    before_u4 = u4_from_dist(step.before_profile.p_by_count)
    after_u4 = u4_from_dist(step.after_profile.p_by_count)
    exact_delta = after_u4 - before_u4

    assert close_one == step.delta_p0
    assert direct_delta == exact_delta
    assert formula_delta == exact_delta

    return StepTax(
        step=step,
        delta_u4=exact_delta,
        close_one=close_one,
        hit_five=hit_five,
        hit_six=hit_six,
        transfer_tax=hit_five + 4 * hit_six,
        hit_by_count=tuple(hit_by_count),
        fixed_by_count=tuple(fixed_by_count),
    )


def audit_run(run: OrderRun) -> RunTaxAudit:
    audit = RunTaxAudit(run=run, steps=tuple(tax_step(step) for step in run.steps))
    assert audit.total_close_one == audit.final_p0
    assert audit.final_u4 - 5 == audit.total_close_one - audit.total_tax
    assert audit.total_tax == 5 + audit.final_p0 - audit.final_u4
    assert audit.final_tail45 == 5 - audit.total_tax
    return audit


def unique_order_runs(row: Row) -> list[OrderRun]:
    dedup: dict[tuple[int, ...], OrderRun] = {}
    for name, order in orders_for(row).items():
        dedup.setdefault(order, run_order(row, name, order))
    return list(dedup.values())


def audit_key(audit: RunTaxAudit) -> tuple[Fraction, int, Fraction, Fraction, int, str]:
    """Order-dependent proof-burden key; final U4 itself is order invariant."""

    return (
        audit.positive_u4_area,
        audit.run.max_support,
        audit.max_positive_delta,
        -audit.tax_before_first_closure,
        audit.last_positive_step,
        audit.run.name,
    )


def tournament_fingerprint(audits: list[RunTaxAudit]) -> str:
    n = len(audits)
    wins = [0 for _ in range(n)]
    edges: set[tuple[int, int]] = set()

    for i, j in combinations(range(n), 2):
        if audit_key(audits[i]) <= audit_key(audits[j]):
            winner, loser = i, j
        else:
            winner, loser = j, i
        wins[winner] += 1
        edges.add((winner, loser))

    cycles = 0
    for i, j, k in combinations(range(n), 3):
        if (
            (i, j) in edges and (j, k) in edges and (k, i) in edges
        ) or (
            (i, k) in edges and (k, j) in edges and (j, i) in edges
        ):
            cycles += 1

    hp_count = 0
    for perm in permutations(range(n)):
        if all((perm[i], perm[i + 1]) in edges for i in range(n - 1)):
            hp_count += 1

    path = sorted(range(n), key=lambda idx: (-wins[idx], audit_key(audits[idx])))
    score_hist = dict(sorted(Counter(wins).items()))
    return (
        f"score_hist={score_hist} directed_3cycles={cycles} HP={hp_count} "
        + "path="
        + " > ".join(audits[idx].run.name for idx in path)
    )


def print_step_tax_ledger(audit: RunTaxAudit, keep: int = 12) -> None:
    print(f"    transfer-tax ledger for {audit.run.name}:")
    prefix_u4 = Fraction(5)
    for idx, step_tax in enumerate(audit.steps[:keep], 1):
        prefix_u4 += step_tax.delta_u4
        print(
            f"      {idx:2d}. add={step_tax.step.added:3d} "
            f"dU4={step_tax.delta_u4} U4={prefix_u4} "
            f"close1={step_tax.close_one} hit5={step_tax.hit_five} "
            f"hit6={step_tax.hit_six} tax={step_tax.transfer_tax} "
            f"hits_by_count={tuple(step_tax.hit_by_count)}"
        )
    if len(audit.steps) > keep:
        print(f"      ... {len(audit.steps) - keep} more steps omitted")


def print_run_summary(audit: RunTaxAudit) -> None:
    run = audit.run
    print(
        f"  {run.name:<21} order={run.order} "
        f"pos_area={audit.positive_u4_area} max_pos={audit.max_positive_delta} "
        f"first_close={audit.first_closure_step} last_pos={audit.last_positive_step} "
        f"preclose_tax={audit.tax_before_first_closure} "
        f"max_support={run.max_support:2d} area_support={run.area_support:3d}"
    )


def print_row_audit(label: str, row: Row) -> RunTaxAudit:
    dist = missed_distribution(row)
    u4 = u4_from_dist(dist)
    cap = CAP.get(len(row))
    print("=" * 90)
    print(f"ROW {label}: E={row}")
    print(
        f"  k={len(row)} span={row[-1]} p0={fmt(dist[0])} "
        f"tail45={fmt(dist[5] + 5 * dist[6])} U4={fmt(u4)} "
        + (f"cap-U4={fmt(cap - u4)}" if cap is not None else "cap-U4=n/a")
    )
    print(f"  missed_dist={[str(x) for x in dist]}")

    audits = [audit_run(run) for run in unique_order_runs(row)]
    audits.sort(key=audit_key)
    for audit in audits:
        print_run_summary(audit)

    print("  Tournament Analysis on insertion schedules")
    print("    vertices: schedules / transfer-tax proof states, not runners")
    print("    observable: smaller positive U4 area, support, max positive jump, then earlier closure control")
    print(f"    {tournament_fingerprint(audits)}")

    best = audits[0]
    print(
        "  best-order identity: "
        f"final_U4-5={best.final_u4 - 5} = close1-tax "
        f"= {best.total_close_one}-{best.total_tax}; "
        f"final_tail45=5-tax={best.final_tail45}"
    )
    print(
        f"  total close1={best.total_close_one} hit5={best.total_hit_five} "
        f"hit6={best.total_hit_six} tax={best.total_tax}"
    )
    print_step_tax_ledger(best)
    return best


def main() -> None:
    print("HYP-2696 / T932 -- LRC14 Bonferroni4 transfer tax")
    print("Arithmetic: exact Fractions over HYP-2691 common sector-wall refinements.")
    print("Preserved predicate: p0(E) <= U4(E), with U4(E)=p0+p5+5p6.")
    print("Assumption challenged: proof vertices are missed-state transitions and schedules, not runners.\n")

    print("EXACT LOCAL LAW")
    print("  For every insertion P -> P union {e}:")
    print("    Delta U4 = mass(1->0) - mass(5->4) - 4*mass(6->5).")
    print("  Positive Bonferroni4 movement is exactly one-missed closure;")
    print("  five/six-missed transitions are a signed tax that lowers U4.\n")

    best_audits: list[tuple[str, RunTaxAudit]] = []
    for label, row in row_bank():
        best_audits.append((label, print_row_audit(label, row)))

    print("=" * 90)
    print("GLOBAL BEST-SCHEDULE SUMMARY")
    for label, audit in best_audits:
        cap = CAP.get(len(audit.run.row))
        slack = cap - audit.final_u4 if cap is not None else None
        print(
            f"  {label:<20} schedule={audit.run.name:<20} "
            f"U4={fmt(audit.final_u4)} "
            + (f"cap-U4={fmt(slack)} " if slack is not None else "")
            + f"pos_area={audit.positive_u4_area} tax={audit.total_tax} "
            f"tail45={audit.final_tail45}"
        )

    print("\nSYNTHESIS")
    print("  The transfer law is exact on every audited row and schedule.")
    print("  Boundary/AP rows can fail the final U4 gate because their high tail")
    print("  survives as finite template mass.  True-wide rows in the bank have")
    print("  enough transfer tax to keep U4 below cap after the same p0 closures.")
    print("  New proof target: classify or bound the possible unpaid 1->0 closures")
    print("  after the 5->4 and 6->5 high-tail tax has been exhausted.")


if __name__ == "__main__":
    main()
