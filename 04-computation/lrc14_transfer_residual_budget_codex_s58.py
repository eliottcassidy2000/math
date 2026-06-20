#!/usr/bin/env python3
"""Residual-budget audit for the HYP-2691 sector-state transfer DP.

THM-554 turns one insertion into

    delta = p0(P union {e}) - p0(P).

The decorrelated value of the same transition is p1(P)/7, because only atoms
where P misses exactly one sector can close, and a decorrelated new speed hits
that needed sector with probability 1/7.  Thus

    residual = delta - p1(P)/7

is the DP-local form of the one-far discrepancy from THM-546.  This script
keeps everything exact and compares each prefix residual with the elementary
component budget (6/49)*V(P)/e, where V(P) is the number of one-missed-sector
runs over the prefix wall partition.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from lrc14_state_transfer_dp_codex_s58 import (
    INNER,
    OrderRun,
    mask_at,
    missed_tuple,
    orders_for,
    row_bank,
    run_order,
    state_profile,
    transfer_step,
)
from lrc14_wide_branch_ridge_codex_s47 import CAP, Row, fmt, missed_distribution, wall_breakpoints


@dataclass(frozen=True)
class ResidualStep:
    added: int
    p1_before: Fraction
    delta: Fraction
    decorrelated: Fraction
    residual: Fraction
    one_missed_runs: int
    bound: Fraction

    @property
    def abs_residual(self) -> Fraction:
        return abs(self.residual)

    @property
    def pressure(self) -> Fraction:
        if self.bound == 0:
            return Fraction(0) if self.residual == 0 else Fraction(10**9)
        return self.abs_residual / self.bound

    @property
    def scaled_abs(self) -> Fraction:
        return self.abs_residual * self.added


@dataclass(frozen=True)
class ResidualRun:
    label: str
    order_name: str
    run: OrderRun
    steps: tuple[ResidualStep, ...]

    @property
    def total_abs(self) -> Fraction:
        return sum((s.abs_residual for s in self.steps), Fraction(0))

    @property
    def total_bound(self) -> Fraction:
        return sum((s.bound for s in self.steps), Fraction(0))

    @property
    def max_pressure(self) -> Fraction:
        return max((s.pressure for s in self.steps), default=Fraction(0))

    @property
    def max_scaled_abs(self) -> Fraction:
        return max((s.scaled_abs for s in self.steps), default=Fraction(0))


def one_missed_run_count(row: Row) -> int:
    row = tuple(sorted(set(row)))
    if not any(row):
        return 0
    d, bps = wall_breakpoints(row)
    l = d // 7
    den2 = 2 * l
    labels: list[int | None] = []
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo:
            continue
        state = missed_tuple(mask_at(row, lo + hi, den2))
        labels.append(state[0] if len(state) == 1 else None)
    if not labels:
        return 0

    total = 0
    for sec in INNER:
        if not any(label == sec for label in labels):
            continue
        prev = labels[-1]
        for label in labels:
            if label == sec and prev != sec:
                total += 1
            prev = label
    return total


def residual_steps(run: OrderRun) -> tuple[ResidualStep, ...]:
    out: list[ResidualStep] = []
    prefix: Row = (0,)
    for added in run.order:
        step = transfer_step(prefix, added)
        prof = state_profile(prefix)
        p1 = prof.p_by_count[1]
        decor = p1 / 7
        residual = step.delta_p0 - decor
        V = one_missed_run_count(prefix)
        bound = Fraction(6 * V, 49 * added)
        # The bound is intentionally the coarse THM-546 component budget.
        assert abs(residual) <= bound
        out.append(
            ResidualStep(
                added=added,
                p1_before=p1,
                delta=step.delta_p0,
                decorrelated=decor,
                residual=residual,
                one_missed_runs=V,
                bound=bound,
            )
        )
        prefix = step.after
    return tuple(out)


def best_runs_for(row: Row) -> list[OrderRun]:
    runs = [run_order(row, name, order) for name, order in orders_for(row).items()]
    dedup: dict[tuple[int, ...], OrderRun] = {}
    for run in runs:
        dedup.setdefault(run.order, run)
    runs = list(dedup.values())
    runs.sort(key=lambda r: (r.max_support, r.area_support, r.max_transitions, r.final_entropy, r.name))
    chosen = [runs[0]]
    inc = next((r for r in runs if r.name == "increasing"), None)
    if inc is not None and inc.order != chosen[0].order:
        chosen.append(inc)
    return chosen


def print_residual_run(rr: ResidualRun) -> None:
    print(
        f"  order={rr.order_name:<21} total_abs={rr.total_abs} ({float(rr.total_abs):.9f}) "
        f"total_bound={rr.total_bound} ({float(rr.total_bound):.9f}) "
        f"max_pressure={rr.max_pressure} ({float(rr.max_pressure):.6f}) "
        f"max_e_abs_residual={rr.max_scaled_abs} ({float(rr.max_scaled_abs):.9f})"
    )
    for idx, step in enumerate(rr.steps, 1):
        if step.p1_before == 0 and step.delta == 0:
            continue
        print(
            f"    {idx:2d}. add={step.added:3d} p1={step.p1_before} "
            f"decor={step.decorrelated} delta={step.delta} residual={step.residual} "
            f"V={step.one_missed_runs:2d} bound={step.bound} pressure={step.pressure}"
        )


def main() -> None:
    print("HYP-2691 residual-budget audit")
    print("residual = exact insertion delta - p1(prefix)/7")
    print("budget = (6/49)*V(prefix)/added_speed, exact Fractions")
    print("V(prefix) counts circular one-missed-sector runs by sector.\n")

    all_residuals: list[ResidualRun] = []
    for label, row in row_bank():
        dist = missed_distribution(row)
        cap = CAP.get(len(row))
        print("=" * 88)
        print(f"ROW {label}: E={row}")
        if cap is not None:
            print(f"  p0={fmt(dist[0])} cap={fmt(cap)} margin={fmt(cap - dist[0])}")
        else:
            print(f"  p0={fmt(dist[0])}")

        for run in best_runs_for(row):
            rr = ResidualRun(
                label=label,
                order_name=run.name,
                run=run,
                steps=residual_steps(run),
            )
            all_residuals.append(rr)
            print_residual_run(rr)
        print()

    print("=" * 88)
    print("GLOBAL RESIDUAL LEADERS")
    for rr in sorted(all_residuals, key=lambda r: (r.max_pressure, r.max_scaled_abs), reverse=True)[:12]:
        print(
            f"  {rr.label:<20} {rr.order_name:<21} "
            f"max_pressure={rr.max_pressure} ({float(rr.max_pressure):.6f}) "
            f"max_e_abs={rr.max_scaled_abs} ({float(rr.max_scaled_abs):.9f}) "
            f"total_abs={rr.total_abs} ({float(rr.total_abs):.9f})"
        )

    print("\nSYNTHESIS")
    print("  1. The DP-local residual is exactly the old one-far discrepancy, but")
    print("     now it appears at every prefix rather than only in a final far step.")
    print("  2. The coarse component budget holds on all audited prefixes.  Its slack")
    print("     is large in many early/high-support steps, suggesting room for a")
    print("     lossy cap-level constant.")
    print("  3. The high-pressure steps are finite structured prefixes.  Those are the")
    print("     candidates to classify by AP/dyadic/cube-root/Ruzsa address before")
    print("     sending the rest to Weyl/BV decorrelation.")


if __name__ == "__main__":
    main()
