#!/usr/bin/env python3
"""HYP-2691 / T927: exact sector-state transfer DP for LRC(14).

The dynamic-programming lens is to keep the missed-sector address until the
last possible moment.  For a prefix P and a new speed e, one wall atom has a
missed-sector state A subset {1,...,6}; after adding e, the new state is either
A or A with exactly one sector removed.  Thus the scalar increment

    p0(P union {e}) - p0(P)

is not a black box: it is the exact mass where P misses exactly one sector and
the new speed lands in that sector.  This script computes those transition
packets over the exact common wall refinement, then asks which runner order
keeps the DP state small.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import log2

from lrc14_wide_branch_ridge_codex_s47 import (
    CAP,
    Row,
    additive_energy,
    fmt,
    missed_distribution,
    normalized,
    squarefree_profile,
    sumset_excess,
    wall_breakpoints,
)


INNER = tuple(range(1, 7))
ALL_INNER = 0b1111110


def missed_tuple(mask: int) -> tuple[int, ...]:
    return tuple(j for j in INNER if not (mask & (1 << j)))


def mask_at(row: Row, midnum: int, den2: int) -> int:
    mask = 0
    for e in row:
        if e:
            mask |= 1 << ((e * midnum // den2) % 7)
    return mask


def entropy(masses: list[Fraction]) -> float:
    total = sum(masses, Fraction(0))
    if not total:
        return 0.0
    out = 0.0
    for mass in masses:
        if mass:
            p = float(mass / total)
            out -= p * log2(p)
    return out


@dataclass(frozen=True)
class StateProfile:
    row: Row
    state_mass: tuple[tuple[tuple[int, ...], Fraction], ...]
    p_by_count: tuple[Fraction, ...]
    support: int
    entropy: float

    @property
    def p0(self) -> Fraction:
        return self.p_by_count[0]


@dataclass(frozen=True)
class TransferStep:
    before: Row
    added: int
    after: Row
    before_profile: StateProfile
    after_profile: StateProfile
    transitions: tuple[tuple[tuple[int, ...], tuple[int, ...], Fraction], ...]
    sector_gain: tuple[Fraction, ...]
    closure_by_missing_count: tuple[Fraction, ...]

    @property
    def delta_p0(self) -> Fraction:
        return self.after_profile.p0 - self.before_profile.p0

    @property
    def transition_count(self) -> int:
        return len(self.transitions)

    @property
    def total_gain(self) -> Fraction:
        return sum(self.sector_gain, Fraction(0))


@dataclass(frozen=True)
class OrderRun:
    name: str
    row: Row
    order: tuple[int, ...]
    steps: tuple[TransferStep, ...]
    final_profile: StateProfile

    @property
    def max_support(self) -> int:
        return max(step.after_profile.support for step in self.steps)

    @property
    def area_support(self) -> int:
        return sum(step.after_profile.support for step in self.steps)

    @property
    def max_transitions(self) -> int:
        return max(step.transition_count for step in self.steps)

    @property
    def final_entropy(self) -> float:
        return self.final_profile.entropy

    @property
    def final_p0(self) -> Fraction:
        return self.final_profile.p0

    @property
    def first_positive_step(self) -> int:
        for idx, step in enumerate(self.steps, 1):
            if step.after_profile.p0:
                return idx
        return 0


def state_profile(row: Row) -> StateProfile:
    row = normalized(row)
    nonzero = [e for e in row if e]
    if not nonzero:
        state = tuple(INNER)
        return StateProfile(
            row=row,
            state_mass=((state, Fraction(1)),),
            p_by_count=(Fraction(0),) * 6 + (Fraction(1),),
            support=1,
            entropy=0.0,
        )

    d, bps = wall_breakpoints(row)
    l = d // 7
    den2 = 2 * l
    states: Counter[tuple[int, ...]] = Counter()
    counts = [Fraction(0) for _ in range(7)]
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo:
            continue
        mass = Fraction(hi - lo, d)
        mask = mask_at(row, lo + hi, den2)
        miss = missed_tuple(mask)
        states[miss] += mass
        counts[len(miss)] += mass

    items = tuple(sorted(states.items(), key=lambda kv: (len(kv[0]), kv[0])))
    return StateProfile(
        row=row,
        state_mass=items,
        p_by_count=tuple(counts),
        support=len(items),
        entropy=entropy([mass for _state, mass in items]),
    )


def transfer_step(before: Row, added: int) -> TransferStep:
    before = normalized(before)
    after = normalized(before + (added,))
    d, bps = wall_breakpoints(after)
    l = d // 7
    den2 = 2 * l
    trans: Counter[tuple[tuple[int, ...], tuple[int, ...]]] = Counter()
    sector_gain = [Fraction(0) for _ in INNER]
    closure = [Fraction(0) for _ in range(7)]

    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo:
            continue
        mass = Fraction(hi - lo, d)
        midnum = lo + hi
        before_mask = mask_at(before, midnum, den2)
        after_mask = before_mask | (1 << ((added * midnum // den2) % 7))
        before_miss = missed_tuple(before_mask)
        after_miss = missed_tuple(after_mask)
        assert set(after_miss) <= set(before_miss)
        assert len(before_miss) - len(after_miss) in (0, 1)
        trans[(before_miss, after_miss)] += mass

        sec = (added * midnum // den2) % 7
        if sec in INNER and not (before_mask & (1 << sec)):
            sector_gain[sec - 1] += mass
        if before_miss and not after_miss:
            closure[len(before_miss)] += mass

    transitions = tuple(
        (a, b, mass)
        for (a, b), mass in sorted(trans.items(), key=lambda kv: (len(kv[0][0]), kv[0][0], kv[0][1]))
    )
    before_profile = state_profile(before)
    after_profile = state_profile(after)
    assert sum(closure, Fraction(0)) == after_profile.p0 - before_profile.p0
    return TransferStep(
        before=before,
        added=added,
        after=after,
        before_profile=before_profile,
        after_profile=after_profile,
        transitions=transitions,
        sector_gain=tuple(sector_gain),
        closure_by_missing_count=tuple(closure),
    )


def two_adic_order(n: int) -> int:
    out = 0
    while n and n % 2 == 0:
        out += 1
        n //= 2
    return out


def greedy_order(row: Row, mode: str) -> tuple[int, ...]:
    prefix: Row = (0,)
    remaining = set(e for e in row if e)
    order: list[int] = []
    while remaining:
        scored = []
        for e in remaining:
            step = transfer_step(prefix, e)
            if mode == "min-support":
                score = (
                    step.after_profile.support,
                    step.after_profile.entropy,
                    -float(step.after_profile.p0),
                    e,
                )
                scored.append((score, e))
            elif mode == "max-gain":
                score = (
                    -step.total_gain,
                    -step.delta_p0,
                    step.after_profile.support,
                    e,
                )
                scored.append((score, e))
            elif mode == "min-transition":
                score = (
                    step.transition_count,
                    step.after_profile.support,
                    step.after_profile.entropy,
                    e,
                )
                scored.append((score, e))
            else:
                raise ValueError(mode)
        scored.sort()
        chosen = scored[0][1]
        order.append(chosen)
        remaining.remove(chosen)
        prefix = normalized(prefix + (chosen,))
    return tuple(order)


def orders_for(row: Row) -> dict[str, tuple[int, ...]]:
    nonzero = tuple(e for e in row if e)
    return {
        "increasing": tuple(sorted(nonzero)),
        "decreasing": tuple(sorted(nonzero, reverse=True)),
        "residue-layer": tuple(sorted(nonzero, key=lambda x: (x % 7, x))),
        "dyadic-tower": tuple(sorted(nonzero, key=lambda x: (-two_adic_order(x), x))),
        "greedy-min-support": greedy_order(row, "min-support"),
        "greedy-max-gain": greedy_order(row, "max-gain"),
        "greedy-min-transition": greedy_order(row, "min-transition"),
    }


def run_order(row: Row, name: str, order: tuple[int, ...]) -> OrderRun:
    prefix: Row = (0,)
    steps: list[TransferStep] = []
    for e in order:
        step = transfer_step(prefix, e)
        steps.append(step)
        prefix = step.after
    return OrderRun(name=name, row=row, order=order, steps=tuple(steps), final_profile=state_profile(prefix))


def top_transitions(step: TransferStep, keep: int = 3) -> list[tuple[tuple[int, ...], tuple[int, ...], Fraction]]:
    return sorted(step.transitions, key=lambda x: (x[2], -len(x[0]), x[0]), reverse=True)[:keep]


def print_order_summary(run: OrderRun) -> None:
    print(
        f"  {run.name:<21} order={run.order} max_support={run.max_support:2d} "
        f"area_support={run.area_support:3d} final_support={run.final_profile.support:2d} "
        f"H={run.final_entropy:.4f} max_trans={run.max_transitions:2d} "
        f"first_p0_step={run.first_positive_step}"
    )


def print_step_ledger(run: OrderRun) -> None:
    print(f"    step ledger for {run.name}:")
    for idx, step in enumerate(run.steps, 1):
        closure = step.closure_by_missing_count
        gain = ",".join(str(x) for x in step.sector_gain)
        print(
            f"      {idx:2d}. add={step.added:3d} support={step.after_profile.support:2d} "
            f"H={step.after_profile.entropy:.4f} p0={step.after_profile.p0} "
            f"delta={step.delta_p0} gain=({gain}) transitions={step.transition_count} "
            f"closure_by_miss={tuple(closure)}"
        )
        tops = "; ".join(f"{a}->{b}:{mass}" for a, b, mass in top_transitions(step))
        print(f"          top transitions: {tops}")


def tournament_order_analysis(runs: list[OrderRun]) -> None:
    """Tournament on order-strategy proof states for one final row."""

    wins = Counter()
    edges: dict[tuple[int, int], int] = {}
    for i, a in enumerate(runs):
        for j, b in enumerate(runs):
            if i >= j:
                continue
            score_a = (
                -a.max_support,
                -a.area_support,
                -a.max_transitions,
                -a.final_entropy,
                a.name,
            )
            score_b = (
                -b.max_support,
                -b.area_support,
                -b.max_transitions,
                -b.final_entropy,
                b.name,
            )
            winner = i if score_a >= score_b else j
            loser = j if winner == i else i
            wins[winner] += 1
            wins.setdefault(loser, wins[loser])
            edges[(winner, loser)] = 1

    cycles = 0
    for i, j, k in combinations(range(len(runs)), 3):
        if (
            (i, j) in edges and (j, k) in edges and (k, i) in edges
        ) or (
            (i, k) in edges and (k, j) in edges and (j, i) in edges
        ):
            cycles += 1
    path = sorted(range(len(runs)), key=lambda idx: (-wins[idx], runs[idx].name))
    print("  Tournament Analysis on DP order strategies")
    print("    vertices: insertion schedules / proof-state addresses, not runners")
    print("    observable: smaller max support, then area, transitions, entropy")
    print(f"    score_hist={dict(sorted(Counter(wins.values()).items()))} directed_3cycles={cycles}")
    print("    Hamiltonian path=" + " > ".join(runs[idx].name for idx in path))


def row_bank() -> list[tuple[str, Row]]:
    rows = [
        ("AP9", tuple(range(9))),
        ("one-gap-top9", (0, 1, 2, 3, 4, 5, 6, 7, 9)),
        ("direct-risk-leader", (0, 4, 6, 8, 10, 12, 14, 15, 16)),
        ("boundary-leader", (0, 2, 4, 6, 8, 10, 12, 14, 15)),
        ("dyadic-block", (0, 1, 2, 4, 8, 12, 16, 20, 24)),
        ("doubled-odd-tail", (0, 1, 2, 4, 8, 14, 26, 34, 38)),
        ("three-cluster-wide", (0, 1, 2, 30, 31, 32, 60, 61, 62)),
        ("AP-triple-phase", (0, 9, 10, 11, 12, 13, 14, 15, 16, 17)),
    ]
    return [(label, normalized(row)) for label, row in rows]


def main() -> None:
    print("HYP-2691 / T927 -- LRC14 exact sector-state transfer DP")
    print("Arithmetic: exact Fractions over common sector-wall refinements.")
    print("Preserved predicate: final p0(E) <= cap_k.")
    print("Destroyed data: actual continuous time after quotienting each wall atom to a missed-sector state.")
    print("Assumption challenged: proof vertices need not be runners; here they are DP states and insertion schedules.\n")

    print("EXACT INSERTION RECURRENCE")
    print("  Adding one speed can remove at most one missed inner sector on any wall atom.")
    print("  Therefore delta p0 is exactly the mass of atoms where the prefix misses one")
    print("  sector and the new speed lands in that sector.  The transition matrix is")
    print("  lower triangular in missed-sector inclusion, with rank-one atom moves.\n")

    all_runs: list[OrderRun] = []
    for label, row in row_bank():
        dist = missed_distribution(row)
        cap = CAP.get(len(row))
        margin = cap - dist[0] if cap is not None else None
        print("=" * 90)
        print(f"ROW {label}: E={row}")
        print(
            f"  k={len(row)} span={row[-1]} p0={fmt(dist[0])} "
            + (f"cap={fmt(cap)} margin={fmt(margin)} " if cap is not None and margin is not None else "")
            + f"missed_dist={[str(x) for x in dist]}"
        )
        print(
            f"  additive: excess={sumset_excess(row)} energy={additive_energy(row)} "
            f"squarefree={squarefree_profile(row)}"
        )
        runs = [run_order(row, name, order) for name, order in orders_for(row).items()]
        # Deduplicate identical orders while preserving the first name.
        dedup: dict[tuple[int, ...], OrderRun] = {}
        for run in runs:
            dedup.setdefault(run.order, run)
        runs = list(dedup.values())
        runs.sort(key=lambda r: (r.max_support, r.area_support, r.max_transitions, r.final_entropy, r.name))
        all_runs.extend(runs)

        for run in runs:
            print_order_summary(run)
        tournament_order_analysis(runs)

        best = runs[0]
        print_step_ledger(best)
        if label in {"direct-risk-leader", "dyadic-block", "three-cluster-wide"} and best.name != "increasing":
            increasing = next((r for r in runs if r.name == "increasing"), None)
            if increasing is not None:
                print_step_ledger(increasing)
        print()

    print("=" * 90)
    print("GLOBAL ORDER-STRATEGY SCORECARD")
    by_name: dict[str, list[OrderRun]] = defaultdict(list)
    for run in all_runs:
        by_name[run.name].append(run)
    for name, runs in sorted(by_name.items()):
        avg_max_support = sum(r.max_support for r in runs) / len(runs)
        avg_area = sum(r.area_support for r in runs) / len(runs)
        avg_trans = sum(r.max_transitions for r in runs) / len(runs)
        print(
            f"  {name:<21} rows={len(runs):2d} avg_max_support={avg_max_support:.2f} "
            f"avg_area={avg_area:.2f} avg_max_trans={avg_trans:.2f}"
        )

    print("\nSYNTHESIS")
    print("  1. The exact DP transition is local: one inserted runner only deletes one")
    print("     sector from the missed-state word on each wall atom.")
    print("  2. Order choice changes the intermediate state explosion even though final")
    print("     p0 is invariant.  This is the LRC analogue of the half-tiling address")
    print("     quotient: choose the address schedule before scalarizing.")
    print("  3. Greedy support/transition schedules expose finite low-state pockets")
    print("     for the resonant rows.  High-support prefixes should be the natural")
    print("     input for Weyl/decorrelation, while low-support prefixes route to")
    print("     AP/dyadic/Ruzsa finite atlases.")
    print("  4. The proof obligation suggested by the DP is a transfer inequality:")
    print("     bound each prefix's one-missed-sector landing mass by either finite")
    print("     address templates or a decorrelation estimate on the transition kernel.")


if __name__ == "__main__":
    main()
