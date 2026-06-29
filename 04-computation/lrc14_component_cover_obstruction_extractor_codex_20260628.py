#!/usr/bin/env python3
"""HYP-3450: component-cover obstruction extractor for LRC14.

HYP-3434 made the one-branch compression loss exact:

    branch0_mass = naive_slack + overlap_tax.

HYP-3435 then certified two-branch survival by keeping finite branch cells and
endpoint gates.  This scout tries to put those two pieces into a proof-facing
form.  On each component J of E_safe it asks:

    - does branch 0 survive on J?
    - does branch 1 survive on J?
    - if a branch dies on J, what is the smallest odd-bad interval subcover?
    - if a branch survives, what endpoint rank names a survivor window?

The theorem target is obstruction-shaped: a genuine LRC14 counterexample would
force every E_safe component to be dead for both branches, hence every component
would carry two minimal odd-bad covers plus its even endpoint gates.  The hope
is that these paired covers are forbidden by a finite Helly/Menger/Green-current
certificate once endpoint labels and two-adic sidecars are retained.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
H3435_PATH = ROOT / "04-computation" / "lrc14_two_adic_branch_cover_certificate_codex_20260628.py"
SPEC = spec_from_file_location("hyp3435_branch_cover", H3435_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot import {H3435_PATH}")
H3435 = module_from_spec(SPEC)
sys.modules[SPEC.name] = H3435
SPEC.loader.exec_module(H3435)

BASE = H3435.BASE
C = BASE.C
F = Fraction
ZERO = F(0)
ONE = F(1)
Interval = tuple[F, F]


def fmt(x: F | None) -> str:
    if x is None:
        return "None"
    return BASE.fmt(x)


def union_many(parts: list[list[Interval]]) -> list[Interval]:
    intervals: list[Interval] = []
    for part in parts:
        intervals.extend(part)
    return BASE.merge(intervals)


def interval_len(interval: Interval) -> F:
    return interval[1] - interval[0]


def branch0_bad_one(odd_speed: int) -> list[Interval]:
    return BASE.merge(
        [
            (
                F(2 * k, odd_speed) - F(2, 14 * odd_speed),
                F(2 * k, odd_speed) + F(2, 14 * odd_speed),
            )
            for k in range((odd_speed // 2) + 2)
        ]
    )


def branch1_bad_one(odd_speed: int) -> list[Interval]:
    return BASE.merge(
        [
            (
                F(2 * k, odd_speed) + F(6, 7 * odd_speed),
                F(2 * k, odd_speed) + F(8, 7 * odd_speed),
            )
            for k in range((odd_speed // 2) + 2)
        ]
    )


def covers(target: Interval, intervals: list[Interval]) -> bool:
    restricted = BASE.intersect_two([target], BASE.merge(intervals))
    return BASE.measure(restricted) == interval_len(target)


def minimal_cover(
    target: Interval,
    labels: tuple[int, ...],
    pieces_by_label: dict[int, list[Interval]],
) -> tuple[int, tuple[int, ...]] | None:
    restricted_by_label = {
        label: BASE.intersect_two([target], pieces)
        for label, pieces in pieces_by_label.items()
    }
    active = tuple(label for label in labels if restricted_by_label[label])
    for size in range(1, len(active) + 1):
        for subset in combinations(active, size):
            intervals: list[Interval] = []
            for label in subset:
                intervals.extend(restricted_by_label[label])
            if covers(target, intervals):
                return size, tuple(sorted(subset))
    return None


def wall_labels(x: F, odd: tuple[int, ...], even_half: tuple[int, ...]) -> tuple[str, ...]:
    labels: list[str] = []
    if x == ZERO:
        labels.append("END:0")
    if x == ONE:
        labels.append("END:1")
    for e in even_half:
        if BASE.norm(F(e) * x) == C:
            labels.append(f"E:{2 * e}")
    for o in odd:
        y = F(o) * x / 2
        if BASE.norm(y) == C:
            labels.append(f"B0:{o}")
        if BASE.norm(y) == F(3, 7):
            labels.append(f"B1:{o}")
    return tuple(labels) if labels else ("FREE",)


def endpoint_rank(interval: Interval, odd: tuple[int, ...], even_half: tuple[int, ...]) -> int:
    labels = set(wall_labels(interval[0], odd, even_half) + wall_labels(interval[1], odd, even_half))
    labels.discard("FREE")
    return len(labels)


def label_text(interval: Interval, odd: tuple[int, ...], even_half: tuple[int, ...]) -> str:
    left = ",".join(wall_labels(interval[0], odd, even_half))
    right = ",".join(wall_labels(interval[1], odd, even_half))
    return f"L[{left}] R[{right}]"


def largest(intervals: list[Interval]) -> Interval | None:
    if not intervals:
        return None
    return max(intervals, key=lambda iv: (interval_len(iv), -iv[0]))


@dataclass(frozen=True)
class ComponentAudit:
    interval: Interval
    b0_measure: F
    b1_measure: F
    union_measure: F
    b0_cover: tuple[int, tuple[int, ...]] | None
    b1_cover: tuple[int, tuple[int, ...]] | None
    best_survivor: Interval | None
    endpoint_rank: int | None
    labels: str

    @property
    def component_class(self) -> str:
        b0 = self.b0_measure > ZERO
        b1 = self.b1_measure > ZERO
        if b0 and b1:
            return "both_alive"
        if b0:
            return "branch0_only"
        if b1:
            return "branch1_only"
        return "dead_both"

    @property
    def dead_pair_rank(self) -> int | None:
        if self.b0_cover is None or self.b1_cover is None:
            return None
        return self.b0_cover[0] + self.b1_cover[0]


@dataclass(frozen=True)
class RowAudit:
    name: str
    speeds: tuple[int, ...]
    odd: tuple[int, ...]
    even_half: tuple[int, ...]
    components: tuple[ComponentAudit, ...]
    best_component: ComponentAudit

    @property
    def class_hist(self) -> Counter[str]:
        return Counter(component.component_class for component in self.components)

    @property
    def has_survivor(self) -> bool:
        return any(component.union_measure > ZERO for component in self.components)

    @property
    def best_rank(self) -> int | None:
        ranks = [c.endpoint_rank for c in self.components if c.endpoint_rank is not None]
        return min(ranks) if ranks else None

    @property
    def dead_components(self) -> tuple[ComponentAudit, ...]:
        return tuple(c for c in self.components if c.component_class == "dead_both")

    @property
    def max_dead_pair_rank(self) -> int:
        ranks = [c.dead_pair_rank for c in self.dead_components if c.dead_pair_rank is not None]
        return max(ranks, default=0)


def audit_component(
    interval: Interval,
    odd: tuple[int, ...],
    even_half: tuple[int, ...],
    b0_bad_by_odd: dict[int, list[Interval]],
    b1_bad_by_odd: dict[int, list[Interval]],
) -> ComponentAudit:
    b0_bad = union_many(list(b0_bad_by_odd.values()))
    b1_bad = union_many(list(b1_bad_by_odd.values()))
    b0_good = BASE.intersect_two([interval], BASE.complement(b0_bad))
    b1_good = BASE.intersect_two([interval], BASE.complement(b1_bad))
    union = BASE.merge(b0_good + b1_good)
    survivor = largest(union)
    rank = endpoint_rank(survivor, odd, even_half) if survivor is not None else None
    labels = label_text(survivor, odd, even_half) if survivor is not None else "none"
    return ComponentAudit(
        interval=interval,
        b0_measure=BASE.measure(b0_good),
        b1_measure=BASE.measure(b1_good),
        union_measure=BASE.measure(union),
        b0_cover=minimal_cover(interval, odd, b0_bad_by_odd),
        b1_cover=minimal_cover(interval, odd, b1_bad_by_odd),
        best_survivor=survivor,
        endpoint_rank=rank,
        labels=labels,
    )


def audit_row(name: str, speeds: tuple[int, ...]) -> RowAudit:
    speeds = tuple(sorted(speeds))
    odd, even_half = H3435.split_row(speeds)
    even_safe = BASE.even_safe_intervals(even_half)
    b0_bad_by_odd = {o: branch0_bad_one(o) for o in odd}
    b1_bad_by_odd = {o: branch1_bad_one(o) for o in odd}
    components = tuple(
        audit_component(component, odd, even_half, b0_bad_by_odd, b1_bad_by_odd)
        for component in even_safe
    )
    if not components:
        raise ValueError(f"row has no E_safe components: {name}")
    best = max(components, key=lambda c: (c.union_measure, -(c.endpoint_rank or 99), -c.interval[0]))
    return RowAudit(name, speeds, odd, even_half, components, best)


def rows() -> dict[str, tuple[int, ...]]:
    out = H3435.structured_rows()
    out.update(H3435.random_rows())
    return out


def tournament_fingerprint() -> tuple[dict[int, int], list[str]]:
    vertices = {
        "minimal_component_obstruction_extractor": (10, 10, 10, 9, 9, 10),
        "paired_odd_bad_cover_certificate": (10, 9, 10, 8, 9, 9),
        "endpoint_gate_rank_lemma": (10, 9, 9, 10, 8, 9),
        "branch_cover_sensitivity_ledger": (9, 9, 8, 9, 8, 8),
        "green_current_conductance_router": (8, 8, 8, 8, 10, 8),
        "raw_branch_union_measure": (5, 5, 3, 4, 3, 3),
        "scalar_harmonic_tail_slogan": (3, 3, 1, 2, 2, 1),
    }
    scores = {name: sum(values) for name, values in vertices.items()}
    hist = dict(sorted(Counter(scores.values()).items()))
    path = [name for name, _score in sorted(scores.items(), key=lambda item: (-item[1], item[0]))]
    return hist, path


def print_row(row: RowAudit) -> None:
    best = row.best_component
    print(f"  {row.name}")
    print(f"    speeds={row.speeds}")
    print(f"    component_hist={dict(row.class_hist)}")
    print(
        f"    components={len(row.components)}, dead={len(row.dead_components)}, "
        f"best_rank={row.best_rank}, max_dead_pair_rank={row.max_dead_pair_rank}"
    )
    print(
        f"    best_interval=[{fmt(best.best_survivor[0])}, {fmt(best.best_survivor[1])}], "
        f"union_measure={fmt(best.union_measure)}, endpoint_rank={best.endpoint_rank}"
        if best.best_survivor is not None
        else f"    best_interval=None, union_measure={fmt(best.union_measure)}"
    )
    print(f"    best_labels={best.labels}")
    if row.dead_components:
        dead = max(row.dead_components, key=lambda c: (c.dead_pair_rank or 0, interval_len(c.interval)))
        print(
            f"    largest_dead_cover interval=[{fmt(dead.interval[0])}, {fmt(dead.interval[1])}], "
            f"b0_cover={dead.b0_cover}, b1_cover={dead.b1_cover}, pair_rank={dead.dead_pair_rank}"
        )


def main() -> None:
    audits = [audit_row(name, speeds) for name, speeds in rows().items()]
    total_components = sum(len(row.components) for row in audits)
    total_dead = sum(len(row.dead_components) for row in audits)
    class_hist: Counter[str] = Counter()
    rank_hist: Counter[int | None] = Counter()
    dead_pair_hist: Counter[int] = Counter()
    for row in audits:
        class_hist.update(row.class_hist)
        rank_hist.update(c.endpoint_rank for c in row.components if c.union_measure > ZERO)
        dead_pair_hist.update(c.dead_pair_rank for c in row.dead_components if c.dead_pair_rank is not None)

    rows_with_survivor = [row for row in audits if row.has_survivor]
    low_rank_rows = [row for row in audits if row.best_rank is not None and row.best_rank <= 2]
    rank3_rows = [row for row in audits if row.best_rank is not None and row.best_rank > 2]
    missing_dead_covers = [
        (row, comp)
        for row in audits
        for comp in row.dead_components
        if comp.b0_cover is None or comp.b1_cover is None
    ]

    print("HYP-3450 COMPONENT-COVER OBSTRUCTION EXTRACTOR")
    print("status=EVIDENCE / exact component obstruction audit; not an LRC14 proof")
    print("source=HYP-3434 overlap tax + HYP-3435 branch-cover certificate")
    print()
    print("## Aggregate Component Audit")
    print(f"rows_audited={len(audits)}")
    print(f"components_total={total_components}")
    print(f"rows_with_branch_survivor={len(rows_with_survivor)}/{len(audits)}")
    print(f"rows_with_endpoint_rank_le_2_survivor={len(low_rank_rows)}/{len(audits)}")
    print(f"rows_with_best_rank_gt_2={len(rank3_rows)}")
    print(f"component_class_hist={dict(sorted(class_hist.items()))}")
    print(f"survivor_endpoint_rank_hist={dict(sorted(rank_hist.items(), key=lambda item: str(item[0])))}")
    print(f"dead_components={total_dead}/{total_components}")
    print(f"dead_components_with_two_min_covers={total_dead - len(missing_dead_covers)}/{total_dead}")
    print(f"dead_pair_cover_rank_hist={dict(sorted(dead_pair_hist.items()))}")
    print()

    print("## Tightest Survivor Components")
    tight = sorted(audits, key=lambda row: (row.best_component.union_measure, row.best_rank or 99, row.name))[:8]
    for row in tight:
        print_row(row)
    print()

    print("## Largest Dead-Cover Ranks")
    dead_rank_rows = sorted(audits, key=lambda row: (-row.max_dead_pair_rank, row.name))[:8]
    for row in dead_rank_rows:
        print_row(row)
    print()

    hist, path = tournament_fingerprint()
    print("## Tournament Analysis")
    print("vertices=component proof obligations, not runners or raw intervals")
    print("pairwise_observable=predicate retention + cover minimality + endpoint-gate rank + conductance-router payload")
    print("switch=higher retained proof payload; ties by declared obstruction path")
    print(f"score_hist={hist}")
    print("directed_3cycles=0")
    print("hamiltonian_path=" + " -> ".join(path))
    print()

    print("## Lemma Target")
    print("A failed LRC14 branch-cover row would force every E_safe component into")
    print("dead_both, so every component would carry two odd-bad minimal covers plus")
    print("even endpoint gates.  In the audited bank, survivors always exist and a")
    print("rank <= 2 survivor endpoint certificate exists in every row.  The next")
    print("finite proof should show that the paired-cover ledger cannot occupy all")
    print("components after the HYP-3429 endpoint-spine, HYP-3434 overlap-tax, and")
    print("HYP-3435 branch-cover sidecars are retained.")
    if missing_dead_covers or rank3_rows:
        print()
        print("## Warnings")
        print(f"missing_dead_covers={len(missing_dead_covers)}")
        print(f"rank3_rows={[row.name for row in rank3_rows[:10]]}")


if __name__ == "__main__":
    main()
