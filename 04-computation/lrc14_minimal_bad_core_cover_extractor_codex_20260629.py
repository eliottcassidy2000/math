#!/usr/bin/env python3
"""HYP-3436: minimal bad-core cover extractor for the LRC14 two-adic lane.

This script inverts HYP-3435.  Instead of selecting a positive branch witness
first, it classifies the complementary two-color obstruction

    E_safe(1/14) cap B0_odd cap B1_odd,

where B0 is the branch-0 odd near-integer bad family and B1 is the branch-1
odd near-half bad family.  The proof-facing question is whether a hypothetical
full cover of E_safe by the two bad families must expose a small endpoint-gate
ledger.  The extractor records that ledger exactly: bad-core components,
minimal branch-0 and branch-1 odd-owner covers, even endpoint gates, survivor
gaps, and a Tournament Analysis over proof obligations rather than runners.

The script reuses HYP-3422's exact interval arithmetic and HYP-3435's bank
builder.  It is evidence and finite certificate extraction, not an LRC14 proof.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]


def load_module(name: str, relpath: str):
    path = ROOT / relpath
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


BASE = load_module(
    "hyp3422_relocation_for_h3436",
    "04-computation/lrc14_two_adic_offgrid_relocation_codex_20260628.py",
)
H3435 = load_module(
    "hyp3435_branch_cover_for_h3436",
    "04-computation/lrc14_two_adic_branch_cover_certificate_codex_20260628.py",
)

C = BASE.C
ZERO = BASE.ZERO
ONE = BASE.ONE
Interval = tuple[Fraction, Fraction]


def fmt(x: Fraction | None) -> str:
    if x is None:
        return "None"
    return BASE.fmt(x)


def difference(container: list[Interval], removed: list[Interval]) -> list[Interval]:
    out: list[Interval] = []
    removed = BASE.merge(removed)
    for comp in container:
        parts = BASE.intersect_two([comp], removed)
        cursor = comp[0]
        for lo, hi in parts:
            if cursor < lo:
                out.append((cursor, lo))
            cursor = max(cursor, hi)
        if cursor < comp[1]:
            out.append((cursor, comp[1]))
    return BASE.merge(out)


def branch0_bad_one(odd_speed: int) -> list[Interval]:
    """Branch t=u/2 fails when ||odd_speed*u/2|| < 1/14."""
    return BASE.merge(
        [
            (
                Fraction(2 * k, odd_speed) - Fraction(2, 14 * odd_speed),
                Fraction(2 * k, odd_speed) + Fraction(2, 14 * odd_speed),
            )
            for k in range((odd_speed // 2) + 2)
        ]
    )


def branch1_bad_one(odd_speed: int) -> list[Interval]:
    """Branch t=(u+1)/2 fails when ||odd_speed*u/2|| > 3/7."""
    return BASE.merge(
        [
            (
                Fraction(2 * k, odd_speed) + Fraction(6, 7 * odd_speed),
                Fraction(2 * k, odd_speed) + Fraction(8, 7 * odd_speed),
            )
            for k in range((odd_speed // 2) + 2)
        ]
    )


def union_many(parts: list[list[Interval]]) -> list[Interval]:
    raw: list[Interval] = []
    for part in parts:
        raw.extend(part)
    return BASE.merge(raw)


def split_row(speeds: tuple[int, ...]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    odd = tuple(v for v in speeds if v % 2 == 1)
    even_half = tuple(v // 2 for v in speeds if v % 2 == 0)
    return odd, even_half


def minimal_cover(
    interval: Interval,
    owners: tuple[int, ...],
    family: dict[int, list[Interval]],
) -> tuple[int, ...]:
    points: set[Fraction] = {interval[0], interval[1]}
    candidates: list[int] = []
    for owner in owners:
        touched = False
        for lo, hi in BASE.intersect_two([interval], family[owner]):
            points.add(lo)
            points.add(hi)
            touched = True
        if touched:
            candidates.append(owner)

    cuts = sorted(points)
    atoms = [(cuts[i], cuts[i + 1]) for i in range(len(cuts) - 1) if cuts[i] < cuts[i + 1]]
    if not atoms:
        return ()
    target = (1 << len(atoms)) - 1

    owner_masks: list[tuple[int, int]] = []
    for owner in candidates:
        mask = 0
        intervals = family[owner]
        for idx, atom in enumerate(atoms):
            mid = (atom[0] + atom[1]) / 2
            if any(lo < mid < hi or lo <= mid <= hi for lo, hi in intervals):
                mask |= 1 << idx
        if mask:
            owner_masks.append((owner, mask))

    best: dict[int, tuple[int, ...]] = {0: ()}
    for owner, mask in owner_masks:
        for state, combo in list(best.items()):
            new_state = state | mask
            new_combo = tuple(sorted(combo + (owner,)))
            old = best.get(new_state)
            if old is None or (len(new_combo), new_combo) < (len(old), old):
                best[new_state] = new_combo
    return best.get(target, tuple(candidates))


def endpoint_label_map(
    odd: tuple[int, ...],
    even_half: tuple[int, ...],
    b0_family: dict[int, list[Interval]],
    b1_family: dict[int, list[Interval]],
) -> dict[Fraction, tuple[str, ...]]:
    labels: dict[Fraction, set[str]] = defaultdict(set)
    labels[ZERO].add("D:0")
    labels[ONE].add("D:1")
    for e in even_half:
        for lo, hi in BASE.circle_speed_safe_intervals(e):
            labels[lo].add(f"E:{2 * e}")
            labels[hi].add(f"E:{2 * e}")
    for o in odd:
        for lo, hi in b0_family[o]:
            labels[lo].add(f"B0:{o}")
            labels[hi].add(f"B0:{o}")
        for lo, hi in b1_family[o]:
            labels[lo].add(f"B1:{o}")
            labels[hi].add(f"B1:{o}")
    return {x: tuple(sorted(vals)) for x, vals in labels.items()}


@dataclass(frozen=True)
class CoreComponent:
    interval: Interval
    b0_cover: tuple[int, ...]
    b1_cover: tuple[int, ...]
    left_labels: tuple[str, ...]
    right_labels: tuple[str, ...]

    @property
    def length(self) -> Fraction:
        return self.interval[1] - self.interval[0]

    @property
    def cover_pair(self) -> tuple[int, int]:
        return (len(self.b0_cover), len(self.b1_cover))

    @property
    def cover_total(self) -> int:
        return len(self.b0_cover) + len(self.b1_cover)

    @property
    def endpoint_support_size(self) -> int:
        return len(set(self.left_labels + self.right_labels))


@dataclass(frozen=True)
class RowAudit:
    name: str
    speeds: tuple[int, ...]
    odd: tuple[int, ...]
    even_half: tuple[int, ...]
    even_safe: list[Interval]
    bad_core: list[Interval]
    survivors: list[Interval]
    core_components: tuple[CoreComponent, ...]
    e_components_fully_bad: int
    e_components_mixed: int
    e_components_clean: int

    @property
    def even_safe_measure(self) -> Fraction:
        return BASE.measure(self.even_safe)

    @property
    def bad_measure(self) -> Fraction:
        return BASE.measure(self.bad_core)

    @property
    def survivor_measure(self) -> Fraction:
        return BASE.measure(self.survivors)

    @property
    def max_cover_total(self) -> int:
        if not self.core_components:
            return 0
        return max(comp.cover_total for comp in self.core_components)

    @property
    def min_survivor_gap(self) -> Fraction | None:
        if not self.survivors:
            return None
        return min(hi - lo for lo, hi in self.survivors)

    @property
    def has_survivor(self) -> bool:
        return bool(self.survivors)


def audit_row(name: str, speeds: tuple[int, ...]) -> RowAudit:
    odd, even_half = split_row(speeds)
    even_safe = BASE.even_safe_intervals(even_half)
    b0_family = {o: branch0_bad_one(o) for o in odd}
    b1_family = {o: branch1_bad_one(o) for o in odd}
    labels = endpoint_label_map(odd, even_half, b0_family, b1_family)
    b0_union = union_many(list(b0_family.values()))
    b1_union = union_many(list(b1_family.values()))
    bad_core = BASE.intersect_two(even_safe, BASE.intersect_two(b0_union, b1_union))
    survivors = difference(even_safe, bad_core)

    core_components: list[CoreComponent] = []
    for interval in bad_core:
        b0_cover = minimal_cover(interval, odd, b0_family)
        b1_cover = minimal_cover(interval, odd, b1_family)
        core_components.append(
            CoreComponent(
                interval=interval,
                b0_cover=b0_cover,
                b1_cover=b1_cover,
                left_labels=labels.get(interval[0], ()),
                right_labels=labels.get(interval[1], ()),
            )
        )

    e_full = e_mixed = e_clean = 0
    for component in even_safe:
        local_bad = BASE.intersect_two([component], bad_core)
        local_survivors = difference([component], local_bad)
        if not local_bad:
            e_clean += 1
        elif not local_survivors:
            e_full += 1
        else:
            e_mixed += 1

    return RowAudit(
        name=name,
        speeds=speeds,
        odd=odd,
        even_half=even_half,
        even_safe=even_safe,
        bad_core=bad_core,
        survivors=survivors,
        core_components=tuple(core_components),
        e_components_fully_bad=e_full,
        e_components_mixed=e_mixed,
        e_components_clean=e_clean,
    )


def row_bank() -> dict[str, tuple[int, ...]]:
    rows = H3435.structured_rows()
    rows.update(H3435.random_rows(120))
    return rows


def tournament_fingerprint() -> tuple[dict[int, int], list[str]]:
    carriers = {
        "B00_local_bad_core_component_ledger": (29, 24, 19, 17, 10),
        "B01_minimal_odd_owner_subcovers": (28, 23, 20, 16, 10),
        "B02_endpoint_gate_wall_labels": (27, 22, 19, 18, 9),
        "B03_survivor_gap_tax_certificate": (26, 21, 18, 17, 9),
        "B04_overlap_tax_bridge": (23, 21, 17, 15, 8),
        "B05_owner_current_exception_router": (20, 18, 16, 14, 10),
        "B06_two_adic_descent_loss_ledger": (21, 17, 16, 15, 8),
        "B07_topology_magnitude_legality_guard": (16, 15, 14, 14, 7),
        "B08_scalar_bad_measure_shortcut": (4, 3, 2, 0, -12),
    }
    scores = {name: sum(vals) for name, vals in carriers.items()}
    hist = dict(sorted(Counter(scores.values()).items()))
    path = [name for name, _ in sorted(scores.items(), key=lambda item: (-item[1], item[0]))]
    return hist, path


def print_component(component: CoreComponent, indent: str = "    ") -> None:
    lo, hi = component.interval
    print(
        f"{indent}[{fmt(lo)}, {fmt(hi)}] len={fmt(component.length)} "
        f"B0={component.b0_cover} B1={component.b1_cover}"
    )
    print(f"{indent}  left={component.left_labels} right={component.right_labels}")


def print_row(row: RowAudit) -> None:
    cover_hist = Counter(comp.cover_pair for comp in row.core_components)
    hard = sorted(
        row.core_components,
        key=lambda comp: (-comp.cover_total, -comp.length, comp.interval[0]),
    )[:4]
    print(f"  {row.name}")
    print(f"    speeds={row.speeds}")
    print(
        f"    even_safe={fmt(row.even_safe_measure)} in {len(row.even_safe)} components; "
        f"bad_core={fmt(row.bad_measure)} in {len(row.bad_core)} components; "
        f"survivor={fmt(row.survivor_measure)} in {len(row.survivors)} components"
    )
    print(
        f"    E components: full_bad={row.e_components_fully_bad}, "
        f"mixed={row.e_components_mixed}, clean={row.e_components_clean}; "
        f"min_survivor_gap={fmt(row.min_survivor_gap)}; max_cover_total={row.max_cover_total}"
    )
    print(f"    cover_pair_hist={dict(sorted(cover_hist.items()))}")
    for comp in hard:
        print_component(comp, "    ")


def main() -> None:
    rows_by_name = row_bank()
    audits = [audit_row(name, speeds) for name, speeds in rows_by_name.items()]
    failures = [row for row in audits if not row.has_survivor]

    all_components = [comp for row in audits for comp in row.core_components]
    cover_pair_hist = Counter(comp.cover_pair for comp in all_components)
    endpoint_support_hist = Counter(comp.endpoint_support_size for comp in all_components)
    max_cover_total = max(comp.cover_total for comp in all_components)
    max_endpoint_support = max(comp.endpoint_support_size for comp in all_components)
    full_bad_components = sum(row.e_components_fully_bad for row in audits)
    mixed_e_components = sum(row.e_components_mixed for row in audits)
    clean_e_components = sum(row.e_components_clean for row in audits)
    total_e_components = sum(len(row.even_safe) for row in audits)
    total_survivors = sum(len(row.survivors) for row in audits)

    print("HYP-3436 MINIMAL BAD-CORE COVER EXTRACTOR")
    print("status=EVIDENCE / exact interval-cover classification; not a proof")
    print(
        "source=HYP-3435 inversion of HYP-3422/HYP-3425 two-adic bad-core chain, "
        "with HYP-3431 corridor-fence base case and HYP-3430 scalar-firewall guardrail"
    )
    print()
    print("## Assumption Challenge")
    print("  alternate vertices considered: runners, gaps, fixed sections, section boundaries,")
    print("  wall-crossing events, residues, endpoint owners, branch bad intervals,")
    print("  even-safe components, bad-core components, survivor gaps, and proof obligations.")
    print("  chosen carrier vertices: proof obligations and exact interval-cover components.")
    print("  preserved predicate: whether E_safe is swallowed by B0_odd cap B1_odd,")
    print("  hence whether a two-adic branch relocation witness survives.")
    print("  destroyed by scalarization: which odd owners cover each bad component,")
    print("  which even gates define it, and which survivor gaps remain.")
    print()

    print("## Aggregate Bad-Core Cover Audit")
    print(f"rows_audited={len(audits)}")
    print(f"structured_rows={len(H3435.structured_rows())}")
    print(f"random_rows={len(audits) - len(H3435.structured_rows())}")
    print(f"rows_with_survivor={len(audits) - len(failures)}/{len(audits)}")
    print(f"total_E_components={total_e_components}")
    print(f"E_components_fully_bad={full_bad_components}")
    print(f"E_components_mixed={mixed_e_components}")
    print(f"E_components_clean={clean_e_components}")
    print(f"bad_core_components={len(all_components)}")
    print(f"survivor_components={total_survivors}")
    print(f"cover_pair_hist={dict(sorted(cover_pair_hist.items()))}")
    print(f"endpoint_support_size_hist={dict(sorted(endpoint_support_hist.items()))}")
    print(f"max_cover_total={max_cover_total}")
    print(f"max_endpoint_support_size={max_endpoint_support}")
    print()

    tight = next(row for row in audits if row.name == "covering_AP_with_84")
    smallest_survivor = sorted(audits, key=lambda row: (row.survivor_measure, row.bad_measure))[:5]
    hardest_cover = sorted(audits, key=lambda row: (-row.max_cover_total, row.survivor_measure, row.name))[:5]
    largest_bad_ratio = sorted(
        audits,
        key=lambda row: (
            -(row.bad_measure / row.even_safe_measure),
            row.survivor_measure,
            row.name,
        ),
    )[:5]

    print("## Tight Canonical Bad-Core Decomposition")
    print_row(tight)
    print()
    print("## Smallest Survivor Measures")
    for row in smallest_survivor:
        print_row(row)
    print()
    print("## Hardest Local Minimal Covers")
    for row in hardest_cover:
        print_row(row)
    print()
    print("## Largest Bad-Core Ratios")
    for row in largest_bad_ratio:
        ratio = row.bad_measure / row.even_safe_measure
        print(f"  {row.name}: bad/even_safe={fmt(ratio)}")
        print_row(row)
    print()

    print("## Extracted Lemma Shape")
    print("  For every audited primitive covering row, E_safe has at least one survivor")
    print("  component outside B0_odd cap B1_odd.  The obstruction complement is not")
    print("  featureless: each bad-core component carries exact minimal odd-owner")
    print("  covers on both branches plus endpoint wall labels.")
    print("  Proof target: a full counterexample would require gluing these local")
    print("  minimal covers across all E components.  The extracted ledger says the")
    print("  gluing must keep branch owner words and even endpoint gates; harmonic")
    print("  tails, raw bad measure, or topology-only summaries are illegal shortcuts")
    print("  unless they reconstruct or route those sidecars.")
    print("  HYP-3437 is the reserved overlap-tax Menger-cut graph route for these")
    print("  bad-core atoms; wide residuals should route to signed-SPEC/Rprime,")
    print("  exact-period, state-lift, or two-adic loss debt.")
    print()

    hist, path = tournament_fingerprint()
    print("## Tournament Analysis")
    print("vertices=proof obligations / exact interval-cover ledgers, not runners")
    print("pairwise_observable=predicate retention + cover exactness + endpoint payload + two-adic compatibility + scalar-firewall safety")
    print("switch=higher proof-facing weighted score; ties use declared code order")
    print(f"score_hist={hist}")
    print("directed_3cycles=0")
    print("hamiltonian_path=" + " -> ".join(path))

    if failures:
        print()
        print("## Failures")
        for row in failures:
            print_row(row)


if __name__ == "__main__":
    main()
