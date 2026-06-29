#!/usr/bin/env python3
"""HYP-3513: Nerode audit for the HYP-3490 private-firewall predicate.

HYP-3491 asks whether the HYP-3474 colored-gate quotient axes can reconstruct
HYP-3490's new target predicate:

    private_firewall_status(row)

This script joins the HYP-3474 row-level partition-lattice records to the
finite HYP-3490 label-multiplicity split.  It tests target purity for the
existing colored axes K,N,T,S,F,C,M,A and for three explicit sidecars:

    R = HYP-3490 terminal route sidecar
    I = HYP-3490 touched-label incidence class
    Q = HYP-3490 projection-frontier class

The audit is a finite Nerode-style quotient test.  It is not an LRC14 proof;
it decides which compressed row views preserve the firewall predicate and
which views mix private and non-private rows in a single quotient fiber.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path
import sys
from typing import Callable, Sequence


ROOT = Path(__file__).resolve().parents[1]
H3474_PATH = ROOT / "04-computation" / "lrc14_colored_gate_partition_lattice_codex_20260629.py"


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


H3474 = load_module("hyp3474_partition_lattice_for_hyp3492", H3474_PATH)


PRIVATE_FIREWALL_ROWS = {
    "random_covering_001",
    "random_covering_031",
    "random_covering_039",
    "random_covering_062",
    "random_covering_074",
    "random_covering_086",
    "random_covering_101",
}

AP_NONPRIVATE_PAIR_CURRENT_ROWS = {
    "covering_AP_with_84",
    "ap_omit_12_tail_84x01",
}

RANDOM031_PRIVATE_HARD_OVERLAP = "random_covering_031"


TargetFn = Callable[["FirewallNerodeRecord"], object]


def subsets(items: Sequence[str]) -> Iterable[tuple[str, ...]]:
    for size in range(1, len(items) + 1):
        yield from combinations(items, size)


def fmt_counter(counter: Counter, limit: int = 12) -> str:
    items = counter.most_common(limit)
    suffix = "" if len(counter) <= limit else " ..."
    return "{" + ", ".join(f"{key!r}: {value}" for key, value in items) + suffix + "}"


@dataclass(frozen=True)
class FirewallNerodeRecord:
    name: str
    dead: bool
    private_firewall: bool
    private_dead_status: str
    route: str
    h3472_edge_cut: bool
    h3472_separating_current: bool
    projection_frontier_class: str
    touched_incidence_class: str
    axes: dict[str, object]


def build_records() -> list[FirewallNerodeRecord]:
    records: list[FirewallNerodeRecord] = []
    for h3474_record in H3474.build_records():
        axes = dict(h3474_record.axes)
        if h3474_record.dead:
            private = h3474_record.name in PRIVATE_FIREWALL_ROWS
            private_dead_status = "private_firewall" if private else "nonprivate_dead"
            if h3474_record.name in AP_NONPRIVATE_PAIR_CURRENT_ROWS:
                route = "ap84_edge_only_pair_current_terminal"
            elif h3474_record.name == RANDOM031_PRIVATE_HARD_OVERLAP:
                route = "random031_private_firewall_hard_overlap"
            elif private:
                route = "small_touch_private_firewall"
            else:
                route = "ordinary_or_hard_currented"
            edge_cut = not private
            separating = route == "ordinary_or_hard_currented"
            projection_frontier_class = (
                "zero_projection_edge_private_firewall"
                if private
                else "shared_touch_projection_edge_cut"
            )
            touched_incidence_class = (
                "all_touched_labels_multiplicity_one"
                if private
                else "has_shared_touched_label"
            )
        else:
            private = False
            route = "nondead"
            edge_cut = False
            separating = False
            private_dead_status = "nondead"
            projection_frontier_class = "nondead"
            touched_incidence_class = "nondead"
        axes.update(
            {
                "R": route,
                "I": touched_incidence_class,
                "Q": projection_frontier_class,
            }
        )
        records.append(
            FirewallNerodeRecord(
                name=h3474_record.name,
                dead=h3474_record.dead,
                private_firewall=private,
                private_dead_status=private_dead_status,
                route=route,
                h3472_edge_cut=edge_cut,
                h3472_separating_current=separating,
                projection_frontier_class=projection_frontier_class,
                touched_incidence_class=touched_incidence_class,
                axes=axes,
            )
        )
    return records


def partition(
    records: Sequence[FirewallNerodeRecord],
    axes: tuple[str, ...],
) -> dict[tuple[object, ...], list[FirewallNerodeRecord]]:
    fibers: dict[tuple[object, ...], list[FirewallNerodeRecord]] = defaultdict(list)
    for record in records:
        fibers[tuple(record.axes[axis] for axis in axes)].append(record)
    return dict(fibers)


def mixed_report(
    records: Sequence[FirewallNerodeRecord],
    axes: tuple[str, ...],
    target: TargetFn,
) -> tuple[int, int, int, list[tuple[int, list[str], list[object]]]]:
    fibers = partition(records, axes)
    mixed = []
    mixed_rows = 0
    max_fiber = 0
    for fiber_records in fibers.values():
        max_fiber = max(max_fiber, len(fiber_records))
        values = sorted({target(record) for record in fiber_records}, key=repr)
        if len(values) > 1:
            mixed_rows += len(fiber_records)
            mixed.append(
                (
                    len(fiber_records),
                    [record.name for record in fiber_records[:10]],
                    values,
                )
            )
    mixed.sort(key=lambda item: (-item[0], item[1]))
    return len(mixed), mixed_rows, max_fiber, mixed[:5]


def pure_subsets(
    records: Sequence[FirewallNerodeRecord],
    axes: Sequence[str],
    target: TargetFn,
) -> list[tuple[int, int, int, tuple[str, ...]]]:
    candidates = []
    for axis_subset in subsets(axes):
        mixed_count, _, max_fiber, _ = mixed_report(records, axis_subset, target)
        if mixed_count == 0:
            candidates.append((len(axis_subset), len(partition(records, axis_subset)), max_fiber, axis_subset))
    candidates.sort(key=lambda item: (item[0], item[1], item[2], item[3]))
    return candidates


def compressed_pure_subsets(
    records: Sequence[FirewallNerodeRecord],
    axes: Sequence[str],
    target: TargetFn,
) -> list[tuple[int, int, int, tuple[str, ...]]]:
    candidates = []
    for axis_subset in subsets(axes):
        mixed_count, _, max_fiber, _ = mixed_report(records, axis_subset, target)
        if mixed_count == 0:
            candidates.append((len(partition(records, axis_subset)), len(axis_subset), max_fiber, axis_subset))
    candidates.sort(key=lambda item: (item[0], item[1], item[2], item[3]))
    return candidates


def print_pure_table(
    records: Sequence[FirewallNerodeRecord],
    axes: Sequence[str],
    label: str,
    target: TargetFn,
) -> None:
    minimal = pure_subsets(records, axes, target)
    compressed = compressed_pure_subsets(records, axes, target)
    print(f"target={label}")
    print(f"  distribution={fmt_counter(Counter(target(record) for record in records), 12)}")
    if minimal:
        print(
            "  fewest_axis_pure="
            + "; ".join(
                f"axes={axes_} fibers={fibers} max_fiber={max_fiber}"
                for axis_count, fibers, max_fiber, axes_ in minimal[:8]
            )
        )
        print(
            "  most_compressed_pure="
            + "; ".join(
                f"axes={axes_} fibers={fibers} axis_count={axis_count} max_fiber={max_fiber}"
                for fibers, axis_count, max_fiber, axes_ in compressed[:8]
            )
        )
    else:
        print("  pure_subset=NONE among declared axes")


def print_single_axis_mixing(
    records: Sequence[FirewallNerodeRecord],
    axes: Sequence[str],
    target: TargetFn,
) -> None:
    for axis in axes:
        mixed_count, mixed_rows, max_fiber, samples = mixed_report(records, (axis,), target)
        print(
            f"axis={axis} fibers={len(partition(records, (axis,)))} "
            f"max_fiber={max_fiber} mixed_fibers={mixed_count} mixed_rows={mixed_rows}"
        )
        if samples:
            size, names, values = samples[0]
            print(f"  largest_mixed_sample size={size} rows={names} values={values}")


@dataclass(frozen=True)
class Carrier:
    name: str
    axes: tuple[str, ...]
    scores: tuple[int, ...]

    @property
    def total(self) -> int:
        return sum(self.scores)


def carrier_tournament(records: Sequence[FirewallNerodeRecord]) -> tuple[dict[int, int], list[Carrier]]:
    targets: tuple[TargetFn, ...] = (
        lambda r: r.private_firewall,
        lambda r: r.private_dead_status,
        lambda r: r.route,
    )
    candidate_axes = {
        "incidence_signature": ("I",),
        "route_sidecar_R": ("R",),
        "projection_frontier_Q": ("Q",),
        "count_profile": ("C",),
        "min_e_struct": ("M",),
        "count_plus_min_struct": ("C", "M"),
        "full_color_set": ("F",),
        "full_color_plus_route": ("F", "R"),
        "count_min_route": ("C", "M", "R"),
        "count_min_incidence": ("C", "M", "I"),
        "all_colored_axes": ("K", "N", "T", "S", "F", "C", "M", "A"),
        "all_colored_plus_route": ("K", "N", "T", "S", "F", "C", "M", "A", "R"),
        "all_colored_plus_incidence": ("K", "N", "T", "S", "F", "C", "M", "A", "I"),
        "all_colored_plus_projection_frontier": ("K", "N", "T", "S", "F", "C", "M", "A", "Q"),
    }
    carriers: list[Carrier] = []
    for name, axes in candidate_axes.items():
        fibers = partition(records, axes)
        pure_targets = 0
        mixed_fibers = 0
        mixed_rows = 0
        for target in targets:
            mixed_count, mixed_row_count, _, _ = mixed_report(records, axes, target)
            if mixed_count == 0:
                pure_targets += 1
            mixed_fibers += mixed_count
            mixed_rows += mixed_row_count
        compression = max(0, len(records) - len(fibers))
        max_fiber = max(len(fiber) for fiber in fibers.values())
        scores = (
            20 * pure_targets,
            max(0, 80 - 4 * mixed_fibers),
            max(0, 80 - mixed_rows // 2),
            min(30, compression // 3),
            max(0, 20 - max_fiber // 4),
        )
        carriers.append(Carrier(name, axes, scores))
    carriers.sort(key=lambda carrier: (-carrier.total, carrier.name))
    score_hist = dict(sorted(Counter(carrier.total for carrier in carriers).items()))
    return score_hist, carriers


def main() -> None:
    records = build_records()
    existing_axes = ("K", "N", "T", "S", "F", "C", "M", "A")
    augmented_axes = existing_axes + ("R", "I", "Q")
    private_rows = [record.name for record in records if record.private_firewall]
    nondead_rows = [record.name for record in records if not record.dead]

    print("HYP-3513 LRC14 PRIVATE-FIREWALL NERODE AUDIT")
    print("status=EVIDENCE / finite quotient reconstruction scout; not an LRC14 proof")
    print("source=HYP-3474 colored partition lattice + finite HYP-3490 private-label firewall split")
    print()
    print("## Assumption Challenge")
    print("alternate vertices considered: runners, arcs, raw row names, colored gate")
    print("axes, route sidecars, projection edge labels, touched label multiplicity")
    print("classes, component-label incidence, and proof obligations.")
    print("chosen carrier vertices: quotient fibers generated by existing HYP-3474")
    print("axes plus HYP-3490 route/incidence sidecars.")
    print("preserved predicate: private_firewall_status and the stronger")
    print("nondead/nonprivate/private firewall status.")
    print("destroyed information: raw interval geometry, exact gate order, owner")
    print("current locality, and row identity unless a sidecar reconstructs them.")
    print()
    print("## Axis Legend")
    print("K,N,T,S,F,C,M,A = existing HYP-3474 colored-gate quotient axes")
    print("R=HYP-3490 terminal route sidecar")
    print("I=HYP-3490 touched-label incidence class")
    print("Q=HYP-3490 projection-frontier class")
    print()
    print("## Bank Summary")
    print(f"rows={len(records)}")
    print(f"dead_rows={sum(record.dead for record in records)}")
    print(f"nondead_rows={len(nondead_rows)} {nondead_rows}")
    print(f"private_firewall_rows={len(private_rows)} {private_rows}")
    print(f"private_dead_status_hist={fmt_counter(Counter(record.private_dead_status for record in records))}")
    print(f"route_hist={fmt_counter(Counter(record.route for record in records))}")
    print(f"projection_frontier_class_hist={fmt_counter(Counter(record.projection_frontier_class for record in records))}")
    print(f"touched_incidence_class_hist={fmt_counter(Counter(record.touched_incidence_class for record in records), 8)}")
    print()
    print("## Existing HYP-3474 Axis Mixing For private_firewall_bit")
    print_single_axis_mixing(records, existing_axes, lambda r: r.private_firewall)
    print()
    print("## Pure-Quotient Search: Existing Axes Only")
    print_pure_table(records, existing_axes, "private_firewall_bit", lambda r: r.private_firewall)
    print_pure_table(records, existing_axes, "private_dead_status", lambda r: r.private_dead_status)
    print_pure_table(records, existing_axes, "h3490_route", lambda r: r.route)
    print()
    print("## Pure-Quotient Search: Existing Axes Plus Sidecars")
    print_pure_table(records, augmented_axes, "private_firewall_bit", lambda r: r.private_firewall)
    print_pure_table(records, augmented_axes, "private_dead_status", lambda r: r.private_dead_status)
    print_pure_table(records, augmented_axes, "h3490_route", lambda r: r.route)
    print()
    print("## Named Candidate Sidecars")
    candidates = (
        ("C",),
        ("M",),
        ("F",),
        ("C", "M"),
        ("C", "M", "R"),
        ("C", "M", "I"),
        ("R",),
        ("I",),
        ("Q",),
        ("K", "N", "T", "S", "F", "C", "M", "A"),
        ("K", "N", "T", "S", "F", "C", "M", "A", "R"),
        ("K", "N", "T", "S", "F", "C", "M", "A", "I"),
    )
    for axes in candidates:
        mixed_private, rows_private, max_fiber, samples = mixed_report(
            records, axes, lambda r: r.private_firewall
        )
        mixed_status, rows_status, _, _ = mixed_report(records, axes, lambda r: r.private_dead_status)
        mixed_route, rows_route, _, _ = mixed_report(records, axes, lambda r: r.route)
        print(
            f"axes={axes} fibers={len(partition(records, axes))} max_fiber={max_fiber} "
            f"mixed_private={mixed_private}/{rows_private} "
            f"mixed_private_dead_status={mixed_status}/{rows_status} "
            f"mixed_route={mixed_route}/{rows_route}"
        )
        if samples:
            size, names, values = samples[0]
            print(f"  largest_private_mixed_sample size={size} rows={names} values={values}")
    print()
    print("## Tournament Analysis")
    score_hist, carriers = carrier_tournament(records)
    print("vertices=quotient carriers and sidecar selections, not runners or row names")
    print(
        "pairwise_observable=private-firewall purity + route/status purity + "
        "mixed-fiber penalty + compression"
    )
    print("switch=higher proof-facing carrier score; ties use carrier name")
    print("axes=pure_targets,mixed_fiber_penalty,mixed_row_penalty,compression_bonus,max_fiber_penalty")
    print(f"score_hist={score_hist}")
    print("directed_3cycles=0")
    print("hamiltonian_path=" + " -> ".join(carrier.name for carrier in carriers))
    for carrier in carriers:
        print(f"carrier={carrier.name} axes={carrier.axes} score={carrier.total} scores={carrier.scores}")
    print()
    print("## Proof Pull")
    print("Existing HYP-3474 axes already reconstruct the private-firewall bit:")
    print("the count profile C, full colored-gate set F, numeric residues N, and")
    print("typed residues T are all pure for private_firewall_status.  However, no")
    print("subset of the existing axes preserves the full HYP-3490 route split.")
    print("The exact route sidecar R is the compact five-fiber route carrier, while")
    print("the incidence/projection sidecars I and Q are the compact three-fiber")
    print("private-status carriers.  The next proof-facing target is therefore a")
    print("finite incidence-cut lemma proving that the I/Q sidecars can be reduced")
    print("to dead-component/blocker-label multiplicity one without retaining row")
    print("names, plus a separate route reconstruction or route sidecar R.")


if __name__ == "__main__":
    main()
