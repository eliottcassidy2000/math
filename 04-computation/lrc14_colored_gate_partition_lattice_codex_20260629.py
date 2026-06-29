#!/usr/bin/env python3
"""HYP-3474: partition-lattice audit for LRC14 colored gate quotients.

HYP-3471 says the current finite bank has a strong colored-gate fact:

    dead_components(row) > 0
      => row has a rank <= 2 E/branch survivor gate.

This script asks the next quotient question.  Which compressed row views are
allowed to forget data while still preserving proof-facing predicates, and
which views mix AP, random, same-branch, cross-branch, and E/branch-only route
types in the same fiber?

This is a finite Nerode-style scout: a quotient is legal for a target predicate
only when every quotient fiber is target-pure.  It is not a proof of LRC14, but
it turns the "quotient may forget only if fibers are pure" guardrail into an
executable object on the HYP-3471 colored-gate bank.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path
import sys
from typing import Callable, Iterable, Sequence


ROOT = Path(__file__).resolve().parents[1]
H3471_PATH = ROOT / "04-computation" / "lrc14_colored_gate_reservoir_codex_20260629.py"


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


H3471 = load_module("hyp3471_colored_gate_reservoir_for_hyp3472", H3471_PATH)


CANONICAL_AP_PALETTE = {
    "B1:7|E:0",
    "E:0|B0:7",
    "E:0|B1:5",
    "B0:5|E:0",
}


def keyset(values: Iterable[object]) -> tuple[object, ...]:
    return tuple(sorted(set(values), key=repr))


def row_numeric_set(row) -> tuple[str, ...]:
    return keyset(H3471.numeric_word(gate) for gate in row.low_rank_gates)  # type: ignore[return-value]


def row_full_set(row) -> tuple[tuple[str, str, str, int, int], ...]:
    return keyset(H3471.full_color_word(gate) for gate in row.low_rank_gates)  # type: ignore[return-value]


def row_counts(row) -> tuple[int, int, int, int, int, int]:
    e_count = sum(1 for gate in row.low_rank_gates if H3471.is_e_branch_gate(gate))
    same_count = sum(1 for gate in row.low_rank_gates if H3471.is_same_branch_gate(gate))
    cross_count = sum(1 for gate in row.low_rank_gates if H3471.is_cross_branch_gate(gate))
    return (
        len(row.low_rank_gates),
        e_count,
        same_count,
        cross_count,
        len(H3471.row_typed_set(row)),
        len(H3471.row_structural_set(row)),
    )


def best_e_branch_gate(row):
    gates = [gate for gate in row.low_rank_gates if H3471.is_e_branch_gate(gate)]
    if not gates:
        return None
    return min(gates, key=lambda g: (g.length, g.row_name, g.component_index, g.interval[0]))


def min_e_typed(row) -> str:
    gate = best_e_branch_gate(row)
    return "none" if gate is None else H3471.typed_word(gate)


def min_e_struct(row) -> tuple[str, str, str, int, int] | str:
    gate = best_e_branch_gate(row)
    return "none" if gate is None else H3471.structural_word(gate)


def min_e_full(row) -> tuple[str, tuple[str, str, str, int, int]] | str:
    gate = best_e_branch_gate(row)
    if gate is None:
        return "none"
    return (H3471.typed_word(gate), H3471.structural_word(gate))


@dataclass(frozen=True)
class RowRecord:
    name: str
    dead: bool
    has_e_branch: bool
    has_same_branch: bool
    has_cross_branch: bool
    only_e_branch: bool
    has_ap84_packet: bool
    min_e_typed: str
    min_e_struct: tuple[str, str, str, int, int] | str
    min_e_full: tuple[str, tuple[str, str, str, int, int]] | str
    axes: dict[str, object]

    @property
    def route_flags(self) -> tuple[bool, bool, bool, bool, bool]:
        return (
            self.dead,
            self.has_ap84_packet,
            self.has_same_branch,
            self.has_cross_branch,
            self.only_e_branch,
        )

    @property
    def theorem_failure(self) -> bool:
        return self.dead and not self.has_e_branch

    @property
    def dead_min_struct_target(self) -> object:
        return "nondead" if not self.dead else self.min_e_struct

    @property
    def dead_min_full_target(self) -> object:
        return "nondead" if not self.dead else self.min_e_full


def build_records() -> list[RowRecord]:
    rows = [H3471.H3453.join_row(name, speeds) for name, speeds in H3471.H3453.H3450.rows().items()]
    records: list[RowRecord] = []
    for row in rows:
        has_e = any(H3471.is_e_branch_gate(gate) for gate in row.low_rank_gates)
        has_same = any(H3471.is_same_branch_gate(gate) for gate in row.low_rank_gates)
        has_cross = any(H3471.is_cross_branch_gate(gate) for gate in row.low_rank_gates)
        only_e = bool(row.low_rank_gates) and all(
            H3471.is_e_branch_gate(gate) for gate in row.low_rank_gates
        )
        has_ap = any(H3471.typed_word(gate) in CANONICAL_AP_PALETTE for gate in row.low_rank_gates)
        axes = {
            "K": H3471.row_kind_set(row),
            "N": row_numeric_set(row),
            "T": H3471.row_typed_set(row),
            "S": H3471.row_structural_set(row),
            "F": row_full_set(row),
            "C": row_counts(row),
            "M": min_e_struct(row),
            "A": has_ap,
        }
        records.append(
            RowRecord(
                name=row.name,
                dead=row.has_dead,
                has_e_branch=has_e,
                has_same_branch=has_same,
                has_cross_branch=has_cross,
                only_e_branch=only_e,
                has_ap84_packet=has_ap,
                min_e_typed=min_e_typed(row),
                min_e_struct=min_e_struct(row),
                min_e_full=min_e_full(row),
                axes=axes,
            )
        )
    return records


def subsets(items: Sequence[str]) -> Iterable[tuple[str, ...]]:
    for size in range(1, len(items) + 1):
        yield from combinations(items, size)


def partition(records: Sequence[RowRecord], axes: tuple[str, ...]) -> dict[tuple[object, ...], list[RowRecord]]:
    fibers: dict[tuple[object, ...], list[RowRecord]] = defaultdict(list)
    for record in records:
        fibers[tuple(record.axes[axis] for axis in axes)].append(record)
    return dict(fibers)


TargetFn = Callable[[RowRecord], object]


def mixed_report(
    records: Sequence[RowRecord],
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
                    [record.name for record in fiber_records[:6]],
                    values[:6],
                )
            )
    mixed.sort(key=lambda item: (-item[0], item[1]))
    return len(mixed), mixed_rows, max_fiber, mixed[:3]


def target_distribution(records: Sequence[RowRecord], target: TargetFn) -> Counter:
    return Counter(target(record) for record in records)


def best_pure_subsets(
    records: Sequence[RowRecord],
    axis_names: Sequence[str],
    target: TargetFn,
    limit: int = 8,
) -> list[tuple[int, int, int, tuple[str, ...]]]:
    candidates: list[tuple[int, int, int, tuple[str, ...]]] = []
    for axes in subsets(axis_names):
        mixed_count, _, max_fiber, _ = mixed_report(records, axes, target)
        if mixed_count == 0:
            fiber_count = len(partition(records, axes))
            candidates.append((len(axes), fiber_count, max_fiber, axes))
    candidates.sort(key=lambda item: (item[0], item[1], item[2], item[3]))
    return candidates[:limit]


def best_compressed_pure_subsets(
    records: Sequence[RowRecord],
    axis_names: Sequence[str],
    target: TargetFn,
    limit: int = 8,
) -> list[tuple[int, int, int, tuple[str, ...]]]:
    candidates: list[tuple[int, int, int, tuple[str, ...]]] = []
    for axes in subsets(axis_names):
        mixed_count, _, max_fiber, _ = mixed_report(records, axes, target)
        if mixed_count == 0:
            fiber_count = len(partition(records, axes))
            candidates.append((fiber_count, len(axes), max_fiber, axes))
    candidates.sort(key=lambda item: (item[0], item[1], item[2], item[3]))
    return candidates[:limit]


@dataclass(frozen=True)
class Carrier:
    name: str
    axes: tuple[str, ...]
    scores: tuple[int, ...]

    @property
    def total(self) -> int:
        return sum(self.scores)


def tournament(records: Sequence[RowRecord]) -> tuple[dict[int, int], list[Carrier]]:
    targets: tuple[TargetFn, ...] = (
        lambda r: r.dead,
        lambda r: r.has_ap84_packet,
        lambda r: r.has_same_branch,
        lambda r: r.has_cross_branch,
        lambda r: r.only_e_branch,
        lambda r: r.route_flags,
    )
    carriers: list[Carrier] = []
    candidate_axes = {
        "endpoint_kind_set": ("K",),
        "numeric_mod14_set": ("N",),
        "typed_mod14_set": ("T",),
        "structural_set": ("S",),
        "full_color_set": ("F",),
        "count_profile": ("C",),
        "min_e_branch_struct": ("M",),
        "ap84_presence_bit": ("A",),
        "typed_plus_min_struct": ("T", "M"),
        "structural_plus_ap_bit": ("S", "A"),
        "full_plus_min_struct": ("F", "M"),
        "all_color_sidecars": ("K", "T", "S", "F", "M", "A"),
    }
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
        compression = max(0, 135 - len(fibers))
        max_fiber = max(len(fiber) for fiber in fibers.values())
        scores = (
            12 * pure_targets,
            max(0, 60 - mixed_fibers),
            max(0, 60 - mixed_rows // 4),
            min(20, compression // 2),
            max(0, 12 - max_fiber // 6),
        )
        carriers.append(Carrier(name, axes, scores))
    carriers.sort(key=lambda carrier: (-carrier.total, carrier.name))
    score_hist = dict(sorted(Counter(carrier.total for carrier in carriers).items()))
    return score_hist, carriers


def fmt_counter(counter: Counter, limit: int = 10) -> str:
    items = counter.most_common(limit)
    suffix = "" if len(counter) <= limit else " ..."
    return "{" + ", ".join(f"{key!r}: {value}" for key, value in items) + suffix + "}"


def print_pure_table(
    records: Sequence[RowRecord],
    axis_names: Sequence[str],
    label: str,
    target: TargetFn,
) -> None:
    minimal = best_pure_subsets(records, axis_names, target)
    compressed = best_compressed_pure_subsets(records, axis_names, target)
    print(f"target={label}")
    print(f"  distribution={fmt_counter(target_distribution(records, target), 12)}")
    if minimal:
        print(
            "  fewest_axis_pure="
            + "; ".join(
                f"axes={axes} fibers={fibers} max_fiber={max_fiber}"
                for axis_count, fibers, max_fiber, axes in minimal[:5]
            )
        )
        print(
            "  most_compressed_pure="
            + "; ".join(
                f"axes={axes} fibers={fibers} max_fiber={max_fiber}"
                for fibers, axis_count, max_fiber, axes in compressed[:5]
            )
        )
    else:
        print("  pure_subset=NONE among declared axes")


def print_duplicate_fibers(records: Sequence[RowRecord], axes: tuple[str, ...], limit: int = 5) -> None:
    fibers = partition(records, axes)
    duplicates = [fiber for fiber in fibers.values() if len(fiber) > 1]
    duplicates.sort(key=lambda fiber: (-len(fiber), [record.name for record in fiber]))
    print(f"axes={axes} fibers={len(fibers)} duplicate_fibers={len(duplicates)}")
    for fiber in duplicates[:limit]:
        print(
            "  size={size} flags={flags} min_struct={min_struct} rows={rows}".format(
                size=len(fiber),
                flags=fiber[0].route_flags,
                min_struct=fiber[0].dead_min_struct_target,
                rows=[record.name for record in fiber[:12]],
            )
        )


def main() -> None:
    records = build_records()
    axis_names = ("K", "N", "T", "S", "F", "C", "M", "A")
    axis_labels = {
        "K": "row endpoint-kind set",
        "N": "row numeric mod-14 residue set",
        "T": "row typed mod-14 residue set",
        "S": "row structural sidecar set",
        "F": "row full colored-gate set",
        "C": "row low-rank gate count profile",
        "M": "minimum E/branch structural gate",
        "A": "canonical AP84 palette presence bit",
    }
    targets: tuple[tuple[str, TargetFn], ...] = (
        ("dead_bit", lambda r: r.dead),
        ("theorem_failure_bit", lambda r: r.theorem_failure),
        ("ap84_packet_bit", lambda r: r.has_ap84_packet),
        ("same_branch_gate_bit", lambda r: r.has_same_branch),
        ("cross_branch_gate_bit", lambda r: r.has_cross_branch),
        ("only_e_branch_gate_bit", lambda r: r.only_e_branch),
        ("route_flags", lambda r: r.route_flags),
        ("dead_min_struct_family", lambda r: r.dead_min_struct_target),
        ("dead_min_full_family", lambda r: r.dead_min_full_target),
    )

    print("HYP-3474 COLORED GATE PARTITION-LATTICE AUDIT")
    print("status=EVIDENCE / finite quotient-legality scout; not an LRC14 proof")
    print("source=HYP-3471 colored gate-reservoir over HYP-3453 survivor-gate/component-cover join")
    print()
    print("## Assumption Challenge")
    print("alternate vertices considered: runners, arcs, raw colors, endpoint kinds,")
    print("numeric residues, typed residues, structural sidecars, row-level color")
    print("sets, minimum E/branch gates, cover components, and proof obligations.")
    print("chosen carrier vertices: row-level quotient fibers in the partition")
    print("lattice generated by eight colored-gate axes.")
    print("preserved predicate: target-purity of proof-facing row predicates inside")
    print("each quotient fiber.")
    print("destroyed information: exact interval geometry, owner-current routing,")
    print("component adjacency, and row identity unless a color/sidecar axis")
    print("explicitly reconstructs it.")
    print()
    print("## Axis Legend")
    for axis in axis_names:
        print(f"{axis}={axis_labels[axis]}")
    print()
    print("## Bank Summary")
    print(f"rows={len(records)}")
    print(f"dead_rows={sum(record.dead for record in records)}")
    print(f"theorem_failures={sum(record.theorem_failure for record in records)}")
    print(f"ap84_packet_rows={sum(record.has_ap84_packet for record in records)}")
    print(f"same_branch_rows={sum(record.has_same_branch for record in records)}")
    print(f"cross_branch_rows={sum(record.has_cross_branch for record in records)}")
    print(f"only_e_branch_rows={sum(record.only_e_branch for record in records)}")
    print()
    print("## Singleton Quotient Mixing")
    for axis in axis_names:
        fibers = partition(records, (axis,))
        max_fiber = max(len(fiber) for fiber in fibers.values())
        route_mixed, route_mixed_rows, _, route_samples = mixed_report(
            records, (axis,), lambda r: r.route_flags
        )
        struct_mixed, struct_mixed_rows, _, _ = mixed_report(
            records, (axis,), lambda r: r.dead_min_struct_target
        )
        print(
            f"axis={axis} fibers={len(fibers)} max_fiber={max_fiber} "
            f"mixed_route_fibers={route_mixed} mixed_route_rows={route_mixed_rows} "
            f"mixed_dead_min_struct_fibers={struct_mixed} "
            f"mixed_dead_min_struct_rows={struct_mixed_rows}"
        )
        if route_samples:
            size, names, values = route_samples[0]
            print(f"  largest_route_mixed_sample size={size} rows={names} values={values}")
    print()
    print("## Pure-Quotient Search")
    for label, target in targets:
        print_pure_table(records, axis_names, label, target)
    print()
    print("## Coarsest Legal Route Quotients")
    route_compressed = best_compressed_pure_subsets(records, axis_names, lambda r: r.route_flags, 12)
    for fibers, axis_count, max_fiber, axes in route_compressed:
        print(f"route_flags_pure axes={axes} fibers={fibers} axis_count={axis_count} max_fiber={max_fiber}")
    print()
    print("## Largest Duplicate Fibers For Best Quotients")
    print_duplicate_fibers(records, ("C",))
    print_duplicate_fibers(records, ("N",))
    print_duplicate_fibers(records, ("T",))
    print_duplicate_fibers(records, ("F",))
    print_duplicate_fibers(records, ("C", "M"))
    print()
    print("## Dead Minimum Gate Family Compression")
    struct_compressed = best_compressed_pure_subsets(
        records, axis_names, lambda r: r.dead_min_struct_target, 12
    )
    for fibers, axis_count, max_fiber, axes in struct_compressed:
        print(
            "dead_min_struct_pure "
            f"axes={axes} fibers={fibers} axis_count={axis_count} max_fiber={max_fiber}"
        )
    full_compressed = best_compressed_pure_subsets(
        records, axis_names, lambda r: r.dead_min_full_target, 8
    )
    for fibers, axis_count, max_fiber, axes in full_compressed:
        print(
            "dead_min_full_pure "
            f"axes={axes} fibers={fibers} axis_count={axis_count} max_fiber={max_fiber}"
        )
    print()
    print("## Theorem-Failure Quotient Note")
    print("theorem_failure_bit is identically false on the HYP-3471 bank.")
    print("Therefore every quotient preserves the audited implication as a fact of")
    print("this bank, but only target-pure quotients preserve the routing data needed")
    print("to turn the bank fact into a labelled packet theorem.")
    print()
    print("## Tournament Analysis")
    score_hist, carriers = tournament(records)
    print("vertices=quotient carriers/sidecar selections, not runners or arcs")
    print("pairwise_observable=route-predicate purity + mixed-fiber penalty + compression")
    print("switch=higher aggregate proof-facing score; ties use carrier name")
    print("axes=pure_targets,mixed_fiber_penalty,mixed_row_penalty,compression_bonus,max_fiber_penalty")
    print(f"score_hist={score_hist}")
    print("directed_3cycles=0")
    print("hamiltonian_path=" + " -> ".join(carrier.name for carrier in carriers))
    for carrier in carriers:
        print(f"carrier={carrier.name} axes={carrier.axes} score={carrier.total} scores={carrier.scores}")
    print()
    print("## Proof Pull")
    print("The computation suggests a packet theorem target:")
    print("  dead-cover obstruction -> minimum E/branch structural gate -> route flags")
    print("where the quotient is allowed to forget only after the partition fiber is")
    print("pure for the relevant route target, or after a named sidecar reconstructs")
    print("the missing coordinate.")
    print("Endpoint-kind and numeric-residue quotients are useful scouts, but they")
    print("mix route flags heavily.  The minimum E/branch structural gate M is the")
    print("smallest single sidecar that sees much of the obstruction skeleton; the")
    print("full row-colored gate set F is the strongest finite fingerprint.")


if __name__ == "__main__":
    main()
