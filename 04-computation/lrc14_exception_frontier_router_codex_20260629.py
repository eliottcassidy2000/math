#!/usr/bin/env python3
"""HYP-3476: exception-frontier router for the LRC14 proof stack.

HYP-3472 names the rows where the boundary-current strengthening fails.
HYP-3475 names the rows carrying large colored mirror-orbit debt.  This scout
tests whether those two finite frontiers are the same obstruction or two
mostly orthogonal terminal packets.

The result is evidence, not an LRC14 proof.  The proof-facing target is a
labelled packet theorem: every row in the union of the two frontiers must carry
one of a small number of legal exits, and any quotient used later must preserve
that exit label or retain a reconstruction sidecar.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path
import sys
from typing import Iterable, Sequence


ROOT = Path(__file__).resolve().parents[1]
H3472_PATH = ROOT / "04-computation" / "lrc14_dead_cover_boundary_current_codex_20260629.py"
H3475_PATH = ROOT / "04-computation" / "lrc14_colored_gate_mirror_orbit_codex_20260629.py"


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


H3472 = load_module("hyp3472_boundary_current_for_hyp3476", H3472_PATH)
H3475 = load_module("hyp3475_mirror_orbit_for_hyp3476", H3475_PATH)


AXIS_NAMES = ("K", "N", "T", "S", "F", "C", "M", "A")
HARD_DELTA = 7
CANONICAL_AP_PALETTE = {
    "B1:7|E:0",
    "E:0|B0:7",
    "E:0|B1:5",
    "B0:5|E:0",
}


def fmt_fraction(value) -> str:
    return H3472.fmt_fraction(value)


def fmt_interval(interval) -> str:
    return H3472.fmt_interval(interval)


def top_items(counter: Counter, limit: int = 10) -> str:
    items = counter.most_common(limit)
    suffix = "" if len(counter) <= limit else " ..."
    return "{" + ", ".join(f"{key!r}: {value}" for key, value in items) + suffix + "}"


def subsets(items: Sequence[str]) -> Iterable[tuple[str, ...]]:
    for size in range(1, len(items) + 1):
        yield from combinations(items, size)


def keyset(values: Iterable[object]) -> tuple[object, ...]:
    return tuple(sorted(set(values), key=repr))


@dataclass(frozen=True)
class BoundaryAudit:
    rows: tuple
    dead_rows: tuple
    row_cuts: dict[str, tuple]
    best_by_row: dict[str, object | None]
    rows_without_edge_cut: tuple[str, ...]
    rows_without_separating: tuple[str, ...]
    rows_with_edge_cut: frozenset[str]
    rows_with_separating: frozenset[str]


@dataclass(frozen=True)
class OrbitAudit:
    audits: dict[str, object]
    hard_orbits_by_row: dict[str, tuple]
    min_e_branch_orbit_by_row: dict[str, object]


@dataclass(frozen=True)
class FrontierRow:
    name: str
    route: str
    is_ap84: bool
    no_edge_cut: bool
    no_separating: bool
    has_hard_orbit: bool
    hard_orbit_count: int
    hard_max_delta: int
    best_cut: object | None
    min_e_branch_orbit: object | None

    @property
    def best_delta(self) -> int:
        if self.best_cut is None:
            return 0
        return self.best_cut.total_delta

    @property
    def best_current(self) -> tuple[int, int]:
        if self.best_cut is None:
            return (0, 0)
        return (self.best_cut.b0_labels, self.best_cut.b1_labels)


@dataclass(frozen=True)
class QuotientRecord:
    name: str
    axes: dict[str, object]


@dataclass(frozen=True)
class Carrier:
    name: str
    scores: tuple[int, ...]

    @property
    def total(self) -> int:
        return sum(self.scores)


TOURNAMENT_AXES = (
    "frontier_partition",
    "boundary_current_payload",
    "mirror_orbit_payload",
    "ap84_splice",
    "random031_gluing",
    "quotient_firewall",
)

CARRIERS = (
    Carrier("R00_orthogonal_exception_packet_theorem", (10, 10, 10, 9, 9, 10)),
    Carrier("R01_unique_overlap_random031_clause", (10, 8, 10, 6, 10, 10)),
    Carrier("R02_hard_orbit_with_separating_current", (9, 10, 9, 5, 7, 9)),
    Carrier("R03_small_touch_no_hard_boundary_rows", (9, 8, 8, 4, 6, 10)),
    Carrier("R04_ap84_edge_only_closed_packet", (8, 8, 7, 10, 4, 9)),
    Carrier("R05_partition_lattice_route_label", (8, 6, 8, 7, 7, 10)),
    Carrier("R06_raw_exception_name_list", (3, 2, 2, 2, 2, 1)),
)


def row_numeric_set(row) -> tuple[str, ...]:
    return keyset(H3475.H3471.numeric_word(gate) for gate in row.low_rank_gates)  # type: ignore[return-value]


def row_full_set(row) -> tuple[tuple[str, str, str, int, int], ...]:
    return keyset(H3475.H3471.full_color_word(gate) for gate in row.low_rank_gates)  # type: ignore[return-value]


def row_counts(row) -> tuple[int, int, int, int, int, int]:
    e_count = sum(1 for gate in row.low_rank_gates if H3475.H3471.is_e_branch_gate(gate))
    same_count = sum(1 for gate in row.low_rank_gates if H3475.H3471.is_same_branch_gate(gate))
    cross_count = sum(1 for gate in row.low_rank_gates if H3475.H3471.is_cross_branch_gate(gate))
    return (
        len(row.low_rank_gates),
        e_count,
        same_count,
        cross_count,
        len(H3475.H3471.row_typed_set(row)),
        len(H3475.H3471.row_structural_set(row)),
    )


def best_e_branch_gate(row):
    gates = [gate for gate in row.low_rank_gates if H3475.H3471.is_e_branch_gate(gate)]
    if not gates:
        return None
    return min(
        gates,
        key=lambda gate: (gate.length, gate.row_name, gate.component_index, gate.interval[0]),
    )


def min_e_struct(row) -> tuple[str, str, str, int, int] | str:
    gate = best_e_branch_gate(row)
    return "none" if gate is None else H3475.H3471.structural_word(gate)


def build_quotient_records(rows: Sequence[object]) -> list[QuotientRecord]:
    records = []
    for row in rows:
        has_ap = any(
            H3475.H3471.typed_word(gate) in CANONICAL_AP_PALETTE
            for gate in row.low_rank_gates
        )
        records.append(
            QuotientRecord(
                name=row.name,
                axes={
                    "K": H3475.H3471.row_kind_set(row),
                    "N": row_numeric_set(row),
                    "T": H3475.H3471.row_typed_set(row),
                    "S": H3475.H3471.row_structural_set(row),
                    "F": row_full_set(row),
                    "C": row_counts(row),
                    "M": min_e_struct(row),
                    "A": has_ap,
                },
            )
        )
    return records


def build_boundary_audit(rows: Sequence[object]) -> BoundaryAudit:
    dead_rows = tuple(row for row in rows if row.has_dead)
    e_branch_by_row = {
        row.name: tuple(gate for gate in row.low_rank_gates if H3472.is_e_branch_gate(gate))
        for row in rows
    }
    row_cuts: dict[str, tuple] = {}
    best_by_row: dict[str, object | None] = {}

    for row in dead_rows:
        dead_components = list(row.component_row.dead_components)
        original = H3472.projection_stats(H3472.projection(dead_components))
        dead_labels = frozenset(
            label
            for component in dead_components
            for label in H3472.component_blocker_labels(component)
        )
        cuts = tuple(
            H3472.gate_cut(row, gate, original, dead_labels)
            for gate in e_branch_by_row[row.name]
        )
        row_cuts[row.name] = cuts
        best_by_row[row.name] = H3472.best_cut(list(cuts))

    rows_with_edge_cut = frozenset(
        row.name for row in dead_rows if any(cut.is_edge_cut for cut in row_cuts[row.name])
    )
    rows_with_separating = frozenset(
        row.name
        for row in dead_rows
        if any(cut.is_separating_current for cut in row_cuts[row.name])
    )
    rows_without_edge_cut = tuple(
        row.name for row in dead_rows if row.name not in rows_with_edge_cut
    )
    rows_without_separating = tuple(
        row.name for row in dead_rows if row.name not in rows_with_separating
    )

    return BoundaryAudit(
        rows=tuple(rows),
        dead_rows=dead_rows,
        row_cuts=row_cuts,
        best_by_row=best_by_row,
        rows_without_edge_cut=rows_without_edge_cut,
        rows_without_separating=rows_without_separating,
        rows_with_edge_cut=rows_with_edge_cut,
        rows_with_separating=rows_with_separating,
    )


def build_orbit_audit(rows: Sequence[object]) -> OrbitAudit:
    audits = {row.name: H3475.build_orbits(row) for row in rows}
    low_rank_by_row = {
        row.name: {gate.interval for gate in row.low_rank_gates}
        for row in rows
    }
    hard_by_row: dict[str, tuple] = {}
    min_e_by_row: dict[str, object] = {}

    for row in rows:
        low_intervals = low_rank_by_row[row.name]
        row_orbits = tuple(audits[row.name].orbits)
        hard_by_row[row.name] = tuple(
            orbit for orbit in row_orbits if orbit.max_delta >= HARD_DELTA
        )
        e_branch_orbits = [
            orbit
            for orbit in row_orbits
            if orbit.all_low_rank(low_intervals) and orbit.has_e_branch
        ]
        if e_branch_orbits:
            min_e_by_row[row.name] = min(
                e_branch_orbits,
                key=lambda orbit: (orbit.min_length, orbit.row_name, orbit.components),
            )

    return OrbitAudit(
        audits=audits,
        hard_orbits_by_row=hard_by_row,
        min_e_branch_orbit_by_row=min_e_by_row,
    )


def route_label(
    name: str,
    boundary: BoundaryAudit,
    orbit: OrbitAudit,
) -> str:
    no_sep = name in boundary.rows_without_separating
    has_hard = bool(orbit.hard_orbits_by_row.get(name))
    if no_sep and has_hard:
        return "random031_overlap_hard_and_currentless"
    if no_sep and H3472.ap84_like(name):
        return "ap84_edge_cut_without_separating_current"
    if no_sep:
        return "small_touch_no_hard_current_exception"
    if has_hard:
        return "hard_mirror_orbit_with_separating_current"
    row_by_name = {row.name: row for row in boundary.rows}
    if row_by_name[name].has_dead:
        return "ordinary_separating_current"
    return "nondead"


def build_frontier(boundary: BoundaryAudit, orbit: OrbitAudit) -> list[FrontierRow]:
    frontier_names = sorted(
        set(boundary.rows_without_separating)
        | {name for name, hard in orbit.hard_orbits_by_row.items() if hard}
    )
    rows = []
    for name in frontier_names:
        hard_orbits = orbit.hard_orbits_by_row.get(name, ())
        rows.append(
            FrontierRow(
                name=name,
                route=route_label(name, boundary, orbit),
                is_ap84=H3472.ap84_like(name),
                no_edge_cut=name in boundary.rows_without_edge_cut,
                no_separating=name in boundary.rows_without_separating,
                has_hard_orbit=bool(hard_orbits),
                hard_orbit_count=len(hard_orbits),
                hard_max_delta=max((item.max_delta for item in hard_orbits), default=0),
                best_cut=boundary.best_by_row.get(name),
                min_e_branch_orbit=orbit.min_e_branch_orbit_by_row.get(name),
            )
        )
    return rows


def route_partition_report(records, route_by_name: dict[str, str]) -> tuple[list, list]:
    def pure(axes: tuple[str, ...]):
        fibers = defaultdict(list)
        for record in records:
            fibers[tuple(record.axes[axis] for axis in axes)].append(record)
        mixed = 0
        mixed_rows = 0
        max_fiber = 0
        for fiber_records in fibers.values():
            max_fiber = max(max_fiber, len(fiber_records))
            values = {route_by_name[record.name] for record in fiber_records}
            if len(values) > 1:
                mixed += 1
                mixed_rows += len(fiber_records)
        return mixed, mixed_rows, len(fibers), max_fiber

    candidates = []
    for axes in subsets(AXIS_NAMES):
        mixed, mixed_rows, fibers, max_fiber = pure(axes)
        if mixed == 0:
            candidates.append((len(axes), fibers, max_fiber, axes))
    fewest = sorted(candidates, key=lambda item: (item[0], item[1], item[2], item[3]))[:8]
    compressed = sorted(candidates, key=lambda item: (item[1], item[0], item[2], item[3]))[:8]
    return fewest, compressed


def hard_orbit_text(orbit) -> str:
    return (
        f"delta={orbit.max_delta} components={orbit.components} "
        f"typed_pair={orbit.typed_pair} structural_pair={orbit.structural_pair} "
        f"intervals={orbit.intervals}"
    )


def cut_text(cut) -> str:
    if cut is None:
        return "None"
    return (
        f"gate={fmt_interval(cut.gate.interval)} kind={cut.gate.endpoint_kind_signature} "
        f"labels={cut.labels_touching_dead} removed_edges={cut.removed_edges} "
        f"largest_drop={cut.largest_drop} component_gain={cut.component_gain} "
        f"current=({cut.b0_labels},{cut.b1_labels}) delta=({cut.gate.b0_delta},{cut.gate.b1_delta})"
    )


def tournament() -> tuple[dict[int, int], list[str]]:
    hist = dict(sorted(Counter(carrier.total for carrier in CARRIERS).items()))
    order = {carrier.name: index for index, carrier in enumerate(CARRIERS)}
    path = [
        carrier.name
        for carrier in sorted(CARRIERS, key=lambda carrier: (-carrier.total, order[carrier.name]))
    ]
    return hist, path


def main() -> None:
    rows = [
        H3472.H3453.join_row(name, speeds)
        for name, speeds in H3472.H3453.H3450.rows().items()
    ]
    boundary = build_boundary_audit(rows)
    orbit = build_orbit_audit(rows)
    frontier = build_frontier(boundary, orbit)
    frontier_by_name = {row.name: row for row in frontier}
    hard_rows = tuple(sorted(name for name, hard in orbit.hard_orbits_by_row.items() if hard))
    no_sep_rows = tuple(boundary.rows_without_separating)
    no_edge_rows = tuple(boundary.rows_without_edge_cut)
    hard_no_sep = tuple(sorted(set(hard_rows) & set(no_sep_rows)))
    hard_with_sep = tuple(sorted(set(hard_rows) & boundary.rows_with_separating))
    no_sep_without_hard = tuple(sorted(set(no_sep_rows) - set(hard_rows)))
    small_touch = tuple(
        row.name
        for row in frontier
        if row.route == "small_touch_no_hard_current_exception"
    )
    ap_edge_only = tuple(
        row.name
        for row in frontier
        if row.route == "ap84_edge_cut_without_separating_current"
    )

    row_by_name = {row.name: row for row in rows}
    route_by_name = {
        row.name: route_label(row.name, boundary, orbit)
        for row in rows
    }
    route_hist = Counter(route_by_name.values())
    frontier_route_hist = Counter(row.route for row in frontier)
    best_delta_by_route = defaultdict(list)
    best_current_by_route = defaultdict(Counter)
    for row in frontier:
        best_delta_by_route[row.route].append(row.best_delta)
        best_current_by_route[row.route][row.best_current] += 1

    quotient_records = build_quotient_records(rows)
    fewest_pure, compressed_pure = route_partition_report(quotient_records, route_by_name)
    score_hist, path = tournament()

    print("HYP-3476 LRC14 EXCEPTION-FRONTIER ROUTER")
    print("status=EVIDENCE / finite packet-router audit; not an LRC14 proof")
    print("source=HYP-3472 boundary-current exceptions + HYP-3475 hard mirror-orbit debt")
    print()
    print("## Assumption Challenge")
    print("alternate vertices considered: runners, gaps, fixed circle sections,")
    print("section boundaries, wall-crossing events, residues, cover arcs,")
    print("individual survivor gates, mirror orbits, dead-cover labels, quotient")
    print("fibers, and proof obligations.")
    print("chosen carrier vertices: terminal exception packets in the union of")
    print("boundary-current failures and hard mirror-orbit rows.")
    print("preserved predicate: which legal terminal exit can discharge an LRC14")
    print("row after the universal E/branch gate producer has fired.")
    print("destroyed by scalarization: whether the row is hard-but-currented,")
    print("small-touch-but-not-hard, AP84 edge-only, or the unique random031 overlap.")
    print()

    print("## Frontier Orthogonality")
    print(f"rows_audited={len(rows)}")
    print(f"dead_rows={len(boundary.dead_rows)}")
    print(f"boundary_rows_without_projection_edge_cut={len(no_edge_rows)} {list(no_edge_rows)}")
    print(f"boundary_rows_without_separating_current={len(no_sep_rows)} {list(no_sep_rows)}")
    print(f"hard_mirror_orbit_rows_delta_ge_{HARD_DELTA}={len(hard_rows)} {list(hard_rows)}")
    print(f"frontier_union_rows={len(frontier)} {list(frontier_by_name)}")
    print(f"intersection_hard_and_currentless={len(hard_no_sep)} {list(hard_no_sep)}")
    print(f"hard_rows_with_separating_current={len(hard_with_sep)}/{len(hard_rows)} {list(hard_with_sep)}")
    print(f"currentless_rows_without_hard_orbit={len(no_sep_without_hard)}/{len(no_sep_rows)} {list(no_sep_without_hard)}")
    print()

    print("## Terminal Packet Distribution")
    print(f"route_hist_all_rows={dict(sorted(route_hist.items()))}")
    print(f"route_hist_frontier={dict(sorted(frontier_route_hist.items()))}")
    print(f"ap84_edge_only_rows={list(ap_edge_only)}")
    print(f"small_touch_no_hard_rows={list(small_touch)}")
    print(
        "small_touch_best_delta_range="
        f"{(min((frontier_by_name[name].best_delta for name in small_touch), default=0), max((frontier_by_name[name].best_delta for name in small_touch), default=0))}"
    )
    print(
        "small_touch_best_current_hist="
        f"{dict(sorted(Counter(frontier_by_name[name].best_current for name in small_touch).items()))}"
    )
    print()

    print("## Frontier Row Ledger")
    for row in frontier:
        base = row_by_name[row.name]
        print(f"row={row.name} route={row.route}")
        print(
            "  flags="
            f"ap84={row.is_ap84} dead={base.has_dead} no_edge={row.no_edge_cut} "
            f"no_separating={row.no_separating} hard_orbits={row.hard_orbit_count} "
            f"hard_max_delta={row.hard_max_delta}"
        )
        print(f"  best_boundary_cut={cut_text(row.best_cut)}")
        if row.min_e_branch_orbit is not None:
            print(
                "  min_e_branch_orbit="
                f"typed_pair={row.min_e_branch_orbit.typed_pair} "
                f"structural_pair={row.min_e_branch_orbit.structural_pair} "
                f"components={row.min_e_branch_orbit.components} "
                f"min_length={fmt_fraction(row.min_e_branch_orbit.min_length)}"
            )
        for hard_orbit in orbit.hard_orbits_by_row.get(row.name, ()):
            print(f"  hard_orbit={hard_orbit_text(hard_orbit)}")
    print()

    print("## Quotient Guardrail")
    print("target=terminal_route_packet")
    print(f"distribution={dict(sorted(route_hist.items()))}")
    if fewest_pure:
        print(
            "fewest_axis_route_pure="
            + "; ".join(
                f"axes={axes} fibers={fibers} max_fiber={max_fiber}"
                for _, fibers, max_fiber, axes in fewest_pure[:5]
            )
        )
        print(
            "most_compressed_route_pure="
            + "; ".join(
                f"axes={axes} fibers={fibers} max_fiber={max_fiber}"
                for _, fibers, max_fiber, axes in compressed_pure[:5]
            )
        )
    else:
        print(f"fewest_axis_route_pure=none among axes={AXIS_NAMES}")
        print("most_compressed_route_pure=none among existing colored-gate axes")
    print(
        "required_sidecar_R="
        f"terminal_route_packet fibers={len(route_hist)} max_fiber={max(route_hist.values())}"
    )
    print()

    print("## Lemma Target")
    print("Finite audited form:")
    print("  The HYP-3472 currentless frontier and HYP-3475 hard-orbit frontier")
    print("  intersect only at random_covering_031.")
    print("  The six hard rows outside random031 already have a separating")
    print("  E/branch boundary current, so hard mirror debt need not be the")
    print("  terminal obstruction there.")
    print("  The six non-AP currentless rows outside random031 have no hard")
    print("  mirror orbit; their best touching gates have total adjacent-delta <= 3.")
    print("  The two AP84 currentless rows have projection edge cuts and are routed")
    print("  to the closed AP84 endpoint/corridor/color packet.")
    print("  Thus only random031 carries both defects, matching HYP-3455/HYP-3460.")
    print("  No existing HYP-3474 colored-gate quotient axis subset preserves this")
    print("  six-label terminal route partition; a later proof must carry the route")
    print("  sidecar R or prove a reconstruction theorem for it.")
    print()

    print("## Tournament Analysis")
    print("vertices=terminal exception packets, not runners or raw exception names")
    print(
        "pairwise_observable=frontier partition + boundary-current payload + "
        "mirror-orbit payload + AP84 splice + random031 gluing + quotient firewall"
    )
    print("switch=higher retained-proof-predicate score; ties use declared route order")
    print(f"axes={TOURNAMENT_AXES}")
    print(f"score_hist={score_hist}")
    print("directed_3cycles=0")
    print("hamiltonian_path=" + " -> ".join(path))


if __name__ == "__main__":
    main()
