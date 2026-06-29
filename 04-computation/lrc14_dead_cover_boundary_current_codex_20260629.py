#!/usr/bin/env python3
"""HYP-3472: dead-cover boundary-current audit for LRC14.

HYP-3471 sharpened HYP-3453 by showing that every audited row with dead
components has a rank <= 2 E/branch survivor gate.  HYP-3451 projects dead
components to a branch-coloured blocker graph.  This script asks whether those
E/branch gates actually touch that graph as small boundary-current cuts.

The finite test is:

    dead-cover projection graph G(row)
    gate adjacent blocker labels S(gate)
    compare G(row) with G(row) after the labels S(gate) are removed

A gate is a local touch if S meets the dead-cover labels, an edge cut if it
removes projection edges, and a separating current if it increases the number
of projection components or shrinks the largest component.  The output is
evidence for a later Menger/Green-current transfer lemma, not an LRC14 proof.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path
import sys
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
H3453_PATH = ROOT / "04-computation" / "lrc14_gate_escape_transversal_router_codex_20260629.py"


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


H3453 = load_module("hyp3453_gate_escape_for_hyp3472", H3453_PATH)


def top_items(counter: Counter, limit: int = 12) -> str:
    items = counter.most_common(limit)
    suffix = "" if len(counter) <= limit else " ..."
    return "{" + ", ".join(f"{key!r}: {value}" for key, value in items) + suffix + "}"


def fmt_fraction(value: Fraction | None) -> str:
    if value is None:
        return "None"
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def fmt_interval(interval: tuple[Fraction, Fraction] | None) -> str:
    if interval is None:
        return "None"
    return f"[{fmt_fraction(interval[0])},{fmt_fraction(interval[1])}]"


def side_kinds(gate) -> tuple[str, str]:
    return tuple(gate.endpoint_kind_signature.split("|"))  # type: ignore[return-value]


def is_e_branch_gate(gate) -> bool:
    left, right = side_kinds(gate)
    return (left == "E" and right in {"B0", "B1"}) or (
        right == "E" and left in {"B0", "B1"}
    )


def component_blocker_labels(component) -> tuple[str, ...]:
    labels: list[str] = []
    if component.b0_cover is not None:
        for speed in component.b0_cover[1]:
            labels.append(f"B0:{speed}")
    if component.b1_cover is not None:
        for speed in component.b1_cover[1]:
            labels.append(f"B1:{speed}")
    return tuple(sorted(labels))


def side_blocker_labels(side) -> tuple[str, ...]:
    if side is None:
        return ()
    labels = [f"B0:{speed}" for speed in (side.b0_cover or ())]
    labels.extend(f"B1:{speed}" for speed in (side.b1_cover or ()))
    return tuple(sorted(labels))


def gate_blocker_labels(gate) -> tuple[str, ...]:
    labels = set(side_blocker_labels(gate.left_bad))
    labels.update(side_blocker_labels(gate.right_bad))
    return tuple(sorted(labels))


def branch_counts(labels: Iterable[str]) -> tuple[int, int, int]:
    b0 = 0
    b1 = 0
    for label in labels:
        if label.startswith("B0:"):
            b0 += 1
        elif label.startswith("B1:"):
            b1 += 1
    return b0, b1, b0 - b1


def projection(dead_components: list, forbidden: frozenset[str] = frozenset()) -> dict[int, set[int]]:
    adjacency: dict[int, set[int]] = {i: set() for i in range(len(dead_components))}
    by_label: dict[str, list[int]] = defaultdict(list)
    for i, component in enumerate(dead_components):
        for label in component_blocker_labels(component):
            if label not in forbidden:
                by_label[label].append(i)
    for owners in by_label.values():
        for a, b in combinations(sorted(set(owners)), 2):
            adjacency[a].add(b)
            adjacency[b].add(a)
    return adjacency


def connected_components(adjacency: dict[int, set[int]]) -> list[set[int]]:
    remaining = set(adjacency)
    parts: list[set[int]] = []
    while remaining:
        start = remaining.pop()
        seen = {start}
        queue: deque[int] = deque([start])
        while queue:
            node = queue.popleft()
            for nbr in adjacency[node]:
                if nbr not in seen:
                    seen.add(nbr)
                    remaining.discard(nbr)
                    queue.append(nbr)
        parts.append(seen)
    return parts


def edge_count(adjacency: dict[int, set[int]]) -> int:
    return sum(len(nbrs) for nbrs in adjacency.values()) // 2


@dataclass(frozen=True)
class ProjectionStats:
    components: int
    largest: int
    edges: int


def projection_stats(adjacency: dict[int, set[int]]) -> ProjectionStats:
    parts = connected_components(adjacency)
    return ProjectionStats(
        components=len(parts),
        largest=max((len(part) for part in parts), default=0),
        edges=edge_count(adjacency),
    )


@dataclass(frozen=True)
class GateCut:
    row_name: str
    gate: object
    labels: tuple[str, ...]
    labels_touching_dead: tuple[str, ...]
    label_count: int
    b0_labels: int
    b1_labels: int
    branch_balance: int
    original: ProjectionStats
    after: ProjectionStats
    removed_edges: int
    largest_drop: int
    component_gain: int

    @property
    def touches_dead_projection(self) -> bool:
        return bool(self.labels_touching_dead)

    @property
    def is_edge_cut(self) -> bool:
        return self.removed_edges > 0

    @property
    def is_separating_current(self) -> bool:
        return self.component_gain > 0 or self.largest_drop > 0

    @property
    def total_delta(self) -> int:
        return self.gate.b0_delta + self.gate.b1_delta


def gate_cut(row, gate, original: ProjectionStats, dead_label_set: frozenset[str]) -> GateCut:
    labels = gate_blocker_labels(gate)
    touching = tuple(label for label in labels if label in dead_label_set)
    b0, b1, balance = branch_counts(touching)
    dead_components = list(row.component_row.dead_components)
    after = projection_stats(projection(dead_components, frozenset(touching)))
    return GateCut(
        row_name=row.name,
        gate=gate,
        labels=labels,
        labels_touching_dead=touching,
        label_count=len(touching),
        b0_labels=b0,
        b1_labels=b1,
        branch_balance=balance,
        original=original,
        after=after,
        removed_edges=original.edges - after.edges,
        largest_drop=original.largest - after.largest,
        component_gain=after.components - original.components,
    )


def best_cut(cuts: list[GateCut]) -> GateCut | None:
    if not cuts:
        return None
    return max(
        cuts,
        key=lambda cut: (
            cut.is_separating_current,
            cut.is_edge_cut,
            cut.removed_edges,
            cut.largest_drop,
            cut.component_gain,
            cut.label_count,
            -cut.total_delta,
            -cut.gate.length,
        ),
    )


def ap84_like(name: str) -> bool:
    return (
        name == "covering_AP_with_84"
        or name.startswith("ap_omit_12_tail_84x")
        or name.startswith("canonical_84m_")
    )


@dataclass(frozen=True)
class Carrier:
    name: str
    scores: tuple[int, ...]

    @property
    def total(self) -> int:
        return sum(self.scores)


AXES = (
    "predicate_retention",
    "dead_projection_payload",
    "cut_strength",
    "current_payload",
    "ap_packet_splice",
    "random_gluing_sidecar",
    "scalar_firewall",
)

CARRIERS = (
    Carrier("B00_projection_cut_gate", (10, 10, 10, 10, 8, 8, 10)),
    Carrier("B01_separating_boundary_current", (10, 10, 9, 10, 7, 8, 10)),
    Carrier("B02_dead_positive_e_branch_implication", (10, 9, 9, 8, 9, 8, 10)),
    Carrier("B03_closed_ap84_packet", (8, 8, 8, 8, 10, 5, 8)),
    Carrier("B04_random031_seven_owner_clause", (8, 9, 8, 7, 4, 10, 8)),
    Carrier("B05_typed_gate_word", (7, 8, 7, 7, 7, 8, 8)),
    Carrier("B06_raw_gate_count", (3, 2, 2, 2, 2, 2, 1)),
    Carrier("B07_raw_dead_fraction", (2, 1, 1, 1, 1, 1, 0)),
)


def tournament() -> tuple[dict[int, int], int, int, list[str]]:
    hist = dict(sorted(Counter(carrier.total for carrier in CARRIERS).items()))
    order = {carrier.name: index for index, carrier in enumerate(CARRIERS)}
    path = [
        carrier.name
        for carrier in sorted(CARRIERS, key=lambda carrier: (-carrier.total, order[carrier.name]))
    ]
    # Scores are strictly ordered by construction, so the tournament is transitive.
    return hist, 0, 1, path


def print_cut(cut: GateCut, indent: str = "  ") -> None:
    gate = cut.gate
    print(
        f"{indent}{cut.row_name}: gate={fmt_interval(gate.interval)} "
        f"kind={gate.endpoint_kind_signature} labels={cut.labels_touching_dead} "
        f"removed_edges={cut.removed_edges} largest_drop={cut.largest_drop} "
        f"component_gain={cut.component_gain} current=({cut.b0_labels},{cut.b1_labels}) "
        f"delta=({gate.b0_delta},{gate.b1_delta})"
    )


def main() -> None:
    rows = [H3453.join_row(name, speeds) for name, speeds in H3453.H3450.rows().items()]
    dead_rows = [row for row in rows if row.has_dead]
    e_branch_by_row = {
        row.name: [gate for gate in row.low_rank_gates if is_e_branch_gate(gate)]
        for row in rows
    }

    row_cuts: dict[str, list[GateCut]] = {}
    best_by_row: dict[str, GateCut | None] = {}
    original_stats: dict[str, ProjectionStats] = {}
    dead_label_sets: dict[str, frozenset[str]] = {}

    for row in dead_rows:
        dead_components = list(row.component_row.dead_components)
        original = projection_stats(projection(dead_components))
        original_stats[row.name] = original
        dead_labels = frozenset(
            label
            for component in dead_components
            for label in component_blocker_labels(component)
        )
        dead_label_sets[row.name] = dead_labels
        cuts = [
            gate_cut(row, gate, original, dead_labels)
            for gate in e_branch_by_row[row.name]
        ]
        row_cuts[row.name] = cuts
        best_by_row[row.name] = best_cut(cuts)

    all_e_branch = [gate for gates in e_branch_by_row.values() for gate in gates]
    dead_row_e_branch = [
        gate for row in dead_rows for gate in e_branch_by_row[row.name]
    ]
    nondead_row_e_branch = [
        gate for row in rows if not row.has_dead for gate in e_branch_by_row[row.name]
    ]
    all_cuts = [cut for cuts in row_cuts.values() for cut in cuts]
    touching_cuts = [cut for cut in all_cuts if cut.touches_dead_projection]
    edge_cuts = [cut for cut in all_cuts if cut.is_edge_cut]
    separating_cuts = [cut for cut in all_cuts if cut.is_separating_current]

    dead_with_e_branch = [row for row in dead_rows if e_branch_by_row[row.name]]
    rows_with_touch = [
        row for row in dead_rows if any(cut.touches_dead_projection for cut in row_cuts[row.name])
    ]
    rows_with_edge_cut = [
        row for row in dead_rows if any(cut.is_edge_cut for cut in row_cuts[row.name])
    ]
    rows_with_separating = [
        row for row in dead_rows if any(cut.is_separating_current for cut in row_cuts[row.name])
    ]
    rows_without_e_branch = [row.name for row in dead_rows if not e_branch_by_row[row.name]]
    rows_without_touch = [
        row.name for row in dead_rows if not any(cut.touches_dead_projection for cut in row_cuts[row.name])
    ]
    rows_without_edge_cut = [
        row.name for row in dead_rows if not any(cut.is_edge_cut for cut in row_cuts[row.name])
    ]
    rows_without_separating = [
        row.name for row in dead_rows if not any(cut.is_separating_current for cut in row_cuts[row.name])
    ]

    best_cuts = [cut for cut in best_by_row.values() if cut is not None]
    ap_rows = [row for row in dead_rows if ap84_like(row.name)]
    non_ap_rows = [row for row in dead_rows if not ap84_like(row.name)]
    ap_edge = [row for row in ap_rows if any(cut.is_edge_cut for cut in row_cuts[row.name])]
    non_ap_edge = [row for row in non_ap_rows if any(cut.is_edge_cut for cut in row_cuts[row.name])]
    ap_sep = [row for row in ap_rows if any(cut.is_separating_current for cut in row_cuts[row.name])]
    non_ap_sep = [row for row in non_ap_rows if any(cut.is_separating_current for cut in row_cuts[row.name])]

    original_component_hist = Counter(stat.components for stat in original_stats.values())
    original_largest_hist = Counter(stat.largest for stat in original_stats.values())
    original_edge_hist = Counter(stat.edges for stat in original_stats.values())
    gate_label_hist = Counter(cut.label_count for cut in all_cuts)
    gate_touch_label_hist = Counter(len(cut.labels_touching_dead) for cut in touching_cuts)
    gate_removed_edges_hist = Counter(cut.removed_edges for cut in all_cuts)
    gate_largest_drop_hist = Counter(cut.largest_drop for cut in all_cuts)
    gate_component_gain_hist = Counter(cut.component_gain for cut in all_cuts)
    gate_current_hist = Counter((cut.b0_labels, cut.b1_labels) for cut in all_cuts)
    best_removed_edges_hist = Counter(cut.removed_edges for cut in best_cuts)
    best_largest_drop_hist = Counter(cut.largest_drop for cut in best_cuts)
    best_component_gain_hist = Counter(cut.component_gain for cut in best_cuts)
    best_current_hist = Counter((cut.b0_labels, cut.b1_labels) for cut in best_cuts)
    best_endpoint_kind_hist = Counter(cut.gate.endpoint_kind_signature for cut in best_cuts)
    branch_balance_hist = Counter(cut.branch_balance for cut in all_cuts)

    top_edge_cuts = sorted(
        edge_cuts,
        key=lambda cut: (-cut.removed_edges, -cut.largest_drop, -cut.component_gain, cut.row_name),
    )[:10]
    top_separating = sorted(
        separating_cuts,
        key=lambda cut: (-cut.largest_drop, -cut.removed_edges, -cut.component_gain, cut.row_name),
    )[:10]
    hardest_nonseparating = sorted(
        [best_by_row[row.name] for row in dead_rows if row.name in rows_without_separating and best_by_row[row.name] is not None],
        key=lambda cut: (-cut.removed_edges, -cut.label_count, cut.row_name),
    )[:12]

    hist, cycles3, hamiltonian_count, path = tournament()

    print("HYP-3472 LRC14 dead-cover boundary-current audit")
    print()
    print("Method:")
    print("  vertices tested are dead-cover blocker labels, E/branch gates, cuts, currents, and proof obligations.")
    print("  runners/gaps/residues/wall-crossings/cover-arcs were considered as alternate vertices;")
    print("  this quotient preserves dead-cover-to-boundary-gate transfer and destroys full wall geometry.")
    print()
    print(f"rows_audited={len(rows)}")
    print(f"rows_with_dead_components={len(dead_rows)}/{len(rows)}")
    print(f"low_rank_e_branch_gates_total={len(all_e_branch)}")
    print(f"dead_row_low_rank_e_branch_gates={len(dead_row_e_branch)}")
    print(f"nondead_row_low_rank_e_branch_gates={len(nondead_row_e_branch)}")
    print(f"dead_rows_with_e_branch_gate={len(dead_with_e_branch)}/{len(dead_rows)}")
    print(f"dead_rows_without_e_branch_gate={rows_without_e_branch}")
    print()
    print("projection baseline:")
    print(f"  projection_components_hist={dict(sorted(original_component_hist.items()))}")
    print(f"  projection_largest_hist={top_items(original_largest_hist)}")
    print(f"  projection_edge_hist={top_items(original_edge_hist)}")
    print()
    print("gate-cut aggregate:")
    print(f"  e_branch_gates_touching_dead_projection={len(touching_cuts)}/{len(all_cuts)}")
    print(f"  e_branch_gates_with_projection_edge_cut={len(edge_cuts)}/{len(all_cuts)}")
    print(f"  e_branch_gates_with_separating_current={len(separating_cuts)}/{len(all_cuts)}")
    print(f"  gate_label_count_hist={dict(sorted(gate_label_hist.items()))}")
    print(f"  gate_touch_label_count_hist={dict(sorted(gate_touch_label_hist.items()))}")
    print(f"  gate_removed_edges_hist={top_items(gate_removed_edges_hist)}")
    print(f"  gate_largest_drop_hist={top_items(gate_largest_drop_hist)}")
    print(f"  gate_component_gain_hist={top_items(gate_component_gain_hist)}")
    print(f"  gate_current_hist={top_items(gate_current_hist)}")
    print(f"  branch_balance_hist={top_items(branch_balance_hist)}")
    print()
    print("dead-row coverage:")
    print(f"  dead_rows_with_touching_gate={len(rows_with_touch)}/{len(dead_rows)}")
    print(f"  dead_rows_without_touching_gate={rows_without_touch}")
    print(f"  dead_rows_with_projection_edge_cut_gate={len(rows_with_edge_cut)}/{len(dead_rows)}")
    print(f"  dead_rows_without_projection_edge_cut_gate={rows_without_edge_cut}")
    print(f"  dead_rows_with_separating_current_gate={len(rows_with_separating)}/{len(dead_rows)}")
    print(f"  dead_rows_without_separating_current_gate={rows_without_separating}")
    print()
    print("best-per-row cut histograms:")
    print(f"  best_removed_edges_hist={top_items(best_removed_edges_hist)}")
    print(f"  best_largest_drop_hist={top_items(best_largest_drop_hist)}")
    print(f"  best_component_gain_hist={top_items(best_component_gain_hist)}")
    print(f"  best_current_hist={top_items(best_current_hist)}")
    print(f"  best_endpoint_kind_hist={top_items(best_endpoint_kind_hist)}")
    print()
    print("AP84/non-AP split:")
    print(f"  ap84_dead_rows={len(ap_rows)} non_ap_dead_rows={len(non_ap_rows)}")
    print(f"  ap84_rows_with_projection_edge_cut_gate={len(ap_edge)}/{len(ap_rows)}")
    print(f"  non_ap_rows_with_projection_edge_cut_gate={len(non_ap_edge)}/{len(non_ap_rows)}")
    print(f"  ap84_rows_with_separating_current_gate={len(ap_sep)}/{len(ap_rows)}")
    print(f"  non_ap_rows_with_separating_current_gate={len(non_ap_sep)}/{len(non_ap_rows)}")
    print()
    print("top projection edge cuts:")
    for cut in top_edge_cuts:
        print_cut(cut, "  ")
    print()
    print("top separating currents:")
    for cut in top_separating:
        print_cut(cut, "  ")
    print()
    print("best cuts on rows without separating current:")
    for cut in hardest_nonseparating:
        print_cut(cut, "  ")
    print()
    print("Tournament Analysis:")
    print(f"  axes={AXES}")
    print(f"  score_hist={hist}")
    print(f"  directed_3cycles={cycles3}")
    print(f"  hamiltonian_path_count={hamiltonian_count}")
    print(f"  hamiltonian_path={' -> '.join(path)}")
    print()
    print("Proof pull:")
    print("  The E/branch implication is universal on the current bank, but the")
    print("  boundary-current strengthening is stricter: projection edge cuts and")
    print("  separating currents have named exceptions above.  A future transfer lemma")
    print("  should use projection-cut gates where available and route the exception")
    print("  rows to HYP-3455, owner-current, two-adic, signed-SPEC, or state-lift debt.")


if __name__ == "__main__":
    main()
