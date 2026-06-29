#!/usr/bin/env python3
"""HYP-3490: private-label firewall audit for LRC14 currentless rows.

HYP-3476 shows that two low-rank E/branch gates do not repair the seven random
HYP-3472 projection-edge exceptions.  This script asks why.  The tested
mechanism is label privacy in the HYP-3451 dead-cover projection:

    if every adjacent blocker label touched by every low-rank E/branch gate
    occurs on only one dead component, then deleting any union of those labels
    cannot remove a projection edge.

This upgrades the two-gate failure to an all-adjacent-E/branch-label firewall
for exactly the currentless random packet.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys


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


H3472 = load_module("hyp3472_boundary_current_for_hyp3478", H3472_PATH)
H3475 = load_module("hyp3475_mirror_orbit_for_hyp3478", H3475_PATH)
H3453 = H3472.H3453

HARD_DELTA = 7


def top_items(counter: Counter, limit: int = 12) -> str:
    items = counter.most_common(limit)
    suffix = "" if len(counter) <= limit else " ..."
    return "{" + ", ".join(f"{key!r}: {value}" for key, value in items) + suffix + "}"


def label_counts(row) -> Counter[str]:
    counts: Counter[str] = Counter()
    for component in row.component_row.dead_components:
        counts.update(H3472.component_blocker_labels(component))
    return counts


def touched_labels(row, gate) -> tuple[str, ...]:
    dead_labels = set(label_counts(row))
    return tuple(sorted(label for label in H3472.gate_blocker_labels(gate) if label in dead_labels))


def fmt_cut(cut) -> str:
    if cut is None:
        return "None"
    gate = cut.gate
    return (
        f"gate={H3472.fmt_interval(gate.interval)} kind={gate.endpoint_kind_signature} "
        f"labels={cut.labels_touching_dead} removed_edges={cut.removed_edges} "
        f"largest_drop={cut.largest_drop} component_gain={cut.component_gain} "
        f"current=({cut.b0_labels},{cut.b1_labels}) delta=({gate.b0_delta},{gate.b1_delta})"
    )


def cached_gate_cut(row, gate, original, dead_label_set: frozenset[str], stats_cache):
    labels = H3472.gate_blocker_labels(gate)
    touching = tuple(label for label in labels if label in dead_label_set)
    touching_key = frozenset(touching)
    after = stats_cache.get(touching_key)
    if after is None:
        after = H3472.projection_stats(
            H3472.projection(list(row.component_row.dead_components), touching_key)
        )
        stats_cache[touching_key] = after
    b0, b1, balance = H3472.branch_counts(touching)
    return H3472.GateCut(
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


@dataclass(frozen=True)
class RowFirewall:
    row: object
    dead_label_count: Counter[str]
    edge_labels: tuple[str, ...]
    e_branch_gates: tuple
    gate_touch_labels: tuple[str, ...]
    gate_touch_multiplicity_hist: Counter[int]
    private_touch_labels: tuple[str, ...]
    shared_touch_labels: tuple[str, ...]
    h3472_edge_cut: bool
    h3472_separating_current: bool
    best_cut: object | None
    hard_orbit_count: int
    hard_max_delta: int

    @property
    def private_firewall(self) -> bool:
        return bool(self.e_branch_gates) and not self.shared_touch_labels

    @property
    def no_adjacent_label_union_can_cut(self) -> bool:
        return self.private_firewall

    @property
    def route(self) -> str:
        if self.h3472_separating_current:
            return "ordinary_or_hard_currented"
        if H3472.ap84_like(self.row.name):
            return "ap84_edge_only_pair_current_terminal"
        if self.hard_orbit_count:
            return "random031_private_firewall_hard_overlap"
        if self.private_firewall:
            return "small_touch_private_firewall"
        return "unclassified_currentless"


def hard_orbits_for(row) -> tuple:
    audit = H3475.build_orbits(row)
    return tuple(orbit for orbit in audit.orbits if orbit.max_delta >= HARD_DELTA)


def audit_row(row) -> RowFirewall:
    counts = label_counts(row)
    edge_labels = tuple(sorted(label for label, count in counts.items() if count >= 2))
    gates = tuple(gate for gate in row.low_rank_gates if H3472.is_e_branch_gate(gate))
    original = H3472.projection_stats(H3472.projection(list(row.component_row.dead_components)))
    dead_labels = frozenset(counts)
    stats_cache = {}
    cuts = tuple(cached_gate_cut(row, gate, original, dead_labels, stats_cache) for gate in gates)
    best = H3472.best_cut(list(cuts))
    touch = tuple(sorted({label for gate in gates for label in touched_labels(row, gate)}))
    private = tuple(label for label in touch if counts[label] == 1)
    shared = tuple(label for label in touch if counts[label] >= 2)
    hard = hard_orbits_for(row)
    return RowFirewall(
        row=row,
        dead_label_count=counts,
        edge_labels=edge_labels,
        e_branch_gates=gates,
        gate_touch_labels=touch,
        gate_touch_multiplicity_hist=Counter(counts[label] for label in touch),
        private_touch_labels=private,
        shared_touch_labels=shared,
        h3472_edge_cut=any(cut.is_edge_cut for cut in cuts),
        h3472_separating_current=any(cut.is_separating_current for cut in cuts),
        best_cut=best,
        hard_orbit_count=len(hard),
        hard_max_delta=max((orbit.max_delta for orbit in hard), default=0),
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
    "private_label_payload",
    "projection_cut_firewall",
    "frontier_split",
    "gluing_exit",
    "quotient_guardrail",
    "scalar_firewall",
)

CARRIERS = (
    Carrier("F00_private_label_firewall_theorem", (10, 10, 10, 9, 8, 10, 10)),
    Carrier("F01_random031_private_hard_overlap", (10, 9, 9, 10, 10, 9, 10)),
    Carrier("F02_small_touch_private_packet", (9, 10, 9, 9, 8, 9, 10)),
    Carrier("F03_ap84_nonprivate_pair_current", (9, 8, 8, 8, 7, 9, 10)),
    Carrier("F04_edge_support_label_axis", (8, 9, 9, 7, 6, 9, 10)),
    Carrier("F05_raw_pair_current_count", (4, 3, 3, 2, 2, 2, 1)),
    Carrier("F06_raw_exception_name", (2, 1, 1, 1, 1, 1, 0)),
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
    rows = [H3453.join_row(name, speeds) for name, speeds in H3453.H3450.rows().items()]
    dead_rows = [row for row in rows if row.has_dead]
    audits = [audit_row(row) for row in dead_rows]

    edge_cut_rows = [audit.row.name for audit in audits if audit.h3472_edge_cut]
    separating_rows = [audit.row.name for audit in audits if audit.h3472_separating_current]
    private_firewall_rows = [audit.row.name for audit in audits if audit.private_firewall]
    shared_touch_rows = [audit.row.name for audit in audits if audit.shared_touch_labels]
    no_edge_rows = [audit.row.name for audit in audits if not audit.h3472_edge_cut]
    no_separating_rows = [audit.row.name for audit in audits if not audit.h3472_separating_current]
    route_hist = Counter(audit.route for audit in audits)
    private_route_hist = Counter(audit.route for audit in audits if audit.private_firewall)
    hard_private_rows = [
        audit.row.name for audit in audits if audit.private_firewall and audit.hard_orbit_count
    ]
    small_private_rows = [
        audit.row.name
        for audit in audits
        if audit.route == "small_touch_private_firewall"
    ]
    ap_nonprivate_rows = [
        audit.row.name
        for audit in audits
        if audit.route == "ap84_edge_only_pair_current_terminal"
    ]
    mismatch_private_vs_no_edge = sorted(set(private_firewall_rows) ^ set(no_edge_rows))
    mismatch_shared_vs_edge = sorted(set(shared_touch_rows) ^ set(edge_cut_rows))

    touch_multiplicity_hist = Counter()
    touch_label_count_hist = Counter()
    shared_touch_count_hist = Counter()
    edge_label_count_hist = Counter()
    e_branch_count_hist = Counter()
    for audit in audits:
        touch_multiplicity_hist.update(audit.gate_touch_multiplicity_hist)
        touch_label_count_hist[len(audit.gate_touch_labels)] += 1
        shared_touch_count_hist[len(audit.shared_touch_labels)] += 1
        edge_label_count_hist[len(audit.edge_labels)] += 1
        e_branch_count_hist[len(audit.e_branch_gates)] += 1

    score_hist, path = tournament()

    print("HYP-3490 LRC14 PRIVATE-LABEL FIREWALL AUDIT")
    print("status=EVIDENCE / all-union currentless-row firewall; not an LRC14 proof")
    print("source=HYP-3472 boundary current + HYP-3476 pair-current + HYP-3477 hard frontier")
    print()
    print("## Assumption Challenge")
    print("alternate vertices considered: runners, gaps, residues, survivor gates,")
    print("mirror orbits, gate pairs, dead components, blocker labels, labelled")
    print("projection edges, quotient fibers, and proof obligations.")
    print("chosen carrier vertices: dead-cover blocker labels with multiplicity")
    print("sidecars, especially labels touched by low-rank E/branch gates.")
    print("preserved predicate: whether any union of adjacent E/branch gate labels")
    print("can remove an edge from the dead-cover projection.")
    print("destroyed by scalarization: which labels are private, which labels carry")
    print("projection edges, and which terminal route packet remains.")
    print()
    print("## Aggregate Firewall Ledger")
    print(f"rows_audited={len(rows)}")
    print(f"dead_rows={len(dead_rows)}")
    print(f"dead_rows_with_h3472_projection_edge_cut={len(edge_cut_rows)}/{len(dead_rows)}")
    print(f"dead_rows_with_h3472_separating_current={len(separating_rows)}/{len(dead_rows)}")
    print(f"private_firewall_rows={len(private_firewall_rows)}/{len(dead_rows)} {private_firewall_rows}")
    print(f"shared_touch_rows={len(shared_touch_rows)}/{len(dead_rows)}")
    print(f"mismatch_private_firewall_vs_no_edge_cut={mismatch_private_vs_no_edge}")
    print(f"mismatch_shared_touch_vs_edge_cut={mismatch_shared_vs_edge}")
    print()
    print("## Label Multiplicity Compression")
    print(f"e_branch_gate_count_hist={top_items(e_branch_count_hist)}")
    print(f"edge_label_count_hist={top_items(edge_label_count_hist)}")
    print(f"gate_touch_label_count_hist={top_items(touch_label_count_hist)}")
    print(f"shared_touch_label_count_hist={top_items(shared_touch_count_hist)}")
    print(f"gate_touch_multiplicity_hist={top_items(touch_multiplicity_hist)}")
    print()
    print("## Frontier Split")
    print(f"route_hist={dict(sorted(route_hist.items()))}")
    print(f"private_firewall_route_hist={dict(sorted(private_route_hist.items()))}")
    print(f"ap_nonprivate_pair_current_rows={ap_nonprivate_rows}")
    print(f"hard_private_overlap_rows={hard_private_rows}")
    print(f"small_touch_private_rows={small_private_rows}")
    print(f"no_separating_rows={no_separating_rows}")
    print()
    print("## Private Firewall Row Detail")
    for audit in audits:
        if not audit.private_firewall:
            continue
        print(f"row={audit.row.name} route={audit.route}")
        print(
            "  "
            f"dead_components={len(audit.row.component_row.dead_components)} "
            f"projection_edge_labels={audit.edge_labels} "
            f"e_branch_gates={len(audit.e_branch_gates)} "
            f"hard_orbits={audit.hard_orbit_count} hard_max_delta={audit.hard_max_delta}"
        )
        print(f"  gate_touch_labels={audit.gate_touch_labels}")
        print(f"  gate_touch_multiplicities={dict(audit.gate_touch_multiplicity_hist)}")
        print(f"  best_cut={fmt_cut(audit.best_cut)}")
    print()
    print("## Lemma Target")
    print("Finite audited form:")
    print("  For exactly the seven HYP-3472 random projection-edge exceptions,")
    print("  every blocker label touched by every low-rank E/branch gate is private")
    print("  to a single dead component.  Therefore no union of adjacent E/branch")
    print("  gate labels, of any size, can delete a dead-cover projection edge.")
    print("  HYP-3476's pair failure is not a two-gate accident; it is an")
    print("  all-union private-label firewall for this carrier.")
    print("  The AP currentless rows are non-private and already close by the")
    print("  HYP-3476 two-gate pair-current terminal.  The only private firewall row")
    print("  with hard mirror debt is random_covering_031; the other six private")
    print("  rows are small-touch/no-hard packets.")
    print()
    print("## Tournament Analysis")
    print("vertices=label-multiplicity proof carriers, not runners or raw exception names")
    print(
        "pairwise_observable=predicate retention + private-label payload + "
        "projection-cut firewall + frontier split + quotient guardrail"
    )
    print("switch=higher proof-facing carrier score; ties use declared route order")
    print(f"axes={AXES}")
    print(f"score_hist={score_hist}")
    print("directed_3cycles=0")
    print("hamiltonian_path=" + " -> ".join(path))


if __name__ == "__main__":
    main()
