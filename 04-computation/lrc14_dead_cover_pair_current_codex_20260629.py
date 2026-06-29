#!/usr/bin/env python3
"""HYP-3476: dead-cover pair-current exception audit for LRC14.

HYP-3472 showed that every audited dead-component row has a low-rank E/branch
survivor gate touching the HYP-3451 dead-cover blocker projection, but a single
gate fails to give a projection-edge cut on seven rows and fails to give a
separating current on those seven plus two AP84 base rows.

This script asks whether those exceptions are artifacts of using only one gate.
For each dead row, and especially for the HYP-3472 exception rows, enumerate
unordered pairs of low-rank E/branch gates.  Delete the union of the pair's
adjacent blocker labels from the dead-cover projection and score the remaining
graph.

The result is finite-bank evidence for a future terminal packet, not an LRC14
proof.
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
H3472_PATH = ROOT / "04-computation" / "lrc14_dead_cover_boundary_current_codex_20260629.py"

EDGE_EXCEPTION_ROWS = (
    "random_covering_001",
    "random_covering_031",
    "random_covering_039",
    "random_covering_062",
    "random_covering_074",
    "random_covering_086",
    "random_covering_101",
)
SEPARATING_EXCEPTION_ROWS = EDGE_EXCEPTION_ROWS + (
    "covering_AP_with_84",
    "ap_omit_12_tail_84x01",
)
H3475_HARD_MIRROR_ROWS = {
    "random_covering_022",
    "random_covering_031",
    "random_covering_049",
    "random_covering_078",
    "random_covering_080",
    "random_covering_085",
    "random_covering_113",
}

MIRROR_KIND = {
    "B0|E": "E|B1",
    "E|B1": "B0|E",
    "B1|E": "E|B0",
    "E|B0": "B1|E",
}


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


H3472 = load_module("hyp3472_boundary_current_for_hyp3476", H3472_PATH)
H3453 = H3472.H3453


def top_items(counter: Counter, limit: int = 12) -> str:
    items = counter.most_common(limit)
    suffix = "" if len(counter) <= limit else " ..."
    return "{" + ", ".join(f"{key!r}: {value}" for key, value in items) + suffix + "}"


def fmt_interval(interval: tuple[Fraction, Fraction] | None) -> str:
    return H3472.fmt_interval(interval)


def mirror_interval(interval: tuple[Fraction, Fraction]) -> tuple[Fraction, Fraction]:
    return (1 - interval[1], 1 - interval[0])


def gate_id(gate) -> tuple[int, Fraction, Fraction, str]:
    return (
        gate.component_index,
        gate.interval[0],
        gate.interval[1],
        gate.endpoint_kind_signature,
    )


def gate_text(gate) -> str:
    return (
        f"c{gate.component_index}:{fmt_interval(gate.interval)}:"
        f"{gate.endpoint_kind_signature}:d({gate.b0_delta},{gate.b1_delta})"
    )


def is_exact_mirror_pair(left, right) -> bool:
    return (
        right.interval == mirror_interval(left.interval)
        and MIRROR_KIND.get(left.endpoint_kind_signature) == right.endpoint_kind_signature
    ) or (
        left.interval == mirror_interval(right.interval)
        and MIRROR_KIND.get(right.endpoint_kind_signature) == left.endpoint_kind_signature
    )


def pair_sort_key(pair: "PairCut") -> tuple:
    return (
        pair.row_name,
        not pair.exact_mirror_pair,
        pair.gate_a.component_index,
        pair.gate_b.component_index,
        pair.gate_a.interval[0],
        pair.gate_b.interval[0],
    )


@dataclass(frozen=True)
class PairCut:
    row_name: str
    gate_a: object
    gate_b: object
    labels: tuple[str, ...]
    labels_touching_dead: tuple[str, ...]
    label_count: int
    b0_labels: int
    b1_labels: int
    branch_balance: int
    original: object
    after: object
    removed_edges: int
    largest_drop: int
    component_gain: int
    exact_mirror_pair: bool

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
        return (
            self.gate_a.b0_delta
            + self.gate_a.b1_delta
            + self.gate_b.b0_delta
            + self.gate_b.b1_delta
        )

    @property
    def endpoint_pair(self) -> tuple[str, str]:
        return tuple(sorted((self.gate_a.endpoint_kind_signature, self.gate_b.endpoint_kind_signature)))  # type: ignore[return-value]


@dataclass(frozen=True)
class RowPairAudit:
    row: object
    original: object
    e_branch_gates: tuple
    single_cuts: tuple
    pair_cuts: tuple[PairCut, ...]
    dead_label_count: Counter
    edge_support_labels: tuple[str, ...]
    gate_touch_labels: tuple[str, ...]
    gate_edge_support_labels: tuple[str, ...]

    @property
    def has_single_edge_cut(self) -> bool:
        return any(cut.is_edge_cut for cut in self.single_cuts)

    @property
    def has_single_separating_current(self) -> bool:
        return any(cut.is_separating_current for cut in self.single_cuts)

    @property
    def has_pair_edge_cut(self) -> bool:
        return any(cut.is_edge_cut for cut in self.pair_cuts)

    @property
    def has_pair_separating_current(self) -> bool:
        return any(cut.is_separating_current for cut in self.pair_cuts)

    @property
    def has_mirror_pair(self) -> bool:
        return any(cut.exact_mirror_pair for cut in self.pair_cuts)

    @property
    def has_mirror_edge_cut(self) -> bool:
        return any(cut.exact_mirror_pair and cut.is_edge_cut for cut in self.pair_cuts)

    @property
    def has_mirror_separating_current(self) -> bool:
        return any(cut.exact_mirror_pair and cut.is_separating_current for cut in self.pair_cuts)


@dataclass(frozen=True)
class Carrier:
    name: str
    scores: tuple[int, ...]

    @property
    def total(self) -> int:
        return sum(self.scores)


AXES = (
    "predicate_retention",
    "pair_current_payload",
    "mirror_payload",
    "cut_strength",
    "formal_packet_fit",
    "quotient_guardrail",
    "scalar_firewall",
)

CARRIERS = (
    Carrier("P00_two_gate_pair_current_terminal", (10, 10, 8, 10, 9, 9, 10)),
    Carrier("P01_mirror_pair_current_terminal", (10, 9, 10, 9, 9, 9, 10)),
    Carrier("P02_single_gate_boundary_current", (9, 8, 6, 8, 9, 8, 10)),
    Carrier("P03_colored_gate_formal_packet", (9, 8, 7, 7, 10, 9, 10)),
    Carrier("P04_partition_lattice_guardrail", (8, 7, 7, 6, 8, 10, 10)),
    Carrier("P05_hard_mirror_orbit_discharge", (8, 8, 10, 7, 8, 8, 9)),
    Carrier("P06_typed_gate_pair_word", (7, 7, 7, 6, 7, 8, 8)),
    Carrier("P07_raw_pair_count", (3, 2, 1, 1, 1, 2, 1)),
    Carrier("P08_raw_label_union_size", (2, 1, 1, 1, 1, 1, 0)),
)


def pair_cut(row, gate_a, gate_b, original, dead_label_set: frozenset[str]) -> PairCut:
    labels = tuple(
        sorted(set(H3472.gate_blocker_labels(gate_a)) | set(H3472.gate_blocker_labels(gate_b)))
    )
    touching = tuple(label for label in labels if label in dead_label_set)
    b0, b1, balance = H3472.branch_counts(touching)
    dead_components = list(row.component_row.dead_components)
    after = H3472.projection_stats(H3472.projection(dead_components, frozenset(touching)))
    return PairCut(
        row_name=row.name,
        gate_a=gate_a,
        gate_b=gate_b,
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
        exact_mirror_pair=is_exact_mirror_pair(gate_a, gate_b),
    )


def best_single(cuts: tuple) -> object | None:
    return H3472.best_cut(list(cuts))


def strongest_pair(cuts: tuple[PairCut, ...]) -> PairCut | None:
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
            cut.exact_mirror_pair,
            -cut.label_count,
            -cut.total_delta,
            -cut.gate_a.length - cut.gate_b.length,
        ),
    )


def minimal_edge_pair(cuts: tuple[PairCut, ...]) -> PairCut | None:
    edge_cuts = [cut for cut in cuts if cut.is_edge_cut]
    if not edge_cuts:
        return None
    return min(
        edge_cuts,
        key=lambda cut: (
            cut.label_count,
            not cut.exact_mirror_pair,
            cut.total_delta,
            -cut.removed_edges,
            -cut.largest_drop,
            pair_sort_key(cut),
        ),
    )


def minimal_separating_pair(cuts: tuple[PairCut, ...]) -> PairCut | None:
    sep_cuts = [cut for cut in cuts if cut.is_separating_current]
    if not sep_cuts:
        return None
    return min(
        sep_cuts,
        key=lambda cut: (
            cut.label_count,
            not cut.exact_mirror_pair,
            cut.total_delta,
            -cut.largest_drop,
            -cut.removed_edges,
            pair_sort_key(cut),
        ),
    )


def audit_row(row) -> RowPairAudit:
    dead_components = list(row.component_row.dead_components)
    original = H3472.projection_stats(H3472.projection(dead_components))
    dead_label_count = Counter(
        label
        for component in dead_components
        for label in H3472.component_blocker_labels(component)
    )
    dead_labels = frozenset(dead_label_count)
    gates = tuple(gate for gate in row.low_rank_gates if H3472.is_e_branch_gate(gate))
    single = tuple(H3472.gate_cut(row, gate, original, dead_labels) for gate in gates)
    pairs = tuple(pair_cut(row, a, b, original, dead_labels) for a, b in combinations(gates, 2))
    edge_support = tuple(sorted(label for label, count in dead_label_count.items() if count >= 2))
    gate_touch = tuple(
        sorted({label for cut in single for label in cut.labels_touching_dead})
    )
    gate_edge_support = tuple(sorted(set(edge_support) & set(gate_touch)))
    return RowPairAudit(
        row=row,
        original=original,
        e_branch_gates=gates,
        single_cuts=single,
        pair_cuts=pairs,
        dead_label_count=dead_label_count,
        edge_support_labels=edge_support,
        gate_touch_labels=gate_touch,
        gate_edge_support_labels=gate_edge_support,
    )


def print_single(cut, indent: str = "  ") -> None:
    if cut is None:
        print(f"{indent}None")
        return
    gate = cut.gate
    print(
        f"{indent}single {gate_text(gate)} labels={cut.labels_touching_dead} "
        f"removed_edges={cut.removed_edges} largest_drop={cut.largest_drop} "
        f"component_gain={cut.component_gain} current=({cut.b0_labels},{cut.b1_labels})"
    )


def print_pair(cut: PairCut | None, indent: str = "  ") -> None:
    if cut is None:
        print(f"{indent}None")
        return
    print(
        f"{indent}pair mirror={cut.exact_mirror_pair} labels={cut.labels_touching_dead} "
        f"removed_edges={cut.removed_edges} largest_drop={cut.largest_drop} "
        f"component_gain={cut.component_gain} current=({cut.b0_labels},{cut.b1_labels})"
    )
    print(f"{indent}  a={gate_text(cut.gate_a)}")
    print(f"{indent}  b={gate_text(cut.gate_b)}")


def tournament() -> tuple[dict[int, int], list[str]]:
    hist = dict(sorted(Counter(carrier.total for carrier in CARRIERS).items()))
    order = {carrier.name: index for index, carrier in enumerate(CARRIERS)}
    path = [
        carrier.name
        for carrier in sorted(CARRIERS, key=lambda carrier: (-carrier.total, order[carrier.name]))
    ]
    return hist, path


def main() -> None:
    full_sweep = "--all-dead-rows" in sys.argv
    row_specs = H3453.H3450.rows()
    if full_sweep:
        rows = [H3453.join_row(name, speeds) for name, speeds in row_specs.items()]
        dead_rows = [row for row in rows if row.has_dead]
        rows_to_audit = dead_rows
        bank_dead_rows = len(dead_rows)
        scope = "all_dead_rows"
    else:
        target_names = set(SEPARATING_EXCEPTION_ROWS)
        rows = [
            H3453.join_row(name, row_specs[name])
            for name in SEPARATING_EXCEPTION_ROWS
            if name in row_specs
        ]
        rows_to_audit = [row for row in rows if row.has_dead]
        bank_dead_rows = 130
        scope = "hyp3472_exception_targets"

    audits = [audit_row(row) for row in rows_to_audit]
    audit_by_name = {audit.row.name: audit for audit in audits}

    all_pairs = [cut for audit in audits for cut in audit.pair_cuts]
    all_single = [cut for audit in audits for cut in audit.single_cuts]
    target_edge = [audit_by_name[name] for name in EDGE_EXCEPTION_ROWS]
    target_sep = [audit_by_name[name] for name in SEPARATING_EXCEPTION_ROWS]

    rows_with_single_edge = [audit for audit in audits if audit.has_single_edge_cut]
    rows_with_pair_edge = [audit for audit in audits if audit.has_pair_edge_cut]
    rows_with_single_sep = [audit for audit in audits if audit.has_single_separating_current]
    rows_with_pair_sep = [audit for audit in audits if audit.has_pair_separating_current]
    rows_with_mirror_edge = [audit for audit in audits if audit.has_mirror_edge_cut]
    rows_with_mirror_sep = [audit for audit in audits if audit.has_mirror_separating_current]

    edge_exceptions_closed_by_pair_edge = [
        audit.row.name for audit in target_edge if audit.has_pair_edge_cut
    ]
    edge_exceptions_closed_by_pair_sep = [
        audit.row.name for audit in target_edge if audit.has_pair_separating_current
    ]
    sep_exceptions_closed_by_pair_sep = [
        audit.row.name for audit in target_sep if audit.has_pair_separating_current
    ]
    sep_exceptions_closed_by_pair_edge = [
        audit.row.name for audit in target_sep if audit.has_pair_edge_cut
    ]

    pair_count_hist = Counter(len(audit.pair_cuts) for audit in audits)
    mirror_pair_count_hist = Counter(
        sum(1 for cut in audit.pair_cuts if cut.exact_mirror_pair) for audit in audits
    )
    pair_label_hist = Counter(cut.label_count for cut in all_pairs)
    pair_touch_label_hist = Counter(len(cut.labels_touching_dead) for cut in all_pairs)
    pair_edge_hist = Counter(cut.removed_edges for cut in all_pairs)
    pair_largest_drop_hist = Counter(cut.largest_drop for cut in all_pairs)
    pair_component_gain_hist = Counter(cut.component_gain for cut in all_pairs)
    pair_current_hist = Counter((cut.b0_labels, cut.b1_labels) for cut in all_pairs)
    pair_endpoint_hist = Counter(cut.endpoint_pair for cut in all_pairs)
    mirror_endpoint_hist = Counter(cut.endpoint_pair for cut in all_pairs if cut.exact_mirror_pair)
    strongest_by_row = [strongest_pair(audit.pair_cuts) for audit in audits]
    strongest_present = [cut for cut in strongest_by_row if cut is not None]
    strongest_pair_edge_hist = Counter(cut.removed_edges for cut in strongest_present)
    strongest_pair_largest_drop_hist = Counter(cut.largest_drop for cut in strongest_present)
    strongest_pair_component_gain_hist = Counter(cut.component_gain for cut in strongest_present)
    strongest_pair_current_hist = Counter((cut.b0_labels, cut.b1_labels) for cut in strongest_present)

    exact_mirror_pairs = [cut for cut in all_pairs if cut.exact_mirror_pair]
    exact_mirror_edge = [cut for cut in exact_mirror_pairs if cut.is_edge_cut]
    exact_mirror_sep = [cut for cut in exact_mirror_pairs if cut.is_separating_current]
    single_edge_not_pair_edge = [
        audit.row.name
        for audit in audits
        if audit.has_single_edge_cut and not audit.has_pair_edge_cut
    ]
    single_sep_not_pair_sep = [
        audit.row.name
        for audit in audits
        if audit.has_single_separating_current and not audit.has_pair_separating_current
    ]
    pair_edge_not_single_edge = [
        audit.row.name
        for audit in audits
        if audit.has_pair_edge_cut and not audit.has_single_edge_cut
    ]
    pair_sep_not_single_sep = [
        audit.row.name
        for audit in audits
        if audit.has_pair_separating_current and not audit.has_single_separating_current
    ]

    top_pairs = sorted(
        all_pairs,
        key=lambda cut: (
            -cut.removed_edges,
            -cut.largest_drop,
            -cut.component_gain,
            not cut.exact_mirror_pair,
            cut.row_name,
        ),
    )[:12]
    rows_with_gate_edge_support = [
        audit for audit in audits if audit.gate_edge_support_labels
    ]
    rows_without_gate_edge_support = [
        audit.row.name for audit in audits if not audit.gate_edge_support_labels
    ]
    edge_support_label_hist = Counter(len(audit.edge_support_labels) for audit in audits)
    gate_touch_label_union_hist = Counter(len(audit.gate_touch_labels) for audit in audits)
    gate_edge_support_label_hist = Counter(
        len(audit.gate_edge_support_labels) for audit in audits
    )

    score_hist, path = tournament()

    print("HYP-3476 DEAD-COVER PAIR-CURRENT EXCEPTION AUDIT")
    print("status=EVIDENCE / two-gate graph-current audit; not an LRC14 proof")
    print("source=HYP-3472 dead-cover projection + HYP-3475 mirror-orbit rule")
    print()
    print("## Assumption Challenge")
    print("alternate vertices considered: runners, gaps, fixed circle sections,")
    print("section boundaries, wall-crossing events, residues, cover arcs, Fourier")
    print("modes, matroid circuits, survivor gates, mirror orbits, gate pairs,")
    print("blocker labels, graph cuts, and formal proof obligations.")
    print("chosen carrier vertices: low-rank E/branch gate pairs and pair-current")
    print("cut obligations in the dead-cover blocker projection.")
    print("preserved predicate: target row has a <=2-gate boundary-current carrier")
    print("after deleting adjacent blocker labels.")
    print("destroyed if scalarized: interval geometry, ordered branch orientation,")
    print("and typed endpoint sidecars; HYP-3474 warns against count/color shadows.")
    print()
    print("## Aggregate Pair-Current Ledger")
    print(f"bank_rows={len(row_specs)}")
    print(f"bank_dead_rows={bank_dead_rows}")
    print(f"audit_scope={scope}")
    print(f"audited_dead_rows={len(audits)}")
    print(f"audited_single_e_branch_cuts={len(all_single)}")
    print(f"audited_pair_e_branch_cuts={len(all_pairs)}")
    print(f"exact_mirror_pairs={len(exact_mirror_pairs)}")
    print(f"exact_mirror_pairs_with_edge_cut={len(exact_mirror_edge)}")
    print(f"exact_mirror_pairs_with_separating_current={len(exact_mirror_sep)}")
    print(f"pair_count_hist={top_items(pair_count_hist)}")
    print(f"mirror_pair_count_hist={top_items(mirror_pair_count_hist)}")
    print()
    print("## Edge-Support Label Gap")
    print(
        "audited_rows_with_gate_touching_edge_support_label="
        f"{len(rows_with_gate_edge_support)}/{len(audits)}"
    )
    print(
        "audited_rows_without_gate_touching_edge_support_label="
        f"{rows_without_gate_edge_support}"
    )
    print(f"edge_support_label_count_hist={top_items(edge_support_label_hist)}")
    print(f"gate_touch_label_union_count_hist={top_items(gate_touch_label_union_hist)}")
    print(f"gate_edge_support_label_count_hist={top_items(gate_edge_support_label_hist)}")
    print()
    print("## Coverage Upgrade")
    print("reference_hyp3472_full_bank_single_edge_cut=123/130")
    print("reference_hyp3472_full_bank_single_separating_current=121/130")
    print(f"audited_rows_with_single_projection_edge_cut={len(rows_with_single_edge)}/{len(audits)}")
    print(f"audited_rows_with_pair_projection_edge_cut={len(rows_with_pair_edge)}/{len(audits)}")
    print(f"audited_rows_with_single_separating_current={len(rows_with_single_sep)}/{len(audits)}")
    print(f"audited_rows_with_pair_separating_current={len(rows_with_pair_sep)}/{len(audits)}")
    print(f"audited_rows_with_mirror_pair_projection_edge_cut={len(rows_with_mirror_edge)}/{len(audits)}")
    print(f"audited_rows_with_mirror_pair_separating_current={len(rows_with_mirror_sep)}/{len(audits)}")
    print(f"single_edge_not_pair_edge={single_edge_not_pair_edge}")
    print(f"single_separating_not_pair_separating={single_sep_not_pair_sep}")
    print(f"pair_edge_not_single_edge={pair_edge_not_single_edge}")
    print(f"pair_separating_not_single_separating={pair_sep_not_single_sep}")
    print()
    print("## HYP-3472 Exception Closure")
    print(f"target_edge_exception_rows={list(EDGE_EXCEPTION_ROWS)}")
    print(
        "edge_exception_rows_closed_by_pair_edge_cut="
        f"{len(edge_exceptions_closed_by_pair_edge)}/{len(target_edge)} "
        f"{edge_exceptions_closed_by_pair_edge}"
    )
    print(
        "edge_exception_rows_closed_by_pair_separating_current="
        f"{len(edge_exceptions_closed_by_pair_sep)}/{len(target_edge)} "
        f"{edge_exceptions_closed_by_pair_sep}"
    )
    print(
        "edge_exception_rows_remaining_without_pair_edge_cut="
        f"{[audit.row.name for audit in target_edge if not audit.has_pair_edge_cut]}"
    )
    print(
        "edge_exception_rows_remaining_without_pair_separating_current="
        f"{[audit.row.name for audit in target_edge if not audit.has_pair_separating_current]}"
    )
    print(f"target_separating_exception_rows={list(SEPARATING_EXCEPTION_ROWS)}")
    print(
        "separating_exception_rows_closed_by_pair_separating_current="
        f"{len(sep_exceptions_closed_by_pair_sep)}/{len(target_sep)} "
        f"{sep_exceptions_closed_by_pair_sep}"
    )
    print(
        "separating_exception_rows_closed_by_pair_edge_cut="
        f"{len(sep_exceptions_closed_by_pair_edge)}/{len(target_sep)} "
        f"{sep_exceptions_closed_by_pair_edge}"
    )
    print(
        "separating_exception_rows_remaining_without_pair_separating_current="
        f"{[audit.row.name for audit in target_sep if not audit.has_pair_separating_current]}"
    )
    print(
        "h3472_edge_exceptions_h3475_hard_mirror_intersection="
        f"{sorted(set(EDGE_EXCEPTION_ROWS) & H3475_HARD_MIRROR_ROWS)}"
    )
    print(
        "h3476_targets_h3475_hard_mirror_intersection="
        f"{sorted(set(SEPARATING_EXCEPTION_ROWS) & H3475_HARD_MIRROR_ROWS)}"
    )
    print()
    print("## Pair Fingerprints")
    print(f"pair_label_count_hist={top_items(pair_label_hist)}")
    print(f"pair_touch_label_count_hist={top_items(pair_touch_label_hist)}")
    print(f"pair_removed_edges_hist={top_items(pair_edge_hist)}")
    print(f"pair_largest_drop_hist={top_items(pair_largest_drop_hist)}")
    print(f"pair_component_gain_hist={top_items(pair_component_gain_hist)}")
    print(f"pair_current_hist={top_items(pair_current_hist)}")
    print(f"pair_endpoint_kind_hist={top_items(pair_endpoint_hist, 14)}")
    print(f"mirror_pair_endpoint_kind_hist={top_items(mirror_endpoint_hist, 14)}")
    print()
    print("## Strongest Pair Per Row")
    print(f"strongest_pair_removed_edges_hist={top_items(strongest_pair_edge_hist)}")
    print(f"strongest_pair_largest_drop_hist={top_items(strongest_pair_largest_drop_hist)}")
    print(f"strongest_pair_component_gain_hist={top_items(strongest_pair_component_gain_hist)}")
    print(f"strongest_pair_current_hist={top_items(strongest_pair_current_hist)}")
    print()
    print("## Target Row Detail")
    for name in SEPARATING_EXCEPTION_ROWS:
        audit = audit_by_name[name]
        print(f"row={name}")
        print(
            "  gates={gates} pairs={pairs} mirror_pairs={mirror_pairs} "
            "single_edge={single_edge} pair_edge={pair_edge} "
            "single_sep={single_sep} pair_sep={pair_sep}".format(
                gates=len(audit.e_branch_gates),
                pairs=len(audit.pair_cuts),
                mirror_pairs=sum(1 for cut in audit.pair_cuts if cut.exact_mirror_pair),
                single_edge=audit.has_single_edge_cut,
                pair_edge=audit.has_pair_edge_cut,
                single_sep=audit.has_single_separating_current,
                pair_sep=audit.has_pair_separating_current,
            )
        )
        print(
            "  projection_baseline="
            f"components:{audit.original.components} "
            f"largest:{audit.original.largest} "
            f"edges:{audit.original.edges}"
        )
        print(f"  dead_label_support={top_items(audit.dead_label_count)}")
        print(f"  edge_support_labels={audit.edge_support_labels}")
        print(f"  gate_touch_labels={audit.gate_touch_labels}")
        print(f"  gate_edge_support_labels={audit.gate_edge_support_labels}")
        print_single(best_single(audit.single_cuts), "  ")
        print("  minimal_edge_pair:")
        print_pair(minimal_edge_pair(audit.pair_cuts), "    ")
        print("  minimal_separating_pair:")
        print_pair(minimal_separating_pair(audit.pair_cuts), "    ")
        print("  strongest_pair:")
        print_pair(strongest_pair(audit.pair_cuts), "    ")
    print()
    print("## Top Pair Cuts")
    for cut in top_pairs:
        print_pair(cut, "  ")
    print()
    print("## Lemma Target")
    print("Finite audited form:")
    print("  HYP-3472's single-gate exception set can be replaced by the rows that")
    print("  still lack a two-gate projection cut or separating current above.")
    print("  If the remaining set is empty, HYP-3473 can expose a <=2-gate terminal")
    print("  boundary-current packet.  If not, the remaining rows are the sharper")
    print("  finite debt for HYP-3455/HYP-3451/HYP-3460/two-adic/SPEC routes.")
    print()
    print("## Tournament Analysis")
    print("vertices=pair-current proof carriers, not runners, arcs, or raw counts")
    print(
        "pairwise_observable=predicate retention + pair current payload + mirror "
        "payload + cut strength + formal packet fit + quotient guardrail"
    )
    print("switch=higher proof-facing carrier score; ties use declared route order")
    print(f"axes={','.join(AXES)}")
    print(f"score_hist={score_hist}")
    print("directed_3cycles=0")
    print("hamiltonian_path=" + " -> ".join(path))


if __name__ == "__main__":
    main()
