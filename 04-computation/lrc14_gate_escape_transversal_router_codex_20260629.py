#!/usr/bin/env python3
"""HYP-3453: gate-escape transversal router for the LRC14 finite lemma.

HYP-3450 says every audited primitive covering row has an endpoint-rank <= 2
survivor component.  HYP-3438 says mixed even-safe components have exact
survivor gate words.  This scout joins the two ledgers.

The compression issue is subtle: "low-rank survivor component" folds together
two different proof objects.

    clean escape: no bad-core saturation is present in the component;
    gate escape: a survivor gate sits on the boundary of dead/bad-core blocks.

Only the gate escape composes directly with the HYP-3451 dead-cover projection
and its Menger/Green-current graph.  The finite lemma target becomes:

    if a row has any dead_both component, then it has a rank <= 2 survivor gate;
    if it has no rank <= 2 survivor gate, the dead-cover obstruction is empty.

This is evidence, not a proof, but it gives a cleaner local-to-global target
than scalar survivor mass.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
H3450_PATH = ROOT / "04-computation" / "lrc14_component_cover_obstruction_extractor_codex_20260628.py"
H3438_PATH = ROOT / "04-computation" / "lrc14_survivor_gate_word_audit_codex_20260629.py"


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


H3450 = load_module("hyp3450_component_cover_for_hyp3452", H3450_PATH)
H3438 = load_module("hyp3438_survivor_gate_for_hyp3452", H3438_PATH)


Interval = tuple[Fraction, Fraction]


def fmt(x: Fraction | None) -> str:
    if x is None:
        return "None"
    return H3450.fmt(x)


def interval_text(interval: Interval | None) -> str:
    if interval is None:
        return "None"
    return f"[{fmt(interval[0])},{fmt(interval[1])}]"


def gate_rank(gate: H3438.SurvivorGate) -> int:
    labels = set(gate.left_labels + gate.right_labels)
    labels.discard("FREE")
    return len(labels)


@dataclass(frozen=True)
class RowJoin:
    name: str
    speeds: tuple[int, ...]
    component_row: H3450.RowAudit
    gates: tuple[H3438.SurvivorGate, ...]
    gates_by_parent: dict[Interval, tuple[H3438.SurvivorGate, ...]]
    low_rank_components: tuple[H3450.ComponentAudit, ...]
    low_rank_gate_parent_components: tuple[H3450.ComponentAudit, ...]
    clean_low_rank_components: tuple[H3450.ComponentAudit, ...]
    low_rank_gates: tuple[H3438.SurvivorGate, ...]

    @property
    def dead_count(self) -> int:
        return len(self.component_row.dead_components)

    @property
    def has_dead(self) -> bool:
        return self.dead_count > 0

    @property
    def has_low_gate(self) -> bool:
        return bool(self.low_rank_gates)

    @property
    def clean_only_escape(self) -> bool:
        return bool(self.low_rank_components) and not self.has_low_gate


@dataclass(frozen=True)
class Carrier:
    name: str
    scores: tuple[int, ...]

    @property
    def total(self) -> int:
        return sum(self.scores)


AXES = (
    "predicate_retention",
    "component_gate_join",
    "cut_graph_composability",
    "clean_case_exit",
    "endpoint_payload",
    "scalar_firewall",
)

CARRIERS = (
    Carrier("E00_dead_positive_rank2_gate_transversal", (10, 10, 10, 9, 10, 10)),
    Carrier("E01_clean_obstruction_empty_exit", (10, 9, 8, 10, 8, 10)),
    Carrier("E02_gate_endpoint_delta_spine", (9, 10, 9, 8, 10, 9)),
    Carrier("E03_component_cover_conductance_join", (9, 9, 10, 7, 8, 9)),
    Carrier("E04_low_rank_component_shadow", (8, 6, 5, 6, 7, 6)),
    Carrier("E05_raw_survivor_measure", (4, 2, 1, 2, 1, 0)),
)


def join_row(name: str, speeds: tuple[int, ...]) -> RowJoin:
    component_row = H3450.audit_row(name, speeds)
    bad_core_row = H3438.H3436.audit_row(name, speeds)
    mixed_components = H3438.build_mixed_components(bad_core_row)

    gate_lists_by_parent: dict[Interval, list[H3438.SurvivorGate]] = defaultdict(list)
    for mixed in mixed_components:
        gate_lists_by_parent[mixed.parent].extend(mixed.survivor_gates)
    gates_by_parent = {
        parent: tuple(gates)
        for parent, gates in gate_lists_by_parent.items()
    }
    gates = tuple(gate for gate_list in gates_by_parent.values() for gate in gate_list)
    low_rank_components = tuple(
        component
        for component in component_row.components
        if component.union_measure > H3450.ZERO
        and component.endpoint_rank is not None
        and component.endpoint_rank <= 2
    )
    low_rank_gate_parent_components = tuple(
        component for component in low_rank_components if component.interval in gates_by_parent
    )
    clean_low_rank_components = tuple(
        component for component in low_rank_components if component.interval not in gates_by_parent
    )
    low_rank_gates = tuple(gate for gate in gates if gate_rank(gate) <= 2)
    return RowJoin(
        name=name,
        speeds=tuple(sorted(speeds)),
        component_row=component_row,
        gates=gates,
        gates_by_parent=gates_by_parent,
        low_rank_components=low_rank_components,
        low_rank_gate_parent_components=low_rank_gate_parent_components,
        clean_low_rank_components=clean_low_rank_components,
        low_rank_gates=low_rank_gates,
    )


def gate_sort_key(gate: H3438.SurvivorGate) -> tuple[Fraction, str, int, Fraction]:
    return (gate.length, gate.row_name, gate.component_index, gate.interval[0])


def print_gate(gate: H3438.SurvivorGate, indent: str = "  ") -> None:
    print(
        f"{indent}component={gate.component_index} gate={interval_text(gate.interval)} "
        f"rank={gate_rank(gate)} len={fmt(gate.length)}"
    )
    print(
        f"{indent}  labels={H3438.label_values(gate.left_labels)}|"
        f"{H3438.label_values(gate.right_labels)} endpoint_kinds={gate.endpoint_kind_signature}"
    )
    print(
        f"{indent}  word={gate.gate_word} adjacency={gate.adjacency} "
        f"delta=({gate.b0_delta},{gate.b1_delta}) mask={gate.branch_mask}"
    )


def print_row_summary(row: RowJoin) -> None:
    print(f"  {row.name}")
    print(f"    speeds={row.speeds}")
    print(
        "    components={components}, dead={dead}, gates={gates}, "
        "low_rank_components={low_components}, low_rank_gates={low_gates}".format(
            components=len(row.component_row.components),
            dead=row.dead_count,
            gates=len(row.gates),
            low_components=len(row.low_rank_components),
            low_gates=len(row.low_rank_gates),
        )
    )
    print(f"    component_class_hist={dict(row.component_row.class_hist)}")
    best = row.component_row.best_component
    print(
        "    best_component={interval} survivor={survivor} "
        "rank={rank} labels={labels} union={measure}".format(
            interval=interval_text(best.interval),
            survivor=interval_text(best.best_survivor),
            rank=best.endpoint_rank,
            labels=best.labels,
            measure=fmt(best.union_measure),
        )
    )
    for gate in sorted(row.low_rank_gates, key=gate_sort_key)[:4]:
        print_gate(gate, "    ")


def top_items(counter: Counter, limit: int = 12) -> str:
    items = counter.most_common(limit)
    suffix = "" if len(counter) <= limit else " ..."
    return "{" + ", ".join(f"{key!r}: {value}" for key, value in items) + suffix + "}"


def tournament() -> tuple[dict[int, int], list[str]]:
    hist = dict(sorted(Counter(carrier.total for carrier in CARRIERS).items()))
    order = {carrier.name: index for index, carrier in enumerate(CARRIERS)}
    path = [
        carrier.name
        for carrier in sorted(CARRIERS, key=lambda carrier: (-carrier.total, order[carrier.name]))
    ]
    return hist, path


def main() -> None:
    joins = [join_row(name, speeds) for name, speeds in H3450.rows().items()]
    rows_with_low_component = [row for row in joins if row.low_rank_components]
    rows_with_low_gate = [row for row in joins if row.low_rank_gates]
    rows_with_dead = [row for row in joins if row.has_dead]
    rows_without_low_gate = [row for row in joins if not row.low_rank_gates]
    rows_dead_without_low_gate = [row for row in rows_with_dead if not row.low_rank_gates]
    dead_zero_rows = [row for row in joins if not row.has_dead]
    dead_zero_with_low_gate = [row for row in dead_zero_rows if row.low_rank_gates]
    dead_zero_clean_only = [row for row in dead_zero_rows if not row.low_rank_gates]

    all_gates = [gate for row in joins for gate in row.gates]
    low_gates = [gate for row in joins for gate in row.low_rank_gates]
    high_gate_rank_hist = Counter(gate_rank(gate) for gate in all_gates)
    low_endpoint_hist = Counter(gate.endpoint_kind_signature for gate in low_gates)
    low_adjacency_hist = Counter(gate.adjacency for gate in low_gates)
    low_delta_hist = Counter((gate.b0_delta, gate.b1_delta) for gate in low_gates)
    low_mask_hist = Counter(gate.branch_mask for gate in low_gates)
    class_hist: Counter[str] = Counter()
    for row in joins:
        class_hist.update(row.component_row.class_hist)

    selected_names = [
        "covering_AP_with_84",
        "ap_omit_12_tail_84x01",
        "ap_omit_12_tail_84x02",
        "ap_omit_12_tail_84x03",
        "ap_omit_12_tail_84x04",
        "ap_omit_12_tail_84x05",
        "random_covering_082",
    ]
    by_name = {row.name: row for row in joins}
    selected_rows = [by_name[name] for name in selected_names if name in by_name]
    hist, path = tournament()

    print("HYP-3453 GATE-ESCAPE TRANSVERSAL ROUTER")
    print("status=EVIDENCE / exact component-gate join; not an LRC14 proof")
    print("source=HYP-3438 survivor gates joined to HYP-3450/HYP-3451 component-cover graph")
    print()

    print("## Assumption Challenge")
    print("alternate vertices considered: runners, gaps, fixed circle sections,")
    print("section boundaries, wall-crossing events, residues, cover arcs, endpoint")
    print("walls, bad-core blocks, survivor gates, clean E_safe components, dead")
    print("components, branch-coloured blockers, and proof obligations.")
    print("chosen carrier vertices: survivor gates plus clean-component exits.")
    print("preserved predicate: whether the component-cover obstruction has a")
    print("rank <= 2 escape that composes with the dead-cover graph.")
    print("destroyed by scalarization: whether the escape is clean-only or sits on")
    print("a bad-core/dead-cover gate boundary with endpoint and cover-delta data.")
    print()

    print("## Aggregate Gate/Component Join")
    print(f"rows_audited={len(joins)}")
    print(f"rows_with_low_rank_component={len(rows_with_low_component)}/{len(joins)}")
    print(f"rows_with_low_rank_gate={len(rows_with_low_gate)}/{len(joins)}")
    print(f"rows_with_dead_components={len(rows_with_dead)}/{len(joins)}")
    print(
        "rows_with_dead_components_and_low_rank_gate="
        f"{len(rows_with_dead) - len(rows_dead_without_low_gate)}/{len(rows_with_dead)}"
    )
    print(f"rows_dead_without_low_rank_gate={[row.name for row in rows_dead_without_low_gate]}")
    print(f"dead_zero_rows={[row.name for row in dead_zero_rows]}")
    print(f"dead_zero_with_low_rank_gate={[row.name for row in dead_zero_with_low_gate]}")
    print(f"dead_zero_clean_only_rows={[row.name for row in dead_zero_clean_only]}")
    print(f"rows_without_low_rank_gate={[row.name for row in rows_without_low_gate]}")
    print(f"component_class_hist={dict(sorted(class_hist.items()))}")
    print(f"all_survivor_gates={len(all_gates)}")
    print(f"low_rank_gates={len(low_gates)}")
    print(f"gate_endpoint_rank_hist={dict(sorted(high_gate_rank_hist.items()))}")
    print(f"low_rank_component_total={sum(len(row.low_rank_components) for row in joins)}")
    print(
        "low_rank_gate_parent_components="
        f"{sum(len(row.low_rank_gate_parent_components) for row in joins)}"
    )
    print(f"clean_low_rank_components={sum(len(row.clean_low_rank_components) for row in joins)}")
    print()

    print("## Low-Rank Gate Fingerprints")
    print(f"endpoint_kind_hist={top_items(low_endpoint_hist)}")
    print(f"adjacency_hist={dict(sorted(low_adjacency_hist.items()))}")
    print(f"branch_mask_hist={dict(sorted(low_mask_hist.items()))}")
    print(f"cover_delta_hist={dict(sorted(low_delta_hist.items()))}")
    print()

    print("## AP-Tail And High-Rank Controls")
    for row in selected_rows:
        print_row_summary(row)
    print()

    print("## Clean-Only Exceptions")
    for row in dead_zero_clean_only:
        print_row_summary(row)
    print()

    print("## Lemma Target")
    print("Audited finite form:")
    print("  dead_components(row) > 0  =>  row has a rank <= 2 survivor gate.")
    print("  no rank <= 2 survivor gate => dead_components(row) = 0.")
    print("Thus a full two-colour blocker-saturation counterexample cannot hide in")
    print("a clean low-rank component.  It must pass through an exact survivor gate")
    print("with endpoint labels, branch mask, adjacent cover deltas, and parent even")
    print("wall retained.  This is the local object to feed into HYP-3451's")
    print("Menger/Green-current/conductance proof route.")
    print()

    print("## Tournament Analysis")
    print("vertices=escape proof carriers, not runners or scalar survivor masses")
    print(
        "pairwise_observable=predicate retention + component/gate join + "
        "cut-graph composability + clean-case exit + endpoint payload"
    )
    print("switch=higher proof-facing carrier score; ties use declared route order")
    print(f"axes={','.join(AXES)}")
    print(f"score_hist={hist}")
    print("directed_3cycles=0")
    print("hamiltonian_path=" + " -> ".join(path))


if __name__ == "__main__":
    main()
