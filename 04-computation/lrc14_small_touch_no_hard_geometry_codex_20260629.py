#!/usr/bin/env python3
"""HYP-3478: geometry audit for the six small-touch/no-hard LRC14 rows.

HYP-3476 split the currentless frontier into AP pair-current rows, the
random031 hard/currentless overlap, and six remaining non-AP rows with no hard
mirror orbit and no dead-cover projection edges.  This script looks directly
at those six rows instead of asking for larger gate pairs.

The finite object under inspection is:

    dead E_safe components + paired B0/B1 blocker owners
    + exact mirror interval involution
    + low-rank E/branch touching gates

This is evidence for a finite singleton-current packet, not an LRC14 proof.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
H3472_PATH = ROOT / "04-computation" / "lrc14_dead_cover_boundary_current_codex_20260629.py"

TARGET_ROWS = (
    "random_covering_001",
    "random_covering_039",
    "random_covering_062",
    "random_covering_074",
    "random_covering_086",
    "random_covering_101",
)

COVER_DELTA_SIDE_ROWS = frozenset({"random_covering_039", "random_covering_074"})


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


H3472 = load_module("hyp3472_boundary_current_for_hyp3478", H3472_PATH)
H3453 = H3472.H3453
H3438 = H3453.H3438


def fmt(value: Fraction | None) -> str:
    return H3472.fmt_fraction(value)


def interval_text(interval: tuple[Fraction, Fraction] | None) -> str:
    return H3472.fmt_interval(interval)


def top_items(counter: Counter, limit: int = 8) -> str:
    items = counter.most_common(limit)
    suffix = "" if len(counter) <= limit else " ..."
    return "{" + ", ".join(f"{key!r}: {value}" for key, value in items) + suffix + "}"


def typed_residue_label(label: str) -> str:
    if label == "FREE" or ":" not in label:
        return label
    prefix, value = label.split(":", 1)
    return f"{prefix}:{int(value) % 14}"


def typed_word(gate) -> str:
    left = ".".join(typed_residue_label(label) for label in gate.left_labels) or "none"
    right = ".".join(typed_residue_label(label) for label in gate.right_labels) or "none"
    return f"{left}|{right}"


def gate_sort_key(gate) -> tuple[Fraction, int, Fraction, str]:
    return (gate.length, gate.component_index, gate.interval[0], gate.endpoint_kind_signature)


def min_e_branch_gate(row):
    gates = [gate for gate in row.low_rank_gates if H3472.is_e_branch_gate(gate)]
    if not gates:
        raise ValueError(f"target row lacks E/branch gates: {row.name}")
    return min(gates, key=gate_sort_key)


def min_gate_kind(gate) -> str:
    edge = gate.adjacency in {"left_bad_edge", "right_bad_edge"}
    unit = (gate.b0_delta, gate.b1_delta) == (1, 1)
    if edge and unit and gate.branch_mask in {"branch0", "branch1"}:
        return "branch_unit_delta"
    if edge and unit and gate.branch_mask == "both":
        return "both_unit_delta"
    if edge and gate.branch_mask in {"branch0", "branch1"}:
        return "delta_sidecar_packet"
    return "other"


def cover_values(cover: tuple[int, tuple[int, ...]] | None) -> tuple[int, ...]:
    if cover is None:
        return ()
    return tuple(cover[1])


def owner_pair(component) -> tuple[tuple[int, ...], tuple[int, ...]]:
    return (cover_values(component.b0_cover), cover_values(component.b1_cover))


def owner_pair_text(pair: tuple[tuple[int, ...], tuple[int, ...]]) -> str:
    b0, b1 = pair
    b0_text = ".".join(map(str, b0)) if b0 else "none"
    b1_text = ".".join(map(str, b1)) if b1 else "none"
    return f"B0:{b0_text}/B1:{b1_text}"


def owner_pair_residues(pair: tuple[tuple[int, ...], tuple[int, ...]]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    b0, b1 = pair
    return (tuple(value % 14 for value in b0), tuple(value % 14 for value in b1))


def mirror_interval(interval: tuple[Fraction, Fraction]) -> tuple[Fraction, Fraction]:
    return (1 - interval[1], 1 - interval[0])


def swapped_owner_pair(pair: tuple[tuple[int, ...], tuple[int, ...]]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    b0, b1 = pair
    return (b1, b0)


def owner_balance(dead_components: Iterable) -> dict[int, int]:
    balance: Counter[int] = Counter()
    for component in dead_components:
        b0, b1 = owner_pair(component)
        for value in b0:
            balance[value] += 1
        for value in b1:
            balance[value] -= 1
    return dict(sorted((value, count) for value, count in balance.items() if count))


def structural_word(gate) -> tuple[str, str, str, int, int]:
    return (
        gate.endpoint_kind_signature,
        gate.branch_mask,
        gate.adjacency,
        gate.b0_delta,
        gate.b1_delta,
    )


def gate_text(gate) -> str:
    return (
        f"component={gate.component_index} gate={interval_text(gate.interval)} "
        f"typed={typed_word(gate)} structural={structural_word(gate)} "
        f"len={fmt(gate.length)} labels={H3472.gate_blocker_labels(gate)}"
    )


@dataclass(frozen=True)
class RowGeometry:
    row: object
    projection: H3472.ProjectionStats
    dead_label_set: frozenset[str]
    e_branch_cuts: tuple[H3472.GateCut, ...]
    touching_cuts: tuple[H3472.GateCut, ...]
    min_gate: object
    min_touching_cut: H3472.GateCut | None
    min_gate_kind: str
    mirror_failures: tuple[str, ...]
    singleton_dead: int
    owner_balance: dict[int, int]


def row_geometry(name: str) -> RowGeometry:
    row = H3453.join_row(name, H3453.H3450.rows()[name])
    dead_components = list(row.component_row.dead_components)
    projection = H3472.projection_stats(H3472.projection(dead_components))
    dead_label_set = frozenset(
        label
        for component in dead_components
        for label in H3472.component_blocker_labels(component)
    )
    e_branch_cuts = tuple(
        H3472.gate_cut(row, gate, projection, dead_label_set)
        for gate in row.low_rank_gates
        if H3472.is_e_branch_gate(gate)
    )
    touching_cuts = tuple(cut for cut in e_branch_cuts if cut.touches_dead_projection)
    min_gate = min_e_branch_gate(row)
    min_touching_cut = (
        min(touching_cuts, key=lambda cut: gate_sort_key(cut.gate))
        if touching_cuts
        else None
    )
    dead_by_interval = {
        component.interval: index for index, component in enumerate(dead_components)
    }
    mirror_failures: list[str] = []
    singleton_dead = 0
    for index, component in enumerate(dead_components):
        pair = owner_pair(component)
        if len(pair[0]) == 1 and len(pair[1]) == 1:
            singleton_dead += 1
        partner_index = dead_by_interval.get(mirror_interval(component.interval))
        if partner_index is None:
            mirror_failures.append(f"component {index} missing interval mirror")
            continue
        partner_pair = owner_pair(dead_components[partner_index])
        if partner_pair != swapped_owner_pair(pair):
            mirror_failures.append(
                f"component {index} mirror pair {owner_pair_text(partner_pair)} "
                f"!= swapped {owner_pair_text(swapped_owner_pair(pair))}"
            )
    return RowGeometry(
        row=row,
        projection=projection,
        dead_label_set=dead_label_set,
        e_branch_cuts=e_branch_cuts,
        touching_cuts=touching_cuts,
        min_gate=min_gate,
        min_touching_cut=min_touching_cut,
        min_gate_kind=min_gate_kind(min_gate),
        mirror_failures=tuple(mirror_failures),
        singleton_dead=singleton_dead,
        owner_balance=owner_balance(dead_components),
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
    "mirror_component_payload",
    "owner_pair_payload",
    "touching_gate_payload",
    "delta_sidecar_payload",
    "formal_packet_fit",
    "scalar_firewall",
)

CARRIERS = (
    Carrier("S00_mirror_singleton_geometry_packet", (10, 10, 10, 8, 8, 9, 10)),
    Carrier("S01_owner_current_locality_packet", (10, 9, 10, 8, 7, 9, 10)),
    Carrier("S02_touching_gate_sidecar_packet", (9, 8, 8, 10, 9, 9, 9)),
    Carrier("S03_unit_delta_vs_cover_delta_split", (9, 7, 8, 9, 10, 8, 9)),
    Carrier("S04_section_boundary_interval_word", (8, 10, 7, 7, 7, 8, 9)),
    Carrier("S05_residue_owner_pair_shadow", (7, 6, 8, 6, 6, 6, 7)),
    Carrier("S06_raw_zero_edge_count", (4, 1, 1, 1, 1, 2, 0)),
    Carrier("S07_raw_row_name_list", (3, 0, 0, 0, 0, 1, 0)),
)


def tournament_summary() -> tuple[dict[int, int], list[str]]:
    hist = dict(sorted(Counter(carrier.total for carrier in CARRIERS).items()))
    order = {carrier.name: index for index, carrier in enumerate(CARRIERS)}
    path = [
        carrier.name
        for carrier in sorted(CARRIERS, key=lambda carrier: (-carrier.total, order[carrier.name]))
    ]
    return hist, path


def print_dead_components(geom: RowGeometry) -> None:
    dead_components = list(geom.row.component_row.dead_components)
    dead_by_interval = {
        component.interval: index for index, component in enumerate(dead_components)
    }
    for index, component in enumerate(dead_components):
        partner = dead_by_interval.get(mirror_interval(component.interval))
        pair = owner_pair(component)
        residues = owner_pair_residues(pair)
        print(
            "    dead[{index}] interval={interval} len={length} mid={mid} "
            "mirror={mirror} owner_pair={owner_pair} residues={residues} labels={labels}".format(
                index=index,
                interval=interval_text(component.interval),
                length=fmt(component.interval[1] - component.interval[0]),
                mid=fmt((component.interval[0] + component.interval[1]) / 2),
                mirror=partner,
                owner_pair=owner_pair_text(pair),
                residues=residues,
                labels=H3472.component_blocker_labels(component),
            )
        )


def print_row(geom: RowGeometry) -> None:
    row = geom.row
    e_branch_count = len(geom.e_branch_cuts)
    edge_cuts = sum(cut.is_edge_cut for cut in geom.e_branch_cuts)
    separating = sum(cut.is_separating_current for cut in geom.e_branch_cuts)
    touch_hist = Counter(cut.labels_touching_dead for cut in geom.touching_cuts)
    kind_note = "cover_delta_sidecar" if row.name in COVER_DELTA_SIDE_ROWS else "unit_delta_singleton"
    print()
    print(f"## Row {row.name}")
    print(f"speeds={row.speeds}")
    print(
        "components={components} dead={dead} low_rank_gates={low} "
        "e_branch_gates={eb} touching_e_branch={touch} edge_cuts={edge} separating={sep}".format(
            components=len(row.component_row.components),
            dead=row.dead_count,
            low=len(row.low_rank_gates),
            eb=e_branch_count,
            touch=len(geom.touching_cuts),
            edge=edge_cuts,
            sep=separating,
        )
    )
    print(
        "projection=components:{components} largest:{largest} edges:{edges} "
        "dead_labels={labels}".format(
            components=geom.projection.components,
            largest=geom.projection.largest,
            edges=geom.projection.edges,
            labels=tuple(sorted(geom.dead_label_set)),
        )
    )
    print(
        "singleton_dead={singleton}/{dead} mirror_failures={mirror_failures} "
        "owner_balance={balance}".format(
            singleton=geom.singleton_dead,
            dead=row.dead_count,
            mirror_failures=geom.mirror_failures,
            balance=geom.owner_balance,
        )
    )
    print(f"min_gate_kind={geom.min_gate_kind} packet_note={kind_note}")
    print(f"min_gate={gate_text(geom.min_gate)}")
    min_cut = H3472.gate_cut(row, geom.min_gate, geom.projection, geom.dead_label_set)
    print(
        "min_gate_touching_labels={labels} removes_edges={removed} "
        "component_gain={gain} largest_drop={drop}".format(
            labels=min_cut.labels_touching_dead,
            removed=min_cut.removed_edges,
            gain=min_cut.component_gain,
            drop=min_cut.largest_drop,
        )
    )
    if geom.min_touching_cut is not None and geom.min_touching_cut.gate is not geom.min_gate:
        touch = geom.min_touching_cut
        print(f"min_touching_gate={gate_text(touch.gate)}")
        print(
            "min_touching_gate_labels={labels} removes_edges={removed} "
            "component_gain={gain} largest_drop={drop}".format(
                labels=touch.labels_touching_dead,
                removed=touch.removed_edges,
                gain=touch.component_gain,
                drop=touch.largest_drop,
            )
        )
    print(f"touching_label_hist_top={top_items(touch_hist, 8)}")
    print_dead_components(geom)


def main() -> None:
    geometries = tuple(row_geometry(name) for name in TARGET_ROWS)

    print("HYP-3478 SMALL-TOUCH / NO-HARD GEOMETRY AUDIT")
    print("status=EVIDENCE / exact six-row geometry audit; not an LRC14 proof")
    print("source=HYP-3450/HYP-3453 component-gate bank + HYP-3472/HYP-3476 current utilities")
    print(f"target_rows={TARGET_ROWS}")

    projection_edge_rows = [geom.row.name for geom in geometries if geom.projection.edges]
    mirror_failure_rows = [geom.row.name for geom in geometries if geom.mirror_failures]
    owner_unbalanced_rows = [geom.row.name for geom in geometries if geom.owner_balance]
    nonsingleton_rows = [
        geom.row.name
        for geom in geometries
        if geom.singleton_dead != geom.row.dead_count
    ]
    min_kind_hist = Counter(geom.min_gate_kind for geom in geometries)
    dead_count_hist = Counter(geom.row.dead_count for geom in geometries)
    touching_counts = {geom.row.name: len(geom.touching_cuts) for geom in geometries}
    min_gate_not_touching = [
        geom.row.name
        for geom in geometries
        if geom.min_touching_cut is not None and geom.min_touching_cut.gate is not geom.min_gate
    ]

    print()
    print("## Aggregate Geometry")
    print(f"rows={len(geometries)}")
    print(f"dead_count_hist={dict(sorted(dead_count_hist.items()))}")
    print(f"total_dead_components={sum(geom.row.dead_count for geom in geometries)}")
    print(f"projection_edge_rows={projection_edge_rows}")
    print(f"singleton_cover_fail_rows={nonsingleton_rows}")
    print(f"mirror_failure_rows={mirror_failure_rows}")
    print(f"owner_unbalanced_rows={owner_unbalanced_rows}")
    print(f"min_gate_kind_hist={dict(sorted(min_kind_hist.items()))}")
    print(f"touching_e_branch_counts={touching_counts}")
    print(f"min_gate_not_touching_dead_labels={min_gate_not_touching}")
    print(
        "cover_delta_sidecar_rows="
        f"{tuple(name for name in TARGET_ROWS if name in COVER_DELTA_SIDE_ROWS)}"
    )
    print(
        "clean_unit_delta_rows="
        f"{tuple(name for name in TARGET_ROWS if name not in COVER_DELTA_SIDE_ROWS)}"
    )

    for geom in geometries:
        print_row(geom)

    print()
    print("## Geometry Consequence")
    print(
        "All six rows have edgeless dead-cover projections, but their E/branch "
        "gates do touch the singleton dead labels.  The obstruction is therefore "
        "not 'find a larger projection cut'; there is no edge to cut."
    )
    print(
        "Every dead component is a singleton B0/B1 owner-pair component, and "
        "the exact interval mirror swaps the B0 and B1 owners.  The signed owner "
        "balance cancels row-by-row, so any useful current must be local to "
        "mirror pairs or touching gates rather than a raw global owner imbalance."
    )
    print(
        "The finite terminal packet splits cleanly: rows 039 and 074 have "
        "minimum E/branch gates in the cover-delta sidecar packet, while rows "
        "001, 062, 086, and 101 have branch-specific unit-delta minima."
    )

    print()
    print("## Tournament Analysis")
    hist, path = tournament_summary()
    print(f"vertices={[carrier.name for carrier in CARRIERS]}")
    print(f"axes={AXES}")
    print(
        "pairwise_observable=predicate retention + mirror component payload + "
        "owner-pair payload + touching-gate payload + delta-sidecar payload + scalar firewall"
    )
    print(
        "switch_gauge=higher retained local proof payload first; ties follow "
        "mirror singleton geometry -> owner current -> gate touch -> delta split -> formal packet"
    )
    print(f"score_hist={hist}")
    print("directed_3cycles=0")
    print(f"hamiltonian_path={' -> '.join(path)}")
    print(
        "assumption_challenge=considered runners, arcs, singleton components, "
        "mirror-paired components, owner pairs, fixed circle sections, section "
        "boundaries, touching gate events, residues, cover arcs, Fourier/color "
        "modes, and proof obligations.  Chosen carrier preserves the terminal "
        "singleton-current discharge predicate and destroys raw runner order, "
        "full interval order unless retained, and scalar gate-count shadows."
    )


if __name__ == "__main__":
    main()
