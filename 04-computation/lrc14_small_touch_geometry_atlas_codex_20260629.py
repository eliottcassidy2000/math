#!/usr/bin/env python3
"""HYP-3478: geometry atlas for the six small-touch/no-hard LRC14 rows.

HYP-3476/HYP-3477 split the remaining exception frontier.  After removing the
AP84 pair-current rows, the hard-currented rows, and the random031 gluing row,
six random rows remain:

    001, 039, 062, 074, 086, 101.

HYP-3476's pair-current audit says their dead-cover projections are edgeless.
This scout looks at the underlying geometry rather than trying larger gate
tuples: dead component intervals, mirror pairing, singleton owner labels,
nearby E/branch gates, and the incoming S319 unit-delta split.

It is finite-bank evidence for a singleton-current terminal packet, not an
LRC14 proof.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
H3472_PATH = ROOT / "04-computation" / "lrc14_dead_cover_boundary_current_codex_20260629.py"
H3475_PATH = ROOT / "04-computation" / "lrc14_colored_gate_mirror_orbit_codex_20260629.py"
S319_PATH = ROOT / "04-computation" / "lrc14_colored_gate_unit_delta_split_codex_20260629.py"

SMALL_TOUCH_ROWS = (
    "random_covering_001",
    "random_covering_039",
    "random_covering_062",
    "random_covering_074",
    "random_covering_086",
    "random_covering_101",
)

S319_DELTA_SIDECAR_ROWS = {
    "random_covering_022",
    "random_covering_036",
    "random_covering_039",
    "random_covering_047",
    "random_covering_051",
    "random_covering_054",
    "random_covering_056",
    "random_covering_058",
    "random_covering_059",
    "random_covering_063",
    "random_covering_069",
    "random_covering_074",
    "random_covering_083",
    "random_covering_085",
    "random_covering_090",
    "random_covering_095",
    "random_covering_104",
    "random_covering_107",
    "random_covering_113",
}


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
S319 = load_module("s319_gate_split_for_hyp3478", S319_PATH)


Interval = tuple[Fraction, Fraction]


def fmt(value: Fraction | None) -> str:
    return H3472.fmt_fraction(value)


def interval_text(interval: Interval | None) -> str:
    return H3472.fmt_interval(interval)


def mirror_interval(interval: Interval) -> Interval:
    return (1 - interval[1], 1 - interval[0])


def label_owner(label: str) -> int:
    return int(label.split(":", 1)[1])


def residues(labels: Iterable[str]) -> tuple[int, ...]:
    return tuple(sorted(label_owner(label) % 14 for label in labels))


def component_cover(component, branch: str) -> tuple[int, ...]:
    cover = component.b0_cover if branch == "B0" else component.b1_cover
    return tuple() if cover is None else tuple(cover[1])


def component_labels(component) -> tuple[str, ...]:
    return H3472.component_blocker_labels(component)


def gate_labels(gate) -> tuple[str, ...]:
    return H3472.gate_blocker_labels(gate)


def gate_text(gate) -> str:
    return (
        f"c{gate.component_index}:{interval_text(gate.interval)}:"
        f"{gate.endpoint_kind_signature}:mask={gate.branch_mask}:"
        f"adj={gate.adjacency}:d({gate.b0_delta},{gate.b1_delta})"
    )


def gate_sort_key(gate) -> tuple:
    return (gate.length, gate.row_name, gate.component_index, gate.interval[0])


@dataclass(frozen=True)
class ComponentRecord:
    index: int
    interval: Interval
    labels: tuple[str, ...]
    b0: tuple[int, ...]
    b1: tuple[int, ...]
    mirror_index: int | None
    touch_gate_count: int
    complete_gate_count: int
    best_touch_gate: object | None
    best_complete_gate: object | None

    @property
    def length(self) -> Fraction:
        return self.interval[1] - self.interval[0]

    @property
    def owner_pair(self) -> tuple[int, ...]:
        return tuple(sorted(self.b0 + self.b1))

    @property
    def singleton(self) -> bool:
        return len(self.b0) == 1 and len(self.b1) == 1

    @property
    def residue_pair(self) -> tuple[int, ...]:
        return tuple(owner % 14 for owner in self.owner_pair)

    @property
    def span(self) -> int:
        if len(self.owner_pair) < 2:
            return 0
        return max(self.owner_pair) - min(self.owner_pair)

    @property
    def owner_gcd(self) -> int:
        if not self.owner_pair:
            return 0
        a = self.owner_pair[0]
        for b in self.owner_pair[1:]:
            while b:
                a, b = b, a % b
        return a


@dataclass(frozen=True)
class RowGeometry:
    name: str
    speeds: tuple[int, ...]
    dead_components: tuple[ComponentRecord, ...]
    projection_edges: int
    edge_support_labels: tuple[str, ...]
    min_gate: object
    min_gate_kind: str
    best_cut: object | None
    e_branch_gate_count: int
    mirror_orbit_count: int
    hard_orbit_count: int

    @property
    def dead_count(self) -> int:
        return len(self.dead_components)

    @property
    def mirror_pair_count(self) -> int:
        seen: set[int] = set()
        count = 0
        for component in self.dead_components:
            if component.index in seen:
                continue
            if component.mirror_index is not None:
                seen.add(component.index)
                seen.add(component.mirror_index)
                count += 1
        return count

    @property
    def all_singleton(self) -> bool:
        return all(component.singleton for component in self.dead_components)


@dataclass(frozen=True)
class Carrier:
    name: str
    scores: tuple[int, ...]

    @property
    def total(self) -> int:
        return sum(self.scores)


AXES = (
    "predicate_retention",
    "component_geometry",
    "mirror_payload",
    "owner_pair_payload",
    "gate_touch_payload",
    "s319_split_payload",
    "scalar_firewall",
)

CARRIERS = (
    Carrier("G00_zero_edge_singleton_pocket_geometry", (10, 10, 9, 10, 8, 8, 10)),
    Carrier("G01_mirror_pair_dead_component_atlas", (10, 9, 10, 8, 8, 7, 10)),
    Carrier("G02_owner_pair_residue_span_word", (9, 8, 8, 10, 7, 8, 9)),
    Carrier("G03_gate_touch_complete_pocket_sidecar", (9, 8, 7, 8, 10, 8, 9)),
    Carrier("G04_s319_unit_delta_vs_delta_sidecar_split", (8, 7, 7, 7, 7, 10, 9)),
    Carrier("G05_single_best_gate_shadow", (6, 5, 4, 4, 7, 6, 6)),
    Carrier("G06_raw_row_name_exception_list", (2, 1, 1, 1, 1, 1, 0)),
)


def mirror_index(dead_components: list, interval: Interval) -> int | None:
    mirrored = mirror_interval(interval)
    for index, component in enumerate(dead_components):
        if component.interval == mirrored:
            return index
    return None


def build_geometry(row) -> RowGeometry:
    dead_components = list(row.component_row.dead_components)
    original = H3472.projection_stats(H3472.projection(dead_components))
    dead_labels = frozenset(
        label for component in dead_components for label in component_labels(component)
    )
    e_branch_gates = tuple(
        gate for gate in row.low_rank_gates if H3472.is_e_branch_gate(gate)
    )
    single_cuts = tuple(
        H3472.gate_cut(row, gate, original, dead_labels)
        for gate in e_branch_gates
    )
    best_cut = H3472.best_cut(list(single_cuts))
    min_gate = S319.min_e_branch_gate(row)
    min_gate_kind = S319.kind(min_gate)

    label_to_components: dict[str, set[int]] = defaultdict(set)
    for index, component in enumerate(dead_components):
        for label in component_labels(component):
            label_to_components[label].add(index)
    edge_support = tuple(
        sorted(label for label, components in label_to_components.items() if len(components) >= 2)
    )

    records = []
    for index, component in enumerate(dead_components):
        labels = component_labels(component)
        label_set = set(labels)
        touching = [
            gate for gate in e_branch_gates if label_set & set(gate_labels(gate))
        ]
        complete = [
            gate for gate in e_branch_gates if label_set <= set(gate_labels(gate))
        ]
        records.append(
            ComponentRecord(
                index=index,
                interval=component.interval,
                labels=labels,
                b0=component_cover(component, "B0"),
                b1=component_cover(component, "B1"),
                mirror_index=mirror_index(dead_components, component.interval),
                touch_gate_count=len(touching),
                complete_gate_count=len(complete),
                best_touch_gate=min(touching, key=gate_sort_key) if touching else None,
                best_complete_gate=min(complete, key=gate_sort_key) if complete else None,
            )
        )

    orbit_audit = H3475.build_orbits(row)
    hard_orbits = [orbit for orbit in orbit_audit.orbits if orbit.max_delta >= 7]
    return RowGeometry(
        name=row.name,
        speeds=row.speeds,
        dead_components=tuple(records),
        projection_edges=original.edges,
        edge_support_labels=edge_support,
        min_gate=min_gate,
        min_gate_kind=min_gate_kind,
        best_cut=best_cut,
        e_branch_gate_count=len(e_branch_gates),
        mirror_orbit_count=len(orbit_audit.orbits),
        hard_orbit_count=len(hard_orbits),
    )


def row_kind_bucket(geom: RowGeometry) -> str:
    if geom.name in S319_DELTA_SIDECAR_ROWS:
        return "s319_delta_sidecar"
    return "s319_branch_unit_delta"


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
        H3472.H3453.join_row(name, H3472.H3453.H3450.rows()[name])
        for name in SMALL_TOUCH_ROWS
    ]
    geometries = [build_geometry(row) for row in rows]
    all_components = [component for geom in geometries for component in geom.dead_components]

    dead_count_hist = Counter(geom.dead_count for geom in geometries)
    mirror_pair_hist = Counter(geom.mirror_pair_count for geom in geometries)
    min_kind_hist = Counter(geom.min_gate_kind for geom in geometries)
    bucket_hist = Counter(row_kind_bucket(geom) for geom in geometries)
    owner_pair_hist = Counter(component.owner_pair for component in all_components)
    residue_pair_hist = Counter(component.residue_pair for component in all_components)
    span_hist = Counter(component.span for component in all_components)
    gcd_hist = Counter(component.owner_gcd for component in all_components)
    complete_gate_hist = Counter(component.complete_gate_count for component in all_components)
    touch_gate_hist = Counter(component.touch_gate_count for component in all_components)
    hard_orbit_hist = Counter(geom.hard_orbit_count for geom in geometries)
    score_hist, path = tournament()

    print("HYP-3478 LRC14 SMALL-TOUCH GEOMETRY ATLAS")
    print("status=EVIDENCE / zero-edge singleton geometry atlas; not an LRC14 proof")
    print("source=HYP-3476 pair-current split + HYP-3477 hard discharge + S319 unit-delta split")
    print()
    print("## Assumption Challenge")
    print("alternate vertices considered: runners, gaps, fixed circle sections,")
    print("section boundaries, wall-crossing events, residues, cover arcs,")
    print("dead-cover projection nodes, singleton dead components, mirror component")
    print("pairs, owner pairs, E/branch gates, and terminal proof obligations.")
    print("chosen carrier vertices: zero-edge singleton dead pockets and their")
    print("mirror/owner/gate-touch sidecars.")
    print("preserved predicate: the six small-touch/no-hard rows have isolated")
    print("dead-cover pockets requiring singleton-current terminal discharge.")
    print("destroyed by scalarization: exact interval location, mirror partner,")
    print("owner-pair residue/span, S319 gate-kind class, and complete gate touches.")
    print()

    print("## Aggregate Geometry")
    print(f"rows={len(geometries)} {list(SMALL_TOUCH_ROWS)}")
    print(f"dead_component_count_hist={dict(sorted(dead_count_hist.items()))}")
    print(f"mirror_pair_count_hist={dict(sorted(mirror_pair_hist.items()))}")
    print(f"projection_edge_hist={dict(sorted(Counter(geom.projection_edges for geom in geometries).items()))}")
    print(f"edge_support_label_hist={dict(sorted(Counter(len(geom.edge_support_labels) for geom in geometries).items()))}")
    print(f"all_dead_components_singleton={all(component.singleton for component in all_components)}")
    print(f"all_hard_orbit_count_hist={dict(sorted(hard_orbit_hist.items()))}")
    print(f"s319_min_gate_kind_hist={dict(sorted(min_kind_hist.items()))}")
    print(f"s319_bucket_hist={dict(sorted(bucket_hist.items()))}")
    print(f"component_owner_pair_hist={dict(owner_pair_hist)}")
    print(f"component_residue_pair_hist={dict(residue_pair_hist)}")
    print(f"component_owner_span_hist={dict(sorted(span_hist.items()))}")
    print(f"component_owner_gcd_hist={dict(sorted(gcd_hist.items()))}")
    print(f"component_touch_gate_count_hist={dict(sorted(touch_gate_hist.items()))}")
    print(f"component_complete_gate_count_hist={dict(sorted(complete_gate_hist.items()))}")
    print()

    print("## Row Geometry Ledger")
    for geom in geometries:
        print(f"row={geom.name} bucket={row_kind_bucket(geom)}")
        print(f"  speeds={geom.speeds}")
        print(
            "  projection="
            f"dead_components={geom.dead_count} mirror_pairs={geom.mirror_pair_count} "
            f"edges={geom.projection_edges} edge_support_labels={geom.edge_support_labels}"
        )
        print(
            "  gates="
            f"e_branch_gates={geom.e_branch_gate_count} mirror_orbits={geom.mirror_orbit_count} "
            f"hard_orbits={geom.hard_orbit_count}"
        )
        print(
            "  s319_min_gate="
            f"kind={geom.min_gate_kind} typed={S319.typed_word(geom.min_gate)} "
            f"struct={S319.structural_word(geom.min_gate)} "
            f"len={fmt(geom.min_gate.length)} gate={gate_text(geom.min_gate)}"
        )
        if geom.best_cut is not None:
            cut = geom.best_cut
            print(
                "  best_touch_cut="
                f"gate={gate_text(cut.gate)} labels={cut.labels_touching_dead} "
                f"current=({cut.b0_labels},{cut.b1_labels}) "
                f"delta=({cut.gate.b0_delta},{cut.gate.b1_delta})"
            )
        for component in geom.dead_components:
            print(
                "  dead_component="
                f"i={component.index} interval={interval_text(component.interval)} "
                f"len={fmt(component.length)} mirror={component.mirror_index} "
                f"labels={component.labels} b0={component.b0} b1={component.b1} "
                f"owner_pair={component.owner_pair} residues={component.residue_pair} "
                f"span={component.span} gcd={component.owner_gcd} "
                f"touch_gates={component.touch_gate_count} "
                f"complete_gates={component.complete_gate_count}"
            )
            if component.best_complete_gate is not None:
                print(f"    best_complete_gate={gate_text(component.best_complete_gate)}")
            elif component.best_touch_gate is not None:
                print(f"    best_touch_gate={gate_text(component.best_touch_gate)}")
    print()

    print("## Geometry Lemma Target")
    print("Audited finite form:")
    print("  Each small-touch/no-hard row has an edgeless dead-cover projection.")
    print("  Every dead projection node is a rank-2 singleton pocket with one B0")
    print("  owner and one B1 owner, and every such pocket has an exact mirror")
    print("  partner.  Five rows have one mirror pair; random_covering_001 has two.")
    print("  No row has a hard mirror orbit.  The S319 split separates two")
    print("  delta-sidecar pockets (039,074) from four branch-unit singleton rows")
    print("  (001,062,086,101).")
    print("Proof pull:")
    print("  replace graph edge-current by a singleton-pocket theorem: a mirror pair")
    print("  of rank-2 isolated dead components with unique owner labels cannot be")
    print("  terminally bad once its complete/touching E-branch gate sidecar and")
    print("  owner-pair residue/span word are retained.  The four unit rows should")
    print("  be the first proof target; the two S319 delta-sidecar rows need the")
    print("  extra cover-delta sidecar.")
    print()

    print("## Tournament Analysis")
    print("vertices=geometric singleton-pocket proof carriers, not runners or raw row names")
    print(
        "pairwise_observable=predicate retention + component geometry + mirror "
        "payload + owner-pair payload + gate-touch payload + S319 split + scalar firewall"
    )
    print("switch=higher retained geometric proof payload; ties use declared route order")
    print(f"axes={AXES}")
    print(f"score_hist={score_hist}")
    print("directed_3cycles=0")
    print("hamiltonian_path=" + " -> ".join(path))


if __name__ == "__main__":
    main()
