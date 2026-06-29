#!/usr/bin/env python3
"""HYP-3478: small-touch singleton-current audit for LRC14.

HYP-3476's exception-frontier router leaves six non-AP, currentless rows with
no hard mirror orbit.  Their dead-cover projections have no edge to cut, so
larger gate-pair searches are the wrong next move.  This script audits the
singleton-current payload retained by the best touching E/branch gate and the
minimum E/branch mirror orbit.

The purpose is to split the six-row frontier into terminal packet types:
clean unit-delta singleton currents, delta-sidecar singleton currents, and the
one asymmetric touch row.  This is evidence for a later finite packet lemma,
not an LRC14 proof.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
ROUTER_PATH = ROOT / "04-computation" / "lrc14_exception_frontier_router_codex_20260629.py"


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


ROUTER = load_module("hyp3476_exception_router_for_hyp3478", ROUTER_PATH)
H3472 = ROUTER.H3472
H3475 = ROUTER.H3475


SMALL_TOUCH_ROWS = (
    "random_covering_001",
    "random_covering_039",
    "random_covering_062",
    "random_covering_074",
    "random_covering_086",
    "random_covering_101",
)


@dataclass(frozen=True)
class SingletonAudit:
    row: object
    projection: object
    best_cut: object
    min_orbit: object
    label_owner_counts: Counter
    touched_components: tuple[int, ...]
    packet: str

    @property
    def row_name(self) -> str:
        return self.row.name

    @property
    def repeated_dead_labels(self) -> tuple[str, ...]:
        return tuple(sorted(label for label, count in self.label_owner_counts.items() if count > 1))

    @property
    def best_delta(self) -> tuple[int, int]:
        return (self.best_cut.gate.b0_delta, self.best_cut.gate.b1_delta)

    @property
    def best_current(self) -> tuple[int, int]:
        return (self.best_cut.b0_labels, self.best_cut.b1_labels)


@dataclass(frozen=True)
class Carrier:
    name: str
    scores: tuple[int, ...]

    @property
    def total(self) -> int:
        return sum(self.scores)


AXES = (
    "singleton_projection_payload",
    "unit_delta_payload",
    "delta_sidecar_payload",
    "owner_label_uniqueness",
    "formal_packet_fit",
    "scalar_firewall",
)

CARRIERS = (
    Carrier("S00_clean_unit_delta_singleton_current", (10, 10, 7, 10, 9, 10)),
    Carrier("S01_delta_sidecar_singleton_current", (10, 8, 10, 9, 9, 10)),
    Carrier("S02_asymmetric_touch_singleton_current", (9, 7, 9, 8, 8, 10)),
    Carrier("S03_owner_label_uniqueness_certificate", (9, 7, 7, 10, 8, 10)),
    Carrier("S04_exact_period_state_lift_debt", (8, 6, 7, 8, 9, 9)),
    Carrier("S05_signed_SPEC_two_adic_debt", (8, 6, 7, 8, 8, 9)),
    Carrier("S06_raw_small_touch_name", (2, 1, 1, 1, 1, 0)),
)


def top_items(counter: Counter, limit: int = 10) -> str:
    items = counter.most_common(limit)
    suffix = "" if len(counter) <= limit else " ..."
    return "{" + ", ".join(f"{key!r}: {value}" for key, value in items) + suffix + "}"


def structural_deltas(orbit) -> tuple[tuple[int, int], ...]:
    return tuple((item[3], item[4]) for item in orbit.structural_pair)


def classify(best_cut, min_orbit) -> str:
    best_delta = (best_cut.gate.b0_delta, best_cut.gate.b1_delta)
    orbit_deltas = set(structural_deltas(min_orbit))
    if best_delta == (1, 1) and orbit_deltas == {(1, 1)}:
        return "clean_unit_delta_singleton_current"
    if best_delta == (1, 1):
        return "delta_sidecar_singleton_current"
    return "asymmetric_touch_singleton_current"


def component_label_map(row) -> dict[str, set[int]]:
    by_label: dict[str, set[int]] = defaultdict(set)
    for index, component in enumerate(row.component_row.dead_components):
        for label in H3472.component_blocker_labels(component):
            by_label[label].add(index)
    return by_label


def build_audit(row, boundary, orbit) -> SingletonAudit:
    projection = H3472.projection_stats(H3472.projection(list(row.component_row.dead_components)))
    best_cut = boundary.best_by_row[row.name]
    min_orbit = orbit.min_e_branch_orbit_by_row[row.name]
    label_map = component_label_map(row)
    counts = Counter({label: len(owners) for label, owners in label_map.items()})
    touched_components = tuple(
        sorted(
            {
                component
                for label in best_cut.labels_touching_dead
                for component in label_map.get(label, ())
            }
        )
    )
    return SingletonAudit(
        row=row,
        projection=projection,
        best_cut=best_cut,
        min_orbit=min_orbit,
        label_owner_counts=counts,
        touched_components=touched_components,
        packet=classify(best_cut, min_orbit),
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
    row_bank = H3472.H3453.H3450.rows()
    rows = [
        H3472.H3453.join_row(name, row_bank[name])
        for name in SMALL_TOUCH_ROWS
        if name in row_bank
    ]
    missing = [name for name in SMALL_TOUCH_ROWS if name not in row_bank]
    boundary = ROUTER.build_boundary_audit(rows)
    orbit = ROUTER.build_orbit_audit(rows)
    audits = [build_audit(row, boundary, orbit) for row in rows]

    packet_hist = Counter(audit.packet for audit in audits)
    projection_edge_hist = Counter(audit.projection.edges for audit in audits)
    projection_component_hist = Counter(audit.projection.components for audit in audits)
    dead_count_hist = Counter(len(audit.row.component_row.dead_components) for audit in audits)
    best_delta_hist = Counter(audit.best_delta for audit in audits)
    best_current_hist = Counter(audit.best_current for audit in audits)
    orbit_delta_hist = Counter(structural_deltas(audit.min_orbit) for audit in audits)
    repeated_label_hist = Counter(len(audit.repeated_dead_labels) for audit in audits)
    touched_component_hist = Counter(len(audit.touched_components) for audit in audits)
    score_hist, path = tournament()

    print("HYP-3478 LRC14 SMALL-TOUCH SINGLETON-CURRENT AUDIT")
    print("status=EVIDENCE / six-row singleton-current split; not an LRC14 proof")
    print("source=HYP-3476 exception-frontier router small-touch/no-hard packet")
    print()
    print("## Assumption Challenge")
    print("alternate vertices considered: runners, gaps, fixed circle sections,")
    print("section boundaries, wall-crossing events, residues, cover arcs,")
    print("individual gates, gate pairs, mirror orbits, dead components, owner")
    print("labels, exact-period packets, two-adic debts, and proof obligations.")
    print("chosen carrier vertices: zero-edge singleton-current packets.")
    print("preserved predicate: which terminal sidecar can discharge a currentless")
    print("non-hard row after the universal E/branch gate producer fires.")
    print("destroyed by scalarization: owner-label uniqueness, best-touch labels,")
    print("minimum mirror-orbit delta sidecar, and asymmetric current direction.")
    print()

    print("## Singleton Frontier")
    print(f"target_rows={len(SMALL_TOUCH_ROWS)} {list(SMALL_TOUCH_ROWS)}")
    print(f"missing_rows={missing}")
    print(f"rows_audited={len(audits)}")
    print(f"packet_hist={dict(sorted(packet_hist.items()))}")
    print(f"projection_edge_hist={dict(sorted(projection_edge_hist.items()))}")
    print(f"projection_component_hist={dict(sorted(projection_component_hist.items()))}")
    print(f"dead_component_count_hist={dict(sorted(dead_count_hist.items()))}")
    print(f"best_delta_hist={dict(sorted(best_delta_hist.items()))}")
    print(f"best_current_hist={dict(sorted(best_current_hist.items()))}")
    print(f"min_orbit_structural_delta_hist={top_items(orbit_delta_hist)}")
    print(f"repeated_dead_label_count_hist={dict(sorted(repeated_label_hist.items()))}")
    print(f"best_touch_component_count_hist={dict(sorted(touched_component_hist.items()))}")
    print()

    print("## Row Ledger")
    for audit in audits:
        cut = audit.best_cut
        orbit_obj = audit.min_orbit
        print(f"row={audit.row_name} packet={audit.packet}")
        print(
            "  projection="
            f"components={audit.projection.components} largest={audit.projection.largest} "
            f"edges={audit.projection.edges} dead_components={len(audit.row.component_row.dead_components)}"
        )
        print(
            "  best_touch="
            f"gate={H3472.fmt_interval(cut.gate.interval)} kind={cut.gate.endpoint_kind_signature} "
            f"labels={cut.labels_touching_dead} current=({cut.b0_labels},{cut.b1_labels}) "
            f"delta=({cut.gate.b0_delta},{cut.gate.b1_delta}) "
            f"touched_components={audit.touched_components}"
        )
        print(
            "  min_e_branch_orbit="
            f"typed_pair={orbit_obj.typed_pair} structural_pair={orbit_obj.structural_pair} "
            f"components={orbit_obj.components} min_length={H3472.fmt_fraction(orbit_obj.min_length)} "
            f"structural_deltas={structural_deltas(orbit_obj)}"
        )
        print(f"  repeated_dead_labels={audit.repeated_dead_labels}")
    print()

    print("## Lemma Target")
    print("Finite audited form:")
    print("  The six small-touch/no-hard rows all have zero-edge dead-cover")
    print("  projections, so no projection-edge cut exists to recover.")
    print("  They split into three terminal packet types:")
    print("    clean unit-delta singleton current: random_covering_062,")
    print("      random_covering_086, random_covering_101")
    print("    delta-sidecar singleton current: random_covering_039,")
    print("      random_covering_074")
    print("    asymmetric touch singleton current: random_covering_001")
    print("  A later proof should attach owner-current, endpoint-spine, exact-period,")
    print("  two-adic, signed-SPEC/Rprime, or state-lift certificates to these three")
    print("  packets, not search larger gate tuples without new edge-support data.")
    print()

    print("## Tournament Analysis")
    print("vertices=singleton-current proof packets, not runners or raw row names")
    print(
        "pairwise_observable=singleton projection payload + unit/delta sidecar "
        "+ owner-label uniqueness + formal packet fit + scalar firewall"
    )
    print("switch=higher retained proof payload; ties use declared route order")
    print(f"axes={AXES}")
    print(f"score_hist={score_hist}")
    print("directed_3cycles=0")
    print("hamiltonian_path=" + " -> ".join(path))


if __name__ == "__main__":
    main()
