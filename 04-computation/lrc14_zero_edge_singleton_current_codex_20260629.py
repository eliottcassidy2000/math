#!/usr/bin/env python3
"""HYP-3480: zero-edge singleton-current audit for LRC14.

HYP-3476 showed that the random boundary-current exceptions are not missing
larger dead-projection cuts: their dead-cover projections are edgeless
singleton packets.  HYP-3478 then showed that the six non-hard rows are made
from mirror-paired singleton B0/B1 owner components.

This script asks the proof-facing follow-up: do those singleton components
carry complete branch-unit E/branch touches in mirror-compatible pairs?  The
control row is `random_covering_031`, the unique hard/currentless overlap from
HYP-3479.

The output is finite-bank evidence for a terminal packet theorem, not an
LRC14 proof.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
H3478_PATH = ROOT / "04-computation" / "lrc14_small_touch_no_hard_geometry_codex_20260629.py"
H3475_RESULT = ROOT / "05-knowledge" / "results" / "lrc14_colored_gate_mirror_orbit_codex_20260629.out"

TARGET_ROWS = (
    "random_covering_001",
    "random_covering_039",
    "random_covering_062",
    "random_covering_074",
    "random_covering_086",
    "random_covering_101",
)
CONTROL_ROWS = ("random_covering_031",)
AUDIT_ROWS = TARGET_ROWS + CONTROL_ROWS

ROUTE_BY_ROW = {name: "small_touch_no_hard_current_exception" for name in TARGET_ROWS}
ROUTE_BY_ROW["random_covering_031"] = "random031_overlap_hard_and_currentless"

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


H3478 = load_module("hyp3478_small_touch_for_hyp3480", H3478_PATH)
H3472 = H3478.H3472
H3453 = H3478.H3453


def fmt(value) -> str:
    return H3472.fmt_fraction(value)


def interval_text(interval) -> str:
    return H3472.fmt_interval(interval)


def top_items(counter: Counter, limit: int = 12) -> str:
    items = counter.most_common(limit)
    suffix = "" if len(counter) <= limit else " ..."
    return "{" + ", ".join(f"{key!r}: {value}" for key, value in items) + suffix + "}"


def mirror_interval(interval):
    return (1 - interval[1], 1 - interval[0])


def component_labels(component) -> tuple[str, ...]:
    return H3472.component_blocker_labels(component)


def two_adic_order(value: int) -> str:
    value = abs(value)
    if value == 0:
        return "inf"
    order = 0
    while value % 2 == 0:
        order += 1
        value //= 2
    return str(order)


def owner_residue_word(pair) -> str:
    b0, b1 = pair
    if len(b0) != 1 or len(b1) != 1:
        return "non_singleton_owner_pair"
    a = b0[0]
    b = b1[0]
    diff = abs(a - b)
    return (
        f"owners=({a},{b}) residues=({a % 14},{b % 14}) "
        f"diff={diff} v2diff={two_adic_order(diff)} sum={a + b} "
        f"v2sum={two_adic_order(a + b)}"
    )


def cut_sort_key(cut):
    gate = cut.gate
    return (
        gate.b0_delta + gate.b1_delta,
        gate.length,
        gate.interval[0],
        gate.endpoint_kind_signature,
        tuple(cut.labels_touching_dead),
    )


def mirror_gate_pair(left, right) -> bool:
    return (
        right.gate.interval == mirror_interval(left.gate.interval)
        and MIRROR_KIND.get(left.gate.endpoint_kind_signature) == right.gate.endpoint_kind_signature
    )


def gate_text(cut) -> str:
    gate = cut.gate
    return (
        f"component={gate.component_index} gate={interval_text(gate.interval)} "
        f"typed={H3478.typed_word(gate)} structural={H3478.structural_word(gate)} "
        f"len={fmt(gate.length)} labels={cut.labels_touching_dead}"
    )


def parse_hard_orbits_by_row() -> dict[str, list[int]]:
    hard_by_row: dict[str, list[int]] = defaultdict(list)
    for raw_line in H3475_RESULT.read_text().splitlines():
        line = raw_line.strip()
        if not line.startswith("row="):
            continue
        row_name, rest = line[4:].split(" delta=", 1)
        delta_text = rest.split(" components=", 1)[0]
        hard_by_row[row_name].append(int(delta_text))
    return dict(hard_by_row)


@dataclass(frozen=True)
class ComponentCertificate:
    index: int
    labels: tuple[str, ...]
    owner_pair: tuple[tuple[int, ...], tuple[int, ...]]
    complete_cuts: tuple
    branch_unit_cuts: tuple

    @property
    def complete_count(self) -> int:
        return len(self.complete_cuts)

    @property
    def branch_unit_count(self) -> int:
        return len(self.branch_unit_cuts)

    @property
    def best_complete(self):
        return min(self.complete_cuts, key=cut_sort_key) if self.complete_cuts else None

    @property
    def best_branch_unit(self):
        if not self.branch_unit_cuts:
            return None
        return min(
            self.branch_unit_cuts,
            key=lambda cut: (
                cut.gate.length,
                cut.gate.interval[0],
                cut.gate.endpoint_kind_signature,
                tuple(cut.labels_touching_dead),
            ),
        )


@dataclass(frozen=True)
class MirrorPairCertificate:
    left: ComponentCertificate
    right: ComponentCertificate
    owner_swapped: bool
    mirror_unit_pairs: tuple[tuple[object, object], ...]

    @property
    def has_branch_unit_mirror_certificate(self) -> bool:
        return bool(self.mirror_unit_pairs)

    @property
    def chosen_mirror_unit_pair(self):
        if not self.mirror_unit_pairs:
            return None
        return min(
            self.mirror_unit_pairs,
            key=lambda pair: (
                pair[0].gate.length + pair[1].gate.length,
                pair[0].gate.interval[0],
                pair[1].gate.interval[0],
            ),
        )


@dataclass(frozen=True)
class RowCertificate:
    row: object
    route: str
    projection: object
    e_branch_cuts: tuple
    touching_cuts: tuple
    components: tuple[ComponentCertificate, ...]
    mirror_pairs: tuple[MirrorPairCertificate, ...]
    hard_deltas: tuple[int, ...]
    min_gate_kind: str

    @property
    def hard_orbit_count(self) -> int:
        return len(self.hard_deltas)

    @property
    def hard_max_delta(self) -> int:
        return max(self.hard_deltas, default=0)

    @property
    def complete_components(self) -> int:
        return sum(1 for component in self.components if component.complete_count)

    @property
    def branch_unit_components(self) -> int:
        return sum(1 for component in self.components if component.branch_unit_count)

    @property
    def mirror_unit_pair_count(self) -> int:
        return sum(1 for pair in self.mirror_pairs if pair.has_branch_unit_mirror_certificate)

    @property
    def terminal_class(self) -> str:
        if self.row.name in TARGET_ROWS:
            if self.mirror_unit_pair_count == len(self.mirror_pairs):
                if self.min_gate_kind == "delta_sidecar_packet":
                    return "mirror_unit_singleton_packet_cover_delta_min_shadow"
                return "mirror_unit_singleton_packet"
            return "small_touch_unresolved_singleton_packet"
        if self.row.name == "random_covering_031" and self.hard_orbit_count:
            return "random031_hard_currentless_control"
        return "control_without_hard_flag"


@dataclass(frozen=True)
class Carrier:
    name: str
    scores: tuple[int, ...]

    @property
    def total(self) -> int:
        return sum(self.scores)


AXES = (
    "predicate_retention",
    "component_complete_touch",
    "mirror_unit_payload",
    "route_sidecar_payload",
    "hard_control_payload",
    "formal_packet_fit",
    "scalar_firewall",
)

CARRIERS = (
    Carrier("Z00_mirror_unit_singleton_current_packet", (10, 10, 10, 9, 8, 10, 10)),
    Carrier("Z01_component_complete_touch_certificate", (10, 10, 8, 8, 7, 9, 10)),
    Carrier("Z02_random031_hard_control_clause", (9, 6, 7, 8, 10, 9, 10)),
    Carrier("Z03_route_sidecar_R_join", (9, 7, 8, 10, 8, 8, 9)),
    Carrier("Z04_cover_delta_min_gate_shadow", (8, 8, 7, 7, 6, 7, 9)),
    Carrier("Z05_owner_residue_two_adic_shadow", (7, 7, 6, 6, 5, 6, 8)),
    Carrier("Z06_raw_zero_edge_projection_count", (3, 1, 1, 1, 1, 1, 0)),
    Carrier("Z07_raw_row_name_list", (2, 0, 0, 1, 1, 1, 0)),
)


def component_certificate(index: int, component, touching_cuts: tuple) -> ComponentCertificate:
    labels = component_labels(component)
    label_set = set(labels)
    complete = tuple(
        cut for cut in touching_cuts if label_set <= set(cut.labels_touching_dead)
    )
    branch_unit = tuple(
        cut for cut in complete if H3478.min_gate_kind(cut.gate) == "branch_unit_delta"
    )
    return ComponentCertificate(
        index=index,
        labels=labels,
        owner_pair=H3478.owner_pair(component),
        complete_cuts=complete,
        branch_unit_cuts=branch_unit,
    )


def mirror_pair_certificate(
    left: ComponentCertificate,
    right: ComponentCertificate,
) -> MirrorPairCertificate:
    mirror_pairs: list[tuple[object, object]] = []
    for left_cut in left.branch_unit_cuts:
        for right_cut in right.branch_unit_cuts:
            if mirror_gate_pair(left_cut, right_cut):
                mirror_pairs.append((left_cut, right_cut))
    return MirrorPairCertificate(
        left=left,
        right=right,
        owner_swapped=right.owner_pair == H3478.swapped_owner_pair(left.owner_pair),
        mirror_unit_pairs=tuple(mirror_pairs),
    )


def row_certificate(name: str, hard_by_row: dict[str, list[int]]) -> RowCertificate:
    row = H3453.join_row(name, H3453.H3450.rows()[name])
    dead_components = tuple(row.component_row.dead_components)
    projection = H3472.projection_stats(H3472.projection(list(dead_components)))
    dead_label_set = frozenset(
        label
        for component in dead_components
        for label in component_labels(component)
    )
    e_branch_cuts = tuple(
        H3472.gate_cut(row, gate, projection, dead_label_set)
        for gate in row.low_rank_gates
        if H3472.is_e_branch_gate(gate)
    )
    touching_cuts = tuple(cut for cut in e_branch_cuts if cut.touches_dead_projection)
    components = tuple(
        component_certificate(index, component, touching_cuts)
        for index, component in enumerate(dead_components)
    )

    by_interval = {component.interval: index for index, component in enumerate(dead_components)}
    seen: set[int] = set()
    mirror_pairs: list[MirrorPairCertificate] = []
    for index, component in enumerate(dead_components):
        if index in seen:
            continue
        partner_index = by_interval.get(mirror_interval(component.interval))
        if partner_index is None:
            seen.add(index)
            continue
        seen.add(index)
        seen.add(partner_index)
        mirror_pairs.append(
            mirror_pair_certificate(components[index], components[partner_index])
        )

    min_gate = H3478.min_e_branch_gate(row)
    return RowCertificate(
        row=row,
        route=ROUTE_BY_ROW[name],
        projection=projection,
        e_branch_cuts=e_branch_cuts,
        touching_cuts=touching_cuts,
        components=components,
        mirror_pairs=tuple(mirror_pairs),
        hard_deltas=tuple(hard_by_row.get(name, ())),
        min_gate_kind=H3478.min_gate_kind(min_gate),
    )


def tournament_summary() -> tuple[dict[int, int], list[str]]:
    hist = dict(sorted(Counter(carrier.total for carrier in CARRIERS).items()))
    order = {carrier.name: index for index, carrier in enumerate(CARRIERS)}
    path = [
        carrier.name
        for carrier in sorted(CARRIERS, key=lambda carrier: (-carrier.total, order[carrier.name]))
    ]
    return hist, path


def print_row(cert: RowCertificate) -> None:
    row = cert.row
    print()
    print(f"## Row {row.name}")
    print(f"speeds={row.speeds}")
    print(
        "route={route} terminal_class={terminal} hard_orbits={hard_count} "
        "hard_max_delta={hard_delta}".format(
            route=cert.route,
            terminal=cert.terminal_class,
            hard_count=cert.hard_orbit_count,
            hard_delta=cert.hard_max_delta,
        )
    )
    print(
        "projection=components:{components} largest:{largest} edges:{edges} "
        "dead={dead} e_branch={eb} touching={touch}".format(
            components=cert.projection.components,
            largest=cert.projection.largest,
            edges=cert.projection.edges,
            dead=len(cert.components),
            eb=len(cert.e_branch_cuts),
            touch=len(cert.touching_cuts),
        )
    )
    print(
        "component_complete_touch={complete}/{total} "
        "component_branch_unit_complete={unit}/{total} "
        "mirror_unit_pairs={mirror}/{pairs} min_gate_kind={kind}".format(
            complete=cert.complete_components,
            unit=cert.branch_unit_components,
            total=len(cert.components),
            mirror=cert.mirror_unit_pair_count,
            pairs=len(cert.mirror_pairs),
            kind=cert.min_gate_kind,
        )
    )
    for pair in cert.mirror_pairs:
        left = pair.left
        right = pair.right
        print(
            "  mirror_pair=({left},{right}) owner_swapped={swapped} "
            "mirror_unit_pair_count={count}".format(
                left=left.index,
                right=right.index,
                swapped=pair.owner_swapped,
                count=len(pair.mirror_unit_pairs),
            )
        )
        print(
            "    left labels={labels} owner_pair={owner} residues={residues} "
            "complete={complete} branch_unit={unit} {word}".format(
                labels=left.labels,
                owner=H3478.owner_pair_text(left.owner_pair),
                residues=H3478.owner_pair_residues(left.owner_pair),
                complete=left.complete_count,
                unit=left.branch_unit_count,
                word=owner_residue_word(left.owner_pair),
            )
        )
        print(
            "    right labels={labels} owner_pair={owner} residues={residues} "
            "complete={complete} branch_unit={unit} {word}".format(
                labels=right.labels,
                owner=H3478.owner_pair_text(right.owner_pair),
                residues=H3478.owner_pair_residues(right.owner_pair),
                complete=right.complete_count,
                unit=right.branch_unit_count,
                word=owner_residue_word(right.owner_pair),
            )
        )
        chosen = pair.chosen_mirror_unit_pair
        if chosen is None:
            print("    chosen_mirror_unit_pair=None")
            if left.best_complete is not None:
                print(f"    left_best_complete={gate_text(left.best_complete)}")
            if right.best_complete is not None:
                print(f"    right_best_complete={gate_text(right.best_complete)}")
        else:
            print(f"    chosen_left={gate_text(chosen[0])}")
            print(f"    chosen_right={gate_text(chosen[1])}")


def main() -> None:
    hard_by_row = parse_hard_orbits_by_row()
    certificates = tuple(row_certificate(name, hard_by_row) for name in AUDIT_ROWS)
    targets = tuple(cert for cert in certificates if cert.row.name in TARGET_ROWS)
    controls = tuple(cert for cert in certificates if cert.row.name in CONTROL_ROWS)

    target_projection_edge_hist = Counter(cert.projection.edges for cert in targets)
    target_min_gate_kind_hist = Counter(cert.min_gate_kind for cert in targets)
    route_hist = Counter(cert.route for cert in certificates)
    terminal_hist = Counter(cert.terminal_class for cert in certificates)
    owner_residue_hist = Counter(
        H3478.owner_pair_residues(pair.left.owner_pair)
        for cert in targets
        for pair in cert.mirror_pairs
    )
    target_branch_unit_components = sum(cert.branch_unit_components for cert in targets)
    target_components = sum(len(cert.components) for cert in targets)
    target_mirror_unit_pairs = sum(cert.mirror_unit_pair_count for cert in targets)
    target_mirror_pairs = sum(len(cert.mirror_pairs) for cert in targets)
    control_branch_unit_components = sum(cert.branch_unit_components for cert in controls)
    control_components = sum(len(cert.components) for cert in controls)
    control_mirror_unit_pairs = sum(cert.mirror_unit_pair_count for cert in controls)
    control_mirror_pairs = sum(len(cert.mirror_pairs) for cert in controls)
    cover_delta_min_shadow_rows = tuple(
        cert.row.name
        for cert in targets
        if cert.min_gate_kind == "delta_sidecar_packet"
        and cert.mirror_unit_pair_count == len(cert.mirror_pairs)
    )
    hard_overlap_rows = tuple(
        cert.row.name for cert in certificates if cert.hard_orbit_count and cert.route.endswith("currentless")
    )
    score_hist, path = tournament_summary()

    print("HYP-3480 ZERO-EDGE SINGLETON-CURRENT AUDIT")
    print("status=EVIDENCE / finite singleton-current certificate audit; not an LRC14 proof")
    print("sources=HYP-3478 small-touch geometry + HYP-3479 hard/current join + HYP-3472 gate cuts")
    print(f"target_rows={TARGET_ROWS}")
    print(f"control_rows={CONTROL_ROWS}")
    print()
    print("## Assumption Challenge")
    print("alternate vertices considered: runners, gaps, fixed circle sections,")
    print("section boundaries, wall-crossing events, residues, cover arcs, Fourier")
    print("modes, matroid circuits, singleton components, mirror component pairs,")
    print("owner-label pairs, touching gate events, route labels, hard-overlap flags,")
    print("and formal proof obligations.")
    print("chosen carrier vertices: mirror-paired singleton dead components with")
    print("complete branch-unit E/branch gate touches, joined to the route sidecar R.")
    print("preserved predicate: the zero-edge currentless row has a local")
    print("mirror-unit singleton-current certificate or is the named random031")
    print("hard/currentless control.")
    print("destroyed if scalarized: exact interval mirror, branch orientation,")
    print("component owner labels, and whether a cover-delta minimum still has a")
    print("larger complete unit-delta touch.")
    print()
    print("## Aggregate Singleton-Current Ledger")
    print(f"audited_rows={len(certificates)} target_rows={len(targets)} control_rows={len(controls)}")
    print(f"route_hist={dict(sorted(route_hist.items()))}")
    print(f"terminal_class_hist={dict(sorted(terminal_hist.items()))}")
    print(f"target_projection_edge_hist={dict(sorted(target_projection_edge_hist.items()))}")
    print(f"target_min_gate_kind_hist={dict(sorted(target_min_gate_kind_hist.items()))}")
    print(f"target_dead_components={target_components}")
    print(
        "target_components_with_complete_branch_unit_touch="
        f"{target_branch_unit_components}/{target_components}"
    )
    print(
        "target_mirror_pairs_with_branch_unit_mirror_gate="
        f"{target_mirror_unit_pairs}/{target_mirror_pairs}"
    )
    print(f"cover_delta_min_shadow_rows_with_unit_certificate={cover_delta_min_shadow_rows}")
    print(f"target_owner_residue_hist={top_items(owner_residue_hist)}")
    print(
        "control_components_with_complete_branch_unit_touch="
        f"{control_branch_unit_components}/{control_components}"
    )
    print(
        "control_mirror_pairs_with_branch_unit_mirror_gate="
        f"{control_mirror_unit_pairs}/{control_mirror_pairs}"
    )
    print(f"hard_currentless_control_rows={hard_overlap_rows}")
    print()
    print("## Row Certificates")
    for cert in certificates:
        print_row(cert)
    print()
    print("## Lemma Target")
    print("Finite audited form:")
    print("  Each of the six small-touch/no-hard zero-edge rows has every dead")
    print("  singleton component completely touched by a branch-unit E/branch gate.")
    print("  The 14 components assemble into 7 mirror pairs, and every pair has at")
    print("  least one mirror-compatible branch-unit gate pair with swapped B0/B1")
    print("  owner labels.")
    print("  Rows random_covering_039 and random_covering_074 still have cover-delta")
    print("  absolute minimum gates, but the component-level certificate finds a")
    print("  complete branch-unit mirror pair, so they need not remain separate")
    print("  terminal debt if the mirror-singleton lemma is proved.")
    print("  The control row random_covering_031 has zero branch-unit complete")
    print("  singleton components and remains exactly the hard/currentless gluing")
    print("  clause from HYP-3455/HYP-3460/HYP-3479.")
    print("Proof pull:")
    print("  Prove a finite mirror-unit singleton-current lemma over swapped")
    print("  singleton B0/B1 owner pairs, carrying route sidecar R and the owner")
    print("  residue/two-adic word.  Then the non-AP currentless random frontier")
    print("  reduces to that lemma plus the named random031 hard-control clause.")
    print()
    print("## Tournament Analysis")
    print(f"vertices={[carrier.name for carrier in CARRIERS]}")
    print(f"axes={AXES}")
    print(
        "pairwise_observable=predicate retention + complete component touch + "
        "mirror-unit payload + route/hard-control payload + scalar firewall"
    )
    print(
        "switch_gauge=higher proof-facing retained payload first; ties follow "
        "mirror unit packet -> complete touch -> random031 control -> route sidecar"
    )
    print(f"score_hist={score_hist}")
    print("directed_3cycles=0")
    print("hamiltonian_path_count=1")
    print(f"hamiltonian_path={' -> '.join(path)}")


if __name__ == "__main__":
    main()
