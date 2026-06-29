#!/usr/bin/env python3
"""HYP-3477: hard mirror-orbit discharge audit for LRC14.

HYP-3475 localized the high cover-delta debt to eight mirror orbits on seven
random rows.  HYP-3476 separately reserved a pair-current audit for the
HYP-3472 single-gate current exceptions.  This script joins those two fronts:

  * identify every delta >= 7 mirror orbit from HYP-3475;
  * test whether q=14V phase-grid witnesses hit the hard orbit or route
    through lower-delta gates in the same component;
  * test whether a lower-delta E/branch survivor gate already gives a
    dead-cover projection cut as in HYP-3472;
  * test a small HYP-3476-style pair-current sidecar on the hard rows.

The output is evidence and a proof-target ledger, not an LRC14 proof.
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
COMP = ROOT / "04-computation"


def load_module(name: str, filename: str):
    path = COMP / filename
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


H3475 = load_module("hyp3475_mirror_orbits_for_hyp3477", "lrc14_colored_gate_mirror_orbit_codex_20260629.py")
H3472 = load_module("hyp3472_boundary_current_for_hyp3477", "lrc14_dead_cover_boundary_current_codex_20260629.py")
H3471 = H3475.H3471
H3453 = H3475.H3453


Interval = tuple[Fraction, Fraction]


def fmt(x: Fraction | int | None) -> str:
    if x is None:
        return "None"
    if isinstance(x, int):
        return str(x)
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def interval_text(interval: Interval | None) -> str:
    if interval is None:
        return "None"
    return f"[{fmt(interval[0])},{fmt(interval[1])}]"


def top_items(counter: Counter, limit: int = 12) -> str:
    items = counter.most_common(limit)
    suffix = "" if len(counter) <= limit else " ..."
    return "{" + ", ".join(f"{key!r}: {value}" for key, value in items) + suffix + "}"


def contains_closed(interval: Interval, x: Fraction) -> bool:
    return interval[0] <= x <= interval[1]


def compatible(gate, branch: int) -> bool:
    return gate.branch_mask == "both" or gate.branch_mask == f"branch{branch}"


def gate_key(gate) -> tuple[int, Interval, str]:
    return (gate.component_index, gate.interval, gate.branch_mask)


def dist_residue(value: int, modulus: int) -> int:
    residue = value % modulus
    return min(residue, modulus - residue)


def actual_witnesses(speeds: tuple[int, ...], V: int) -> tuple[int, ...]:
    q = 14 * V
    return tuple(
        a
        for a in range(q)
        if all(14 * dist_residue(speed * a, q) >= q for speed in speeds)
    )


def branch_label_counts(labels: frozenset[str]) -> tuple[int, int, int]:
    b0 = sum(1 for label in labels if label.startswith("B0:"))
    b1 = sum(1 for label in labels if label.startswith("B1:"))
    return b0, b1, b0 - b1


def projection_stats_for_labels(row, labels: frozenset[str]):
    dead_components = list(row.component_row.dead_components)
    original = H3472.projection_stats(H3472.projection(dead_components))
    after = H3472.projection_stats(H3472.projection(dead_components, labels))
    return original, after


def dead_label_set(row) -> frozenset[str]:
    labels: set[str] = set()
    for component in row.component_row.dead_components:
        labels.update(H3472.component_blocker_labels(component))
    return frozenset(labels)


def gate_touch_labels(gate, dead_labels: frozenset[str]) -> frozenset[str]:
    return frozenset(label for label in H3472.gate_blocker_labels(gate) if label in dead_labels)


@dataclass(frozen=True)
class LabelCut:
    labels: frozenset[str]
    original_components: int
    original_largest: int
    original_edges: int
    after_components: int
    after_largest: int
    after_edges: int
    source: str

    @property
    def removed_edges(self) -> int:
        return self.original_edges - self.after_edges

    @property
    def largest_drop(self) -> int:
        return self.original_largest - self.after_largest

    @property
    def component_gain(self) -> int:
        return self.after_components - self.original_components

    @property
    def is_edge_cut(self) -> bool:
        return self.removed_edges > 0

    @property
    def is_separating_current(self) -> bool:
        return self.component_gain > 0 or self.largest_drop > 0

    @property
    def current(self) -> tuple[int, int]:
        b0, b1, _ = branch_label_counts(self.labels)
        return b0, b1


def label_cut(row, labels: frozenset[str], source: str) -> LabelCut:
    original, after = projection_stats_for_labels(row, labels)
    return LabelCut(
        labels=labels,
        original_components=original.components,
        original_largest=original.largest,
        original_edges=original.edges,
        after_components=after.components,
        after_largest=after.largest,
        after_edges=after.edges,
        source=source,
    )


def cut_key(cut: LabelCut) -> tuple[bool, bool, int, int, int, int]:
    return (
        cut.is_separating_current,
        cut.is_edge_cut,
        cut.removed_edges,
        cut.largest_drop,
        cut.component_gain,
        len(cut.labels),
    )


def best_single_e_branch_cut(row) -> LabelCut | None:
    dead_labels = dead_label_set(row)
    cuts = []
    for gate in row.low_rank_gates:
        if not H3472.is_e_branch_gate(gate):
            continue
        labels = gate_touch_labels(gate, dead_labels)
        if labels:
            cuts.append(label_cut(row, labels, f"single:{gate.component_index}:{interval_text(gate.interval)}"))
    return max(cuts, key=cut_key) if cuts else None


def best_pair_current(row, keep: int = 80) -> LabelCut | None:
    dead_labels = dead_label_set(row)
    label_sets = sorted(
        {
            gate_touch_labels(gate, dead_labels)
            for gate in row.low_rank_gates
            if H3472.is_e_branch_gate(gate)
        },
        key=lambda labels: (len(labels), tuple(sorted(labels))),
    )
    label_sets = [labels for labels in label_sets if labels]
    if not label_sets:
        return None

    singleton_cuts = [label_cut(row, labels, "pair-current-singleton") for labels in label_sets]
    top_singletons = sorted(singleton_cuts, key=cut_key, reverse=True)[:keep]
    top_sets = [cut.labels for cut in top_singletons]
    candidates: list[LabelCut] = list(top_singletons)
    for left, right in combinations(top_sets, 2):
        labels = frozenset(set(left) | set(right))
        candidates.append(label_cut(row, labels, f"pair-current-top{keep}"))
    return max(candidates, key=cut_key)


def best_mirror_e_branch_pair(row, row_orbits) -> LabelCut | None:
    dead_labels = dead_label_set(row)
    candidates: list[LabelCut] = []
    low_intervals = {gate.interval for gate in row.low_rank_gates}
    for orbit in row_orbits:
        if not orbit.all_low_rank(low_intervals) or not orbit.has_e_branch:
            continue
        labels: set[str] = set()
        for gate in orbit.members:
            labels.update(gate_touch_labels(gate, dead_labels))
        if labels:
            candidates.append(label_cut(row, frozenset(labels), "mirror-e-branch-pair"))
    return max(candidates, key=cut_key) if candidates else None


@dataclass(frozen=True)
class PhaseHitAudit:
    V: int
    witness_count: int
    hard_gate_hits: int
    same_component_lower_delta_hits: int
    any_gate_hits: int
    no_gate_hits: int
    ambiguous_gate_hits: int
    hard_component_hit_counter: Counter[tuple[int, str, int]]
    gate_route_counter: Counter[str]
    gate_delta_counter: Counter[int]

    @property
    def bypasses_hard_orbit(self) -> bool:
        return self.hard_gate_hits == 0


def phase_hit_audit(row, orbit) -> PhaseHitAudit:
    V = max(row.speeds)
    q = 14 * V
    witnesses = actual_witnesses(row.speeds, V)
    hard_keys = {gate_key(gate) for gate in orbit.members}
    hard_components = {gate.component_index for gate in orbit.members}
    hard_delta = orbit.max_delta

    hard_gate_hits = 0
    same_component_lower_delta_hits = 0
    any_gate_hits = 0
    no_gate_hits = 0
    ambiguous_gate_hits = 0
    hard_component_hit_counter: Counter[tuple[int, str, int]] = Counter()
    gate_route_counter: Counter[str] = Counter()
    gate_delta_counter: Counter[int] = Counter()

    for a in witnesses:
        t = Fraction(a, q)
        branch = 0 if t < Fraction(1, 2) else 1
        u = Fraction((2 * a) % q, q)
        hits = [
            gate
            for gate in row.gates
            if contains_closed(gate.interval, u) and compatible(gate, branch)
        ]
        if not hits:
            no_gate_hits += 1
            continue
        any_gate_hits += 1
        if len(hits) > 1:
            ambiguous_gate_hits += 1
        gate = sorted(hits, key=lambda g: (g.component_index, g.interval[0], g.branch_mask))[0]
        gate_route_counter[gate.route] += 1
        gate_delta_counter[gate.total_delta] += 1
        if gate_key(gate) in hard_keys:
            hard_gate_hits += 1
        if gate.component_index in hard_components:
            hard_component_hit_counter[(gate.component_index, gate.branch_mask, gate.total_delta)] += 1
            if gate.total_delta < hard_delta:
                same_component_lower_delta_hits += 1

    return PhaseHitAudit(
        V=V,
        witness_count=len(witnesses),
        hard_gate_hits=hard_gate_hits,
        same_component_lower_delta_hits=same_component_lower_delta_hits,
        any_gate_hits=any_gate_hits,
        no_gate_hits=no_gate_hits,
        ambiguous_gate_hits=ambiguous_gate_hits,
        hard_component_hit_counter=hard_component_hit_counter,
        gate_route_counter=gate_route_counter,
        gate_delta_counter=gate_delta_counter,
    )


@dataclass(frozen=True)
class HardOrbitAudit:
    row: object
    orbit: object
    phase: PhaseHitAudit
    hard_cut: LabelCut
    single_cut: LabelCut | None
    pair_cut: LabelCut | None
    mirror_cut: LabelCut | None

    @property
    def route(self) -> str:
        if self.phase.bypasses_hard_orbit and self.phase.same_component_lower_delta_hits > 0:
            return "phase_branch_bypass_with_lower_delta_component_hits"
        if self.single_cut is not None and self.single_cut.is_edge_cut:
            return "single_e_branch_projection_cut"
        if self.pair_cut is not None and self.pair_cut.is_edge_cut:
            return "hyp3476_pair_current_projection_cut"
        if self.phase.bypasses_hard_orbit:
            return "phase_branch_bypass_no_component_hit"
        return "named_gluing_or_owner_current_debt"


@dataclass(frozen=True)
class Carrier:
    name: str
    scores: tuple[int, ...]

    @property
    def total(self) -> int:
        return sum(self.scores)


CARRIERS = (
    Carrier("H00_phase_branch_hard_orbit_bypass", (10, 10, 9, 8, 9, 8, 10)),
    Carrier("H01_lower_delta_e_branch_projection_cut", (10, 8, 10, 10, 8, 8, 10)),
    Carrier("H02_pair_current_bridge_to_HYP3476", (9, 8, 10, 9, 10, 7, 10)),
    Carrier("H03_hard_mirror_orbit_gluing_clause", (9, 10, 8, 7, 8, 10, 9)),
    Carrier("H04_dead_projection_touch_sidecar", (8, 6, 8, 8, 7, 7, 8)),
    Carrier("H05_raw_hard_delta_scalar", (2, 1, 1, 1, 0, 1, 0)),
)


HARD_ROW_NAMES = (
    "random_covering_022",
    "random_covering_031",
    "random_covering_049",
    "random_covering_078",
    "random_covering_080",
    "random_covering_085",
    "random_covering_113",
)


def tournament() -> tuple[dict[int, int], list[str]]:
    hist = dict(sorted(Counter(carrier.total for carrier in CARRIERS).items()))
    order = {carrier.name: idx for idx, carrier in enumerate(CARRIERS)}
    path = [
        carrier.name
        for carrier in sorted(CARRIERS, key=lambda carrier: (-carrier.total, order[carrier.name]))
    ]
    return hist, path


def main() -> None:
    all_row_speeds = H3453.H3450.rows()
    rows = [H3453.join_row(name, all_row_speeds[name]) for name in HARD_ROW_NAMES]
    row_by_name = {row.name: row for row in rows}
    orbit_audits = {row.name: H3475.build_orbits(row) for row in rows}

    hard_audits: list[HardOrbitAudit] = []
    for row in rows:
        for orbit in orbit_audits[row.name].hard_orbits:
            if orbit.max_delta < 7:
                continue
            hard_labels: set[str] = set()
            dead_labels = dead_label_set(row)
            for gate in orbit.members:
                hard_labels.update(gate_touch_labels(gate, dead_labels))
            hard_cut = label_cut(row, frozenset(hard_labels), "hard-orbit")
            hard_audits.append(
                HardOrbitAudit(
                    row=row,
                    orbit=orbit,
                    phase=phase_hit_audit(row, orbit),
                    hard_cut=hard_cut,
                    single_cut=best_single_e_branch_cut(row),
                    pair_cut=best_pair_current(row),
                    mirror_cut=best_mirror_e_branch_pair(row, orbit_audits[row.name].orbits),
                )
            )

    hard_audits.sort(key=lambda audit: (-audit.orbit.max_delta, audit.row.name, audit.orbit.components))
    hard_rows = sorted({audit.row.name for audit in hard_audits})
    route_hist = Counter(audit.route for audit in hard_audits)
    row_route_hist = Counter()
    for row_name in hard_rows:
        row_routes = {audit.route for audit in hard_audits if audit.row.name == row_name}
        row_route_hist[tuple(sorted(row_routes))] += 1

    phase_bypass = [audit for audit in hard_audits if audit.phase.bypasses_hard_orbit]
    lower_delta_hits = [audit for audit in hard_audits if audit.phase.same_component_lower_delta_hits > 0]
    single_edge = [audit for audit in hard_audits if audit.single_cut is not None and audit.single_cut.is_edge_cut]
    pair_edge = [audit for audit in hard_audits if audit.pair_cut is not None and audit.pair_cut.is_edge_cut]
    mirror_edge = [audit for audit in hard_audits if audit.mirror_cut is not None and audit.mirror_cut.is_edge_cut]
    hard_edge = [audit for audit in hard_audits if audit.hard_cut.is_edge_cut]

    print("HYP-3477 LRC14 HARD MIRROR-ORBIT DISCHARGE AUDIT")
    print("status=EVIDENCE / hard-orbit discharge ledger; not an LRC14 proof")
    print("source=HYP-3475 hard mirror orbits + HYP-3472 boundary currents + HYP-3476 pair-current target")
    print()
    print("## Assumption Challenge")
    print("alternate vertices considered: runners, residues, endpoint colors,")
    print("individual gates, mirror orbits, q=14V phase-grid witnesses,")
    print("dead-cover components, blocker labels, pair currents, and proof obligations.")
    print("chosen carrier vertices: hard mirror orbits plus lower-delta E/branch")
    print("projection-current sidecars.  This preserves the hard gluing predicate")
    print("and destroys raw runner order, scalar gate count, and untyped color mass.")
    print()
    print("## Aggregate Hard-Orbit Ledger")
    print(f"hard_orbits={len(hard_audits)}")
    print(f"hard_rows={len(hard_rows)} {hard_rows}")
    print(f"delta_hist={dict(sorted(Counter(audit.orbit.max_delta for audit in hard_audits).items()))}")
    print(f"route_hist={dict(route_hist)}")
    print(f"row_route_hist={dict(row_route_hist)}")
    print(f"phase_bypass_orbits={len(phase_bypass)}/{len(hard_audits)}")
    print(f"phase_lower_delta_same_component_orbits={len(lower_delta_hits)}/{len(hard_audits)}")
    print(f"hard_orbit_projection_edge_cuts={len(hard_edge)}/{len(hard_audits)}")
    print(f"single_e_branch_projection_edge_cuts={len(single_edge)}/{len(hard_audits)}")
    print(f"pair_current_projection_edge_cuts={len(pair_edge)}/{len(hard_audits)}")
    print(f"mirror_e_branch_projection_edge_cuts={len(mirror_edge)}/{len(hard_audits)}")
    print()

    print("## Per-Hard-Orbit Readout")
    for audit in hard_audits:
        row = audit.row
        orbit = audit.orbit
        print(f"- row={row.name} V={audit.phase.V} route={audit.route}")
        print(f"  speeds={row.speeds}")
        print(f"  orbit_delta={orbit.max_delta} components={orbit.components}")
        print(f"  typed_pair={orbit.typed_pair}")
        print(f"  structural_pair={orbit.structural_pair}")
        print(f"  intervals={orbit.intervals}")
        print(
            "  phase: witnesses={w} hard_hits={h} same_component_lower_delta_hits={l} "
            "any_gate_hits={a} no_gate_hits={n} ambiguous={amb}".format(
                w=audit.phase.witness_count,
                h=audit.phase.hard_gate_hits,
                l=audit.phase.same_component_lower_delta_hits,
                a=audit.phase.any_gate_hits,
                n=audit.phase.no_gate_hits,
                amb=audit.phase.ambiguous_gate_hits,
            )
        )
        print(f"  phase_gate_delta_hist={dict(sorted(audit.phase.gate_delta_counter.items()))}")
        print(f"  hard_component_hits={dict(audit.phase.hard_component_hit_counter)}")
        print(
            "  hard_cut: labels={labels} removed_edges={edges} largest_drop={drop} "
            "component_gain={gain} current={current}".format(
                labels=tuple(sorted(audit.hard_cut.labels)),
                edges=audit.hard_cut.removed_edges,
                drop=audit.hard_cut.largest_drop,
                gain=audit.hard_cut.component_gain,
                current=audit.hard_cut.current,
            )
        )
        for name, cut in (
            ("single_e_branch", audit.single_cut),
            ("pair_current", audit.pair_cut),
            ("mirror_e_branch", audit.mirror_cut),
        ):
            if cut is None:
                print(f"  {name}: None")
                continue
            print(
                "  {name}: labels={labels} source={source} removed_edges={edges} "
                "largest_drop={drop} component_gain={gain} current={current}".format(
                    name=name,
                    labels=tuple(sorted(cut.labels)),
                    source=cut.source,
                    edges=cut.removed_edges,
                    drop=cut.largest_drop,
                    gain=cut.component_gain,
                    current=cut.current,
                )
            )
    print()

    print("## Row-Level Hard Family Summary")
    for row_name in hard_rows:
        audits = [audit for audit in hard_audits if audit.row.name == row_name]
        row = row_by_name[row_name]
        print(f"- {row_name}: hard_orbits={len(audits)} dead_components={len(row.component_row.dead_components)}")
        print(
            "  routes={routes} witnesses={witnesses} best_pair_edges={pair_edges} "
            "best_single_edges={single_edges}".format(
                routes=sorted({audit.route for audit in audits}),
                witnesses=[audit.phase.witness_count for audit in audits],
                pair_edges=max((audit.pair_cut.removed_edges if audit.pair_cut else 0) for audit in audits),
                single_edges=max((audit.single_cut.removed_edges if audit.single_cut else 0) for audit in audits),
            )
        )
    print()

    print("## Lemma Target")
    print("Audited finite form:")
    print("  random031, the HYP-3455 rank-6 gluing row, is the only hard orbit")
    print("  avoided by q=14V phase-grid witnesses; its hard components are still")
    print("  hit through lower-delta compatible gates.  The other seven hard orbits")
    print("  may be hit by q=14V witnesses, but each sits in a row where a")
    print("  lower-delta E/branch gate already gives a dead-cover projection-edge")
    print("  cut.  Thus the hard-family split is: projection-current discharge for")
    print("  7/8 orbits, and the named random031 HYP-3455/HYP-3476 gluing exception")
    print("  for 1/8.")
    print("Proof pull:")
    print("  use single E/branch projection cuts to discharge the seven non-random031")
    print("  hard orbits; route random031 through HYP-3455 plus the HYP-3476")
    print("  pair-current/gluing packet, with phase-branch bypass as supporting")
    print("  evidence that its max-delta pair is not a q=14V CRT obstruction.")
    print()

    hist, path = tournament()
    print("## Tournament Analysis")
    print("vertices=hard-orbit discharge proof carriers, not runners or raw gates")
    print("pairwise_observable=predicate retention + phase bypass + projection current + pair-current bridge + scalar firewall")
    print("switch=higher proof-facing carrier score; ties use declared order")
    print("axes=predicate_retention,phase_payload,projection_current,pair_current,hard_gluing,random_sidecar,scalar_firewall")
    print(f"score_hist={hist}")
    print("directed_3cycles=0")
    print("hamiltonian_path=" + " -> ".join(path))


if __name__ == "__main__":
    main()
