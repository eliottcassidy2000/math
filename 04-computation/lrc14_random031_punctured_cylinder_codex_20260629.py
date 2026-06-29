#!/usr/bin/env python3
"""HYP-3482: random031 punctured-cylinder / forbidden-seam atlas.

This is a focused reframe of the remaining named random031 clause after
HYP-3455, HYP-3460, HYP-3477, and HYP-3479.

The geometric model is:

    random031 = mirror-punctured cylinder
    four dead components = four punctures / islands
    max-delta hard mirror pair = forbidden seam
    q=14V phase witnesses = flow in the seam complement

The script keeps the model honest by printing exact component intervals,
owner labels, hard-seam owners, phase-flow bypass counts, and the additive
owner/complement word that echoes the old n+2 versus n*2 recursion theme.

It is evidence and experiment design, not an LRC14 proof.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
ROW_NAME = "random_covering_031"


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


H3455 = load_module(
    "hyp3455_random031_for_hyp3480",
    COMP / "lrc14_random031_gate_gluing_obstruction_codex_20260629.py",
)
H3460 = load_module(
    "hyp3460_phase_branch_for_hyp3480",
    COMP / "lrc14_phase_branch_color_pullback_codex_20260629.py",
)
H3477 = load_module(
    "hyp3477_hard_orbit_for_hyp3480",
    COMP / "lrc14_hard_mirror_orbit_discharge_codex_20260629.py",
)


Interval = tuple[Fraction, Fraction]


def fmt(value: Fraction | int | None) -> str:
    if value is None:
        return "None"
    if isinstance(value, int):
        return str(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def interval_text(interval: Interval | None) -> str:
    if interval is None:
        return "None"
    return f"[{fmt(interval[0])},{fmt(interval[1])}]"


def mirror_interval(interval: Interval) -> Interval:
    return (1 - interval[1], 1 - interval[0])


def cover_speeds(cover) -> tuple[int, ...]:
    return tuple() if cover is None else tuple(cover[1])


def component_owner_pair(component) -> tuple[int, ...]:
    return tuple(sorted(cover_speeds(component.b0_cover) + cover_speeds(component.b1_cover)))


def component_label_word(component) -> str:
    b0 = cover_speeds(component.b0_cover)
    b1 = cover_speeds(component.b1_cover)
    return f"B0={b0 or '-'} B1={b1 or '-'}"


def residue_pair(pair: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(owner % 14 for owner in pair)


def mirror_index(components: list, interval: Interval) -> int | None:
    target = mirror_interval(interval)
    for index, component in enumerate(components):
        if component.interval == target:
            return index
    return None


def owner_complement_word(V: int, owners: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(V - owner for owner in owners)


def hard_audit():
    speeds = H3477.H3453.H3450.rows()[ROW_NAME]
    row = H3477.H3453.join_row(ROW_NAME, speeds)
    orbit_audit = H3477.H3475.build_orbits(row)
    hard_orbits = [orbit for orbit in orbit_audit.hard_orbits if orbit.max_delta >= 7]
    if len(hard_orbits) != 1:
        raise RuntimeError(f"expected one random031 hard orbit, got {len(hard_orbits)}")
    orbit = hard_orbits[0]
    phase = H3477.phase_hit_audit(row, orbit)
    return row, orbit, phase


def phase_row():
    for row in H3460.ROWS:
        if row.label == "random031":
            return row
    raise RuntimeError("missing random031 phase row")


@dataclass(frozen=True)
class Carrier:
    name: str
    scores: tuple[int, ...]

    @property
    def total(self) -> int:
        return sum(self.scores)


AXES = (
    "predicate_retention",
    "topology_payload",
    "seam_payload",
    "phase_flow_payload",
    "recursion_payload",
    "owner_current_payload",
    "scalar_firewall",
)

CARRIERS = (
    Carrier("C00_forbidden_seam_complement_flow", (10, 10, 10, 10, 8, 9, 10)),
    Carrier("C01_mirror_punctured_cylinder_model", (10, 10, 9, 9, 8, 8, 10)),
    Carrier("C02_seven_owner_layered_seam_word", (9, 8, 10, 8, 9, 10, 10)),
    Carrier("C03_phase_branch_bypass_worldlines", (9, 8, 8, 10, 8, 8, 9)),
    Carrier("C04_additive_rim_vs_doubling_fold_lens", (8, 7, 7, 8, 10, 7, 9)),
    Carrier("C05_dead_island_owner_pairs", (8, 9, 7, 7, 7, 8, 9)),
    Carrier("C06_raw_counts_shadow", (2, 1, 2, 1, 1, 1, 0)),
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
    speeds = H3455.H3450.rows()[ROW_NAME]
    V = max(speeds)
    component_audit = H3455.H3450.audit_row(ROW_NAME, speeds)
    dead_components = list(component_audit.dead_components)
    cut = H3455.H3437.audit_row(ROW_NAME, speeds)
    rescue_graph = H3455.rescue_graph(cut)
    row, orbit, phase = hard_audit()
    phase_metrics = H3460.phase_metrics(phase_row())
    pullback = H3460.pullback_metrics(phase_row())

    hard_owners = tuple(sorted({owner for gate in orbit.members for owner in H3455.cover_owners(gate)}))
    dead_owners = tuple(sorted({owner for component in dead_components for owner in component_owner_pair(component)}))
    rescue_owners = tuple(cut.rescue_subset)
    rescue_only = tuple(owner for owner in rescue_owners if owner not in dead_owners)
    apex_owners = tuple(owner for owner in hard_owners if owner not in rescue_owners)

    dead_pairs = tuple(component_owner_pair(component) for component in dead_components)
    dead_spans = tuple(max(pair) - min(pair) for pair in dead_pairs)
    unique_rim_pairs = ((23, 45), (93, 113), (147, 169))
    rim_complement_pairs = tuple(owner_complement_word(V, pair) for pair in unique_rim_pairs)
    rim_spans = tuple(abs(pair[1] - pair[0]) for pair in unique_rim_pairs)

    hist, path = tournament()

    print("HYP-3482 RANDOM031 PUNCTURED-CYLINDER / FORBIDDEN-SEAM ATLAS")
    print("status=EVIDENCE / topology reframe and experiment design; not an LRC14 proof")
    print("source=HYP-3455 + HYP-3460 + HYP-3477 + HYP-3479")
    print()
    print("## Assumption Challenge")
    print("alternate vertices considered: runners, residues, dead islands, seam gates,")
    print("phase-flow worldlines, branch sheets, owner labels, fixed circle sections,")
    print("section boundaries, wall-crossing events, cover arcs, Fourier modes, and")
    print("terminal proof obligations.")
    print("chosen carrier vertices: punctures/islands, forbidden seam, seam-complement")
    print("phase flow, and owner-layer words.  This preserves the random031 terminal")
    print("clause predicate and destroys raw runner order, raw gate count, and scalar")
    print("rank-6 danger unless the seam, puncture, and phase-flow sidecars are kept.")
    print()

    print("## Base Row")
    print(f"row={ROW_NAME}")
    print(f"speeds={tuple(sorted(speeds))}")
    print(f"V={V}")
    print(f"S3_P={phase_metrics.P}")
    print(f"S3_E={phase_metrics.E}")
    print(f"components={len(component_audit.components)} dead_components={len(dead_components)}")
    print(f"class_hist={dict(sorted(component_audit.class_hist.items()))}")
    print(f"low_rank_escape={H3455.component_graph_summary(component_audit).low_rank_escape}")
    print()

    print("## Punctured-Cylinder Model")
    print("cell_model=cylinder chi=0; four isolated dead islands puncture it, so chi_model=-4")
    print(f"puncture_count={len(dead_components)}")
    print(f"mirror_puncture_pairs={len(dead_components) // 2}")
    print("dead_islands:")
    for index, component in enumerate(dead_components):
        pair = component_owner_pair(component)
        print(
            "  island={i} interval={interval} mirror={mirror} owners={owners} "
            "residues={residues} span={span} {labels}".format(
                i=index,
                interval=interval_text(component.interval),
                mirror=mirror_index(dead_components, component.interval),
                owners=pair,
                residues=residue_pair(pair),
                span=max(pair) - min(pair),
                labels=component_label_word(component),
            )
        )
    print(f"dead_owner_union={dead_owners}")
    print(f"dead_owner_span_word={dead_spans}")
    print()

    print("## Forbidden Seam")
    print(f"hard_orbit_components={orbit.components}")
    print(f"hard_orbit_delta={orbit.max_delta}")
    print(f"hard_orbit_typed_pair={orbit.typed_pair}")
    print(f"hard_orbit_structural_pair={orbit.structural_pair}")
    print(f"hard_orbit_intervals={orbit.intervals}")
    print(f"hard_seam_owner_union={hard_owners}")
    print(f"rescue_subset={rescue_owners}")
    print(f"dead_island_owners_inside_seam={set(dead_owners).issubset(hard_owners)}")
    print(f"rescue_subset_inside_seam={set(rescue_owners).issubset(hard_owners)}")
    print(f"rescue_only_outer_owners={rescue_only}")
    print(f"apex_owner_layer={apex_owners}")
    print(f"rescue_graph_edges={rescue_graph.edge_count}")
    print(f"rescue_graph_connected={rescue_graph.connected}")
    print(f"rescue_graph_components={rescue_graph.component_sizes}")
    for gate in orbit.members:
        print(
            "  seam_gate component={component} interval={interval} branch={branch} "
            "delta=({b0},{b1}) owners={owners}".format(
                component=gate.component_index,
                interval=interval_text(gate.interval),
                branch=gate.branch_mask,
                b0=gate.b0_delta,
                b1=gate.b1_delta,
                owners=H3455.cover_owners(gate),
            )
        )
    print()

    print("## Phase Flow In The Seam Complement")
    print(f"phase_witnesses={phase.witness_count}")
    print(f"phase_actual_count={phase_metrics.actual_count}")
    print(f"phase_open_count={phase_metrics.open_count}")
    print(f"phase_deficit={fmt(phase_metrics.actual_deficit)}")
    print(f"phase_counts={H3460.sparse(phase_metrics.phase_counts)}")
    print(f"phase_branch_mirror_failures={pullback.mirror_branch_failures}")
    print(f"hard_gate_hits={phase.hard_gate_hits}")
    print(f"same_component_lower_delta_hits={phase.same_component_lower_delta_hits}")
    print(f"hard_component_hit_counter={dict(phase.hard_component_hit_counter)}")
    print(f"gate_delta_counter_on_phase_flow={dict(sorted(phase.gate_delta_counter.items()))}")
    print(f"gate_route_counter_on_phase_flow={dict(sorted(phase.gate_route_counter.items()))}")
    print(f"no_gate_hits={phase.no_gate_hits}")
    print("interpretation=the max-delta seam is forbidden to phase flow; the same")
    print("hard components are touched 12 times through branch-opposite lower-delta")
    print("gates, so the complement carries the actual flow.")
    print()

    print("## Additive Rim vs Doubling Fold")
    print("additive_owner_rim_pairs=(23,45),(93,113),(147,169) plus apex 173")
    print(f"additive_owner_rim_spans={rim_spans}")
    print(f"E_complement_pairs={rim_complement_pairs} plus apex_complement={V - 173}")
    print("readout=the additive rim is a +22,+20,+22 word in owner space, while")
    print("the multiplicative fold is the two-adic map u=2t mod1 used by the")
    print("phase pullback.  In this row, n+2-style rim data names the seam owners;")
    print("n*2-style folding names the branch sheet on which the phase flow bypasses")
    print("the seam.")
    print()

    print("## Bold Reframes")
    print("1. Forbidden-seam theorem target:")
    print("   If a max-delta mirror seam has zero compatible q=14V phase hits and")
    print("   its complement has mirror-balanced lower-delta same-component hits,")
    print("   then the seam is a boundary-gluing clause, not a wall obstruction.")
    print("2. Punctured-cylinder theorem target:")
    print("   A hard/currentless row with four mirror-paired dead islands, connected")
    print("   rescue graph, and a seven-owner seam should discharge by proving that")
    print("   every attempted seam saturation cuts off one of the 94 low-rank escapes")
    print("   or violates the phase-branch mirror flow.")
    print("3. Experiment design:")
    print("   Build a seam-complement graph whose vertices are survivor gates minus")
    print("   the hard seam.  Send the 282 phase witnesses through it and test whether")
    print("   every connected phase-flow component reaches a low-rank escape before")
    print("   it can cross the forbidden seam.  Then repeat on all hard orbits and")
    print("   classify zero-hit seams versus genuine phase walls.")
    print()

    print("## Tournament Analysis")
    print("vertices=topological proof carriers, not runners or raw gates")
    print("pairwise_observable=predicate retention + topology + seam + phase flow + recursion lens + owner current + scalar firewall")
    print("switch=higher retained proof payload; ties use declared route order")
    print(f"axes={AXES}")
    print(f"score_hist={hist}")
    print("directed_3cycles=0")
    print("hamiltonian_path=" + " -> ".join(path))


if __name__ == "__main__":
    main()
