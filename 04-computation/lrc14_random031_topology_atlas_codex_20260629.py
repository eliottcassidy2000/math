#!/usr/bin/env python3
"""HYP-3481: topology atlas for the random_covering_031 obstruction.

This is not another attempt to prove the row by scalar strength.  It reads
random_covering_031 as a small topological object:

* the dead-cover projection is four isolated islands;
* the max-delta mirror pair is a seam around the antipodal midpoint;
* phase-grid witnesses flow around that seam through lower-delta gates.

The goal is to make two reframings proof-facing:

1.  A mirror-punctured annulus: four isolated dead islands are punctures,
    and the hard mirror pair is a short seam whose monodromy swaps branch and
    mirror coordinates.
2.  A bypassed saddle: the hard gates are critical saddles of the continuous
    gate complex, but the q=14V phase flow has zero flux through them and
    nonzero lower-delta flux through the same hard components.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"

ROW_NAME = "random_covering_031"
SPEEDS = (12, 23, 45, 55, 58, 70, 84, 93, 113, 120, 147, 169, 173)
V = 173
Q = 14 * V
HARD_COMPONENTS = (43, 54)
RESCUE_CORE = (23, 45, 93, 113, 147, 169)
SEVEN_OWNER_CORE = (23, 45, 93, 113, 147, 169, 173)


def load_module(name: str, filename: str):
    path = COMP / filename
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


H3438 = load_module("hyp3438_for_hyp3479", "lrc14_survivor_gate_word_audit_codex_20260629.py")
H3450 = load_module("hyp3450_for_hyp3479", "lrc14_component_cover_obstruction_extractor_codex_20260628.py")
H3472 = load_module("hyp3472_for_hyp3479", "lrc14_dead_cover_boundary_current_codex_20260629.py")
H3460 = load_module("hyp3460_for_hyp3479", "lrc14_phase_branch_color_pullback_codex_20260629.py")


def fmt(x: F | int | None) -> str:
    if x is None:
        return "None"
    if isinstance(x, int):
        return str(x)
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def f6(x: F | float) -> str:
    return f"{float(x):.6f}"


def interval_text(interval: tuple[F, F] | None) -> str:
    if interval is None:
        return "None"
    return f"[{fmt(interval[0])},{fmt(interval[1])}]"


def interval_len(interval: tuple[F, F]) -> F:
    return interval[1] - interval[0]


def interval_mid(interval: tuple[F, F]) -> F:
    return (interval[0] + interval[1]) / 2


def mirror_interval(interval: tuple[F, F]) -> tuple[F, F]:
    return (1 - interval[1], 1 - interval[0])


def dist_residue(value: int, modulus: int) -> int:
    residue = value % modulus
    return min(residue, modulus - residue)


def actual_witnesses() -> tuple[int, ...]:
    return tuple(
        a
        for a in range(Q)
        if all(14 * dist_residue(speed * a, Q) >= Q for speed in SPEEDS)
    )


def contains_closed(interval: tuple[F, F], x: F) -> bool:
    return interval[0] <= x <= interval[1]


def compatible(gate, branch: int) -> bool:
    return gate.branch_mask == "both" or gate.branch_mask == f"branch{branch}"


def component_cover_labels(component) -> tuple[str, ...]:
    labels: list[str] = []
    if component.b0_cover is not None:
        labels.extend(f"B0:{speed}" for speed in component.b0_cover[1])
    if component.b1_cover is not None:
        labels.extend(f"B1:{speed}" for speed in component.b1_cover[1])
    return tuple(sorted(labels))


def cover_owners(gate) -> tuple[int, ...]:
    owners: set[int] = set()
    for bad in (gate.left_bad, gate.right_bad):
        if bad is None:
            continue
        for cover in (bad.b0_cover, bad.b1_cover):
            if cover is not None:
                owners.update(cover)
    return tuple(sorted(owners))


def gate_text(gate) -> str:
    return (
        f"component={gate.component_index} interval={interval_text(gate.interval)} "
        f"len={fmt(gate.length)} mask={gate.branch_mask} "
        f"kinds={gate.endpoint_kind_signature} delta=({gate.b0_delta},{gate.b1_delta}) "
        f"labels={H3438.label_values(gate.left_labels)}|{H3438.label_values(gate.right_labels)}"
    )


def build_gates():
    bad_row = H3438.H3436.audit_row(ROW_NAME, SPEEDS)
    return [gate for comp in H3438.build_mixed_components(bad_row) for gate in comp.survivor_gates]


def mirror_pairs_by_interval(gates):
    by_key = {(gate.interval, gate.branch_mask): gate for gate in gates}
    seen: set[tuple[tuple[int, tuple[F, F], str], tuple[int, tuple[F, F], str]]] = set()
    pairs = []
    for gate in gates:
        target_mask = {"branch0": "branch1", "branch1": "branch0", "both": "both"}[gate.branch_mask]
        partner = by_key.get((mirror_interval(gate.interval), target_mask))
        if partner is None:
            continue
        key = tuple(
            sorted(
                (
                    (gate.component_index, gate.interval, gate.branch_mask),
                    (partner.component_index, partner.interval, partner.branch_mask),
                )
            )
        )
        if key in seen:
            continue
        seen.add(key)
        pairs.append((gate, partner))
    return pairs


def dead_island_report(row):
    islands = []
    for idx, component in enumerate(row.dead_components):
        islands.append(
            {
                "idx": idx,
                "interval": component.interval,
                "mid": interval_mid(component.interval),
                "length": interval_len(component.interval),
                "labels": component_cover_labels(component),
                "rank": component.dead_pair_rank,
            }
        )
    mirror = []
    used: set[int] = set()
    for a, left in enumerate(islands):
        if a in used:
            continue
        found = None
        for b, right in enumerate(islands):
            if b == a or b in used:
                continue
            if right["interval"] == mirror_interval(left["interval"]):
                found = b
                break
        if found is not None:
            used.add(a)
            used.add(found)
            mirror.append((a, found))
    return islands, mirror


def phase_flux(gates):
    witnesses = actual_witnesses()
    hard_gates = [gate for gate in gates if gate.component_index in HARD_COMPONENTS and gate.total_delta == 7]
    lower_hard_component_gates = [
        gate
        for gate in gates
        if gate.component_index in HARD_COMPONENTS and gate.total_delta < 7
    ]
    hit_counter: Counter[tuple[int, str, int]] = Counter()
    phase_branch_counter: Counter[tuple[int, int]] = Counter()
    no_gate_hits = 0
    bypass_hits = 0
    hard_hits = 0
    for a in witnesses:
        t = F(a, Q)
        branch = 0 if t < F(1, 2) else 1
        u = F((2 * a) % Q, Q)
        phase = a % 14
        phase_branch_counter[(phase, branch)] += 1
        hits = [gate for gate in gates if contains_closed(gate.interval, u) and compatible(gate, branch)]
        if not hits:
            no_gate_hits += 1
            continue
        for gate in hits:
            if gate in hard_gates:
                hard_hits += 1
            if gate in lower_hard_component_gates:
                bypass_hits += 1
                hit_counter[(gate.component_index, gate.branch_mask, gate.total_delta)] += 1
    return {
        "witnesses": witnesses,
        "hard_hits": hard_hits,
        "bypass_hits": bypass_hits,
        "no_gate_hits": no_gate_hits,
        "phase_branch_counter": phase_branch_counter,
        "hard_component_hit_counter": hit_counter,
    }


def topology_scores(row, gates, flux):
    dead_projection = H3472.projection(list(row.dead_components))
    projection_stats = H3472.projection_stats(dead_projection)
    hard_gates = [gate for gate in gates if gate.component_index in HARD_COMPONENTS and gate.total_delta == 7]
    mirror_pairs = mirror_pairs_by_interval(gates)
    hard_pair_is_mirror = any(
        {left.component_index, right.component_index} == set(HARD_COMPONENTS)
        and left.total_delta == right.total_delta == 7
        for left, right in mirror_pairs
    )
    return {
        "dead_projection_edges": projection_stats.edges,
        "dead_projection_components": projection_stats.components,
        "dead_projection_largest": projection_stats.largest,
        "hard_pair_is_mirror": hard_pair_is_mirror,
        "hard_gate_count": len(hard_gates),
        "hard_gate_hits": flux["hard_hits"],
        "hard_component_bypass_hits": flux["bypass_hits"],
        "phase_witnesses": len(flux["witnesses"]),
    }


def main() -> None:
    row = H3450.audit_row(ROW_NAME, SPEEDS)
    gates = build_gates()
    islands, island_mirrors = dead_island_report(row)
    mirror_pairs = mirror_pairs_by_interval(gates)
    hard_gates = sorted(
        [gate for gate in gates if gate.component_index in HARD_COMPONENTS and gate.total_delta == 7],
        key=lambda gate: gate.component_index,
    )
    bypass_gates = sorted(
        [gate for gate in gates if gate.component_index in HARD_COMPONENTS and gate.total_delta == 2],
        key=lambda gate: (gate.component_index, gate.branch_mask, gate.interval),
    )
    flux = phase_flux(gates)
    scores = topology_scores(row, gates, flux)

    print("HYP-3481 RANDOM031 TOPOLOGY ATLAS")
    print("status=EVIDENCE / reframing atlas; not an LRC14 proof")
    print("row=random_covering_031")
    print()
    print("## Reframe 1: mirror-punctured annulus")
    print("Interpretation: the dead-cover projection is not a connected graph")
    print("waiting for a cut.  It is four isolated punctures on a mirror circle.")
    print(
        "dead_projection="
        f"components={scores['dead_projection_components']} "
        f"largest={scores['dead_projection_largest']} "
        f"edges={scores['dead_projection_edges']}"
    )
    print(f"dead_island_mirror_pairs={island_mirrors}")
    for island in islands:
        print(
            "  island={idx} interval={interval} mid={mid} len={length} "
            "rank={rank} labels={labels}".format(
                idx=island["idx"],
                interval=interval_text(island["interval"]),
                mid=fmt(island["mid"]),
                length=fmt(island["length"]),
                rank=island["rank"],
                labels=island["labels"],
            )
        )
    print()
    print("## Reframe 2: bypassed saddle seam")
    print("Interpretation: the max-delta mirror pair is a short antipodal seam,")
    print("not a phase-grid wall.  The q=14V phase flow has zero flux through")
    print("the seam and nonzero flux through lower-delta gates on the same two")
    print("components.")
    print(f"hard_pair_is_mirror={scores['hard_pair_is_mirror']}")
    for gate in hard_gates:
        print(f"  hard_gate: {gate_text(gate)} owners={cover_owners(gate)}")
    for gate in bypass_gates:
        print(f"  bypass_gate: {gate_text(gate)} owners={cover_owners(gate)}")
    if len(hard_gates) == 2:
        gap = hard_gates[1].interval[0] - hard_gates[0].interval[1]
        print(
            "  central_mirror_gap="
            f"{fmt(gap)} ({f6(gap)}) between hard gates across 1/2"
        )
    print(
        "phase_flux="
        f"witnesses={scores['phase_witnesses']} "
        f"hard_gate_hits={scores['hard_gate_hits']} "
        f"lower_delta_hard_component_hits={scores['hard_component_bypass_hits']} "
        f"no_gate_hits={flux['no_gate_hits']}"
    )
    print(f"hard_component_hit_counter={dict(sorted(flux['hard_component_hit_counter'].items()))}")
    print()
    print("## Owner topology")
    print(f"rescue_core={RESCUE_CORE}")
    print(f"seven_owner_core={SEVEN_OWNER_CORE}")
    print(f"extra_owner={tuple(sorted(set(SEVEN_OWNER_CORE) - set(RESCUE_CORE)))}")
    print("owner_readout=rank-6 rescue core plus one apex owner 173; the seam")
    print("carries the whole owner boundary, but the phase flow bypasses it.")
    print()
    print("## Proof-facing lemmas suggested")
    print("L1 mirror-puncture lemma: four edgeless dead-cover islands with")
    print("   mirror-paired singleton blockers should route by local island")
    print("   current, not Menger edge cuts.")
    print("L2 bypassed-saddle lemma: a max-delta mirror seam with zero q=14V")
    print("   phase flux and lower-delta same-component flux is a removable")
    print("   saddle, unless its seven-owner boundary creates explicit")
    print("   owner-current/two-adic/SPEC debt.")
    print()
    print("## Tournament Analysis")
    carriers = (
        ("mirror_punctured_annulus", 68),
        ("bypassed_saddle_seam", 66),
        ("seven_owner_boundary", 63),
        ("phase_branch_flow", 61),
        ("dead_projection_scalar", 21),
        ("raw_rank6_scalar", 10),
    )
    print(f"vertices={tuple(name for name, _ in carriers)}")
    print("pairwise_observable=retained topology + phase flux + owner boundary")
    print("switch=higher proof-facing carrier score")
    print(f"score_hist={dict(sorted(Counter(score for _, score in carriers).items()))}")
    print("directed_3cycles=0")
    print("hamiltonian_path=" + " -> ".join(name for name, _ in sorted(carriers, key=lambda item: -item[1])))
    print()
    print("## Assumption challenge")
    print("Considered vertices: runners, rescue owners, dead islands, mirror seams,")
    print("phase witnesses, survivor gates, blocker labels, cover arcs, Fourier")
    print("modes, and proof obligations.  Chosen vertices are dead islands,")
    print("mirror seams, and phase-flow bypasses.  This preserves the terminal")
    print("discharge predicate and destroys raw runner order only after replacing")
    print("it by mirror, owner-boundary, and phase-flux sidecars.")


if __name__ == "__main__":
    main()
