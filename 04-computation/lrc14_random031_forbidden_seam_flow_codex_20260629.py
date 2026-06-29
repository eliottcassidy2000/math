#!/usr/bin/env python3
"""HYP-3484: random031 as forbidden seam and mirror-punctured cylinder.

This scout is deliberately geometric.  HYP-3455 isolated the rank-6
random_covering_031 obstruction as a seven-owner max-delta mirror gate pair,
while HYP-3460/HYP-3477 showed that the q=14V phase witnesses avoid those
max-delta gates and touch the same components only through lower-delta
opposite-branch gates.

Here the working picture is:

    four isolated dead components  = mirror-paired punctures,
    max-delta mirror gate pair     = forbidden seam,
    q=14V phase witnesses          = flow on the seam complement.

The script recomputes that packet, performs two gate-surgery tests, and records
the carry/doubling coordinates suggested by the older n+2 versus n*2 recursion
thread.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
ROW_NAME = "random_covering_031"


def load_module(name: str, filename: str):
    path = COMP / filename
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


H3455 = load_module(
    "hyp3455_random031_for_hyp3480",
    "lrc14_random031_gate_gluing_obstruction_codex_20260629.py",
)
H3460 = load_module(
    "hyp3460_phase_pullback_for_hyp3480",
    "lrc14_phase_branch_color_pullback_codex_20260629.py",
)

H3438 = H3455.H3438
H3450 = H3455.H3450


def fmt(x: F | int | None) -> str:
    if x is None:
        return "None"
    if isinstance(x, int):
        return str(x)
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def interval_text(interval: tuple[F, F] | None) -> str:
    if interval is None:
        return "None"
    return f"[{fmt(interval[0])},{fmt(interval[1])}]"


def sparse(counter: Counter) -> dict:
    return {key: counter[key] for key in sorted(counter) if counter[key]}


def contains_closed(interval: tuple[F, F], x: F) -> bool:
    return interval[0] <= x <= interval[1]


def compatible(gate, branch: int) -> bool:
    return gate.branch_mask == "both" or gate.branch_mask == f"branch{branch}"


def gate_key(gate) -> tuple[int, tuple[F, F], str]:
    return (gate.component_index, gate.interval, gate.branch_mask)


def gate_sort_key(gate) -> tuple[int, F, str]:
    return (gate.component_index, gate.interval[0], gate.branch_mask)


def labels_for_component(component) -> tuple[str, ...]:
    labels: list[str] = []
    if component.b0_cover is not None:
        labels.extend(f"B0:{owner}" for owner in component.b0_cover[1])
    if component.b1_cover is not None:
        labels.extend(f"B1:{owner}" for owner in component.b1_cover[1])
    return tuple(sorted(labels))


def owners_for_component(component) -> tuple[int, ...]:
    owners: set[int] = set()
    if component.b0_cover is not None:
        owners.update(component.b0_cover[1])
    if component.b1_cover is not None:
        owners.update(component.b1_cover[1])
    return tuple(sorted(owners))


def speed_from_label(label: str) -> int | None:
    if ":" not in label:
        return None
    head, value = label.split(":", 1)
    if head not in {"B0", "B1", "E"}:
        return None
    try:
        return int(value)
    except ValueError:
        return None


def gate_endpoint_speeds(gate) -> tuple[int, ...]:
    speeds = []
    for label in gate.left_labels + gate.right_labels:
        speed = speed_from_label(label)
        if speed is not None:
            speeds.append(speed)
    return tuple(sorted(set(speeds)))


def gate_word(gate) -> str:
    return (
        f"c{gate.component_index}:{interval_text(gate.interval)}:"
        f"{H3438.label_kinds(gate.left_labels)}|{H3438.label_kinds(gate.right_labels)}:"
        f"mask={gate.branch_mask}:adj={gate.adjacency}:"
        f"d({gate.b0_delta},{gate.b1_delta}):route={gate.route}:"
        f"labels={H3438.label_values(gate.left_labels)}|{H3438.label_values(gate.right_labels)}:"
        f"owners={H3455.cover_owners(gate)}"
    )


@dataclass(frozen=True)
class PhaseHit:
    a: int
    phase: int
    branch: int
    u: F
    gate_key: tuple[int, tuple[F, F], str]
    component: int
    gate_branch: str
    total_delta: int
    route: str


@dataclass(frozen=True)
class Surgery:
    any_gate_hits: int
    no_gate_hits: int
    ambiguous_gate_hits: int
    route_counter: Counter[str]
    delta_counter: Counter[int]


@dataclass(frozen=True)
class Carrier:
    name: str
    scores: tuple[int, ...]

    @property
    def total(self) -> int:
        return sum(self.scores)


CARRIERS = (
    Carrier("F00_forbidden_seam_complement_flow", (10, 10, 10, 10, 9, 10, 10)),
    Carrier("F01_lower_delta_mirror_bypass_channel", (10, 10, 9, 9, 10, 9, 10)),
    Carrier("F02_mirror_punctured_cylinder", (10, 9, 8, 10, 8, 8, 10)),
    Carrier("F03_c27_carry_vs_phase_doubling_split", (9, 8, 10, 8, 9, 10, 10)),
    Carrier("F04_seven_owner_gluing_clause", (9, 10, 10, 8, 8, 9, 10)),
    Carrier("F05_small_touch_singleton_pocket_shadow", (8, 8, 8, 9, 8, 7, 10)),
    Carrier("F06_raw_wall_or_delta_scalar", (2, 1, 1, 0, 1, 0, 0)),
)


def tournament() -> tuple[dict[int, int], list[str]]:
    hist = dict(sorted(Counter(carrier.total for carrier in CARRIERS).items()))
    order = {carrier.name: idx for idx, carrier in enumerate(CARRIERS)}
    path = [
        carrier.name
        for carrier in sorted(CARRIERS, key=lambda carrier: (-carrier.total, order[carrier.name]))
    ]
    return hist, path


def phase_hits_for_gates(
    witnesses: tuple[int, ...],
    gates: list,
    q: int,
    removed: set[tuple[int, tuple[F, F], str]] | None = None,
) -> tuple[Surgery, list[PhaseHit]]:
    removed = removed or set()
    any_gate_hits = 0
    no_gate_hits = 0
    ambiguous_gate_hits = 0
    route_counter: Counter[str] = Counter()
    delta_counter: Counter[int] = Counter()
    hits_out: list[PhaseHit] = []
    active_gates = [gate for gate in gates if gate_key(gate) not in removed]

    for a in witnesses:
        t = F(a, q)
        branch = 0 if t < F(1, 2) else 1
        u = F((2 * a) % q, q)
        hits = [
            gate
            for gate in active_gates
            if contains_closed(gate.interval, u) and compatible(gate, branch)
        ]
        if not hits:
            no_gate_hits += 1
            continue
        any_gate_hits += 1
        if len(hits) > 1:
            ambiguous_gate_hits += 1
        gate = sorted(hits, key=gate_sort_key)[0]
        route_counter[gate.route] += 1
        delta_counter[gate.total_delta] += 1
        hits_out.append(
            PhaseHit(
                a=a,
                phase=a % 14,
                branch=branch,
                u=u,
                gate_key=gate_key(gate),
                component=gate.component_index,
                gate_branch=gate.branch_mask,
                total_delta=gate.total_delta,
                route=gate.route,
            )
        )

    return (
        Surgery(
            any_gate_hits=any_gate_hits,
            no_gate_hits=no_gate_hits,
            ambiguous_gate_hits=ambiguous_gate_hits,
            route_counter=route_counter,
            delta_counter=delta_counter,
        ),
        hits_out,
    )


def dead_projection_edges(row) -> int:
    projection = H3455.dead_projection(row)
    return sum(len(neighbors) for neighbors in projection.values()) // 2


def mirror_partner_index(interval: tuple[F, F], components: list) -> int | None:
    target = H3455.mirror_interval(interval)
    for idx, component in enumerate(components):
        if component.interval == target:
            return idx
    return None


def main() -> None:
    speeds = H3450.rows()[ROW_NAME]
    component_row = H3450.audit_row(ROW_NAME, speeds)
    bad_row = H3438.H3436.audit_row(ROW_NAME, speeds)
    gates = [gate for component in H3438.build_mixed_components(bad_row) for gate in component.survivor_gates]
    max_delta = max(gate.total_delta for gate in gates)
    seam_gates = sorted((gate for gate in gates if gate.total_delta == max_delta), key=gate_sort_key)
    seam_keys = {gate_key(gate) for gate in seam_gates}
    seam_components = {gate.component_index for gate in seam_gates}
    seam_owners = tuple(sorted({owner for gate in seam_gates for owner in H3455.cover_owners(gate)}))

    phase_row = next(row for row in H3460.ROWS if row.label == "random031")
    phase_metrics = H3460.phase_metrics(phase_row)
    pullback = H3460.pullback_metrics(phase_row)
    witnesses = H3460.actual_witnesses(phase_row.speeds, phase_row.V)
    q = 14 * phase_row.V

    base_surgery, base_hits = phase_hits_for_gates(witnesses, gates, q)
    seam_hits = [hit for hit in base_hits if hit.gate_key in seam_keys]
    bypass_hits = [
        hit
        for hit in base_hits
        if hit.component in seam_components and hit.gate_key not in seam_keys and hit.total_delta < max_delta
    ]
    bypass_keys = {hit.gate_key for hit in bypass_hits}
    seam_deleted, _ = phase_hits_for_gates(witnesses, gates, q, seam_keys)
    bypass_deleted, _ = phase_hits_for_gates(witnesses, gates, q, bypass_keys)
    seam_and_bypass_deleted, _ = phase_hits_for_gates(witnesses, gates, q, seam_keys | bypass_keys)

    components = list(component_row.components)
    dead_components = [(idx, component) for idx, component in enumerate(components) if component.component_class == "dead_both"]
    dead_mirror_pairs = {
        tuple(sorted((idx, mirror_partner_index(component.interval, components))))
        for idx, component in dead_components
        if mirror_partner_index(component.interval, components) is not None
    }

    seam_endpoint_speeds = tuple(sorted({speed for gate in seam_gates for speed in gate_endpoint_speeds(gate)}))
    seam_owner_mod14 = tuple((owner, owner % 14) for owner in seam_owners)
    seam_owner_mod27 = tuple((owner, owner % 27) for owner in seam_owners)
    seam_owner_mod7 = tuple((owner, owner % 7) for owner in seam_owners)
    c27_duplicate_residues = {
        residue: count
        for residue, count in Counter(owner % 27 for owner in seam_owners).items()
        if count > 1
    }
    triangular_identity_lhs = 8 * (14 * 13 // 2) + 1
    triangular_identity_rhs = (2 * 14 - 1) ** 2
    phase_doubling_hist = Counter((phase, (2 * phase) % 14) for phase in phase_metrics.phase_counts)
    witness_gap_hist = Counter(
        (witnesses[(idx + 1) % len(witnesses)] - witnesses[idx]) % q
        for idx in range(len(witnesses))
    )

    hist, path = tournament()

    print("HYP-3484 RANDOM031 FORBIDDEN-SEAM FLOW GEOMETRY")
    print("=" * 78)
    print("status=EVIDENCE / geometric-topological reframing; not an LRC14 proof")
    print("source=HYP-3455 random031 gluing + HYP-3460 phase pullback + HYP-3477 hard-orbit discharge")
    print()

    print("## Assumption Challenge")
    print("alternate vertices considered: runners, gaps, fixed circle sections,")
    print("section boundaries, wall-crossing events, residues, cover arcs,")
    print("Fourier modes, matroid circuits, survivor gates, dead islands,")
    print("phase-flow hits, carry residues, and proof obligations.")
    print("chosen carrier vertices: four dead punctures, the max-delta mirror seam,")
    print("the lower-delta mirror bypass channel, C=27 carry coordinates, and")
    print("the phase-flow complement.  This preserves the LRC branch relocation")
    print("predicate plus the q=14V CRT witness predicate, while destroying raw")
    print("runner order and raw wall-count scalars.")
    print()

    print("## Mirror-Punctured Cylinder")
    print(f"row={ROW_NAME}")
    print(f"speeds={tuple(sorted(speeds))}")
    print(f"components={len(component_row.components)}")
    print(f"component_class_hist={dict(sorted(component_row.class_hist.items()))}")
    print(f"dead_islands={len(dead_components)}")
    print(f"dead_mirror_pairs={tuple(sorted(dead_mirror_pairs))}")
    print(f"dead_projection_edges={dead_projection_edges(component_row)}")
    print(f"projection_components={len(H3455.connected_components(H3455.dead_projection(component_row)))}")
    for idx, component in dead_components:
        print(
            "  puncture=i={idx} interval={interval} mirror={mirror} "
            "labels={labels} owners={owners} owner_mod14={mod14} owner_mod27={mod27}".format(
                idx=idx,
                interval=interval_text(component.interval),
                mirror=mirror_partner_index(component.interval, components),
                labels=labels_for_component(component),
                owners=owners_for_component(component),
                mod14=tuple(owner % 14 for owner in owners_for_component(component)),
                mod27=tuple(owner % 27 for owner in owners_for_component(component)),
            )
        )
    print()

    print("## Forbidden Seam")
    print(f"survivor_gates={len(gates)}")
    print(f"max_delta={max_delta}")
    print(f"seam_gate_count={len(seam_gates)}")
    print(f"seam_components={tuple(sorted(seam_components))}")
    print(f"seam_owner_union={seam_owners}")
    print(f"seam_owner_count={len(seam_owners)}")
    print(f"seam_owner_mod14={seam_owner_mod14}")
    print(f"seam_owner_mod27={seam_owner_mod27}")
    print(f"seam_owner_mod7={seam_owner_mod7}")
    print(f"c27_duplicate_residues={dict(sorted(c27_duplicate_residues.items()))}")
    print(f"seam_endpoint_speeds={seam_endpoint_speeds}")
    print(f"seam_endpoint_residues_mod14={tuple((speed, speed % 14) for speed in seam_endpoint_speeds)}")
    print("seam_gates:")
    for gate in seam_gates:
        print(f"  {gate_word(gate)}")
    print()

    print("## Phase Flow On The Complement")
    print(f"V={phase_row.V}")
    print(f"q=14V={q}")
    print(f"P={phase_metrics.P}")
    print(f"E={phase_metrics.E}")
    print(f"actual_witnesses={len(witnesses)}")
    print(f"phase_counts={sparse(phase_metrics.phase_counts)}")
    print(f"phase_branch_counts={sparse(pullback.phase_branch_counts)}")
    print(f"component_class_hits={dict(sorted(pullback.component_class_counts.items()))}")
    print(f"gate_route_counts={dict(sorted(base_surgery.route_counter.items()))}")
    print(f"gate_total_delta_counts={dict(sorted(base_surgery.delta_counter.items()))}")
    print(f"any_gate_hits={base_surgery.any_gate_hits}")
    print(f"no_gate_hits={base_surgery.no_gate_hits}")
    print(f"ambiguous_gate_hits={base_surgery.ambiguous_gate_hits}")
    print(f"hard_seam_phase_hits={len(seam_hits)}")
    print(f"lower_delta_same_component_bypass_hits={len(bypass_hits)}")
    print(
        "bypass_component_branch_delta_hist="
        f"{dict(sorted(Counter((hit.component, hit.gate_branch, hit.total_delta) for hit in bypass_hits).items()))}"
    )
    print(f"bypass_phase_hist={dict(sorted(Counter(hit.phase for hit in bypass_hits).items()))}")
    print(f"bypass_branch_hist={dict(sorted(Counter(hit.branch for hit in bypass_hits).items()))}")
    print(f"bypass_route_hist={dict(sorted(Counter(hit.route for hit in bypass_hits).items()))}")
    print(f"unique_bypass_gate_count={len(bypass_keys)}")
    print("bypass_gates:")
    gate_by_key = {gate_key(gate): gate for gate in gates}
    for key in sorted(bypass_keys, key=lambda item: (item[0], item[1][0], item[2])):
        print(f"  {gate_word(gate_by_key[key])}")
    print(f"witness_gap_hist_top={witness_gap_hist.most_common(12)}")
    print(f"phase_doubling_edges={dict(sorted(phase_doubling_hist.items()))}")
    print()

    print("## Gate-Surgery Experiments")
    print("seam deletion is the forbidden-seam test: if it is really forbidden,")
    print("phase routing should be unchanged when the max-delta seam is removed.")
    print(
        "  remove_seam: any_gate_hits={a} no_gate_hits={n} ambiguous={m} "
        "delta_any={da} delta_no={dn}".format(
            a=seam_deleted.any_gate_hits,
            n=seam_deleted.no_gate_hits,
            m=seam_deleted.ambiguous_gate_hits,
            da=seam_deleted.any_gate_hits - base_surgery.any_gate_hits,
            dn=seam_deleted.no_gate_hits - base_surgery.no_gate_hits,
        )
    )
    print("bypass deletion is the complement-flow test: the lower-delta mirror")
    print("channel should carry exactly the same-component phase traffic.")
    print(
        "  remove_bypass: any_gate_hits={a} no_gate_hits={n} ambiguous={m} "
        "lost_any={lost} gained_no={gain}".format(
            a=bypass_deleted.any_gate_hits,
            n=bypass_deleted.no_gate_hits,
            m=bypass_deleted.ambiguous_gate_hits,
            lost=base_surgery.any_gate_hits - bypass_deleted.any_gate_hits,
            gain=bypass_deleted.no_gate_hits - base_surgery.no_gate_hits,
        )
    )
    print(
        "  remove_seam_and_bypass: any_gate_hits={a} no_gate_hits={n} ambiguous={m}".format(
            a=seam_and_bypass_deleted.any_gate_hits,
            n=seam_and_bypass_deleted.no_gate_hits,
            m=seam_and_bypass_deleted.ambiguous_gate_hits,
        )
    )
    print()

    print("## n+2 Versus n*2 Recursion Coordinates")
    print(f"n=14 C=2n-1={2 * 14 - 1}")
    print(f"triangular_bridge_check=8*C(n,2)+1={triangular_identity_lhs}")
    print(f"(2n-1)^2={triangular_identity_rhs}")
    print(f"triangular_identity_holds={triangular_identity_lhs == triangular_identity_rhs}")
    print("additive_n_plus_2_readout=dead punctures and C=27 carry seam")
    print("multiplicative_n_times_2_readout=phase flow under u=2t mod 1")
    print("prediction=random031 is the interface defect: the +2/carry face sees")
    print("a seven-owner max-delta seam, while the *2/doubling face routes around")
    print("that seam and only pays a 12-hit lower-delta bypass toll.")
    print()

    print("## Prediction And Experiment Design")
    print("P1 seam-deletion invariance: deleting the two max-delta seam gates")
    print("should leave q=14V phase routing unchanged, as observed above.")
    print("P2 bypass-channel sensitivity: deleting the lower-delta same-component")
    print("bypass gates should move exactly 12 witnesses from gate-hit to no-gate")
    print("unless a wider gate basis adds a secondary channel.")
    print("P3 carry-lift stress: perturb or lift seam owners by C=27 classes;")
    print("preserving the endpoint mod-14 pair (1,9) should keep the seam/bypass")
    print("split more often than preserving raw owner size or raw delta.")
    print("P4 puncture filling: compare HYP-3478 singleton pockets against")
    print("random031's four punctures; isolated dead islands should be harmless")
    print("until a seven-owner seam connects their surrounding flow components.")
    print()

    print("## Tournament Analysis")
    print("vertices=geometric/topological proof carriers, not runners or arcs")
    print("pairwise_observable=puncture isolation + forbidden-seam avoidance + lower-delta bypass + C27 carry/doubling split + scalar firewall")
    print("switch=higher retained LRC predicate and phase-flow payload; ties use declared carrier order")
    print("axes=(branch_predicate,phase_exactness,seam_charge,puncture_topology,bypass_payload,carry_payload,scalar_firewall)")
    print(f"score_hist={hist}")
    print("directed_3cycles=0")
    print("hamiltonian_path=" + " -> ".join(path))
    print()

    print("## Proof Pull")
    print("Treat random031's hard pair as a forbidden seam, not a wall to cross.")
    print("The complement carries all q=14V phase flow; the max-delta seam carries")
    print("the seven-owner gluing charge but receives zero compatible phase hits.")
    print("The only phase contact with the seam components is the lower-delta")
    print("mirror bypass, exactly 6 hits on each side.  A useful named clause")
    print("should therefore glue: mirror-punctured cylinder + seven-owner seam")
    print("+ lower-delta bypass + C=27 carry/doubling split.")


if __name__ == "__main__":
    main()
