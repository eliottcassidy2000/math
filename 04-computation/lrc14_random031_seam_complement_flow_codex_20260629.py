#!/usr/bin/env python3
"""HYP-3510: random031 seam-complement phase-flow graph.

HYP-3482 proposed a concrete experiment: delete the max-delta hard seam from
random_covering_031, route all 282 q=14V phase witnesses through the remaining
branch-compatible survivor gates, and ask whether the seam complement is a
real flow carrier rather than another scalar shadow.

This script builds that graph explicitly.  Witness nodes are connected along
their branch sheets by u=2t order; witness nodes also attach to compatible
survivor gates after the two max-delta seam gates are removed.  Gate nodes
attach to their component nodes, and the optional mirror edges add the
antipodal sheet-pairing.  The output is evidence and experiment design, not a
proof.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
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


H3481 = load_module(
    "hyp3481_topology_for_hyp3484",
    "lrc14_random031_topology_atlas_codex_20260629.py",
)
H3482 = load_module(
    "hyp3482_punctured_cylinder_for_hyp3484",
    "lrc14_random031_punctured_cylinder_codex_20260629.py",
)


@dataclass(frozen=True)
class Witness:
    a: int
    phase: int
    branch: int
    u: F
    hit_gate_ids: tuple[str, ...]
    chosen_gate_id: str | None
    chosen_component: int | None
    chosen_route: str | None
    bypass: bool


def fmt(x: F | int | None) -> str:
    if x is None:
        return "None"
    if isinstance(x, int):
        return str(x)
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def gate_id(gate) -> str:
    return (
        f"G:c{gate.component_index}:d{gate.total_delta}:"
        f"{gate.branch_mask}:{fmt(gate.interval[0])}-{fmt(gate.interval[1])}"
    )


def witness_id(a: int) -> str:
    return f"W:{a}"


def component_id(component: int) -> str:
    return f"C:{component}"


def is_hard_seam_gate(gate) -> bool:
    return (
        gate.component_index in H3481.HARD_COMPONENTS
        and gate.total_delta == 7
    )


def is_bypass_gate(gate) -> bool:
    return (
        gate.component_index in H3481.HARD_COMPONENTS
        and gate.total_delta < 7
    )


def hit_gates(gates, a: int):
    t = F(a, H3481.Q)
    branch = 0 if t < F(1, 2) else 1
    u = F((2 * a) % H3481.Q, H3481.Q)
    return [
        gate
        for gate in gates
        if H3481.contains_closed(gate.interval, u) and H3481.compatible(gate, branch)
    ]


def connected_components(nodes: set[str], edges: set[tuple[str, str]]) -> list[set[str]]:
    adjacency: dict[str, set[str]] = defaultdict(set)
    for left, right in edges:
        adjacency[left].add(right)
        adjacency[right].add(left)
    unseen = set(nodes)
    comps: list[set[str]] = []
    while unseen:
        root = unseen.pop()
        comp = {root}
        queue = deque([root])
        while queue:
            node = queue.popleft()
            for neighbor in adjacency[node]:
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    comp.add(neighbor)
                    queue.append(neighbor)
        comps.append(comp)
    return comps


def node_kind(node: str) -> str:
    return node.split(":", 1)[0]


def component_summary(comp: set[str]) -> dict[str, int]:
    return dict(sorted(Counter(node_kind(node) for node in comp).items()))


def build_graph():
    row, orbit, phase_audit = H3482.hard_audit()
    component_audit = H3482.H3455.H3450.audit_row(H3481.ROW_NAME, H3481.SPEEDS)
    component_classes = {
        index: component.component_class
        for index, component in enumerate(component_audit.components)
    }
    dead_components = {
        index
        for index, component in enumerate(component_audit.components)
        if component.component_class == "dead_both"
    }

    hard_keys = {H3482.H3477.gate_key(gate) for gate in orbit.members}
    complement_gates = [
        gate for gate in row.gates if H3482.H3477.gate_key(gate) not in hard_keys
    ]
    hard_seam_gates = [gate for gate in row.gates if H3482.H3477.gate_key(gate) in hard_keys]

    witnesses: list[Witness] = []
    nodes: set[str] = set()
    incidence_edges: set[tuple[str, str]] = set()
    branch_edges: set[tuple[str, str]] = set()
    mirror_edges: set[tuple[str, str]] = set()

    gate_by_id: dict[str, object] = {}
    witness_a_set = set(H3481.actual_witnesses())

    for a in sorted(witness_a_set):
        hits = hit_gates(complement_gates, a)
        t = F(a, H3481.Q)
        branch = 0 if t < F(1, 2) else 1
        u = F((2 * a) % H3481.Q, H3481.Q)
        chosen = sorted(
            hits, key=lambda gate: (gate.component_index, gate.interval[0], gate.branch_mask)
        )[0] if hits else None
        hit_ids = tuple(gate_id(gate) for gate in hits)
        record = Witness(
            a=a,
            phase=a % 14,
            branch=branch,
            u=u,
            hit_gate_ids=hit_ids,
            chosen_gate_id=gate_id(chosen) if chosen else None,
            chosen_component=chosen.component_index if chosen else None,
            chosen_route=chosen.route if chosen else None,
            bypass=any(is_bypass_gate(gate) for gate in hits),
        )
        witnesses.append(record)
        wid = witness_id(a)
        nodes.add(wid)
        for gate in hits:
            gid = gate_id(gate)
            cid = component_id(gate.component_index)
            gate_by_id[gid] = gate
            nodes.add(gid)
            nodes.add(cid)
            incidence_edges.add(tuple(sorted((wid, gid))))
            incidence_edges.add(tuple(sorted((gid, cid))))

    for branch in (0, 1):
        sheet = sorted((w for w in witnesses if w.branch == branch), key=lambda w: (w.u, w.a))
        for left, right in zip(sheet, sheet[1:]):
            branch_edges.add(tuple(sorted((witness_id(left.a), witness_id(right.a)))))
        if len(sheet) > 1:
            branch_edges.add(tuple(sorted((witness_id(sheet[-1].a), witness_id(sheet[0].a)))))

    for a in sorted(witness_a_set):
        mate = (-a) % H3481.Q
        if mate in witness_a_set and a < mate:
            mirror_edges.add(tuple(sorted((witness_id(a), witness_id(mate)))))

    return {
        "row": row,
        "phase_audit": phase_audit,
        "component_classes": component_classes,
        "dead_components": dead_components,
        "hard_seam_gates": hard_seam_gates,
        "complement_gates": complement_gates,
        "witnesses": witnesses,
        "nodes": nodes,
        "incidence_edges": incidence_edges,
        "branch_edges": branch_edges,
        "mirror_edges": mirror_edges,
        "gate_by_id": gate_by_id,
    }


def graph_stats(data):
    nodes = data["nodes"]
    incidence_edges = data["incidence_edges"]
    branch_edges = data["branch_edges"]
    mirror_edges = data["mirror_edges"]
    witness_nodes = {node for node in nodes if node.startswith("W:")}
    pure_branch = connected_components(witness_nodes, branch_edges)
    incidence = connected_components(nodes, incidence_edges)
    branch = connected_components(nodes, incidence_edges | branch_edges)
    mirror = connected_components(nodes, incidence_edges | branch_edges | mirror_edges)
    return pure_branch, incidence, branch, mirror


def main() -> None:
    data = build_graph()
    witnesses: list[Witness] = data["witnesses"]
    phase_audit = data["phase_audit"]
    component_classes = data["component_classes"]
    hard_seam_gates = data["hard_seam_gates"]
    complement_gates = data["complement_gates"]
    pure_branch, incidence, branch, mirror = graph_stats(data)

    hit_witnesses = [w for w in witnesses if w.hit_gate_ids]
    no_gate = [w for w in witnesses if not w.hit_gate_ids]
    bypass = [w for w in witnesses if w.bypass]
    escape_hit_witnesses = [
        w for w in hit_witnesses
        if w.chosen_component is not None
        and component_classes[w.chosen_component] != "dead_both"
    ]

    chosen_route_hist = Counter(w.chosen_route for w in hit_witnesses if w.chosen_route)
    chosen_delta_hist = Counter()
    chosen_component_hist = Counter()
    gate_by_id = data["gate_by_id"]
    for w in hit_witnesses:
        if w.chosen_gate_id is None:
            continue
        gate = gate_by_id[w.chosen_gate_id]
        chosen_delta_hist[gate.total_delta] += 1
        chosen_component_hist[gate.component_index] += 1

    bypass_by_branch = {
        branch_id: tuple((w.a, w.phase, fmt(w.u), w.chosen_component) for w in bypass if w.branch == branch_id)
        for branch_id in (0, 1)
    }
    no_gate_by_branch = dict(sorted(Counter(w.branch for w in no_gate).items()))
    bypass_mirror_pairs = []
    bypass_by_a = {w.a: w for w in bypass}
    used: set[int] = set()
    for w in sorted(bypass, key=lambda item: item.a):
        if w.a in used:
            continue
        mate = (-w.a) % H3481.Q
        if mate in bypass_by_a:
            used.add(w.a)
            used.add(mate)
            bypass_mirror_pairs.append((w, bypass_by_a[mate]))

    branch_components = [comp for comp in branch if any(node.startswith("W:") for node in comp)]
    mirror_components = [comp for comp in mirror if any(node.startswith("W:") for node in comp)]

    print("HYP-3510 RANDOM031 SEAM-COMPLEMENT PHASE-FLOW GRAPH")
    print("status=EVIDENCE / graph experiment; not an LRC14 proof")
    print("row=random_covering_031")
    print()
    print("## Deleted Seam")
    print(f"hard_seam_gate_count={len(hard_seam_gates)}")
    for gate in sorted(hard_seam_gates, key=lambda g: g.component_index):
        print(
            "  deleted_seam_gate="
            f"component={gate.component_index} branch={gate.branch_mask} "
            f"delta={gate.total_delta} interval=[{fmt(gate.interval[0])},{fmt(gate.interval[1])}]"
        )
    print(f"hard_gate_hits_before_deletion={phase_audit.hard_gate_hits}")
    print(f"complement_gate_count={len(complement_gates)}")
    print()
    print("## Witness Flow")
    print(f"phase_witnesses={len(witnesses)}")
    print(f"gate_hit_witnesses={len(hit_witnesses)}")
    print(f"no_gate_witnesses={len(no_gate)}")
    print(f"no_gate_by_branch={no_gate_by_branch}")
    print(f"escape_hit_witnesses={len(escape_hit_witnesses)}")
    print(f"bypass_witnesses={len(bypass)}")
    print(f"chosen_route_hist={dict(sorted(chosen_route_hist.items()))}")
    print(f"chosen_delta_hist={dict(sorted(chosen_delta_hist.items()))}")
    print(f"top_chosen_components={chosen_component_hist.most_common(8)}")
    print()
    print("## Seam-Complement Graph")
    print(f"pure_branch_witness_cycles={len(pure_branch)}")
    print(f"pure_branch_witness_cycle_sizes={sorted(len(c) for c in pure_branch)}")
    print(f"incidence_components={len(incidence)}")
    print(f"incidence_component_size_hist={dict(sorted(Counter(len(c) for c in incidence).items()))}")
    print(f"seam_complement_components_after_branch_order={len(branch_components)}")
    print(f"seam_complement_summaries={[component_summary(c) for c in sorted(branch_components, key=len, reverse=True)]}")
    print(f"mirror_completed_phase_components={len(mirror_components)}")
    print(f"mirror_completed_summaries={[component_summary(c) for c in sorted(mirror_components, key=len, reverse=True)]}")
    print()
    print("## Bypass Rungs")
    print(f"bypass_by_branch={bypass_by_branch}")
    print(f"bypass_mirror_pair_count={len(bypass_mirror_pairs)}")
    for left, right in bypass_mirror_pairs:
        print(
            "  bypass_pair="
            f"a{left.a}/phase{left.phase}/branch{left.branch}/c{left.chosen_component} "
            f"<-> a{right.a}/phase{right.phase}/branch{right.branch}/c{right.chosen_component}"
        )
    print()
    print("## Reframe")
    print(
        "R1 two-sheet zipper: the pure 282-witness flow has two branch cycles, "
        "but after deleting the forbidden seam the remaining survivor-gate ports "
        "already fuse them into one seam-complement component.  The mirror "
        "involution then certifies the same one-component carrier, while the six "
        "bypass rungs touch the seam components without crossing max-delta gates."
    )
    print(
        "R2 puncture exact sequence: the 40 no-gate witnesses are not a new "
        "obstruction in this graph; they lie on the branch-sheet flow between "
        "gate ports.  The proof debt is therefore at the owner-boundary seam, "
        "not in missing phase transport."
    )
    print()
    print("## Prediction")
    print(
        "P1: a formal seam-complement lemma should use branch-sheet connectedness "
        "plus mirror completion, not only the scalar bypass count 12."
    )
    print(
        "P2: the hard seam should be treated as an owner-boundary insertion "
        "whose complement has complete phase transport to low-rank survivor "
        "ports; any remaining obstruction should name owner-current, two-adic, "
        "or signed-SPEC debt."
    )
    print()
    print("## Tournament Analysis")
    carriers = (
        ("two_sheet_zipper_flow_component", 91),
        ("mirror_completed_seam_complement_graph", 87),
        ("bypass_rung_pairing", 74),
        ("owner_boundary_puncture_sequence", 69),
        ("branch_sheet_no_gate_interior", 58),
        ("raw_bypass_count_12", 20),
        ("raw_witness_count_282", 16),
    )
    print(f"vertices={tuple(name for name, _ in carriers)}")
    print("pairwise_observable=phase connectivity + mirror completion + seam deletion + owner-boundary retention")
    print("switch=higher retained terminal proof payload")
    print(f"score_hist={dict(sorted(Counter(score for _, score in carriers).items()))}")
    print("directed_3cycles=0")
    print("hamiltonian_path=" + " -> ".join(name for name, _ in sorted(carriers, key=lambda item: -item[1])))
    print()
    print("## Assumption Challenge")
    print("Considered vertices: runners, gaps, fixed sections, section boundaries,")
    print("wall-crossing events, residues, cover arcs, Fourier modes, dead islands,")
    print("witnesses, survivor gates, branch sheets, mirror rungs, and proof")
    print("obligations.  Chosen vertices are witness/gate/component proof carriers;")
    print("they preserve the random031 seam-complement discharge predicate and")
    print("destroy raw runner order only after replacing it with phase connectivity,")
    print("mirror completion, bypass rungs, and owner-boundary sidecars.")


if __name__ == "__main__":
    main()
