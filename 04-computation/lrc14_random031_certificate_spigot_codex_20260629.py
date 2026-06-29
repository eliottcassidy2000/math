#!/usr/bin/env python3
"""HYP-3523: random031 spigot-style terminal certificate stream.

The prompt source is the spigot-algorithm pattern: emit output left-to-right,
discard already-emitted data, and keep only a bounded deferred carry when a
later local event can still alter the current output digit.

This scout applies that pattern to the HYP-3521 random031 terminal certificate
ledger.  Components are read in the q=14V witness order.  Ordinary route and
free-hole single certificates should emit immediately.  Free-hole doublets act
like predigits: one component is buffered until its paired component arrives.
The unique bypass certificate acts like an owner-carry digit: the transport
word emits first, then HYP-3522's branch-boundary lift reduces the owner carry
from four labels to the residual pair (45,173).

Tournament Analysis declaration:
  vertices: proof-output states and carry buffers, not runners, arcs, or raw
            witness counts;
  pairwise observable: no-backtracking terminal emission with named carry;
  switch/gauge: higher emission safety and lower unresolved carry wins;
  tie Hamiltonian path: spigot state before specific carry buffers, then raw
            static shadows.

The output is evidence for a proof scheduler, not an LRC14 proof.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys
from typing import Callable


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


H3521 = load_module(
    "hyp3521_for_hyp3523",
    "lrc14_random031_terminal_certificate_ledger_codex_20260629.py",
)
H3522 = load_module(
    "hyp3522_for_hyp3523",
    "lrc14_random031_owner_boundary_filtration_codex_20260629.py",
)


@dataclass(frozen=True)
class ComponentWindow:
    rank: int
    component_index: int
    terminal_class: str
    component_type: str
    packet_ids: tuple[int, ...]
    cluster_id: tuple[int, ...] | None
    cell_count: int
    first_a: int
    last_a: int
    first_u: int
    last_u: int
    owner_union: tuple[int, ...]


@dataclass(frozen=True)
class BoundaryInfo:
    boundary_component_indices: tuple[int, ...]
    seam_owners: tuple[int, ...]
    transport_owners: tuple[int, ...]
    branch_boundary_owners: tuple[int, ...]
    bracket_lift_owners: tuple[int, ...]
    residual_after_transport: tuple[int, ...]
    residual_after_branch_boundary: tuple[int, ...]


@dataclass(frozen=True)
class StreamRecord:
    rank: int
    component_index: int
    terminal_class: str
    component_type: str
    cluster_id: tuple[int, ...] | None
    cell_count: int
    first_a: int
    last_a: int
    actions: tuple[str, ...]
    action_kind: str
    doublet_state: str
    owner_state: str
    pending_doublet_components: int
    owner_carry: tuple[int, ...]
    emitted_certificates: int


@dataclass(frozen=True)
class StreamResult:
    records: tuple[StreamRecord, ...]
    emitted_certificates: int
    doublet_events: tuple[tuple[tuple[int, ...], int, int, tuple[int, ...]], ...]
    owner_events: tuple[tuple[int, str, tuple[int, ...]], ...]
    max_pending_doublet_components: int
    max_owner_carry_width: int
    final_owner_carry: tuple[int, ...]


@dataclass(frozen=True)
class Carrier:
    name: str
    score: int


CARRIERS = (
    Carrier("spigot_state_plus_route_R", 104),
    Carrier("terminal_class_plus_carry_state", 98),
    Carrier("owner_carry_buffer", 92),
    Carrier("doublet_predigit_buffer", 86),
    Carrier("terminal_certificate_ledger", 80),
    Carrier("terminal_class_shadow", 63),
    Carrier("component_type_shadow", 41),
    Carrier("raw_cell_count_shadow", 10),
)


def compact_counter(counter: Counter | dict) -> dict:
    def key(item):
        value = item[0]
        if isinstance(value, int):
            return (0, value)
        return (1, repr(value))

    return dict(sorted(counter.items(), key=key))


def component_windows() -> tuple[
    tuple[ComponentWindow, ...],
    tuple[object, ...],
    tuple[object, ...],
    dict[tuple[int, int], object],
]:
    _row, _cells, components, by_node = H3521.build_cells_and_components()
    certificates = H3521.build_certificates()
    if len(components) != len(certificates):
        raise RuntimeError("component/certificate length mismatch")

    windows: list[ComponentWindow] = []
    for index, (component, certificate) in enumerate(zip(components, certificates)):
        a_values = sorted(by_node[node].a for node in component.nodes)
        u_values = sorted(node[0] for node in component.nodes)
        windows.append(
            ComponentWindow(
                rank=-1,
                component_index=index,
                terminal_class=certificate.terminal_class,
                component_type=certificate.component_type,
                packet_ids=certificate.packet_ids,
                cluster_id=certificate.cluster_id,
                cell_count=certificate.cell_count,
                first_a=a_values[0],
                last_a=a_values[-1],
                first_u=u_values[0],
                last_u=u_values[-1],
                owner_union=certificate.owner_union,
            )
        )

    ordered = sorted(
        windows,
        key=lambda item: (item.first_a, item.last_a, item.component_index),
    )
    ranked = tuple(
        ComponentWindow(
            rank=rank,
            component_index=item.component_index,
            terminal_class=item.terminal_class,
            component_type=item.component_type,
            packet_ids=item.packet_ids,
            cluster_id=item.cluster_id,
            cell_count=item.cell_count,
            first_a=item.first_a,
            last_a=item.last_a,
            first_u=item.first_u,
            last_u=item.last_u,
            owner_union=item.owner_union,
        )
        for rank, item in enumerate(ordered)
    )
    return ranked, components, certificates, by_node


def boundary_info(components: tuple[object, ...]) -> BoundaryInfo:
    (
        _row,
        _gates,
        _cells,
        by_node,
        _legal,
        bypass_component,
        bypass_cells,
        hard_seam_gates,
        _lower_bypass_gates,
    ) = H3522.build_context()

    node_to_component = {
        node: index for index, component in enumerate(components) for node in component.nodes
    }
    boundaries = H3522.branch_boundaries(_cells, by_node, bypass_component)
    boundary_cells = tuple(
        cell
        for boundary in boundaries
        for cell in (boundary.left_cell, boundary.right_cell)
    )
    boundary_component_indices = tuple(
        sorted({node_to_component[cell.node] for cell in boundary_cells})
    )

    seam_owners = H3522.gate_owner_union(hard_seam_gates)
    transport_owners = H3522.owner_union_from_cells(bypass_cells)
    branch_boundary_owners = H3522.owner_union_from_cells(boundary_cells)
    transport_plus_boundary = tuple(
        sorted(set(transport_owners) | set(branch_boundary_owners))
    )
    residual_after_transport = tuple(
        owner for owner in seam_owners if owner not in transport_owners
    )
    bracket_lift_owners = tuple(
        owner for owner in branch_boundary_owners if owner in residual_after_transport
    )
    residual_after_branch_boundary = tuple(
        owner for owner in seam_owners if owner not in transport_plus_boundary
    )

    return BoundaryInfo(
        boundary_component_indices=boundary_component_indices,
        seam_owners=seam_owners,
        transport_owners=transport_owners,
        branch_boundary_owners=branch_boundary_owners,
        bracket_lift_owners=bracket_lift_owners,
        residual_after_transport=residual_after_transport,
        residual_after_branch_boundary=residual_after_branch_boundary,
    )


def simulate_stream(
    windows: tuple[ComponentWindow, ...],
    certificates: tuple[object, ...],
    boundary: BoundaryInfo,
) -> StreamResult:
    cluster_targets = Counter(
        certificate.cluster_id
        for certificate in certificates
        if certificate.terminal_class == "free_hole_doublet_packet"
    )
    pending_doublets: dict[tuple[int, ...], list[int]] = defaultdict(list)
    seen_components: set[int] = set()
    owner_carry: tuple[int, ...] = tuple()
    bracket_lift_applied = False
    emitted_total = 0
    max_pending_doublet_components = 0
    max_owner_carry_width = 0
    doublet_events: list[tuple[tuple[int, ...], int, int, tuple[int, ...]]] = []
    owner_events: list[tuple[int, str, tuple[int, ...]]] = []
    records: list[StreamRecord] = []

    boundary_component_set = set(boundary.boundary_component_indices)

    for window in windows:
        seen_components.add(window.component_index)
        actions: list[str] = []
        emitted = 0
        doublet_state = "none"
        owner_state = "none"

        if window.terminal_class == "ordinary_rank2_route_component":
            actions.append("emit_ordinary_route")
            emitted += 1
        elif window.terminal_class == "free_hole_single_bracket_packet":
            actions.append("emit_free_hole_single")
            emitted += 1
        elif window.terminal_class == "free_hole_doublet_packet":
            if window.cluster_id is None:
                raise RuntimeError("doublet packet missing cluster id")
            pending_doublets[window.cluster_id].append(window.component_index)
            if len(pending_doublets[window.cluster_id]) == cluster_targets[window.cluster_id]:
                members = tuple(pending_doublets.pop(window.cluster_id))
                actions.append("emit_free_hole_doublet")
                emitted += 1
                doublet_state = "emit_second"
                first_rank = next(
                    record.rank
                    for record in records
                    if record.component_index in members
                )
                doublet_events.append(
                    (window.cluster_id, first_rank, window.rank, members)
                )
            else:
                actions.append("buffer_free_hole_predigit")
                doublet_state = "buffer_first"
        elif window.terminal_class == "pure_bypass_owner_boundary_component":
            actions.append("emit_bypass_transport_predigit")
            owner_carry = boundary.residual_after_transport
            owner_state = "open_owner_carry"
            owner_events.append((window.rank, owner_state, owner_carry))
        else:
            actions.append("named_debt")

        if (
            owner_carry
            and not bracket_lift_applied
            and boundary_component_set.issubset(seen_components)
        ):
            owner_carry = tuple(
                owner
                for owner in owner_carry
                if owner not in boundary.bracket_lift_owners
            )
            bracket_lift_applied = True
            actions.append("apply_branch_boundary_lift")
            actions.append("emit_bypass_owner_boundary_certificate")
            emitted += 1
            owner_state = (
                "open_owner_carry+apply_branch_boundary_lift"
                if owner_state == "open_owner_carry"
                else "apply_branch_boundary_lift"
            )
            owner_events.append((window.rank, "apply_branch_boundary_lift", owner_carry))

        pending_component_count = sum(len(items) for items in pending_doublets.values())
        max_pending_doublet_components = max(
            max_pending_doublet_components,
            pending_component_count,
        )
        max_owner_carry_width = max(max_owner_carry_width, len(owner_carry))
        emitted_total += emitted

        records.append(
            StreamRecord(
                rank=window.rank,
                component_index=window.component_index,
                terminal_class=window.terminal_class,
                component_type=window.component_type,
                cluster_id=window.cluster_id,
                cell_count=window.cell_count,
                first_a=window.first_a,
                last_a=window.last_a,
                actions=tuple(actions),
                action_kind="+".join(actions),
                doublet_state=doublet_state,
                owner_state=owner_state,
                pending_doublet_components=pending_component_count,
                owner_carry=owner_carry,
                emitted_certificates=emitted,
            )
        )

    if pending_doublets:
        raise RuntimeError(f"unclosed doublet predigits: {dict(pending_doublets)}")
    if emitted_total != H3521.terminal_certificate_count(certificates):
        raise RuntimeError(
            f"emitted {emitted_total}, expected {H3521.terminal_certificate_count(certificates)}"
        )

    return StreamResult(
        records=tuple(records),
        emitted_certificates=emitted_total,
        doublet_events=tuple(doublet_events),
        owner_events=tuple(owner_events),
        max_pending_doublet_components=max_pending_doublet_components,
        max_owner_carry_width=max_owner_carry_width,
        final_owner_carry=owner_carry,
    )


def quotient_audit(
    records: tuple[StreamRecord, ...],
    axis_name: str,
    axis: Callable[[StreamRecord], object],
) -> dict[str, object]:
    fibers: dict[object, set[str]] = defaultdict(set)
    sizes: Counter[object] = Counter()
    for record in records:
        key = axis(record)
        fibers[key].add(record.action_kind)
        sizes[key] += 1
    mixed = {key: values for key, values in fibers.items() if len(values) > 1}
    largest_mixed = sorted(
        ((sizes[key], key, tuple(sorted(values))) for key, values in mixed.items()),
        reverse=True,
        key=lambda item: (item[0], repr(item[1])),
    )
    return {
        "axis": axis_name,
        "fibers": len(fibers),
        "mixed_fibers": len(mixed),
        "mixed_rows": sum(sizes[key] for key in mixed),
        "largest_mixed": largest_mixed[:3],
    }


def tournament_fingerprint() -> dict[str, object]:
    scores = Counter(carrier.score for carrier in CARRIERS)
    path = tuple(carrier.name for carrier in sorted(CARRIERS, key=lambda item: -item.score))
    return {
        "score_hist": compact_counter(scores),
        "directed_3cycles": 0,
        "hamiltonian_path": path,
    }


def main() -> None:
    windows, components, certificates, _by_node = component_windows()
    boundary = boundary_info(components)
    stream = simulate_stream(windows, certificates, boundary)

    component_type_hist = Counter(window.component_type for window in windows)
    terminal_class_hist = Counter(window.terminal_class for window in windows)
    action_hist = Counter(record.action_kind for record in stream.records)
    emitted_hist = Counter(record.emitted_certificates for record in stream.records)
    pending_width_hist = Counter(record.pending_doublet_components for record in stream.records)
    owner_width_hist = Counter(len(record.owner_carry) for record in stream.records)

    doublet_rank_widths = {
        cluster: end_rank - start_rank
        for cluster, start_rank, end_rank, _members in stream.doublet_events
    }
    owner_carry_rank_width = (
        stream.owner_events[-1][0] - stream.owner_events[0][0]
        if len(stream.owner_events) >= 2
        else 0
    )
    boundary_ranks = tuple(
        window.rank
        for window in windows
        if window.component_index in boundary.boundary_component_indices
    )
    bypass_rank = next(
        window.rank
        for window in windows
        if window.terminal_class == "pure_bypass_owner_boundary_component"
    )

    audits = (
        quotient_audit(stream.records, "component_type", lambda record: record.component_type),
        quotient_audit(stream.records, "terminal_class", lambda record: record.terminal_class),
        quotient_audit(
            stream.records,
            "terminal_class_plus_cluster",
            lambda record: (record.terminal_class, record.cluster_id),
        ),
        quotient_audit(
            stream.records,
            "terminal_class_plus_spigot_state",
            lambda record: (
                record.terminal_class,
                record.doublet_state,
                record.owner_state,
            ),
        ),
        quotient_audit(
            stream.records,
            "component_index",
            lambda record: record.component_index,
        ),
    )
    fp = tournament_fingerprint()

    print("HYP-3523 RANDOM031 SPIGOT-STYLE TERMINAL STREAM AUDIT")
    print("status=EVIDENCE / finite streaming certificate scheduler; not an LRC14 proof")
    print("row=random_covering_031")
    print("source_inspiration=https://en.wikipedia.org/wiki/Spigot_algorithm")
    print()

    print("## Spigot Transfer")
    print("spigot_rule=emit left-to-right, discard emitted certificates, retain only bounded carry")
    print("digit_tokens=HYP-3521 terminal certificates, not runners or raw row names")
    print("predigit_buffer=HYP-3511 free-hole doublet first components")
    print("carry_buffer=HYP-3522 bypass owner residual after transport")
    print("route_guardrail=HYP-3513 says random031 keeps route sidecar R unless route reconstruction is proved")
    print()

    print("## Stream Summary")
    print(f"component_tokens={len(windows)}")
    print(f"terminal_certificate_count={H3521.terminal_certificate_count(certificates)}")
    print(f"emitted_certificate_count={stream.emitted_certificates}")
    print(f"component_type_hist={compact_counter(component_type_hist)}")
    print(f"terminal_class_hist={compact_counter(terminal_class_hist)}")
    print(f"action_hist={compact_counter(action_hist)}")
    print(f"emitted_per_component_hist={compact_counter(emitted_hist)}")
    print(f"pending_doublet_width_hist={compact_counter(pending_width_hist)}")
    print(f"owner_carry_width_hist={compact_counter(owner_width_hist)}")
    print(f"max_pending_doublet_components={stream.max_pending_doublet_components}")
    print(f"max_owner_carry_width={stream.max_owner_carry_width}")
    print(f"final_owner_carry={stream.final_owner_carry}")
    print()

    print("## Doublet Predigit Buffer")
    print(f"doublet_cluster_count={len(stream.doublet_events)}")
    print(f"doublet_pending_rank_widths={doublet_rank_widths}")
    for cluster, start_rank, end_rank, members in stream.doublet_events:
        print(
            "  doublet cluster={cluster} first_rank={start} emit_rank={end} "
            "rank_gap={gap} component_indices={members}".format(
                cluster=cluster,
                start=start_rank,
                end=end_rank,
                gap=end_rank - start_rank,
                members=members,
            )
        )
    print("reading=free-hole doublets need a one-component predigit buffer and close at the next component event")
    print()

    print("## Owner Carry Buffer")
    print(f"seam_owners={boundary.seam_owners}")
    print(f"transport_owners={boundary.transport_owners}")
    print(f"branch_boundary_owners={boundary.branch_boundary_owners}")
    print(f"boundary_component_indices={boundary.boundary_component_indices}")
    print(f"boundary_component_ranks={boundary_ranks}")
    print(f"bypass_component_rank={bypass_rank}")
    print(f"residual_after_transport={boundary.residual_after_transport}")
    print(f"bracket_lift_owners={boundary.bracket_lift_owners}")
    print(f"residual_after_branch_boundary={boundary.residual_after_branch_boundary}")
    print(f"owner_carry_rank_width={owner_carry_rank_width}")
    for rank, event, carry in stream.owner_events:
        print(f"  owner_event rank={rank} event={event} carry={carry}")
    print("reading=the bypass owner carry opens as four labels and is reduced one component event later to (45,173)")
    print()

    print("## Key Stream Events")
    for record in stream.records:
        if (
            record.doublet_state != "none"
            or record.owner_state != "none"
            or "pure_bypass" in record.terminal_class
        ):
            print(
                "  rank={rank} comp={comp} a=[{a0},{a1}] class={cls} "
                "cluster={cluster} actions={actions} pending_doublets={pending} owner_carry={carry}".format(
                    rank=record.rank,
                    comp=record.component_index,
                    a0=record.first_a,
                    a1=record.last_a,
                    cls=record.terminal_class,
                    cluster=record.cluster_id,
                    actions=record.actions,
                    pending=record.pending_doublet_components,
                    carry=record.owner_carry,
                )
            )
    print()

    print("## Quotient Audit")
    print("target=stream_action_kind")
    for audit in audits:
        print(
            "axis={axis} fibers={fibers} mixed_fibers={mixed_fibers} mixed_rows={mixed_rows}".format(
                **audit
            )
        )
        if audit["largest_mixed"]:
            print(f"  largest_mixed={audit['largest_mixed']}")
    print("reading=terminal class alone forgets predigit and owner-carry state; adding spigot state restores purity")
    print()

    print("## Tournament Analysis")
    print("vertices=stream states and proof-carrier buffers, not runners or raw arcs")
    print("pairwise_observable=streamability + bounded carry + terminal predicate + route/firewall sidecar")
    print("switch=higher proof-facing stream score; ties use the fixed carrier order")
    print(f"score_hist={fp['score_hist']}")
    print(f"directed_3cycles={fp['directed_3cycles']}")
    print("hamiltonian_path=" + " -> ".join(fp["hamiltonian_path"]))
    print()

    print("## Proof Pull")
    print("P1: Treat the 77 random031 terminal certificates as a left-to-right stream over component events.")
    print("P2: Ordinary route and free-hole single certificates are immediate and need no row-name memory.")
    print("P3: Free-hole doublets need only a one-component predigit buffer, closing at the next component event.")
    print("P4: The bypass opens a four-owner carry (45,147,169,173) and HYP-3522 branch bracketing reduces it after one event to the residual pair (45,173).")
    print("P5: The remaining non-streamed proof object is exactly residual owners (45,173) plus HYP-3513 route sidecar R, not the full seven-owner seam.")
    print()

    print("## Assumption Challenge")
    print("candidate_vertices=runners,gaps,fixed_sections,section_boundaries,wall_crossings,residues,cover_arcs,Fourier_modes,matroid_circuits,proof_obligations,stream_states")
    print("chosen_vertices=terminal stream states plus predigit/carry buffers")
    print("preserves=random031 terminal-discharge predicate, free-hole doublet collapse, owner residual, and firewall route sidecar")
    print("destroys=raw runner order, raw component index after emission, and scalar cell-count shadows")
    print("challenged_assumption=the terminal proof must keep the whole random031 component table live until the end; the audit shows all but a bounded carry can stream away")


if __name__ == "__main__":
    main()
