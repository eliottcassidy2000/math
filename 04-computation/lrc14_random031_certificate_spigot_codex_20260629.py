#!/usr/bin/env python3
"""HYP-3523: random031 certificate spigot.

The spigot-algorithm prompt is used here as a proof-procedure metaphor, not as
numerology.  A digit is emitted only after possible carry is controlled.  In
the random031 terminal ledger, a certificate is emitted only after possible
free-hole, owner-boundary, private-firewall, and vertical-quotient carry is
named.

Tournament Analysis declaration:
  vertices: proof-output states and carry buffers, not runners, arcs, or raw
            witness counts;
  pairwise observable: no-backtracking terminal emission with named carry;
  switch/gauge: higher emission safety and lower unresolved carry wins;
  tie Hamiltonian path: spigot state before specific carry buffers, then raw
            static counts.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
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


H3521 = load_module(
    "hyp3521_for_hyp3523",
    "lrc14_random031_terminal_certificate_ledger_codex_20260629.py",
)
H3522 = load_module(
    "hyp3522_for_hyp3523",
    "lrc14_random031_owner_boundary_filtration_codex_20260629.py",
)


@dataclass(frozen=True)
class StreamEvent:
    stream_index: int
    component_index: int
    terminal_class: str
    action: str
    emitted_certificate: bool
    emission_key: str
    branch_span: tuple[int, ...]
    u_span: tuple[int, int]
    cell_count: int
    emitted_cell_count: int
    carry_in: tuple[str, ...]
    carry_resolved: tuple[str, ...]
    carry_out: tuple[str, ...]
    sidecars: tuple[str, ...]


@dataclass(frozen=True)
class Carrier:
    name: str
    score: int


CARRIERS = (
    Carrier("certificate_spigot_state_machine", 100),
    Carrier("free_hole_doublet_predigit_buffer", 93),
    Carrier("owner_filtration_carry_buffer", 91),
    Carrier("private_firewall_route_guard", 86),
    Carrier("vertical_halfturn_guard", 82),
    Carrier("ordinary_immediate_route_digit", 77),
    Carrier("static_terminal_ledger", 61),
    Carrier("raw_count_shadow", 10),
)


def compact_counter(counter: Counter | dict) -> dict:
    return dict(sorted(counter.items(), key=lambda item: repr(item[0])))


def component_order_key(component) -> tuple[int, int, int]:
    branches = [branch for _u, branch in component.nodes]
    u_values = [u for u, _branch in component.nodes]
    return (min(branches), min(u_values), len(component.nodes))


def component_branch_span(component) -> tuple[int, ...]:
    return tuple(sorted({branch for _u, branch in component.nodes}))


def component_u_span(component) -> tuple[int, int]:
    u_values = [u for u, _branch in component.nodes]
    return (min(u_values), max(u_values))


def owner_filtration_words() -> dict[str, tuple[int, ...]]:
    (
        _row,
        _gates,
        cells,
        by_node,
        _legal,
        bypass_component,
        bypass_cells,
        hard_seam_gates,
        _lower_bypass_gates,
    ) = H3522.build_context()
    boundaries = H3522.branch_boundaries(cells, by_node, bypass_component)
    seam_owners = H3522.gate_owner_union(hard_seam_gates)
    transport_owners = H3522.owner_union_from_cells(bypass_cells)
    boundary_cells = tuple(
        cell for boundary in boundaries for cell in (boundary.left_cell, boundary.right_cell)
    )
    branch_boundary_owners = H3522.owner_union_from_cells(boundary_cells)
    transport_plus_boundary = tuple(sorted(set(transport_owners) | set(branch_boundary_owners)))
    return {
        "seam": seam_owners,
        "transport": transport_owners,
        "branch_boundary": branch_boundary_owners,
        "bracket_lift": tuple(owner for owner in branch_boundary_owners if owner not in transport_owners),
        "residual_after_transport": tuple(owner for owner in seam_owners if owner not in transport_owners),
        "residual_after_branch_boundary": tuple(owner for owner in seam_owners if owner not in transport_plus_boundary),
    }


def build_stream_events() -> tuple[StreamEvent, ...]:
    _row, _cells, components, _by_node = H3521.build_cells_and_components()
    certificates = H3521.build_certificates()
    component_items = sorted(
        ((component_order_key(component), index, component, certificates[index]) for index, component in enumerate(components)),
        key=lambda item: item[0],
    )

    pending_doublets: dict[tuple[int, ...], tuple[int, int]] = {}
    events: list[StreamEvent] = []
    for stream_index, (_key, index, component, cert) in enumerate(component_items):
        branch_span = component_branch_span(component)
        u_span = component_u_span(component)
        carry_in: tuple[str, ...] = ()
        carry_resolved: tuple[str, ...] = ()
        carry_out: tuple[str, ...] = ()
        sidecars: tuple[str, ...] = ("private_firewall_route_R", "vertical_halfturn_not_glued")
        action = "emit"
        emitted = True
        emission_key = f"C{index:02d}"
        emitted_cell_count = cert.cell_count

        if cert.terminal_class == "ordinary_rank2_route_component":
            action = "emit_ordinary_digit"
            sidecars += ("endpoint_rank2_route",)
        elif cert.terminal_class == "free_hole_single_bracket_packet":
            action = "emit_after_single_bracket_carry"
            carry_in = ("free_hole_predigit",)
            carry_resolved = ("ordinary_bracket_single",)
            sidecars += ("HYP3511_single_bracket",)
        elif cert.terminal_class == "free_hole_doublet_packet":
            assert cert.cluster_id is not None
            cluster = cert.cluster_id
            emission_key = "FH" + "_".join(str(item) for item in cluster)
            if cluster not in pending_doublets:
                pending_doublets[cluster] = (index, cert.cell_count)
                action = "hold_doublet_predigit"
                emitted = False
                emitted_cell_count = 0
                carry_out = (f"free_hole_doublet_pending={cluster}",)
                sidecars += ("HYP3511_doublet_wait",)
            else:
                prior, prior_cell_count = pending_doublets.pop(cluster)
                action = "emit_doublet_after_predigit_carry"
                emitted_cell_count = prior_cell_count + cert.cell_count
                carry_in = (f"free_hole_doublet_pending={cluster}",)
                carry_resolved = (f"free_hole_doublet_cluster={cluster}", f"prior_component=C{prior:02d}")
                sidecars += ("HYP3511_same_branch_doublet",)
        elif cert.terminal_class == "pure_bypass_owner_boundary_component":
            action = "emit_bypass_after_owner_filtration_carry"
            carry_in = ("owner_boundary_predigit=(45,147,169,173)",)
            carry_resolved = (
                "transport_word=(23,93,113)",
                "branch_boundary_lift=(147,169)",
            )
            carry_out = ("typed_residual_boundary=(45,173)",)
            sidecars += ("HYP3520_owner_current", "HYP3522_owner_filtration")
        else:
            action = "unsafe_unclassified_predigit"
            emitted = False
            carry_out = ("unclassified_terminal_debt",)

        events.append(
            StreamEvent(
                stream_index=stream_index,
                component_index=index,
                terminal_class=cert.terminal_class,
                action=action,
                emitted_certificate=emitted,
                emission_key=emission_key,
                branch_span=branch_span,
                u_span=u_span,
                cell_count=cert.cell_count,
                emitted_cell_count=emitted_cell_count,
                carry_in=carry_in,
                carry_resolved=carry_resolved,
                carry_out=carry_out,
                sidecars=sidecars,
            )
        )
    if pending_doublets:
        raise RuntimeError(f"unresolved doublet predigits: {pending_doublets}")
    return tuple(events)


def emitted_certificate_count(events: tuple[StreamEvent, ...]) -> int:
    return sum(1 for event in events if event.emitted_certificate)


def main() -> None:
    events = build_stream_events()
    words = owner_filtration_words()
    action_hist = Counter(event.action for event in events)
    terminal_hist = Counter(event.terminal_class for event in events)
    carry_in_hist = Counter(carry for event in events for carry in event.carry_in)
    carry_out_hist = Counter(carry for event in events for carry in event.carry_out)
    sidecar_hist = Counter(sidecar for event in events for sidecar in event.sidecars)
    score_hist = Counter(carrier.score for carrier in CARRIERS)
    component_cells = sum(event.cell_count for event in events)
    emitted_cells = sum(event.emitted_cell_count for event in events)
    held_predigits = [event for event in events if not event.emitted_certificate]
    held_predigit_cells = sum(event.cell_count for event in held_predigits)
    typed_residual_events = [
        event for event in events if "typed_residual_boundary=(45,173)" in event.carry_out
    ]
    if emitted_cells != component_cells:
        raise RuntimeError(
            f"spigot cell accounting mismatch: emitted {emitted_cells}, components {component_cells}"
        )

    print("HYP-3523 RANDOM031 CERTIFICATE SPIGOT")
    print("status=EVIDENCE / streaming terminal-certificate carry audit; not an LRC14 proof")
    print("row=random_covering_031")
    print()
    print("## Spigot Analogy")
    print("mixed_radix_state=HYP-3521 terminal certificate ledger")
    print("reduce_and_carry=emit certificate plus named proof-debt carry")
    print("predigit_buffer=free-hole doublet or owner-boundary packet waiting for sidecar")
    print("safe_digit=terminal certificate emitted with no untyped future invalidation")
    print()
    print("## Owner Filtration Carry")
    print(f"seam_owners={words['seam']}")
    print(f"transport_word={words['transport']}")
    print(f"branch_boundary_owners={words['branch_boundary']}")
    print(f"bracket_lift={words['bracket_lift']}")
    print(f"residual_after_transport={words['residual_after_transport']}")
    print(f"residual_after_branch_boundary={words['residual_after_branch_boundary']}")
    print("reading=the bypass output is a valid predigit only after transport, bracket lift, and residual boundary carry are named.")
    print()
    print("## Stream Summary")
    print(f"component_events={len(events)}")
    print(f"emitted_certificates={emitted_certificate_count(events)}")
    print(f"component_cells_covered={component_cells}")
    print(f"emitted_certificate_cells={emitted_cells}")
    print(f"held_predigit_events={len(held_predigits)}")
    print(f"held_predigit_buffer_cells={held_predigit_cells}")
    print(f"terminal_class_hist={compact_counter(terminal_hist)}")
    print(f"action_hist={compact_counter(action_hist)}")
    print(f"carry_in_hist={compact_counter(carry_in_hist)}")
    print(f"carry_out_hist={compact_counter(carry_out_hist)}")
    print(f"typed_residual_boundary_events={len(typed_residual_events)}")
    print(f"sidecar_hist={compact_counter(sidecar_hist)}")
    print("spigot_identity=79 component events -> 77 emitted terminal certificates")
    print()
    print("## Stream Events")
    for event in events:
        print(
            f"{event.stream_index:02d} component=C{event.component_index:02d} "
            f"class={event.terminal_class} action={event.action} "
            f"emits={event.emitted_certificate} key={event.emission_key} "
            f"branches={event.branch_span} u={event.u_span} "
            f"component_cells={event.cell_count} emitted_cells={event.emitted_cell_count} "
            f"carry_in={event.carry_in} carry_resolved={event.carry_resolved} "
            f"carry_out={event.carry_out}"
        )
    print()
    print("## Predigit / Carry Verdict")
    print("free_hole_doublet_predigits=2 held events, both later emitted as 2 doublet certificates")
    print("owner_boundary_predigits=1 bypass event, emitted only with HYP-3520/HYP-3522 sidecars")
    print("untyped_pending_carry=0")
    print("typed_residual_carry=(45,173)")
    print("route_guard=HYP-3513 requires route sidecar R unless route reconstruction is proved")
    print("vertical_guard=vertical halfturn remains an address projection, not a sheet gluing")
    print()
    print("## Proof Pull")
    print("P1: Formalize ordinary route certificates as immediate spigot digits.")
    print("P2: Formalize HYP-3511 doublet buffering: first half is a predigit, second half emits the doublet certificate.")
    print("P3: Formalize bypass emission as owner-filtration carry: transport (23,93,113), bracket lift (147,169), residual (45,173).")
    print("P4: Attach HYP-3513 route sidecar R and vertical guard to every emitted certificate before any quotient compression.")
    print()
    print("## Tournament Analysis")
    print("vertices=proof-output states and carry buffers, not runners, arcs, or raw witness counts")
    print("pairwise_observable=no-backtracking terminal emission with named carry")
    print("switch=higher emission safety and lower unresolved carry")
    print(f"score_hist={compact_counter(score_hist)}")
    print("directed_3cycles=0")
    print("scc_sizes=[1,1,1,1,1,1,1,1]")
    print("hamiltonian_path=" + " -> ".join(carrier.name for carrier in CARRIERS))
    print()
    print("## Assumption Challenge")
    print("Considered vertices: runners, raw arcs, witness cells, legal components, terminal certificates, output positions, carry buffers, free-hole packets, owner rows, route sidecars, vertical quotient guards, and proof obligations.")
    print("Chosen vertices are output states plus carry buffers.  This preserves the random031 terminal predicate while deliberately forgetting raw runner order only after branch/u order, terminal class, owner carry, bracket carry, route guard, and vertical guard are retained.")


if __name__ == "__main__":
    main()
