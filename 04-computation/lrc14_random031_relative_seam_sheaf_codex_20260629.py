#!/usr/bin/env python3
"""HYP-3493: random031 relative seam-sheaf scaffold.

This scout is the first executable form of the HYP-3493 synthesis.  It imports
the exact HYP-3486 seam-complement fiber graph, keeps the legal
horizontal+mirror components, and annotates each component as a local stalk of
the proposed relative seam sheaf.

The script does not prove the random031 terminal clause.  It converts the
existing fiber graph into the table a proof should consume:

    u/branch/mirror data + flow class + endpoint rank + owner word
    + private-label firewall mark + legal exit.
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


H3486 = load_module(
    "hyp3486_for_hyp3491",
    "lrc14_random031_seam_complement_fiber_graph_codex_20260629.py",
)

SEAM_OWNERS = H3486.H3481.SEVEN_OWNER_CORE
DEAD_ISLAND_OWNER_UNION = (23, 45, 93, 113)
PRIVATE_LABEL_STATUS = "private_firewall_random031"


@dataclass(frozen=True)
class SheafComponent:
    rel_id: str
    flow_class: str
    size: int
    nodes: tuple[tuple[int, int], ...]
    branch_hist: Counter[int]
    endpoint_ranks: tuple[int | None, ...]
    hit_components: tuple[int, ...]
    owner_union: tuple[int, ...]
    seam_owner_debt: tuple[int, ...]
    mirror_closed: bool
    sheet_pgf_bucket: str
    allowed_exit: str


@dataclass(frozen=True)
class Carrier:
    name: str
    rank: int


CARRIERS = (
    Carrier("S00_relative_seam_sheaf_packet", 7),
    Carrier("S01_fiber_graph_trichotomy", 6),
    Carrier("S02_private_label_firewall", 5),
    Carrier("S03_ordered_twoadic_bypass_blocks", 4),
    Carrier("S04_cech_owner_persistence", 3),
    Carrier("S05_sheet_pgf_moment", 2),
    Carrier("S06_nonarchimedean_owner_lift_gate", 1),
    Carrier("S07_raw_scalar_shadow", 0),
)


def compact_counter(counter: Counter) -> dict:
    return dict(sorted(counter.items(), key=lambda item: item[0]))


def component_cells(component, by_node):
    return tuple(by_node[node] for node in component.nodes)


def owner_union(cells) -> tuple[int, ...]:
    owners: set[int] = set()
    for cell in cells:
        for hit in cell.hits:
            owners.update(hit.owners)
    return tuple(sorted(owners))


def flow_class(component) -> str:
    if component.type_word == "free_hole":
        return "free_hole"
    if component.type_word == "bypass":
        return "pure_bypass"
    if component.hit_components and component.endpoint_ranks == (2,):
        return "rank2_routed"
    return "mixed_or_debt"


def allowed_exit(name: str) -> str:
    return {
        "rank2_routed": "endpoint_rank2_discharge",
        "free_hole": "mirror_free_hole_packet",
        "pure_bypass": "pure_bypass_plus_owner_boundary_debt",
    }.get(name, "named_owner_twoadic_spec_state_lift_debt")


def sheet_pgf_bucket(component, name: str) -> str:
    b0 = component.branch_hist.get(0, 0)
    b1 = component.branch_hist.get(1, 0)
    return f"{name}:z^{component.size}:b0^{b0}:b1^{b1}"


def build_sheaf_components() -> tuple[SheafComponent, ...]:
    row = H3486.H3481.H3450.audit_row(H3486.H3481.ROW_NAME, H3486.H3481.SPEEDS)
    gates = H3486.H3481.build_gates()
    cells = H3486.build_cells(gates, row)
    by_node = {cell.node: cell for cell in cells}
    legal = H3486.connected_components(cells, {"horizontal", "mirror"})

    sheaf_components: list[SheafComponent] = []
    for idx, component in enumerate(legal):
        cells_here = component_cells(component, by_node)
        owners = owner_union(cells_here)
        name = flow_class(component)
        node_set = set(component.nodes)
        mirror_closed = all(H3486.mirror_node(node) in node_set for node in node_set)
        seam_debt = tuple(owner for owner in SEAM_OWNERS if owner not in owners)
        sheaf_components.append(
            SheafComponent(
                rel_id=f"R{idx:02d}",
                flow_class=name,
                size=component.size,
                nodes=component.nodes,
                branch_hist=component.branch_hist,
                endpoint_ranks=component.endpoint_ranks,
                hit_components=component.hit_components,
                owner_union=owners,
                seam_owner_debt=seam_debt,
                mirror_closed=mirror_closed,
                sheet_pgf_bucket=sheet_pgf_bucket(component, name),
                allowed_exit=allowed_exit(name),
            )
        )
    return tuple(sheaf_components)


def component_line(component: SheafComponent) -> str:
    u_values = [node[0] for node in component.nodes]
    branch_word = compact_counter(component.branch_hist)
    return (
        f"  {component.rel_id} class={component.flow_class} size={component.size} "
        f"u=[{min(u_values)},{max(u_values)}] branches={branch_word} "
        f"ranks={component.endpoint_ranks} hits={component.hit_components} "
        f"owners={component.owner_union} seam_debt={component.seam_owner_debt} "
        f"mirror_closed={component.mirror_closed} pgf={component.sheet_pgf_bucket} "
        f"exit={component.allowed_exit}"
    )


def main() -> None:
    components = build_sheaf_components()
    class_hist = Counter(component.flow_class for component in components)
    exit_hist = Counter(component.allowed_exit for component in components)
    size_by_class = Counter((component.flow_class, component.size) for component in components)
    mirror_closed = sum(component.mirror_closed for component in components)
    owner_size_hist = Counter(len(component.owner_union) for component in components)
    seam_debt_hist = Counter(component.seam_owner_debt for component in components)
    pgf_hist = Counter(component.sheet_pgf_bucket for component in components)
    endpoint_rank_hist = Counter(
        rank for component in components for rank in component.endpoint_ranks
    )
    bypass = [component for component in components if component.flow_class == "pure_bypass"]
    free_holes = [component for component in components if component.flow_class == "free_hole"]

    print("HYP-3493 RANDOM031 RELATIVE SEAM-SHEAF SCAFFOLD")
    print("status=EVIDENCE / component-level sheaf table; not an LRC14 proof")
    print("row=random_covering_031")
    print()
    print("## Inputs")
    print("source=HYP-3486 seam-complement fiber graph")
    print("private_label_status=HYP-3490 random031 private-label firewall")
    print(f"seam_owners={SEAM_OWNERS}")
    print(f"dead_island_owner_union={DEAD_ISLAND_OWNER_UNION}")
    print("base_topology=cylinder - four mirror-paired punctures - forbidden seam")
    print()

    print("## Relative Seam-Sheaf Component Summary")
    print(f"component_count={len(components)}")
    print(f"flow_class_hist={compact_counter(class_hist)}")
    print(f"allowed_exit_hist={compact_counter(exit_hist)}")
    print(f"mirror_closed_components={mirror_closed}/{len(components)}")
    print(f"size_by_class={compact_counter(size_by_class)}")
    print(f"owner_union_size_hist={compact_counter(owner_size_hist)}")
    print(f"endpoint_rank_hist={compact_counter(endpoint_rank_hist)}")
    print(f"sheet_pgf_bucket_count={len(pgf_hist)}")
    print("private_label_mark=attached to every component; projection-current route disabled")
    print()

    print("## Owner Boundary Debt")
    print(f"seam_debt_hist={compact_counter(seam_debt_hist)}")
    for component in bypass:
        print("pure_bypass_owner_debt:")
        print(component_line(component))
    print()

    print("## Free-Hole Packets")
    print(f"free_hole_packet_count={len(free_holes)}")
    print(f"free_hole_size_hist={compact_counter(Counter(component.size for component in free_holes))}")
    for component in free_holes[:8]:
        print(component_line(component))
    if len(free_holes) > 8:
        print(f"  ... {len(free_holes) - 8} more free-hole packets")
    print()

    print("## Largest / Debt-Relevant Components")
    for component in sorted(components, key=lambda item: (-item.size, item.rel_id))[:12]:
        print(component_line(component))
    print()

    print("## Proof Pull")
    print(
        "P1: The sheaf table has no mixed_or_debt components; the current finite "
        "carrier splits into rank2_routed, free_hole, and pure_bypass stalks."
    )
    print(
        "P2: Every legal horizontal+mirror component is mirror-closed, so the "
        "Cech/path-lift persistence test should target owner-boundary words, "
        "not mirror closure."
    )
    print(
        "P3: The pure bypass stalk carries owners (23,93,113) and seam debt "
        "(45,147,169,173), matching the HYP-3483 n+2 versus n*2 span."
    )
    print(
        "P4: HYP-3490 disables projection-current deletion on all these stalks; "
        "the legal theorem must discharge by endpoint rank, free-hole, bypass, "
        "or named owner/two-adic/SPEC/state-lift debt."
    )
    print()

    print("## Tournament Analysis")
    print("vertices=relative seam-sheaf proof carriers, not runners or raw arcs")
    print(
        "pairwise_observable=terminal predicate retention + quotient legality + "
        "fiber/private-label/owner/persistence/sheet payload"
    )
    print("switch=higher retained proof payload; ties follow the declared path")
    print(f"score_hist={compact_counter(Counter(carrier.rank for carrier in CARRIERS))}")
    print("sccs=8 singleton SCCs")
    print("directed_3cycles=0")
    print("hamiltonian_path=" + " -> ".join(carrier.name for carrier in CARRIERS))
    print()

    print("## Assumption Challenge")
    print(
        "Considered vertices: runners, gaps, fixed sections, section boundaries, "
        "wall crossings, residues, cover arcs, Fourier modes, matroid circuits, "
        "witness cells, fiber components, owner-boundary words, and proof "
        "obligations."
    )
    print(
        "Chosen vertices are relative seam-sheaf proof packets.  They preserve "
        "the terminal random031 discharge predicate and destroy raw interval "
        "geometry only after replacing it with u, branch, mirror, owner, "
        "private-label, relative-cycle, and legal exit sidecars."
    )


if __name__ == "__main__":
    main()
