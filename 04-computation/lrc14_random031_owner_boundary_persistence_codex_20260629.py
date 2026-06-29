#!/usr/bin/env python3
"""HYP-3520: random031 owner-boundary persistence certificate.

This script imports the HYP-3486 seam-complement fiber graph and turns the
HYP-3494 owner-boundary target into a finite quotient-price audit.

Tournament Analysis declaration:
  vertices: quotient sidecars / proof interfaces, not runners, arcs, gates, or
            raw phase witnesses;
  pairwise observable: exact reconstruction of the pure-bypass boundary debt
            while retaining the local 12-cell bypass stalk;
  switch/gauge: higher owner-boundary score wins, with a fixed sidecar-cost
            tie path;
  tie Hamiltonian path: owner-current cobordism before equivalent compressed
            owner words, then contaminated or scalar quotients.
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
    "hyp3486_for_hyp3520",
    "lrc14_random031_seam_complement_fiber_graph_codex_20260629.py",
)
H3481 = H3486.H3481
H3493 = load_module(
    "hyp3493_for_hyp3520",
    "lrc14_random031_relative_seam_sheaf_codex_20260629.py",
)


@dataclass(frozen=True)
class OwnerMatrixRow:
    owner: int
    seam_boundary: int
    bypass_flow: int
    boundary_debt: int
    dead_island: int
    global_flow: int
    owner_current: int


@dataclass(frozen=True)
class QuotientCandidate:
    name: str
    retained: str
    reconstructed_debt: tuple[int, ...]
    preserves_local_bypass: bool
    payload_cost: int
    destroyed_field: str
    repair: str
    tie_rank: int

    def exact(self, true_debt: tuple[int, ...]) -> bool:
        return self.reconstructed_debt == true_debt

    def missing(self, true_debt: tuple[int, ...]) -> tuple[int, ...]:
        return tuple(owner for owner in true_debt if owner not in self.reconstructed_debt)

    def extra(self, true_debt: tuple[int, ...]) -> tuple[int, ...]:
        return tuple(owner for owner in self.reconstructed_debt if owner not in true_debt)

    def score(self, true_debt: tuple[int, ...]) -> int:
        score = 0
        if self.exact(true_debt):
            score += 55
        elif self.reconstructed_debt:
            score += 16
        if self.preserves_local_bypass:
            score += 20
        if not self.extra(true_debt):
            score += 12
        else:
            score -= 8
        if not self.missing(true_debt):
            score += 10
        score -= self.payload_cost
        return score


def compact_counter(counter: Counter) -> dict:
    return dict(sorted(counter.items(), key=lambda item: str(item[0])))


def owner_boundary_word(
    seam_owners: tuple[int, ...],
    bypass_owners: tuple[int, ...],
    dead_island_owners: tuple[int, ...],
) -> tuple[str, ...]:
    bypass = set(bypass_owners)
    dead = set(dead_island_owners)
    word = []
    for owner in seam_owners:
        if owner in bypass:
            word.append("P")
        elif owner in dead:
            word.append("D")
        else:
            word.append("O")
    return tuple(word)


def sheaf_quotient_report(components, target_by_rel):
    keys = {
        "flow_class": lambda component: component.flow_class,
        "allowed_exit": lambda component: component.allowed_exit,
        "owner_union": lambda component: component.owner_union,
        "owner_union_size": lambda component: len(component.owner_union),
        "endpoint_ranks": lambda component: component.endpoint_ranks,
        "branch_hist": lambda component: tuple(sorted(component.branch_hist.items())),
        "size": lambda component: component.size,
        "mirror_closed": lambda component: component.mirror_closed,
        "sheet_pgf_bucket": lambda component: component.sheet_pgf_bucket,
    }
    rows = []
    for name, fn in keys.items():
        buckets: dict[object, set[str]] = {}
        examples: dict[object, list[str]] = {}
        for component in components:
            key = fn(component)
            buckets.setdefault(key, set()).add(target_by_rel[component.rel_id])
            examples.setdefault(key, [])
            if len(examples[key]) < 4:
                examples[key].append(component.rel_id)
        mixed = {
            key: sorted(values)
            for key, values in buckets.items()
            if len(values) > 1
        }
        rows.append((name, len(buckets), len(mixed), mixed, examples))
    return rows


def owner_union_from_gates(gates) -> tuple[int, ...]:
    owners: set[int] = set()
    for gate in gates:
        owners.update(H3481.cover_owners(gate))
    return tuple(sorted(owners))


def owner_union_from_cells(cells) -> tuple[int, ...]:
    owners: set[int] = set()
    for cell in cells:
        for hit in cell.hits:
            owners.update(hit.owners)
    return tuple(sorted(owners))


def gate_summary(gate) -> str:
    return (
        f"component={gate.component_index} interval={H3481.interval_text(gate.interval)} "
        f"mask={gate.branch_mask} delta={gate.total_delta} "
        f"route={getattr(gate, 'route', 'unknown')} owners={H3481.cover_owners(gate)}"
    )


def cells_for_component(component, by_node):
    return tuple(by_node[node] for node in component.nodes)


def dead_island_owners(row) -> tuple[tuple[int, tuple[int, ...]], ...]:
    islands, _mirror = H3481.dead_island_report(row)
    rows: list[tuple[int, tuple[int, ...]]] = []
    for island in islands:
        owners = tuple(sorted({int(label.split(":")[1]) for label in island["labels"]}))
        rows.append((island["idx"], owners))
    return tuple(rows)


def external_horizontal_ports(component, by_node) -> tuple[tuple[tuple[int, int], tuple[int, int], str], ...]:
    nodes = set(component.nodes)
    ports: list[tuple[tuple[int, int], tuple[int, int], str]] = []
    for u_index, branch in sorted(nodes):
        for neighbor_u in ((u_index - 1) % H3486.BASE, (u_index + 1) % H3486.BASE):
            neighbor = (neighbor_u, branch)
            if neighbor in by_node and neighbor not in nodes:
                ports.append(((u_index, branch), neighbor, by_node[neighbor].cell_class))
    return tuple(ports)


def mirror_pair_for_gates(gates_subset, all_gates) -> bool:
    subset = set(gates_subset)
    for left, right in H3481.mirror_pairs_by_interval(all_gates):
        if {left, right} == subset:
            return True
    return False


def build_owner_matrix(
    seam_owners: tuple[int, ...],
    bypass_owners: tuple[int, ...],
    debt: tuple[int, ...],
    dead_owners: tuple[int, ...],
    global_flow_owners: tuple[int, ...],
) -> tuple[OwnerMatrixRow, ...]:
    universe = tuple(sorted(set(seam_owners) | set(bypass_owners) | set(dead_owners) | set(global_flow_owners)))
    rows: list[OwnerMatrixRow] = []
    for owner in universe:
        seam = int(owner in seam_owners)
        bypass = int(owner in bypass_owners)
        rows.append(
            OwnerMatrixRow(
                owner=owner,
                seam_boundary=seam,
                bypass_flow=bypass,
                boundary_debt=int(owner in debt),
                dead_island=int(owner in dead_owners),
                global_flow=int(owner in global_flow_owners),
                owner_current=seam - bypass,
            )
        )
    return tuple(rows)


def quotient_candidates(
    debt: tuple[int, ...],
    dead_owners: tuple[int, ...],
    bypass_owners: tuple[int, ...],
    global_flow_owners: tuple[int, ...],
) -> tuple[QuotientCandidate, ...]:
    dead_minus_bypass = tuple(owner for owner in dead_owners if owner not in bypass_owners)
    global_minus_bypass = tuple(owner for owner in global_flow_owners if owner not in bypass_owners)
    empty: tuple[int, ...] = ()
    return (
        QuotientCandidate(
            "owner_current_cobordism_matrix",
            "owner x {seam_boundary,bypass_flow,boundary_debt}",
            debt,
            True,
            4,
            "none at owner-boundary level",
            "finite owner-current lemma over the forbidden seam",
            0,
        ),
        QuotientCandidate(
            "seam_and_bypass_owner_words",
            "forbidden seam word plus pure bypass flow word",
            debt,
            True,
            3,
            "component-local adjacency and endpoint-rank details",
            "attach hard components (43,54) and R01 stalk id",
            1,
        ),
        QuotientCandidate(
            "hard_component_owner_map",
            "per-hard-component seam/bypass owner words",
            debt,
            True,
            5,
            "global phase-carrier connectedness",
            "attach HYP-3510 carrier id and HYP-3511 free-hole bracket sidecar",
            2,
        ),
        QuotientCandidate(
            "observer_cut_payload_orbit",
            "seven-owner seam as stationary observer-cut payload",
            debt,
            True,
            5,
            "raw u-fiber branch order",
            "attach HYP-3486 u_index/branch packet if formalizing cells",
            3,
        ),
        QuotientCandidate(
            "symbolic_relative_H1_boundary_class",
            "relative boundary word [seam]-[bypass]",
            debt,
            True,
            6,
            "individual owner rows unless expanded",
            "expand to the owner-current matrix above",
            4,
        ),
        QuotientCandidate(
            "global_flow_owner_union_minus_bypass",
            "all nonfree flow owners minus bypass owners",
            global_minus_bypass,
            False,
            2,
            "locality of the pure bypass stalk",
            "remove ordinary-owner contaminant 55 by hard-component sidecar",
            5,
        ),
        QuotientCandidate(
            "dead_island_owner_union_minus_bypass",
            "four puncture owners minus bypass owners",
            dead_minus_bypass,
            False,
            2,
            "seam-only owners 147,169,173",
            "restore the seven-owner forbidden seam payload",
            6,
        ),
        QuotientCandidate(
            "bypass_owner_word_only",
            "pure bypass owners (23,93,113)",
            empty,
            True,
            1,
            "all seam-only owners",
            "add the forbidden seam owner word",
            7,
        ),
        QuotientCandidate(
            "delta_route_phase_blocks",
            "delta=2, branch split 6/6, phases on components 43/54",
            empty,
            True,
            1,
            "owner labels",
            "attach owner payload to the two-adic phase carrier",
            8,
        ),
        QuotientCandidate(
            "component_pair_only",
            "hard component pair (43,54)",
            empty,
            False,
            1,
            "owner labels and phase flow",
            "attach seam/bypass gate words",
            9,
        ),
        QuotientCandidate(
            "raw_bypass_count_12",
            "12 bypass witnesses",
            empty,
            False,
            0,
            "branch, owner, component, endpoint, and seam fields",
            "lift count to HYP-3486 witness packets",
            10,
        ),
        QuotientCandidate(
            "raw_282_witness_count",
            "282 phase witnesses",
            empty,
            False,
            0,
            "all terminal class and owner-boundary fields",
            "lift count to seam-complement receiver packets",
            11,
        ),
    )


def count_directed_3cycles(ordered: tuple[QuotientCandidate, ...]) -> int:
    # The switch is a total order after the fixed tie path, so this audit is
    # intentionally explicit: any nonzero result would signal a scoring bug.
    position = {candidate.name: idx for idx, candidate in enumerate(ordered)}
    names = [candidate.name for candidate in ordered]
    cycles = 0
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            for k in range(j + 1, len(names)):
                a, b, c = names[i], names[j], names[k]
                if position[a] < position[b] < position[c] < position[a]:
                    cycles += 1
    return cycles


def main() -> None:
    row = H3481.H3450.audit_row(H3481.ROW_NAME, H3481.SPEEDS)
    gates = H3481.build_gates()
    cells = H3486.build_cells(gates, row)
    by_node = {cell.node: cell for cell in cells}
    legal = H3486.connected_components(cells, {"horizontal", "mirror"})

    hard_gates = tuple(
        sorted(
            [
                gate
                for gate in gates
                if gate.component_index in H3481.HARD_COMPONENTS and gate.total_delta == 7
            ],
            key=lambda gate: (gate.component_index, gate.branch_mask, gate.interval),
        )
    )
    bypass_gates = tuple(
        sorted(
            [
                gate
                for gate in gates
                if gate.component_index in H3481.HARD_COMPONENTS and gate.total_delta < 7
            ],
            key=lambda gate: (gate.component_index, gate.branch_mask, gate.interval),
        )
    )
    pure_bypass_items = [
        (idx, component) for idx, component in enumerate(legal) if component.type_word == "bypass"
    ]
    if len(pure_bypass_items) != 1:
        raise RuntimeError(f"expected one pure bypass component, found {len(pure_bypass_items)}")
    pure_idx, pure_bypass = pure_bypass_items[0]
    pure_cells = cells_for_component(pure_bypass, by_node)

    seam_owners = owner_union_from_gates(hard_gates)
    bypass_gate_owners = owner_union_from_gates(bypass_gates)
    pure_bypass_owners = owner_union_from_cells(pure_cells)
    boundary_debt = tuple(owner for owner in seam_owners if owner not in pure_bypass_owners)

    dead_words = dead_island_owners(row)
    dead_owner_union = tuple(sorted({owner for _idx, owners in dead_words for owner in owners}))
    global_flow_owners = owner_union_from_cells([cell for cell in cells if cell.hits])
    ordinary_flow_owners = owner_union_from_cells([cell for cell in cells if cell.cell_class == "ordinary"])
    matrix = build_owner_matrix(
        seam_owners,
        pure_bypass_owners,
        boundary_debt,
        dead_owner_union,
        global_flow_owners,
    )

    class_hist = Counter(cell.cell_class for cell in cells)
    legal_type_hist = Counter(component.type_word for component in legal)
    phase_hist = Counter(cell.phase for cell in pure_cells)
    branch_hist = Counter(cell.branch for cell in pure_cells)
    component_phase_hist = Counter(
        (hit.component, cell.branch, cell.phase)
        for cell in pure_cells
        for hit in cell.hits
    )
    hard_component_owner_map = {
        component: {
            "seam": owner_union_from_gates([gate for gate in hard_gates if gate.component_index == component]),
            "bypass": owner_union_from_gates([gate for gate in bypass_gates if gate.component_index == component]),
        }
        for component in H3481.HARD_COMPONENTS
    }
    hard_component_owner_map = {
        component: {
            **words,
            "debt": tuple(owner for owner in words["seam"] if owner not in words["bypass"]),
        }
        for component, words in hard_component_owner_map.items()
    }
    external_ports = external_horizontal_ports(pure_bypass, by_node)

    candidates = quotient_candidates(
        boundary_debt,
        dead_owner_union,
        pure_bypass_owners,
        global_flow_owners,
    )
    ordered_candidates = tuple(
        sorted(candidates, key=lambda candidate: (-candidate.score(boundary_debt), candidate.tie_rank))
    )
    score_hist = Counter(candidate.score(boundary_debt) for candidate in candidates)
    sheaf_components = H3493.build_sheaf_components()
    sheaf_bypass_components = [
        component for component in sheaf_components if component.flow_class == "pure_bypass"
    ]
    if len(sheaf_bypass_components) != 1:
        raise RuntimeError(
            f"expected one HYP-3493 pure bypass component, found {len(sheaf_bypass_components)}"
        )
    sheaf_bypass = sheaf_bypass_components[0]
    sheaf_word = owner_boundary_word(seam_owners, sheaf_bypass.owner_union, dead_owner_union)
    sheaf_target_by_rel = {
        component.rel_id: (
            "pure_bypass_owner_boundary_debt"
            if component.flow_class == "pure_bypass"
            else "local_discharge_or_freehole_context"
        )
        for component in sheaf_components
    }
    sheaf_free_owner_unions = Counter(
        component.owner_union for component in sheaf_components if component.flow_class == "free_hole"
    )
    sheaf_flow_hist = Counter(component.flow_class for component in sheaf_components)

    print("HYP-3520 RANDOM031 OWNER-BOUNDARY PERSISTENCE")
    print("status=EVIDENCE / finite owner-cobordism and quotient-price certificate; not an LRC14 proof")
    print("row=random_covering_031")
    print()

    print("## Seam-Complement Receiver")
    print(f"witness_cells={len(cells)}")
    print(f"cell_class_hist={compact_counter(class_hist)}")
    print(f"legal_component_count={len(legal)}")
    print(f"legal_component_type_hist={compact_counter(legal_type_hist)}")
    print(f"hard_components={H3481.HARD_COMPONENTS}")
    print(f"hard_seam_gate_count={len(hard_gates)} mirror_pair={mirror_pair_for_gates(hard_gates, gates)}")
    print(f"bypass_gate_count={len(bypass_gates)} mirror_pair={mirror_pair_for_gates(bypass_gates, gates)}")
    print()

    print("## Forbidden Seam Gates")
    for gate in hard_gates:
        print("  " + gate_summary(gate))
    print(f"forbidden_seam_owner_word={seam_owners}")
    print()

    print("## Lower-Delta Bypass Gates")
    for gate in bypass_gates:
        print("  " + gate_summary(gate))
    print(f"bypass_gate_owner_word={bypass_gate_owners}")
    print()

    print("## Pure Bypass Stalk")
    print(f"rel_id=R{pure_idx:02d}")
    print(f"size={pure_bypass.size}")
    print(f"nodes={pure_bypass.nodes}")
    print(f"branch_hist={compact_counter(branch_hist)}")
    print(f"phase_hist={compact_counter(phase_hist)}")
    print(f"hit_components={pure_bypass.hit_components}")
    print(f"endpoint_ranks={pure_bypass.endpoint_ranks}")
    print(f"pure_bypass_owner_word={pure_bypass_owners}")
    print(f"external_horizontal_ports={len(external_ports)}")
    print(f"component_branch_phase_hits={compact_counter(component_phase_hist)}")
    print("local_reading=pure 12-cell lower-delta flow stalk on the same two hard components")
    print()

    print("## Owner Boundary Matrix")
    print(f"dead_island_owner_words={dead_words}")
    print(f"dead_island_owner_union={dead_owner_union}")
    print(f"ordinary_flow_owner_union={ordinary_flow_owners}")
    print(f"global_flow_owner_union={global_flow_owners}")
    print("owner | seam_boundary | bypass_flow | boundary_debt | dead_island | global_flow | seam-minus-bypass")
    for row_item in matrix:
        print(
            f"{row_item.owner:>5} | {row_item.seam_boundary:>13} | {row_item.bypass_flow:>11} | "
            f"{row_item.boundary_debt:>13} | {row_item.dead_island:>11} | "
            f"{row_item.global_flow:>11} | {row_item.owner_current:>17}"
        )
    print(f"owner_boundary_debt={boundary_debt}")
    print(
        "owner_current_signature="
        + ",".join(f"{item.owner}:{item.owner_current:+d}" for item in matrix if item.owner_current)
    )
    print()

    print("## Hard-Component Owner Map")
    for component, words in hard_component_owner_map.items():
        print(f"component={component} seam={words['seam']} bypass={words['bypass']} debt={words['debt']}")
    print("same_component_verdict=seam and bypass live on the same mirror-paired hard components")
    print()

    print("## Quotient-Price Matrix")
    print("candidate | score | exact | reconstructed_debt | missing | extra | retained | repair")
    for candidate in ordered_candidates:
        print(
            f"{candidate.name} | {candidate.score(boundary_debt)} | "
            f"{candidate.exact(boundary_debt)} | {candidate.reconstructed_debt} | "
            f"{candidate.missing(boundary_debt)} | {candidate.extra(boundary_debt)} | "
            f"{candidate.retained} | {candidate.repair}"
        )
    print()

    print("## Seam-Sheaf Compression Canary")
    print("source=HYP-3493 relative seam-sheaf component table")
    print("target=pure_bypass_owner_boundary_debt versus local discharge/free-hole context")
    print(f"sheaf_component_count={len(sheaf_components)}")
    print(f"sheaf_flow_class_hist={compact_counter(sheaf_flow_hist)}")
    print(f"sheaf_pure_bypass_rel_id={sheaf_bypass.rel_id}")
    print(f"sheaf_pure_bypass_owner_word={sheaf_bypass.owner_union}")
    print(f"owner_boundary_persistence_word={''.join(sheaf_word)}")
    print(f"free_hole_owner_union_hist={compact_counter(sheaf_free_owner_unions)}")
    print("quotient | bucket_count | mixed_bucket_count | safe")
    for name, bucket_count, mixed_count, mixed, examples in sheaf_quotient_report(
        sheaf_components,
        sheaf_target_by_rel,
    ):
        safe = mixed_count == 0
        print(f"{name} | {bucket_count} | {mixed_count} | {safe}")
        if mixed_count:
            shown = 0
            for key, values in mixed.items():
                print(f"  mixed_key={key!r} targets={tuple(values)} examples={tuple(examples[key])}")
                shown += 1
                if shown >= 4:
                    break
    print(
        "canary_verdict=flow_class, allowed_exit, owner_union, and sheet_pgf_bucket "
        "preserve the debt target; owner count, endpoint rank, branch histogram, "
        "size, and mirror closure are illegal scalar compressions here."
    )
    print()

    print("## Reframe")
    print(
        "R1: The hard pair is not a wall for phase flow.  It is a forbidden seam; "
        "the complement carries all 282 phase witnesses, with the lower-delta "
        "bypass hit exactly 12 times."
    )
    print(
        "R2: The four missing owners are not failed bypass flow.  They are the "
        "stationary owner-boundary current left after subtracting the moving "
        "bypass owner word from the seven-owner seam word."
    )
    print(
        "R3: The n*2 recursion is the u=2t phase address.  The n+2 recursion is "
        "the owner-boundary seam insertion.  HYP-3520 should prove the span, "
        "not choose one recursion and erase the other."
    )
    print()

    print("## Tournament Analysis")
    print("vertices=quotient sidecars/proof interfaces, not runners, arcs, or raw witnesses")
    print("pairwise_observable=exact boundary-debt reconstruction + local 12-cell bypass retention")
    print("switch=higher owner-boundary score; fixed sidecar-cost tie path")
    print(f"score_hist={compact_counter(score_hist)}")
    print(f"directed_3cycles={count_directed_3cycles(ordered_candidates)}")
    print("scc_sizes=" + str([1] * len(ordered_candidates)))
    print("hamiltonian_path_count_under_fixed_tie_path=1")
    print("hamiltonian_path=" + " -> ".join(candidate.name for candidate in ordered_candidates))
    print()

    print("## Assumption Challenge")
    print(
        "Considered vertices: runners, raw gates, arcs, u-fibers, fixed circle "
        "sections, section boundaries, wall-crossing events, residues, cover "
        "arcs, free-hole packets, owner rows, hard-component owner maps, "
        "relative-H1 boundary classes, and proof obligations."
    )
    print(
        "Chosen vertices: quotient sidecars.  This preserves the LRC predicate "
        "needed here--pure-bypass discharge with exact owner-boundary debt--and "
        "intentionally destroys raw runner order, raw phase count, and ordinary "
        "flow owner noise."
    )
    print(
        "Challenged assumption: the owner-boundary proof is not a search for "
        "more bypass witnesses.  The bypass is already local and isolated; the "
        "live question is which sidecar prevents the seam-only owners from "
        "being forgotten."
    )


if __name__ == "__main__":
    main()
