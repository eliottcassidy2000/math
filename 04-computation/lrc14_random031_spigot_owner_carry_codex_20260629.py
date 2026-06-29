#!/usr/bin/env python3
"""HYP-3523: random031 spigot owner-carry audit.

This script uses the spigot-algorithm prompt as a proof-design analogy:
emit only the owner certificates that are locally stable, keep the remaining
owner debt as a bounded carry state, and reject quotients whose "digit" is
ambiguous across terminal classes.

Tournament Analysis declaration:
  vertices: streaming proof interfaces / observer states, not runners or arcs;
  pairwise observable: stable emitted owner word, residual carry width, and
            quotient safety against the HYP-3520/HYP-3522 terminal predicate;
  switch/gauge: safe emission before smaller residual carry before lower
            sidecar cost;
  tie Hamiltonian path: exact spigot state, filtration words, persistence word,
            then scalar shadows.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
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


H3520 = load_module(
    "hyp3520_for_hyp3523",
    "lrc14_random031_owner_boundary_persistence_codex_20260629.py",
)
H3522 = load_module(
    "hyp3522_for_hyp3523",
    "lrc14_random031_owner_boundary_filtration_codex_20260629.py",
)


@dataclass(frozen=True)
class EmissionStage:
    name: str
    observer: str
    emitted: tuple[int, ...]
    incoming_carry: tuple[int, ...]
    outgoing_carry: tuple[int, ...]
    stable: bool
    certificate: str
    destroyed: str
    proof_pull: str


@dataclass(frozen=True)
class StreamInterface:
    name: str
    emitted: tuple[int, ...]
    residual: tuple[int, ...]
    stable: bool
    contaminant: tuple[int, ...]
    sidecar_cost: int
    destroys: str
    repair: str
    tie_rank: int

    def score(self) -> int:
        score = 0
        if self.stable:
            score += 70
        score += 18 if self.emitted else 0
        score += 18 if self.residual == (45, 173) else 0
        score += max(0, 10 - len(self.residual))
        score -= 12 * len(self.contaminant)
        score -= self.sidecar_cost
        return score


def compact_counter(counter: Counter) -> dict:
    try:
        return dict(sorted(counter.items(), key=lambda item: item[0]))
    except TypeError:
        return dict(sorted(counter.items(), key=lambda item: str(item[0])))


def difference(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    right_set = set(right)
    return tuple(owner for owner in left if owner not in right_set)


def union_words(*words: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted({owner for word in words for owner in word}))


def boundary_word_by_branch(boundaries) -> dict[int, tuple[int, ...]]:
    out = {}
    for boundary in boundaries:
        left = H3522.owner_union_from_hits(boundary.left_cell.hits)
        right = H3522.owner_union_from_hits(boundary.right_cell.hits)
        out[boundary.branch] = union_words(left, right)
    return out


def boundary_intersections(boundaries) -> dict[int, tuple[int, ...]]:
    out = {}
    for boundary in boundaries:
        left = set(H3522.owner_union_from_hits(boundary.left_cell.hits))
        right = set(H3522.owner_union_from_hits(boundary.right_cell.hits))
        out[boundary.branch] = tuple(sorted(left & right))
    return out


def directed_3cycle_count(ordered: tuple[StreamInterface, ...]) -> int:
    position = {item.name: index for index, item in enumerate(ordered)}
    cycles = 0
    for a, b, c in combinations(ordered, 3):
        ab = a if position[a.name] < position[b.name] else b
        bc = b if position[b.name] < position[c.name] else c
        ca = c if position[c.name] < position[a.name] else a
        if ab is a and bc is b and ca is c:
            cycles += 1
        if ab is b and bc is c and ca is a:
            cycles += 1
    return cycles


def build_context():
    (
        row,
        gates,
        cells,
        by_node,
        legal,
        bypass_component,
        bypass_cells,
        hard_seam_gates,
        lower_bypass_gates,
    ) = H3522.build_context()
    boundaries = H3522.branch_boundaries(cells, by_node, bypass_component)
    mirror = H3522.mirror_pairs(by_node, bypass_component)
    persistence = H3520.persistence_summary()
    dead_owners = tuple(
        sorted(
            {
                owner
                for component in row.dead_components
                for owner in H3522.component_owner_pair(component)
            }
        )
    )
    return {
        "row": row,
        "gates": gates,
        "cells": cells,
        "by_node": by_node,
        "legal": legal,
        "bypass_component": bypass_component,
        "bypass_cells": bypass_cells,
        "hard_seam_gates": hard_seam_gates,
        "lower_bypass_gates": lower_bypass_gates,
        "boundaries": boundaries,
        "mirror": mirror,
        "persistence": persistence,
        "dead_owners": dead_owners,
    }


def build_stream(ctx) -> dict:
    seam = H3522.gate_owner_union(ctx["hard_seam_gates"])
    transport = H3522.owner_union_from_cells(ctx["bypass_cells"])
    lower_bypass = H3522.gate_owner_union(ctx["lower_bypass_gates"])
    boundary_cells = tuple(
        cell
        for boundary in ctx["boundaries"]
        for cell in (boundary.left_cell, boundary.right_cell)
    )
    branch_boundary = H3522.owner_union_from_cells(boundary_cells)
    branch_by_branch = boundary_word_by_branch(ctx["boundaries"])
    branch_intersections = boundary_intersections(ctx["boundaries"])
    branch_lift = difference(branch_boundary, transport)
    transport_plus_branch = union_words(transport, branch_boundary)
    residual_after_transport = difference(seam, transport)
    residual_after_branch = difference(seam, transport_plus_branch)
    transport_only = difference(transport, branch_boundary)
    mirror_owner_persistent = all(
        H3522.owner_union_from_hits(left.hits) == H3522.owner_union_from_hits(right.hits)
        for left, right in ctx["mirror"]
    )
    bypass_owner_words = Counter(
        H3522.owner_union_from_hits(cell.hits) for cell in ctx["bypass_cells"]
    )
    lower_matches_transport = lower_bypass == transport
    branch_word_constant = len(set(branch_by_branch.values())) == 1
    branch_intersection_constant = len(set(branch_intersections.values())) == 1
    free_stats = H3522.free_hole_bracket_stats()
    persistence_word = "".join(
        H3520.owner_boundary_word(seam, transport, ctx["dead_owners"])
    )
    global_flow = H3520.owner_union_from_cells([cell for cell in ctx["cells"] if cell.hits])
    global_minus_bypass = difference(global_flow, transport)

    stages = (
        EmissionStage(
            "S0_load_forbidden_seam",
            "max_delta_seam_owner_word",
            (),
            (),
            seam,
            True,
            "two mirror seam gates on components (43,54)",
            "phase order until u/branch sidecar is attached",
            "establish the seven-owner initial payload",
        ),
        EmissionStage(
            "S1_emit_transport_word",
            "lower_delta_bypass_stalk",
            transport,
            seam,
            residual_after_transport,
            lower_matches_transport and mirror_owner_persistent,
            "all 12 bypass cells and six mirror pairs carry (23,93,113)",
            "seam-only owners unless carry is retained",
            "prove transport-word constancy",
        ),
        EmissionStage(
            "S2_emit_branch_bracket_lift",
            "ordinary_neighbors_of_bypass_interval",
            branch_lift,
            residual_after_transport,
            residual_after_branch,
            branch_word_constant and branch_intersection_constant,
            "both branch sheets have boundary word (23,93,147,169) and intersection (169,)",
            "raw branch order after the bracket sidecar is emitted",
            "prove branch-boundary bracket lift",
        ),
        EmissionStage(
            "S3_emit_free_hole_cap_certificates",
            "HYP-3511_bracketed_free_holes",
            (),
            residual_after_branch,
            residual_after_branch,
            (
                free_stats["packet_count"] == 14
                and free_stats["bracketed_count"] + free_stats["half_open_count"] == 14
                and ctx["persistence"]["unresolved"] == []
            ),
            "14 free-hole packets are bracketed as 10 singles plus 2 half-open doublets",
            "free-hole cell order after bracket/doublet id is retained",
            "splice HYP-3511 into the owner-boundary terminal dispatch",
        ),
        EmissionStage(
            "S4_defer_residual_pair",
            "puncture_apex_owner_carry",
            (),
            residual_after_branch,
            residual_after_branch,
            residual_after_branch == (45, 173),
            "residual carry is the two-owner pair (45,173)",
            "terminal closure until the HYP-3490/HYP-3513 route sidecar is proved",
            "prove the residual two-owner boundary lemma",
        ),
    )

    interfaces = (
        StreamInterface(
            "exact_spigot_owner_state",
            union_words(transport, branch_lift),
            residual_after_branch,
            True,
            (),
            4,
            "nothing terminal; keeps transport, lift, residual, and route hook",
            "formalize as finite-state owner-carry transducer",
            0,
        ),
        StreamInterface(
            "transport_branch_residual_words",
            union_words(transport, branch_lift),
            residual_after_branch,
            True,
            (),
            5,
            "raw u-order and individual cell ids",
            "attach HYP-3486 cell packet only inside local lemmas",
            1,
        ),
        StreamInterface(
            "owner_boundary_persistence_word",
            tuple(),
            residual_after_branch,
            True,
            (),
            4,
            "which owner was emitted in which stage",
            "expand PDPPOOO by HYP-3522 filtration when proving residual",
            2,
        ),
        StreamInterface(
            "seam_and_bypass_owner_words",
            transport,
            residual_after_transport,
            True,
            (),
            3,
            "branch-boundary lift (147,169)",
            "add ordinary-neighbor bracket sidecar",
            3,
        ),
        StreamInterface(
            "flow_class_owner_union_sheet_pgf",
            transport,
            residual_after_transport,
            True,
            (),
            4,
            "owner-current order inside the pure bypass",
            "attach HYP-3522 filtration before residual lemma",
            4,
        ),
        StreamInterface(
            "global_flow_minus_bypass",
            transport,
            global_minus_bypass,
            False,
            tuple(owner for owner in global_minus_bypass if owner not in seam),
            2,
            "pure-bypass locality",
            "remove contaminant 55 with hard-component owner map",
            5,
        ),
        StreamInterface(
            "bypass_owner_word_only",
            transport,
            (),
            False,
            (),
            1,
            "all stationary owner-boundary carry",
            "add forbidden seam owner word",
            6,
        ),
        StreamInterface(
            "endpoint_rank_shadow",
            (),
            (),
            False,
            (),
            1,
            "owner labels and terminal class",
            "replace with owner_presence_word or seam_debt_word",
            7,
        ),
        StreamInterface(
            "raw_owner_count_shadow",
            (),
            (),
            False,
            (),
            0,
            "which three owners and whether they are discharged",
            "replace count by owner_union",
            8,
        ),
    )
    ordered_interfaces = tuple(
        sorted(interfaces, key=lambda item: (-item.score(), item.tie_rank))
    )

    return {
        "seam": seam,
        "transport": transport,
        "lower_bypass": lower_bypass,
        "branch_boundary": branch_boundary,
        "branch_by_branch": branch_by_branch,
        "branch_intersections": branch_intersections,
        "branch_lift": branch_lift,
        "transport_only": transport_only,
        "residual_after_transport": residual_after_transport,
        "residual_after_branch": residual_after_branch,
        "transport_plus_branch": transport_plus_branch,
        "mirror_owner_persistent": mirror_owner_persistent,
        "bypass_owner_words": bypass_owner_words,
        "lower_matches_transport": lower_matches_transport,
        "branch_word_constant": branch_word_constant,
        "branch_intersection_constant": branch_intersection_constant,
        "persistence_word": persistence_word,
        "global_flow": global_flow,
        "global_minus_bypass": global_minus_bypass,
        "stages": stages,
        "interfaces": interfaces,
        "ordered_interfaces": ordered_interfaces,
    }


def print_stage(stage: EmissionStage) -> None:
    print(
        f"{stage.name} observer={stage.observer} emitted={stage.emitted} "
        f"incoming_carry={stage.incoming_carry} outgoing_carry={stage.outgoing_carry} "
        f"stable={stage.stable}"
    )
    print(f"  certificate={stage.certificate}")
    print(f"  destroys={stage.destroyed}")
    print(f"  proof_pull={stage.proof_pull}")


def main() -> None:
    ctx = build_context()
    stream = build_stream(ctx)
    stages = stream["stages"]
    ordered_interfaces = stream["ordered_interfaces"]
    carry_widths = tuple(len(stage.outgoing_carry) for stage in stages)
    emitted_widths = tuple(len(stage.emitted) for stage in stages)
    strict_after_load = all(
        len(stages[index].outgoing_carry) <= len(stages[index - 1].outgoing_carry)
        for index in range(1, len(stages))
    )
    terminal_residual_closed = False

    print("HYP-3523 RANDOM031 SPIGOT OWNER-CARRY AUDIT")
    print("status=EVIDENCE / streaming proof-interface scout; not an LRC14 proof")
    print("row=random_covering_031")
    print("source_inspiration=https://en.wikipedia.org/wiki/Spigot_algorithm")
    print(
        "spigot_reading=emit only locally stable owner certificates; keep the "
        "unemitted owner-boundary debt as bounded carry."
    )
    print()

    print("## Input Certificates")
    print(f"seam_owner_word={stream['seam']}")
    print(f"transport_owner_word={stream['transport']}")
    print(f"lower_bypass_owner_word={stream['lower_bypass']}")
    print(f"branch_boundary_owner_word={stream['branch_boundary']}")
    print(f"branch_boundary_by_branch={stream['branch_by_branch']}")
    print(f"branch_boundary_intersections={stream['branch_intersections']}")
    print(f"branch_lift={stream['branch_lift']}")
    print(f"transport_only={stream['transport_only']}")
    print(f"residual_after_transport={stream['residual_after_transport']}")
    print(f"residual_after_branch={stream['residual_after_branch']}")
    print(f"owner_boundary_persistence_word={stream['persistence_word']}")
    print(f"mirror_owner_persistent={stream['mirror_owner_persistent']}")
    print(f"bypass_owner_word_hist={compact_counter(stream['bypass_owner_words'])}")
    print(f"lower_bypass_matches_transport={stream['lower_matches_transport']}")
    print(f"branch_word_constant={stream['branch_word_constant']}")
    print(f"branch_intersection_constant={stream['branch_intersection_constant']}")
    print()

    print("## Spigot Emission Schedule")
    for stage in stages:
        print_stage(stage)
    print(f"carry_widths={carry_widths}")
    print(f"emitted_widths={emitted_widths}")
    print(f"carry_monotone_after_load={strict_after_load}")
    print(f"bounded_carry_max_after_transport={max(carry_widths[1:])}")
    print(f"terminal_residual_closed={terminal_residual_closed}")
    print(
        "reading=HYP-3522 converts the HYP-3520 four-owner carry into a "
        "two-owner residual carry; it does not close that residual by itself."
    )
    print()

    print("## Binary Relations")
    print(
        "emit_precedence="
        "S0_load_forbidden_seam < S1_emit_transport_word < "
        "S2_emit_branch_bracket_lift < S3_emit_free_hole_cap_certificates < "
        "S4_defer_residual_pair"
    )
    print(
        "carry_refinement="
        f"{stream['seam']} -> {stream['residual_after_transport']} -> "
        f"{stream['residual_after_branch']}"
    )
    print(
        "safe_forget_relation="
        + str(tuple(ctx["persistence"]["clean_quotients"]))
    )
    print(
        "forbidden_forget_relation="
        + str(tuple(ctx["persistence"]["lossy_quotients"]))
    )
    print(
        "contaminated_global_flow_minus_bypass="
        + str(stream["global_minus_bypass"])
    )
    print()

    print("## Stream Interface Tournament")
    print("interface | score | stable | emitted | residual | contaminant | destroys | repair")
    for interface in ordered_interfaces:
        print(
            f"{interface.name} | {interface.score()} | {interface.stable} | "
            f"{interface.emitted} | {interface.residual} | {interface.contaminant} | "
            f"{interface.destroys} | {interface.repair}"
        )
    print(
        "score_hist="
        + str(compact_counter(Counter(interface.score() for interface in stream["interfaces"])))
    )
    print(f"directed_3cycles={directed_3cycle_count(ordered_interfaces)}")
    print(
        "hamiltonian_path="
        + " -> ".join(interface.name for interface in ordered_interfaces)
    )
    print()

    print("## Proof Pull")
    print(
        "P1: State the random031 pure-bypass proof as a bounded owner-carry "
        "transducer: seam payload in, transport digit out, branch-lift digit "
        "out, residual carry (45,173) retained."
    )
    print(
        "P2: Treat the spigot safety predicate as quotient constancy: a digit "
        "may be emitted only when all observers in its fiber agree."
    )
    print(
        "P3: Raw owner count, endpoint rank, component size, and mirror closure "
        "are illegal emitters because HYP-3520 shows they mix terminal classes."
    )
    print(
        "P4: The next terminal lemma should be phrased as 'no legal residual "
        "two-owner carry survives the route sidecar R', not as a search for "
        "more bypass hits."
    )
    print()

    print("## Assumption Challenge")
    print(
        "Considered vertices: runners, arcs, hard gates, bypass cells, branch "
        "boundaries, owner labels, free-hole packets, u-fibers, carry states, "
        "quotient observers, and proof obligations."
    )
    print(
        "Chosen vertices are streaming proof interfaces.  This preserves the "
        "random031 terminal owner predicate and intentionally destroys raw "
        "phase order only after the emitted owner word and residual carry have "
        "been recorded."
    )
    print(
        "Challenged assumption: the pure bypass is not a wall-crossing event "
        "that must carry all seven owners.  It is a finite-state emitter whose "
        "stable output is (23,93,113)+(147,169), with carry (45,173)."
    )


if __name__ == "__main__":
    main()
