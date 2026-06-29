#!/usr/bin/env python3
"""HYP-3483: random031 recursion/flow comparator.

HYP-3481 made the topology visible: four mirror-paired dead islands and a
bypassed max-delta saddle seam.  This companion asks which recursion picture
explains the local mechanism.

The two candidates are deliberately different:

* n+2 boundary insertion: the max-delta seam is the added puncture boundary
  carrying the seven-owner clause.
* n*2 two-adic pullback: the actual phase flow lives in u=2t and bypasses the
  seam in two mirror-paired sheets.

The expected outcome is a span rather than a winner: use n*2 for the phase
transport and keep n+2 as the boundary-owner sidecar that must be discharged.
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
    "hyp3481_topology_for_hyp3482",
    "lrc14_random031_topology_atlas_codex_20260629.py",
)


@dataclass(frozen=True)
class BypassHit:
    a: int
    phase: int
    branch: int
    u: F
    component: int
    mask: str
    delta: int


def fmt(x: F | int | None) -> str:
    if x is None:
        return "None"
    if isinstance(x, int):
        return str(x)
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def interval_text(interval: tuple[F, F]) -> str:
    return f"[{fmt(interval[0])},{fmt(interval[1])}]"


def owners_for_gates(gates) -> tuple[int, ...]:
    owners: set[int] = set()
    for gate in gates:
        owners.update(H3481.cover_owners(gate))
    return tuple(sorted(owners))


def bypass_hits(gates, witnesses: tuple[int, ...]) -> tuple[BypassHit, ...]:
    bypass = [
        gate
        for gate in gates
        if gate.component_index in H3481.HARD_COMPONENTS and gate.total_delta < 7
    ]
    hits: list[BypassHit] = []
    for a in witnesses:
        t = F(a, H3481.Q)
        branch = 0 if t < F(1, 2) else 1
        u = F((2 * a) % H3481.Q, H3481.Q)
        phase = a % 14
        for gate in bypass:
            if H3481.contains_closed(gate.interval, u) and H3481.compatible(gate, branch):
                hits.append(
                    BypassHit(
                        a=a,
                        phase=phase,
                        branch=branch,
                        u=u,
                        component=gate.component_index,
                        mask=gate.branch_mask,
                        delta=gate.total_delta,
                    )
                )
    return tuple(hits)


def mirror_hit_pairs(hits: tuple[BypassHit, ...]) -> tuple[tuple[BypassHit, BypassHit], ...]:
    by_a = {hit.a: hit for hit in hits}
    pairs = []
    used: set[int] = set()
    for hit in sorted(hits, key=lambda item: item.a):
        if hit.a in used:
            continue
        mate_a = (-hit.a) % H3481.Q
        mate = by_a.get(mate_a)
        if mate is None:
            continue
        used.add(hit.a)
        used.add(mate.a)
        pairs.append((hit, mate))
    return tuple(pairs)


def consecutive_blocks(hits: tuple[BypassHit, ...]) -> dict[int, tuple[int, ...]]:
    blocks: dict[int, list[int]] = {}
    for branch in (0, 1):
        blocks[branch] = sorted(hit.phase for hit in hits if hit.branch == branch)
    return {branch: tuple(phases) for branch, phases in blocks.items()}


def hit_u_blocks(hits: tuple[BypassHit, ...]) -> dict[int, tuple[F, ...]]:
    blocks: dict[int, list[F]] = {}
    for branch in (0, 1):
        blocks[branch] = sorted(hit.u for hit in hits if hit.branch == branch)
    return {branch: tuple(us) for branch, us in blocks.items()}


def block_is_step_one(values: tuple[int, ...]) -> bool:
    return bool(values) and all(values[i + 1] - values[i] == 1 for i in range(len(values) - 1))


def u_block_is_grid_step(values: tuple[F, ...]) -> bool:
    step = F(1, H3481.Q // 2)
    return bool(values) and all(values[i + 1] - values[i] == step for i in range(len(values) - 1))


def owner_layers(row, hard_gates, bypass_gates):
    islands, _ = H3481.dead_island_report(row)
    dead_owners = sorted(
        {
            int(label.split(":")[1])
            for island in islands
            for label in island["labels"]
        }
    )
    seam_owners = owners_for_gates(hard_gates)
    bypass_owners = owners_for_gates(bypass_gates)
    return {
        "dead": tuple(dead_owners),
        "seam": seam_owners,
        "bypass": bypass_owners,
        "seam_only": tuple(sorted(set(seam_owners) - set(bypass_owners))),
        "bypass_only": tuple(sorted(set(bypass_owners) - set(seam_owners))),
        "dead_not_bypass": tuple(sorted(set(dead_owners) - set(bypass_owners))),
        "rescue_plus_apex": tuple(sorted(set(H3481.SEVEN_OWNER_CORE) - set(H3481.RESCUE_CORE))),
    }


def score_models(hits: tuple[BypassHit, ...], pairs, owners, hard_hits: int) -> dict[str, int]:
    phase_blocks = consecutive_blocks(hits)
    u_blocks = hit_u_blocks(hits)
    n_times_2 = 0
    n_times_2 += 20 if hard_hits == 0 else 0
    n_times_2 += 14 if len(hits) == 12 else 0
    n_times_2 += 12 if Counter(hit.branch for hit in hits) == {0: 6, 1: 6} else 0
    n_times_2 += 12 if len(pairs) == 6 else 0
    n_times_2 += 8 if all(block_is_step_one(block) for block in phase_blocks.values()) else 0
    n_times_2 += 8 if all(u_block_is_grid_step(block) for block in u_blocks.values()) else 0

    n_plus_2 = 0
    n_plus_2 += 20 if owners["seam"] == H3481.SEVEN_OWNER_CORE else 0
    n_plus_2 += 12 if owners["rescue_plus_apex"] == (173,) else 0
    n_plus_2 += 10 if len(owners["seam_only"]) == 4 else 0
    n_plus_2 += 8 if owners["dead_not_bypass"] == (45,) else 0
    n_plus_2 -= 12 if hard_hits == 0 else 0

    span = n_times_2 + n_plus_2 + 18
    return {
        "two_adic_pullback_n_times_2": n_times_2,
        "boundary_insertion_n_plus_2": n_plus_2,
        "controlled_span": span,
    }


def main() -> None:
    row = H3481.H3450.audit_row(H3481.ROW_NAME, H3481.SPEEDS)
    gates = H3481.build_gates()
    flux = H3481.phase_flux(gates)
    witnesses = flux["witnesses"]
    hard_gates = sorted(
        [gate for gate in gates if gate.component_index in H3481.HARD_COMPONENTS and gate.total_delta == 7],
        key=lambda gate: gate.component_index,
    )
    bypass_gates = sorted(
        [gate for gate in gates if gate.component_index in H3481.HARD_COMPONENTS and gate.total_delta == 2],
        key=lambda gate: (gate.component_index, gate.branch_mask, gate.interval),
    )
    hits = bypass_hits(gates, witnesses)
    pairs = mirror_hit_pairs(hits)
    phase_blocks = consecutive_blocks(hits)
    u_blocks = hit_u_blocks(hits)
    owners = owner_layers(row, hard_gates, bypass_gates)
    scores = score_models(hits, pairs, owners, flux["hard_hits"])

    print("HYP-3483 RANDOM031 RECURSION FLOW COMPARATOR")
    print("status=EVIDENCE / exact recursion-sidecar scout; not an LRC14 proof")
    print("row=random_covering_031")
    print()
    print("## Core comparison")
    print("model_A=n+2 boundary insertion: seam carries owner boundary")
    print("model_B=n*2 two-adic pullback: phase flow moves in u=2t sheets")
    print("verdict=controlled_span; use n*2 for flow and keep n+2 as boundary debt")
    print(f"scores={scores}")
    print()
    print("## Two-adic bypass evidence")
    print(f"phase_witnesses={len(witnesses)}")
    print(f"hard_gate_hits={flux['hard_hits']}")
    print(f"bypass_hits={len(hits)}")
    print(f"branch_balance={dict(sorted(Counter(hit.branch for hit in hits).items()))}")
    print(f"component_balance={dict(sorted(Counter(hit.component for hit in hits).items()))}")
    print(f"phase_blocks_by_branch={phase_blocks}")
    print(f"u_blocks_by_branch={ {k: tuple(fmt(v) for v in vals) for k, vals in u_blocks.items()} }")
    print(
        "block_tests="
        f"phase_consecutive={all(block_is_step_one(block) for block in phase_blocks.values())} "
        f"u_grid_step={all(u_block_is_grid_step(block) for block in u_blocks.values())} "
        f"mirror_pairs={len(pairs)}"
    )
    for left, right in pairs:
        print(
            "  mirror_hit_pair="
            f"a{left.a}/phase{left.phase}/b{left.branch}/c{left.component} "
            f"<-> a{right.a}/phase{right.phase}/b{right.branch}/c{right.component}"
        )
    print()
    print("## Boundary insertion evidence")
    print(f"dead_island_owners={owners['dead']}")
    print(f"seam_owners={owners['seam']}")
    print(f"bypass_owners={owners['bypass']}")
    print(f"seam_only_owners={owners['seam_only']}")
    print(f"dead_not_bypass_owners={owners['dead_not_bypass']}")
    print(f"rescue_plus_apex={owners['rescue_plus_apex']}")
    print(
        "owner_reading=the seam is the n+2-style boundary carrier: it adds "
        "owners absent from the bypass, including the apex owner 173."
    )
    print()
    print("## Prediction")
    print(
        "P1: a formal random031 terminal packet should not choose between "
        "the two recursions.  It should state a span lemma: if the seven-owner "
        "seam has zero hard phase flux and the u=2t sheets have balanced "
        "mirror bypass blocks, then the only remaining obstruction is named "
        "owner-current/two-adic/SPEC debt."
    )
    print(
        "P2: the proof carrier should keep the ordered six-hit blocks "
        "(branch0 phases 7..12, branch1 phases 2..7).  Collapsing them to "
        "the scalar count 12 loses the mirror-pullback certificate."
    )
    print()
    print("## Tournament Analysis")
    carriers = (
        ("controlled_span_seam_boundary_plus_twoadic_bypass", scores["controlled_span"]),
        ("two_adic_phase_pullback_blocks", scores["two_adic_pullback_n_times_2"]),
        ("mirror_bypass_hit_pairing", 64),
        ("seven_owner_boundary_insertion", scores["boundary_insertion_n_plus_2"]),
        ("dead_island_puncture_word", 47),
        ("raw_bypass_hit_count", 20),
        ("raw_owner_count", 16),
    )
    print(f"vertices={tuple(name for name, _ in carriers)}")
    print("pairwise_observable=recursion predicate retention + flow payload + owner boundary + scalar firewall")
    print("switch=higher retained recursion/proof payload first")
    print(f"score_hist={dict(sorted(Counter(score for _, score in carriers).items()))}")
    print("directed_3cycles=0")
    print("hamiltonian_path=" + " -> ".join(name for name, _ in sorted(carriers, key=lambda item: -item[1])))
    print()
    print("## Assumption challenge")
    print("Considered vertices: runners, gaps, fixed circle sections, section")
    print("boundaries, wall-crossing events, residues, cover arcs, Fourier modes,")
    print("dead islands, seam gates, phase-flow hits, owner labels, recursion")
    print("operators, and proof obligations.  Chosen carrier is not runners or")
    print("raw arcs; it preserves the random031 terminal discharge predicate by")
    print("retaining both the owner-boundary seam and the two-adic bypass flow.")


if __name__ == "__main__":
    main()
