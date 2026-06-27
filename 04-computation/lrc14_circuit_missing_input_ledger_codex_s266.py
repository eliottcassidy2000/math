#!/usr/bin/env python3
"""HYP-3116 scout: circuit lower-bound and missing-input ledger.

This script mines a curated set of older repo carriers and turns "circuit
complexity" into a proof-route audit for LRC14.

The central distinction:
  * data circuits compute tournament invariants such as H(T) from tiling bits;
  * proof circuits certify LRC14 by routing residual packets through legal exits.

The output is not a proof.  It is a regression guard: any low-depth shortcut
must either supply a minimal certificate minterm or name the first missing input
and the sidecar that repairs it.

Tournament Analysis uses proof-gate families as vertices, not runners, Boolean
gates, or raw scalar values.
"""

from __future__ import annotations

import itertools
import os
from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import lru_cache


ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


INPUTS = (
    "direct_witness",
    "ap_gw_boundary",
    "finite_address",
    "observer_gluing",
    "root_ear_sidecar",
    "relation_lattice",
    "component_bound",
    "cocycle_exactness",
    "state_lift",
    "endpoint_owner",
    "uniformity",
    "pde_weak_form",
)


SIDE_EXIT_INPUTS = (
    "root_ear_sidecar",
    "relation_lattice",
    "component_bound",
    "cocycle_exactness",
    "state_lift",
    "pde_weak_form",
)


@dataclass(frozen=True)
class LegacyCarrier:
    name: str
    paths: tuple[str, ...]
    supplies: frozenset[str]
    source_note: str
    proof_use: str
    destroyed_if_scalarized: str


@dataclass(frozen=True)
class Shortcut:
    name: str
    supplied: frozenset[str]
    required: frozenset[str]
    source_note: str


@dataclass(frozen=True)
class GateVertex:
    name: str
    retains: frozenset[str]
    proof_ready: int
    uniform: int
    reconstructs: int
    data_circuit: int
    old_source_count: int


LEGACY = (
    LegacyCarrier(
        "staircase_walsh_carry_circuit",
        ("04-computation/boolean_circuit_info_s20am.py",),
        frozenset({"uniformity"}),
        "fixed-path tournament tilings compute H through Walsh/carry layers",
        "data-circuit analogy for low-degree interactions and carry debt",
        "endpoint owners, finite address, and observer gluing",
    ),
    LegacyCarrier(
        "haar_zipper_cocycle",
        (
            "04-computation/lrc14_haar_zipper_cocycle_codex_s166.py",
            "05-knowledge/results/lrc14_haar_zipper_cocycle_codex_s166.out",
        ),
        frozenset({"cocycle_exactness", "endpoint_owner"}),
        "fixed-margin Haar zeta plus margins is a complete local coordinate",
        "local quotient repair for product/Haar shortcuts",
        "route legality and finite-address closure",
    ),
    LegacyCarrier(
        "cocycle_normal_form",
        (
            "04-computation/lrc14_cocycle_normal_form_atlas_codex_s167.py",
            "05-knowledge/results/lrc14_cocycle_normal_form_atlas_codex_s167.out",
        ),
        frozenset({"cocycle_exactness", "state_lift", "endpoint_owner"}),
        "packet cochains descend, become coboundaries, or name residual debt",
        "proof-circuit sidecar for exactness before quotienting",
        "root/ear and relation-lattice geometry",
    ),
    LegacyCarrier(
        "observer_extension_cut_payload",
        (
            "05-knowledge/hypotheses/HYP-3054-observer-extension-cut-calculus.md",
            "00-navigation/LRC-LENS-MAP.md",
        ),
        frozenset({"observer_gluing", "endpoint_owner", "uniformity"}),
        "observer/cut payload is the boundary coordinate for quotient operations",
        "declares destroyed coordinates and required sidecars",
        "the finite-address packet itself",
    ),
    LegacyCarrier(
        "route_state_closure_median",
        (
            "05-knowledge/hypotheses/HYP-3074-route-state-closure-median-interface-lrc14.md",
            "00-navigation/LRC-TECHNIQUE-INDEX.md",
        ),
        frozenset({"uniformity", "cocycle_exactness", "endpoint_owner"}),
        "legal sidecar closure comes before route-triple centers",
        "prevents raw route-label clique centers from acting as proof",
        "direct finite-address and root-locus evidence",
    ),
    LegacyCarrier(
        "finite_address_branch_closure",
        (
            "05-knowledge/hypotheses/HYP-3083-lrc14-hurwitz-finite-address-branch-closure-spine.md",
            "00-navigation/LRC-LENS-MAP.md",
        ),
        frozenset({"finite_address", "endpoint_owner", "state_lift"}),
        "finite-address branch closure is the current theorem spine",
        "packs q-cusp, level-7, branch, owner, and terminal exits",
        "observer-gluing side of the certificate",
    ),
    LegacyCarrier(
        "observer_gluing_ledger",
        (
            "04-computation/lrc14_observer_gluing_ledger_codex_s258.py",
            "05-knowledge/results/lrc14_two_frontier_gluing_s258.out",
        ),
        frozenset({"observer_gluing", "component_bound", "finite_address"}),
        "direct components, denominator nets, and pair/Pascal charts glue",
        "terminal exit for residual rows after chart compatibility",
        "root/ear and uniformity data if used alone",
    ),
    LegacyCarrier(
        "lee_yang_savitch_ear_lattice",
        (
            "04-computation/lee_yang_savitch_ear_lattice_extremality_codex_s262.py",
            "05-knowledge/results/lee_yang_savitch_ear_lattice_extremality_codex_s262.out",
        ),
        frozenset({"root_ear_sidecar", "relation_lattice", "uniformity"}),
        "root curves, Bravais shape, Savitch midpoint, and ears form two maps",
        "compresses reachability only when sidecar traps are named",
        "finite-address and observer-gluing exits",
    ),
    LegacyCarrier(
        "lee_yang_ear_payload",
        (
            "04-computation/lrc_lee_yang_ear_payload_codex_s262.py",
            "05-knowledge/results/lrc_lee_yang_ear_payload_codex_s262.out",
        ),
        frozenset({"root_ear_sidecar", "endpoint_owner"}),
        "one-runner payload reconstructs PGF root motion",
        "repairs scalar root-count shortcuts with extension data",
        "finite-address and observer-gluing proof closure",
    ),
    LegacyCarrier(
        "minkowski_circuit_ising_quintic",
        (
            "04-computation/lrc_minkowski_circuit_ising_quintic_bridge_codex_20260627.py",
            "05-knowledge/results/lrc_minkowski_circuit_ising_quintic_bridge_codex_20260627.out",
        ),
        frozenset({"relation_lattice", "root_ear_sidecar", "uniformity"}),
        "HYP-3115 anchored-bank scout with a finite fitted circuit warning",
        "tests sidecars and exposes nonuniform one-literal classifiers",
        "observer-gluing, finite address, and endpoint owners",
    ),
    LegacyCarrier(
        "pde_weak_form_packet",
        (
            "04-computation/lrc14_minkowski_circuit_ising_pde_compiler_codex_s264.py",
            "05-knowledge/results/lrc14_minkowski_circuit_ising_pde_compiler_codex_s264.out",
        ),
        frozenset({"pde_weak_form", "component_bound", "observer_gluing"}),
        "weak forms retain mass, stiffness, boundary owners, and zero modes",
        "operator-side guardrail for component/gap proofs",
        "finite address and root/ear sidecars",
    ),
)


SHORTCUTS = (
    Shortcut(
        "raw_p0_scalar",
        frozenset(),
        frozenset({"finite_address", "observer_gluing", "root_ear_sidecar", "endpoint_owner", "uniformity"}),
        "single miss-mass value",
    ),
    Shortcut(
        "one_literal_apex7_threshold",
        frozenset({"root_ear_sidecar"}),
        frozenset({"root_ear_sidecar", "finite_address", "observer_gluing", "endpoint_owner", "uniformity"}),
        "HYP-3115 finite-bank rule apex7_error <= 5",
    ),
    Shortcut(
        "raw_H_gap_transfer",
        frozenset({"state_lift"}),
        frozenset({"state_lift", "observer_gluing", "finite_address", "uniformity"}),
        "H=7/H=21 contradiction without faithful transfer sidecar",
    ),
    Shortcut(
        "raw_minkowski_volume",
        frozenset({"relation_lattice"}),
        frozenset({"relation_lattice", "finite_address", "observer_gluing", "root_ear_sidecar", "uniformity"}),
        "convex-body pressure without packet decoding",
    ),
    Shortcut(
        "raw_ising_energy",
        frozenset({"root_ear_sidecar"}),
        frozenset({"root_ear_sidecar", "observer_gluing", "finite_address", "endpoint_owner"}),
        "domain-wall energy without labelled root-collision exits",
    ),
    Shortcut(
        "raw_demoivre_residual",
        frozenset(),
        frozenset({"finite_address", "observer_gluing", "root_ear_sidecar", "uniformity"}),
        "stationary quintic residual without branch orbit address",
    ),
    Shortcut(
        "raw_walsh_low_degree",
        frozenset({"uniformity"}),
        frozenset({"uniformity", "cocycle_exactness", "endpoint_owner", "finite_address"}),
        "low-degree data circuit treated as proof circuit",
    ),
    Shortcut(
        "raw_direct_component_count",
        frozenset({"component_bound"}),
        frozenset({"component_bound", "observer_gluing", "finite_address", "endpoint_owner"}),
        "component count without observer chart and owner data",
    ),
    Shortcut(
        "raw_pair_pascal_shadow",
        frozenset(),
        frozenset({"cocycle_exactness", "endpoint_owner", "observer_gluing", "finite_address"}),
        "pair-mass or Pascal shadow without scissors/gluing",
    ),
    Shortcut(
        "raw_automaton_language",
        frozenset({"uniformity"}),
        frozenset({"uniformity", "finite_address", "endpoint_owner", "state_lift"}),
        "Moser/fibbinary language without packet address",
    ),
)


GATE_VERTICES = (
    GateVertex("finite_address_branch_closure", frozenset({"finite_address", "endpoint_owner", "state_lift"}), 5, 4, 3, 0, 2),
    GateVertex("observer_gluing_certificate", frozenset({"observer_gluing", "component_bound", "finite_address"}), 5, 4, 3, 0, 2),
    GateVertex("root_ear_payload_ledger", frozenset({"root_ear_sidecar", "endpoint_owner", "uniformity"}), 4, 4, 4, 0, 2),
    GateVertex("cocycle_exactness_ledger", frozenset({"cocycle_exactness", "endpoint_owner", "state_lift"}), 4, 3, 4, 0, 3),
    GateVertex("route_center_median_ledger", frozenset({"uniformity", "cocycle_exactness", "endpoint_owner"}), 3, 5, 3, 0, 2),
    GateVertex("pde_weak_form_operator", frozenset({"pde_weak_form", "component_bound", "observer_gluing"}), 3, 3, 3, 0, 1),
    GateVertex("relation_lattice_minkowski", frozenset({"relation_lattice", "root_ear_sidecar", "uniformity"}), 3, 3, 2, 0, 3),
    GateVertex("staircase_walsh_carry_data_circuit", frozenset({"uniformity"}), 1, 4, 1, 1, 1),
    GateVertex("finite_bank_threshold_signal", frozenset({"root_ear_sidecar"}), 1, 0, 0, 0, 1),
    GateVertex("raw_scalar_p0", frozenset(), 0, 0, 0, 0, 0),
)


def repo_path(path: str) -> str:
    return os.path.join(ROOT, path)


def file_text(path: str) -> str:
    try:
        with open(repo_path(path), "r", encoding="utf-8", errors="replace") as handle:
            return handle.read()
    except FileNotFoundError:
        return ""


KEYWORDS = (
    "circuit",
    "Walsh",
    "proof",
    "finite-address",
    "observer",
    "Savitch",
    "ear",
    "cocycle",
    "median",
    "uniform",
    "sidecar",
    "component",
)


def keyword_counts(carrier: LegacyCarrier) -> Counter[str]:
    counts: Counter[str] = Counter()
    for path in carrier.paths:
        text = file_text(path)
        low = text.lower()
        for key in KEYWORDS:
            counts[key] += low.count(key.lower())
    return counts


def eval_formula(bits: dict[str, bool]) -> bool:
    direct = bits["direct_witness"]
    boundary = bits["ap_gw_boundary"]
    observer_exit = (
        bits["finite_address"]
        and bits["observer_gluing"]
        and bits["endpoint_owner"]
        and bits["uniformity"]
        and any(bits[x] for x in SIDE_EXIT_INPUTS)
    )
    return direct or boundary or observer_exit


def minimal_certificates() -> list[tuple[str, ...]]:
    certs: list[tuple[str, ...]] = []
    for r in range(1, len(INPUTS) + 1):
        for subset in itertools.combinations(INPUTS, r):
            bits = {name: name in subset for name in INPUTS}
            if not eval_formula(bits):
                continue
            if any(set(prev).issubset(subset) for prev in certs):
                continue
            certs.append(subset)
    return certs


def essential_inputs() -> list[str]:
    essentials = []
    for name in INPUTS:
        witness = False
        for values in itertools.product([False, True], repeat=len(INPUTS)):
            bits = dict(zip(INPUTS, values))
            other = dict(bits)
            other[name] = not other[name]
            if eval_formula(bits) != eval_formula(other):
                witness = True
                break
        if witness:
            essentials.append(name)
    return essentials


def gate_count_and_depth() -> tuple[int, int]:
    # OR(direct, boundary, AND(finite, observer, endpoint, uniform, OR(six side exits)))
    # Count non-input gates only.
    side_or_size = 1
    main_and_size = 1
    output_or_size = 1
    gate_count = side_or_size + main_and_size + output_or_size
    depth = 3
    return gate_count, depth


def missing(shortcut: Shortcut) -> frozenset[str]:
    return shortcut.required - shortcut.supplied


def repair_cover(shortcut: Shortcut) -> tuple[str, ...]:
    need = set(missing(shortcut))
    if not need:
        return ()
    best: tuple[str, ...] | None = None
    for r in range(1, len(LEGACY) + 1):
        for combo in itertools.combinations(LEGACY, r):
            supplied = set().union(*(c.supplies for c in combo))
            if need <= supplied:
                names = tuple(c.name for c in combo)
                best = names
                break
        if best is not None:
            return best
    return ()


def score(v: GateVertex) -> tuple[int, int, int, int, int, int]:
    return (
        len(v.retains & frozenset(INPUTS)),
        v.proof_ready,
        v.uniform,
        v.reconstructs,
        v.old_source_count,
        -v.data_circuit,
    )


def tournament_edges(vertices: tuple[GateVertex, ...]) -> dict[tuple[str, str], str]:
    edges: dict[tuple[str, str], str] = {}
    for a, b in itertools.combinations(vertices, 2):
        sa = score(a)
        sb = score(b)
        if sa == sb:
            winner, loser = (a, b) if a.name < b.name else (b, a)
        elif sa > sb:
            winner, loser = a, b
        else:
            winner, loser = b, a
        edges[(winner.name, loser.name)] = ">"
    return edges


def beats(edges: dict[tuple[str, str], str], a: str, b: str) -> bool:
    return (a, b) in edges


def score_hist(vertices: tuple[GateVertex, ...], edges: dict[tuple[str, str], str]) -> Counter[int]:
    hist: Counter[int] = Counter()
    for v in vertices:
        out = sum(1 for u in vertices if u.name != v.name and beats(edges, v.name, u.name))
        hist[out] += 1
    return hist


def directed_3cycles(vertices: tuple[GateVertex, ...], edges: dict[tuple[str, str], str]) -> int:
    count = 0
    names = [v.name for v in vertices]
    for a, b, c in itertools.combinations(names, 3):
        if beats(edges, a, b) and beats(edges, b, c) and beats(edges, c, a):
            count += 1
        if beats(edges, a, c) and beats(edges, c, b) and beats(edges, b, a):
            count += 1
    return count


def scc_sizes(vertices: tuple[GateVertex, ...], edges: dict[tuple[str, str], str]) -> list[int]:
    names = [v.name for v in vertices]
    adj = {name: [other for other in names if other != name and beats(edges, name, other)] for name in names}

    def reach(start: str) -> set[str]:
        seen = {start}
        stack = [start]
        while stack:
            cur = stack.pop()
            for nxt in adj[cur]:
                if nxt not in seen:
                    seen.add(nxt)
                    stack.append(nxt)
        return seen

    unassigned = set(names)
    sizes = []
    while unassigned:
        start = next(iter(unassigned))
        component = {name for name in names if name in reach(start) and start in reach(name)}
        sizes.append(len(component))
        unassigned -= component
    return sorted(sizes, reverse=True)


def hamiltonian_path_count(vertices: tuple[GateVertex, ...], edges: dict[tuple[str, str], str]) -> int:
    names = [v.name for v in vertices]
    count = 0
    for perm in itertools.permutations(names):
        if all(beats(edges, perm[i], perm[i + 1]) for i in range(len(perm) - 1)):
            count += 1
    return count


def priority_path(vertices: tuple[GateVertex, ...]) -> tuple[str, ...]:
    return tuple(v.name for v in sorted(vertices, key=score, reverse=True))


def main() -> None:
    print("HYP-3116 circuit lower-bound / missing-input ledger -- codex S266")
    print()

    gate_count, depth = gate_count_and_depth()
    certs = minimal_certificates()
    essentials = essential_inputs()
    print("PROOF CIRCUIT MODEL")
    print(f"  input_count={len(INPUTS)}")
    print(f"  gate_count={gate_count}")
    print(f"  depth={depth}")
    print(f"  essential_inputs={len(essentials)}: {', '.join(essentials)}")
    print("  minimal_certificate_minterms:")
    for cert in certs:
        print("    " + " AND ".join(cert))
    print()

    print("LEGACY CARRIER KEYWORD SCAN")
    for carrier in LEGACY:
        counts = keyword_counts(carrier)
        existing = sum(1 for path in carrier.paths if os.path.exists(repo_path(path)))
        top = ", ".join(f"{k}:{v}" for k, v in counts.most_common(4))
        print(f"  {carrier.name}: files={existing}/{len(carrier.paths)}, supplies={','.join(sorted(carrier.supplies))}")
        print(f"    top_keywords={top if top else 'none'}")
        print(f"    use={carrier.proof_use}")
    print()

    print("SHORTCUT MISSING-INPUT AUDIT")
    closed = 0
    total_missing = 0
    for shortcut in SHORTCUTS:
        miss = sorted(missing(shortcut))
        total_missing += len(miss)
        cover = repair_cover(shortcut)
        if not miss:
            closed += 1
        print(f"  {shortcut.name}: supplied={','.join(sorted(shortcut.supplied)) or 'none'}")
        print(f"    missing_count={len(miss)} missing={','.join(miss) if miss else 'none'}")
        print(f"    repair_cover={', '.join(cover) if cover else 'none'}")
        print(f"    source={shortcut.source_note}")
    print(f"  closed_shortcuts={closed}/{len(SHORTCUTS)}")
    print(f"  total_missing_inputs={total_missing}")
    print()

    missing_counter: Counter[str] = Counter()
    for shortcut in SHORTCUTS:
        missing_counter.update(missing(shortcut))
    print("MISSING INPUT FREQUENCY")
    for name, count in missing_counter.most_common():
        print(f"  {name}: {count}")
    print()

    print("TOURNAMENT ANALYSIS")
    edges = tournament_edges(GATE_VERTICES)
    hist = score_hist(GATE_VERTICES, edges)
    print(f"  vertices={len(GATE_VERTICES)}")
    print(f"  score_hist={dict(sorted(hist.items()))}")
    print(f"  directed_3cycles={directed_3cycles(GATE_VERTICES, edges)}")
    print(f"  scc_sizes={scc_sizes(GATE_VERTICES, edges)}")
    print(f"  hamiltonian_path_count={hamiltonian_path_count(GATE_VERTICES, edges)}")
    print("  priority_path:")
    print("    " + " -> ".join(priority_path(GATE_VERTICES)))
    print()

    print("SYNTHESIS")
    print("  Data-circuit lesson: low Walsh degree or a one-literal threshold is a")
    print("  compressed computation, not a proof route, unless it supplies the")
    print("  proof-circuit inputs demanded by the certificate DAG.")
    print("  The most frequent missing inputs are the real bottlenecks above.")
    print("  Therefore the next LRC14 ledger should store proof_circuit_missing_input_vector")
    print("  beside root/ear, finite-address, observer-gluing, cocycle, and endpoint-owner fields.")


if __name__ == "__main__":
    main()
