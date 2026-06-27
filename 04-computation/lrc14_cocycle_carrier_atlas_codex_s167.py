#!/usr/bin/env python3
"""Cocycle carrier atlas for LRC14 proof technology.

This is a synthesis script, not a proof search.  It records the repo's current
LRC14 proof routes as cocycle carriers: each carrier has a base, a fiber, a
local signed defect, a closedness/gluing law, and exits that make quotienting
legal.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "05-knowledge/results/lrc14_cocycle_carrier_atlas_codex_s167.out"


PROOF_FEATURES = [
    "predicate",
    "fiber",
    "local",
    "closedness",
    "coboundary",
    "boundary_atom",
    "descent",
    "named_residual",
    "formal",
    "lrc14",
    "cross_domain",
]

EXPLORATION_WEIGHTS = {
    "predicate": 2,
    "fiber": 2,
    "local": 2,
    "closedness": 3,
    "coboundary": 1,
    "boundary_atom": 1,
    "descent": 2,
    "named_residual": 2,
    "formal": 1,
    "lrc14": 2,
    "cross_domain": 4,
}


@dataclass(frozen=True)
class Carrier:
    name: str
    base: str
    vertices: str
    cochain: str
    closedness: str
    exits: str
    lrc_pull: str
    destroys: str
    anchors: str
    scores: Dict[str, int]

    @property
    def proof_score(self) -> int:
        return sum(self.scores[f] for f in PROOF_FEATURES)

    @property
    def exploration_score(self) -> int:
        return sum(EXPLORATION_WEIGHTS[f] * self.scores[f] for f in PROOF_FEATURES)


CARRIERS: List[Carrier] = [
    Carrier(
        "labelled_packet_master_cocycle",
        "packet base P(S)",
        "qdiv/Farey/Haar/endpoint/state/certificate fibers",
        "total lost-coordinate obstruction omega_P",
        "all local route defects must agree on overlaps",
        "strict witness, AP/GW atom, state-lift residual, or exact certificate",
        "the unifying theorem target for all LRC14 quotients",
        "nothing if all coordinates stay labelled; too large for direct proof",
        "HYP-2976, HYP-2987, HYP-2990, LRC-TECHNIQUE-INDEX",
        dict(predicate=3, fiber=3, local=2, closedness=3, coboundary=2,
             boundary_atom=3, descent=3, named_residual=3, formal=2, lrc14=3,
             cross_domain=2),
    ),
    Carrier(
        "haar_zipper_zeta",
        "2 x 2 fixed-margin packet table",
        "row/column margin fibers",
        "zeta=T00-T01-T10+T11",
        "fixed-margin moves change zeta in gcd-4 steps",
        "color cancellation, boundary stop, handoff, descent, F7 debt",
        "prevents row/column margins from erasing mixed Haar sign",
        "global endpoint owners and exact-period labels",
        "HYP-2991, HYP-2989, HYP-2595",
        dict(predicate=3, fiber=3, local=3, closedness=3, coboundary=2,
             boundary_atom=2, descent=2, named_residual=3, formal=3, lrc14=3,
             cross_domain=2),
    ),
    Carrier(
        "endpoint_credit_winding",
        "danger-arc endpoint transition graph",
        "active endpoint-owner arcs",
        "K(a,b)=14rs(R(a)-L(b))",
        "unit-winding cycle is exactly an open cover",
        "potential Phi, AP/GW closed cycle, or K33/state lift",
        "turns counterexamples into positive endpoint circulations",
        "raw Haar mass and runner identity away from endpoint graph",
        "HYP-2970, HYP-2965, HYP-2986",
        dict(predicate=3, fiber=3, local=3, closedness=3, coboundary=3,
             boundary_atom=3, descent=1, named_residual=3, formal=2, lrc14=3,
             cross_domain=2),
    ),
    Carrier(
        "carry_crt_cocycle",
        "Res_27 lift tower v=r+27k",
        "carry fibers and owner certificates",
        "k mod 14 controls parity, apex divisibility, and pair sums",
        "floor lift requires globally coherent carry alignment",
        "owner/Cprime repair, apex tax, scalar lift, or strictness",
        "makes parity, multiples of 14, and pair pinches one object",
        "continuous endpoint order unless owners are reattached",
        "HYP-2230, HYP-2241, HYP-2166",
        dict(predicate=3, fiber=3, local=3, closedness=2, coboundary=2,
             boundary_atom=2, descent=3, named_residual=2, formal=2, lrc14=3,
             cross_domain=2),
    ),
    Carrier(
        "owner_deletion_derivative",
        "visible Res_27/gcd shell quotient",
        "deleted-speed cards paired with private-owner flags",
        "delta_owner(v)=obligations uncovered by deleting v",
        "nonzero local carry must change paired owner derivative",
        "private-owner bit, carry support, or strict tax",
        "repairs visible quotient leaks without retaining full carry",
        "phase order and raw speed identity",
        "HYP-2241, HYP-2237, HYP-2230",
        dict(predicate=3, fiber=3, local=3, closedness=2, coboundary=2,
             boundary_atom=1, descent=3, named_residual=2, formal=2, lrc14=3,
             cross_domain=2),
    ),
    Carrier(
        "ramanujan_exact_period_trace",
        "exact denominator q and primitive q-th roots",
        "primitive phase packets",
        "c_q(a-b) as exact-period character trace",
        "orthogonality closes over exact-period fibers",
        "primitive projector, endpoint-owner repair, or squarefree residual",
        "repairs divisor and squarefree quotients with period labels",
        "endpoint ownership and strict-open topology if used alone",
        "HYP-2979, HYP-2978",
        dict(predicate=2, fiber=3, local=2, closedness=3, coboundary=3,
             boundary_atom=2, descent=2, named_residual=2, formal=3, lrc14=3,
             cross_domain=3),
    ),
    Carrier(
        "fejer_toeplitz_moment",
        "Fourier modes of C_S(t)-1",
        "Toeplitz moment cones",
        "quadratic Fejer form v* M v",
        "nonnegative cover forces PSD moments",
        "negative quadratic certificate or interval proof payload",
        "dual-certificate route for positive-open packets",
        "packet family unless center/degree/atom bank are retained",
        "HYP-2974, HYP-2981",
        dict(predicate=3, fiber=2, local=2, closedness=3, coboundary=3,
             boundary_atom=2, descent=2, named_residual=2, formal=3, lrc14=3,
             cross_domain=2),
    ),
    Carrier(
        "boundary_moment_multichart",
        "exact-period denominator charts",
        "boundary-moment feasible regions",
        "chart-to-chart moment defect",
        "one covered chart must glue with other charts before obstruction",
        "positive Haar-open chart, AP/GW equality, or residual chart debt",
        "prevents one all-covered denominator chart from being overread",
        "single-chart scalar obstruction",
        "HYP-2969, HYP-2971",
        dict(predicate=2, fiber=3, local=2, closedness=3, coboundary=2,
             boundary_atom=2, descent=2, named_residual=3, formal=2, lrc14=3,
             cross_domain=2),
    ),
    Carrier(
        "product_rule_tiling",
        "Haar rectangles / fixed-path Walsh staircase tiles",
        "2D address-retained product fibers",
        "product class: zero, atom, owner strip, handoff, descent, residual",
        "multiplication is exact before strip-count scalarization",
        "orthogonality, owner reconstruction, dual annihilation, or residual",
        "makes tiling shadows legal only after product-rule descent",
        "endpoint owners if product classes are not packet-labelled",
        "HYP-2989, HYP-2988, THM-351",
        dict(predicate=2, fiber=3, local=3, closedness=3, coboundary=2,
             boundary_atom=2, descent=3, named_residual=2, formal=3, lrc14=3,
             cross_domain=3),
    ),
    Carrier(
        "farey_k33_determinant",
        "exact M=p/q Farey node",
        "unit-excess and factor-incidence fibers",
        "determinant/excess e=14p-q and K_{p,q} incidence depth",
        "branch labels glue only with exact scale retained",
        "p=1 AP/GW, p=2 C27, p>=3 K33/Fejer/state lift",
        "keeps product/mediant analogies honest",
        "endpoint geometry if product value is used alone",
        "HYP-2932, HYP-2934, HYP-2945",
        dict(predicate=2, fiber=3, local=2, closedness=2, coboundary=1,
             boundary_atom=2, descent=3, named_residual=3, formal=2, lrc14=3,
             cross_domain=3),
    ),
    Carrier(
        "c27_unital_transfer",
        "C=27 antipodal shell / q=3 unital block",
        "hole/double transfer fibers",
        "local pair-completion transfer defect",
        "global row must glue branch-local completions",
        "petal, GW transfer, K33 branch, or state-lift debt",
        "explains hidden q=3 transfer not visible to scalar Haar mass",
        "global realizability if unital chart is isolated",
        "HYP-2937, HYP-2942, HYP-2940",
        dict(predicate=2, fiber=3, local=2, closedness=2, coboundary=1,
             boundary_atom=3, descent=3, named_residual=3, formal=1, lrc14=3,
             cross_domain=3),
    ),
    Carrier(
        "root_packet_boundary",
        "type-A root-sign chain complex",
        "open walks and closed root packets",
        "sum of roots = endpoint boundary; closed packets have zero boundary",
        "cycles close exactly when endpoint charge vanishes",
        "packet compatibility, incidence rank, or homology residue",
        "provides the cleanest cochain language for tournaments/OCF",
        "LRC arithmetic labels unless imported as proof carrier",
        "RootPackets.lean, HYP-1814, path-homology notes",
        dict(predicate=1, fiber=2, local=3, closedness=3, coboundary=3,
             boundary_atom=1, descent=2, named_residual=2, formal=3, lrc14=1,
             cross_domain=3),
    ),
    Carrier(
        "metagraph_transfer_cocycle",
        "operation/metagraph state graph",
        "class-to-class extension strips",
        "successor defect or transfer-matrix curvature",
        "walks glue when transfer defects telescope",
        "state-vector refinement, support residue, or SCC debt",
        "turns recursion/mode changes into path-integral constraints",
        "exact LRC endpoint arithmetic unless labels travel with edges",
        "HYP-1835, metagraph transfer-chain notes",
        dict(predicate=1, fiber=2, local=2, closedness=2, coboundary=2,
             boundary_atom=1, descent=3, named_residual=2, formal=2, lrc14=2,
             cross_domain=3),
    ),
    Carrier(
        "sequence_shadow_difference",
        "n+2 / n*2 / companion sequence recursions",
        "fixed, merged, nonfixed, q-shadow sequence fibers",
        "finite difference under add/multiply mode changes",
        "companion shadows must commute before next-term scalarization",
        "state split, derivative repair, or discarded analogy",
        "keeps series recurrences as transport data rather than numerology",
        "LRC predicate if no packet/fiber labels are attached",
        "S633, sequence-shadow reflections, LRC operation-metagraph work",
        dict(predicate=1, fiber=2, local=2, closedness=2, coboundary=2,
             boundary_atom=0, descent=3, named_residual=1, formal=1, lrc14=1,
             cross_domain=3),
    ),
    Carrier(
        "octahedral_hodge_current",
        "octahedral/K4 face-current support",
        "curl/divergence current fibers",
        "discrete curl/divergence defect",
        "closed current requires zero boundary after tail labels",
        "Hodge split, tail repair, or residual current",
        "cross-domain model for state-lift F7 as unpaired current",
        "LRC exact-period and endpoint labels unless explicitly attached",
        "HYP-2887, HYP-2990",
        dict(predicate=1, fiber=2, local=2, closedness=3, coboundary=2,
             boundary_atom=1, descent=2, named_residual=3, formal=2, lrc14=1,
             cross_domain=3),
    ),
    Carrier(
        "ocf_activity_coimage",
        "odd-cycle compatibility graph",
        "activity/coimage packet fibers",
        "hard-core activity residue",
        "closed packet gas after compatibility constraints",
        "coimage sector, path homology residue, or support repair",
        "guards against reading H or Omega as raw scalar counts",
        "LRC packet labels unless used as residual model",
        "HYP-2618, THM-081, HYP-2990",
        dict(predicate=1, fiber=2, local=2, closedness=2, coboundary=2,
             boundary_atom=1, descent=2, named_residual=3, formal=2, lrc14=1,
             cross_domain=3),
    ),
]


def compare(a: Carrier, b: Carrier, weighted: bool = False) -> int:
    if weighted:
        sa = a.exploration_score
        sb = b.exploration_score
        if sa != sb:
            return 1 if sa > sb else -1
    wins_a = 0
    wins_b = 0
    for f in PROOF_FEATURES:
        if a.scores[f] > b.scores[f]:
            wins_a += 1
        elif b.scores[f] > a.scores[f]:
            wins_b += 1
    if wins_a != wins_b:
        return 1 if wins_a > wins_b else -1
    if a.proof_score != b.proof_score:
        return 1 if a.proof_score > b.proof_score else -1
    # Stable tie path: earlier list position wins.
    return 1 if CARRIERS.index(a) < CARRIERS.index(b) else -1


def build_edges(weighted: bool = False) -> List[List[bool]]:
    n = len(CARRIERS)
    edges = [[False] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        c = compare(CARRIERS[i], CARRIERS[j], weighted=weighted)
        if c > 0:
            edges[i][j] = True
        else:
            edges[j][i] = True
    return edges


def score_hist(edges: List[List[bool]]) -> Dict[int, int]:
    hist: Dict[int, int] = {}
    for row in edges:
        s = sum(row)
        hist[s] = hist.get(s, 0) + 1
    return dict(sorted(hist.items()))


def directed_3cycles(edges: List[List[bool]]) -> int:
    total = 0
    for i, j, k in combinations(range(len(edges)), 3):
        if edges[i][j] and edges[j][k] and edges[k][i]:
            total += 1
        if edges[i][k] and edges[k][j] and edges[j][i]:
            total += 1
    return total


def scc_sizes(edges: List[List[bool]]) -> List[int]:
    n = len(edges)
    index = 0
    stack: List[int] = []
    on_stack = [False] * n
    indices = [-1] * n
    low = [0] * n
    sizes: List[int] = []

    def strongconnect(v: int) -> None:
        nonlocal index
        indices[v] = index
        low[v] = index
        index += 1
        stack.append(v)
        on_stack[v] = True
        for w, has_edge in enumerate(edges[v]):
            if not has_edge:
                continue
            if indices[w] == -1:
                strongconnect(w)
                low[v] = min(low[v], low[w])
            elif on_stack[w]:
                low[v] = min(low[v], indices[w])
        if low[v] == indices[v]:
            size = 0
            while True:
                w = stack.pop()
                on_stack[w] = False
                size += 1
                if w == v:
                    break
            sizes.append(size)

    for v in range(n):
        if indices[v] == -1:
            strongconnect(v)
    return sizes


def hamiltonian_path_count(edges: List[List[bool]]) -> int:
    n = len(edges)
    dp: Dict[Tuple[int, int], int] = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1 << n):
        for v in range(n):
            count = dp.get((mask, v), 0)
            if not count:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if edges[v][w]:
                    dp[(mask | (1 << w), w)] = dp.get((mask | (1 << w), w), 0) + count
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


def canonical_path(edges: List[List[bool]]) -> List[str]:
    remaining = set(range(len(edges)))
    path: List[int] = []
    while remaining:
        best = max(
            remaining,
            key=lambda i: (
                sum(1 for j in remaining if i != j and edges[i][j]),
                CARRIERS[i].proof_score,
                -i,
            ),
        )
        path.append(best)
        remaining.remove(best)
    return [CARRIERS[i].name for i in path]


def edge_flips(a: List[List[bool]], b: List[List[bool]]) -> int:
    flips = 0
    n = len(a)
    for i, j in combinations(range(n), 2):
        if a[i][j] != b[i][j]:
            flips += 1
    return flips


def format_table(rows: Iterable[Tuple[str, str, str, str]]) -> str:
    lines = []
    for name, base, cochain, exits in rows:
        lines.append(f"- {name}")
        lines.append(f"  base: {base}")
        lines.append(f"  cochain: {cochain}")
        lines.append(f"  exits: {exits}")
    return "\n".join(lines)


def main() -> None:
    proof_edges = build_edges(weighted=False)
    exploration_edges = build_edges(weighted=True)

    lines: List[str] = []
    lines.append("LRC14 COCYCLE CARRIER ATLAS")
    lines.append("=" * 72)
    lines.append("")
    lines.append("Executive readout:")
    lines.append(
        "A useful LRC cocycle is a signed local obstruction to forgetting a "
        "coordinate.  It lives on a labelled packet fiber, closes when local "
        "certificates glue, and becomes harmless only as a coboundary/potential, "
        "a dual-annihilated mode, a known boundary atom, a descended family, or "
        "a named residual sector."
    )
    lines.append("")
    lines.append("Carrier dictionary:")
    lines.append(
        format_table(
            (c.name, c.base, c.cochain, c.exits)
            for c in CARRIERS
        )
    )
    lines.append("")
    lines.append("Proof-gauge Tournament Analysis:")
    lines.append(f"vertices={len(CARRIERS)}")
    lines.append(f"score_hist={score_hist(proof_edges)}")
    lines.append(f"directed_3cycles={directed_3cycles(proof_edges)}")
    lines.append(f"SCC_sizes={scc_sizes(proof_edges)}")
    lines.append(f"Hamiltonian_path_count={hamiltonian_path_count(proof_edges)}")
    lines.append("canonical_path=")
    for name in canonical_path(proof_edges):
        lines.append(f"  > {name}")
    lines.append("")
    lines.append("Exploration-gauge comparison:")
    lines.append(f"score_hist={score_hist(exploration_edges)}")
    lines.append(f"directed_3cycles={directed_3cycles(exploration_edges)}")
    lines.append(f"SCC_sizes={scc_sizes(exploration_edges)}")
    lines.append(f"Hamiltonian_path_count={hamiltonian_path_count(exploration_edges)}")
    lines.append(f"edge_flips_vs_proof_gauge={edge_flips(proof_edges, exploration_edges)}")
    lines.append("")
    lines.append("Top cocycle proof obligations:")
    obligations = [
        (
            "Packet cocycle theorem",
            "Define omega_Q for every quotient Q on P(S)-fibers and prove each "
            "omega_Q is exact, dual-annihilated, descended, AP/GW, or residualized.",
        ),
        (
            "Endpoint-credit potential backend",
            "Turn HYP-2970's no-positive-winding condition into exact interval "
            "certificates for packet families.",
        ),
        (
            "Haar zeta packet grid",
            "Compute zeta/product-rule classes on HYP-2963 packets, not only "
            "abstract 2 x 2 tables.",
        ),
        (
            "Carry-owner no-leak theorem",
            "Prove local/global carry cocycles either reconstruct from paired "
            "owner derivatives or pay strict loneliness tax.",
        ),
        (
            "F7 definition",
            "Define F7 as a non-exact, non-annihilated cocycle class after all "
            "known teeth fail; specify its coefficient ring and state-lift map.",
        ),
    ]
    for name, body in obligations:
        lines.append(f"- {name}: {body}")
    lines.append("")
    lines.append("Cocycle packet theorem skeleton:")
    lines.append(
        "Let P(S) retain qdiv, exact M/Farey data, Haar/Baire topology, endpoint "
        "owners, exact-period labels, Fejer certificate data, C27/K33/state "
        "labels, and smoothing route.  For any quotient Q:P(S)->Y, construct a "
        "fiber cochain omega_Q measuring the coordinate Q forgets.  Q is "
        "theorem-safe only if the LRC predicate is constant on Y-fibers or "
        "omega_Q is reconstructed, exact as a coboundary/potential, annihilated "
        "by a dual certificate, descended to a smaller packet family, identified "
        "as AP/GW boundary equality, or emitted as a named F7/THM-572 residual."
    )
    lines.append("")
    lines.append("Most important change of viewpoint:")
    lines.append(
        "Tournaments are not only built on runners.  For cocycle work the better "
        "vertices are fibers, proof carriers, endpoint arrows, Haar squares, "
        "carry derivatives, exact-period modes, transfer states, and residual "
        "obligations.  Edges compare which carrier preserves the cocycle needed "
        "for the next gluing step."
    )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
