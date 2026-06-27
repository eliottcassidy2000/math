#!/usr/bin/env python3
"""S166: zipper theorem pattern atlas for the current LRC14 proof stack.

This is a synthesis scout, not an LRC14 proof.

The prompt asks to pull recent agents together through discrepancy theory, the
two-dimensional Haar product rule, tournament tilings, and "more zipper
theorems."  The shared structure across the recent repo lanes is:

    local certificate side A
    local certificate side B
    a labelled interface where they meet
    equality stops
    a named residual when they do not meet

An admissible quotient may forget a coordinate only when the other side of the
zipper reconstructs it, orthogonality/boundary atoms annihilate it, or the
forgotten coordinate is emitted as a residual proof obligation.

Tournament Analysis:
  vertices: zipper proof patterns, not runners.
  pairwise observable: which pattern retains more theorem-bearing payload:
      LRC predicate, local-to-global gluing, boundary stops, residual naming,
      formal checkability, cross-domain transfer, anti-scalar guardrails, and
      computable next tests.
  switch/gauge: dimensionwise majority over that retention vector, then
      weighted score, then a fixed Hamiltonian tie path.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations


RETENTION_KEYS = (
    "predicate",
    "local_gluing",
    "boundary_stop",
    "residual_name",
    "formal_check",
    "cross_transfer",
    "anti_scalar",
    "computable_next",
)

WEIGHTS = {
    "predicate": 6,
    "local_gluing": 5,
    "boundary_stop": 5,
    "residual_name": 6,
    "formal_check": 4,
    "cross_transfer": 3,
    "anti_scalar": 5,
    "computable_next": 4,
}


@dataclass(frozen=True)
class ZipperPattern:
    name: str
    domain: str
    left_teeth: str
    right_teeth: str
    boundary_stop: str
    residual_atom: str
    quotient_warning: str
    proof_arrow: str
    related: tuple[str, ...]
    retention: dict[str, int]

    def weighted_score(self) -> int:
        return sum(WEIGHTS[k] * self.retention[k] for k in RETENTION_KEYS)


PATTERNS = [
    ZipperPattern(
        name="haar_fourier_product",
        domain="discrepancy / Haar product / tournament tiling",
        left_teeth="row-column or endpoint-wall fibers",
        right_teeth="mixed Haar products, fixed-margin switches, Fejer modes",
        boundary_stop="same-tile AP/GW boundary atom",
        residual_atom="named mixed-product F7 harmonic coefficient",
        quotient_warning="row and column margins forget the checkerboard sign",
        proof_arrow="replace raw component count K by independent labelled mixed switches",
        related=("HYP-2992", "HYP-2989", "HYP-2595", "HYP-2594"),
        retention=dict(
            predicate=3,
            local_gluing=3,
            boundary_stop=2,
            residual_name=3,
            formal_check=2,
            cross_transfer=3,
            anti_scalar=3,
            computable_next=3,
        ),
    ),
    ZipperPattern(
        name="fejer_interval_packet",
        domain="Fourier-Toeplitz dual / interval certificates",
        left_teeth="exact safe component, rational center, packet fiber",
        right_teeth="divisor-curried Fejer atom bank and negative interval bound",
        boundary_stop="AP/GW PSD-blind equality atom",
        residual_atom="family-compression or formal-interval backend debt",
        quotient_warning="floating Fejer signs are useless after packet labels are dropped",
        proof_arrow="compress HYP-2963 rows into reusable packet-family certificates",
        related=("HYP-2981", "HYP-2974", "HYP-2963"),
        retention=dict(
            predicate=3,
            local_gluing=2,
            boundary_stop=3,
            residual_name=2,
            formal_check=3,
            cross_transfer=2,
            anti_scalar=3,
            computable_next=3,
        ),
    ),
    ZipperPattern(
        name="tope_cocircuit_wall",
        domain="oriented endpoint arrangement",
        left_teeth="open all-safe topes",
        right_teeth="all-safe boundary cocircuits with owner labels",
        boundary_stop="zero-dimensional AP/GW owner-sum wall",
        residual_atom="no-tope/no-cocircuit forbidden wall packet",
        quotient_warning="safe Haar mass blurs positive open witnesses with closed equality",
        proof_arrow="classify forbidden wall packets by owner parity, wall slide, K33, or F7",
        related=("HYP-2986", "HYP-2951", "HYP-2948"),
        retention=dict(
            predicate=3,
            local_gluing=2,
            boundary_stop=3,
            residual_name=3,
            formal_check=3,
            cross_transfer=2,
            anti_scalar=3,
            computable_next=2,
        ),
    ),
    ZipperPattern(
        name="exposure_poset_kernel",
        domain="proof-channel exposure audit",
        left_teeth="row witnesses and packet labels",
        right_teeth="q-witness, open-Haar, Fejer, K33/C27, pressure channels",
        boundary_stop="AP/GW taut boundary channel",
        residual_atom="UNEXPOSED_SOURCE_KERNEL",
        quotient_warning="a missing numerical witness is not a counterexample until all channels fail",
        proof_arrow="prove the no-hidden-exposure-kernel lemma familywise",
        related=("HYP-2988", "HYP-2987", "HYP-2981"),
        retention=dict(
            predicate=3,
            local_gluing=2,
            boundary_stop=2,
            residual_name=3,
            formal_check=3,
            cross_transfer=2,
            anti_scalar=3,
            computable_next=3,
        ),
    ),
    ZipperPattern(
        name="ramanujan_exact_period",
        domain="Ramanujan sums / primitive phases / divisor guardrails",
        left_teeth="primitive q-period phase packets",
        right_teeth="endpoint-owner sums, Fejer/Ramanujan projections, unit shells",
        boundary_stop="q=14 primitive trace of AP/GW owner cancellation",
        residual_atom="prime-power or repeated-prime packet side channel",
        quotient_warning="squarefree or gcd scalars erase q=25,27,36,63,84,98,168,280,4312",
        proof_arrow="tensor primitive-period projectors with Haar owner grids",
        related=("HYP-2979", "HYP-2978", "HYP-2982"),
        retention=dict(
            predicate=3,
            local_gluing=2,
            boundary_stop=2,
            residual_name=3,
            formal_check=2,
            cross_transfer=3,
            anti_scalar=3,
            computable_next=2,
        ),
    ),
    ZipperPattern(
        name="smoothing_kaczynski_policy",
        domain="analytic smoothing / Kaczynski-Abel boundary",
        left_teeth="kernel or smoothing policy",
        right_teeth="endpoint owner, exact period, far-approach boundary clock",
        boundary_stop="boundary-defect atom emitted by a failed homotopy",
        residual_atom="state-lift obligation from failed admissible policies",
        quotient_warning="large-sieve and Selberg weights are middle packets, not final certificates",
        proof_arrow="turn every smoothing choice into a declared packet-family policy",
        related=("HYP-2985", "HYP-2984", "HYP-2983", "HYP-2982"),
        retention=dict(
            predicate=2,
            local_gluing=3,
            boundary_stop=2,
            residual_name=3,
            formal_check=2,
            cross_transfer=3,
            anti_scalar=3,
            computable_next=2,
        ),
    ),
    ZipperPattern(
        name="fixed_margin_johnson",
        domain="fixed-margin reductions / Johnson harmonic sectors",
        left_teeth="fixed margin or packet fibers",
        right_teeth="low-row core, count sector, Johnson harmonic sector",
        boundary_stop="finite core comparison",
        residual_atom="named F7 Johnson-harmonic packet",
        quotient_warning="count sectors are safe only after packet margins are fixed",
        proof_arrow="define F7 as a preserved harmonic residual, not an unknown bucket",
        related=("HYP-2987", "HYP-2963", "HYP-2956", "arXiv:2606.22636"),
        retention=dict(
            predicate=3,
            local_gluing=3,
            boundary_stop=2,
            residual_name=3,
            formal_check=2,
            cross_transfer=3,
            anti_scalar=2,
            computable_next=2,
        ),
    ),
    ZipperPattern(
        name="apex_sheaf_gluing",
        domain="apex-lift sheaf / Vitali handoff",
        left_teeth="local cheap-pair sections over mod-7 or mod-27 charts",
        right_teeth="endpoint anchors, small-pair pinches, measure-positive covers",
        boundary_stop="AP/V* tie-wall section",
        residual_atom="failed gluing that must create positive measure",
        quotient_warning="round-class identity forgets the local section labels",
        proof_arrow="prove local sections glue or obstruction positivity fires",
        related=("HYP-2101", "HYP-2104", "THM-397", "THM-398"),
        retention=dict(
            predicate=2,
            local_gluing=3,
            boundary_stop=3,
            residual_name=2,
            formal_check=2,
            cross_transfer=2,
            anti_scalar=3,
            computable_next=2,
        ),
    ),
    ZipperPattern(
        name="convolution_irreducibility_lift",
        domain="irreducible polynomials / coefficient tilings",
        left_teeth="visible coefficient boundary totals",
        right_teeth="hidden integral factor grid and value-factor token split",
        boundary_stop="fixed-divisor or Cohn prime place-value certificate",
        residual_atom="no nontrivial convolution lift",
        quotient_warning="coefficient profile forgets the hidden support/factor grid",
        proof_arrow="port no-integral-lift certificates back to LRC packet ledgers and code supports",
        related=("HYP-2452", "HYP-2450", "HYP-2430"),
        retention=dict(
            predicate=2,
            local_gluing=3,
            boundary_stop=2,
            residual_name=3,
            formal_check=3,
            cross_transfer=3,
            anti_scalar=3,
            computable_next=2,
        ),
    ),
    ZipperPattern(
        name="unit_distance_cyclotomic_norm",
        domain="unit-distance carriers / cyclotomic norm shells",
        left_teeth="finite additive patches or distance graph constraints",
        right_teeth="cyclotomic/norm shell and lattice-unit compatibility",
        boundary_stop="known lattice/equilateral/cyclotomic equality carrier",
        residual_atom="non-product configuration with retained norm labels",
        quotient_warning="Euclidean distance alone forgets algebraic unit and shell data",
        proof_arrow="use norm-shell labels as a quotient guardrail for geometry transfers",
        related=("unit-distance", "cyclotomic", "HYP-2978", "HYP-2450"),
        retention=dict(
            predicate=2,
            local_gluing=2,
            boundary_stop=2,
            residual_name=2,
            formal_check=1,
            cross_transfer=3,
            anti_scalar=2,
            computable_next=1,
        ),
    ),
]


TIE_PATH = [p.name for p in PATTERNS]


def majority_margin(a: ZipperPattern, b: ZipperPattern) -> int:
    return sum(
        (a.retention[k] > b.retention[k]) - (a.retention[k] < b.retention[k])
        for k in RETENTION_KEYS
    )


def orient(a: ZipperPattern, b: ZipperPattern) -> tuple[str, str]:
    margin = majority_margin(a, b)
    if margin > 0:
        return a.name, b.name
    if margin < 0:
        return b.name, a.name

    sa, sb = a.weighted_score(), b.weighted_score()
    if sa > sb:
        return a.name, b.name
    if sb > sa:
        return b.name, a.name

    return (
        (a.name, b.name)
        if TIE_PATH.index(a.name) < TIE_PATH.index(b.name)
        else (b.name, a.name)
    )


def tournament_edges() -> dict[tuple[str, str], str]:
    edges: dict[tuple[str, str], str] = {}
    for a, b in combinations(PATTERNS, 2):
        winner, loser = orient(a, b)
        edges[(winner, loser)] = f"margin={majority_margin(a, b)}"
    return edges


def score_hist(edges: dict[tuple[str, str], str]) -> tuple[dict[int, int], dict[str, int]]:
    scores = {p.name: 0 for p in PATTERNS}
    for winner, _ in edges:
        scores[winner] += 1
    hist: dict[int, int] = {}
    for score in scores.values():
        hist[score] = hist.get(score, 0) + 1
    return dict(sorted(hist.items())), scores


def directed_3cycles(edges: dict[tuple[str, str], str]) -> int:
    names = [p.name for p in PATTERNS]
    edge_set = {(a, b) for a, b in edges}
    total = 0
    for a, b, c in combinations(names, 3):
        if (
            ((a, b) in edge_set and (b, c) in edge_set and (c, a) in edge_set)
            or ((a, c) in edge_set and (c, b) in edge_set and (b, a) in edge_set)
        ):
            total += 1
    return total


def scc_sizes(edges: dict[tuple[str, str], str]) -> list[int]:
    names = [p.name for p in PATTERNS]
    out = {n: set() for n in names}
    rev = {n: set() for n in names}
    for a, b in edges:
        out[a].add(b)
        rev[b].add(a)

    seen: set[str] = set()
    order: list[str] = []

    def dfs(v: str) -> None:
        seen.add(v)
        for w in out[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for n in names:
        if n not in seen:
            dfs(n)

    seen.clear()
    sizes: list[int] = []

    def rdfs(v: str, bag: list[str]) -> None:
        seen.add(v)
        bag.append(v)
        for w in rev[v]:
            if w not in seen:
                rdfs(w, bag)

    for n in reversed(order):
        if n not in seen:
            bag: list[str] = []
            rdfs(n, bag)
            sizes.append(len(bag))
    return sorted(sizes, reverse=True)


def hamiltonian_path_count(edges: dict[tuple[str, str], str]) -> int:
    names = [p.name for p in PATTERNS]
    idx = {name: i for i, name in enumerate(names)}
    out_mask = [0] * len(names)
    for a, b in edges:
        out_mask[idx[a]] |= 1 << idx[b]

    full = (1 << len(names)) - 1
    dp: dict[tuple[int, int], int] = {}
    for i in range(len(names)):
        dp[(1 << i, i)] = 1
    for mask in range(1 << len(names)):
        for last in range(len(names)):
            count = dp.get((mask, last), 0)
            if not count:
                continue
            candidates = out_mask[last] & ~mask
            while candidates:
                bit = candidates & -candidates
                nxt = bit.bit_length() - 1
                dp[(mask | bit, nxt)] = dp.get((mask | bit, nxt), 0) + count
                candidates ^= bit
    return sum(dp.get((full, i), 0) for i in range(len(names)))


def theorem_schema() -> str:
    return """Zipper theorem schema:
  For every primitive packet P in a classified family F:
    1. compute the left teeth L(P) and right teeth R(P) with labels attached;
    2. if L(P) and R(P) meet compatibly, close the certificate;
    3. if P lands on a declared stop, record the equality/boundary atom;
    4. otherwise emit a residual whose labels are strong enough for the next theorem.

  A quotient may forget coordinate z only if z is reconstructed by the opposite
  tooth, annihilated by orthogonality or a boundary stop, or named in the residual.
"""


def main() -> None:
    edges = tournament_edges()
    hist, scores = score_hist(edges)

    print("S166 ZIPPER THEOREM PATTERN ATLAS")
    print("status: synthesis scout / proof-interface map, not an LRC14 proof")
    print("namespace repair: HYP-2988=exposure, HYP-2989=minimal Haar square, HYP-2990=abstract zipper atlas, HYP-2991=local Haar zipper cocycle, HYP-2993=concrete LRC zipper pattern atlas, HYP-2992=Haar tile atlas")
    print()
    print(theorem_schema())

    print("PATTERNS")
    for i, p in enumerate(PATTERNS, 1):
        print(f"{i:02d}. {p.name}")
        print(f"    domain: {p.domain}")
        print(f"    left: {p.left_teeth}")
        print(f"    right: {p.right_teeth}")
        print(f"    stop: {p.boundary_stop}")
        print(f"    residual: {p.residual_atom}")
        print(f"    quotient warning: {p.quotient_warning}")
        print(f"    proof arrow: {p.proof_arrow}")
        print(f"    related: {', '.join(p.related)}")
        print(f"    retention_score: {p.weighted_score()} vector={tuple(p.retention[k] for k in RETENTION_KEYS)}")
    print()

    print("TOURNAMENT ANALYSIS")
    print(f"vertices={len(PATTERNS)}")
    print("pairwise_observable=retained theorem payload across predicate, local gluing, stops, residuals, formal checkability, transfer, anti-scalar guardrail, computable next test")
    print("switch_gauge=dimensionwise majority; weighted score; fixed Hamiltonian tie path")
    print("tie_path=" + " > ".join(TIE_PATH))
    print(f"score_hist={hist}")
    print(f"directed_3cycles={directed_3cycles(edges)}")
    print(f"SCC_sizes={scc_sizes(edges)}")
    print(f"Hamiltonian_path_count={hamiltonian_path_count(edges)}")
    print("scores=" + ", ".join(f"{name}:{scores[name]}" for name in sorted(scores)))
    print()

    print("HIGH-LEVERAGE NEXT MOVES")
    moves = [
        "Haar-Fejer compression: group HYP-2963 packet rows by mixed Haar switch signatures before interval certificate generation.",
        "Ramanujan-Haar tensor zipper: attach primitive-period c_q labels to each endpoint-wall/Haar rectangle cell.",
        "No-hidden-kernel theorem split: prove exposure_poset_kernel has no residual after HYP-2992 typed coefficients and HYP-2981 interval certificates are attached.",
        "F7 definition: use fixed_margin_johnson to define a residual by preserved harmonic sector rather than by failure of known tests.",
        "Convolution-lift transfer: treat LRC blocker ledgers like coefficient boundary totals and search for hidden support-lift infeasibility certificates.",
        "Apex/Vitali transfer: use local section gluing as the measure-theoretic analogue of the Haar same-tile boundary stop.",
    ]
    for i, move in enumerate(moves, 1):
        print(f"{i}. {move}")

    print()
    print("CORE HYPOTHESIS")
    print(
        "The LRC14 proof stack is converging on a labelled packet zipper: "
        "positive-open Haar/topes and negative Fejer/Toeplitz duals are the two "
        "sides; AP/GW are stops; K33/C27/petal/covering/F7 are named residual "
        "routes.  The missing proof should show that no primitive packet can "
        "remain after both sides close with quotient-safe labels attached."
    )


if __name__ == "__main__":
    main()
