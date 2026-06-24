#!/usr/bin/env python3
"""S165: Abstract zipper theorem atlas.

This is an exploratory proof-technology artifact, not a theorem proof.  The
prompt asks for more "zipper theorems" across past topics.  The common pattern
in the recent LRC14 work, unit-distance spines, Farey-product/K33 packets,
octahedral currents, root packets, boundary moments, OCF coimages, and
good-cut/SCC support is not an average or scalar analogy.  It is a controlled
kernel rule:

    a quotient may forget a coordinate only when the target predicate is
    constant on fibers, the coordinate is reconstructible from retained
    labels, a dual certificate annihilates the forgotten direction, or the
    direction is routed to a named residual sector.

Each "zipper candidate" below records two toothed carriers, a slider that
identifies their fibers, the predicate preserved by the join, and the labels
that cannot be forgotten.  Tournament Analysis compares proof carriers rather
than runners.  The pairwise observable is majority comparison of retention
features:

    predicate, fiber labels, two-sided transfer, declared kernel, finite
    checkability, family compression, formalizability, anti-scalar guardrail,
    and cross-domain signal.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations


VECTOR_KEYS = (
    "predicate",
    "fiber_label",
    "two_sided_transfer",
    "kernel_declared",
    "finite_check",
    "family_compression",
    "formalizable",
    "anti_scalar_guard",
    "cross_domain_signal",
)


@dataclass(frozen=True)
class ZipperCandidate:
    name: str
    prior: str
    domains: tuple[str, ...]
    left_teeth: str
    right_teeth: str
    slider: str
    preserved_predicate: str
    forbidden_forgetting: tuple[str, ...]
    open_arrow: str
    theorem_template: str
    vector: tuple[int, ...]


CANDIDATES = [
    ZipperCandidate(
        name="lrc14_certificate_handoff_zipper",
        prior="HYP-2987/T1071",
        domains=("LRC14", "Fejer", "Ramanujan", "state-lift"),
        left_teeth="labelled source packets with qdiv/Farey/Haar/endpoints",
        right_teeth="Fejer, Ramanujan, endpoint, twist, moment, and state-lift certificates",
        slider="packet fiber P(S)",
        preserved_predicate="strict witness, AP/GW equality, or forbidden THM-572 state lift",
        forbidden_forgetting=(
            "exact M=p/q and e=14p-q",
            "strict-open versus boundary-only status",
            "endpoint owners and packet route",
            "formal certificate payload",
        ),
        open_arrow="prove O1-O6: source-kernel, interval backend, family compression, smoothing, state lift, F7",
        theorem_template="labelled packet zipper theorem",
        vector=(5, 5, 5, 4, 4, 4, 4, 5, 5),
    ),
    ZipperCandidate(
        name="kernel_tope_smoothing_zipper",
        prior="HYP-2984/HYP-2985/HYP-2986",
        domains=("LRC14", "kernel homotopy", "oriented topes", "smoothing"),
        left_teeth="open safe topes and boundary cocircuits",
        right_teeth="smoothing policies and Fejer/Ramanujan/Kaczynski exits",
        slider="regular-open component or named boundary-defect atom",
        preserved_predicate="open-stable packet certificate or AP/GW boundary atom",
        forbidden_forgetting=(
            "support radius",
            "boundary owner pair sums",
            "smoothing clock",
            "approach class",
        ),
        open_arrow="lift named-row homotopy/topes into HYP-2963 packet families",
        theorem_template="admissible kernel-homotopy zipper",
        vector=(5, 5, 4, 5, 4, 3, 4, 5, 4),
    ),
    ZipperCandidate(
        name="octahedral_hodge_current_zipper",
        prior="HYP-2887",
        domains=("octahedral graph", "Hodge", "LRC14 currents"),
        left_teeth="vertex divergence on the L(K4) octahedral support",
        right_teeth="face curls plus Abel/Cauchy tail defects",
        slider="octahedral current class after low-height wall deletion",
        preserved_predicate="realizability of a repeated-packet lift or named curl obstruction",
        forbidden_forgetting=(
            "divergence source labels",
            "face-curl owner labels",
            "tail estimate used to kill residual flow",
        ),
        open_arrow="prove the octahedral Hodge estimate after wall deletion",
        theorem_template="divergence-curl zipper theorem",
        vector=(4, 4, 5, 5, 3, 3, 3, 4, 5),
    ),
    ZipperCandidate(
        name="c27_unital_pair_completion_zipper",
        prior="HYP-2942",
        domains=("C27 shell", "unital design", "AP/GW", "K33"),
        left_teeth="AP/GW calibrated pair-completion slots",
        right_teeth="branch-local K33, H12, D9, and petal completion packets",
        slider="unit pair slot with retained branch label",
        preserved_predicate="unique branch-local four-point completion or state-lift debt",
        forbidden_forgetting=(
            "q=3 unital branch",
            "AP/GW calibration",
            "near/K33 versus petal route",
        ),
        open_arrow="turn positive branch-local lift into a global state-lift discharge",
        theorem_template="pair-completion zipper theorem",
        vector=(4, 5, 4, 4, 5, 3, 3, 5, 4),
    ),
    ZipperCandidate(
        name="farey_product_k33_zipper",
        prior="HYP-2932/HYP-2934/HYP-2945",
        domains=("Farey", "complete bipartite graphs", "perfect products", "K33"),
        left_teeth="exact Farey branch e=14p-q and denominator shell",
        right_teeth="K_{p,q} incidence blow-up and product/factor ledger",
        slider="unit-excess route index p after exact M is retained",
        preserved_predicate="p=2 C27 petal/two-block branch or p>=3 K33 incidence packet",
        forbidden_forgetting=(
            "denominator q",
            "summand shell C=27",
            "factor-fiber address",
            "minor label",
        ),
        open_arrow="prove every remaining q=14 non-AP/GW atom enters C27 or K33 packet lane",
        theorem_template="incidence-product zipper theorem",
        vector=(4, 5, 3, 3, 4, 4, 3, 5, 5),
    ),
    ZipperCandidate(
        name="boundary_moment_multichart_zipper",
        prior="HYP-2969",
        domains=("LRC14", "boundary moments", "fixed-margin packets", "gK8"),
        left_teeth="adaptive exact-period packet boundary charts",
        right_teeth="missed-depth sector vector and L_y/gK8 moment readouts",
        slider="labelled packet over multiple denominator charts",
        preserved_predicate="positive Haar/moment image, K33 state-lift debt, or named harmonic sector",
        forbidden_forgetting=(
            "chart identity",
            "packet label",
            "missed-depth sector",
            "K33/TournamentStateLift debt",
        ),
        open_arrow="replace one-chart all-covered tests by a multi-chart feasible-region theorem",
        theorem_template="boundary-moment bridge zipper",
        vector=(5, 5, 4, 4, 3, 3, 3, 5, 4),
    ),
    ZipperCandidate(
        name="shell1_root_packet_mouth_zipper",
        prior="HYP-1815/HYP-2664",
        domains=("LRC14", "root packets", "shell-1", "comb pruning"),
        left_teeth="shell-1 carry gate on the tower {1,2,4,8}",
        right_teeth="root-packet incidence rank and mouth-owner ledger",
        slider="shell-full addressed packet with holes/tails/carry retained",
        preserved_predicate="finite AP-tail residue burden after carry gate",
        forbidden_forgetting=(
            "which shell-1 bit is deleted",
            "root-packet incidence rank",
            "old-mouth owner",
        ),
        open_arrow="prove shell-1 deletion theorem and then enumerate only shell-full packets",
        theorem_template="carry-gate/root-rank zipper theorem",
        vector=(4, 5, 4, 4, 4, 3, 2, 5, 4),
    ),
    ZipperCandidate(
        name="unit_distance_spine_ear_zipper",
        prior="HYP-2620/THM-408",
        domains=("unit distance", "Hamiltonian spines", "endpoint ears"),
        left_teeth="traceable unit graph and endpoint masks",
        right_teeth="endpoint-compatible deletion ears and geometric cocyclic extensions",
        slider="spine endpoint section over vertex deletion",
        preserved_predicate="unit Hamiltonian spine survives extension unless geometry forbids the attachment",
        forbidden_forgetting=(
            "endpoint option set",
            "deleted-ear identity",
            "unit-cocyclic realizability",
        ),
        open_arrow="classify degree-4/5 unit-cocyclic attachments over exact n=21 cores",
        theorem_template="endpoint-ear zipper theorem",
        vector=(4, 4, 5, 4, 5, 3, 3, 4, 5),
    ),
    ZipperCandidate(
        name="ocf_activity_coimage_zipper",
        prior="HYP-2618",
        domains=("OCF", "Condorcet", "hard-core activity", "LRC coimages"),
        left_teeth="hard-core activity-2 partition of Omega(T)",
        right_teeth="compatible packet address and signed partition evaluation",
        slider="coimage carrier with retained compatibility labels",
        preserved_predicate="forbidden H-value/coimage signal without pretending it is raw noise stability",
        forbidden_forgetting=(
            "compatibility relation",
            "activity parameter",
            "signed versus unsigned partition function",
        ),
        open_arrow="define the exact LRC packet coimage evaluation, not just its mass shadow",
        theorem_template="activity/coimage zipper theorem",
        vector=(3, 4, 4, 3, 4, 3, 2, 4, 5),
    ),
    ZipperCandidate(
        name="good_cut_scc_support_zipper",
        prior="THM-354/INV-237",
        domains=("staircase tilings", "tournaments", "SCC support", "good cuts"),
        left_teeth="base-path interval good-cut count",
        right_teeth="strong-component condensation boundaries",
        slider="support residue g_P(T)=n-#SCC(T)",
        preserved_predicate="strong connectivity/top good-cut bucket across the tiling-to-tournament map",
        forbidden_forgetting=(
            "chosen Hamiltonian base path",
            "component-boundary order",
            "range parity of quotient transport",
        ),
        open_arrow="use SCC support residue as a quotient-transport coordinate in LRC packet atlases",
        theorem_template="support-residue zipper theorem",
        vector=(4, 4, 4, 4, 5, 3, 4, 4, 5),
    ),
    ZipperCandidate(
        name="raw_scalar_shadow",
        prior="negative control",
        domains=("all domains",),
        left_teeth="one number",
        right_teeth="another number",
        slider="unlabelled equality or monotone comparison",
        preserved_predicate="usually none unless a separate theorem supplies the kernel",
        forbidden_forgetting=(
            "route",
            "fiber",
            "endpoint owner",
            "topology",
            "certificate payload",
        ),
        open_arrow="use only as a diagnostic after the zipper labels have been retained",
        theorem_template="anti-template",
        vector=(1, 1, 1, 1, 1, 1, 1, 1, 1),
    ),
]

TIE_PATH = [candidate.name for candidate in CANDIDATES]

THEOREM_TEMPLATES = [
    (
        "Labelled packet zipper lemma",
        "If two carriers L and R map to the same labelled packet base P, the "
        "target predicate is constant on P-fibers, and every kernel coordinate "
        "is reconstructed, annihilated, or routed to a named residual sector, "
        "then local L- and R-certificates glue to a certificate over the fiber "
        "product L x_P R.",
    ),
    (
        "No-free-slider lemma",
        "A proposed zipper fails exactly where the slider forgets a load-bearing "
        "label.  A failure certificate should name one of four defects: predicate "
        "mixing on a fiber, unreconstructed coordinate, unannihilated kernel, or "
        "unnamed residual sector.",
    ),
    (
        "Alternating gauntlet theorem",
        "Run left teeth and right teeth alternately.  If every failed tooth emits "
        "a smaller named packet or a boundary-moment/state-lift debt, and the "
        "debt order is well founded, the gauntlet terminates in a certificate, "
        "a known equality atom, or a forbidden lift.",
    ),
    (
        "Harmonic residual theorem",
        "The last bucket cannot be an anonymous exception.  It must be a named "
        "representation, homology, cocircuit, curl, Johnson-harmonic, or "
        "state-lift sector with an explicit predicate.",
    ),
]


def candidate_index() -> dict[str, int]:
    return {candidate.name: i for i, candidate in enumerate(CANDIDATES)}


def beats(a: ZipperCandidate, b: ZipperCandidate) -> bool:
    wins_a = sum(1 for x, y in zip(a.vector, b.vector) if x > y)
    wins_b = sum(1 for x, y in zip(a.vector, b.vector) if x < y)
    if wins_a != wins_b:
        return wins_a > wins_b
    return candidate_index()[a.name] < candidate_index()[b.name]


def tournament_edges() -> dict[str, set[str]]:
    edges = {candidate.name: set() for candidate in CANDIDATES}
    for a, b in combinations(CANDIDATES, 2):
        if beats(a, b):
            edges[a.name].add(b.name)
        else:
            edges[b.name].add(a.name)
    return edges


def score_hist(edges: dict[str, set[str]]) -> dict[int, int]:
    hist: dict[int, int] = {}
    for outs in edges.values():
        hist[len(outs)] = hist.get(len(outs), 0) + 1
    return dict(sorted(hist.items()))


def directed_3cycles(edges: dict[str, set[str]]) -> list[tuple[str, str, str]]:
    cycles: list[tuple[str, str, str]] = []
    names = [candidate.name for candidate in CANDIDATES]
    for a, b, c in combinations(names, 3):
        if b in edges[a] and c in edges[b] and a in edges[c]:
            cycles.append((a, b, c))
        elif a in edges[b] and c in edges[a] and b in edges[c]:
            cycles.append((b, a, c))
    return cycles


def sccs(edges: dict[str, set[str]]) -> list[list[str]]:
    names = [candidate.name for candidate in CANDIDATES]
    index = 0
    stack: list[str] = []
    on_stack: set[str] = set()
    idx: dict[str, int] = {}
    low: dict[str, int] = {}
    comps: list[list[str]] = []

    def visit(v: str) -> None:
        nonlocal index
        idx[v] = index
        low[v] = index
        index += 1
        stack.append(v)
        on_stack.add(v)
        for w in sorted(edges[v], key=candidate_index().get):
            if w not in idx:
                visit(w)
                low[v] = min(low[v], low[w])
            elif w in on_stack:
                low[v] = min(low[v], idx[w])
        if low[v] == idx[v]:
            comp: list[str] = []
            while True:
                w = stack.pop()
                on_stack.remove(w)
                comp.append(w)
                if w == v:
                    break
            comps.append(sorted(comp, key=candidate_index().get))

    for name in names:
        if name not in idx:
            visit(name)
    return comps


def hamiltonian_path_count(edges: dict[str, set[str]]) -> int:
    names = [candidate.name for candidate in CANDIDATES]
    n = len(names)

    @lru_cache(None)
    def dp(mask: int, last: int) -> int:
        if mask == (1 << last):
            return 1
        prev_mask = mask ^ (1 << last)
        total = 0
        for prev in range(n):
            if prev_mask & (1 << prev) and names[last] in edges[names[prev]]:
                total += dp(prev_mask, prev)
        return total

    full = (1 << n) - 1
    return sum(dp(full, last) for last in range(n))


def sorted_scores(edges: dict[str, set[str]]) -> list[tuple[str, int]]:
    return sorted(
        ((candidate.name, len(edges[candidate.name])) for candidate in CANDIDATES),
        key=lambda item: (-item[1], candidate_index()[item[0]]),
    )


def main() -> None:
    edges = tournament_edges()
    cycles = directed_3cycles(edges)
    comps = sccs(edges)
    scores = sorted_scores(edges)

    print("S165 abstract zipper theorem atlas")
    print("=" * 42)
    print()
    print("retention vector keys:")
    for i, key in enumerate(VECTOR_KEYS, 1):
        print(f"  {i}. {key}")
    print()

    print("candidate zipper packets:")
    for candidate in CANDIDATES:
        print(f"- {candidate.name} [{candidate.prior}]")
        print(f"  domains: {', '.join(candidate.domains)}")
        print(f"  left teeth: {candidate.left_teeth}")
        print(f"  right teeth: {candidate.right_teeth}")
        print(f"  slider: {candidate.slider}")
        print(f"  preserves: {candidate.preserved_predicate}")
        print(f"  forbidden forgetting: {', '.join(candidate.forbidden_forgetting)}")
        print(f"  open arrow: {candidate.open_arrow}")
        print(f"  template: {candidate.theorem_template}")
        print(f"  vector: {candidate.vector}")
    print()

    print("tournament fingerprint:")
    print(f"  score_hist={score_hist(edges)}")
    print(f"  directed_3cycles={len(cycles)}")
    print(f"  SCC_sizes={[len(comp) for comp in comps]}")
    print(f"  Hamiltonian_path_count={hamiltonian_path_count(edges)}")
    print("  top score order:")
    for name, score in scores:
        print(f"    {score:2d}  {name}")
    if cycles:
        print("  sample directed 3-cycles:")
        for cycle in cycles[:8]:
            print("    " + " -> ".join(cycle) + f" -> {cycle[0]}")
    print()

    print("theorem templates:")
    for title, body in THEOREM_TEMPLATES:
        print(f"- {title}: {body}")
    print()

    print("LRC14 transfer readout:")
    print(
        "  The strongest next zipper is not another scalar invariant.  It is a "
        "labelled-packet theorem saying that the Fejer interval, Ramanujan "
        "exact-period, endpoint bridge, tope/cocircuit, boundary-moment, and "
        "state-lift carriers are all projections of one packet object."
    )
    print(
        "  The finite obstruction target is therefore a controlled residual: "
        "F7 must name the harmonic/state-lift sector left after every zipper "
        "tooth has either certified a strict interval or emitted AP/GW equality."
    )
    print(
        "  The cross-domain warning is stable: when a quotient forgets endpoint "
        "ears, shell bits, C27 branch labels, curl/divergence owners, or exact "
        "periods, the theorem does not become false by mystery.  It becomes "
        "untyped."
    )


if __name__ == "__main__":
    main()
