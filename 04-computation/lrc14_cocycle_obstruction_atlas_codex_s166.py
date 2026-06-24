#!/usr/bin/env python3
"""S166: cocycle obstruction atlas for the current LRC14 proof stack.

This is a proof-interface computation, not a proof of LRC14.

HYP-2991 identified the local fixed-margin Haar coordinate

    zeta(T) = T00 - T01 - T10 + T11

as a genuine cocycle-like obstruction: margins alone forget the switch
direction.  HYP-2994 lifts that viewpoint from one square to the whole proof
stack.  The useful abstraction is a small cochain ledger:

    C0: packet labels, owner potentials, sections, exact-period residues
    C1: handoff arrows, endpoint transfers, smoothing gauges, source pullbacks
    C2: Haar switches, tope curls, moment curls, color squares, state-lift faces

A proof quotient is admissible only if its obstruction is an exact coboundary,
is annihilated by a dual certificate or stop, is periodic/torsion with labels
attached, or is emitted as a named residual such as F7/THM-572.

Tournament Analysis:
  vertices: cocycle carriers / proof obligations, not runners.
  pairwise observable: retained predicate payload and obstruction control.
  switch/gauge: dimensionwise majority over the retention vector; ties use
      weighted score and then a declared Hamiltonian path.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations


RETENTION_KEYS = (
    "predicate",
    "fiber",
    "exactness",
    "closedness",
    "torsion_period",
    "residual_name",
    "formal_check",
    "anti_scalar",
    "computable_next",
)

WEIGHTS = {
    "predicate": 7,
    "fiber": 6,
    "exactness": 6,
    "closedness": 5,
    "torsion_period": 4,
    "residual_name": 7,
    "formal_check": 5,
    "anti_scalar": 6,
    "computable_next": 4,
}

PREDICATE_LABELS = (
    "lrc_predicate",
    "packet_fiber",
    "endpoint_owner",
    "exact_scale",
    "phase_period",
    "boundary_topology",
    "mixed_haar_sign",
    "dual_certificate",
    "family_compression",
    "state_lift_route",
    "formal_check",
    "residual_name",
)


@dataclass(frozen=True)
class CocycleCarrier:
    name: str
    kind: str
    cochain_level: str
    domain: str
    local_object: str
    differential_test: str
    exactness_route: str
    residual_if_not_exact: str
    anchors: tuple[str, ...]
    preserves: frozenset[str]
    retention: dict[str, int]

    def weighted_score(self) -> int:
        return sum(WEIGHTS[key] * self.retention[key] for key in RETENTION_KEYS)

    def destroys(self) -> tuple[str, ...]:
        return tuple(label for label in PREDICATE_LABELS if label not in self.preserves)


def fs(*items: str) -> frozenset[str]:
    return frozenset(items)


CARRIERS = [
    CocycleCarrier(
        name="fixed_margin_haar_zeta",
        kind="closed_2_cocycle",
        cochain_level="C2",
        domain="2x2 fixed-margin Haar/tournament switch",
        local_object="zeta(T)=T00-T01-T10+T11",
        differential_test="fixed row/column margins make the switch direction invisible unless zeta is retained",
        exactness_route="exact only after a packet potential reconstructs the mixed sign or a dual discrepancy cancels it",
        residual_if_not_exact="unpaired mixed Haar coefficient -> F7/THM-572",
        anchors=("HYP-2991", "HYP-2989", "HYP-2595"),
        preserves=fs(
            "lrc_predicate",
            "packet_fiber",
            "mixed_haar_sign",
            "boundary_topology",
            "dual_certificate",
            "formal_check",
            "residual_name",
        ),
        retention=dict(
            predicate=3,
            fiber=3,
            exactness=2,
            closedness=3,
            torsion_period=1,
            residual_name=3,
            formal_check=3,
            anti_scalar=3,
            computable_next=3,
        ),
    ),
    CocycleCarrier(
        name="haar_tile_interaction_curl",
        kind="typed_2_cocycle",
        cochain_level="C2",
        domain="dyadic Haar rectangle product atlas",
        local_object="orthogonal zero / same tile / owner strip / cross handoff / nested refinement",
        differential_test="product signs balance globally; labelled packet asymmetry is the curl to detect",
        exactness_route="owner strips and nested refinements become exact after the handoff arrow names owner/scale",
        residual_if_not_exact="zero-open packet with no typed coefficient -> same-tile stop or state-lift atom",
        anchors=("HYP-2992", "HYP-2993", "HYP-2988"),
        preserves=fs(
            "lrc_predicate",
            "packet_fiber",
            "endpoint_owner",
            "boundary_topology",
            "mixed_haar_sign",
            "family_compression",
            "residual_name",
        ),
        retention=dict(
            predicate=3,
            fiber=3,
            exactness=2,
            closedness=3,
            torsion_period=1,
            residual_name=3,
            formal_check=2,
            anti_scalar=3,
            computable_next=3,
        ),
    ),
    CocycleCarrier(
        name="certificate_handoff_delta",
        kind="1_coboundary_atlas",
        cochain_level="C1",
        domain="HYP-2987 zipper handoff arrows",
        local_object="delta(packet certificate potential) along proof carriers",
        differential_test="each handoff must preserve the LRC predicate and name any coordinate it drops",
        exactness_route="a compatible certificate arrow is a coboundary between packet potentials",
        residual_if_not_exact="open arrow O1-O6, especially state-lift/F7 definition debt",
        anchors=("HYP-2987", "HYP-2990", "HYP-2993"),
        preserves=fs(
            "lrc_predicate",
            "packet_fiber",
            "endpoint_owner",
            "exact_scale",
            "phase_period",
            "dual_certificate",
            "family_compression",
            "state_lift_route",
            "formal_check",
            "residual_name",
        ),
        retention=dict(
            predicate=3,
            fiber=3,
            exactness=3,
            closedness=2,
            torsion_period=2,
            residual_name=3,
            formal_check=3,
            anti_scalar=3,
            computable_next=3,
        ),
    ),
    CocycleCarrier(
        name="tope_cocircuit_boundary_class",
        kind="boundary_coboundary_or_stop",
        cochain_level="C1/C2",
        domain="oriented endpoint arrangement on R/Z",
        local_object="open tope boundary and zero-dimensional cocircuit",
        differential_test="boundary of a safe open cell is AP/GW-style owner-labelled equality",
        exactness_route="open topes are exact witness sections; cocircuits are declared stops",
        residual_if_not_exact="no-tope/no-cocircuit forbidden wall packet",
        anchors=("HYP-2986", "HYP-2951", "HYP-2948"),
        preserves=fs(
            "lrc_predicate",
            "endpoint_owner",
            "boundary_topology",
            "dual_certificate",
            "formal_check",
            "residual_name",
        ),
        retention=dict(
            predicate=3,
            fiber=2,
            exactness=3,
            closedness=3,
            torsion_period=1,
            residual_name=3,
            formal_check=3,
            anti_scalar=3,
            computable_next=2,
        ),
    ),
    CocycleCarrier(
        name="exposure_kernel_cech_class",
        kind="cech_1_cocycle",
        cochain_level="C1",
        domain="proof-channel exposure poset",
        local_object="local exposure covers glued over packet families",
        differential_test="all local exposure sections must agree on overlaps or emit UNEXPOSED_SOURCE_KERNEL",
        exactness_route="familywise no-hidden-kernel theorem makes the Cech class vanish",
        residual_if_not_exact="UNEXPOSED_SOURCE_KERNEL",
        anchors=("HYP-2988", "HYP-2963", "HYP-2981"),
        preserves=fs(
            "lrc_predicate",
            "packet_fiber",
            "dual_certificate",
            "family_compression",
            "formal_check",
            "residual_name",
        ),
        retention=dict(
            predicate=3,
            fiber=3,
            exactness=2,
            closedness=3,
            torsion_period=1,
            residual_name=3,
            formal_check=3,
            anti_scalar=3,
            computable_next=3,
        ),
    ),
    CocycleCarrier(
        name="fejer_interval_dual_coboundary",
        kind="dual_coboundary",
        cochain_level="C1",
        domain="Fejer/Toeplitz interval certificates",
        local_object="negative interval certificate as a dual potential difference",
        differential_test="exact center, atom bank, and packet key must match the primal Haar/tope packet",
        exactness_route="interval backend proves the dual certificate is exact on a packet family",
        residual_if_not_exact="formal interval backend or family-compression debt",
        anchors=("HYP-2981", "HYP-2974", "HYP-2963"),
        preserves=fs(
            "lrc_predicate",
            "packet_fiber",
            "exact_scale",
            "dual_certificate",
            "family_compression",
            "formal_check",
            "residual_name",
        ),
        retention=dict(
            predicate=3,
            fiber=3,
            exactness=3,
            closedness=2,
            torsion_period=1,
            residual_name=2,
            formal_check=3,
            anti_scalar=3,
            computable_next=3,
        ),
    ),
    CocycleCarrier(
        name="ramanujan_period_torsor",
        kind="torsion_1_cocycle",
        cochain_level="C1",
        domain="Ramanujan exact-period packets",
        local_object="primitive-period phase torsor with c_q labels",
        differential_test="squarefree/gcd quotients must not erase prime-power exact periods",
        exactness_route="exact-period projector splits phase torsion before smoothing or Fejer handoff",
        residual_if_not_exact="prime-power side channel for q=25,27,36,63,84,98,168,280,4312",
        anchors=("HYP-2979", "HYP-2978", "HYP-2982"),
        preserves=fs(
            "lrc_predicate",
            "packet_fiber",
            "endpoint_owner",
            "phase_period",
            "dual_certificate",
            "formal_check",
            "residual_name",
        ),
        retention=dict(
            predicate=3,
            fiber=3,
            exactness=2,
            closedness=2,
            torsion_period=3,
            residual_name=3,
            formal_check=2,
            anti_scalar=3,
            computable_next=2,
        ),
    ),
    CocycleCarrier(
        name="smoothing_homotopy_gauge",
        kind="gauge_1_cochain_with_boundary",
        cochain_level="C1",
        domain="kernel homotopy and admissible smoothing",
        local_object="kernel deformation gauge plus boundary-defect atom",
        differential_test="homotopy is safe only inside strict open support or after boundary defect is named",
        exactness_route="open support radius makes the gauge exact on a packet; zero-open emits a stop",
        residual_if_not_exact="failed smoothing policy -> state-lift obligation",
        anchors=("HYP-2984", "HYP-2985", "HYP-2983"),
        preserves=fs(
            "lrc_predicate",
            "packet_fiber",
            "endpoint_owner",
            "exact_scale",
            "boundary_topology",
            "state_lift_route",
            "residual_name",
        ),
        retention=dict(
            predicate=3,
            fiber=2,
            exactness=3,
            closedness=2,
            torsion_period=2,
            residual_name=3,
            formal_check=2,
            anti_scalar=3,
            computable_next=2,
        ),
    ),
    CocycleCarrier(
        name="boundary_moment_curl",
        kind="moment_2_cocycle",
        cochain_level="C2",
        domain="boundary-moment / danger-count ledgers",
        local_object="curl of endpoint moments around covering packets",
        differential_test="scalar moment shadows must agree with endpoint-owner and source-spectrum labels",
        exactness_route="moment dual closes only after packet ownership is attached",
        residual_if_not_exact="covering boundary-moment kernel or F7 harmonic sector",
        anchors=("HYP-2969", "HYP-2956", "HYP-2954"),
        preserves=fs(
            "lrc_predicate",
            "packet_fiber",
            "endpoint_owner",
            "boundary_topology",
            "family_compression",
            "state_lift_route",
            "residual_name",
        ),
        retention=dict(
            predicate=3,
            fiber=2,
            exactness=2,
            closedness=3,
            torsion_period=1,
            residual_name=3,
            formal_check=2,
            anti_scalar=3,
            computable_next=2,
        ),
    ),
    CocycleCarrier(
        name="color_resonance_square",
        kind="global_2_cocycle",
        cochain_level="C2",
        domain="colored discrepancy / CRT resonance",
        local_object="color-compatible mixed-mode square",
        differential_test="mixed Haar signs must pair with compatible color/CRT residue data",
        exactness_route="global resonance cancels a local mixed cocycle by signed pairing",
        residual_if_not_exact="unpaired color-incompatible mixed mode",
        anchors=("HYP-2595", "HYP-2594", "HYP-2991"),
        preserves=fs(
            "lrc_predicate",
            "packet_fiber",
            "phase_period",
            "mixed_haar_sign",
            "dual_certificate",
            "residual_name",
        ),
        retention=dict(
            predicate=3,
            fiber=2,
            exactness=2,
            closedness=3,
            torsion_period=2,
            residual_name=3,
            formal_check=2,
            anti_scalar=3,
            computable_next=2,
        ),
    ),
    CocycleCarrier(
        name="source_spectrum_pullback",
        kind="observer_1_cocycle",
        cochain_level="C1",
        domain="observer-source tournament movie",
        local_object="phase-indexed source pullback over Farey binding nodes",
        differential_test="source must survive pullback through phase, exact M, endpoint, and packet labels",
        exactness_route="q-witness or positive source interval makes the pullback exact",
        residual_if_not_exact="adaptive qdiv>14 source-kernel/state-lift debt",
        anchors=("HYP-2953", "HYP-2486", "THM-523"),
        preserves=fs(
            "lrc_predicate",
            "packet_fiber",
            "endpoint_owner",
            "exact_scale",
            "phase_period",
            "boundary_topology",
            "state_lift_route",
            "residual_name",
        ),
        retention=dict(
            predicate=3,
            fiber=2,
            exactness=2,
            closedness=2,
            torsion_period=2,
            residual_name=3,
            formal_check=2,
            anti_scalar=3,
            computable_next=2,
        ),
    ),
    CocycleCarrier(
        name="apex_sheaf_gluing_class",
        kind="cech_1_cocycle",
        cochain_level="C1",
        domain="apex-lift certificate sheaf",
        local_object="local sections over mod-7/mod-27 charts",
        differential_test="sections glue on overlaps or failed gluing must force positive measure",
        exactness_route="local cheap-pair sections agree after endpoint-owner labels are kept",
        residual_if_not_exact="Vitali/AP-alignment gluing obstruction",
        anchors=("HYP-2101", "HYP-2104", "THM-397", "THM-398"),
        preserves=fs(
            "lrc_predicate",
            "endpoint_owner",
            "boundary_topology",
            "state_lift_route",
            "residual_name",
        ),
        retention=dict(
            predicate=2,
            fiber=2,
            exactness=2,
            closedness=3,
            torsion_period=1,
            residual_name=2,
            formal_check=2,
            anti_scalar=3,
            computable_next=1,
        ),
    ),
    CocycleCarrier(
        name="octahedral_hodge_current",
        kind="harmonic_2_residual",
        cochain_level="H2",
        domain="octahedral/Clebsch/K33 state-lift current",
        local_object="divergence-free curl class on state-lift faces",
        differential_test="if local packet exits fail, the residual must carry a harmonic state-lift address",
        exactness_route="exact only after explicit HYP-2908/THM-572 state lift is built",
        residual_if_not_exact="named harmonic F7 sector",
        anchors=("HYP-2887", "HYP-2908", "THM-572"),
        preserves=fs(
            "lrc_predicate",
            "packet_fiber",
            "boundary_topology",
            "family_compression",
            "state_lift_route",
            "formal_check",
            "residual_name",
        ),
        retention=dict(
            predicate=3,
            fiber=2,
            exactness=1,
            closedness=3,
            torsion_period=2,
            residual_name=3,
            formal_check=2,
            anti_scalar=3,
            computable_next=2,
        ),
    ),
    CocycleCarrier(
        name="fixed_margin_johnson_F7",
        kind="named_harmonic_residual",
        cochain_level="H2",
        domain="fixed-margin count sector plus Johnson-like harmonic sector",
        local_object="packet-preserved harmonic residual after count-sector reduction",
        differential_test="count sector may close while harmonic sector remains",
        exactness_route="define the harmonic coordinate and prove it vanishes or constructs state lift",
        residual_if_not_exact="F7 as a named Johnson-harmonic packet",
        anchors=("HYP-2993", "HYP-2987", "HYP-2956", "arXiv:2606.22636"),
        preserves=fs(
            "lrc_predicate",
            "packet_fiber",
            "family_compression",
            "state_lift_route",
            "formal_check",
            "residual_name",
        ),
        retention=dict(
            predicate=3,
            fiber=3,
            exactness=1,
            closedness=3,
            torsion_period=2,
            residual_name=3,
            formal_check=2,
            anti_scalar=3,
            computable_next=2,
        ),
    ),
    CocycleCarrier(
        name="raw_scalar_shadow",
        kind="negative_control",
        cochain_level="C0-shadow",
        domain="raw K, row/column margins, qdiv, density, or component count",
        local_object="visible scalar after quotienting all packet labels",
        differential_test="fails whenever two fibers share the scalar but differ by zeta, owner, period, or state route",
        exactness_route="none; must be lifted to a labelled carrier before theorem use",
        residual_if_not_exact="false quotient / hidden cocycle",
        anchors=("HYP-2990", "HYP-2991", "HYP-2978"),
        preserves=fs("lrc_predicate"),
        retention=dict(
            predicate=1,
            fiber=0,
            exactness=0,
            closedness=0,
            torsion_period=0,
            residual_name=0,
            formal_check=1,
            anti_scalar=0,
            computable_next=1,
        ),
    ),
]

TIE_PATH = [carrier.name for carrier in CARRIERS]


COCHAIN_LEDGER = {
    "C0": [
        "packet labels",
        "endpoint-owner potentials",
        "safe-section choices",
        "Fejer rational centers",
        "exact-period residues",
        "state-lift route names",
    ],
    "C1": [
        "certificate handoff arrows",
        "endpoint owner transfers",
        "smoothing/kernel gauges",
        "Ramanujan phase shifts",
        "source-spectrum pullbacks",
        "local section restrictions",
    ],
    "C2": [
        "fixed-margin Haar zeta squares",
        "Haar tile curls",
        "tope/cocircuit boundary curls",
        "color-resonance squares",
        "boundary moment curls",
        "octahedral/K33 state-lift faces",
    ],
    "H2": [
        "unpaired mixed Haar mode",
        "no-hidden-kernel survivor",
        "Johnson-harmonic F7 sector",
        "THM-572 state-lift obstruction",
    ],
}


OBSTRUCTION_FACES = [
    (
        "fixed_margin_square",
        "closed local 2-cocycle",
        "retain zeta or pair by color-compatible discrepancy",
        "HYP-2991",
    ),
    (
        "haar_fejer_dual_square",
        "primal/dual coboundary check",
        "match typed Haar packet to interval Fejer atom bank",
        "HYP-2993/HYP-2981",
    ),
    (
        "endpoint_tope_boundary",
        "boundary coboundary or declared stop",
        "open tope gives witness; AP/GW cocircuit is a stop",
        "HYP-2986",
    ),
    (
        "period_smoothing_square",
        "torsion gauge",
        "keep exact-period labels before smoothing or Selberg quotient",
        "HYP-2979/HYP-2985",
    ),
    (
        "source_exposure_cech_face",
        "Cech 1-cocycle",
        "local exposure sections glue or emit UNEXPOSED_SOURCE_KERNEL",
        "HYP-2988/HYP-2953",
    ),
    (
        "octahedral_harmonic_face",
        "harmonic residual",
        "failed local exits must construct HYP-2908/THM-572 state lift",
        "HYP-2887/THM-572",
    ),
]


def majority_margin(a: CocycleCarrier, b: CocycleCarrier) -> int:
    return sum(
        (a.retention[k] > b.retention[k]) - (a.retention[k] < b.retention[k])
        for k in RETENTION_KEYS
    )


def orient(a: CocycleCarrier, b: CocycleCarrier) -> tuple[str, str]:
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


def tournament_edges() -> dict[tuple[str, str], int]:
    edges: dict[tuple[str, str], int] = {}
    for a, b in combinations(CARRIERS, 2):
        winner, loser = orient(a, b)
        edges[(winner, loser)] = majority_margin(a, b)
    return edges


def score_hist(edges: dict[tuple[str, str], int]) -> tuple[dict[int, int], dict[str, int]]:
    scores = {carrier.name: 0 for carrier in CARRIERS}
    for winner, _ in edges:
        scores[winner] += 1
    hist = Counter(scores.values())
    return dict(sorted(hist.items())), dict(sorted(scores.items()))


def directed_3cycles(edges: dict[tuple[str, str], int]) -> int:
    names = [carrier.name for carrier in CARRIERS]
    edge_set = set(edges)
    total = 0
    for a, b, c in combinations(names, 3):
        if (
            ((a, b) in edge_set and (b, c) in edge_set and (c, a) in edge_set)
            or ((a, c) in edge_set and (c, b) in edge_set and (b, a) in edge_set)
        ):
            total += 1
    return total


def scc_sizes(edges: dict[tuple[str, str], int]) -> list[int]:
    names = [carrier.name for carrier in CARRIERS]
    out = {name: set() for name in names}
    rev = {name: set() for name in names}
    for winner, loser in edges:
        out[winner].add(loser)
        rev[loser].add(winner)

    seen: set[str] = set()
    order: list[str] = []

    def dfs(v: str) -> None:
        seen.add(v)
        for w in out[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for name in names:
        if name not in seen:
            dfs(name)

    seen.clear()
    sizes: list[int] = []

    def rdfs(v: str, bag: list[str]) -> None:
        seen.add(v)
        bag.append(v)
        for w in rev[v]:
            if w not in seen:
                rdfs(w, bag)

    for name in reversed(order):
        if name not in seen:
            bag: list[str] = []
            rdfs(name, bag)
            sizes.append(len(bag))
    return sorted(sizes, reverse=True)


def hamiltonian_path_count(edges: dict[tuple[str, str], int]) -> int:
    names = [carrier.name for carrier in CARRIERS]
    index = {name: i for i, name in enumerate(names)}
    out_mask = [0] * len(names)
    for winner, loser in edges:
        out_mask[index[winner]] |= 1 << index[loser]

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


def label_coverage() -> dict[str, list[str]]:
    coverage: dict[str, list[str]] = defaultdict(list)
    for carrier in CARRIERS:
        for label in carrier.preserves:
            coverage[label].append(carrier.name)
    return dict(sorted(coverage.items()))


def kind_counts() -> dict[str, int]:
    return dict(sorted(Counter(carrier.kind for carrier in CARRIERS).items()))


def level_counts() -> dict[str, int]:
    return dict(sorted(Counter(carrier.cochain_level for carrier in CARRIERS).items()))


def print_wrapped_list(prefix: str, values: tuple[str, ...] | list[str]) -> None:
    print(prefix + ", ".join(values))


def main() -> None:
    edges = tournament_edges()
    hist, scores = score_hist(edges)
    coverage = label_coverage()

    print("S166 LRC14 COCYCLE OBSTRUCTION ATLAS")
    print("status: proof-interface atlas / not an LRC14 proof")
    print("namespace: HYP-2994 / T1077")
    print("upstream signal: HYP-2992=Haar tile atlas, HYP-2993=zipper pattern atlas; this pass is the cochain/cocycle lift")
    print()

    print("COCHAIN LEDGER")
    for level, items in COCHAIN_LEDGER.items():
        print(f"{level}: " + "; ".join(items))
    print()

    print("OBSTRUCTION FACE DICTIONARY")
    for name, nature, discharge, anchor in OBSTRUCTION_FACES:
        print(f"- {name}: {nature}; discharge={discharge}; anchor={anchor}")
    print()

    print("CARRIER ATLAS")
    for i, carrier in enumerate(CARRIERS, 1):
        print(f"{i:02d}. {carrier.name}")
        print(f"    kind={carrier.kind} level={carrier.cochain_level}")
        print(f"    domain={carrier.domain}")
        print(f"    local_object={carrier.local_object}")
        print(f"    differential_test={carrier.differential_test}")
        print(f"    exactness_route={carrier.exactness_route}")
        print(f"    residual_if_not_exact={carrier.residual_if_not_exact}")
        print(f"    anchors={', '.join(carrier.anchors)}")
        print(f"    preserves={', '.join(sorted(carrier.preserves))}")
        print(f"    destroys={', '.join(carrier.destroys())}")
        print(
            "    retention="
            + str(tuple(carrier.retention[key] for key in RETENTION_KEYS))
            + f" weighted_score={carrier.weighted_score()}"
        )
    print()

    print("CLASSIFICATION COUNTS")
    print(f"kinds={kind_counts()}")
    print(f"cochain_levels={level_counts()}")
    weak_labels = {label: names for label, names in coverage.items() if len(names) <= 4}
    print("sparsely_preserved_labels=" + str({k: len(v) for k, v in weak_labels.items()}))
    print()

    print("TOURNAMENT ANALYSIS")
    print(f"vertices={len(CARRIERS)}")
    print("pairwise_observable=retained LRC predicate payload plus obstruction-control vector")
    print("switch_gauge=dimensionwise majority over predicate/fiber/exactness/closedness/torsion/residual/formal/anti-scalar/computable fields")
    print("tie_path=" + " > ".join(TIE_PATH))
    print(f"score_hist={hist}")
    print(f"directed_3cycles={directed_3cycles(edges)}")
    print(f"SCC_sizes={scc_sizes(edges)}")
    print(f"Hamiltonian_path_count={hamiltonian_path_count(edges)}")
    print("scores=" + ", ".join(f"{name}:{score}" for name, score in scores.items()))
    leaders = sorted(scores, key=lambda name: (-scores[name], -next(c.weighted_score() for c in CARRIERS if c.name == name), name))
    print("leader_path_by_score=" + " > ".join(leaders))
    print()

    print("COCYCLE READOUT")
    readout = [
        "Exact coboundaries are not optional bookkeeping: they are the only safe way to forget packet potentials.",
        "Closed but non-exact local squares are useful only when a dual certificate pairs them or a named residual receives them.",
        "Torsion/period classes are real proof data; squarefree and gcd quotients erase live prime-power packets.",
        "F7 should be defined as a preserved harmonic residual sector, not as the complement of known tests.",
        "The raw scalar shadow is the negative control: it preserves the slogan but destroys nearly every obstruction coordinate.",
    ]
    for i, item in enumerate(readout, 1):
        print(f"{i}. {item}")

    print()
    print("NEXT COMPUTABLE TESTS")
    tests = [
        "Attach zeta signatures to HYP-2963 packet rows and count independent color-compatible mixed modes.",
        "Tensor Ramanujan c_q labels with Haar tile classes to catch q=25/27/36/63/84/98/168/280/4312 side channels.",
        "Group Fejer interval certificates by identical cocycle signature before interval generation.",
        "Build a no-hidden-kernel Cech check: every exposure cover overlap glues or emits a named residual.",
        "Define an F7 record with fields (packet fiber, mixed cocycle, harmonic sector, state-lift route, failed exits).",
    ]
    for i, item in enumerate(tests, 1):
        print(f"{i}. {item}")


if __name__ == "__main__":
    main()
