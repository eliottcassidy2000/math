#!/usr/bin/env python3
"""Cocycle normal-form atlas for the current LRC14 proof stack.

This is a proof-interface computation, not an LRC proof.  It turns recent
packet, Haar, Farey, C27, Fourier/PSD, tope, and tournament methods into a
single cocycle vocabulary and then runs the repository's Tournament Analysis
protocol on those cocycle carriers.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


RETENTION_KEYS = (
    "predicate",
    "base_fiber",
    "gauge_invariance",
    "coboundary_test",
    "dual_annihilator",
    "local_to_global",
    "formalizable",
    "residual_named",
    "anti_scalar_guard",
)


@dataclass(frozen=True)
class Carrier:
    name: str
    degree: str
    base: str
    cochain: str
    cocycle_equation: str
    coboundary_or_exit: str
    lrc_pull: str
    anchors: tuple[str, ...]
    vector: tuple[int, ...]
    markers: tuple[str, ...]
    first_nonzero_family: str


CARRIERS = (
    Carrier(
        name="labelled_packet_total_cocycle",
        degree="mixed total class",
        base="LRC proof-packet complex: rows, endpoint cells, exact M/Farey node, proof exits",
        cochain="all retained labels as one packet cochain",
        cocycle_equation="delta(packet)=0 on every declared commuting proof square",
        coboundary_or_exit="fiber-constant, reconstructed, dual-annihilated, or routed to named residual",
        lrc_pull="normal form: every usable quotient must say where each forgotten direction went",
        anchors=("HYP-2996", "HYP-2995", "HYP-2994", "HYP-2992", "HYP-2990", "HYP-2991", "HYP-2987", "LTI-147"),
        vector=(5, 5, 5, 5, 5, 5, 5, 5, 5),
        markers=("labelled packet", "normal form", "zipper", "residual sector"),
        first_nonzero_family="not a residual family; this is the bookkeeping object",
    ),
    Carrier(
        name="haar_zipper_2cocycle",
        degree="2-cocycle / square curl",
        base="2 x 2 Haar or fixed-margin tournament square",
        cochain="zeta(T)=T00-T01-T10+T11",
        cocycle_equation="fixed row/column margins imply only the mixed curl can move",
        coboundary_or_exit="orthogonal cancellation, boundary stop, owner strip, cross handoff, nested descent, F7",
        lrc_pull="row/column margins are unsafe unless zeta is retained or discharged",
        anchors=("HYP-2991", "HYP-2989", "HYP-2595"),
        vector=(5, 5, 4, 5, 4, 5, 4, 5, 5),
        markers=("zeta", "Haar", "fixed-margin", "mixed Haar", "2 x 2"),
        first_nonzero_family="F_zeta: unpaired mixed Haar curl after all zipper teeth fail",
    ),
    Carrier(
        name="endpoint_owner_boundary_cocycle",
        degree="1-cocycle / boundary current",
        base="regular-open circle cells cut by danger endpoints",
        cochain="signed endpoint owner current on the tope boundary graph",
        cocycle_equation="boundary owner current sums to zero around closed equality atoms",
        coboundary_or_exit="AP/GW boundary cocircuit, open tope, or named wall residual",
        lrc_pull="separates positive strict witnesses from closed endpoint-only atoms",
        anchors=("HYP-2986", "HYP-2984", "HYP-2949"),
        vector=(5, 5, 4, 5, 3, 4, 5, 5, 5),
        markers=("endpoint", "owner", "boundary cocircuit", "tope", "regular-open"),
        first_nonzero_family="F_wall: no open tope and no AP/GW boundary cocircuit",
    ),
    Carrier(
        name="farey_excess_mediant_1cocycle",
        degree="1-cocycle / determinant defect",
        base="Farey edge path adjacent to 1/14",
        cochain="e(p/q)=14p-q and det((p,q),(p',q'))",
        cocycle_equation="mediant addition preserves primitive edge data only with determinant labels",
        coboundary_or_exit="e=0 AP/GW equality, e=1 unit-excess child, e>1 non-tight branch",
        lrc_pull="binding scale is q plus Farey excess, not product or numerator scalar alone",
        anchors=("HYP-2930", "HYP-2931", "HYP-2932", "HYP-2945"),
        vector=(5, 5, 5, 4, 3, 5, 5, 5, 5),
        markers=("Farey", "mediant", "e=14p-q", "3/41", "unit-excess"),
        first_nonzero_family="F_farey: first nonzero excess class beyond the 1/14 boundary",
    ),
    Carrier(
        name="c27_carry_lift_1cocycle",
        degree="1-cocycle / covering carry",
        base="C27 antipodal shell cover over residues",
        cochain="carry k in v=r+27k, plus gcd shell route",
        cocycle_equation="lift changes must preserve shell route and unit/nonunit monodromy",
        coboundary_or_exit="unit petal discharge, GW D3 branch, K33 D9 state-lift debt",
        lrc_pull="residue shadows lie unless the carry and gcd-shell route are retained",
        anchors=("HYP-2937", "HYP-2942", "HYP-2171"),
        vector=(5, 5, 5, 4, 3, 4, 4, 5, 5),
        markers=("C27", "carry", "shell", "D9", "H12", "unital"),
        first_nonzero_family="F_carry: nonunit C27 monodromy, usually K33/D9",
    ),
    Carrier(
        name="fejer_toeplitz_dual_coboundary",
        degree="dual coboundary / PSD annihilator",
        base="Fourier moment cone of the danger multiplicity C_S(t)-1",
        cochain="finite Toeplitz moment functional or Fejer trigonometric square",
        cocycle_equation="covering would force every Toeplitz moment matrix to be PSD",
        coboundary_or_exit="negative Fejer/Toeplitz certificate is a dual coboundary for a safe interval",
        lrc_pull="low harmonic PSD failure is a certificate, not just a numerical shadow",
        anchors=("HYP-2974", "HYP-2981", "HYP-2987"),
        vector=(5, 4, 4, 5, 5, 4, 5, 4, 5),
        markers=("Fejer", "Toeplitz", "PSD", "Fourier", "moment"),
        first_nonzero_family="F_psd: harmonic moment cone fails to certify or fails to descend",
    ),
    Carrier(
        name="ramanujan_exact_period_character_cocycle",
        degree="1-cocycle / exact-period character",
        base="denominator and divisor lattice of primitive periods",
        cochain="Ramanujan exact-period packet c_q(a-b) with Mobius inversion labels",
        cocycle_equation="period projectors commute only when prime-power and endpoint labels are kept",
        coboundary_or_exit="exact-period projector, prime-power side channel, or smoothing handoff",
        lrc_pull="divisor quotients are unsafe without endpoint, qdiv, and prime-power labels",
        anchors=("HYP-2978", "HYP-2979", "HYP-2982"),
        vector=(4, 4, 5, 4, 4, 4, 5, 4, 5),
        markers=("Ramanujan", "exact-period", "divisor", "Mobius", "prime-power"),
        first_nonzero_family="F_period: exact-period character survives divisor quotient",
    ),
    Carrier(
        name="tope_cocircuit_wall_cocycle",
        degree="oriented matroid cocycle",
        base="hyperplane arrangement of danger endpoints",
        cochain="tope sign vector and boundary cocircuit owner pair",
        cocycle_equation="crossing a wall changes one sign and records the cocircuit",
        coboundary_or_exit="open all-safe tope, AP/GW equality cocircuit, or forbidden wall packet",
        lrc_pull="closed boundary debt and open witness mass must not be scalarized together",
        anchors=("HYP-2986", "HYP-2948", "HYP-2951"),
        vector=(5, 5, 4, 5, 3, 4, 5, 5, 5),
        markers=("tope", "cocircuit", "wall", "endpoint", "owner pair"),
        first_nonzero_family="F_cocircuit: wall packet not discharged by open-tope or equality atom",
    ),
    Carrier(
        name="tournament_path_h1_cocycle",
        degree="1-cocycle modulo coboundaries",
        base="tournament path-homology or proof-carrier tournament",
        cochain="edge comparison weights; triangle curl records nontransitive pressure",
        cocycle_equation="transitive triples kill curl; directed 3-cycles expose H^1/H^2 pressure",
        coboundary_or_exit="score potential coboundary, SCC support residue, or named cyclic obstruction",
        lrc_pull="tournament iso classes help only after the quotient states what labels survived",
        anchors=("beta1_upper_bound", "HYP-2924", "HYP-2990"),
        vector=(4, 3, 4, 5, 2, 4, 4, 4, 4),
        markers=("tournament", "cocycle", "H^1", "SCC", "directed 3"),
        first_nonzero_family="F_tournament: non-potential tournament pressure class",
    ),
    Carrier(
        name="boundary_moment_multichart_cocycle",
        degree="chart-transition cocycle",
        base="adaptive exact-period boundary-moment charts",
        cochain="missed-sector depth vector across denominator charts",
        cocycle_equation="one all-covered chart is harmless unless transition data also closes",
        coboundary_or_exit="multi-chart feasible region, positive Haar image, or K33/state-lift debt",
        lrc_pull="covering residuals need chart identity and transition labels",
        anchors=("HYP-2969", "HYP-2965", "HYP-2968"),
        vector=(4, 5, 4, 3, 4, 4, 3, 5, 5),
        markers=("boundary-moment", "multi-chart", "missed", "gK8", "L_y"),
        first_nonzero_family="F_chart: boundary-moment transition class fails to close",
    ),
    Carrier(
        name="state_lift_obstruction_class",
        degree="named residual cohomology class",
        base="TournamentStateLift / THM-572 residual sector",
        cochain="packet value agreeing with forbidden H=7 lift",
        cocycle_equation="all lower exits failed and the residual must be represented in the lift",
        coboundary_or_exit="construct forbidden state lift or split the residual into earlier cocycles",
        lrc_pull="F7 is not an exception bucket; it must be a named obstruction class",
        anchors=("THM-572", "HYP-2908", "HYP-2963", "LTI-039", "HYP-2987", "HYP-2991"),
        vector=(5, 4, 5, 3, 3, 5, 5, 5, 5),
        markers=("state-lift", "F7", "THM-572", "H=7", "TournamentStateLift"),
        first_nonzero_family="F_lift: irreducible named state-lift obstruction",
    ),
    Carrier(
        name="curried_section_derivative_cocycle",
        degree="1-cocycle / section derivative",
        base="curried proof maps after freezing runners, endpoints, modes, or packet coordinates",
        cochain="difference between adjacent sections of a partially applied proof function",
        cocycle_equation="section derivatives must commute around frozen-coordinate squares",
        coboundary_or_exit="reconstruct forgotten coordinate, prove section independence, or expose defect",
        lrc_pull="curried functions are useful exactly when their section derivative is controlled",
        anchors=("HYP-2990", "LRC technique index", "curried-functions thread"),
        vector=(4, 4, 5, 4, 3, 4, 3, 4, 5),
        markers=("curried", "section", "derivative", "coordinate", "quotient"),
        first_nonzero_family="F_section: frozen-coordinate derivative refuses to commute",
    ),
    Carrier(
        name="raw_scalar_shadow",
        degree="0-shadow / negative control",
        base="unlabelled scalar summaries",
        cochain="one number such as residue class count, component count, or product value",
        cocycle_equation="none without an external theorem identifying the kernel",
        coboundary_or_exit="diagnostic only after cocycle channels are already discharged",
        lrc_pull="explains why residue, product, component, and iso-class shadows keep lying",
        anchors=("HYP-2990", "HYP-2991", "HYP-2978"),
        vector=(2, 1, 1, 1, 1, 1, 1, 1, 1),
        markers=("raw scalar", "shadow", "component count", "product value", "residue"),
        first_nonzero_family="not a family; this is the anti-template",
    ),
)


DOCUMENTS = (
    "00-navigation/LRC-TECHNIQUE-INDEX.md",
    "comms/POKE-COORDINATION.md",
    "05-knowledge/hypotheses/HYP-2996-lrc14-residual-section-packet-grid-verification.md",
    "05-knowledge/hypotheses/HYP-2995-lrc14-cocycle-carrier-atlas.md",
    "05-knowledge/hypotheses/HYP-2994-lrc14-cocycle-obstruction-atlas.md",
    "05-knowledge/hypotheses/HYP-2992-lrc14-cocycle-sheaf-exactness.md",
    "05-knowledge/hypotheses/HYP-2993-lrc14-zipper-theorem-pattern-atlas.md",
    "05-knowledge/hypotheses/HYP-2991-lrc14-haar-zipper-cocycle.md",
    "05-knowledge/hypotheses/HYP-2990-abstract-zipper-theorem-atlas.md",
    "05-knowledge/hypotheses/HYP-2989-lrc14-haar-product-discrepancy-tournament-tiling.md",
    "05-knowledge/hypotheses/HYP-2987-lrc14-certificate-handoff-atlas.md",
    "05-knowledge/hypotheses/HYP-2986-lrc14-tope-wall-cocircuit.md",
    "05-knowledge/hypotheses/HYP-2985-lrc14-admissible-smoothing-dispatcher.md",
    "05-knowledge/hypotheses/HYP-2978-lrc14-ramanujan-divisor-quotient-guardrails.md",
    "05-knowledge/hypotheses/HYP-2974-lrc14-fourier-toeplitz-dual-scout.md",
    "05-knowledge/hypotheses/HYP-2969-lrc14-boundary-moment-packet-ledger.md",
    "05-knowledge/hypotheses/HYP-2963-lrc14-labelled-packet-counterexample-audit.md",
    "05-knowledge/hypotheses/HYP-2942-lrc14-c27-unital-block-lift.md",
    "05-knowledge/hypotheses/HYP-2937-lrc14-c27-shell-transfer-spectrum.md",
    "05-knowledge/hypotheses/HYP-2930-lrc14-farey-mediant-tournament-interface.md",
    "05-knowledge/results/beta1_upper_bound.out",
    "05-knowledge/results/heisenberg_frustration_bridge.out",
)


def text_for(path: str) -> str:
    file_path = ROOT / path
    if not file_path.exists():
        return ""
    return file_path.read_text(encoding="utf-8", errors="ignore")


def marker_hits(text: str, carrier: Carrier) -> int:
    lower = text.lower()
    return sum(lower.count(marker.lower()) for marker in carrier.markers)


def winner(i: int, j: int) -> int:
    """Return the winning carrier index under majority retention comparison.

    Pairwise observable: compare all coordinates in RETENTION_KEYS.
    Switch/gauge: orient toward the carrier with more strictly larger
    coordinates.  Tie Hamiltonian path: CARRIERS order.
    """
    left = CARRIERS[i].vector
    right = CARRIERS[j].vector
    left_wins = sum(a > b for a, b in zip(left, right))
    right_wins = sum(b > a for a, b in zip(left, right))
    if left_wins > right_wins:
        return i
    if right_wins > left_wins:
        return j
    return min(i, j)


def build_edges() -> dict[int, set[int]]:
    edges: dict[int, set[int]] = {i: set() for i in range(len(CARRIERS))}
    for i, j in combinations(range(len(CARRIERS)), 2):
        w = winner(i, j)
        loser = j if w == i else i
        edges[w].add(loser)
    return edges


def count_3cycles(edges: dict[int, set[int]]) -> list[tuple[str, str, str]]:
    cycles = []
    for i, j, k in combinations(range(len(CARRIERS)), 3):
        triples = ((i, j, k), (i, k, j), (j, i, k), (j, k, i), (k, i, j), (k, j, i))
        for a, b, c in triples:
            if b in edges[a] and c in edges[b] and a in edges[c]:
                cycles.append((CARRIERS[a].name, CARRIERS[b].name, CARRIERS[c].name))
                break
    return cycles


def scc_sizes(edges: dict[int, set[int]]) -> list[int]:
    n = len(CARRIERS)
    reverse = {i: set() for i in range(n)}
    for src, dsts in edges.items():
        for dst in dsts:
            reverse[dst].add(src)

    seen: set[int] = set()
    order: list[int] = []

    def dfs(v: int) -> None:
        seen.add(v)
        for nxt in edges[v]:
            if nxt not in seen:
                dfs(nxt)
        order.append(v)

    for v in range(n):
        if v not in seen:
            dfs(v)

    seen.clear()
    sizes: list[int] = []

    def rdfs(v: int, bucket: list[int]) -> None:
        seen.add(v)
        bucket.append(v)
        for nxt in reverse[v]:
            if nxt not in seen:
                rdfs(nxt, bucket)

    for v in reversed(order):
        if v not in seen:
            bucket: list[int] = []
            rdfs(v, bucket)
            sizes.append(len(bucket))
    return sorted(sizes, reverse=True)


def hamiltonian_path_count(edges: dict[int, set[int]]) -> int:
    n = len(CARRIERS)
    dp: dict[tuple[int, int], int] = {}
    for i in range(n):
        dp[(1 << i, i)] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp.get((mask, last), 0)
            if not count:
                continue
            for nxt in edges[last]:
                if not (mask & (1 << nxt)):
                    dp[(mask | (1 << nxt), nxt)] = dp.get((mask | (1 << nxt), nxt), 0) + count
    full = (1 << n) - 1
    return sum(dp.get((full, last), 0) for last in range(n))


def scan_documents() -> tuple[dict[str, int], Counter[tuple[str, str]], list[str]]:
    docs = [(path, text_for(path)) for path in DOCUMENTS]
    present_docs = [path for path, text in docs if text]
    hits: dict[str, int] = {}
    per_doc_presence: dict[str, set[str]] = {}
    for path, text in docs:
        if not text:
            continue
        carriers_here = set()
        for carrier in CARRIERS:
            count = marker_hits(text, carrier)
            hits[carrier.name] = hits.get(carrier.name, 0) + count
            if count:
                carriers_here.add(carrier.name)
        per_doc_presence[path] = carriers_here

    co = Counter()
    for carriers_here in per_doc_presence.values():
        for a, b in combinations(sorted(carriers_here), 2):
            co[(a, b)] += 1
    return hits, co, present_docs


def print_header(title: str) -> None:
    print()
    print(title)
    print("-" * len(title))


def main() -> None:
    print("S167 cocycle normal-form atlas for LRC14")
    print("=" * 54)

    print_header("[1] Exact definitions")
    print("Proof packet complex P:")
    print("  vertices: labelled proof states, not automatically runners")
    print("  edges: admissible quotient, handoff, wall-crossing, or comparison moves")
    print("  2-cells: commuting proof squares, Haar squares, tournament triples, chart overlaps")
    print("Coefficient groups:")
    print("  integers/rationals for counts and moments; signs for topes; characters for periods")
    print("Cochain:")
    print("  a k-observable a assigns a value to k-cells of P")
    print("Cocycle condition:")
    print("  delta a = 0 on all cells whose boundary preserves the target LRC predicate")
    print("Coboundary/certificate condition:")
    print("  a = delta phi on a quotient fiber, or a dual certificate annihilates a")
    print("Sound quotient condition:")
    print("  every forgotten direction is fiber-constant, reconstructed, dual-annihilated,")
    print("  or mapped to a named residual group such as F7/THM-572.")

    print_header("[2] Cocycle carriers")
    for carrier in CARRIERS:
        anchor_text = ", ".join(carrier.anchors)
        vector_text = ", ".join(f"{key}={value}" for key, value in zip(RETENTION_KEYS, carrier.vector))
        print(f"- {carrier.name} [{carrier.degree}]")
        print(f"  base: {carrier.base}")
        print(f"  cochain: {carrier.cochain}")
        print(f"  delta law: {carrier.cocycle_equation}")
        print(f"  exit: {carrier.coboundary_or_exit}")
        print(f"  LRC pull: {carrier.lrc_pull}")
        print(f"  first nonzero family: {carrier.first_nonzero_family}")
        print(f"  anchors: {anchor_text}")
        print(f"  retention vector: {vector_text}")

    edges = build_edges()
    scores = {CARRIERS[i].name: len(edges[i]) for i in range(len(CARRIERS))}
    score_hist = Counter(scores.values())
    cycles = count_3cycles(edges)

    print_header("[3] Tournament Analysis")
    print("vertices: cocycle channels / proof obligations, not runners")
    print("pairwise observable: coordinate-wise retention comparison")
    print("switch/gauge: majority of strictly larger retention coordinates")
    print("tie Hamiltonian path:")
    print("  " + " > ".join(carrier.name for carrier in CARRIERS))
    for name, score in sorted(scores.items(), key=lambda item: (-item[1], item[0])):
        carrier = next(c for c in CARRIERS if c.name == name)
        print(f"  {name:44s} score={score:2d} vector={carrier.vector}")
    print(f"score_hist={dict(sorted(score_hist.items()))}")
    print(f"directed_3cycles={len(cycles)}")
    if cycles:
        print("sample directed 3-cycles:")
        for a, b, c in cycles[:8]:
            print(f"  {a} -> {b} -> {c} -> {a}")
    print(f"SCC_sizes={scc_sizes(edges)}")
    print(f"Hamiltonian_path_count={hamiltonian_path_count(edges)}")

    hits, co, present_docs = scan_documents()
    print_header("[4] Repository co-occurrence scan")
    print(f"documents requested={len(DOCUMENTS)} present={len(present_docs)}")
    print("documents scanned:")
    for path in present_docs:
        print(f"  {path}")
    print("carrier marker hits:")
    for name, count in sorted(hits.items(), key=lambda item: (-item[1], item[0])):
        print(f"  {name:44s} {count}")
    print("strongest co-occurrences by document count:")
    for (a, b), count in co.most_common(15):
        print(f"  {a:42s} <-> {b:42s} docs={count}")

    print_header("[5] AP/GW cocycle closure profile")
    closure_rows = (
        ("qdiv/Farey", "e=14p-q is zero at the tight 1/14 boundary"),
        ("endpoint/tope", "regular-open safe mass is zero; only boundary cocircuit debt remains"),
        ("Haar zipper", "mixed zeta is stopped at a boundary atom, not a free slider"),
        ("C27/Jacobsthal", "AP is base tiling; GW is the unique gated D3 acceleration"),
        ("Fejer/Toeplitz", "PSD-failure certificates vanish because there is no positive safe interval"),
        ("tournament", "raw iso class is insufficient; labelled optimum/Farey packet is the useful class"),
        ("state lift", "no F7/THM-572 residual is invoked for AP/GW"),
    )
    for label, statement in closure_rows:
        print(f"  {label:18s}: {statement}")

    print_header("[6] Counterexample family normal form")
    print("A primitive LRC14 counterexample would need a first nonzero class:")
    for carrier in CARRIERS:
        if carrier.name in {"labelled_packet_total_cocycle", "raw_scalar_shadow"}:
            continue
        print(f"  {carrier.first_nonzero_family}")
    print("If every listed class is a coboundary, dual-annihilated, or exits to AP/GW,")
    print("the packet has no remaining typed place to hide; the only legal last bucket")
    print("is the named state_lift_obstruction_class.")

    print_header("[7] Assumption challenge")
    print("alternate vertex sets considered:")
    print("  runners, arcs, fixed circle sections, section boundaries, endpoint walls,")
    print("  residues, cover arcs, Fourier modes, Haar tiles, divisor packets,")
    print("  matroid cocircuits, chart overlaps, proof obligations, and cocycle channels.")
    print("chosen vertices:")
    print("  cocycle channels / proof obligations.")
    print("preserved predicate:")
    print("  strict witness, AP/GW equality, positive certificate, or named state-lift residual.")
    print("destroyed information:")
    print("  raw runner identity and continuous scalar order are intentionally destroyed;")
    print("  exact M/Farey, endpoint owners, carry, Haar curl, period labels, and residual")
    print("  class are retained because they are the current load-bearing cocycles.")
    print("challenged assumption:")
    print("  Tournament Analysis vertices do not need to be objects of the original problem;")
    print("  here the useful tournament is over obstruction classes.")

    print_header("[8] Cocycle normal-form theorem target")
    print("Candidate theorem:")
    print("  For every reduced LRC14 packet, the retained proof data form a cocycle in")
    print("  the packet complex.  If the packet has no strict safe interval, then each")
    print("  non-boundary cocycle is either a coboundary on the quotient fiber, killed")
    print("  by a dual certificate, transferred through a typed zipper tooth, or mapped")
    print("  to the named F7/THM-572 obstruction class.  AP and GW are exactly the")
    print("  zero-open boundary packets where all channels close before F7.")
    print("Next finite test:")
    print("  build HYP-2963 packet-level cocycle ledgers and record the first nonzero")
    print("  class for every low-frontier row; no row should land in raw_scalar_shadow.")


if __name__ == "__main__":
    main()
