#!/usr/bin/env python3
"""
tournament_contradiction_grammar_codex_s260.py

Build a programmatic grammar around the project's tournament proof-by-
contradiction pattern.  The seed pattern is:

    forced encoding -> tournament invariant -> forbidden value/shape -> contradiction

The script imports the S31ah certificate engine, integrates the concurrent
HYP-3099 tournament-as-proof applications, adds a wider menu of tournament
proof moves, and scores those moves against current repo targets: LRC14
observer gluing, the H=7/H=21 route, level-7 descent, HYP-2963 packet
normalization, and the sexy-prime sidecar question.

This is a research routing artifact.  It does not claim a new theorem.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
import itertools
import re
from typing import Iterable

import tournament_certificate_engine_kps as kps


STOPWORDS = {
    "a",
    "an",
    "and",
    "as",
    "at",
    "be",
    "by",
    "for",
    "from",
    "if",
    "in",
    "into",
    "is",
    "it",
    "of",
    "on",
    "or",
    "the",
    "to",
    "with",
}


@dataclass(frozen=True)
class Technique:
    name: str
    family: str
    mode: str
    certificate: str
    preserves: tuple[str, ...]
    destroys: tuple[str, ...]
    vertices: tuple[str, ...]
    trigger_keywords: tuple[str, ...]
    applications: tuple[str, ...]
    next_experiment: str


@dataclass(frozen=True)
class Target:
    name: str
    key: str
    domain: str
    claim_type: str
    encoded_as: tuple[str, ...]
    desired_signal: tuple[str, ...]
    sidecars: tuple[str, ...]
    notes: str


def tokens(parts: Iterable[str]) -> set[str]:
    text = " ".join(parts).lower()
    return {w for w in re.findall(r"[a-z0-9_]+", text) if w not in STOPWORDS}


def make_techniques() -> list[Technique]:
    return [
        Technique(
            "H forbidden-value certificate",
            "forbidden_value",
            "contradiction",
            "If a forced tournament shadow has H even, H=7, or H=21, the encoding is impossible.",
            ("Hamiltonian path count", "odd-cycle conflict independence polynomial", "SCC product factorization"),
            ("geometry of the original LRC row", "which cycle witnesses produced the value"),
            ("tournament isomorphism classes", "Omega components", "proof obligations"),
            ("H", "7", "21", "odd", "Omega", "SCC", "Hamiltonian", "contradiction"),
            ("h21_closure", "level7_state_lift", "observer_gluing", "loose_digraph_guardrail"),
            "Feed every terminal state-lift quotient through H_spectrum_certificate before promoting it.",
        ),
        Technique(
            "Redei odd-path parity",
            "forbidden_value",
            "contradiction",
            "The Hamiltonian path count of a tournament is always odd.",
            ("Hamiltonian paths", "linear extensions of forced comparison data"),
            ("local metric sizes", "endpoint owners"),
            ("comparison states", "proof-route choices", "sieve channels"),
            ("Redei", "parity", "odd", "Hamiltonian", "path", "count"),
            ("sexy_prime_sidecar", "observer_gluing", "h21_closure"),
            "Use it as a parity sanity check whenever a count is asserted to be a tournament path count.",
        ),
        Technique(
            "Landau score majorization",
            "score_polytope",
            "nonrealizability",
            "A proposed win-score vector must satisfy Landau prefix inequalities.",
            ("pairwise comparison margins", "score histograms", "degree budgets"),
            ("cycle structure", "metric scale"),
            ("proof obligations", "residue channels", "sidecar columns"),
            ("score", "Landau", "majorization", "degree", "win", "profile"),
            ("observer_gluing", "pascal_cap", "sexy_prime_sidecar", "a000568_defect"),
            "Turn sidecar dominance claims into score sequences and reject non-Landau shadows.",
        ),
        Technique(
            "Omega forbidden-shape miner",
            "Omega_realizability",
            "contradiction",
            "A graph forced to be Omega(T) must be in the tournament odd-cycle conflict spectrum.",
            ("odd-cycle sharing", "independence polynomial", "local obstruction shape"),
            ("absolute runner scale", "which original arc created a cycle"),
            ("Omega vertices", "cycle witnesses", "conflict pairs"),
            ("Omega", "conflict", "cycle", "P4", "K3", "realizable", "shape"),
            ("h21_closure", "level7_state_lift", "a000568_defect", "loose_digraph_guardrail"),
            "Extend THM-200/202 by enumerating small forbidden Omega components and matching repo residuals.",
        ),
        Technique(
            "Alpha-vector Newton gate",
            "independence_polynomial",
            "sidecar obstruction",
            "Small Omega independence polynomials observed in tournaments obey strong log-concavity patterns.",
            ("alpha-vector", "independence polynomial", "packing count"),
            ("cycle labels", "metric normalization"),
            ("independent cycle packs", "Omega components"),
            ("alpha", "Newton", "log", "concavity", "independence", "polynomial"),
            ("h21_closure", "pascal_cap", "symbolic_overlap"),
            "Mine alpha-vectors up to the feasible enumeration limit and treat persistent failures as candidates.",
        ),
        Technique(
            "Odd-cycle census gap miner",
            "cycle_spectrum",
            "nonrealizability",
            "The tuple of directed 3-, 5-, and 7-cycle counts has an achieved spectrum with gaps.",
            ("directed odd-cycle counts", "cycle census", "local cyclicity"),
            ("which pairwise certificate caused the cycle",),
            ("small tournaments", "winding residues", "observer charts"),
            ("cycle", "census", "c3", "c5", "c7", "winding", "residue"),
            ("coarse_winding", "fine_mod7_winding", "a000568_defect"),
            "Use it on fine mod-p winding tournaments after resolving antipodal ties.",
        ),
        Technique(
            "Odd automorphism obstruction",
            "symmetry",
            "contradiction",
            "Finite tournament automorphism groups have odd order, so forced involutive symmetry is suspect.",
            ("visible automorphism parity", "observer-cut orbit ledger"),
            ("metric slack", "nonvisible symmetries"),
            ("rooted perspectives", "observer cuts", "quotient fibers"),
            ("automorphism", "involution", "orbit", "symmetry", "observer", "cut"),
            ("a000568_defect", "observer_gluing", "packet_normalizer"),
            "Search quotient fibers for fake order-2 symmetries created by forgotten observer cuts.",
        ),
        Technique(
            "Skew-spectrum interlacing",
            "linear_algebra",
            "sidecar obstruction",
            "Skew-adjacency spectra of tournaments are constrained and interlace under deletion.",
            ("linear spectral signature", "deletion interlacing", "regularity alarms"),
            ("combinatorial witness labels",),
            ("adjacency matrices", "state-lift matrices", "observer columns"),
            ("spectrum", "eigenvalue", "skew", "matrix", "interlacing", "regular"),
            ("matrix_observability", "fine_mod7_winding", "packet_normalizer"),
            "Attach skew-spectrum rows to HYP-2963 quotient columns as a cheap non-scalar sidecar.",
        ),
        Technique(
            "H-max rigidity",
            "extremal",
            "rigidity",
            "If a quotient reaches an H-extreme, it should inherit regular/Paley-like constraints.",
            ("extremal H", "regularity", "Paley/BIBD comparison"),
            ("non-extremal residual geometry",),
            ("regular tournaments", "fine prime residue charts"),
            ("H", "max", "regular", "Paley", "BIBD", "extreme", "prime"),
            ("fine_mod7_winding", "level7_state_lift", "coarse_winding"),
            "Use mod-7 fine charts to see whether a residual is forced toward Paley T7 or away from it.",
        ),
        Technique(
            "SCC product descent",
            "factorization",
            "contradiction",
            "H factors over strongly connected components, so forbidden prime factors can localize debt.",
            ("SCC factorization", "component product", "local debt split"),
            ("arcs between SCCs", "exact geometry inside components"),
            ("strong components", "proof-route factors", "packet subledgers"),
            ("SCC", "product", "factor", "component", "7", "21", "debt"),
            ("h21_closure", "level7_state_lift", "packet_normalizer"),
            "Split any terminal H value by SCC before treating it as a single global obstruction.",
        ),
        Technique(
            "Trienerment lift for ties",
            "encoding_repair",
            "repair",
            "When a coarse winding map creates ties, lift to ternary/tie-aware or fine-scale data before H tests.",
            ("tie status", "fine-scale residue", "nondegenerate comparison"),
            ("coarse scalar H", "raw mod-14 shadow"),
            ("tied pairs", "residue packets", "fine mod-p charts"),
            ("tie", "trienerment", "winding", "fine", "mod7", "coarse", "antipodal"),
            ("coarse_winding", "fine_mod7_winding", "observer_gluing"),
            "Replace coarse mod-14 H by a tie-aware chart before using forbidden-value certificates.",
        ),
        Technique(
            "Proof-obligation dominance tournament",
            "meta_tournament",
            "proof_routing",
            "Treat proof obligations as vertices and orient by retained predicate plus discharge cost.",
            ("LRC predicate retention", "destroyed-coordinate ledger", "terminal exit"),
            ("raw theorem-route comfort",),
            ("proof obligations", "observer charts", "sidecar bundles"),
            ("proof", "obligation", "observer", "gluing", "sidecar", "terminal", "exit"),
            ("observer_gluing", "packet_normalizer", "pascal_cap", "sexy_prime_sidecar"),
            "Use tournament fingerprints to expose cycles among proof routes, not just scalar orderings.",
        ),
        Technique(
            "Median/Helly route-center tournament",
            "route_geometry",
            "repair",
            "Route triples should have legal centers after enough sidecars are retained.",
            ("median center", "route compatibility", "sidecar closure"),
            ("raw route labels",),
            ("route states", "certificate sidecars", "discharge hubs"),
            ("median", "Helly", "route", "center", "triple", "sidecar"),
            ("packet_normalizer", "observer_gluing", "branch_closure"),
            "Compare new tournament certificates against the existing median-center packet interface.",
        ),
        Technique(
            "No-naked-bridge orientation",
            "graph_orientation",
            "contradiction_or_repair",
            "A proof branch graph must not contract a load-bearing bridge without a protected reverse path.",
            ("strong orientation", "bridge protection", "reverse verification"),
            ("scalar branch count",),
            ("branch graph edges", "residual exits", "certificate routes"),
            ("bridge", "Robbins", "orientation", "branch", "protected", "reverse"),
            ("packet_normalizer", "residual_capacitor", "observer_gluing"),
            "Apply to any tournament-derived contradiction before deleting a sidecar edge.",
        ),
        Technique(
            "Metagraph nonembedding",
            "metagraph",
            "nonrealizability",
            "A claimed transition system must embed in the tournament metagraph with the right parity and flips.",
            ("edge-flip parity", "metagraph walk", "even-graph dual"),
            ("actual metric values",),
            ("tournament classes", "flip edges", "state transitions"),
            ("metagraph", "edge", "flip", "parity", "even", "walk", "dual"),
            ("a000568_defect", "packet_normalizer", "loose_digraph_guardrail"),
            "Check whether a residual transition ladder is an actual tournament-class walk.",
        ),
        Technique(
            "Transfer-matrix symmetry",
            "linear_algebra",
            "contradiction",
            "A transfer matrix forced to be tournament-symmetric cannot realize an asymmetric target table.",
            ("matrix symmetry", "two-boundary transfer count"),
            ("which internal witness realizes a transfer",),
            ("transfer states", "boundary labels", "matrix entries"),
            ("transfer", "matrix", "symmetry", "asymmetry", "boundary"),
            ("matrix_observability", "packet_normalizer", "observer_gluing"),
            "Use as a low-cost filter for observer-chart overlap tables.",
        ),
        Technique(
            "Path-homology rank sidecar",
            "homology",
            "sidecar obstruction",
            "Tournament-derived incidence can carry boundary rank/torsion constraints.",
            ("boundary rank", "cycle residue", "incidence torsion"),
            ("metric threshold",),
            ("incidence rows", "cycle residues", "boundary images"),
            ("path", "homology", "rank", "torsion", "boundary", "cycle"),
            ("residual_capacitor", "packet_normalizer", "a000568_defect"),
            "Attach only after a metric or route certificate exists; rank alone is not LRC currency.",
        ),
        Technique(
            "Residue-sieve tournament",
            "number_theory_sidecar",
            "sidecar_obstruction",
            "Local residue channels can be oriented by which channel preserves more forbidden-residue data.",
            ("finite residue word", "sieve weight", "local obstruction"),
            ("parity-breaking distribution", "global primality"),
            ("prime gap rays", "residue classes", "sieve channels"),
            ("residue", "sieve", "prime", "sexy", "gap", "parity", "distribution"),
            ("sexy_prime_sidecar", "level7_state_lift", "fine_mod7_winding"),
            "Use for sexy primes as admissibility and bookkeeping, not as a lower-bound sieve proof.",
        ),
        Technique(
            "Proof-by-survival miner",
            "negative_information",
            "necessary_condition",
            "If no certificate fires, the survivor profile is a useful necessary condition.",
            ("surviving invariant profile", "certificate failure ledger"),
            ("sufficiency",),
            ("certificate columns", "target residuals", "negative examples"),
            ("survival", "necessary", "condition", "profile", "failure", "ledger"),
            ("observer_gluing", "h21_closure", "packet_normalizer", "sexy_prime_sidecar"),
            "Record no-hit profiles as sidecar requirements instead of discarding failed attacks.",
        ),
        Technique(
            "Functorial certificate pullback",
            "category",
            "proof_routing",
            "A certificate can be pulled back along an encoding only if the predicate is fiber-constant.",
            ("predicate preservation", "fiber-constant quotient", "pullback discipline"),
            ("coordinates not represented in the functor",),
            ("quotient maps", "fiber ledgers", "certificate functors"),
            ("functor", "pullback", "fiber", "constant", "quotient", "predicate"),
            ("observer_gluing", "packet_normalizer", "coarse_winding", "sexy_prime_sidecar"),
            "Make every tournament analogy state the pullback map before using the certificate.",
        ),
        Technique(
            "Random gap miner",
            "discovery",
            "candidate_generator",
            "Sample larger tournaments for invariant gaps that survive beyond exact enumeration.",
            ("candidate gaps", "spectrum onset", "persistent forbidden values"),
            ("proof of permanence",),
            ("sampled tournaments", "invariant values", "gap tables"),
            ("random", "gap", "spectrum", "sample", "persistent", "discover"),
            ("h21_closure", "fine_mod7_winding", "packet_normalizer"),
            "Generate leads, then hand permanent candidates to exact structural proof.",
        ),
        Technique(
            "Tropical score-polytope certificate",
            "tropical",
            "nonrealizability",
            "Piecewise-linear score and slack data can expose impossible dominance regions.",
            ("slack margins", "normal fan", "score polytope face"),
            ("exact cycle witnesses",),
            ("normal-fan cells", "score regions", "active bottlenecks"),
            ("tropical", "slack", "normal", "fan", "score", "polytope", "bottleneck"),
            ("pascal_cap", "residual_capacitor", "packet_normalizer"),
            "Pull active-bottleneck normal-fan data into tournament score inequalities.",
        ),
        Technique(
            "Automaton-state tournament",
            "automata",
            "proof_routing",
            "Orient automatic states by route-purity and sidecar cost, then inspect SCCs and cycles.",
            ("automatic fiber purity", "state transition debt", "route-sidecar cost"),
            ("exact magnitude unless sidecar attached",),
            ("automaton states", "fiber classes", "repair teeth"),
            ("automaton", "state", "fiber", "route", "transition", "purity"),
            ("packet_normalizer", "observer_gluing", "level7_state_lift"),
            "Use after the HYP-2963 packet bank emits exact residual state words.",
        ),
        Technique(
            "Matrix observability tournament",
            "observability",
            "proof_routing",
            "Columns of a sidecar matrix become vertices, oriented by which separates more live fibers.",
            ("observable column rank", "fiber separation", "hidden coordinate audit"),
            ("nonlinear interactions unless columns are paired",),
            ("sidecar columns", "payload matrices", "quotient fibers"),
            ("matrix", "observability", "column", "fiber", "separation", "hidden"),
            ("matrix_observability", "packet_normalizer", "observer_gluing"),
            "Treat tournament certificates as columns in the existing controlled-forgetting matrix.",
        ),
    ]


def make_targets() -> list[Target]:
    return [
        Target(
            "LRC14 observer-gluing residual",
            "observer_gluing",
            "lrc14",
            "proof-frontier coverage",
            ("observer charts", "proof obligations", "terminal exits"),
            ("predicate retention", "chart overlap", "normalized arc", "CRT lift", "sidecar"),
            ("endpoint owner", "active binder", "cap defect", "branch K33", "level7"),
            "Current S258/S259 frontier: direct witness chart and Pascal chart must glue.",
        ),
        Target(
            "Coarse mod-14 winding degeneracy",
            "coarse_winding",
            "lrc14",
            "guardrail",
            ("winding tournament", "tied antipodal pairs", "coarse residues"),
            ("tie repair", "fine scale", "nondegenerate comparison"),
            ("mod7 lift", "dyadic lift", "trienerment"),
            "Coarse H is not theorem currency once antipodal ties appear.",
        ),
        Target(
            "Fine mod-7 winding and Paley comparison",
            "fine_mod7_winding",
            "lrc14",
            "candidate obstruction",
            ("fine residue tournament", "prime field chart", "Paley T7 comparison"),
            ("regularity", "H max", "cycle census", "skew spectrum"),
            ("CRT level7 status", "residue owners"),
            "The c=7 lift is the live field-like chart after 14=2*7 breaks the prime-field shortcut.",
        ),
        Target(
            "H=21 permanent-gap closure",
            "h21_closure",
            "tournament",
            "structural contradiction",
            ("Omega component", "alpha-vector", "H value", "SCC factors"),
            ("H=21 forbidden", "Omega realizability", "cycle spectrum"),
            ("Moon pancyclicity", "component product", "small graph candidates"),
            "Continue the H=7 method by ruling out all candidate I(Omega,2)=21 components.",
        ),
        Target(
            "Level-7 state-lift residual",
            "level7_state_lift",
            "lrc14",
            "state-lift contradiction",
            ("F7 residual sectors", "K33 lift", "CRT c=7 descent"),
            ("H=7 or H=21 terminal", "SCC factor", "mod7 sidecar"),
            ("c=2 lift", "covering residual count", "endpoint owner"),
            "If a terminal state-lift emits a complete tournament with forbidden H, the residual dies.",
        ),
        Target(
            "Pascal cap and symbolic overlap ledger",
            "pascal_cap",
            "lrc14",
            "moment/cap proof route",
            ("pair mass", "inclusion-exclusion overlap", "cap defect"),
            ("score majorization", "alpha-vector", "normal fan"),
            ("S3/S4/Perron", "endpoint sectors", "symbolic overlap"),
            "Pair mass is exact through j<=3; j=4,5 need higher-order deviation closure.",
        ),
        Target(
            "HYP-2963 packet normalizer",
            "packet_normalizer",
            "lrc14",
            "global proof interface",
            ("finite-address packets", "sidecar matrix", "route exits"),
            ("fiber purity", "bridge protection", "median center", "observability"),
            ("exact M", "endpoint topology", "Haar", "Fejer", "K33", "q-cusp"),
            "The final proof spine needs route-pure packet emission or named residual debt.",
        ),
        Target(
            "A000568 first-defect observer cuts",
            "a000568_defect",
            "tournament",
            "controlled-forgetting guardrail",
            ("rooted perspectives", "observer cuts", "incident words"),
            ("automorphism parity", "metagraph embedding", "Omega shape"),
            ("edge-sector deck", "rectangle residues", "deletion parent"),
            "First class-count failures are missing observer payload, not raw enumeration noise.",
        ),
        Target(
            "Residual capacitor and bridge graph",
            "residual_capacitor",
            "lrc14",
            "repair/obstruction",
            ("residual branch graph", "capacitor cuts", "endpoint owner strips"),
            ("strong orientation", "no naked bridge", "normal fan"),
            ("Haar zeta", "primitive deck", "route certificate"),
            "Bridge contraction is legal only after protected reverse paths or terminal exits are named.",
        ),
        Target(
            "Sexy prime fixed-gap sidecar",
            "sexy_prime_sidecar",
            "number theory",
            "sieve bookkeeping",
            ("fixed gap ray", "midpoint residues", "sieve channels"),
            ("residue admissibility", "parity debt", "distribution debt"),
            ("prime-power sidecars", "Hardy-Littlewood weight", "terminal exits"),
            "Tournaments help audit local residue channels but do not break sieve parity.",
        ),
        Target(
            "Loose digraph H=7 guardrail",
            "loose_digraph_guardrail",
            "tournament",
            "analogy failure",
            ("loose directed graph", "incomplete comparison shadow", "Omega-like graph"),
            ("complete tournament condition", "H forbidden only after valid pullback"),
            ("missing edges", "tie path", "encoding functor"),
            "Arbitrary digraph H=7 examples do not refute THM-200; completeness is load-bearing.",
        ),
        Target(
            "Matrix observability sidecar",
            "matrix_observability",
            "lrc14",
            "quotient audit",
            ("payload matrix", "sidecar columns", "fiber-separation rows"),
            ("transfer symmetry", "skew spectrum", "column rank"),
            ("observer cut orbit", "hidden coordinate", "dual annihilator"),
            "A sidecar matrix can say which tournament certificate column actually separates the live fiber.",
        ),
    ]


def technique_score(tech: Technique, target: Target) -> int:
    tech_words = tokens(
        (
            tech.name,
            tech.family,
            tech.mode,
            tech.certificate,
            tech.next_experiment,
            *tech.preserves,
            *tech.vertices,
            *tech.trigger_keywords,
            *tech.applications,
        )
    )
    target_words = tokens(
        (
            target.name,
            target.key,
            target.domain,
            target.claim_type,
            target.notes,
            *target.encoded_as,
            *target.desired_signal,
            *target.sidecars,
        )
    )
    overlap = tech_words & target_words
    score = len(overlap)
    score += 4 if target.key in tech.applications else 0
    score += 2 if tech.family in target_words else 0
    if tech.mode in {"contradiction", "nonrealizability"} and target.claim_type in {
        "candidate obstruction",
        "state-lift contradiction",
        "structural contradiction",
    }:
        score += 3
    if tech.mode in {"repair", "proof_routing"} and target.claim_type in {
        "proof-frontier coverage",
        "global proof interface",
        "guardrail",
    }:
        score += 3
    if "parity" in tech_words and "parity" in target_words:
        score += 3
    return score


def score_table(techniques: list[Technique], targets: list[Target]) -> list[tuple[int, str, str]]:
    rows = []
    for tech in techniques:
        for target in targets:
            rows.append((technique_score(tech, target), tech.name, target.name))
    return sorted(rows, key=lambda row: (-row[0], row[1], row[2]))


def landau_samples() -> list[tuple[str, tuple[int, ...], bool]]:
    samples = [
        ("regular_13_score", (6,) * 13, True),
        ("Paley_7_score", (3,) * 7, True),
        ("impossible_low_prefix_7", (0, 0, 0, 6, 6, 6, 6), False),
        ("too_many_sources_6", (0, 0, 0, 5, 5, 5), False),
        ("transitive_6", tuple(range(6)), True),
    ]
    return [(name, seq, kps.is_landau(seq)) for name, seq, _ in samples]


def enumerate_small_spectra(max_n: int = 5) -> dict[str, object]:
    cumulative_h: set[int] = set()
    by_n: dict[int, dict[str, object]] = {}
    alpha_examples: Counter[tuple[int, ...]] = Counter()
    cycle_census: Counter[tuple[int, tuple[tuple[int, int], ...]]] = Counter()
    even_ham = 0
    total = 0
    for n in range(1, max_n + 1):
        h_values: set[int] = set()
        hams: Counter[int] = Counter()
        scores_seen: Counter[tuple[int, ...]] = Counter()
        for adj in kps.all_tournaments(n):
            total += 1
            h_value = kps.H_value(adj)
            h_values.add(h_value)
            cumulative_h.add(h_value)
            ham = kps.ham_path_count(adj)
            hams[ham] += 1
            if ham % 2 == 0:
                even_ham += 1
            alpha_examples[tuple(kps.alpha_vector(adj))] += 1
            cycle_census[(n, tuple(sorted(kps.odd_cycle_counts(adj).items())))] += 1
            scores_seen[tuple(kps.scores(adj))] += 1
        by_n[n] = {
            "count": 2 ** (n * (n - 1) // 2),
            "H_values": sorted(h_values),
            "ham_path_values": sorted(hams),
            "score_sequences": len(scores_seen),
        }
    max_h = max(cumulative_h)
    odd_gaps = [v for v in range(1, max_h + 1, 2) if v not in cumulative_h]
    return {
        "by_n": by_n,
        "cumulative_H": sorted(cumulative_h),
        "odd_gaps_to_max": odd_gaps,
        "alpha_common": alpha_examples.most_common(8),
        "cycle_census_common": cycle_census.most_common(8),
        "even_hamiltonian_path_count_hits": even_ham,
        "total_tournaments": total,
    }


def certificate_probes() -> list[tuple[int, str]]:
    values = [0, 1, 2, 3, 5, 7, 8, 9, 14, 21, 35, 63, 189]
    out = []
    for value in values:
        reason = kps.H_spectrum_certificate(value)
        out.append((value, reason or "not rejected by current permanent-gap table"))
    return out


def dominance_tournament(
    techniques: list[Technique], targets: list[Target], selected_names: list[str]
) -> tuple[list[Technique], list[list[int]]]:
    selected = [tech for tech in techniques if tech.name in selected_names]
    n = len(selected)
    adj = [[0] * n for _ in range(n)]
    scores_by_target = {
        (tech.name, target.key): technique_score(tech, target)
        for tech in selected
        for target in targets
    }
    for i, j in itertools.combinations(range(n), 2):
        a = selected[i]
        b = selected[j]
        a_wins = 0
        b_wins = 0
        for target in targets:
            sa = scores_by_target[(a.name, target.key)]
            sb = scores_by_target[(b.name, target.key)]
            if sa > sb:
                a_wins += 1
            elif sb > sa:
                b_wins += 1
        if a_wins > b_wins:
            adj[i][j] = 1
        elif b_wins > a_wins:
            adj[j][i] = 1
        else:
            # Tie gauge: fewer destroyed coordinates, then name order.
            a_key = (len(a.destroys), a.name)
            b_key = (len(b.destroys), b.name)
            if a_key <= b_key:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
    return selected, adj


def redei_path(adj: list[list[int]]) -> list[int]:
    path: list[int] = []
    for v in range(len(adj)):
        inserted = False
        for idx, w in enumerate(path):
            if adj[v][w]:
                path.insert(idx, v)
                inserted = True
                break
        if not inserted:
            path.append(v)
    return path


def tournament_fingerprint(adj: list[list[int]]) -> dict[str, object]:
    score_hist = Counter(sum(row) for row in adj)
    triples = 0
    for a, b, c in itertools.combinations(range(len(adj)), 3):
        edges = adj[a][b] + adj[b][c] + adj[c][a]
        if edges in {0, 3}:
            triples += 1
    return {
        "score_histogram": dict(sorted(score_hist.items())),
        "directed_3_cycles": triples,
        "SCCs": [sorted(comp) for comp in kps.sccs(adj)],
        "Hamiltonian_path_count": kps.ham_path_count(adj),
        "Redei_path": redei_path(adj),
    }


def application_ledger() -> list[tuple[str, str, str]]:
    return [
        (
            "LRC14 observer gluing",
            "Use proof-obligation dominance, functorial pullback, and survival profiles before any H contradiction.",
            "A forbidden H hit is terminal only after the chart overlap map proves the tournament pullback preserves the LRC predicate.",
        ),
        (
            "Coarse mod-14 winding",
            "Use trienerment/fine-scale repair rather than raw H.",
            "Antipodal ties make the coarse tournament a guardrail; move to c=7/c=2 lift sidecars.",
        ),
        (
            "Fine mod-7 residual",
            "Compare against Paley/regular H-max, score, cycle census, and skew spectrum.",
            "If the residual is forced into Paley-like regularity, exploit rigidity; if not, use the failed regularity as sidecar data.",
        ),
        (
            "H=21 closure",
            "Combine H forbidden-value, Omega-shape mining, alpha-vector, and SCC product descent.",
            "The next exact job is not another scalar search but an Omega-component realizability exclusion, with K_10 the clean single-component package and P4 already excluded.",
        ),
        (
            "K33/F7 state lift",
            "Run H=7/H=21 only on complete tournament certificates, then localize by SCC factors.",
            "Loose digraph shadows may have H=7; completeness and pullback validity are load-bearing.",
        ),
        (
            "Apex-7 bridge",
            "Retire the slogan apex-7 equals forbidden-H=7; keep the order-2 antipodal / odd-automorphism lever.",
            "The apex antipodal matching is triangle-free, structurally opposite to Omega=K3, so the H=7 match is a coincidence.",
        ),
        (
            "Pascal cap route",
            "Use score/tropical/polytope and alpha-vector gates as higher-order deviation auditors.",
            "Pair mass is a legal equinumerous shadow through j<=3 but needs S3/S4/Perron or symbolic overlap debt at j=4,5.",
        ),
        (
            "HYP-2963 normalizer",
            "Treat tournament certificates as sidecar columns and rank them by fiber separation.",
            "A no-hit profile still records which sidecars are necessary for global packet emission and for the Lean LRCBleedingEdgeFrontier wrapper.",
        ),
        (
            "A000568 defect",
            "Use odd automorphism, metagraph nonembedding, and observer-cut payloads.",
            "The first defect is missing rooted/incident observer data, not raw class-count noise.",
        ),
        (
            "Residual capacitors",
            "Apply no-naked-bridge orientation before contracting proof branch graphs.",
            "Any tournament-derived contradiction must survive protected reverse verification or name bridge debt.",
        ),
        (
            "Sexy prime conjecture",
            "Use residue-sieve tournaments for admissibility ledgers and local-channel audits.",
            "This does not break parity or prove fixed-gap infinitude; it clarifies what a legal sieve sidecar would need to retain.",
        ),
    ]


def print_section(title: str) -> None:
    print()
    print(f"=== {title} ===")


def main() -> None:
    techniques = make_techniques()
    targets = make_targets()
    rows = score_table(techniques, targets)
    family_counts = Counter(tech.family for tech in techniques)
    target_best: dict[str, list[tuple[int, str]]] = defaultdict(list)
    for score, tech, target in rows:
        target_best[target].append((score, tech))

    print("Tournament contradiction grammar -- codex-2026-06-27-S260")
    print("Integrated source: tournament_certificate_engine_kps.py / S31ah")
    print("Rebased incoming: HYP-3099/S65 applications and S31ah single-component H ladder")
    print("Formal sink: TournamentH7.LRCBleedingEdgeFrontier packages the observer/pascal/witness ledgers.")
    print("Incoming constraints: H=7 <-> Omega K_3 gap; H=21 <-> Omega K_10 gap;")
    print("  cap exchange is non-transitive but bounded; baby-Hodge holes are c5/spectral;")
    print("  apex-7 equals H=7 is a coincidence, order-2 antipodal symmetry is the usable lever.")
    print(f"Technique count: {len(techniques)}")
    print(f"Target count: {len(targets)}")
    print(f"Technique families: {dict(sorted(family_counts.items()))}")

    print_section("H-spectrum certificate probes")
    for value, reason in certificate_probes():
        print(f"H={value}: {reason}")

    print_section("Landau sanity probes")
    for name, seq, result in landau_samples():
        print(f"{name}: sequence={seq} landau={result}")

    print_section("Small exact tournament spectra")
    spectra = enumerate_small_spectra(5)
    print(f"Enumerated tournaments: {spectra['total_tournaments']}")
    print(f"Even Hamiltonian path count violations: {spectra['even_hamiltonian_path_count_hits']}")
    for n, data in spectra["by_n"].items():
        print(
            f"n={n}: count={data['count']} H_values={data['H_values']} "
            f"ham_path_values={data['ham_path_values']} score_sequences={data['score_sequences']}"
        )
    print(f"Cumulative H values through n=5: {spectra['cumulative_H']}")
    print(f"Odd H gaps through max exact value: {spectra['odd_gaps_to_max'][:40]}")
    print("Common alpha-vectors:")
    for alpha, count in spectra["alpha_common"]:
        print(f"  {alpha}: {count}")
    print("Common cycle census entries:")
    for census, count in spectra["cycle_census_common"]:
        print(f"  {census}: {count}")

    print_section("Top scored technique-target applications")
    for score, tech, target in rows[:36]:
        print(f"{score:02d} | {tech} -> {target}")

    print_section("Best moves by current target")
    for target in sorted(target_best):
        best = target_best[target][:4]
        best_text = "; ".join(f"{score}:{tech}" for score, tech in best)
        print(f"{target}: {best_text}")

    selected_names = []
    for _, tech, _ in rows:
        if tech not in selected_names:
            selected_names.append(tech)
        if len(selected_names) == 8:
            break
    selected, adj = dominance_tournament(techniques, targets, selected_names)
    fingerprint = tournament_fingerprint(adj)

    print_section("Tournament Analysis of the generated technique frontier")
    print("Vertices: certificate families / proof obligations, not runners.")
    print("Pairwise observable: number of target ledgers on which one technique scores higher.")
    print("Switch/gauge: orient A->B if A wins more target ledgers; tie by fewer destroyed coordinates.")
    print("Selected vertices:")
    for idx, tech in enumerate(selected):
        print(f"  {idx}: {tech.name}")
    print(f"Score histogram: {fingerprint['score_histogram']}")
    print(f"Directed 3-cycles: {fingerprint['directed_3_cycles']}")
    print(f"SCCs: {fingerprint['SCCs']}")
    print(f"Hamiltonian path count: {fingerprint['Hamiltonian_path_count']}")
    path_names = [selected[idx].name for idx in fingerprint["Redei_path"]]
    print(f"Tie Hamiltonian path: {' > '.join(path_names)}")

    print_section("Application ledger")
    for target, move, readout in application_ledger():
        print(f"{target}:")
        print(f"  move: {move}")
        print(f"  readout: {readout}")

    print_section("Assumption challenge")
    print("Considered vertex sets: runners, arcs, residues, Omega components, score sequences,")
    print("proof obligations, observer charts, sidecar columns, quotient maps, finite sieve channels,")
    print("automaton states, rooted perspectives, and branch-graph edges.")
    print("Chosen default for this pass: certificate functors and target proof obligations.")
    print("Preserved predicate: whether a pulled-back tournament certificate can force a contradiction,")
    print("repair a degenerate encoding, or record a necessary sidecar for LRC14/sexy-prime ledgers.")
    print("Destroyed information: concrete runner geometry, exact metric scale, and distributional")
    print("number-theory strength unless the target sidecars explicitly retain them.")
    print("Challenged assumption: tournaments are not only H-value contradiction machines;")
    print("they can be score-polytopes, sidecar matrices, route medians, bridge guards,")
    print("residue-channel auditors, and no-hit necessary-condition miners.")


if __name__ == "__main__":
    main()
