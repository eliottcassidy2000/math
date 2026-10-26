#!/usr/bin/env python3
"""
Broad repo concept atlas for the LRC14 endgame.

This is a synthesis script, not a search proof.  It gathers a deliberately wide
cross-section of prior repo mechanisms and scores how strongly each one speaks
to the current LRC14 proof interface:

    q/Farey branch
    C=27 typed shell
    Kpq/K33 incidence
    visible-vs-hidden relation mass
    HYP-2908 state-lift packets
    finite/analytic split

The main output is a carrier tournament whose vertices are mechanism families,
not runners.  The point is to make the "almost random" archaeology reusable:
which old ideas should be imported into the bleeding-edge LRC14 route, and in
which order?
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
import re
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class Concept:
    name: str
    refs: tuple[str, ...]
    family: str
    summary: str
    lrc14_connection: str
    score: tuple[int, int, int, int, int, int]
    tags: tuple[str, ...]

    @property
    def total(self) -> int:
        return sum(self.score)


SCORE_FIELDS = (
    "branch_retention",
    "typed_visibility",
    "state_lift_fit",
    "finite_certifiability",
    "anti_scalar_guard",
    "cross_problem_signal",
)


CONCEPTS: tuple[Concept, ...] = (
    Concept(
        name="Farey/q binding scale",
        refs=("HYP-2930", "HYP-2931", "T1027"),
        family="current LRC14",
        summary="Exact value M=p/q has gap (14p-q)/(14q); q remains theorem scale.",
        lrc14_connection="Primary address for distinguishing tight, unit-excess, and nonunit branches.",
        score=(5, 3, 3, 5, 5, 3),
        tags=("farey", "q", "branch"),
    ),
    Concept(
        name="C=27 shell and Yoneda coimage",
        refs=("HYP-2161", "HYP-2934", "T1030"),
        family="coimage/resonance",
        summary="At n=14, C=2n-1=27 is the conservative shell probe candidate.",
        lrc14_connection="Controls the p=2 petal/two-block branch and unit/nonunit shell holes.",
        score=(5, 5, 4, 4, 5, 4),
        tags=("C27", "coimage", "yoneda", "shell"),
    ),
    Concept(
        name="Bigraded summand/multiplicand signature",
        refs=("HYP-2935", "T1031"),
        family="current LRC14",
        summary="Addition gives relation shells; multiplication tests observer visibility.",
        lrc14_connection="Separates visible folds, hidden balanced collisions, C27 shells, and Kpq rank after the Farey branch is known.",
        score=(5, 5, 4, 5, 5, 4),
        tags=("summand", "multiplicand", "visible", "hidden"),
    ),
    Concept(
        name="Kpq/K33 incidence wall",
        refs=("HYP-2932", "HYP-2934", "T1028"),
        family="graph/incidence",
        summary="p/q has an incidence blow-up K_{p,q}; p>=3 is the first K33-capable lane.",
        lrc14_connection="Routes p>=3 unit-excess rows to finite three-owner obstruction packets.",
        score=(4, 4, 5, 5, 4, 4),
        tags=("Kpq", "K33", "incidence"),
    ),
    Concept(
        name="Octahedron L(K4) support-six current",
        refs=("HYP-2887", "HYP-2933", "T1002"),
        family="packet/Hodge",
        summary="Repeated support-six packets form an octahedral current/curl carrier with cycle rank 7.",
        lrc14_connection="Candidate carrier for signed support-six relation mass before HYP-2908 state lift.",
        score=(3, 5, 5, 4, 4, 4),
        tags=("octahedron", "support-six", "hodge"),
    ),
    Concept(
        name="Clebsch/halved-cube covariance",
        refs=("HYP-2891", "HYP-2892", "HYP-2933"),
        family="packet/Hodge",
        summary="Folded residual masks preserve pair covariance/cut data while destroying scalar depth.",
        lrc14_connection="Good for organizing folded residual masks after low-depth AP/Fejer rows are removed.",
        score=(3, 4, 4, 4, 5, 4),
        tags=("Clebsch", "halfcube", "covariance"),
    ),
    Concept(
        name="Depth-corrected Bonferroni parity",
        refs=("HYP-2903", "T1017"),
        family="measure/overlap",
        summary="High Newton packets are signed by missing-depth parity, not just order.",
        lrc14_connection="Controls the multi-large branch where finite moments and naive Bonferroni-3 fail.",
        score=(3, 5, 4, 3, 5, 4),
        tags=("bonferroni", "helly", "depth"),
    ),
    Concept(
        name="Vitali wall / anti-Poisson coimage",
        refs=("S604", "THM-406", "HYP-2152"),
        family="measure/overlap",
        summary="p0 collapse is all-orders overlap cancellation; finite moments can be blind.",
        lrc14_connection="Explains why AP, GW, C27 petals, and hidden-fold rows must not be scalarized by energy alone.",
        score=(4, 5, 4, 3, 5, 5),
        tags=("vitali", "anti-poisson", "coimage"),
    ),
    Concept(
        name="Fejer/additive-energy AP-facing carrier",
        refs=("HYP-2889", "HYP-2890", "HYP-+2873"),
        family="Fourier/additive",
        summary="Additive energy is a spectral fourth moment and AP-facing Fejer shadow, not a monotone invariant.",
        lrc14_connection="Useful only after visible/hidden folds and branch labels are retained.",
        score=(2, 3, 2, 4, 5, 4),
        tags=("fejer", "energy", "Fourier"),
    ),
    Concept(
        name="Exact-period Mobius/totient packet ledger",
        refs=("HYP-2886", "HYP-2899", "HYP-2900"),
        family="number-theory packets",
        summary="Exact denominators live on divisor packets with phi capacity and CRT curvature.",
        lrc14_connection="Analytic branch after finite C27/K33 routing; prevents finite denominator basis fantasies.",
        score=(4, 4, 3, 4, 5, 4),
        tags=("totient", "Mobius", "CRT", "exact-period"),
    ),
    Concept(
        name="Boundary-state induction / one-large peeler",
        refs=("HYP-2904", "HYP-2905", "HYP-2906"),
        family="induction",
        summary="Large-speed induction needs safe-set component budget or one witness interval.",
        lrc14_connection="Shrinks the proof to top-balanced bounded/finite-core atoms before the carrier split.",
        score=(3, 3, 3, 5, 4, 3),
        tags=("induction", "boundary-state", "peeler"),
    ),
    Concept(
        name="Forbidden-H7 state lift",
        refs=("HYP-2908", "THM-572", "T1025"),
        family="tournament endpoint",
        summary="A bad LRC atom is dead if it constructs a tournament packet with H=7.",
        lrc14_connection="Terminal contradiction for three-owner/K33 or support-six packet branches.",
        score=(4, 4, 5, 5, 4, 4),
        tags=("H7", "state-lift", "OCF"),
    ),
    Concept(
        name="Tournament deck derivative repair",
        refs=("HYP-2236", "S660"),
        family="projection repair",
        summary="Deck collisions are fixed by pairing each deleted card with its boundary derivative.",
        lrc14_connection="Template for repairing LRC projections: card plus owner/visibility label, not global scalar add-ons.",
        score=(3, 4, 4, 4, 5, 5),
        tags=("deck", "derivative", "boundary"),
    ),
    Concept(
        name="Unit-distance Hamiltonian flop",
        refs=("HYP-2202", "T1026"),
        family="unit distance",
        summary="Geometric unit paths persist while canonical tiling-order HPs can flop at n=7.",
        lrc14_connection="Warns that a guaranteed tournament path may be an order artifact unless its side-channel support is labelled.",
        score=(2, 4, 3, 4, 5, 5),
        tags=("unit-distance", "hamiltonian", "flop"),
    ),
    Concept(
        name="Unit-distance impairment spectroscopy",
        refs=("HYP-2194", "S622"),
        family="unit distance",
        summary="Drop direction/gain side channels and watch exact small-size methods fail.",
        lrc14_connection="Suggests deliberately impairing C27 shell lanes or K33 owners to find load-bearing proof channels.",
        score=(2, 5, 3, 4, 5, 5),
        tags=("unit-distance", "impairment", "Helly"),
    ),
    Concept(
        name="Cauldron adversarial residue schedules",
        refs=("HYP-2190", "HYP-2191", "HYP-2192"),
        family="games/additive",
        summary="Base cauldrons are Schur-frontier hard; adversarial variants add correlated residue streams.",
        lrc14_connection="Concrete toy model for density-blind two-block/correlated-residue failures.",
        score=(2, 4, 2, 4, 5, 5),
        tags=("cauldron", "Schur", "residue"),
    ),
    Concept(
        name="Summand graph / Zeckendorf / forbidden 7,21",
        refs=("summand-graph-fermat-zeckendorf", "THM-079", "THM-200"),
        family="additive/tournament",
        summary="Summand graph phase transitions encode triangular/tournament thresholds and forbidden H values.",
        lrc14_connection="Backbone for reading C27 and p+q as additive shell data rather than raw magnitude.",
        score=(3, 4, 3, 3, 4, 5),
        tags=("summand", "Zeckendorf", "forbidden-H"),
    ),
    Concept(
        name="OCF as hard-core partition function",
        refs=("THM-002", "HYP-1800", "support-residue-calculus"),
        family="tournament endpoint",
        summary="H(T)=I(Omega(T),2); forbidden values are conflict-packet realizability gaps.",
        lrc14_connection="Gives the packet category for HYP-2908; scalar H is only the evaluation after packet realization.",
        score=(3, 4, 5, 4, 5, 5),
        tags=("OCF", "partition-function", "noise"),
    ),
    Concept(
        name="Collatz/two-block correlated residue",
        refs=("lrc-collatz-famous-problems-not-pvsnp", "collatz-iterated-log-tower"),
        family="residue analogy",
        summary="Density/drift heuristics are blind to exact correlated residue cycles.",
        lrc14_connection="Analogy for why LRC residuals need labelled congruence packets, not averaged density alone.",
        score=(2, 3, 2, 3, 5, 5),
        tags=("Collatz", "two-block", "residue"),
    ),
    Concept(
        name="p-curvature/support gate",
        refs=("p_curvature_lrc14_support_gate", "grothendieck-katz-p-curvature-and-lrc14-ledgers"),
        family="arithmetic gate",
        summary="Local prime gates can certify support structure without carrying full analytic content.",
        lrc14_connection="Potential exact-period/local-prime filter for residual branches after branch labels are known.",
        score=(3, 3, 3, 3, 4, 4),
        tags=("p-curvature", "support", "prime"),
    ),
    Concept(
        name="Pollock/Sierpinski carry-pair lift",
        refs=("HYP-2491", "HYP-2497"),
        family="additive basis",
        summary="Single residues were too weak; lifted defect pairs revealed the dyadic carry structure.",
        lrc14_connection="Another projection-repair template: lift from scalar residue to paired carry/owner defects.",
        score=(2, 3, 3, 4, 5, 5),
        tags=("Pollock", "carry", "defect-pair"),
    ),
)


THEME_PATTERNS = {
    "farey": r"Farey|mediant|p/q|p\+q|p\*q",
    "C27": r"C=27|2n-1|Res_27|27=3\^3|27=3",
    "K33": r"K33|K_\{3,3\}|K_3,3|K\{3,3\}|K_\{p,q\}|Kpq",
    "visible_hidden": r"visible|hidden|fold|balanced collision|pair-sum",
    "coimage": r"coimage|Yoneda|anti-Poisson|Vitali|p_0",
    "Helly": r"Helly|Bonferroni|depth",
    "OCF": r"OCF|I\(Omega|I\(\\Omega|Hamiltonian path|forbidden H|H=7|H=21",
    "unit_distance": r"unit.distance|Moser|Hamiltonian.*unit|impairment",
    "cauldron": r"cauldron|Schur|last-boil|adversarial",
    "exact_period": r"totient|Mobius|CRT|exact-period|denominator wall",
    "octa_clebsch": r"octahedron|Clebsch|halved cube|support-six|Hodge",
}


SCAN_FILES = (
    "00-navigation/TANGENTS.md",
    "00-navigation/CONCEPT-MAP.md",
    "05-knowledge/hypotheses/INDEX.md",
    "07-reflections/anti-poisson-coimage-atlas-s604.md",
    "07-reflections/lrc-coimage-yoneda-2nm1-cancellation-s605.md",
    "07-reflections/summand-graph-fermat-zeckendorf.md",
    "05-knowledge/hypotheses/HYP-2935-lrc14-bigraded-relation-signature.md",
    "05-knowledge/hypotheses/HYP-2934-lrc14-summand-multiplicand-farey-bridge.md",
    "05-knowledge/hypotheses/HYP-2933-lrc14-farey-graph-pz-carrier-synthesis.md",
)


def read_file(rel: str) -> str:
    path = ROOT / rel
    try:
        return path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        return path.read_text(encoding="utf-8", errors="replace")


def theme_counts() -> dict[str, int]:
    counts: dict[str, int] = {}
    corpus = "\n".join(read_file(rel) for rel in SCAN_FILES if (ROOT / rel).exists())
    for theme, pattern in THEME_PATTERNS.items():
        counts[theme] = len(re.findall(pattern, corpus, flags=re.IGNORECASE))
    return dict(sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])))


def family_summary(concepts: Iterable[Concept]) -> list[tuple[str, int, int, tuple[str, ...]]]:
    buckets: dict[str, list[Concept]] = defaultdict(list)
    for concept in concepts:
        buckets[concept.family].append(concept)
    rows = []
    for family, vals in buckets.items():
        rows.append(
            (
                family,
                len(vals),
                sum(v.total for v in vals),
                tuple(v.name for v in sorted(vals, key=lambda c: (-c.total, c.name))[:3]),
            )
        )
    return sorted(rows, key=lambda row: (-row[2], row[0]))


def role_tournament(concepts: tuple[Concept, ...]) -> tuple[dict[int, set[int]], Counter[int], int, tuple[int, ...], int]:
    edges: dict[int, set[int]] = {i: set() for i in range(len(concepts))}
    for i, a in enumerate(concepts):
        for j, b in enumerate(concepts):
            if i == j:
                continue
            key_a = (a.score, a.total, -len(a.name), a.name)
            key_b = (b.score, b.total, -len(b.name), b.name)
            if key_a > key_b:
                edges[i].add(j)

    score_hist = Counter(len(edges[i]) for i in edges)
    c3 = 0
    n = len(concepts)
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                if j in edges[i] and k in edges[j] and i in edges[k]:
                    c3 += 1
                if k in edges[i] and j in edges[k] and i in edges[j]:
                    c3 += 1
    hp = tuple(sorted(range(n), key=lambda i: (concepts[i].score, concepts[i].total, concepts[i].name), reverse=True))
    flips_vs_maturity = 0
    maturity_order = {name: rank for rank, name in enumerate([
        "OCF as hard-core partition function",
        "Summand graph / Zeckendorf / forbidden 7,21",
        "Fejer/additive-energy AP-facing carrier",
        "Boundary-state induction / one-large peeler",
        "Forbidden-H7 state lift",
        "C=27 shell and Yoneda coimage",
        "Farey/q binding scale",
        "Kpq/K33 incidence wall",
        "Bigraded summand/multiplicand signature",
    ])}
    for i in range(n):
        for j in range(i + 1, n):
            a_old = maturity_order.get(concepts[i].name, 999)
            b_old = maturity_order.get(concepts[j].name, 999)
            old_pref_i = a_old < b_old or (a_old == b_old and concepts[i].name < concepts[j].name)
            role_pref_i = j in edges[i]
            if old_pref_i != role_pref_i:
                flips_vs_maturity += 1
    return edges, score_hist, c3, hp, flips_vs_maturity


def pairwise_synergies(concepts: tuple[Concept, ...]) -> list[tuple[int, str, str, tuple[str, ...]]]:
    rows = []
    for i, a in enumerate(concepts):
        for j in range(i + 1, len(concepts)):
            b = concepts[j]
            shared = tuple(sorted(set(a.tags).intersection(b.tags)))
            branch_bridge = int(("C27" in a.tags or "C27" in b.tags) and ("K33" in a.tags or "K33" in b.tags))
            score = (
                len(shared) * 3
                + min(a.score[0] + b.score[0], 9)
                + min(a.score[1] + b.score[1], 9)
                + min(a.score[2] + b.score[2], 9)
                + branch_bridge * 7
            )
            if score >= 24:
                rows.append((score, a.name, b.name, shared))
    return sorted(rows, reverse=True)[:12]


def main() -> None:
    print("S135 BROAD CONCEPT ATLAS FOR LRC14")
    print("=" * 78)
    print("Scored dimensions:", ", ".join(SCORE_FIELDS))
    print()

    print("[1] Repo motif counts from selected archaeology files")
    for theme, count in theme_counts().items():
        print(f"  {theme:16s} {count:5d}")
    print()

    print("[2] Top mechanism carriers for the current LRC14 proof interface")
    for idx, concept in enumerate(sorted(CONCEPTS, key=lambda c: (c.score, c.total, c.name), reverse=True), 1):
        if idx > 14:
            break
        score_str = "/".join(str(x) for x in concept.score)
        refs = ",".join(concept.refs[:3])
        print(f"  {idx:2d}. {concept.name:44s} total={concept.total:2d} score={score_str:11s} refs={refs}")
        print(f"      {concept.lrc14_connection}")
    print()

    print("[3] Family-level readout")
    for family, count, total, names in family_summary(CONCEPTS):
        print(f"  {family:22s} concepts={count:2d} total={total:3d} leaders={'; '.join(names)}")
    print()

    print("[4] High-value cross-thread synergies")
    for score, a, b, shared in pairwise_synergies(CONCEPTS):
        shared_str = ",".join(shared) if shared else "-"
        print(f"  score={score:2d} :: {a}  <->  {b}  shared={shared_str}")
    print()

    print("[5] Tournament Analysis on mechanism families")
    edges, score_hist, c3, hp, flips = role_tournament(CONCEPTS)
    print("  vertices considered/challenged:")
    print("    runners, arcs, speed rows, exact fractions, shell probes, cover arcs,")
    print("    relation packets, graph carriers, finite moments, OCF packets, and proof obligations.")
    print("  chosen vertices: repo mechanism families / proof carriers.")
    print("  pair observable:")
    print("    lexicographic role vector (branch retention, typed visibility, state-lift fit,")
    print("    finite certifiability, anti-scalar guard, cross-problem signal).")
    print("  switch/gauge: larger role vector wins; ties use total score and stable name order.")
    print(f"  fingerprint: score_hist={dict(sorted(score_hist.items()))} c3={c3} scc={(1,) * len(CONCEPTS)} hp=1")
    print("  first Hamiltonian path:")
    print("    " + " > ".join(CONCEPTS[i].name for i in hp[:10]))
    print(f"  old-maturity-only gauge would flip {flips} carrier pairs.")
    print()

    print("[6] Proof-facing synthesis")
    print("  The broad repo view reinforces one rule: keep labels until they have done")
    print("  their proof job.  The LRC14 endgame should not ask whether a row has high")
    print("  additive energy, many folds, a nice tournament, or a large product in the")
    print("  abstract.  It should ask which conservative packet is still live after")
    print("  the Farey/q branch is known.")
    print()
    print("  Proposed route order:")
    print("    (1) q/Farey branch and excess e=14p-q;")
    print("    (2) C=27 unit/nonunit shell probe for p=2 rows;")
    print("    (3) Kpq/K33 owner-incidence probe for p>=3 rows;")
    print("    (4) visible-vs-hidden summand/multiplicand relation signature;")
    print("    (5) octahedral/Clebsch/Hodge packet if support-six mass remains;")
    print("    (6) HYP-2908 tournament-state lift or exact-period analytic branch.")
    print()
    print("  Novel next lemma shape:")
    print("    If a top-balanced LRC14 residual survives the current peeler and is not")
    print("    AP/GW, then its conservative probe atlas contains either a C27 p=2")
    print("    shell defect with finite petal/lift discharge, or a sign-visible")
    print("    three-owner incidence packet whose connected OCF packet has activity-two")
    print("    value 7.  The scalar estimates enter only after this typed packet address.")


if __name__ == "__main__":
    main()
