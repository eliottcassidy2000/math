#!/usr/bin/env python3
"""LRC14 cocycle-sheaf unification pass.

This is a proof-interface synthesis, not a proof of LRC14.

The previous Haar-product work isolated the local 2-by-2 fixed-margin
cocycle

    zeta(T) = T00 - T01 - T10 + T11.

This pass treats the rest of the LRC14 stack in the same language.  A proof
route is useful when it names a cocycle, gives the base/fiber on which it
lives, and declares how the cocycle becomes exact, is annihilated by a dual
certificate, restricts to a smaller family, or remains as named residual
cohomology.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class Carrier:
    name: str
    kind: str
    cocycle: str
    base: str
    differential: str
    exact_or_annihilated: str
    residual: str
    destroys: str
    anchors: tuple[str, ...]
    vector: tuple[int, ...]


DIMENSIONS = (
    "local_coordinate",
    "packet_base",
    "boundary_operator",
    "annihilator",
    "restriction_handoff",
    "endpoint_owner",
    "exact_period",
    "formal_check",
    "named_residual",
)


CARRIERS: tuple[Carrier, ...] = (
    Carrier(
        name="haar_zipper_zeta",
        kind="local mixed Haar cocycle",
        cocycle="zeta(T)=T00-T01-T10+T11 on a 2x2 fixed-margin table",
        base="row/column margin fiber, then labelled LRC packet fiber",
        differential="fixed-margin switch changes zeta by 4 while margins stay fixed",
        exact_or_annihilated="color-compatible discrepancy pairs opposite zeta signs or orthogonality kills disjoint rectangles",
        residual="unpaired zeta is F7 / THM-572 state-lift debt",
        destroys="raw row/column margins and raw component count K",
        anchors=("HYP-2991", "HYP-2989", "HYP-2595"),
        vector=(5, 5, 5, 5, 5, 4, 2, 4, 5),
    ),
    Carrier(
        name="certificate_handoff_debt",
        kind="global proof-arrow cocycle",
        cocycle="debt cochain on q-witness, Fejer, Ramanujan, endpoint, moment, smoothing, and state-lift exits",
        base="HYP-2963 labelled packet families",
        differential="failed tooth emits a handoff obligation to the next proof clock",
        exact_or_annihilated="all six atlas arrows O1-O6 close or route to a named smaller packet",
        residual="SOURCE-SPECTRUM-UNKNOWN / F7",
        destroys="which proof method originally found the certificate",
        anchors=("HYP-2987", "HYP-2990", "THM-572"),
        vector=(4, 5, 5, 4, 5, 5, 4, 5, 5),
    ),
    Carrier(
        name="exposure_kernel_class",
        kind="unexposed-source cocycle",
        cocycle="row survives every exposure channel and therefore represents a kernel class",
        base="certificate-channel poset over AP-neighborhood and hard-frontier rows",
        differential="each exposure edge kills one obstruction coordinate",
        exact_or_annihilated="Fejer PSD dual, positive Haar mass, q-witness, C27/K33/petal label, or moment handoff exposes the class",
        residual="UNEXPOSED_SOURCE_KERNEL",
        destroys="fine geometry of the exposing certificate after exposure is recorded",
        anchors=("HYP-2988", "HYP-2974", "HYP-2963"),
        vector=(4, 5, 4, 5, 4, 4, 4, 5, 5),
    ),
    Carrier(
        name="tope_cocircuit_boundary_current",
        kind="oriented-matroid cocycle",
        cocycle="endpoint-owner boundary current on danger-arrangement walls",
        base="cyclic tope/cocircuit arrangement cut by (14k +/- 1)/(14v)",
        differential="crossing a wall changes the open tope and records an owner-labelled cocircuit",
        exact_or_annihilated="open all-safe tope is a witness; AP/GW zero-sum endpoints are boundary atoms",
        residual="no-tope/no-cocircuit forbidden wall packet",
        destroys="continuous runner identity inside a fixed endpoint arrangement",
        anchors=("HYP-2986", "HYP-2975", "HYP-2970"),
        vector=(4, 5, 5, 3, 4, 5, 2, 4, 4),
    ),
    Carrier(
        name="ramanujan_exact_period_cocycle",
        kind="character/exact-period cocycle",
        cocycle="primitive phase kernel c_q(a-b) on Z/qZ packet functions",
        base="exact denominator and unit-orbit packet fibers",
        differential="restriction from full phase data to exact-period primitive components",
        exact_or_annihilated="Ramanujan orthogonality annihilates off-period leakage when endpoint/state labels are retained",
        residual="prime-power or repeated-prime packet missed by squarefree quotients",
        destroys="endpoint ownership and open-vs-boundary status if used alone",
        anchors=("HYP-2979", "HYP-2978", "HYP-2982"),
        vector=(3, 5, 4, 5, 4, 3, 5, 4, 4),
    ),
    Carrier(
        name="fejer_toeplitz_dual_cochain",
        kind="dual-certificate cochain",
        cocycle="danger multiplicity F_S(t)=C_S(t)-1 paired with Fejer/Toeplitz test functions",
        base="safe interval packet with exact center, degree, atom bank, and route label",
        differential="negative trigonometric-square pairing is a Farkas/PSD coboundary certificate",
        exact_or_annihilated="formal interval Fejer certificate proves positive safe mass",
        residual="PSD-through-degree packet: AP/GW equality or state-lift-owned sector",
        destroys="combinatorial owner labels unless the packet manifest keeps them",
        anchors=("HYP-2974", "HYP-2981", "HYP-2988"),
        vector=(3, 4, 4, 5, 4, 3, 3, 5, 4),
    ),
    Carrier(
        name="smoothing_boundary_cocycle",
        kind="analytic deformation cocycle",
        cocycle="exceptional boundary/approach class created by changing kernels or smoothing functions",
        base="labelled packet plus kernel support radius and Kaczynski approach class",
        differential="kernel homotopy boundary defect or smoothing dispatcher arrow",
        exact_or_annihilated="open support radius, Abel decay, explicit formula, or large-sieve precondition discharges the defect",
        residual="true-wide resonant packet with named approach/state-lift debt",
        destroys="prime-power and endpoint data if reduced to squarefree scalar weights",
        anchors=("HYP-2985", "HYP-2984", "HYP-2983", "HYP-2982"),
        vector=(3, 4, 5, 4, 5, 4, 4, 4, 4),
    ),
    Carrier(
        name="pair_tension_coboundary",
        kind="literal runner-pair coboundary",
        cocycle="w_ij=v_j-v_i with w_ij+w_jk+w_ki=0 on every triple",
        base="complete graph on runners or pair-cell tournament",
        differential="delta v is a 1-coboundary, hence its triangle coboundary vanishes",
        exact_or_annihilated="observer-source restrictions and Sidon-defect labels decide which pair tensions survive",
        residual="nonfree second-order pair-tournament class",
        destroys="individual speed names after passing to pair differences",
        anchors=("HYP-2027", "HYP-1976", "HYP-2486"),
        vector=(5, 3, 5, 3, 3, 4, 2, 3, 4),
    ),
    Carrier(
        name="carry_residue_cocycle",
        kind="CRT/carry side-information cocycle",
        cocycle="carry in v=r+27k, plus owner route and Cprime tick window",
        base="Res_27 proof atom quotient with primitive floor atom",
        differential="changing lifts changes carry while residue shadow can remain fixed",
        exact_or_annihilated="conditional-independence theorem would make raw integer presentation irrelevant given carry/owner/window data",
        residual="mixed floor-vs-strict shadow route in one residue fiber",
        destroys="integer presentation and high address entropy",
        anchors=("HYP-2171", "HYP-2167", "HYP-2168"),
        vector=(3, 4, 3, 3, 4, 4, 3, 4, 4),
    ),
    Carrier(
        name="path_homology_witness_cocycle",
        kind="literal tournament homology cocycle",
        cocycle="H1 witness cocycle z normalized by <z, bad boundary>=1 after tournament flips",
        base="tournament path-homology chain complex and deleted-vertex restrictions",
        differential="bad transitive-triple boundary creates or kills the witness class",
        exact_or_annihilated="restriction/cokernel tests show exactly one missing direction in bad vertices",
        residual="nonzero homology class rather than scalar tournament statistic",
        destroys="LRC packet labels unless transported through a state-lift functor",
        anchors=("HYP-362", "HYP-366", "HYP-372"),
        vector=(5, 2, 5, 4, 3, 2, 2, 5, 5),
    ),
    Carrier(
        name="raw_scalar_shadow",
        kind="negative control",
        cocycle="none; scalar shadow after the cochain address has been forgotten",
        base="component counts, qdiv-only values, strip counts, or tournament isomorphism class",
        differential="undefined because the base has already quotiented the switch direction",
        exact_or_annihilated="valid only after a previous carrier proves fiber constancy or annihilation",
        residual="mixed fibers masquerading as evidence",
        destroys="local cocycle, endpoint owner, exact period, and proof route",
        anchors=("HYP-2594", "HYP-2978", "HYP-2990"),
        vector=(1, 1, 1, 1, 1, 1, 1, 2, 1),
    ),
)


DOCS = (
    "05-knowledge/hypotheses/HYP-2991-lrc14-haar-zipper-cocycle.md",
    "05-knowledge/hypotheses/HYP-2990-abstract-zipper-theorem-atlas.md",
    "05-knowledge/hypotheses/HYP-2989-lrc14-haar-product-discrepancy-tournament-tiling.md",
    "05-knowledge/hypotheses/HYP-2988-lrc14-exposure-poset-proof-pass.md",
    "05-knowledge/hypotheses/HYP-2987-lrc14-certificate-handoff-atlas.md",
    "05-knowledge/hypotheses/HYP-2986-lrc14-tope-wall-cocircuit.md",
    "05-knowledge/hypotheses/HYP-2985-lrc14-admissible-smoothing-dispatcher.md",
    "05-knowledge/hypotheses/HYP-2984-lrc14-kernel-homotopy-boundary-defect.md",
    "05-knowledge/hypotheses/HYP-2982-lrc14-analytic-sieve-packet-weights.md",
    "05-knowledge/hypotheses/HYP-2979-lrc14-ramanujan-exact-period-projector.md",
    "05-knowledge/hypotheses/HYP-2978-lrc14-ramanujan-divisor-quotient-guardrails.md",
    "05-knowledge/hypotheses/HYP-2974-lrc14-fourier-toeplitz-dual-scout.md",
    "05-knowledge/hypotheses/HYP-2171-lrc-n14-information-bottleneck.md",
    "05-knowledge/hypotheses/INDEX.md",
)


TERMS = (
    "cocycle",
    "coboundary",
    "cohomology",
    "boundary",
    "cocircuit",
    "exact",
    "residual",
    "quotient",
    "fiber",
    "state-lift",
    "Ramanujan",
    "Fejer",
    "Haar",
    "endpoint",
    "smoothing",
    "carry",
    "tension",
)


def read_docs() -> tuple[dict[str, str], list[str]]:
    found: dict[str, str] = {}
    missing: list[str] = []
    for rel in DOCS:
        path = ROOT / rel
        if path.exists():
            found[rel] = path.read_text(encoding="utf-8", errors="ignore")
        else:
            missing.append(rel)
    return found, missing


def concept_counts(docs: dict[str, str]) -> tuple[Counter[str], Counter[tuple[str, str]]]:
    hits: Counter[str] = Counter()
    co: Counter[tuple[str, str]] = Counter()
    for text in docs.values():
        lower = text.lower()
        present = []
        for term in TERMS:
            count = lower.count(term.lower())
            if count:
                hits[term] += count
                present.append(term)
        for a, b in combinations(sorted(set(present)), 2):
            co[(a, b)] += 1
    return hits, co


def orient(a: int, b: int) -> bool:
    """Return True when carrier a dominates carrier b."""
    av = CARRIERS[a].vector
    bv = CARRIERS[b].vector
    a_wins = sum(x > y for x, y in zip(av, bv))
    b_wins = sum(y > x for x, y in zip(av, bv))
    if a_wins != b_wins:
        return a_wins > b_wins
    adot = sum(av)
    bdot = sum(bv)
    if adot != bdot:
        return adot > bdot
    return a < b


def tournament() -> dict[int, set[int]]:
    adj = {i: set() for i in range(len(CARRIERS))}
    for i, j in combinations(range(len(CARRIERS)), 2):
        if orient(i, j):
            adj[i].add(j)
        else:
            adj[j].add(i)
    return adj


def directed_3cycles(adj: dict[int, set[int]]) -> int:
    total = 0
    for a, b, c in combinations(adj, 3):
        if b in adj[a] and c in adj[b] and a in adj[c]:
            total += 1
        if c in adj[a] and b in adj[c] and a in adj[b]:
            total += 1
    return total


def scc_sizes(adj: dict[int, set[int]]) -> list[int]:
    radj = {v: set() for v in adj}
    for v, outs in adj.items():
        for w in outs:
            radj[w].add(v)

    def reach(graph: dict[int, set[int]], start: int) -> set[int]:
        seen = {start}
        stack = [start]
        while stack:
            v = stack.pop()
            for w in graph[v]:
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
        return seen

    remaining = set(adj)
    sizes = []
    while remaining:
        start = next(iter(remaining))
        comp = reach(adj, start) & reach(radj, start)
        sizes.append(len(comp))
        remaining -= comp
    return sorted(sizes, reverse=True)


def hamiltonian_path_count(adj: dict[int, set[int]]) -> int:
    n = len(adj)
    dp: dict[tuple[int, int], int] = {}
    for i in range(n):
        dp[(1 << i, i)] = 1
    for mask in range(1 << n):
        for last in range(n):
            ways = dp.get((mask, last), 0)
            if not ways:
                continue
            for nxt in adj[last]:
                bit = 1 << nxt
                if mask & bit:
                    continue
                dp[(mask | bit, nxt)] = dp.get((mask | bit, nxt), 0) + ways
    full = (1 << n) - 1
    return sum(dp.get((full, i), 0) for i in range(n))


def transitive_path(adj: dict[int, set[int]]) -> list[str]:
    scores = {i: len(adj[i]) for i in adj}
    return [CARRIERS[i].name for i in sorted(adj, key=lambda i: (-scores[i], i))]


def print_table() -> None:
    print("Cocycle carrier table")
    print("---------------------")
    for carrier in CARRIERS:
        anchors = ", ".join(carrier.anchors)
        print(f"* {carrier.name} [{carrier.kind}]")
        print(f"  cocycle: {carrier.cocycle}")
        print(f"  base: {carrier.base}")
        print(f"  differential: {carrier.differential}")
        print(f"  exact/annihilated: {carrier.exact_or_annihilated}")
        print(f"  residual: {carrier.residual}")
        print(f"  destroys: {carrier.destroys}")
        print(f"  anchors: {anchors}")
        print(f"  vector: {dict(zip(DIMENSIONS, carrier.vector))}")
        print()


def main() -> None:
    docs, missing = read_docs()
    hits, co = concept_counts(docs)
    adj = tournament()
    scores = {i: len(adj[i]) for i in adj}
    score_hist = Counter(scores.values())
    path = transitive_path(adj)

    print("S167 LRC14 cocycle-sheaf unification")
    print("=" * 72)
    print()
    print("Thesis")
    print("------")
    print("Every useful LRC14 proof carrier can be read as one of four things:")
    print("  1. a cocycle coordinate that a quotient must retain;")
    print("  2. a coboundary/exactness statement that cancels it;")
    print("  3. a restriction or zipper handoff to a smaller packet/family;")
    print("  4. a named cohomology class, the residual F7/THM-572 debt.")
    print()
    print("The local Haar zeta is the model example, but not the whole story.")
    print("Endpoint currents, exact-period Ramanujan kernels, Fejer dual pairings,")
    print("smoothing boundary defects, carry side information, pair tensions, and")
    print("path-homology witnesses all fit the same controlled-kernel grammar.")
    print()

    print("[1] Recent-document cocycle vocabulary")
    print(f"  documents scanned: {len(docs)}")
    if missing:
        print("  missing documents:")
        for rel in missing:
            print(f"    {rel}")
    print("  term hits:")
    for term, count in hits.most_common():
        print(f"    {term:12s} {count:5d}")
    print("  strongest co-occurrences:")
    for (a, b), count in co.most_common(14):
        print(f"    {a:12s} <-> {b:12s} docs={count}")
    print()

    print("[2] Cocycle equations / proof meanings")
    equations = [
        (
            "fixed-margin Haar",
            "zeta(T+lambda*[[1,-1],[-1,1]]) = zeta(T)+4 lambda",
            "margins are 0-data; zeta is the forgotten switch coordinate",
        ),
        (
            "pair tension",
            "w_ij=v_j-v_i, so w_ij+w_jk+w_ki=0",
            "pair tournaments are not free: their labels are coboundaries",
        ),
        (
            "endpoint current",
            "AP/GW taut owner pairs have sum 0 mod 14",
            "boundary equality is a cocircuit/current, not scalar zero",
        ),
        (
            "Ramanujan projector",
            "c_q(a-b)=sum_{d|gcd(q,a-b)} d*mu(q/d)",
            "exact-period components are character cocycles on Z/qZ",
        ),
        (
            "Fejer dual",
            "<C_S-1, |P|^2> < 0",
            "a negative pairing is a dual coboundary certificate for safety",
        ),
        (
            "zipper handoff",
            "d(packet debt)=sum failed teeth routed to named exits",
            "proof gluing is exactness of the handoff complex at C1",
        ),
    ]
    for name, equation, meaning in equations:
        print(f"  {name:20s}: {equation}")
        print(f"  {'':20s}  {meaning}")
    print()

    print("[3] Tournament Analysis over cocycle carriers")
    print("  vertices: cocycle carriers / proof obligations, not runners")
    print("  pairwise observable: retained coordinates in")
    print("    " + ", ".join(DIMENSIONS))
    print("  switch/gauge: majority comparison of retention vector")
    print("  tie Hamiltonian path: " + " > ".join(carrier.name for carrier in CARRIERS))
    print()
    for i in sorted(adj, key=lambda j: (-scores[j], j)):
        carrier = CARRIERS[i]
        print(f"  {carrier.name:34s} score={scores[i]:2d} vector={carrier.vector}")
    print(f"  score_hist={dict(sorted(score_hist.items()))}")
    print(f"  directed_3cycles={directed_3cycles(adj)}")
    print(f"  scc_sizes={scc_sizes(adj)}")
    print(f"  hamiltonian_path_count={hamiltonian_path_count(adj)}")
    print("  canonical_path=" + " > ".join(path))
    print()

    print("[4] Candidate exactness theorem")
    print("  Define a labelled packet cochain complex:")
    print("    C0 = packet fibers with exact M/qdiv, open-boundary state, owners,")
    print("         exact-period labels, and proof-route labels.")
    print("    C1 = local emitted cocycles: zeta switches, endpoint currents,")
    print("         Ramanujan phases, Fejer dual debts, smoothing defects, carries,")
    print("         pair tensions, and handoff obligations.")
    print("    C2 = incompatibilities between exits: unnamed F7 classes,")
    print("         state-lift debt, or quotient fibers mixing the LRC predicate.")
    print("  LRC14 would follow from exactness at C1 for primitive non-AP/GW rows:")
    print("    ker(d1 on labelled packets) = im(d0 from known certificates)")
    print("  with AP/GW as declared boundary cohomology and THM-572/F7 as the only")
    print("  named residual summand.  In plain terms: every cocycle emitted by a")
    print("  possible counterexample is cancelled, exposed, restricted, or named.")
    print()

    print("[5] Assumption challenge")
    print("  Alternate vertices considered: runners, arcs, dyadic rectangles,")
    print("  endpoint walls, residues, exact-period modes, Fejer atoms, smoothing")
    print("  policies, pair tensions, carry lifts, path-homology classes, and proof")
    print("  obligations.")
    print("  Chosen vertices: cocycle carriers, because they preserve the predicate")
    print("  of interest under quotienting.")
    print("  Preserved: local switch sign, packet base, endpoint/current labels,")
    print("  exact-period/dual labels, handoff arrows, and named residual sectors.")
    print("  Destroyed: raw runner names, scalar component counts, unmarked")
    print("  tournament isomorphism classes, and squarefree/strip shadows without")
    print("  a reconstruction or annihilation certificate.")


if __name__ == "__main__":
    main()
