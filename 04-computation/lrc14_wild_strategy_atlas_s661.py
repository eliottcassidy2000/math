#!/usr/bin/env python3
"""S661: wild proof-strategy atlas for the LRC n=14 seam.

This is not an exact verifier.  It is a structured brainstorming artifact for
HYP-2237: take the current LRC n=14 quotient tower seriously, then import the
repo's strongest side-channel methods (deck derivatives, Kakeya pinned fibers,
union-closed pressure, sheaf gluing, pincer/rigidity, and automata) as possible
ways to prove the remaining carry-owner no-leak theorem.

Tournament Analysis:
  Vertices are proof strategies / carriers, not runners.
  The pairwise observable is a majority vote over proof leverage,
  near-term actionability, side-channel retention, repo support, wild transfer
  value, computational cheapness, risk inverse, and theorem specificity.
  The tie Hamiltonian path is the listed strategy order.

Assumption challenge:
  Candidate vertices considered: runners, speeds, gaps, fixed circle sections,
  wall-crossing events, residues, cover arcs, endpoint owners, carry cocycles,
  deletion derivatives, sheaf charts, automaton states, proof obligations, and
  strategies.  This script chooses strategies because the user's request is a
  proof-route search.  The preserved predicate is "could this route prove the
  LRC n=14 no-leak/carry-owner seam?"  The destroyed data are exact speed-row
  identities and exact phase geometry, which must be reattached by any route
  that graduates from atlas to theorem.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, replace
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCAN_DIRS = [
    "00-navigation",
    "01-canon",
    "04-computation",
    "05-knowledge",
    "07-reflections",
]
SCAN_EXTS = {".md", ".py", ".lean", ".txt"}
MAX_READ_BYTES = 900_000
SKIP_PARTS = {".git", ".lake", "__pycache__", ".pytest_cache", "node_modules"}


@dataclass(frozen=True)
class Strategy:
    name: str
    thesis: str
    keywords: tuple[str, ...]
    proof_leverage: int
    near_term: int
    side_retention: int
    wild_transfer: int
    cheapness: int
    risk_inverse: int
    theorem_specificity: int
    next_experiment: str
    danger: str
    repo_hits: int = 0
    repo_files: int = 0
    top_files: tuple[tuple[str, int], ...] = ()


BASE_STRATEGIES: tuple[Strategy, ...] = (
    Strategy(
        name="carry_deletion_derivative",
        thesis=(
            "Reconstruct the carry cocycle from local deletions/perturbations: "
            "if all one-runner and two-runner derivatives are strict unless "
            "they match scalar AP/Vstar, global floor carries are forced."
        ),
        keywords=("carry", "deletion", "derivative", "loss profile", "Res_27", "HYP-2167"),
        proof_leverage=10,
        near_term=9,
        side_retention=10,
        wild_transfer=7,
        cheapness=7,
        risk_inverse=8,
        theorem_specificity=10,
        next_experiment=(
            "For each HYP-2164 floor/near-floor atom, delete one coordinate or "
            "toggle one carry and record exact changes in M, active wall pairs, "
            "gcd-shell mass, n-clock residues, and owner route."
        ),
        danger="Derivative data may not be visible from the quotient without a marking theorem.",
    ),
    Strategy(
        name="owner_concurrency_jackknife",
        thesis=(
            "Use the Kakeya/Falconer lesson: scalar support saturates early, so "
            "hold odd-wall/pair-sum support fixed and jackknife owner, pinned, "
            "and concurrency labels until the false floors reveal the missing channel."
        ),
        keywords=("Kakeya", "Falconer", "pinned", "concurrency", "owner", "jackknife", "no-leak"),
        proof_leverage=10,
        near_term=8,
        side_retention=10,
        wild_transfer=9,
        cheapness=8,
        risk_inverse=8,
        theorem_specificity=9,
        next_experiment=(
            "Build an LRC analogue of the F_5 line-family audit: same odd wall "
            "and C=27 support, varying owner/carry packets, then measure which "
            "packets decide floor versus strict."
        ),
        danger="A finite-field metaphor can overfit unless every retained label has an LRC predicate.",
    ),
    Strategy(
        name="apex_sheaf_gluing",
        thesis=(
            "Treat cheap-pair certificates as local sections over unit, nonunit, "
            "and apex charts; prove every failed gluing creates endpoint-owner "
            "positive measure."
        ),
        keywords=("sheaf", "gluing", "apex", "cheap pair", "endpoint owner", "section"),
        proof_leverage=9,
        near_term=8,
        side_retention=10,
        wild_transfer=8,
        cheapness=7,
        risk_inverse=7,
        theorem_specificity=9,
        next_experiment=(
            "Extend S579's section table from one-lifts to carry derivatives: "
            "AP, Vstar, gcd-3/gcd-9 moves, and minimal apex carries with "
            "restriction maps between charts."
        ),
        danger="Sheaf language must remain a finite table of certificates, not abstraction fog.",
    ),
    Strategy(
        name="three_state_middle_automaton",
        thesis=(
            "Refine tournament edges into L/M/R automata.  A counterexample must "
            "keep too many proof cells terminal-M; prove no all-middle circuit "
            "survives owner, CRT, and endpoint peeling."
        ),
        keywords=("three-state", "automaton", "L/M/R", "middle", "endpoint", "residue automaton"),
        proof_leverage=9,
        near_term=7,
        side_retention=9,
        wild_transfer=9,
        cheapness=8,
        risk_inverse=7,
        theorem_specificity=8,
        next_experiment=(
            "Compile owner/carry obligations into L/M/R states and build the "
            "middle graph for the S609/S654 residual; extract closed-middle "
            "circuits or certify none."
        ),
        danger="Automata can silently reintroduce the global CRT if states are too coarse.",
    ),
    Strategy(
        name="two_block_helly_extractor",
        thesis=(
            "Do not solve the full large-owner CRT; extract singleton/pair "
            "determinant Helly contradictions for every bounded residual."
        ),
        keywords=("two-block", "determinant", "Helly", "CRT", "large-owner", "singleton"),
        proof_leverage=9,
        near_term=8,
        side_retention=8,
        wild_transfer=7,
        cheapness=9,
        risk_inverse=8,
        theorem_specificity=9,
        next_experiment=(
            "Route S654 apex-carry rows that survive cheap exits into the "
            "S599 determinant language and request minimal empty subfamilies."
        ),
        danger="This mainly attacks the Cprime/multiple branch, not every carry floor shadow.",
    ),
    Strategy(
        name="pincer_grip_ledger",
        thesis=(
            "Pair-sum pinches are the bottom jaw; shields, anchors, and endpoint "
            "owners are the top jaw.  Prove every nonmeeting pincer leaves a "
            "labelled escape ledger."
        ),
        keywords=("pincer", "pinch", "shield", "anchor", "endpoint", "observer", "grip"),
        proof_leverage=8,
        near_term=8,
        side_retention=9,
        wild_transfer=8,
        cheapness=8,
        risk_inverse=8,
        theorem_specificity=8,
        next_experiment=(
            "For each HYP-2164 strict row type, store the first pincer meeting "
            "or the exact escape label: D|v shield, endpoint anchor, Cprime, "
            "or positive-measure owner."
        ),
        danger="May only repackage existing pair-pinch certificates unless it exports new escape labels.",
    ),
    Strategy(
        name="union_closed_obligation_pressure",
        thesis=(
            "Model proof obligations as a closure family.  Frequency-like scalar "
            "counts collide; set pressure may be the missing repair that decides "
            "which obligation closes under owner/carry union."
        ),
        keywords=("union-closed", "set pressure", "pressure", "proof obligations", "closure"),
        proof_leverage=7,
        near_term=7,
        side_retention=8,
        wild_transfer=9,
        cheapness=8,
        risk_inverse=6,
        theorem_specificity=6,
        next_experiment=(
            "Create the family of obligation packets for AP/Vstar/carry moves "
            "and compute pressure(A)=sum_B |A union B|; look for scalar twins "
            "separated by pressure."
        ),
        danger="The closure predicate must be defined carefully or pressure becomes decorative.",
    ),
    Strategy(
        name="forcing_absoluteness_audit",
        thesis=(
            "Treat the least-positive Res_27 section as a constructible model and "
            "carry/owner lifts as generic extensions; prove floor classification "
            "is absolute exactly for the proposed sufficient statistic."
        ),
        keywords=("forcing", "absoluteness", "scalar twins", "sufficient statistic", "generic"),
        proof_leverage=8,
        near_term=7,
        side_retention=10,
        wild_transfer=8,
        cheapness=7,
        risk_inverse=7,
        theorem_specificity=8,
        next_experiment=(
            "Generate scalar twins with identical Res_27 proof atoms but varied "
            "carry/owner packets; mark which predicates are absolute under the "
            "lift and which change."
        ),
        danger="The set-theory analogy is method-only; it cannot replace the finite lift theorem.",
    ),
    Strategy(
        name="unit_distance_impairment_transfer",
        thesis=(
            "Borrow impairment spectroscopy: deliberately damage LRC channels "
            "one at a time and price released shallow residues or false floors."
        ),
        keywords=("impairment", "damage", "direction-drop", "channel ledger", "unit-distance", "jackknife"),
        proof_leverage=7,
        near_term=8,
        side_retention=8,
        wild_transfer=8,
        cheapness=9,
        risk_inverse=8,
        theorem_specificity=7,
        next_experiment=(
            "Extend S624 from witness-orbit deletion to carry-owner deletion: "
            "remove one channel from AP/Vstar/2AP and record released residues, "
            "redundancy prices, and new cheap pairs."
        ),
        danger="A damage ledger gives necessity/pricing before sufficiency.",
    ),
    Strategy(
        name="orbit_boundary_rigidity",
        thesis=(
            "Split n=14 into unit sheet, gcd-boundary sheet, and apex chart; "
            "prove every boundary monodromy class is owned by endpoint, pincer, "
            "or CRT labels."
        ),
        keywords=("orbit", "rigidity", "monodromy", "boundary", "gcd", "apex", "unit sheet"),
        proof_leverage=8,
        near_term=7,
        side_retention=9,
        wild_transfer=8,
        cheapness=7,
        risk_inverse=7,
        theorem_specificity=8,
        next_experiment=(
            "Compute the monodromy of active wall-pair labels around the "
            "G=<2,-1> shell fold and mark where unit/gcd/apex labels mix."
        ),
        danger="Orbit purity must be checked at the observer-coupled label level.",
    ),
    Strategy(
        name="depth_polynomial_root_atlas",
        thesis=(
            "Use THM-406 directly: p0=0 is a root of the depth polynomial.  "
            "Compare AP, Vstar, carry moves, and strict controls by all-order "
            "overlap signatures instead of low moments."
        ),
        keywords=("depth polynomial", "p_0", "anti-Poisson", "coimage", "overlap", "THM-406"),
        proof_leverage=8,
        near_term=6,
        side_retention=8,
        wild_transfer=8,
        cheapness=5,
        risk_inverse=6,
        theorem_specificity=7,
        next_experiment=(
            "For the smallest carry derivative rows, compute the depth "
            "generating polynomial at delta=1/14 and locate the first "
            "coefficient/order that separates floor from strict."
        ),
        danger="All-order overlap data can become expensive without a compression lemma.",
    ),
    Strategy(
        name="matroid_fold_circuit_no_leak",
        thesis=(
            "Treat visible folds and hidden 4-term packets as matroid circuits; "
            "prove circuit-free rows have discrepancy margin and circuit rows "
            "route through denominator shields."
        ),
        keywords=("matroid", "circuit", "fold", "3-term", "4-term", "denominator shield"),
        proof_leverage=7,
        near_term=6,
        side_retention=8,
        wild_transfer=8,
        cheapness=6,
        risk_inverse=6,
        theorem_specificity=7,
        next_experiment=(
            "Add a circuit signature to the fixed odd-wall scan: visible folds, "
            "hidden sums, denominator shields, and whether each circuit is "
            "observer-coupled or observer-blind."
        ),
        danger="Prior sessions showed 4-term richness alone is safe; the circuit type must be observer-coupled.",
    ),
    Strategy(
        name="tropical_minplus_safe_polytope",
        thesis=(
            "View M(V) as a tropical maximin over pair-sum wall normals; prove "
            "non-AP/Vstar cells have a positive tropical slack certificate."
        ),
        keywords=("tropical", "min-plus", "polytope", "maximin", "wall", "slack"),
        proof_leverage=6,
        near_term=5,
        side_retention=7,
        wild_transfer=9,
        cheapness=5,
        risk_inverse=5,
        theorem_specificity=6,
        next_experiment=(
            "Build a tiny tropical chart for AP, Vstar, unit-shift AP, and "
            "single carry moves: active normals, slack cone, and degeneracy rank."
        ),
        danger="May produce beautiful geometry without an arithmetic lift/CRT bound.",
    ),
    Strategy(
        name="fourier_major_minor_shell_split",
        thesis=(
            "Split safe measure into major C=27 shell packets and minor lift "
            "oscillations, mirroring circle-method logic; prove major packets "
            "dominate unless the row is AP/Vstar."
        ),
        keywords=("Fourier", "major", "minor", "shell", "character", "Paley", "spectral"),
        proof_leverage=6,
        near_term=5,
        side_retention=7,
        wild_transfer=8,
        cheapness=4,
        risk_inverse=5,
        theorem_specificity=6,
        next_experiment=(
            "Compute C=27 character sums for AP/Vstar/carry moves and test "
            "whether the gcd-3/gcd-9 packets have a spectral separator."
        ),
        danger="Analytic estimates may be much weaker than existing finite certificates.",
    ),
    Strategy(
        name="rooted_perspective_source_quotient",
        thesis=(
            "Use observer-source rigidity: quotient source/sink perspectives only "
            "after threshold labels and incident words are retained."
        ),
        keywords=("rooted perspective", "source", "sink", "observer", "threshold", "A000568", "Burnside"),
        proof_leverage=6,
        near_term=6,
        side_retention=8,
        wild_transfer=7,
        cheapness=7,
        risk_inverse=7,
        theorem_specificity=6,
        next_experiment=(
            "Decorate source/sink rooted perspective classes with odd-wall, "
            "endpoint-owner, and carry labels; check whether decorated fibers "
            "become pure on known n=14 residuals."
        ),
        danger="Unmarked A000568-style quotients are known to leak LRC predicates.",
    ),
)


def iter_text_files() -> list[Path]:
    files: list[Path] = []
    for dirname in SCAN_DIRS:
        base = ROOT / dirname
        if not base.exists():
            continue
        for path in base.rglob("*"):
            if any(part in SKIP_PARTS for part in path.parts):
                continue
            if path.is_file() and path.suffix.lower() in SCAN_EXTS:
                try:
                    if path.stat().st_size <= MAX_READ_BYTES:
                        files.append(path)
                except OSError:
                    pass
    return sorted(files)


def scan_all_strategies(strategies: tuple[Strategy, ...], files: list[Path]) -> tuple[Strategy, ...]:
    """Read each file once and update every strategy's keyword counters."""
    per_file = {s.name: Counter() for s in strategies}
    totals = Counter()
    keyword_table = {s.name: tuple(k.lower() for k in s.keywords) for s in strategies}

    for path in files:
        try:
            text = path.read_text(encoding="utf-8", errors="ignore").lower()
        except OSError:
            continue
        rel = path.relative_to(ROOT).as_posix()
        for s in strategies:
            hits = sum(text.count(k) for k in keyword_table[s.name])
            if hits:
                per_file[s.name][rel] = hits
                totals[s.name] += hits

    scanned: list[Strategy] = []
    for s in strategies:
        top = tuple(per_file[s.name].most_common(4))
        scanned.append(
            replace(
                s,
                repo_hits=totals[s.name],
                repo_files=len(per_file[s.name]),
                top_files=top,
            )
        )
    return tuple(scanned)


def repo_support_score(strategy: Strategy, all_strategies: tuple[Strategy, ...]) -> int:
    max_hits = max(s.repo_hits for s in all_strategies) or 1
    max_files = max(s.repo_files for s in all_strategies) or 1
    hit_score = round(10 * strategy.repo_hits / max_hits)
    file_score = round(10 * strategy.repo_files / max_files)
    return max(1, round((hit_score + file_score) / 2))


def profile(strategy: Strategy, all_strategies: tuple[Strategy, ...]) -> dict[str, int]:
    return {
        "proof_leverage": strategy.proof_leverage,
        "near_term": strategy.near_term,
        "side_retention": strategy.side_retention,
        "repo_support": repo_support_score(strategy, all_strategies),
        "wild_transfer": strategy.wild_transfer,
        "cheapness": strategy.cheapness,
        "risk_inverse": strategy.risk_inverse,
        "theorem_specificity": strategy.theorem_specificity,
    }


def majority_winner(a: Strategy, b: Strategy, all_strategies: tuple[Strategy, ...]) -> str:
    pa = profile(a, all_strategies)
    pb = profile(b, all_strategies)
    votes_a = 0
    votes_b = 0
    for key in pa:
        if pa[key] > pb[key]:
            votes_a += 1
        elif pb[key] > pa[key]:
            votes_b += 1
    if votes_a > votes_b:
        return a.name
    if votes_b > votes_a:
        return b.name
    return a.name


def build_tournament(strategies: tuple[Strategy, ...]) -> dict[str, set[str]]:
    adj = {s.name: set() for s in strategies}
    for a, b in combinations(strategies, 2):
        winner = majority_winner(a, b, strategies)
        loser = b.name if winner == a.name else a.name
        adj[winner].add(loser)
    return adj


def directed_3cycles(adj: dict[str, set[str]]) -> int:
    nodes = sorted(adj)
    total = 0
    for a, b, c in combinations(nodes, 3):
        if b in adj[a] and c in adj[b] and a in adj[c]:
            total += 1
        elif c in adj[a] and b in adj[c] and a in adj[b]:
            total += 1
    return total


def scc_sizes(adj: dict[str, set[str]]) -> list[int]:
    nodes = sorted(adj)
    rev = {u: set() for u in nodes}
    for u in nodes:
        for v in adj[u]:
            rev[v].add(u)

    def reach(start: str, graph: dict[str, set[str]]) -> set[str]:
        seen = {start}
        stack = [start]
        while stack:
            u = stack.pop()
            for v in graph[u]:
                if v not in seen:
                    seen.add(v)
                    stack.append(v)
        return seen

    remaining = set(nodes)
    sizes: list[int] = []
    while remaining:
        u = min(remaining)
        comp = reach(u, adj) & reach(u, rev)
        sizes.append(len(comp))
        remaining -= comp
    return sorted(sizes, reverse=True)


def hamiltonian_path_count(adj: dict[str, set[str]]) -> int:
    nodes = sorted(adj)
    n = len(nodes)
    index = {node: i for i, node in enumerate(nodes)}
    out_masks = [0] * n
    for u in nodes:
        ui = index[u]
        for v in adj[u]:
            out_masks[ui] |= 1 << index[v]
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            val = dp[mask][last]
            if not val:
                continue
            remaining = (~mask) & ((1 << n) - 1)
            nxts = out_masks[last] & remaining
            while nxts:
                bit = nxts & -nxts
                j = bit.bit_length() - 1
                dp[mask | bit][j] += val
                nxts -= bit
    return sum(dp[-1])


def rank_strategies(strategies: tuple[Strategy, ...]) -> list[tuple[int, Strategy, dict[str, int]]]:
    rows: list[tuple[int, Strategy, dict[str, int]]] = []
    for s in strategies:
        p = profile(s, strategies)
        score = (
            3 * p["proof_leverage"]
            + 3 * p["theorem_specificity"]
            + 2 * p["side_retention"]
            + 2 * p["near_term"]
            + p["repo_support"]
            + p["wild_transfer"]
            + p["cheapness"]
            + p["risk_inverse"]
        )
        rows.append((score, s, p))
    return sorted(rows, key=lambda row: (-row[0], row[1].name))


def print_search_ledger(strategies: tuple[Strategy, ...], files: list[Path]) -> None:
    print("A. Repo search ledger")
    print(f"  scanned_files={len(files)} dirs={','.join(SCAN_DIRS)} max_read_bytes={MAX_READ_BYTES}")
    print("  top keyword-supported strategies:")
    for s in sorted(strategies, key=lambda x: (-x.repo_files, -x.repo_hits, x.name))[:8]:
        print(f"    {s.name}: hits={s.repo_hits}, files={s.repo_files}")
        for rel, hits in s.top_files[:2]:
            print(f"      {hits:4d}  {rel}")
    print()


def print_ranked_atlas(strategies: tuple[Strategy, ...]) -> None:
    print("B. Strategy atlas ranking")
    for rank, (score, s, p) in enumerate(rank_strategies(strategies), start=1):
        print(
            f"  {rank:2d}. {s.name:34s} score={score:3d} "
            f"lev={p['proof_leverage']} spec={p['theorem_specificity']} "
            f"side={p['side_retention']} near={p['near_term']} "
            f"repo={p['repo_support']} wild={p['wild_transfer']}"
        )
        print(f"      thesis: {s.thesis}")
        print(f"      next:   {s.next_experiment}")
        print(f"      risk:   {s.danger}")
    print()


def print_tournament(strategies: tuple[Strategy, ...]) -> None:
    adj = build_tournament(strategies)
    scores = Counter(len(adj[s.name]) for s in strategies)
    print("C. Tournament Analysis over proof strategies")
    print("  vertices are proof strategies / carriers, not runners")
    print("  pairwise gauges: proof leverage, specificity, side retention, near-term action, repo support, wild transfer, cheapness, risk inverse")
    print(f"  score_hist={dict(sorted(scores.items()))}")
    print(f"  directed_3cycles={directed_3cycles(adj)}")
    print(f"  scc_sizes={scc_sizes(adj)}")
    print(f"  hamiltonian_paths={hamiltonian_path_count(adj)}")
    leaders = sorted(strategies, key=lambda s: (-len(adj[s.name]), s.name))[:6]
    print("  majority leaders:")
    for s in leaders:
        print(f"    {s.name}: outscore={len(adj[s.name])}")
    print()


def print_synthesis(strategies: tuple[Strategy, ...]) -> None:
    ranked = [s for _, s, _ in rank_strategies(strategies)]
    print("D. S661 synthesis")
    print("  strongest near-term theorem shape:")
    print(
        "    no-leak owner derivative theorem = fixed odd wall + C=27 gcd shell + "
        "carry/owner/deletion derivatives force AP, Vstar, 2AP, or strict looseness"
    )
    print("  four concrete proof moves:")
    print(f"    1. {ranked[0].name}: {ranked[0].next_experiment}")
    print(f"    2. {ranked[1].name}: {ranked[1].next_experiment}")
    print(f"    3. {ranked[2].name}: {ranked[2].next_experiment}")
    print(f"    4. {ranked[3].name}: {ranked[3].next_experiment}")
    print("  wacky but plausible reserve routes:")
    for s in ranked[-4:]:
        print(f"    - {s.name}: {s.thesis}")
    print()


def main() -> None:
    files = iter_text_files()
    scanned = scan_all_strategies(BASE_STRATEGIES, files)
    print("S661 LRC14 wild strategy atlas")
    print("=" * 40)
    print_search_ledger(scanned, files)
    print_ranked_atlas(scanned)
    print_tournament(scanned)
    print_synthesis(scanned)


if __name__ == "__main__":
    main()
