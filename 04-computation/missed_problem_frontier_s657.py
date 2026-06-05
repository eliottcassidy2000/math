#!/usr/bin/env python3
"""
S657: missed important problem frontier scout.

This is not a solver and not a claim that the repo is near solving famous
problems.  It is a triage atlas: compare important outside problems against the
repo's carrier methods, then rank targets where a small finite/restricted
result looks realistically reachable.

Tournament Analysis:
  Vertices are external problem clusters.
  Pairwise observable is a tuple:
    (repo_leverage, near_term_actionability, undercoverage, external_importance,
     source_confidence, -risk_of_raw_fame_chasing).
  The tie Hamiltonian path is the ranked list of candidate frontiers.

Assumption challenge:
  Candidate vertices considered included named conjectures, proof techniques,
  invariant families, finite toy models, and local carrier programs.  This
  script chooses problem clusters because the preserved predicate is "could the
  repo make a concrete stride here?"  It destroys fine-grained literature
  status and should therefore be followed by a focused literature check before
  any proof claim.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
SCAN_DIRS = [
    "00-navigation",
    "01-canon",
    "03-artifacts",
    "05-knowledge/hypotheses",
    "05-knowledge/variables",
    "07-reflections",
]
MAX_SCAN_BYTES = 300_000


@dataclass(frozen=True)
class Problem:
    slug: str
    label: str
    keywords: tuple[str, ...]
    external_importance: int
    repo_leverage: int
    near_term_actionability: int
    source_confidence: int
    risk_of_raw_fame_chasing: int
    source: str
    carrier_match: str
    first_experiment: str


PROBLEMS = [
    Problem(
        "graph_reconstruction",
        "Graph reconstruction / deck invariants",
        ("graph reconstruction", "reconstructible", "deck", "vertex-deleted"),
        5,
        5,
        5,
        5,
        2,
        "RIT reconstruction-number data and the classical Ulam graph reconstruction conjecture",
        "repo already has deletion profiles, endpoint transfer, old-projection kill/near-kill ledgers",
        "build a deck-to-carrier table for tournaments: which of H, OCF alpha, beta, SCC, score, and deletion-loss profiles are reconstructible from cards",
    ),
    Problem(
        "hadwiger_nelson",
        "Hadwiger-Nelson / chromatic unit-distance graphs",
        ("Hadwiger", "Nelson", "chromatic number of the plane", "unit-distance graph"),
        5,
        5,
        5,
        5,
        2,
        "MathWorld and Polymath16 records: plane chromatic number is 5, 6, or 7",
        "unit-distance work already has Moser carriers, traceable spines, direction support, and chi questions",
        "compute color-critical side-channel fingerprints for Moser/Sawin-style finite carriers: direction support, traceability, independence ratio, and tournament flip gauges",
    ),
    Problem(
        "kakeya_falconer",
        "Kakeya / Falconer distance-incidence geometry",
        ("Kakeya", "Falconer", "distance set", "decoupling", "projection theory"),
        5,
        5,
        4,
        5,
        3,
        "IPAM Kakeya overview and finite-field Kakeya references",
        "direction tubes, finite-field lines, unit-distance energy, LRC forbidden arcs, and Paley/finite-field carriers all align",
        "make a finite-field Kakeya/Falconer toy atlas over F_p^2: direction coverage, distance fibers, Paley character orientations, and tournament fingerprints",
    ),
    Problem(
        "sunflower_extractors",
        "Sunflower / robust sunflower / extractor set systems",
        ("sunflower", "robust sunflower", "extractor", "Erdos-Rado"),
        4,
        5,
        4,
        5,
        2,
        "Theory of Computing robust sunflower/extractor paper",
        "OCF is already a disjoint odd-cycle packing partition function; sunflower cores are a natural missing hypergraph side channel",
        "scan odd-cycle support hypergraphs for sunflower cores and robust-sunflower thresholds; compare H=63 exact-core rows with THM-025 near-kill rows",
    ),
    Problem(
        "union_closed",
        "Union-closed sets / Frankl frequency carriers",
        ("union-closed", "Frankl", "union closed", "closure system"),
        4,
        4,
        5,
        4,
        2,
        "arXiv:2208.03803 notes on union-closed set frequencies",
        "closure under union matches cauldron finite-sums, divisor/aliquot fixed carriers, and observer-source frequency predicates",
        "enumerate small closure systems as tournaments whose vertices are elements or sets; rank carriers by majority-frequency witness retention",
    ),
    Problem(
        "erdos_hajnal",
        "Erdos-Hajnal for tournaments and blow-ups",
        ("Erdos-Hajnal", "Erdos Hajnal", "transitive subtournament", "blow-up"),
        5,
        4,
        4,
        5,
        2,
        "Open Problem Garden EH page and tournament EH papers",
        "THM410 and S652 square-node substitution are literal blow-up/tournament carriers; H-mass can be a soft transitivity statistic",
        "test H-mass and OCF alpha spectra on H-free tournament families and lexicographic blow-ups before searching for exact transitive-subtournament bounds",
    ),
    Problem(
        "rota_ryser_rainbow",
        "Rota basis / Ryser-Brualdi-Stein rainbow transversals",
        ("Rota", "Ryser", "Brualdi", "Latin square", "rainbow", "transversal"),
        4,
        4,
        3,
        4,
        3,
        "SIAM Rota basis paper and Latin/rainbow transversal problem pages",
        "matroid circuits, endpoint-transfer matrices, cauldron schedules, and Latin/rainbow parity are underused in the repo",
        "build a rainbow-basis obstruction tournament where vertices are color/base obligations; report SCCs, Hall deficits, and Hamiltonian transversals",
    ),
    Problem(
        "caccetta_haggkvist",
        "Caccetta-Haggkvist short directed cycles",
        ("Caccetta", "Haggkvist", "Häggkvist", "short directed cycle", "minimum outdegree"),
        5,
        5,
        3,
        5,
        2,
        "arXiv:math/0605646 AIM-linked Caccetta-Haggkvist survey",
        "repo already has CH residue probes, LRC clock-return analogues, and directed-cycle obstruction language",
        "consolidate a CH-vs-LRC return-residue atlas and test r=3 extremal oriented graphs through deletion/owner side channels",
    ),
    Problem(
        "erdos_distinct_distances",
        "Erdos distinct distances and pinned distances",
        ("distinct distances", "pinned distance", "Guth", "Katz", "Szem", "Trotter"),
        5,
        4,
        3,
        5,
        3,
        "Guth-Katz Annals distinct distances theorem and active pinned-distance variants",
        "unit-distance energy work sees the opposite extremal face: many repeated distances rather than many distinct distances",
        "run dual energy ledgers on finite point sets: unit-distance support, distinct-distance support, pinned-distance distribution, and tournament flips between them",
    ),
    Problem(
        "graceful_tree",
        "Graceful tree labeling / difference-set carriers",
        ("graceful", "Rosa", "Kotzig", "tree labeling"),
        3,
        3,
        4,
        3,
        2,
        "standard graph-theory open problem lists",
        "difference labels resemble LRC shell witnesses and pair-sum/discriminant carriers, but tree structure is not central here yet",
        "prototype a graceful-labeling carrier for caterpillars and spiders using pair-difference tournaments and endpoint-transfer checks",
    ),
    Problem(
        "clay_control",
        "Clay Millennium controls",
        ("Riemann", "Hodge", "Navier", "Yang-Mills", "Birch", "Swinnerton", "P vs NP"),
        5,
        2,
        1,
        5,
        5,
        "Clay Mathematics Institute official list",
        "too broad for near-term repo progress except as analogy/control rows",
        "use only as negative controls: if a candidate ranks below Clay controls, it is probably not actionable",
    ),
]


def repo_text() -> str:
    chunks: list[str] = []
    for rel in SCAN_DIRS:
        base = ROOT / rel
        if not base.exists():
            continue
        for path in base.rglob("*"):
            if path.suffix.lower() not in {".md", ".py", ".lean", ".txt", ".html"}:
                continue
            try:
                if path.stat().st_size > MAX_SCAN_BYTES:
                    continue
            except OSError:
                continue
            try:
                chunks.append(path.read_text(errors="ignore"))
            except OSError:
                pass
    return "\n".join(chunks)


def count_hits(text: str, keywords: tuple[str, ...]) -> int:
    total = 0
    for keyword in keywords:
        total += len(re.findall(re.escape(keyword), text, flags=re.IGNORECASE))
    return total


def undercoverage_score(hits: int) -> int:
    if hits == 0:
        return 5
    if hits <= 5:
        return 5
    if hits <= 20:
        return 4
    if hits <= 80:
        return 3
    if hits <= 250:
        return 2
    return 1


def metric(problem: Problem, hits: int) -> tuple[int, int, int, int, int, int]:
    return (
        problem.repo_leverage,
        problem.near_term_actionability,
        undercoverage_score(hits),
        problem.external_importance,
        problem.source_confidence,
        -problem.risk_of_raw_fame_chasing,
    )


def route_tournament(rows: list[tuple[Problem, int]]) -> dict[str, object]:
    n = len(rows)
    wins = [0] * n
    edges: set[tuple[int, int]] = set()
    for i, j in combinations(range(n), 2):
        mi = metric(*rows[i])
        mj = metric(*rows[j])
        if mi > mj:
            winner = i
        elif mj > mi:
            winner = j
        else:
            winner = min(i, j)
        loser = j if winner == i else i
        edges.add((winner, loser))
        wins[winner] += 1

    c3 = 0
    for i, j, k in combinations(range(n), 3):
        if (i, j) in edges and (j, k) in edges and (k, i) in edges:
            c3 += 1
        if (j, i) in edges and (k, j) in edges and (i, k) in edges:
            c3 += 1

    adj = [[False] * n for _ in range(n)]
    for u, v in edges:
        adj[u][v] = True

    def reach(start: int, graph: list[list[bool]]) -> set[int]:
        seen = {start}
        stack = [start]
        while stack:
            u = stack.pop()
            for v, ok in enumerate(graph[u]):
                if ok and v not in seen:
                    seen.add(v)
                    stack.append(v)
        return seen

    radj = [[adj[j][i] for j in range(n)] for i in range(n)]
    remaining = set(range(n))
    scc_sizes: list[int] = []
    while remaining:
        start = min(remaining)
        comp = reach(start, adj) & reach(start, radj)
        scc_sizes.append(len(comp))
        remaining -= comp
    scc_sizes.sort(reverse=True)

    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last in range(n):
            val = dp[mask][last]
            if not val:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += val
    hpaths = sum(dp[-1])
    order = sorted(range(n), key=lambda i: (-wins[i], rows[i][0].slug))
    return {
        "top_order": [rows[i][0].slug for i in order],
        "score_hist": {score: wins.count(score) for score in sorted(set(wins))},
        "directed_3cycles": c3,
        "scc_sizes": scc_sizes,
        "hamiltonian_paths": hpaths,
    }


def main() -> None:
    text = repo_text()
    rows = [(problem, count_hits(text, problem.keywords)) for problem in PROBLEMS]
    rows.sort(key=lambda row: (metric(*row), -row[1], row[0].slug), reverse=True)

    print("S657 missed important problem frontier carrier atlas")
    print("=" * 72)
    print()

    print("A. Ranked frontier candidates")
    for rank, (problem, hits) in enumerate(rows, start=1):
        print(f"{rank:2d}. {problem.slug}")
        print(f"    label={problem.label}")
        print(f"    repo_hits={hits} undercoverage={undercoverage_score(hits)}")
        print(f"    metric={metric(problem, hits)}")
        print(f"    source={problem.source}")
        print(f"    carrier_match={problem.carrier_match}")
        print(f"    first_experiment={problem.first_experiment}")
    print()

    print("B. Tournament Analysis")
    fp = route_tournament(rows)
    for key, value in fp.items():
        print(f"  {key}={value}")
    print()

    print("C. Practical interpretation")
    print("  The best missed targets are not the most famous targets.  They are")
    print("  problems where the repo already has a carrier but has not built the")
    print("  correct finite toy model: decks for reconstruction, color-critical")
    print("  unit-distance graphs, Kakeya/Falconer finite-field incidence ledgers,")
    print("  sunflower cores inside OCF hypergraphs, and union-closed frequency")
    print("  carriers.  Clay-scale problems remain useful as negative controls,")
    print("  not as first attack surfaces.")


if __name__ == "__main__":
    main()
