#!/usr/bin/env python3
"""
S656: Continuum Hypothesis as a cardinal-shadow carrier atlas.

This script is deliberately an atlas, not a CH solver.  CH is independent of
ZFC: Godel's constructible universe gives a CH model, and Cohen forcing gives
models of not-CH, assuming the usual consistency background.  The repo-useful
question is: what does CH teach our finite quotient programs?

Answer: cardinality is a scalar shadow.  The proof object is the retained
model/forcing/generic side channel, just as LRC needs carry/owner labels and
tournament equinumerosity needs packet fibers.

Tournament Analysis:
  Vertices are proof routes / CH lessons, not sets or reals.
  Pairwise observable is a tuple of side-channel retention, repo transfer
  value, foundational centrality, and risk of misleading scalar analogy.

Assumption challenge:
  Candidate vertices considered: reals, subsets of N, ordinals, cardinals,
  Boolean algebras, forcing notions, generic filters, inner models, quotient
  maps, proof obligations, and repo carriers.  This script chooses proof-route
  carriers because they preserve the independence/absoluteness predicate; raw
  cardinal values destroy the relevant model information.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations


def finite_powerset_shadow(limit: int = 9) -> list[dict[str, int]]:
    rows = []
    for n in range(limit + 1):
        power = 2**n
        rows.append(
            {
                "n": n,
                "pow": power,
                "finite_successor": n + 1,
                "overshoot": power - (n + 1),
                "intermediate_sizes": max(0, power - n - 1),
            }
        )
    return rows


FOUNDATION_ROWS = [
    {
        "world": "Constructible universe L",
        "method": "inner model",
        "continuum": "aleph_1",
        "CH": "yes",
        "retained_channel": "definability / canonical well-order of reals",
        "repo_echo": "least-positive section with enough labels retained",
    },
    {
        "world": "Cohen-style forcing extension",
        "method": "generic extension",
        "continuum": "can be aleph_2 or larger",
        "CH": "no",
        "retained_channel": "generic reals added while ordinals are preserved",
        "repo_echo": "lift/carry fiber changes scalar floor behavior",
    },
    {
        "world": "Forcing-axiom worlds",
        "method": "axiom-extension program",
        "continuum": "often aleph_2",
        "CH": "no",
        "retained_channel": "maximality/compatibility side conditions",
        "repo_echo": "owner/Cprime certificate saturates residual cases",
    },
    {
        "world": "ZFC alone",
        "method": "base theory",
        "continuum": "undetermined",
        "CH": "independent",
        "retained_channel": "not enough model data to decide",
        "repo_echo": "raw quotient before sufficient statistic",
    },
]


REPO_BRIDGES = [
    {
        "thread": "HYP-2187 equinumerosity/equidecomposability",
        "scalar_shadow": "same number of objects",
        "missing_middle": "same count need not preserve H/beta/packet fibers",
        "retained_channel": "(H, beta1, odd-cycle packets)",
    },
    {
        "thread": "HYP-2185 infinite Go",
        "scalar_shadow": "finite cutoff count Z(1)",
        "missing_middle": "ordinal value/pole order invisible after scalarization",
        "retained_channel": "boundary fuel depth omega*r",
    },
    {
        "thread": "HYP-2171 LRC information bottleneck",
        "scalar_shadow": "Res_27 row address",
        "missing_middle": "floor predicate not determined by row address alone",
        "retained_channel": "owner route + carry cocycle + Cprime window",
    },
    {
        "thread": "HYP-2230 parity/carry",
        "scalar_shadow": "least-positive residue r",
        "missing_middle": "parity and apex obstruction live in carry k",
        "retained_channel": "v=r+27k, parity/apex/pair-sum congruences",
    },
    {
        "thread": "HYP-2231 sin/cos wall",
        "scalar_shadow": "(14,21) diagonal echo",
        "missing_middle": "active LRC wall is off-diagonal odd complement pairs",
        "retained_channel": "odd wall pairs + even slack + C=27 gcd shells",
    },
    {
        "thread": "OCF / H gaps",
        "scalar_shadow": "H=I(Omega,2)",
        "missing_middle": "same H can split by packet fibers; gaps can be structural",
        "retained_channel": "independence polynomial coefficients / SCC packets",
    },
]


@dataclass(frozen=True)
class Route:
    name: str
    side_channel: int
    repo_transfer: int
    foundation: int
    risk_low: int
    note: str

    @property
    def score_tuple(self) -> tuple[int, int, int, int]:
        return (self.side_channel, self.repo_transfer, self.foundation, self.risk_low)


ROUTES = [
    Route(
        "forcing/generic side-channel",
        5,
        5,
        5,
        4,
        "CH changes by adding generic reals while preserving ordinal scaffolding.",
    ),
    Route(
        "inner-model canonical section",
        5,
        4,
        5,
        4,
        "L gives a CH section with definable well-order side data.",
    ),
    Route(
        "absoluteness audit",
        4,
        5,
        5,
        5,
        "Ask which predicates survive forcing/quotienting unchanged.",
    ),
    Route(
        "equinumerosity-vs-fiber bridge",
        4,
        5,
        4,
        5,
        "Same count is weak; retained fibers make it proof-useful.",
    ),
    Route(
        "ordinal boundary-state transfer",
        4,
        4,
        4,
        5,
        "Open-game ordinal value remembers pole order/boundary fuel.",
    ),
    Route(
        "LRC sufficient-statistic program",
        4,
        5,
        3,
        4,
        "Carry/owner/Cprime labels are the finite analogue of model side data.",
    ),
    Route(
        "raw cardinal numerology",
        1,
        1,
        3,
        1,
        "Bare 2^aleph_0 slogans are too lossy for repo transfer.",
    ),
]


def route_tournament(routes: list[Route]) -> dict[str, object]:
    n = len(routes)
    wins = [0] * n
    edges: dict[tuple[int, int], int] = {}
    for i, j in combinations(range(n), 2):
        a = routes[i].score_tuple
        b = routes[j].score_tuple
        if a > b:
            winner = i
        elif b > a:
            winner = j
        else:
            winner = min(i, j)
        loser = j if winner == i else i
        edges[(winner, loser)] = 1
        wins[winner] += 1

    c3 = 0
    for i, j, k in combinations(range(n), 3):
        eij = (i, j) in edges
        ejk = (j, k) in edges
        eki = (k, i) in edges
        if (eij and ejk and eki) or ((j, i) in edges and (k, j) in edges and (i, k) in edges):
            c3 += 1

    # Count Hamiltonian paths by DP over the small route tournament.
    adj = [[False] * n for _ in range(n)]
    for (u, v) in edges:
        adj[u][v] = True

    def reachable_from(start: int, graph: list[list[bool]]) -> set[int]:
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
        comp = reachable_from(start, adj) & reachable_from(start, radj)
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
    hpaths = sum(dp[(1 << n) - 1])
    order = sorted(range(n), key=lambda i: (-wins[i], i))
    return {
        "top_order": [routes[i].name for i in order],
        "score_hist": {score: wins.count(score) for score in sorted(set(wins))},
        "directed_3cycles": c3,
        "scc_sizes": scc_sizes,
        "hamiltonian_paths": hpaths,
    }


def main() -> None:
    print("S656 Continuum Hypothesis carrier atlas")
    print("=" * 72)
    print()

    print("A. Baseline facts, typed")
    print("  CH: 2^aleph_0 = aleph_1?")
    print("  Godel: if ZFC is consistent, ZFC + CH is consistent via L.")
    print("  Cohen: if ZFC is consistent, ZFC + not-CH is consistent via forcing.")
    print("  Importance: CH made independence, forcing, and model side-channels unavoidable.")
    print("  Repo caution: this is set-theory CH, not Caccetta-Haggkvist CH.")
    print()

    print("B. Finite powerset shadow: why raw 2^n is too coarse")
    print("  n  2^n  n+1  overshoot  intermediate finite sizes")
    for row in finite_powerset_shadow():
        print(
            f"  {row['n']:1d}  {row['pow']:3d}  {row['finite_successor']:3d}"
            f"  {row['overshoot']:9d}  {row['intermediate_sizes']:25d}"
        )
    print("  Finite analogy is not CH; it is a warning that powerset growth is a scalar jump.")
    print()

    print("C. Model/forcing carrier table")
    for row in FOUNDATION_ROWS:
        print(f"  {row['world']}:")
        print(f"    method={row['method']}")
        print(f"    continuum={row['continuum']}  CH={row['CH']}")
        print(f"    retained_channel={row['retained_channel']}")
        print(f"    repo_echo={row['repo_echo']}")
    print()

    print("D. Repo bridge table")
    for row in REPO_BRIDGES:
        print(f"  {row['thread']}:")
        print(f"    scalar_shadow={row['scalar_shadow']}")
        print(f"    missing_middle={row['missing_middle']}")
        print(f"    retained_channel={row['retained_channel']}")
    print()

    print("E. Tournament Analysis over CH transfer routes")
    fp = route_tournament(ROUTES)
    for key, value in fp.items():
        print(f"  {key}={value}")
    print()

    print("Conclusion")
    print("  CH is important here because it is the archetype of a quotient being")
    print("  too weak: the scalar cardinal 2^aleph_0 does not determine the universe")
    print("  of reals without a model/generic side channel.  The repo analogue is to")
    print("  stop treating counts, residues, H-values, or wall scalars as complete")
    print("  proof objects until the predicate-preserving fiber is named.")


if __name__ == "__main__":
    main()
