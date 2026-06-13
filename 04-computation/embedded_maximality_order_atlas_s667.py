#!/usr/bin/env python3
"""S667: embedded maximality atlas.

Prompt:

    merge in the concept of embedded maximality.
    think of the rational numbers with the ordering less than.
    extend beyond that.

The working idea is that "maximal" is not an intrinsic scalar predicate.  It is
relative to an ambient embedding and a class of allowed extensions:

    maximal(object, ambient, allowed_extensions).

The rationals with < are the clean toy.  A finite embedded chain has a maximum
inside itself, but the ambient dense order realizes a point above it and points
in every open cut.  Thus the finite maximum is not ambient maximal.  Dedekind
cuts add the missing address coordinate.

Tournament Analysis:
  Vertices are embedded-maximality lanes, not runners or arcs.  Pairwise
  observable is (exact_evidence, ambient_extendability, address_need,
  derivative_leverage, LRC_transfer, actionability).  The switch is majority,
  ties use the curated order as Hamiltonian path.

Assumption challenge:
  Candidate vertices considered included rationals, cuts, finite chains,
  runners, endpoint cores, owner obligations, graph cards, point sets, bases,
  forcing models, and proof obligations.  The script uses lanes because the
  preserved predicate is "which ambient address makes a local maximum stable,
  and which extension destroys it?"
"""

from __future__ import annotations

from collections import Counter, deque
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations


@dataclass(frozen=True)
class Lane:
    name: str
    family: str
    local_maximum: str
    ambient_embedding: str
    missing_address: str
    destructive_extension: str
    derivative_readout: str
    exact_evidence: int
    ambient_extendability: int
    address_need: int
    derivative_leverage: int
    lrc_transfer: int
    actionability: int

    @property
    def vector(self) -> tuple[int, int, int, int, int, int]:
        return (
            self.exact_evidence,
            self.ambient_extendability,
            self.address_need,
            self.derivative_leverage,
            self.lrc_transfer,
            self.actionability,
        )

    @property
    def total(self) -> int:
        return sum(self.vector)


LANES = [
    Lane(
        "Q,< finite-chain maximum",
        "order toy",
        "a finite embedded chain has a last element",
        "countable dense linear order without endpoints",
        "which open cut / endpoint boundary the chain occupies",
        "insert a rational above the finite maximum or inside any gap",
        "extension-count ledger over n+1 cuts",
        5,
        5,
        5,
        4,
        4,
        5,
    ),
    Lane(
        "Dedekind cut completion",
        "order toy",
        "finite denominator search has a best lower approximation",
        "Q embedded in its completion R",
        "cut address: lower set with no rational maximum",
        "raise the best approximation inside the same irrational cut",
        "cut-improvement sequence / endpoint-loss sum",
        5,
        5,
        5,
        4,
        4,
        5,
    ),
    Lane(
        "LRC14 owner-private stability",
        "LRC exact",
        "AP/Vstar/2AP are floor atoms in the visible Res_27 shell",
        "local +27 carry space over the same visible shell",
        "paired owner-private deletion flag",
        "nonzero local carry preserving the visible shell",
        "owner-deletion loss ledger; S666 mixed fibers drop to zero",
        5,
        5,
        5,
        5,
        5,
        5,
    ),
    Lane(
        "Endpoint protection core",
        "LRC endpoint",
        "an all-protected endpoint set is a maximal obstruction core",
        "endpoint incidence graph inside the quotient denominator",
        "protector owner for each endpoint",
        "peel an endpoint with no remaining private protector",
        "core-peeling derivative / exposed endpoint sum",
        4,
        5,
        5,
        5,
        5,
        4,
    ),
    Lane(
        "Tournament endpoint private child",
        "tournaments",
        "private endpoint child is a pivot-maximal witness",
        "endpoint-extension matrix n -> n+1",
        "parent owner row and child column address",
        "complement merging creates non-private SC collision columns",
        "rank/private-pivot deletion ledger",
        4,
        5,
        5,
        5,
        4,
        4,
    ),
    Lane(
        "Graph reconstruction rooted deck",
        "graphs",
        "deck-derived scalar/subgraph counts look maximal",
        "graph embedded in its rooted/deleted-card deck",
        "deleted vertex / card-to-boundary address",
        "unrooted deck collision with different rooted assignment",
        "rooted Kelly/Kocay constraint sums",
        4,
        4,
        5,
        5,
        3,
        4,
    ),
    Lane(
        "Unit-distance frontier owner",
        "unit distance",
        "edge/direction support can be extremal in a fixed point set",
        "point set embedded in extension/deletion frontier",
        "ear/frontier owner and unit-spine endpoint address",
        "add/delete a point while preserving coarse direction support",
        "frontier edge-loss and spine-loss sums",
        3,
        5,
        4,
        4,
        4,
        4,
    ),
    Lane(
        "A000568 marked observer fiber",
        "tournaments/LRC",
        "unmarked isomorphism class is maximal as a quotient object",
        "class embedded in rooted/observer perspectives",
        "observer root and endpoint threshold address",
        "change root while preserving unmarked class",
        "observer-source deletion sums",
        4,
        4,
        5,
        4,
        4,
        4,
    ),
    Lane(
        "Rado/Fraisse extension property",
        "homogeneous structures",
        "finite pattern is maximal inside its induced substructure",
        "universal homogeneous limit",
        "one-point extension type",
        "realize every compatible extension type",
        "extension-type count / age derivative",
        3,
        5,
        5,
        3,
        3,
        3,
    ),
    Lane(
        "Matroid basis activity",
        "matroids",
        "a basis is maximal independent inside the matroid",
        "basis embedded in an ordered matroid / activity expansion",
        "edge/order activity address",
        "exchange an active element without changing coarse rank",
        "deletion-contraction activity sum",
        3,
        4,
        5,
        5,
        3,
        4,
    ),
    Lane(
        "Forcing model maximality",
        "logic",
        "a cardinal or truth value is maximal inside one model",
        "model embedded in forcing/generic extension",
        "generic filter / forcing-name address",
        "force an extension where the scalar statement changes",
        "truth-value ledger over model boundary",
        3,
        5,
        5,
        3,
        3,
        3,
    ),
]


def dense_order_lab() -> dict[str, object]:
    base = [Fraction(0), Fraction(1), Fraction(2)]
    cut_witnesses = [
        Fraction(-1),
        Fraction(1, 2),
        Fraction(3, 2),
        Fraction(3),
    ]
    cuts = [
        "(-infty,0)",
        "(0,1)",
        "(1,2)",
        "(2,infty)",
    ]
    return {
        "base_chain": base,
        "finite_maximum": max(base),
        "one_point_extension_cuts": list(zip(cuts, cut_witnesses)),
        "ambient_maximum_after_extension": max(base + [Fraction(3)]),
        "lesson": "finite maximality is destroyed by the ambient DLO extension property",
    }


def sqrt2_cut_lab(base_den_limit: int = 12, search_den_limit: int = 200) -> dict[str, object]:
    lower = []
    for q in range(1, base_den_limit + 1):
        for p in range(0, 2 * q + 1):
            r = Fraction(p, q)
            if r * r < 2:
                lower.append(r)
    finite_best = max(lower)
    better = None
    for q in range(base_den_limit + 1, search_den_limit + 1):
        for p in range(0, 2 * q + 1):
            r = Fraction(p, q)
            if finite_best < r and r * r < 2:
                better = r
                break
        if better is not None:
            break
    return {
        "cut": "{q in Q : q^2 < 2}",
        "finite_denominator_limit": base_den_limit,
        "finite_best_lower": finite_best,
        "strictly_better_lower_with_larger_denominator": better,
        "finite_best_square": finite_best * finite_best,
        "better_square": better * better if better else None,
        "lesson": "finite search has a maximum; the embedded Dedekind cut has no rational maximum",
    }


def ambient_boundary_lab() -> dict[str, object]:
    point = Fraction(1)
    subset = [Fraction(0), Fraction(1, 2), point]
    ambients = {
        "finite subset {0,1/2,1}": point == max(subset),
        "closed interval [0,1]_Q": True,
        "open interval (0,1)_Q": False,
        "whole Q": False,
    }
    return {
        "point": point,
        "ambient_maximality": ambients,
        "lesson": "the same rational is maximal or not depending on the embedding boundary",
    }


def s666_projection_lab() -> dict[str, object]:
    return {
        "visible_shadow": {"groups": 3, "mixed": 3, "max_bucket": 378},
        "visible+cheap_pair": {"groups": 112, "mixed": 3, "max_bucket": 106},
        "visible+owner_cover_count": {"groups": 976, "mixed": 2, "max_bucket": 2},
        "visible+owner_private_flag": {"groups": 1067, "mixed": 0, "max_bucket": 2},
        "visible+owner_private_count": {"groups": 1134, "mixed": 0, "max_bucket": 1},
        "lesson": "embedded maximality appears when the visible maximum is tested in the owner-deletion ambient",
    }


def tournament(lanes: list[Lane]) -> dict[str, object]:
    n = len(lanes)
    adj = [[0] * n for _ in range(n)]
    out = [0] * n
    for i, j in combinations(range(n), 2):
        vi, vj = lanes[i].vector, lanes[j].vector
        wi = sum(a > b for a, b in zip(vi, vj))
        wj = sum(a < b for a, b in zip(vi, vj))
        winner, loser = (i, j) if wi > wj or (wi == wj and i < j) else (j, i)
        adj[winner][loser] = 1
        out[winner] += 1

    c3 = 0
    for a, b, c in combinations(range(n), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (
            adj[a][c] and adj[c][b] and adj[b][a]
        ):
            c3 += 1

    radj = [[i for i in range(n) if adj[i][j]] for j in range(n)]
    seen = [False] * n
    order: list[int] = []
    for start in range(n):
        if seen[start]:
            continue
        stack = [(start, False)]
        while stack:
            v, done = stack.pop()
            if done:
                order.append(v)
                continue
            if seen[v]:
                continue
            seen[v] = True
            stack.append((v, True))
            for w in range(n):
                if adj[v][w] and not seen[w]:
                    stack.append((w, False))

    seen = [False] * n
    sccs: list[list[str]] = []
    for start in reversed(order):
        if seen[start]:
            continue
        comp = []
        q = deque([start])
        seen[start] = True
        while q:
            v = q.popleft()
            comp.append(lanes[v].name)
            for w in radj[v]:
                if not seen[w]:
                    seen[w] = True
                    q.append(w)
        sccs.append(comp)

    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            val = dp[mask][last]
            if not val:
                continue
            for nxt in range(n):
                if not (mask & (1 << nxt)) and adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += val

    return {
        "score_hist": dict(sorted(Counter(out).items())),
        "directed_3cycles": c3,
        "scc_sizes": sorted([len(c) for c in sccs], reverse=True),
        "hamiltonian_paths": sum(dp[-1]),
        "top_order": [lanes[i].name for i in sorted(range(n), key=lambda k: (-out[k], -lanes[k].total, k))],
    }


def main() -> None:
    print("=" * 78)
    print("S667 embedded maximality order atlas")
    print("=" * 78)
    print()
    print("Principle")
    print("  maximality is not unary: maximal(object, ambient embedding, allowed extensions)")
    print("  Q,< is the toy: finite chains have maxima; dense ambient extensions destroy them.")
    print()

    print("A. Dense-order finite-chain lab")
    dlo = dense_order_lab()
    for key, value in dlo.items():
        print(f"  {key}={value}")
    print()

    print("B. Dedekind-cut address lab")
    cut = sqrt2_cut_lab()
    for key, value in cut.items():
        print(f"  {key}={value}")
    print()

    print("C. Same point, different ambient boundary")
    amb = ambient_boundary_lab()
    for key, value in amb.items():
        print(f"  {key}={value}")
    print()

    print("D. S666 embedded repair readout")
    for key, value in s666_projection_lab().items():
        print(f"  {key}={value}")
    print()

    print("E. Embedded-maximality lanes")
    print(f"{'lane':<36} {'family':<22} {'total':>5} vector")
    for lane in sorted(LANES, key=lambda lane: (-lane.total, lane.name)):
        print(f"{lane.name:<36} {lane.family:<22} {lane.total:5d} {lane.vector}")
    print()

    print("F. Lane details")
    for lane in LANES:
        print(f"  {lane.name}")
        print(f"    local maximum:       {lane.local_maximum}")
        print(f"    ambient embedding:   {lane.ambient_embedding}")
        print(f"    missing address:     {lane.missing_address}")
        print(f"    destructive extension:{lane.destructive_extension}")
        print(f"    derivative readout:  {lane.derivative_readout}")
    print()

    fp = tournament(LANES)
    print("G. Tournament Analysis")
    print("  vertices=embedded-maximality lanes")
    print("  observable=(exact, ambient extension, address need, derivative, LRC transfer, action)")
    print(f"  score_hist={fp['score_hist']}")
    print(f"  directed_3cycles={fp['directed_3cycles']}")
    print(f"  scc_sizes={fp['scc_sizes']}")
    print(f"  hamiltonian_paths={fp['hamiltonian_paths']}")
    print("  top_order:")
    for name in fp["top_order"]:
        print(f"    {name}")
    print()

    print("H. Synthesis")
    print("  Embedded maximality turns 'maximal' into a boundary-sensitive predicate.")
    print("  Q,< says finite maxima are real inside a finite image but unstable in the")
    print("  dense ambient order; Dedekind completion names the missing cut address.")
    print("  For LRC14, S666 suggests the analogous ambient is the owner-deletion")
    print("  ledger: visible Res_27 maxima remain stable only when their private-owner")
    print("  addresses survive every allowed carry extension.")


if __name__ == "__main__":
    main()
