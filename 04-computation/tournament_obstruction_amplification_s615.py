#!/usr/bin/env python3
"""
tournament_obstruction_amplification_s615.py

S615 scout for the user's proposed bridge:

    forbidden tournament H-values 7/21
    -> LRC 2n-1 resonance at 21/27
    -> unit-distance n^1.014 arithmetic amplification
    -> Collatz two-block prime-2/prime-3 residuals.

This is a synthesis scout, not a proof.  Its main negative result is important:
the public tournament Hamiltonian-path literature has no known "1.014" exponent
matching Sawin's unit-distance exponent.  So the shared object should be treated
as a carrier-amplification probe, not as an established numerical identity.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from math import gcd


def factor(n: int) -> dict[int, int]:
    out: dict[int, int] = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def phi3(x: int) -> int:
    return x * x + x + 1


def q_for_lrc(n: int) -> int:
    return 2 * n - 1


def vstar_tight_on_lattice(n: int) -> bool:
    q = q_for_lrc(n)
    return q % 3 == 0


def dist_mod(r: int, q: int) -> int:
    r %= q
    return min(r, q - r)


def vstar(n: int) -> list[int]:
    return sorted([x for x in range(1, n) if x != n - 2] + [2 * n - 4])


def lattice_cover_radius(n: int) -> int:
    q = q_for_lrc(n)
    speeds = vstar(n)
    return max(min(dist_mod(v * m, q) for v in speeds) for m in range(1, q))


@dataclass(frozen=True)
class Route:
    name: str
    evidence: int
    cross_domain: int
    actionability: int
    side_channels: int
    exponent_relevance: int
    risk: int


ROUTES = (
    Route("OCF forbidden H=7/21 as low-order obstruction", 5, 5, 4, 4, 2, 1),
    Route("LRC gcd(3,2n-1) lattice law plus carry conservativity", 5, 5, 5, 5, 3, 1),
    Route("unit-distance CM/class-field carrier gives n^1.014", 5, 4, 3, 5, 5, 2),
    Route("Collatz two-block Baker gap as residual template", 4, 4, 3, 4, 2, 2),
    Route("Paley/max-H tournament asymptotics as amplifier", 3, 3, 3, 3, 2, 3),
    Route("search for permanent H-gaps beyond 7 and 21", 3, 3, 4, 2, 1, 3),
    Route("raw numeric equality: tournament exponent = 1.014", 1, 3, 2, 1, 5, 5),
    Route("raw quotient forgetting without side channels", 2, 2, 2, 0, 1, 5),
)


def route_votes(a: Route, b: Route) -> tuple[int, int]:
    criteria = [
        (a.evidence > b.evidence, b.evidence > a.evidence),
        (a.cross_domain > b.cross_domain, b.cross_domain > a.cross_domain),
        (a.actionability > b.actionability, b.actionability > a.actionability),
        (a.side_channels > b.side_channels, b.side_channels > a.side_channels),
        (a.exponent_relevance > b.exponent_relevance, b.exponent_relevance > a.exponent_relevance),
        (a.risk < b.risk, b.risk < a.risk),
    ]
    av = sum(1 for x, y in criteria if x and not y)
    bv = sum(1 for x, y in criteria if y and not x)
    return av, bv


def route_tournament(routes: tuple[Route, ...]) -> list[list[int]]:
    n = len(routes)
    adj = [[0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        iv, jv = route_votes(routes[i], routes[j])
        if iv > jv or (iv == jv and i < j):
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return adj


def score_hist(adj: list[list[int]]) -> dict[int, int]:
    hist: dict[int, int] = {}
    for row in adj:
        d = sum(row)
        hist[d] = hist.get(d, 0) + 1
    return dict(sorted(hist.items()))


def directed_triangles(adj: list[list[int]]) -> int:
    total = 0
    for i, j, k in combinations(range(len(adj)), 3):
        if adj[i][j] and adj[j][k] and adj[k][i]:
            total += 1
        if adj[i][k] and adj[k][j] and adj[j][i]:
            total += 1
    return total


def scc_sizes(adj: list[list[int]]) -> list[int]:
    n = len(adj)

    def reach(starts: list[int], forward: bool) -> set[int]:
        seen = set(starts)
        stack = list(starts)
        while stack:
            v = stack.pop()
            for w in range(n):
                edge = adj[v][w] if forward else adj[w][v]
                if edge and w not in seen:
                    seen.add(w)
                    stack.append(w)
        return seen

    left = set(range(n))
    sizes = []
    while left:
        v = next(iter(left))
        comp = reach([v], True) & reach([v], False)
        sizes.append(len(comp))
        left -= comp
    return sorted(sizes, reverse=True)


def hamiltonian_paths(adj: list[list[int]]) -> int:
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last in range(n):
            val = dp[mask][last]
            if not val:
                continue
            for nxt in range(n):
                if not (mask >> nxt) & 1 and adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += val
    return sum(dp[-1])


def print_arithmetic_ledger() -> None:
    print("Forbidden-H and LRC resonance ledger")
    h7 = phi3(2)
    h21 = 3 * h7
    print(f"  Phi_3(2) = 2^2 + 2 + 1 = {h7}")
    print(f"  3 * Phi_3(2) = {h21}")
    print("  Repo theorem status: H(T)=7 and H(T)=21 are permanent tournament gaps;")
    print("  HYP-2180 signal: strong-component multiplicativity explains the 7/21 gap.")
    print("  HYP-2179 signal: robust finite-spectrum/round-LRC propagation is now recorded.")
    print("  S612b LRC signal: tight AP single-swaps split into exceptional")
    print("    reflection sporadics and the infinite mirror-retaining doubling family.")
    print()
    print("  LRC C=2n-1 windows around the current frontier:")
    print(f"  {'n':>3} {'C':>3} {'factor(C)':>14} {'3|C':>5} {'7|C':>5} {'V* lattice':>12} {'cover radius':>12}")
    for n in range(8, 15):
        q = q_for_lrc(n)
        f = "*".join(f"{p}^{e}" if e > 1 else str(p) for p, e in factor(q).items())
        tight = "tight" if n % 2 == 0 and vstar_tight_on_lattice(n) else "loose/na"
        radius = lattice_cover_radius(n) if n % 2 == 0 else "-"
        print(
            f"  {n:>3} {q:>3} {f:>14} {str(q % 3 == 0):>5} "
            f"{str(q % 7 == 0):>5} {tight:>12} {str(radius):>12}"
        )
    print()


def print_exponent_ledger() -> None:
    print("Exponent ledger")
    print("  Unit distance: Sawin gives planar sets with more than n^1.014 unit pairs.")
    print("  Carrier: high-degree algebraic number fields, small discriminant, many small-norm primes.")
    print("  Tournament max-H literature: max H(T) is Theta(n!/2^(n-1));")
    print("    Alon's upper-bound slack is polynomial n^(3/2), and random regular")
    print("    tournaments give constant-factor e asymptotics over the random baseline.")
    print("  Strong-tournament minimum HP literature has exponential base 5^(1/3),")
    print("    again not 1.014.")
    print("  S615 conclusion: the exact 1.014 is not yet a tournament theorem.")
    print("    Treat it as a target carrier-amplification exponent to search for in")
    print("    proof-obligation tournaments or arithmetic tournament carriers.")
    print()


def print_route_tournament() -> None:
    print("Tournament Analysis on transfer routes")
    adj = route_tournament(ROUTES)
    scores = [sum(row) for row in adj]
    for route, score in sorted(zip(ROUTES, scores), key=lambda x: (-x[1], x[0].name)):
        print(
            f"  score={score} route={route.name}; "
            f"features=(evidence={route.evidence}, cross={route.cross_domain}, "
            f"action={route.actionability}, side={route.side_channels}, "
            f"exp={route.exponent_relevance}, risk={route.risk})"
        )
    print(f"  score histogram: {score_hist(adj)}")
    print(f"  directed 3-cycles: {directed_triangles(adj)}")
    print(f"  SCC sizes: {scc_sizes(adj)}")
    print(f"  Hamiltonian path count: {hamiltonian_paths(adj)}")
    print(
        "  Reading: the winning route is not raw numeric matching.  It is the "
        "side-channel carrier route: OCF/H-gaps, LRC gcd shells, unit-distance "
        "arithmetic towers, and Collatz two-block gaps are all retained-carrier "
        "amplification or obstruction problems."
    )
    print()


def main() -> None:
    print("==== S615 tournament obstruction / amplification scout ====")
    print_arithmetic_ledger()
    print_exponent_ledger()
    print_route_tournament()
    print("Hypothesis HYP-2181:")
    print("  H=7/21 are the visible low-order forbidden faces of the tournament OCF.")
    print("  LRC C=21/27 and the S612b swap split propagate it into shell arithmetic.")
    print("  Unit-distance n^1.014 is the current arithmetic-carrier amplification model.")
    print("  Collatz supplies the two-block gap template.")
    print("  The exact shared exponent is open; the shared structure is the retained carrier.")


if __name__ == "__main__":
    main()
