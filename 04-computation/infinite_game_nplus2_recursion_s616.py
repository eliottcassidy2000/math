#!/usr/bin/env python3
"""
infinite_game_nplus2_recursion_s616.py

S616 scout: use Hamkins-style infinite-game/ordinal-value thinking to sharpen
the repo's n+2 recursion programs.

External prompt context: the linked X post did not expose text to the web
fetcher, but public Hamkins pages verify the relevant background: infinite
games, infinite Go, and finite positions with ordinal values such as omega^2+k.

Mathematical point:
  * S575 says even LRC n has the clean polygon/tournament ladder and n -> n+2
    preserves it.
  * HYP-1881/HYP-1891 say LRC denominator recursion has two natural modes:
    odd-root column motion b -> b+2 and row motion n -> 2n.
  * HYP-2177 says the V* doubling sporadic is tight iff 3 | (2n-1).

The new observation is that the clean n+2 ladder carries a finite automaton:
under n -> n+2, C=2n-1 moves by C -> C+4, so C mod 3 advances by +1.
Therefore the V* family is not a vague infinite family; it is alive exactly on
one state of a period-3 n+2 automaton.  The full row/column recursion has
ordinal proof rank omega^2, just like an infinite-game value hierarchy.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations


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


def factor_str(n: int) -> str:
    return "*".join(
        f"{p}^{e}" if e > 1 else str(p) for p, e in factor(n).items()
    )


def canonical_residue(r: int, modulus: int) -> int:
    r %= modulus
    return min(r, modulus - r)


def doubling_reflection_orbits(modulus: int) -> list[list[int]]:
    """Orbits on antipodal residue shells under <2,-1> for odd modulus."""
    universe = set(range(1, (modulus + 1) // 2))
    orbits: list[list[int]] = []
    while universe:
        start = min(universe)
        orbit = set()
        cur = start
        while cur not in orbit:
            orbit.add(cur)
            cur = canonical_residue(2 * cur, modulus)
        # Multiplication by 2 is invertible for odd modulus, but the reflected
        # canonical path can enter the same cycle from every member.
        closed = set(orbit)
        changed = True
        while changed:
            changed = False
            for x in list(closed):
                y = canonical_residue(2 * x, modulus)
                if y not in closed:
                    closed.add(y)
                    changed = True
        orbits.append(sorted(closed))
        universe -= closed
    return sorted(orbits, key=lambda xs: (min(xs), len(xs)))


@dataclass(frozen=True)
class CleanLadderState:
    n_lrc: int
    runners: int
    c_modulus: int
    c_factor: str
    c_mod_3: int
    vstar_alive: bool
    shell_orbits: int
    orbit_sizes: tuple[int, ...]
    ordinal_token: str


def clean_ladder_state(n_lrc: int) -> CleanLadderState:
    c = 2 * n_lrc - 1
    orbits = doubling_reflection_orbits(c)
    phase = c % 3
    column_index = (n_lrc - 4) // 2
    return CleanLadderState(
        n_lrc=n_lrc,
        runners=n_lrc - 1,
        c_modulus=c,
        c_factor=factor_str(c),
        c_mod_3=phase,
        vstar_alive=(phase == 0),
        shell_orbits=len(orbits),
        orbit_sizes=tuple(len(o) for o in orbits),
        ordinal_token=f"omega*0+{column_index}",
    )


@dataclass(frozen=True)
class DenominatorRank:
    denominator: int
    odd_root: int
    row_depth: int
    column_index: int
    rank_token: str


def denominator_rank(n: int) -> DenominatorRank:
    row = 0
    odd = n
    while odd % 2 == 0:
        row += 1
        odd //= 2
    col = (odd - 1) // 2
    # Lexicographic proof rank: finish rows inside each odd-root column, then
    # move to the next column.  This is an omega^2-shaped recursion skeleton.
    return DenominatorRank(n, odd, row, col, f"omega*{col}+{row}")


@dataclass(frozen=True)
class Route:
    name: str
    preserves_lrc_predicate: int
    infinite_rank_value: int
    finite_automaton_value: int
    side_channel_retention: int
    proof_actionability: int
    risk: int

    @property
    def score(self) -> int:
        return (
            2 * self.preserves_lrc_predicate
            + self.infinite_rank_value
            + self.finite_automaton_value
            + self.side_channel_retention
            + self.proof_actionability
            - 2 * self.risk
        )


ROUTES = (
    Route("omega^2 two-mode proof rank", 5, 5, 2, 5, 5, 1),
    Route("period-3 V* n+2 automaton", 5, 3, 5, 4, 5, 1),
    Route("clean dihedral polygon ladder", 4, 4, 3, 4, 4, 1),
    Route("<2,-1> shell-orbit ledger", 5, 3, 3, 5, 4, 1),
    Route("unit-distance/infinite-game carrier analogy", 3, 5, 2, 4, 3, 2),
    Route("raw next-n brute force", 2, 1, 1, 1, 2, 4),
    Route("unseen X-post text as theorem", 0, 0, 0, 0, 0, 5),
)


def route_tournament(routes: tuple[Route, ...]) -> list[list[int]]:
    n = len(routes)
    adj = [[0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        a, b = routes[i], routes[j]
        key_a = (a.score, -a.risk, a.preserves_lrc_predicate, a.name)
        key_b = (b.score, -b.risk, b.preserves_lrc_predicate, b.name)
        if key_a >= key_b:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return adj


def directed_triangles(adj: list[list[int]]) -> int:
    total = 0
    for a, b, c in combinations(range(len(adj)), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (
            adj[a][c] and adj[c][b] and adj[b][a]
        ):
            total += 1
    return total


def hamiltonian_path_count(adj: list[list[int]]) -> int:
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last in range(n):
            if not dp[mask][last]:
                continue
            for nxt in range(n):
                if (mask >> nxt) & 1:
                    continue
                if adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += dp[mask][last]
    return sum(dp[-1])


def print_clean_ladder() -> None:
    print("Clean even-LRC n+2 ladder")
    print("n  m  C=2n-1  factor(C)  Cmod3  V*_alive  shell_orbits  orbit_sizes  rank")
    states = [clean_ladder_state(n) for n in range(4, 31, 2)]
    for s in states:
        print(
            f"{s.n_lrc:2d} {s.runners:2d} {s.c_modulus:7d} {s.c_factor:9s} "
            f"{s.c_mod_3:5d} {str(s.vstar_alive):>9s} "
            f"{s.shell_orbits:12d} {str(s.orbit_sizes):>16s} {s.ordinal_token}"
        )
    alive = [s.n_lrc for s in states if s.vstar_alive]
    print()
    print(f"V* alive columns through n=30: {alive}")
    print("Automaton law: n -> n+2 sends C -> C+4, hence C mod 3 advances by +1.")
    print("Therefore V* is alive every third clean step: n = 8,14,20,26,...")
    print("Equivalently n == 2 mod 6 on the even ladder.")
    print()


def print_two_mode_rank() -> None:
    print("Two-mode denominator rank, HYP-1881/HYP-1891")
    print("N  odd_root  row_depth  column  ordinal_rank")
    for n in [7, 14, 16, 18, 21, 27, 28, 32, 36, 40, 54]:
        r = denominator_rank(n)
        print(
            f"{r.denominator:2d} {r.odd_root:9d} {r.row_depth:10d} "
            f"{r.column_index:7d} {r.rank_token}"
        )
    print()
    print("Reading: odd-root motion b -> b+2 is the outer omega ladder;")
    print("row motion N -> 2N is the inner finite/omega ladder.  A proof that")
    print("discharges rows inside columns and then moves columns has omega^2 shape.")
    print("This is the precise infinite-game import: use ordinal rank to prevent")
    print("an infinite proof game from cycling through larger finite boards.")
    print()


def print_route_tournament() -> None:
    print("Tournament Analysis over n+2 proof routes")
    adj = route_tournament(ROUTES)
    scores = [sum(row) for row in adj]
    hist: dict[int, int] = {}
    for score in scores:
        hist[score] = hist.get(score, 0) + 1
    for idx, route in sorted(
        enumerate(ROUTES),
        key=lambda item: (scores[item[0]], item[1].score, -item[1].risk),
        reverse=True,
    ):
        print(
            f"score={scores[idx]} route={route.name}; "
            f"features=(pred={route.preserves_lrc_predicate}, "
            f"rank={route.infinite_rank_value}, automaton={route.finite_automaton_value}, "
            f"side={route.side_channel_retention}, action={route.proof_actionability}, "
            f"risk={route.risk})"
        )
    print(f"score histogram: {dict(sorted(hist.items()))}")
    print(f"directed 3-cycles: {directed_triangles(adj)}")
    print(f"Hamiltonian path count: {hamiltonian_path_count(adj)}")
    print()


def main() -> None:
    print("S616 infinite-game n+2 recursion scout")
    print("Source limitation: X status text was inaccessible; public Hamkins pages")
    print("confirm the relevant infinite-games / infinite-Go / ordinal-value context.")
    print()
    print_clean_ladder()
    print_two_mode_rank()
    print_route_tournament()
    print("Math improvements")
    print("1. V* is a period-3 automaton under the clean n+2 ladder, not a")
    print("   miscellaneous infinite sporadic family.")
    print("2. HYP-1881/HYP-1891 become an omega^2 proof-rank program: columns")
    print("   are odd roots, rows are dyadic depth.")
    print("3. The next LRC n=14 transfer should prove a finite automaton theorem")
    print("   for the n+2 column phases plus a rank-decreasing row discharge.")
    print("4. Infinite-game thinking says not to seek one finite invariant only;")
    print("   seek a well-founded rank on proof states and a periodic automaton")
    print("   for the families that remain alive forever.")


if __name__ == "__main__":
    main()
