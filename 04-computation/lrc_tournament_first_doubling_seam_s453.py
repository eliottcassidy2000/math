#!/usr/bin/env python3
"""
lrc_tournament_first_doubling_seam_s453.py

codex-2026-05-31 S453

Probe the first-doubling seam shared by tournament recursion and LRC.

Every integer is n = 2^r * odd.  The row step odd -> 2*odd is not an
ordinary doubling:

* matching pairs gain one extra pair from the odd unmatched vertex;
* unit residues satisfy U(2m) ~= U(m), not two copies of U(m);
* the SC tournament blowup forces a universal near-regular twin score;
* LRC first-even rows open quotient/nonunit room without new unit witnesses.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from itertools import combinations
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
S411 = SourceFileLoader(
    "lrc_column_row_modes_s411",
    str(ROOT / "04-computation" / "lrc_column_row_modes_s411.py"),
).load_module()


def fmt_frac(value: Fraction | None) -> str:
    if value is None:
        return "-"
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def hist(values: list[int] | tuple[int, ...]) -> str:
    parts = []
    for value in sorted(set(values)):
        count = sum(1 for x in values if x == value)
        parts.append(f"{value}^{count}" if count > 1 else str(value))
    return " ".join(parts)


def vp2(value: int) -> int:
    out = 0
    while value and value % 2 == 0:
        out += 1
        value //= 2
    return out


def odd_core(value: int) -> int:
    return value >> vp2(value)


def phi(n: int) -> int:
    return sum(1 for a in range(1, n + 1) if gcd(a, n) == 1)


def units(n: int) -> tuple[int, ...]:
    return tuple(a for a in range(1, n) if gcd(a, n) == 1)


def odd_lift_to_2m(a: int, m: int) -> int:
    """The unique odd lift of a mod m to a unit mod 2m."""
    return a if a % 2 == 1 else a + m


def pairs(n: int) -> int:
    return n // 2


def walsh_degree(n: int) -> int:
    return 2 * ((n - 1) // 2)


def transitive(n: int) -> list[list[bool]]:
    return [[i < j for j in range(n)] for i in range(n)]


def cyclic_regular(n: int) -> list[list[bool]]:
    assert n % 2 == 1
    half = (n - 1) // 2
    adj = [[False] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        if (j - i) % n <= half:
            adj[i][j] = True
        else:
            adj[j][i] = True
    return adj


def sc_blowup(adj: list[list[bool]]) -> list[list[bool]]:
    n = len(adj)
    out = [[False] * (2 * n) for _ in range(2 * n)]

    def idx(v: int, bit: int) -> int:
        return 2 * v + bit

    for v in range(n):
        out[idx(v, 0)][idx(v, 1)] = True
    for u, v in combinations(range(n), 2):
        if adj[u][v]:
            winner, loser = u, v
        else:
            winner, loser = v, u
        out[idx(winner, 0)][idx(loser, 0)] = True
        out[idx(winner, 1)][idx(loser, 1)] = True
        out[idx(loser, 0)][idx(winner, 1)] = True
        out[idx(loser, 1)][idx(winner, 0)] = True
    return out


def lex_blowup(adj: list[list[bool]]) -> list[list[bool]]:
    n = len(adj)
    out = [[False] * (2 * n) for _ in range(2 * n)]

    def idx(v: int, bit: int) -> int:
        return 2 * v + bit

    for v in range(n):
        out[idx(v, 0)][idx(v, 1)] = True
    for u, v in combinations(range(n), 2):
        winner, loser = (u, v) if adj[u][v] else (v, u)
        for a in (0, 1):
            for b in (0, 1):
                out[idx(winner, a)][idx(loser, b)] = True
    return out


def scores(adj: list[list[bool]]) -> tuple[int, ...]:
    return tuple(sum(1 for value in row if value) for row in adj)


def hp_count(adj: list[list[bool]]) -> int:
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp[mask][last]
            if not count:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += count
    return sum(dp[-1])


def all_tournaments(n: int):
    edges = list(combinations(range(n), 2))
    for mask in range(1 << len(edges)):
        adj = [[False] * n for _ in range(n)]
        for bit, (i, j) in enumerate(edges):
            if mask & (1 << bit):
                adj[j][i] = True
            else:
                adj[i][j] = True
        yield adj


def print_arithmetic_seam() -> None:
    print("FIRST DOUBLING SEAM: ARITHMETIC LEDGER")
    print("=" * 96)
    print("odd m -> 2m is a parity-sheet selection, not a normal 2-copy lift.")
    print()
    print(
        "  m  2m  pairs defect  phi defect  later phi defect  "
        "unit density drop  Walsh missing"
    )
    for m in (1, 3, 5, 7, 9, 11, 15, 21):
        child = 2 * m
        pair_defect = pairs(child) - 2 * pairs(m)
        phi_defect = 2 * phi(m) - phi(child)
        later_phi_defect = 2 * phi(child) - phi(2 * child)
        density_drop = Fraction(phi(child), child) / Fraction(phi(m), m)
        walsh_missing = (child - 1) - walsh_degree(child)
        print(
            f"  {m:>2} {child:>3} {pair_defect:>13} {phi_defect:>11} "
            f"{later_phi_defect:>17} {fmt_frac(density_drop):>18} "
            f"{walsh_missing:>14}"
        )
    print()
    print("Critical exact facts")
    print("- pairs(2m) = 2*pairs(m) + 1 for odd m: the unmatched vertex gains a twin.")
    print("- phi(2m) = phi(m) for odd m: the unit boundary does not double.")
    print("- phi(4m) = 2*phi(2m): later row steps are honest doublings.")
    print("- Walsh degree at even 2m misses the top size degree by exactly 1.")


def print_unit_lifts() -> None:
    print()
    print("UNIT LIFTS THROUGH THE FIRST PARITY SHEET")
    print("=" * 96)
    for m in (7, 9, 15):
        lifts = [(a, odd_lift_to_2m(a, m)) for a in units(m)]
        print(f"  m={m:<2} U(m)->U(2m): {lifts}")
    print()
    print("Interpretation: reduction U(2m)->U(m) is bijective at the first seam.")
    print("The second sheet exists arithmetically, but exactly one lift is odd and survives.")


def print_tournament_seam() -> None:
    print()
    print("SC TOURNAMENT BLOWUP: TWIN SCORE COLLAPSE")
    print("=" * 96)
    print("The SC blowup keeps T and T^op data, but forces universal twin scores.")
    print()
    print("  base        n  score(T)        H(T)  score(SC)       H(SC)  H(Lex)")
    examples = [
        ("transitive", 3, transitive(3)),
        ("cyclic", 3, cyclic_regular(3)),
        ("transitive", 5, transitive(5)),
        ("cyclic", 5, cyclic_regular(5)),
    ]
    for name, n, adj in examples:
        sc = sc_blowup(adj)
        lex = lex_blowup(adj)
        print(
            f"  {name:<10} {n:>2}  {hist(scores(adj)):<14} {hp_count(adj):>5} "
            f"{hist(scores(sc)):<15} {hp_count(sc):>7} {hp_count(lex):>7}"
        )

    print()
    print("All labelled tournaments: SC blowup concentrates H while lex blowup magnifies variation.")
    print("  n  H(T) range     H(SC) range       H(Lex) range      SC score")
    for n in (3, 5):
        h_values = []
        sc_values = []
        lex_values = []
        score_shapes = set()
        for adj in all_tournaments(n):
            h_values.append(hp_count(adj))
            sc = sc_blowup(adj)
            lex = lex_blowup(adj)
            sc_values.append(hp_count(sc))
            lex_values.append(hp_count(lex))
            score_shapes.add(hist(scores(sc)))
        print(
            f"  {n:>1}  {min(h_values):>4}..{max(h_values):<6} "
            f"{min(sc_values):>8}..{max(sc_values):<8} "
            f"{min(lex_values):>8}..{max(lex_values):<8} "
            f"{next(iter(score_shapes))}"
        )


def print_lrc_first_even_rows() -> None:
    print()
    print("LRC FIRST-EVEN ROWS AND ROW-PARENT LADDERS")
    print("=" * 96)
    print("For even n, lpd(n)=n/2: the quotient ladder is the row parent.")
    print()
    print("  n   state      lpd  gap/th      debt  product   meaning")
    for n in (6, 10, 14, 16, 18, 22, 26):
        row = S411.best_lpd_ladder(n)
        state = f"2^{vp2(n)}*{odd_core(n)}"
        meaning = "pure dyadic column" if odd_core(n) == 1 else f"first child of {odd_core(n)}"
        if vp2(n) > 1 and odd_core(n) != 1:
            meaning = f"deeper child of {odd_core(n)}"
        print(
            f"  {n:>2}  {state:<9} {row.lpd:>3} {fmt_frac(row.arch_gap):>8} "
            f"{row.debt:>7} {fmt_frac(row.product):>8}  {meaning}"
        )


def print_synthesis() -> None:
    print()
    print("S453 SYNTHESIS")
    print("=" * 96)
    print("The first doubling is a parity seam, not a generic row step.")
    print("Tournament side: an odd unmatched vertex gains a twin, and SC blowup")
    print("turns arbitrary score data into the universal near-regular split n/n-1.")
    print("LRC side: unit endpoints do not double; U(2m) is the odd-lift copy of U(m).")
    print("The new room is nonunit quotient room, so first-even denominators such as")
    print("14 and 18 should be attacked by showing every use of that room exports")
    print("endpoint debt into the row-parent/product-depth ledger.")


def main() -> None:
    print_arithmetic_seam()
    print_unit_lifts()
    print_tournament_seam()
    print_lrc_first_even_rows()
    print_synthesis()


if __name__ == "__main__":
    main()
