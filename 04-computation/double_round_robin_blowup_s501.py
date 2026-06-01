#!/usr/bin/env python3
"""
double_round_robin_blowup_s501.py

codex-2026-06-01 S501c

Formalize the "double round robin" meaning of doubling tournament vertices.

The useful object is not the lexicographic blowup.  It is a 2-fiber lift:
each old vertex v becomes v_0 -> v_1, and every old pair {u,v} becomes a
2x2 block.  A double-round-robin block is one where the base winner wins one
perfect matching of the K_{2,2} block and the base loser wins the complementary
perfect matching.

The matching choice is a voltage bit sigma_uv.  The all-zero voltage is the
SC blowup already used in the repo.  Sheet flips at vertices gauge-transform
sigma; triangle parities are the gauge-invariant content.
"""

from __future__ import annotations

from collections import Counter
from itertools import combinations, product


def idx(v: int, bit: int) -> int:
    return 2 * v + bit


def pairs(n: int) -> list[tuple[int, int]]:
    return list(combinations(range(n), 2))


def transitive(n: int) -> list[list[bool]]:
    return [[i < j for j in range(n)] for i in range(n)]


def cyclic_regular(n: int) -> list[list[bool]]:
    assert n % 2 == 1
    half = (n - 1) // 2
    adj = [[False] * n for _ in range(n)]
    for i, j in pairs(n):
        if (j - i) % n <= half:
            adj[i][j] = True
        else:
            adj[j][i] = True
    return adj


def edge_key(u: int, v: int) -> tuple[int, int]:
    return (u, v) if u < v else (v, u)


def get_sigma(sigma: dict[tuple[int, int], int], u: int, v: int) -> int:
    return sigma[edge_key(u, v)]


def voltage_from_mask(n: int, mask: int) -> dict[tuple[int, int], int]:
    sigma: dict[tuple[int, int], int] = {}
    for bit, edge in enumerate(pairs(n)):
        sigma[edge] = (mask >> bit) & 1
    return sigma


def fundamental_voltage(n: int, mask: int) -> dict[tuple[int, int], int]:
    """Gauge-normalized voltage with all root edges 0 and free edges among 1..n-1."""
    sigma = {edge: 0 for edge in pairs(n)}
    free = list(combinations(range(1, n), 2))
    for bit, edge in enumerate(free):
        sigma[edge] = (mask >> bit) & 1
    return sigma


def double_round_robin_lift(
    base: list[list[bool]], sigma: dict[tuple[int, int], int]
) -> list[list[bool]]:
    """Signed SC/double-round-robin 2-fiber lift."""
    n = len(base)
    out = [[False] * (2 * n) for _ in range(2 * n)]

    for v in range(n):
        out[idx(v, 0)][idx(v, 1)] = True

    for u, v in pairs(n):
        if base[u][v]:
            winner, loser = u, v
        else:
            winner, loser = v, u
        s = get_sigma(sigma, u, v)
        for r in (0, 1):
            out[idx(winner, r)][idx(loser, r ^ s)] = True
            out[idx(loser, r)][idx(winner, r ^ s ^ 1)] = True

    return out


def lex_blowup(base: list[list[bool]]) -> list[list[bool]]:
    n = len(base)
    out = [[False] * (2 * n) for _ in range(2 * n)]
    for v in range(n):
        out[idx(v, 0)][idx(v, 1)] = True
    for u, v in pairs(n):
        winner, loser = (u, v) if base[u][v] else (v, u)
        for a in (0, 1):
            for b in (0, 1):
                out[idx(winner, a)][idx(loser, b)] = True
    return out


def is_tournament(adj: list[list[bool]]) -> bool:
    n = len(adj)
    for i in range(n):
        if adj[i][i]:
            return False
    for i, j in pairs(n):
        if adj[i][j] == adj[j][i]:
            return False
    return True


def scores(adj: list[list[bool]]) -> tuple[int, ...]:
    return tuple(sum(row) for row in adj)


def count_hp(adj: list[list[bool]]) -> int:
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
                if mask & (1 << nxt):
                    continue
                if adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += val
    return sum(dp[-1])


def block_wins(adj: list[list[bool]], u: int, v: int) -> frozenset[tuple[int, int]]:
    """Cells (a,b) where u_a beats v_b."""
    cells = set()
    for a in (0, 1):
        for b in (0, 1):
            if adj[idx(u, a)][idx(v, b)]:
                cells.add((a, b))
    return frozenset(cells)


def is_perfect_matching(cells: frozenset[tuple[int, int]]) -> bool:
    if len(cells) != 2:
        return False
    return (
        sorted(a for a, _ in cells) == [0, 1]
        and sorted(b for _, b in cells) == [0, 1]
    )


def aggregate_block_counts(adj: list[list[bool]], n: int) -> Counter[int]:
    counts = Counter()
    for u, v in pairs(n):
        counts[len(block_wins(adj, u, v))] += 1
    return counts


def quotient_team_wins(adj: list[list[bool]], n: int) -> tuple[int, ...]:
    """Wins by each doubled pair against other doubled pairs; internal arc ignored."""
    wins = [0] * n
    for u, v in pairs(n):
        uv = len(block_wins(adj, u, v))
        wins[u] += uv
        wins[v] += 4 - uv
    return tuple(wins)


def triangle_parity(sigma: dict[tuple[int, int], int], a: int, b: int, c: int) -> int:
    return get_sigma(sigma, a, b) ^ get_sigma(sigma, b, c) ^ get_sigma(sigma, a, c)


def gauge_normalize(n: int, sigma: dict[tuple[int, int], int], root: int = 0) -> tuple[int, ...]:
    """Normalize by sheet flips so all root-edge voltages are zero."""
    tau = [0] * n
    for v in range(n):
        if v != root:
            tau[v] = get_sigma(sigma, root, v)
    normalized = []
    for i, j in pairs(n):
        normalized.append(get_sigma(sigma, i, j) ^ tau[i] ^ tau[j])
    return tuple(normalized)


def fundamental_parities(n: int, sigma: dict[tuple[int, int], int]) -> tuple[int, ...]:
    return tuple(triangle_parity(sigma, 0, i, j) for i, j in combinations(range(1, n), 2))


def print_block_taxonomy() -> None:
    print("2x2 BLOCK TAXONOMY")
    print("=" * 80)
    all_cells = [(0, 0), (0, 1), (1, 0), (1, 1)]
    by_count: Counter[int] = Counter()
    balanced = []
    matching = []
    star = []
    for choices in product((0, 1), repeat=4):
        cells = frozenset(cell for bit, cell in zip(choices, all_cells) if bit)
        by_count[len(cells)] += 1
        if len(cells) == 2:
            balanced.append(cells)
            if is_perfect_matching(cells):
                matching.append(cells)
            else:
                star.append(cells)
    print(f"all orientations of a 2x2 fiber block: {sum(by_count.values())}")
    print(f"by number of wins for the first fiber: {dict(sorted(by_count.items()))}")
    print(f"balanced 2-2 blocks: {len(balanced)}")
    print(f"double-round-robin matching blocks: {len(matching)}")
    print(f"balanced row/column-star blocks, not fiber-schedule balanced: {len(star)}")
    print("matching blocks:")
    for cells in matching:
        print(f"  {sorted(cells)}")
    print()
    print("Lex blowup uses 4-0 or 0-4 blocks; it magnifies the old result.")
    print("SC/voltage blowups use matching 2-2 blocks; they store the old result as a parity.")


def print_voltage_counts() -> None:
    print()
    print("VOLTAGE/GAUGE COUNT")
    print("=" * 80)
    print("  n   edges  voltages  sheet-flip gauge classes  fundamental triangle bits")
    for n in range(2, 8):
        m = n * (n - 1) // 2
        classes = 2 ** (m - n + 1)
        fund = (n - 1) * (n - 2) // 2
        print(f"  {n:>1} {m:>7} {2**m:>9} {classes:>26} {fund:>28}")

    print()
    print("Gauge normalization check")
    for n in range(2, 7):
        normalized = set()
        parity_vectors = set()
        for mask in range(1 << (n * (n - 1) // 2)):
            sigma = voltage_from_mask(n, mask)
            normalized.add(gauge_normalize(n, sigma))
            parity_vectors.add(fundamental_parities(n, sigma))
        target = 2 ** ((n - 1) * (n - 2) // 2)
        print(
            f"  n={n}: normalized={len(normalized)} parity_vectors={len(parity_vectors)} "
            f"target={target}"
        )


def verify_lift_family() -> None:
    print()
    print("LIFT FAMILY VERIFICATION")
    print("=" * 80)
    examples = [
        ("transitive-3", transitive(3)),
        ("cyclic-3", cyclic_regular(3)),
        ("transitive-4", transitive(4)),
        ("transitive-5", transitive(5)),
        ("cyclic-5", cyclic_regular(5)),
    ]
    for name, base in examples:
        n = len(base)
        free_bits = (n - 1) * (n - 2) // 2
        score_hist = Counter()
        team_hist = Counter()
        h_values = []
        all_tourn = True
        all_blocks_matching = True
        for mask in range(1 << free_bits):
            sigma = fundamental_voltage(n, mask)
            lift = double_round_robin_lift(base, sigma)
            all_tourn = all_tourn and is_tournament(lift)
            for u, v in pairs(n):
                all_blocks_matching = all_blocks_matching and is_perfect_matching(
                    block_wins(lift, u, v)
                )
            score_hist[tuple(sorted(scores(lift)))] += 1
            team_hist[quotient_team_wins(lift, n)] += 1
            if 2 * n <= 10:
                h_values.append(count_hp(lift))
        h_summary = "-"
        if h_values:
            h_counter = Counter(h_values)
            h_summary = (
                f"min={min(h_values)} max={max(h_values)} distinct={len(h_counter)} "
                f"top={sorted(h_counter.items())[:5]}"
            )
        print(f"  {name}: gauges={1 << free_bits}")
        print(f"    tournament={all_tourn} matching_blocks={all_blocks_matching}")
        print(f"    score_hist={dict(score_hist)}")
        print(f"    quotient_team_wins={dict(team_hist)}")
        print(f"    H over gauges: {h_summary}")

    print()
    print("Lex contrast")
    for name, base in [("transitive-5", transitive(5)), ("cyclic-5", cyclic_regular(5))]:
        n = len(base)
        lex = lex_blowup(base)
        print(f"  {name}: block win counts={dict(aggregate_block_counts(lex, n))}")
        print(f"    score={tuple(sorted(scores(lex)))}")
        print(f"    quotient_team_wins={quotient_team_wins(lex, n)}")


def print_sc_rule() -> None:
    print()
    print("FORMAL RULE")
    print("=" * 80)
    print("For base arc u -> v and voltage sigma_uv in F2:")
    print("  u_r -> v_{r+sigma_uv}")
    print("  v_r -> u_{r+sigma_uv+1}")
    print("for r=0,1.  Internal arcs are v_0 -> v_1.")
    print()
    print("sigma=0 for every edge gives the repo's SC blowup:")
    print("  same sheet follows T; opposite sheet follows T^op.")
    print()
    print("A sheet flip tau_v at each vertex sends")
    print("  sigma_uv -> sigma_uv + tau_u + tau_v.")
    print("The invariant data are triangle parities")
    print("  p_abc = sigma_ab + sigma_bc + sigma_ac.")


def main() -> None:
    print_block_taxonomy()
    print_voltage_counts()
    print_sc_rule()
    verify_lift_family()
    print()
    print("SYNTHESIS")
    print("=" * 80)
    print("Doubling vertices turns a tournament edge into a fiberwise double round robin.")
    print("The old winner is not stored as a quotient score; every doubled pair ties 2-2.")
    print("It is stored in which perfect matching of the 2x2 block carries the old edge.")
    print("The hidden matching choices form a signed 2-lift; modulo sheet flips there are")
    print("2^C(n-1,2) choices, the same cube dimension as fixed-base tournament tilings.")


if __name__ == "__main__":
    main()
