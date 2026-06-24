#!/usr/bin/env python3
"""S165: Haar product discrepancy / tournament tiling synthesis.

This is a proof-interface scout, not a proof.  The user asked to synthesize
recent LRC14 agent work with discrepancy theory and the two-dimensional Haar
product rule.

The executable point is tiny but load-bearing:

    h_I(x) h_J(y) = h_{I x J}(x,y)

On the two children of each interval this is the 2-by-2 sign matrix

    [[+1, -1],
     [-1, +1]].

That same matrix is the elementary fixed-margin switch for 2-by-2 contingency
tables, the checkerboard atom behind mixed discrepancy, and the square tile
hidden inside the tournament tiling model.  Row/column margins forget its sign;
the mixed Haar coefficient remembers it.

Tournament Analysis is run on proof carriers rather than runners.  The
pairwise observable is retention of:

    LRC predicate, row/column margins, mixed Haar sign, fixed-margin switch,
    color resonance, Fejer/dual compatibility, tope/cocircuit labels,
    quotient guardrail.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from itertools import combinations, permutations
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
RESULT = REPO / "05-knowledge" / "results" / "lrc14_haar_product_discrepancy_tiling_codex_s165.out"


HAAR_X = (1, -1)
HAAR_Y = (1, -1)
HAAR_PRODUCT = tuple(tuple(x * y for y in HAAR_Y) for x in HAAR_X)
SWITCH_VECTOR = ((1, -1), (-1, 1))


TABLES = {
    "diagonal_packet": ((1, 0), (0, 1)),
    "anti_diagonal_packet": ((0, 1), (1, 0)),
    "uniform_packet": ((1, 1), (1, 1)),
    "left_column_packet": ((1, 0), (1, 0)),
    "top_row_packet": ((1, 1), (0, 0)),
}


def row_sums(table: tuple[tuple[int, int], tuple[int, int]]) -> tuple[int, int]:
    return tuple(sum(row) for row in table)


def col_sums(table: tuple[tuple[int, int], tuple[int, int]]) -> tuple[int, int]:
    return tuple(table[0][j] + table[1][j] for j in range(2))


def total(table: tuple[tuple[int, int], tuple[int, int]]) -> int:
    return sum(sum(row) for row in table)


def haar_coeff(table: tuple[tuple[int, int], tuple[int, int]]) -> int:
    return sum(table[i][j] * HAAR_PRODUCT[i][j] for i in range(2) for j in range(2))


def add_tables(
    a: tuple[tuple[int, int], tuple[int, int]],
    b: tuple[tuple[int, int], tuple[int, int]],
) -> tuple[tuple[int, int], tuple[int, int]]:
    return tuple(tuple(a[i][j] + b[i][j] for j in range(2)) for i in range(2))  # type: ignore[return-value]


def scale_table(k: int, table: tuple[tuple[int, int], tuple[int, int]]) -> tuple[tuple[int, int], tuple[int, int]]:
    return tuple(tuple(k * table[i][j] for j in range(2)) for i in range(2))  # type: ignore[return-value]


def all_tables_with_total(n: int) -> list[tuple[tuple[int, int], tuple[int, int]]]:
    out = []
    for a in range(n + 1):
        for b in range(n - a + 1):
            for c in range(n - a - b + 1):
                d = n - a - b - c
                out.append(((a, b), (c, d)))
    return out


@dataclass(frozen=True)
class Carrier:
    name: str
    vector: tuple[int, ...]
    note: str


CARRIERS = [
    Carrier(
        "labelled_haar_square_packet",
        (5, 5, 5, 5, 5, 4, 5, 5),
        "keeps margins, mixed sign, switch direction, resonance, and wall labels",
    ),
    Carrier(
        "fixed_margin_switch_cocircuit",
        (5, 5, 5, 5, 3, 3, 5, 5),
        "the 2x2 Markov move [[1,-1],[-1,1]] with wall-crossing orientation",
    ),
    Carrier(
        "mixed_haar_discrepancy",
        (4, 4, 5, 3, 4, 5, 4, 4),
        "detects checkerboard imbalance after row/column margins vanish",
    ),
    Carrier(
        "colored_resonance_congruence",
        (4, 3, 5, 3, 5, 3, 3, 4),
        "HYP-2595 compatible color sum condition for surviving Fourier/Haar modes",
    ),
    Carrier(
        "tope_cocircuit_wall_label",
        (4, 3, 4, 4, 3, 3, 5, 4),
        "HYP-2986 open-tope versus boundary-cocircuit distinction",
    ),
    Carrier(
        "row_column_margin_shadow",
        (3, 5, 1, 1, 2, 2, 2, 2),
        "fixed margins only; cannot distinguish diagonal from anti-diagonal",
    ),
    Carrier(
        "raw_component_count_K",
        (2, 1, 1, 1, 1, 1, 1, 1),
        "HYP-2594 continuous component count; safe but very lossy in tails",
    ),
]

TIE_PATH = [carrier.name for carrier in CARRIERS]


def beats(a: Carrier, b: Carrier) -> bool:
    aw = sum(1 for x, y in zip(a.vector, b.vector) if x > y)
    bw = sum(1 for x, y in zip(a.vector, b.vector) if x < y)
    if aw != bw:
        return aw > bw
    return TIE_PATH.index(a.name) < TIE_PATH.index(b.name)


def adjacency() -> dict[str, set[str]]:
    adj = {carrier.name: set() for carrier in CARRIERS}
    for a, b in combinations(CARRIERS, 2):
        if beats(a, b):
            adj[a.name].add(b.name)
        else:
            adj[b.name].add(a.name)
    return adj


def directed_3cycles(adj: dict[str, set[str]]) -> int:
    count = 0
    names = list(adj)
    for a, b, c in combinations(names, 3):
        edges = {
            (a, b): b in adj[a],
            (b, a): a in adj[b],
            (a, c): c in adj[a],
            (c, a): a in adj[c],
            (b, c): c in adj[b],
            (c, b): b in adj[c],
        }
        if (
            edges[(a, b)] and edges[(b, c)] and edges[(c, a)]
        ) or (
            edges[(a, c)] and edges[(c, b)] and edges[(b, a)]
        ):
            count += 1
    return count


def score_hist(adj: dict[str, set[str]]) -> dict[int, int]:
    return dict(sorted(Counter(len(v) for v in adj.values()).items()))


def hamiltonian_paths(adj: dict[str, set[str]]) -> list[tuple[str, ...]]:
    names = list(adj)
    paths = []
    for order in permutations(names):
        if all(order[i + 1] in adj[order[i]] for i in range(len(order) - 1)):
            paths.append(order)
    return paths


def scc_sizes(adj: dict[str, set[str]]) -> list[int]:
    names = list(adj)
    radj = {name: set() for name in names}
    for a, outs in adj.items():
        for b in outs:
            radj[b].add(a)

    def reach(graph: dict[str, set[str]], start: str) -> set[str]:
        seen = {start}
        q = deque([start])
        while q:
            cur = q.popleft()
            for nxt in graph[cur]:
                if nxt not in seen:
                    seen.add(nxt)
                    q.append(nxt)
        return seen

    remaining = set(names)
    sizes = []
    while remaining:
        start = next(iter(remaining))
        comp = reach(adj, start) & reach(radj, start)
        sizes.append(len(comp))
        remaining -= comp
    return sorted(sizes, reverse=True)


def table_text(table: tuple[tuple[int, int], tuple[int, int]]) -> str:
    return f"[[{table[0][0]},{table[0][1]}],[{table[1][0]},{table[1][1]}]]"


def main() -> None:
    lines: list[str] = []
    lines.append("S165 Haar product discrepancy / tournament tiling synthesis")
    lines.append("=" * 72)
    lines.append("")

    lines.append("2D Haar product atom")
    lines.append(f"h_x={HAAR_X} h_y={HAAR_Y}")
    lines.append(f"h_x tensor h_y = {table_text(HAAR_PRODUCT)}")
    lines.append(f"fixed-margin switch vector = {table_text(SWITCH_VECTOR)}")
    lines.append(
        "identity: the Haar product sign matrix is exactly the 2x2 fixed-margin switch"
    )
    lines.append("")

    lines.append("Packet table readout")
    for name, table in TABLES.items():
        lines.append(
            f"{name:22s} table={table_text(table):13s} rows={row_sums(table)} "
            f"cols={col_sums(table)} total={total(table)} Hxy={haar_coeff(table)}"
        )
    lines.append("")

    groups: dict[tuple[tuple[int, int], tuple[int, int], int], list[tuple[str, int]]] = defaultdict(list)
    for name, table in TABLES.items():
        groups[(row_sums(table), col_sums(table), total(table))].append((name, haar_coeff(table)))

    lines.append("Quotient collision under row/column margins")
    for key, vals in groups.items():
        coeffs = sorted({v for _, v in vals})
        if len(vals) > 1 and len(coeffs) > 1:
            lines.append(f"margin_key={key} has mixed Hxy values: {vals}")
    lines.append(
        "minimal collision: diagonal_packet and anti_diagonal_packet have identical margins "
        "but opposite mixed Haar coefficients."
    )
    lines.append("")

    diag = TABLES["diagonal_packet"]
    anti = TABLES["anti_diagonal_packet"]
    switched = add_tables(anti, SWITCH_VECTOR)
    lines.append("Minimal wall-crossing / Markov switch")
    lines.append(f"anti_diagonal + switch = {table_text(switched)}")
    lines.append(f"equals diagonal: {switched == diag}")
    lines.append(f"Hxy jump = {haar_coeff(diag) - haar_coeff(anti)}")
    lines.append(
        "This is the smallest fixed-margin move; in the tope language it is a wall "
        "crossing, and in the tiling language it is the oriented square tile."
    )
    lines.append("")

    lines.append("Exhaustive 2x2 fixed-margin stress for totals <= 6")
    for n in range(2, 7):
        collision_count = 0
        opposite_count = 0
        for table in all_tables_with_total(n):
            key = (row_sums(table), col_sums(table), total(table))
            mates = [
                other
                for other in all_tables_with_total(n)
                if (row_sums(other), col_sums(other), total(other)) == key
                and other != table
            ]
            if mates:
                collision_count += 1
                if any(haar_coeff(other) != haar_coeff(table) for other in mates):
                    opposite_count += 1
        lines.append(
            f"total={n}: tables_with_margin_mates={collision_count}, "
            f"mates_changing_Hxy={opposite_count}"
        )
    lines.append(
        "Conclusion: margins are useful fixed fibers, but the mixed product sign is "
        "a separate quotient guardrail even in the smallest model."
    )
    lines.append("")

    adj = adjacency()
    paths = hamiltonian_paths(adj)
    lines.append("Tournament Analysis on proof carriers")
    lines.append(
        "pairwise observable: retention of LRC predicate, margins, mixed Haar sign, "
        "fixed-margin switch, color resonance, Fejer compatibility, tope/cocircuit "
        "labels, quotient guardrail"
    )
    lines.append(f"score_hist={score_hist(adj)}")
    lines.append(f"directed_3cycles={directed_3cycles(adj)}")
    lines.append(f"SCC_sizes={scc_sizes(adj)}")
    lines.append(f"hamiltonian_paths={len(paths)}")
    if paths:
        lines.append("canonical_path=" + " > ".join(paths[0]))
    lines.append("")

    lines.append("Synthesis")
    lines.append(
        "HYP-2594 counted continuous components K. HYP-2595 showed only color-compatible "
        "resonances survive. The Haar product square explains why: row/column boundary "
        "data are margins, while the surviving discrepancy is the mixed product coefficient."
    )
    lines.append(
        "HYP-2986's tope/cocircuit packets are the geometric version of the same switch: "
        "an open tope is a positive cell, AP/GW are boundary atoms, and a true residual "
        "would have to survive every mixed-product switch without exposing a witness."
    )
    lines.append(
        "Tournament tiling is therefore not a metaphor bolted on afterward. The fixed "
        "Hamiltonian path chooses the row/column observer; the Haar product tile records "
        "the orientation sign that the row/column quotient would otherwise forget."
    )

    text = "\n".join(lines) + "\n"
    print(text, end="")


if __name__ == "__main__":
    main()
