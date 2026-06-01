#!/usr/bin/env python3
"""
Symbolic-dynamics probes for the Lonely Runner Conjecture.

codex-2026-06-01-S541

The LRC clock for a fixed primitive speed set is coded as a periodic word.
Open chambers between section-boundary events carry one of four observer
danger symbols:

  G = both observer-adjacent danger sections empty;
  L = only the left observer-adjacent section occupied;
  R = only the right observer-adjacent section occupied;
  B = both observer-adjacent sections occupied.

Boundary events carry W if the closed LRC inequality holds exactly at the wall,
and "." otherwise.  Thus the compactified target alphabet is {G, W}.  LRC says
every primitive speed set has a compactified target in this periodic word.

Tournament Analysis declaration
-------------------------------
Vertices:
  symbolic states/return-word letters, not runners.

Pairwise observable:
  for two non-target symbols a,b, compare which symbol appears first more often
  in target-to-target return words.

Switch/gauge:
  orient a -> b when a wins the first-occurrence comparison; ties use the fixed
  symbol order "." < B < L < R.

Tie Hamiltonian path:
  the fixed symbol order above.

Fingerprints:
  block-complexity p(L), target gap, return-word counts, bad subshift SCCs,
  return-order score histograms, directed 3-cycles, SCCs, and Hamiltonian-path
  counts.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction
from functools import reduce
from itertools import combinations, permutations
from math import gcd, isinf


ZERO = Fraction(0)
ONE = Fraction(1)
TARGETS = {"G", "W"}
TIE_ORDER = {".": 0, "B": 1, "L": 2, "R": 3}


def frac(x: Fraction) -> Fraction:
    return x - Fraction(x.numerator // x.denominator)


def dist_to_integer(x: Fraction) -> Fraction:
    y = frac(x)
    return min(y, ONE - y)


def primitive_speed_sets(n: int, max_speed: int):
    for speeds in combinations(range(1, max_speed + 1), n - 1):
        if reduce(gcd, speeds) == 1:
            yield speeds


def section_walls(speeds: tuple[int, ...], n: int) -> list[Fraction]:
    walls = {ZERO}
    for v in speeds:
        for k in range(n * v):
            walls.add(Fraction(k, n * v))
    return sorted(walls)


def occupancy(speeds: tuple[int, ...], n: int, t: Fraction) -> tuple[int, ...]:
    counts = [0] * n
    for v in speeds:
        counts[int(n * frac(Fraction(v) * t)) % n] += 1
    return tuple(counts)


def danger_symbol(speeds: tuple[int, ...], n: int, t: Fraction) -> str:
    counts = occupancy(speeds, n, t)
    left = counts[-1] > 0
    right = counts[0] > 0
    if not left and not right:
        return "G"
    if left and right:
        return "B"
    return "L" if left else "R"


def wall_symbol(speeds: tuple[int, ...], n: int, t: Fraction) -> str:
    threshold = Fraction(1, n)
    if all(dist_to_integer(Fraction(v) * t) >= threshold for v in speeds):
        return "W"
    return "."


def compact_word(speeds: tuple[int, ...], n: int) -> tuple[str, ...]:
    """Return the wall/chamber periodic word over '.', W, G, L, R, B."""
    walls = section_walls(speeds, n)
    word: list[str] = []
    for i, wall in enumerate(walls):
        next_wall = walls[i + 1] if i + 1 < len(walls) else ONE
        if next_wall <= wall:
            continue
        word.append(wall_symbol(speeds, n, wall))
        word.append(danger_symbol(speeds, n, (wall + next_wall) / 2))
    return tuple(word)


def cyclic_blocks(word: tuple[str, ...], length: int) -> set[tuple[str, ...]]:
    if not word:
        return set()
    m = len(word)
    return {tuple(word[(i + j) % m] for j in range(length)) for i in range(m)}


def target_positions(word: tuple[str, ...]) -> list[int]:
    return [i for i, symbol in enumerate(word) if symbol in TARGETS]


def longest_target_gap(word: tuple[str, ...]) -> int | None:
    positions = target_positions(word)
    if not positions:
        return None
    m = len(word)
    gaps = []
    for a, b in zip(positions, positions[1:] + [positions[0] + m]):
        gaps.append(b - a - 1)
    return max(gaps)


def return_words(word: tuple[str, ...]) -> list[tuple[str, ...]]:
    positions = target_positions(word)
    if not positions:
        return []
    m = len(word)
    out = []
    for a, b in zip(positions, positions[1:] + [positions[0] + m]):
        out.append(tuple(word[(a + j) % m] for j in range(1, b - a)))
    return out


def transition_graph(word: tuple[str, ...], exclude_targets: bool = False):
    graph: dict[str, set[str]] = defaultdict(set)
    if not word:
        return graph
    m = len(word)
    for i, a in enumerate(word):
        b = word[(i + 1) % m]
        if exclude_targets and (a in TARGETS or b in TARGETS):
            continue
        graph[a].add(b)
        graph.setdefault(b, set())
    return graph


def strongly_connected_components(graph: dict[str, set[str]]) -> list[set[str]]:
    nodes = set(graph)
    for outs in graph.values():
        nodes.update(outs)

    order: list[str] = []
    seen: set[str] = set()

    def dfs(v: str):
        seen.add(v)
        for w in graph.get(v, ()):
            if w not in seen:
                dfs(w)
        order.append(v)

    for v in nodes:
        if v not in seen:
            dfs(v)

    rev: dict[str, set[str]] = defaultdict(set)
    for v in nodes:
        rev.setdefault(v, set())
    for v, outs in graph.items():
        for w in outs:
            rev[w].add(v)

    comps: list[set[str]] = []
    seen.clear()

    def rdfs(v: str, comp: set[str]):
        seen.add(v)
        comp.add(v)
        for w in rev.get(v, ()):
            if w not in seen:
                rdfs(w, comp)

    for v in reversed(order):
        if v not in seen:
            comp: set[str] = set()
            rdfs(v, comp)
            comps.append(comp)
    return comps


def has_bad_subshift_cycle(word: tuple[str, ...]) -> bool:
    graph = transition_graph(word, exclude_targets=True)
    for comp in strongly_connected_components(graph):
        if len(comp) > 1:
            return True
        v = next(iter(comp))
        if v in graph.get(v, set()):
            return True
    return False


def first_occurrence_tournament(words: list[tuple[str, ...]]):
    letters = sorted(
        {s for w in words for s in w if s not in TARGETS},
        key=lambda s: TIE_ORDER.get(s, 99),
    )
    n = len(letters)
    if n == 0:
        return letters, tuple()
    index = {s: i for i, s in enumerate(letters)}
    adj = [[0] * n for _ in range(n)]
    for i, a in enumerate(letters):
        for j, b in enumerate(letters):
            if i == j:
                continue
            wins = 0
            losses = 0
            for word in words:
                pa = word.index(a) if a in word else None
                pb = word.index(b) if b in word else None
                if pa is None and pb is None:
                    continue
                if pb is None or (pa is not None and pa < pb):
                    wins += 1
                elif pa is None or pb < pa:
                    losses += 1
            if wins > losses:
                adj[i][j] = 1
            elif wins == losses and TIE_ORDER.get(a, 99) < TIE_ORDER.get(b, 99):
                adj[i][j] = 1
    return letters, tuple(tuple(row) for row in adj)


def score_histogram(adj: tuple[tuple[int, ...], ...]) -> tuple[int, ...]:
    return tuple(sorted(sum(row) for row in adj))


def directed_triangles(adj: tuple[tuple[int, ...], ...]) -> int:
    n = len(adj)
    count = 0
    for a, b, c in combinations(range(n), 3):
        edges = adj[a][b] + adj[b][c] + adj[c][a]
        if edges in (0, 3):
            count += 1
    return count


def tournament_scc_count(adj: tuple[tuple[int, ...], ...]) -> int:
    graph = {
        str(i): {str(j) for j, bit in enumerate(row) if bit}
        for i, row in enumerate(adj)
    }
    return len(strongly_connected_components(graph))


def hamiltonian_path_count(adj: tuple[tuple[int, ...], ...]) -> int:
    n = len(adj)
    if n == 0:
        return 1
    total = 0
    for p in permutations(range(n)):
        if all(adj[p[i]][p[i + 1]] for i in range(n - 1)):
            total += 1
    return total


def analyze_speed_set(speeds: tuple[int, ...], n: int) -> dict:
    word = compact_word(speeds, n)
    targets = target_positions(word)
    rwords = return_words(word)
    letters, adj = first_occurrence_tournament(rwords)
    open_hit = "G" in word
    wall_hit = "W" in word
    return {
        "speeds": speeds,
        "word_len": len(word),
        "open_hit": open_hit,
        "wall_hit": wall_hit,
        "target_count": len(targets),
        "max_gap": longest_target_gap(word),
        "p": tuple(len(cyclic_blocks(word, k)) for k in range(1, 5)),
        "return_types": len(set(rwords)),
        "bad_cycle": has_bad_subshift_cycle(word),
        "letters": tuple(letters),
        "score_hist": score_histogram(adj),
        "triangles": directed_triangles(adj),
        "scc": tournament_scc_count(adj) if adj else 0,
        "hp": hamiltonian_path_count(adj),
        "word_sample": "".join(word[:80]),
    }


def summarize(n: int, max_speed: int) -> list[dict]:
    rows = [analyze_speed_set(speeds, n) for speeds in primitive_speed_sets(n, max_speed)]
    print(f"\n=== n={n}, max_speed={max_speed}, primitive_sets={len(rows)} ===")
    open_hits = sum(1 for row in rows if row["open_hit"])
    wall_only = sum(1 for row in rows if (not row["open_hit"]) and row["wall_hit"])
    missing = [row for row in rows if not row["open_hit"] and not row["wall_hit"]]
    print(f"target coverage: open={open_hits}, wall_only={wall_only}, missing={len(missing)}")

    if rows:
        avg_p = []
        for k in range(4):
            avg_p.append(sum(row["p"][k] for row in rows) / len(rows))
        finite_gaps = [row["max_gap"] for row in rows if row["max_gap"] is not None]
        print(
            "avg block complexity p(1..4): "
            + ", ".join(f"{x:.2f}" for x in avg_p)
        )
        print(
            f"max target gap={max(finite_gaps) if finite_gaps else 'none'}, "
            f"avg return types={sum(row['return_types'] for row in rows) / len(rows):.2f}, "
            f"bad-subshift-cycle sets={sum(1 for row in rows if row['bad_cycle'])}/{len(rows)}"
        )

    by_score = Counter(row["score_hist"] for row in rows)
    by_tri = Counter(row["triangles"] for row in rows)
    by_scc = Counter(row["scc"] for row in rows)
    by_hp = Counter(row["hp"] for row in rows)
    print(f"return-order score histograms: {dict(sorted(by_score.items()))}")
    print(f"return-order directed 3-cycles: {dict(sorted(by_tri.items()))}")
    print(f"return-order SCC counts: {dict(sorted(by_scc.items()))}")
    print(f"return-order Hamiltonian-path counts: {dict(sorted(by_hp.items()))}")

    for row in rows:
        if not row["open_hit"] or row["speeds"] == tuple(range(1, n)):
            label = "AP/wall" if row["speeds"] == tuple(range(1, n)) else "no-open"
            print(
                f"example {label} speeds={row['speeds']} "
                f"open={row['open_hit']} wall={row['wall_hit']} "
                f"targets={row['target_count']} max_gap={row['max_gap']} "
                f"p={row['p']} return_types={row['return_types']} "
                f"letters={row['letters']} score={row['score_hist']} "
                f"tri={row['triangles']} hp={row['hp']}"
            )
            print(f"  compact word prefix: {row['word_sample']}")
            break
    return rows


def main() -> None:
    configs = [(4, 10), (5, 8), (6, 7), (7, 8)]
    all_rows = []
    for n, max_speed in configs:
        all_rows.extend((n, row) for row in summarize(n, max_speed))

    missing = [(n, row["speeds"]) for n, row in all_rows if not row["open_hit"] and not row["wall_hit"]]
    print("\n=== compactified target-free candidates ===")
    if missing:
        for n, speeds in missing:
            print(f"missing target: n={n}, speeds={speeds}")
    else:
        print("none in the bounded scan")

    print("\n=== interpretation ===")
    print(
        "The open word can avoid G on the AP wall, but the compactified wall "
        "symbol W repairs it.  A genuine counterexample would be a periodic "
        "arithmetic word avoiding both G and W.  The bad-subshift cycle counts "
        "show that coarse symbolic factors often admit spurious target-free "
        "cycles; the proof problem is to rule out arithmetic realizability of "
        "those cycles or refine the alphabet until they disappear."
    )


if __name__ == "__main__":
    main()
