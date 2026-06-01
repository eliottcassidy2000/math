#!/usr/bin/env python3
"""Fixed-n LRC source menus as tournament isomorphism-class atlases.

codex-2026-06-01 S525

For total LRC denominator n, a lonely/source state places the m=n-1 moving
runners inside the safe arc [1/n, 1-1/n], of length L=1-2/n.  Ordering those
points in the arc turns the moving runners into a half-turn tournament on m
vertices:

    i -> j for i<j iff x_j - x_i < 1/2.

Thus every backward edge, relative to arc order, is a pair whose separation is
larger than 1/2.  The backward pairs form a Ferrers filter in the triangular
pair poset.  This script enumerates those Ferrers interval tournaments and
compares their isomorphism classes with the open source classes reached by
bounded speed clocks.

Tournament Analysis declaration:
    vertices: fixed total denominators n=4..9
    pairwise observable: source-menu profile
        (open iso count, feasible Ferrers signatures, A000568 base size,
         max Hamiltonian-path count, number of H-values)
    switch/gauge: lexicographic comparison of the profile vector
    tie Hamiltonian path: increasing n

Stored output:
    05-knowledge/results/lrc_fixed_n_tournament_class_atlas_s525.out
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations
from math import comb, gcd


Adj = tuple[tuple[int, ...], ...]
Signature = tuple[int, ...]


A000568_TABLE = {
    1: 1,
    2: 1,
    3: 2,
    4: 4,
    5: 12,
    6: 56,
    7: 456,
    8: 6880,
    9: 191536,
    10: 9733056,
}


def a000568(m: int) -> int:
    return A000568_TABLE[m]


def catalan(m: int) -> int:
    return comb(2 * m, m) // (m + 1)


def ferrers_signatures(m: int):
    """Yield monotone threshold signatures for m ordered arc points.

    c_i is the first index j such that x_j-x_i > 1/2.  If no such j exists,
    c_i=m.  Monotonicity of differences forces c_0 <= c_1 <= ... <= c_{m-2}.
    """

    def rec(i: int, previous: int, cur: list[int]):
        if i == m - 1:
            yield tuple(cur)
            return
        for value in range(max(previous, i + 1), m + 1):
            yield from rec(i + 1, value, cur + [value])

    yield from rec(0, 1, [])


def signature_adjacency(signature: Signature) -> Adj:
    m = len(signature) + 1
    adj = [[0] * m for _ in range(m)]
    for i, j in combinations(range(m), 2):
        if j >= signature[i]:
            adj[j][i] = 1
        else:
            adj[i][j] = 1
    return tuple(tuple(row) for row in adj)


def feasible_at_length(signature: Signature, length: float, eps: float = 1e-8) -> bool:
    """Difference-constraint feasibility for an open source chamber.

    Variables are ordered positions x_0<...<x_{m-1}.  The system is translation
    invariant, so feasibility is detected by absence of a negative cycle.
    """

    m = len(signature) + 1
    edges: list[tuple[int, int, float]] = []

    def add(u: int, v: int, bound: float) -> None:
        # x_v - x_u <= bound
        edges.append((u, v, bound))

    for i in range(m - 1):
        add(i + 1, i, -eps)  # x_{i+1}-x_i >= eps
    add(0, m - 1, length)  # total span <= safe arc length

    for i, j in combinations(range(m), 2):
        if j >= signature[i]:
            add(j, i, -0.5 - eps)  # x_j-x_i > 1/2
        else:
            add(i, j, 0.5 - eps)  # x_j-x_i < 1/2

    dist = [0.0] * m
    for _ in range(m):
        changed = False
        for u, v, bound in edges:
            candidate = dist[u] + bound
            if dist[v] > candidate + 1e-12:
                dist[v] = candidate
                changed = True
        if not changed:
            return True
    return False


@lru_cache(maxsize=None)
def perms(m: int) -> tuple[tuple[int, ...], ...]:
    return tuple(permutations(range(m)))


@lru_cache(maxsize=None)
def canon(adj: Adj) -> tuple[int, ...]:
    m = len(adj)
    best: tuple[int, ...] | None = None
    for p in perms(m):
        flat = tuple(adj[p[i]][p[j]] for i in range(m) for j in range(m) if i != j)
        if best is None or flat < best:
            best = flat
    assert best is not None
    return best


def hamiltonian_paths(adj: Adj) -> int:
    m = len(adj)
    full = (1 << m) - 1
    dp = [[0] * m for _ in range(1 << m)]
    for start in range(m):
        dp[1 << start][start] = 1
    for mask in range(1 << m):
        for last in range(m):
            cur = dp[mask][last]
            if not cur:
                continue
            for nxt in range(m):
                if not ((mask >> nxt) & 1) and adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += cur
    return sum(dp[full])


def score_sequence(adj: Adj) -> tuple[int, ...]:
    return tuple(sorted(sum(row) for row in adj))


def c3_count(adj: Adj) -> int:
    count = 0
    for a, b, c in combinations(range(len(adj)), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (
            adj[a][c] and adj[c][b] and adj[b][a]
        ):
            count += 1
    return count


def open_interval_menu(total_n: int) -> dict[tuple[int, ...], dict]:
    m = total_n - 1
    length = 1 - 2 / total_n
    classes: dict[tuple[int, ...], dict] = {}
    for signature in ferrers_signatures(m):
        if not feasible_at_length(signature, length):
            continue
        adj = signature_adjacency(signature)
        key = canon(adj)
        if key not in classes:
            classes[key] = {
                "adj": adj,
                "signatures": [],
                "H": hamiltonian_paths(adj),
                "score": score_sequence(adj),
                "c3": c3_count(adj),
            }
        classes[key]["signatures"].append(signature)
    return classes


def primitive(combo: tuple[int, ...]) -> bool:
    g = 0
    for value in combo:
        g = gcd(g, value)
    return g == 1


def frac(value: Fraction) -> Fraction:
    return value - (value.numerator // value.denominator)


def dist0(value: Fraction) -> Fraction:
    f = frac(value)
    return min(f, 1 - f)


def speed_walls(speeds: tuple[int, ...], total_n: int) -> list[Fraction]:
    walls = {Fraction(0), Fraction(1)}
    threshold = Fraction(1, total_n)
    for speed in speeds[1:]:
        for k in range(speed):
            walls.add(frac((Fraction(k) + threshold) / speed))
            walls.add(frac((Fraction(k) - threshold) / speed))
    for i, j in combinations(range(1, total_n), 2):
        d = abs(speeds[i] - speeds[j])
        if d == 0:
            continue
        for k in range(2 * d):
            walls.add(Fraction(k, 2 * d))
    return sorted(walls)


def open_speed_menu(total_n: int, max_speed: int) -> dict[tuple[int, ...], dict]:
    m = total_n - 1
    threshold = Fraction(1, total_n)
    classes: dict[tuple[int, ...], dict] = {}
    total_sets = 0
    source_sets = 0
    for combo in combinations(range(1, max_speed + 1), m):
        if not primitive(combo):
            continue
        speeds = (0,) + combo
        total_sets += 1
        reached = False
        walls = speed_walls(speeds, total_n)
        for left, right in zip(walls, walls[1:]):
            if not left < right:
                continue
            t = (left + right) / 2
            if any(dist0(Fraction(speed) * t) < threshold for speed in speeds[1:]):
                continue
            reached = True
            adj = [[0] * m for _ in range(m)]
            for a, b in combinations(range(m), 2):
                phase = frac(Fraction(speeds[a + 1] - speeds[b + 1]) * t)
                if 0 < phase < Fraction(1, 2):
                    adj[a][b] = 1
                else:
                    adj[b][a] = 1
            adj_t = tuple(tuple(row) for row in adj)
            key = canon(adj_t)
            if key not in classes:
                classes[key] = {
                    "adj": adj_t,
                    "H": hamiltonian_paths(adj_t),
                    "score": score_sequence(adj_t),
                    "first_speeds": speeds,
                }
        if reached:
            source_sets += 1
    classes["_meta"] = {"total_sets": total_sets, "source_sets": source_sets}  # type: ignore[index]
    return classes


def class_spectrum(classes: dict[tuple[int, ...], dict]) -> tuple[tuple[int, tuple[int, ...], int, int], ...]:
    rows = []
    for data in classes.values():
        rows.append((data["H"], data["score"], data["c3"], len(data["signatures"])))
    return tuple(sorted(rows))


def profile_vector(total_n: int, classes: dict[tuple[int, ...], dict]) -> tuple:
    h_values = {data["H"] for data in classes.values()}
    feasible = sum(len(data["signatures"]) for data in classes.values())
    return (
        len(classes),
        feasible,
        a000568(total_n - 1),
        max(h_values),
        len(h_values),
        total_n,
    )


def orient_tournament(rows: list[tuple[int, dict[tuple[int, ...], dict]]]) -> list[list[bool]]:
    m = len(rows)
    adj = [[False] * m for _ in range(m)]
    profiles = [profile_vector(total_n, classes) for total_n, classes in rows]
    for i in range(m):
        for j in range(i + 1, m):
            if profiles[i] > profiles[j] or (profiles[i] == profiles[j] and i < j):
                adj[i][j] = True
            else:
                adj[j][i] = True
    return adj


def score_hist(adj: list[list[bool]]) -> tuple[tuple[int, int], ...]:
    return tuple(sorted(Counter(sum(row) for row in adj).items()))


def tournament_c3(adj: list[list[bool]]) -> int:
    count = 0
    for a, b, c in combinations(range(len(adj)), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (
            adj[a][c] and adj[c][b] and adj[b][a]
        ):
            count += 1
    return count


def scc_sizes(adj: list[list[bool]]) -> tuple[int, ...]:
    n = len(adj)
    graph = [[j for j, edge in enumerate(adj[i]) if edge] for i in range(n)]
    rgraph = [[] for _ in range(n)]
    for i, row in enumerate(graph):
        for j in row:
            rgraph[j].append(i)

    seen: set[int] = set()
    order: list[int] = []

    def dfs(v: int) -> None:
        seen.add(v)
        for w in graph[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for v in range(n):
        if v not in seen:
            dfs(v)

    seen.clear()
    sizes = []

    def rdfs(v: int) -> int:
        seen.add(v)
        size = 1
        for w in rgraph[v]:
            if w not in seen:
                size += rdfs(w)
        return size

    for v in reversed(order):
        if v not in seen:
            sizes.append(rdfs(v))
    return tuple(sorted(sizes, reverse=True))


def hamiltonian_path_count_bool(adj: list[list[bool]]) -> int:
    n = len(adj)
    full = (1 << n) - 1
    dp: dict[tuple[int, int], int] = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1 << n):
        for last in range(n):
            cur = dp.get((mask, last), 0)
            if not cur:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if adj[last][nxt]:
                    dp[(mask | (1 << nxt), nxt)] = dp.get((mask | (1 << nxt), nxt), 0) + cur
    return sum(dp.get((full, v), 0) for v in range(n))


def main() -> None:
    print("LRC fixed-n tournament class atlas -- codex-2026-06-01 S525")
    print()
    print("PART A: open source menus as Ferrers interval tournaments")
    print("=" * 72)
    rows: list[tuple[int, dict[tuple[int, ...], dict]]] = []
    for total_n in range(4, 10):
        classes = open_interval_menu(total_n)
        rows.append((total_n, classes))
        feasible = sum(len(data["signatures"]) for data in classes.values())
        h_values = sorted({data["H"] for data in classes.values()})
        m = total_n - 1
        print(
            "n={} m={} L=1-2/n={:.6f} A000568(m)={} Catalan={} feasible_signatures={} open_iso={} H={}".format(
                total_n,
                m,
                1 - 2 / total_n,
                a000568(m),
                catalan(m),
                feasible,
                len(classes),
                h_values,
            )
        )
    print()

    print("PART B: bounded speed-clock cross-check of the open menu")
    print("=" * 72)
    for total_n, max_speed in [(4, 14), (5, 12), (6, 11), (7, 10), (8, 9)]:
        geometric = set(open_interval_menu(total_n))
        speed_classes = open_speed_menu(total_n, max_speed)
        meta = speed_classes.pop("_meta")  # type: ignore[arg-type]
        speed_keys = set(speed_classes)
        print(
            "n={} max_speed={} primitive_sets={} source_sets={} speed_open_iso={} geometric_open_iso={} match={}".format(
                total_n,
                max_speed,
                meta["total_sets"],
                meta["source_sets"],
                len(speed_keys),
                len(geometric),
                speed_keys == geometric,
            )
        )
    print()
    print("Closed/wall compactification note:")
    print("  S512/S520 source menus include equality-wall witnesses.  The open")
    print("  Ferrers tournament menu gives 1,2,4,6,10 classes for n=4..8, while")
    print("  stored closed/wall speed menus report 2,2,6,6,>=12.  The surplus")
    print("  classes are boundary/tie-path data, not raw open A000568 classes.")
    print()

    print("PART C: class spectra")
    print("=" * 72)
    for total_n, classes in rows:
        print(f"n={total_n} open source class spectrum:")
        for H, score, c3, multiplicity in class_spectrum(classes):
            print(f"  H={H:<4} score={score} c3={c3:<2} ferrers_signatures={multiplicity}")
        print()

    print("PART D: Tournament Analysis across fixed-n profiles")
    print("=" * 72)
    profile_tournament = orient_tournament(rows)
    print(f"score_hist={score_hist(profile_tournament)}")
    print(f"c3={tournament_c3(profile_tournament)}")
    print(f"SCCs={scc_sizes(profile_tournament)}")
    print(f"Hamiltonian_paths={hamiltonian_path_count_bool(profile_tournament)}")
    print("top order by score:")
    scores = [(sum(profile_tournament[i]), rows[i][0], profile_vector(*rows[i])) for i in range(len(rows))]
    for score, total_n, profile in sorted(scores, reverse=True):
        print(f"  score={score} n={total_n} profile={profile}")
    print()

    print("SYNTHESIS")
    print("=" * 72)
    print("1. Fixed total n turns LRC into source reachability in a marked")
    print("   tournament on n vertices; after deleting the source observer, the")
    print("   target lives in tournament classes on m=n-1 vertices.")
    print("2. In an open source cell those m vertices lie in one arc of length")
    print("   1-2/n.  Relative to their arc order, backward edges form a Ferrers")
    print("   filter, and the fixed safe-arc length cuts Catalan many Ferrers")
    print("   signatures down to a much smaller interval-feasible menu.")
    print("3. Bounded speed clocks hit exactly the same open isomorphism classes")
    print("   for n=4..8.  Wall witnesses explain the extra closed-menu classes")
    print("   from S512/S520, so the true LRC object is an open tournament menu")
    print("   plus a compactified boundary fiber.")
    print("4. The problem at a particular n is therefore: prove every primitive")
    print("   arithmetic clock reaches this tiny Ferrers interval menu, or its")
    print("   compactified wall boundary, inside the A000568(n-1) class universe.")


if __name__ == "__main__":
    main()
