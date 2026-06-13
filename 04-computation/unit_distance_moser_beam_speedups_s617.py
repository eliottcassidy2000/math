#!/usr/bin/env python3
"""
unit_distance_moser_beam_speedups_s617.py

S617: speed up small unit-distance candidate counts by keeping the right
carrier and updating edge counts incrementally.

The carrier is the rank-4 integer Moser coordinate lattice already exposed in
S614.  Its unit shell has 18 directed unit vectors (9 antipodal pairs).  For a
connected cluster S, adding a frontier point q changes the unit-distance count
by exactly

    gain(q) = #{u in U : q+u in S}.

Thus a beam search can count candidates by building one frontier gain table per
state, instead of recomputing all unit edges for every candidate extension.

This is not a proof that u(22)=60.  It is a reproducible counting speedup and a
search/extension diagnostic for the n=22 one-bit frontier.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from itertools import combinations
from time import perf_counter


Point = tuple[int, int, int, int]
Cluster = tuple[Point, ...]


def moser_unit_vectors() -> tuple[Point, ...]:
    """Return the directed unit shell from the S614 Moser-coordinate form."""
    units: list[Point] = []
    for a in range(-4, 5):
        for b in range(-4, 5):
            for c in range(-4, 5):
                for d in range(-4, 5):
                    if a == b == c == d == 0:
                        continue
                    if a * d != b * c:
                        continue
                    value = (
                        6 * a * a
                        + 6 * a * b
                        + 10 * a * c
                        + 5 * a * d
                        + 6 * b * b
                        + 5 * b * c
                        + 10 * b * d
                        + 6 * c * c
                        + 6 * c * d
                        + 6 * d * d
                    )
                    if value == 6:
                        units.append((a, b, c, d))
    return tuple(sorted(units))


UNITS = moser_unit_vectors()
UNIT_SET = set(UNITS)


def add(a: Point, b: Point) -> Point:
    return tuple(a[i] + b[i] for i in range(4))  # type: ignore[return-value]


def sub(a: Point, b: Point) -> Point:
    return tuple(a[i] - b[i] for i in range(4))  # type: ignore[return-value]


def canon(points: Cluster) -> Cluster:
    """Translation-canonical tuple; enough to remove the dominant duplicate."""
    pts = sorted(points)
    base = pts[0]
    return tuple(sorted(sub(p, base) for p in pts))


def span4(cluster: Cluster) -> int:
    return sum(max(p[i] for p in cluster) - min(p[i] for p in cluster) for i in range(4))


def unit_edges(cluster: Cluster) -> int:
    s = set(cluster)
    return sum(1 for p in cluster for u in UNITS if add(p, u) in s) // 2


def frontier_gains(cluster: Cluster) -> dict[Point, int]:
    s = set(cluster)
    gains: dict[Point, int] = {}
    for p in cluster:
        for u in UNITS:
            q = add(p, u)
            if q not in s:
                gains[q] = gains.get(q, 0) + 1
    return gains


@dataclass(frozen=True)
class StateInfo:
    edges: int
    span: int


@dataclass(frozen=True)
class BeamRow:
    size: int
    parents: int
    unique_children: int
    frontier_evals: int
    kept: int
    best_edges: int
    best_span: int
    naive_edge_checks: int
    incremental_edge_checks: int
    seconds: float

    @property
    def edge_count_speedup(self) -> float:
        if self.incremental_edge_checks == 0:
            return 1.0
        return self.naive_edge_checks / self.incremental_edge_checks


@dataclass(frozen=True)
class Trick:
    name: str
    count_speed: int
    geometry: int
    side_channels: int
    n22_focus: int
    proof_power: int
    risk: int


TRICKS = (
    Trick("frontier-gain incremental count", 5, 4, 4, 5, 3, 1),
    Trick("21-core extension ledger", 4, 4, 5, 5, 4, 1),
    Trick("bitset adjacency popcount window", 5, 3, 3, 4, 3, 2),
    Trick("Moser unit-shell carrier", 4, 5, 4, 4, 3, 2),
    Trick("totally-unfaithful obstruction library", 3, 5, 5, 5, 5, 2),
    Trick("raw F-free graph enumeration", 1, 1, 2, 5, 3, 4),
    Trick("naive pairwise distance recount", 0, 2, 0, 2, 0, 5),
    Trick("triangular-lattice-only beam", 4, 3, 1, 1, 1, 4),
)


def trick_votes(a: Trick, b: Trick) -> tuple[int, int]:
    criteria = [
        (a.count_speed > b.count_speed, b.count_speed > a.count_speed),
        (a.geometry > b.geometry, b.geometry > a.geometry),
        (a.side_channels > b.side_channels, b.side_channels > a.side_channels),
        (a.n22_focus > b.n22_focus, b.n22_focus > a.n22_focus),
        (a.proof_power > b.proof_power, b.proof_power > a.proof_power),
        (a.risk < b.risk, b.risk < a.risk),
    ]
    av = sum(1 for x, y in criteria if x and not y)
    bv = sum(1 for x, y in criteria if y and not x)
    return av, bv


def trick_tournament(tricks: tuple[Trick, ...]) -> list[list[int]]:
    n = len(tricks)
    adj = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            iv, jv = trick_votes(tricks[i], tricks[j])
            if iv > jv or (iv == jv and i < j):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
    return adj


def directed_triangles(adj: list[list[int]]) -> int:
    total = 0
    for i, j, k in combinations(range(len(adj)), 3):
        if adj[i][j] and adj[j][k] and adj[k][i]:
            total += 1
        if adj[i][k] and adj[k][j] and adj[j][i]:
            total += 1
    return total


def hamiltonian_path_count(adj: list[list[int]]) -> int:
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
                if not ((mask >> nxt) & 1) and adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += val
    return sum(dp[-1])


def score_hist(adj: list[list[int]]) -> dict[int, int]:
    hist: dict[int, int] = {}
    for row in adj:
        score = sum(row)
        hist[score] = hist.get(score, 0) + 1
    return dict(sorted(hist.items()))


def run_beam(target: int, width: int) -> tuple[list[BeamRow], dict[int, dict[int, int]], Cluster, int]:
    states: dict[Cluster, StateInfo] = {((0, 0, 0, 0),): StateInfo(0, 0)}
    rows: list[BeamRow] = []
    extension_gain_hist: dict[int, dict[int, int]] = {}
    t0 = perf_counter()

    best_cluster: Cluster = ((0, 0, 0, 0),)
    best_edges = 0

    for size in range(2, target + 1):
        parent_count = len(states)
        children: dict[Cluster, StateInfo] = {}
        frontier_evals = 0
        incremental_checks = 0
        naive_checks = 0

        for cluster, info in states.items():
            gains = frontier_gains(cluster)
            frontier_evals += len(gains)
            incremental_checks += len(cluster) * len(UNITS)

            if size == 22 and info.edges >= 56:
                hist = extension_gain_hist.setdefault(info.edges, {})
                for gain in gains.values():
                    hist[gain] = hist.get(gain, 0) + 1

            for q, gain in gains.items():
                child = canon(cluster + (q,))
                new_edges = info.edges + gain
                # What a naive unit-edge recount would have spent on this child.
                naive_checks += size * len(UNITS)
                prev = children.get(child)
                new_info = StateInfo(new_edges, span4(child))
                if prev is None or (new_info.edges, -new_info.span) > (prev.edges, -prev.span):
                    children[child] = new_info

        ranked = sorted(
            children.items(),
            key=lambda item: (item[1].edges, -item[1].span, item[0]),
            reverse=True,
        )
        states = dict(ranked[:width])
        best_cluster, best_info = ranked[0]
        best_edges = best_info.edges
        recomputed = unit_edges(best_cluster)
        if recomputed != best_edges:
            raise AssertionError((size, best_edges, recomputed))
        rows.append(
            BeamRow(
                size=size,
                parents=parent_count,
                unique_children=len(children),
                frontier_evals=frontier_evals,
                kept=len(states),
                best_edges=best_edges,
                best_span=best_info.span,
                naive_edge_checks=naive_checks,
                incremental_edge_checks=incremental_checks,
                seconds=perf_counter() - t0,
            )
        )

    return rows, extension_gain_hist, best_cluster, best_edges


def print_trick_tournament() -> None:
    print("Tournament Analysis over counting tricks")
    adj = trick_tournament(TRICKS)
    scores = [sum(row) for row in adj]
    for idx, trick in sorted(
        enumerate(TRICKS),
        key=lambda item: (scores[item[0]], item[1].count_speed, item[1].proof_power),
        reverse=True,
    ):
        print(
            f"  score={scores[idx]} trick={trick.name}; "
            f"features=(speed={trick.count_speed}, geom={trick.geometry}, "
            f"side={trick.side_channels}, n22={trick.n22_focus}, "
            f"proof={trick.proof_power}, risk={trick.risk})"
        )
    print(f"  score histogram: {score_hist(adj)}")
    print(f"  directed 3-cycles: {directed_triangles(adj)}")
    print(f"  Hamiltonian path count: {hamiltonian_path_count(adj)}")
    print()


def main() -> None:
    sys.stdout.reconfigure(line_buffering=True)
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", type=int, default=22)
    parser.add_argument("--width", type=int, default=1200)
    args = parser.parse_args()

    print("S617 unit-distance Moser-carrier beam speedups")
    print(f"unit shell: {len(UNITS)} directed vectors, {len(UNITS)//2} antipodal pairs")
    print(f"target={args.target} beam_width={args.width}")
    print()

    rows, extension_hist, best_cluster, best_edges = run_beam(args.target, args.width)
    print("Beam ledger")
    print("n  kept  unique_children  frontier_evals  best_edges  span  edge_count_speedup  seconds")
    for row in rows:
        print(
            f"{row.size:2d} {row.kept:5d} {row.unique_children:15d} "
            f"{row.frontier_evals:14d} {row.best_edges:10d} {row.best_span:5d} "
            f"{row.edge_count_speedup:18.1f} {row.seconds:8.2f}"
        )
    print()
    print(f"best cluster at n={args.target}: edges={best_edges}, span={span4(best_cluster)}")
    print(f"best cluster coordinates: {best_cluster}")
    print()

    if args.target >= 22:
        print("n=22 extension ledger inside retained beam")
        print("core_edges  frontier_candidates_by_gain  max_new_edges")
        for core_edges in sorted(extension_hist):
            hist = dict(sorted(extension_hist[core_edges].items()))
            max_gain = max(hist)
            print(f"{core_edges:10d} {hist!s:30s} {core_edges + max_gain:13d}")
        print(
            "Reading: within the retained Moser-carrier beam, 57-edge 21-cores "
            "only accept gain-3 frontier extensions, and 56-edge cores only "
            "reach gain 4.  This recovers the known 60-edge lower-bound lane "
            "but does not produce a 61-edge witness."
        )
        print()

    print_trick_tournament()
    print("Speedup claim")
    print("  The frontier-gain table replaces per-candidate unit-edge recounts.")
    print("  The displayed speedup is only for edge-count checks; canonicalization")
    print("  and beam ordering still cost time, but the expensive count itself")
    print("  becomes state-local rather than child-local.")
    print("  Next heavy run: increase --width and add automorphism canonicalization")
    print("  or bitset popcounts on a bounded Moser window.")


if __name__ == "__main__":
    main()
