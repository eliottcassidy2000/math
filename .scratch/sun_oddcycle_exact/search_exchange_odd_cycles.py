#!/usr/bin/env python3
"""Finite exact search for odd cycles in Sun two-role exchange graphs."""

from __future__ import annotations

import argparse
import itertools
import math
from collections import defaultdict, deque


LOWERS = (2, 3, 5, 7)
DEGREES = (2, 4, 6, 8)


def values(degree: int, lower: int, limit: int) -> list[tuple[int, int]]:
    out = []
    n = lower
    while (v := math.comb(n, degree)) <= limit:
        out.append((n, v))
        n += 1
    return out


def shared_count(a: tuple[int, ...], b: tuple[int, ...]) -> int:
    return sum(x == y for x, y in zip(a, b))


def odd_cycle(reps: list[tuple[int, ...]]) -> list[int] | None:
    adj = [set() for _ in reps]
    for i, j in itertools.combinations(range(len(reps)), 2):
        if shared_count(reps[i], reps[j]) == 2:
            adj[i].add(j)
            adj[j].add(i)
    color = {}
    for root in range(len(adj)):
        if root in color:
            continue
        color[root] = 0
        parent = {root: -1}
        todo = deque([root])
        while todo:
            v = todo.popleft()
            for u in sorted(adj[v]):
                if u not in color:
                    color[u] = color[v] ^ 1
                    parent[u] = v
                    todo.append(u)
                elif color[u] == color[v]:
                    pv = []
                    x = v
                    while x != -1:
                        pv.append(x)
                        x = parent[x]
                    pu = []
                    x = u
                    while x != -1:
                        pu.append(x)
                        x = parent[x]
                    pos = {x: i for i, x in enumerate(pv)}
                    meet = next(x for x in pu if x in pos)
                    return pv[: pos[meet] + 1] + list(reversed(pu[: pu.index(meet)])) + [v]
    return None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=10_000)
    parser.add_argument("--positive-higher", action="store_true")
    args = parser.parse_args()
    lowers = (2, 4, 6, 8) if args.positive_higher else LOWERS
    atoms = [values(d, lower, args.limit) for d, lower in zip(DEGREES, lowers)]
    fibres: dict[int, list[tuple[int, int, int, int]]] = defaultdict(list)
    tuples = 0
    for (w, a), (x, b), (y, c), (z, d) in itertools.product(*atoms):
        n = a + b + c + d
        if n <= args.limit:
            fibres[n].append((w, x, y, z))
            tuples += 1
    for n in sorted(fibres):
        witness = odd_cycle(fibres[n])
        if witness is not None:
            print(
                f"limit={args.limit} positive_higher={args.positive_higher} "
                f"cumulative_tuples={tuples} first_n={n}"
            )
            print(f"fibre_size={len(fibres[n])} closed_walk_length={len(witness)-1}")
            for i in witness:
                print("representation=" + ",".join(map(str, fibres[n][i])))
            return
    print(
        f"limit={args.limit} positive_higher={args.positive_higher} "
        f"cumulative_tuples={tuples} no_odd_exchange_cycle=TRUE"
    )


if __name__ == "__main__":
    main()
