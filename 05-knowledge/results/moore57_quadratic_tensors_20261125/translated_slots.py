"""Translate a single certified sheet assignment into six simultaneously saturated sheets."""
from collections import defaultdict
from pathlib import Path
import sys

from three_stars import A, B, IDENTITY, PIVOT, SLOTS, TRADES, make_chart, routes

CENTERS = (PIVOT, A, B)
EXTRA = ((1, 50, 3), (1, 50, 6))
SHIFTS = (0, 8, 16, 24, 40, 48)


def read_transversal(path):
    lines = Path(path).read_text().splitlines()
    count = int(lines[0].split()[1])
    trades = tuple(tuple(map(int, line.split())) for line in lines[1:count + 1])
    marker, h, x = lines[count + 1].split()
    assert marker == "selected"
    choices = [tuple(map(int, line.split())) for line in lines[count + 2:count + 55]]
    return trades, [(int(h), int(x)), *choices]


def contacts(chart, h, x):
    return tuple(chart[h, j][x] for j in CENTERS)


def placement(chart, reference, shifts):
    lookup = {contacts(chart, h, y): (h, y) for h in SLOTS for y in IDENTITY}
    assert len(lookup) == 3024
    stars = {}
    for d in shifts:
        source, *others = [lookup[tuple((v + d) % 56 for v in edge)] for edge in reference]
        b, x = source
        assert sorted([b] + [h for h, y in others]) == SLOTS
        for j in CENTERS:
            existing = [r[x] for r in routes(chart, b, j)]
            new = [chart[h, j][y] for h, y in others]
            assert sorted(existing + new) == list(IDENTITY)
        stars[source] = others
    return stars


def check_union(chart, stars):
    graph = defaultdict(set)
    for (i, j), perm in chart.items():
        for x, y in enumerate(perm):
            graph[i, x].add((j, y))
    additions = {tuple(sorted((source, target))) for source, targets in stars.items()
                 for target in targets}
    for u, v in additions:
        graph[u].add(v)
        graph[v].add(u)
    for u, v in additions:
        for r, s in ((u, v), (v, u)):
            assert sum(branch == s[0] for branch, sheet in graph[r]) == 1
        assert graph[u].isdisjoint(graph[v])
        for w in graph[u] - {v}:
            assert graph[w].isdisjoint(graph[v] - {u})
    return len(additions)


def lift(seed_file):
    extra, selection = read_transversal(seed_file)
    fresh = TRADES + extra
    original = make_chart(fresh)
    triple_list = [contacts(original, h, x) for h, x in selection]
    traded = make_chart(fresh + EXTRA)
    stars = placement(traded, triple_list, SHIFTS)
    assert set(stars) == {(1, 0), (50, 8), (50, 16), (50, 24), (1, 40), (1, 48)}
    print("Centers", sorted(stars), "distinct new edges", check_union(traded, stars))
    return traded, stars


if __name__ == "__main__":
    lift(sys.argv[1])
