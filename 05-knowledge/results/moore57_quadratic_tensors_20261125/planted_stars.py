"""Replay two partial Moore graphs on the actual 3,250 vertices."""
from collections import Counter
from itertools import combinations
from pathlib import Path

from three_stars import IDENTITY, SLOTS, make_chart, routes, symmetric_group_certificate
from translated_slots import contacts, lift, placement, read_transversal

HERE = Path(__file__).resolve().parent


def full_graph(chart, stars):
    def leaf(branch, sheet):
        return 58 + 56 * branch + sheet

    edges = [set() for _ in range(3250)]

    def join(x, y):
        assert x != y
        edges[x].add(y)
        edges[y].add(x)

    for branch in range(57):
        join(0, 1 + branch)
        for sheet in IDENTITY:
            join(1 + branch, leaf(branch, sheet))
    for (i, j), row in chart.items():
        if i < j:
            for x, y in enumerate(row):
                join(leaf(i, x), leaf(j, y))
    for (b, x), choices in stars.items():
        for h, y in choices:
            join(leaf(b, x), leaf(h, y))

    for i in range(57):
        for s in IDENTITY:
            neighbors = edges[leaf(i, s)]
            branches = [(v - 58) // 56 for v in neighbors if v >= 58]
            assert i not in branches and len(set(branches)) == len(branches)
    for b, x in stars:
        assert len(edges[leaf(b, x)]) == 57
    for branch in set(range(57)) - set(SLOTS):
        assert all(len(edges[leaf(branch, s)]) == 57 for s in IDENTITY)
    for v, neighbors in enumerate(edges):
        reached = set(neighbors) | {v}
        for w in neighbors:
            for u in edges[w]:
                if u != v:
                    assert u not in reached, (v, w, u)
                    reached.add(u)
    degrees = Counter(map(len, edges))
    assert max(degrees) == 57
    size = sum(d * n for d, n in degrees.items()) // 2
    print("Edges", size, "degree distribution", sorted(degrees.items()))
    return size, degrees[57]


def orbit_chart(filename, step, repairs=()):
    trades, selected = read_transversal(HERE / filename)
    chart = make_chart(trades)
    for row in chart.values():
        assert all(row[(s + step) % 56] == (row[s] + step) % 56 for s in IDENTITY)
    triples = [contacts(chart, h, s) for h, s in selected]
    chart = make_chart(trades + repairs)
    stars = placement(chart, triples, range(0, 56, step))
    assert set(stars) == {(selected[0][0], (selected[0][1] + d) % 56) for d in range(0, 56, step)}
    return chart, stars


if __name__ == '__main__':
    six, stars = lift(HERE / 'planted_site.txt')
    words = {(0, 28): (0, 43, 49), (0, 56): (1, 1, 5), (28, 56): (2, 28)}
    for pair in combinations((0, 28, 56), 2):
        print("Symmetric group", pair, symmetric_group_certificate(routes(six, *pair), words[pair]))
    assert full_graph(six, stars) == (12807, 232)
    seven, stars = orbit_chart('period_seven_site.txt', 8, ((10, 32, 0),))
    words = {(0, 28): (1, 24, 31), (0, 56): (10, 14), (28, 56): (11, 23)}
    for pair, word in words.items():
        print("Symmetric group", pair, symmetric_group_certificate(routes(seven, *pair), word))
    assert full_graph(seven, stars) == (12860, 233)
    fourteen, stars = orbit_chart('period_fourteen_site.txt', 4, ((10, 45, 2), (27, 49, 0)))
    words = {(0, 28): (27, 44, 38), (0, 56): (6, 45, 27), (28, 56): (8, 28, 45)}
    for pair, word in words.items():
        print("Symmetric group", pair, symmetric_group_certificate(routes(fourteen, *pair), word))
    assert full_graph(fourteen, stars) == (13231, 240)
