"""Check the midpoint tournaments and count all sheet-symmetric Hamilton cycles."""
import sys
from itertools import combinations
from pathlib import Path

from square36 import graph_edges, square

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "chanpak_ternary_buckets_20261122"))
from target_flip import flip_and_check


def check_midpoints():
    targets = set()
    orders = set()
    for line in (HERE / "midpoints.rows").read_text().splitlines():
        midpoint, u, v, *rows = map(int, line.split())
        n = len(rows)
        assert midpoint not in targets and 0 <= u < n and 0 <= v < n
        assert all(0 <= row < 1 << n and not (row >> i & 1)
                   for i, row in enumerate(rows))
        assert all((rows[i] >> j & 1) != (rows[j] >> i & 1)
                   for i, j in combinations(range(n), 2))
        assert rows[u] >> v & 1
        flip_and_check(midpoint, (u, v, rows))
        targets.add(midpoint)
        orders.add(n)
    assert targets == set(range(100, 3601, 100)) and max(orders) == 10
    return sorted(orders)


def paths_between(neighbors, start, end):
    full = (1 << len(neighbors)) - 1
    # Bit i represents vertex i+1. Do not consume the terminal early.
    path = [start]

    def extend(free):
        if free == 1 << (end - 1):
            if neighbors[path[-1] - 1] & free:
                yield tuple(path + [end])
            return
        available = neighbors[path[-1] - 1] & free & ~(1 << (end - 1))
        while available:
            bit = available & -available
            available ^= bit
            path.append(bit.bit_length())
            yield from extend(free ^ bit)
            path.pop()

    yield from extend(full ^ (1 << (start - 1)))


def lifts():
    n = 36
    neighbors = [0] * n
    for a, b in graph_edges(n):
        neighbors[a - 1] |= 1 << (b - 1)
        neighbors[b - 1] |= 1 << (a - 1)
    rails = [a for a in range(1, n + 1) if square(2*a)]
    labels = {100*a + s for a in range(1, n + 1) for s in (1, -1)}
    edges = sum(square(a+b) for a, b in combinations(labels, 2))
    assert rails == [2, 8, 18, 32] and edges == 118
    counts = []
    for a, b in combinations(rails, 2):
        count = 0
        for path in paths_between(neighbors, a, b):
            first = [100*v + (-1)**i for i, v in enumerate(path)]
            other = [100*v - (-1)**i for i, v in enumerate(path)]
            cycle = first + other[::-1]
            assert set(cycle) == labels and len(cycle) == len(labels)
            assert all(square(x+y) for x, y in zip(cycle, cycle[1:] + cycle[:1]))
            count += 1
        counts.append(count)
    assert counts == [20, 18, 42, 81, 66, 162]
    return counts


def signed_examples():
    edges = []
    for x in range(1, 37, 2):
        for sign in (1, -1):
            numerator = 3*x + sign
            y = numerator // (numerator & -numerator)
            if y <= 36 and square(x+y):
                edges.append((x, sign, y))
    assert edges == [(3, -1, 1), (9, 1, 7)]
    for x, j in ((2277, 4), (36437, 8)):
        assert 3*x + 1 == (1 << j)*427 and square(x+427)


if __name__ == "__main__":
    print("Midpoint certificates on orders", check_midpoints())
    counts = lifts()
    print("Hamilton paths per pair of vertical edges", counts, "total", sum(counts))
    signed_examples()
