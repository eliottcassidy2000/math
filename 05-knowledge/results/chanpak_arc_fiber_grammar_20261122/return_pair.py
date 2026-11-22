"""Ordered-slot transfer and least equal-height arcs seen by a triangle return."""
from collections import Counter
from functools import cache
from itertools import combinations
from math import comb
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "chanpak_gap_star_ports_20261122"))
from replay import flipped, subpaths


def slot_reader(rows):
    """Compute slot values from labelled induced path counts and ordered covers."""
    h = subpaths(rows)
    h[0] = 1
    full = len(h) - 1

    @cache
    def covers(mask, slots):
        if not slots:
            return int(mask == 0)
        if slots == 1:
            return h[mask]
        total = 0
        part = mask
        while True:
            total += h[part] * covers(mask ^ part, slots - 1)
            if not part:
                return total
            part = (part - 1) & mask

    return lambda slots: covers(full, slots)


def return_pair(rows):
    """Orient S -> c -> b -> S; the original vertices keep their labels."""
    n = len(rows)
    return tuple(row | (1 << n) for row in rows) + (1 << (n + 1), (1 << n) - 1)


def small_catalog(order):
    path = HERE.parents[2] / "session212" / "cycle-levels" / f"classes{order}.rows"
    for line in path.read_text().splitlines():
        yield tuple(map(int, line.split()))


def four_blocks(a, b, reverse=False):
    """Label the singletons first, followed by A and B."""
    m, n = len(a), len(b)
    left = ((1 << m) - 1) << 2
    right = ((1 << n) - 1) << (m + 2)
    return ((left | (0 if reverse else 2), right | int(reverse)) +
            tuple((row << 2) | right | 2 for row in a) +
            tuple((row << (m + 2)) | 1 for row in b))


def determinant_change(a, b, j):
    return j * (a(j + 1) * b(j) - a(j) * b(j + 1)
                + a(j) * b(j - 1) - a(j - 1) * b(j))


def host_reader(a, b, reverse):
    """Count absolute covers in the square using compressed module visits."""
    base = four_blocks((0,), (0,), reverse)
    singles = (0, 1)
    slots = Counter()

    def extend(last, mask, left, right):
        slots[mask, left, right] += 1
        for vertex in range(4):
            if last >= 0 and not (base[last] >> vertex & 1):
                continue
            if vertex in singles and mask >> vertex & 1:
                continue
            extend(vertex, mask | ((1 << vertex) if vertex in singles else 0),
                   left + (vertex == 2), right + (vertex == 3))

    extend(-1, 0, 0, 0)

    @cache
    def nonempty(module, k):
        return sum((-1)**(k-i) * comb(k, i) * module(i) for i in range(k+1))

    def read(j):
        current = Counter({(0, 0, 0): 1})
        for _ in range(j):
            following = Counter()
            for (mask, x, y), count in current.items():
                for (added, u, v), ways in slots.items():
                    if not mask & added:
                        following[mask | added, x + u, y + v] += count * ways
            current = following
        return sum(count * nonempty(a, x) * nonempty(b, y)
                   for (mask, x, y), count in current.items() if mask == 3)

    return read


def delayed_flip():
    for left, right in [((2, 4, 1), (2, 0)), ((2, 0), (0,))]:
        a, b = slot_reader(left), slot_reader(right)
        graphs = [four_blocks(left, right, reverse) for reverse in (False, True)]
        assert flipped(graphs[0], 0, 1) == graphs[1]
        actual = [slot_reader(graph) for graph in graphs]
        compressed = [host_reader(a, b, reverse) for reverse in (False, True)]
        for j in range(1, 4):
            assert [f(j) for f in actual] == [f(j) for f in compressed]
            assert actual[1](j) - actual[0](j) == determinant_change(a, b, j)

    # The two order-15 modules are ordinal sums of the archived order-7 and -8 factors.
    factors = (((8,9,3,100,15,23,55), (144,5,161,7,78,27,47,122)),
               ((32,97,19,7,11,28,61), (112,97,67,7,14,28,56,127)))
    def ordinal_reader(rows):
        x, y = map(slot_reader, rows)
        return cache(lambda j: x(j) * y(j))
    a, b = map(ordinal_reader, factors)
    assert [a(j) for j in range(1, 4)] == [b(j) for j in range(1, 4)] == [17325,5105184,317017071]
    assert (a(4), b(4)) == (8535609600,8540125440)
    old, new = [host_reader(a, b, reverse) for reverse in (False, True)]
    for j in range(1, 4):
        assert new(j) - old(j) == determinant_change(a, b, j)
    q = [[f(j) for j in range(1, 4)] for f in (old, new)]
    # H(R^2)=2Q3+5Q2-Q1.
    print("32-vertex square (H,Q2,Q3)", *q)
    second = tuple(2*z + 5*y - x for x, y, z in q)
    assert second == (1644123578453985193299, 1644114988863765765459)
    print("second-wrap counts", *second)


def main():
    # The two independent counts use labelled path DP before and after the return.
    checks = 0
    for n in range(1, 6):
        for rows in small_catalog(n):
            before = slot_reader(rows)
            after = slot_reader(return_pair(rows))
            for j in range(1, 5):
                assert after(j) == j * (before(j + 1) + j * before(j) - before(j - 1))
            for i, j in combinations(range(n), 2):
                other = slot_reader(flipped(rows, i, j))
                if before(1) == other(1):
                    # At these orders values through 2 and monic parity determine Q.
                    assert before(2) == other(2)
                checks += 1
    assert checks == 151
    print("All 151 isomorphism-class/arc choices through order five checked")

    for rows, edge, expected in (
        ((6, 12, 24, 33, 11, 23), (0, 3), (37, 376, 352, 413, 389)),
        ((2, 12, 73, 49, 39, 71, 27), (3, 6), (91, 1256, 1244, 1347, 1335)),
    ):
        twin = flipped(rows, *edge)
        a, b = slot_reader(rows), slot_reader(twin)
        wa, wb = return_pair(rows), return_pair(twin)
        assert flipped(wa, *edge) == wb
        assert a(1) == b(1)
        assert (a(1), a(2), b(2), slot_reader(wa)(1), slot_reader(wb)(1)) == expected
        print("source order", len(rows), "marked arc", edge, "H,Q2,Q2',returned H,H'", expected)


if __name__ == "__main__":
    main()
    delayed_flip()
