"""Ordered-slot transfer and least equal-height arcs seen by a triangle return."""
from functools import cache
from itertools import combinations
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
