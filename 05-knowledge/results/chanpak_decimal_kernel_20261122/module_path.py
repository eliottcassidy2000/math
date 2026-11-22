"""Check the small and perfect-number examples for the mod-three path criterion."""

from itertools import combinations
from math import prod
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "chanpak_gap_star_ports_20261122"))
from replay import flipped, subpaths
from e31_certificates import ordinal


def height_and_factor(rows):
    counts = subpaths(rows)
    counts[0] = 1
    full = len(counts) - 1
    q = sum(count * counts[full ^ subset] for subset, count in enumerate(counts))
    return counts[full], counts[full] + q


def bad_proper_modules(rows):
    bad = []
    for size in range(2, len(rows)):
        for vertices in combinations(range(len(rows)), size):
            mask = sum(1 << v for v in vertices)
            if any(row & mask not in (0, mask) for v, row in enumerate(rows) if v not in vertices):
                continue
            sub = tuple(sum(((rows[v] >> u) & 1) << i for i, u in enumerate(vertices))
                        for v in vertices)
            if height_and_factor(sub)[0] % 3 == 0:
                bad.append(vertices)
    return bad


if __name__ == "__main__":
    # The inner directed triangle beats vertex 3, which beats 4, which beats the triangle.
    small = (10, 12, 9, 16, 7)
    edited = flipped(small, 0, 1)
    assert bad_proper_modules(small) == [(0, 1, 2)]
    assert height_and_factor(small)[0] == 15 and height_and_factor(edited)[0] == 9
    tri = (2, 4, 1)
    pair = ordinal(tri, tri)
    assert subpaths(pair)[-1] == 9 and subpaths(flipped(pair, 1, 4))[-1] == 17
    inner = ordinal((0,), pair)
    # inner -> a -> b -> inner
    a, b = len(inner), len(inner) + 1
    source = tuple(row | (1 << a) for row in inner) + (1 << b, (1 << a) - 1)
    assert subpaths(source)[-1] == 297
    assert all(subpaths(flipped(source, u, v))[-1] % 3 == 0
               for u, v in combinations(range(len(source)), 2))
    assert subpaths(flipped(flipped(source, 2, 5), a, b))[-1] == 17
    arms = [(2, 8, 35, 69, 15, 27, 55), (6, 12, 72, 129, 15, 31, 59, 119),
            (10, 12, 65, 164, 15, 23, 59, 119), (136, 37, 73, 82, 39, 141, 51, 94),
            (26, 28, 145, 100, 104, 71, 135, 59)]
    data = [height_and_factor(rows) for rows in arms]
    assert data == [(49, 861), (219, 3963), (257, 4681), (563, 8275), (595, 8723)]
    long_factor = (1 << 61) + 1
    lower = long_factor * prod(factor for _, factor in data)
    assert lower == (1 << 60) * ((1 << 61) - 1) - 1
    h, factor = height_and_factor(flipped(arms[1], 4, 6))
    upper_candidate = lower // data[1][1] * factor
    assert h == 215 and upper_candidate % 3 == 0
    print("Before and after the internal edit:", lower, upper_candidate)
