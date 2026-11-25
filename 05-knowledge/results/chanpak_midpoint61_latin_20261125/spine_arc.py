"""Literal Hamilton subset DP for the one-apex spine-flip classification."""
from functools import lru_cache
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "chanpak_gap_star_ports_20261122"))
from replay import flipped, profile, star_rows, subpaths


def spine(n):
    return tuple((1 << n) - (2 << v) for v in range(n))


@lru_cache(maxsize=1000)
def paths(rows):
    return subpaths(rows)[-1]


def midpoint(rows, i, j):
    before = paths(rows)
    after = paths(flipped(rows, i, j))
    return (before + after) // 2 if abs(before - after) == 2 else None


def allowed(center, interior):
    if interior == 0:
        return center in (2, 4, 10)
    if interior == 2:
        return center in (10, 18)
    if interior != 1:
        return False
    return center in (2, 4) or any(
        center in (6 * (1 << j), 9 * (1 << j) - 2)
        for j in range(1, center.bit_length() + 1)
    )


def small_tables():
    for m, expected, expected_zero_masks in (
        (1, {2, 4, 12, 16, 24}, [0, 3]), (2, {10, 18}, [0])
    ):
        n = m + 2
        core = spine(n)
        edited = flipped(core, 0, n - 1)
        centers = set()
        zero_masks = []
        for mask in range(1 << n):
            old, new = profile(core, mask), profile(edited, mask)
            if sum(new[2]) == sum(old[2]):
                zero_masks.append(mask)
            for right in range(4):
                rows = star_rows(core, mask, 0, right)
                value = midpoint(rows, 0, n - 1)
                if value is not None:
                    centers.add(value)
        assert centers == expected and zero_masks == expected_zero_masks
        print("interior", m, "small centers", sorted(centers), "zero masks", zero_masks)


def all_short_boards():
    total = units = 0
    for n in range(2, 8):
        core = spine(n)
        for mask in range(1 << n):
            rows = star_rows(core, mask, 0, 0)
            for i in range(n):
                for j in range(i + 1, n):
                    value = midpoint(rows, i, j)
                    if value is not None:
                        assert allowed(value, j - i - 1), (n, mask, i, j, value)
                        units += 1
                    total += 1
    print("literal word/arc choices", total, "unit gaps", units)


def family_examples():
    for j in range(1, 10):
        for bits, center in (("110" + "0" * j, 6 * (1 << j)),
                             ("1001" + "0" * j, 9 * (1 << j) - 2)):
            mask = sum(1 << i for i, bit in enumerate(bits) if bit == "1")
            rows = star_rows(spine(len(bits)), mask, 0, 0)
            assert paths(rows) == center - 1
            assert paths(flipped(rows, 0, 2)) == center + 1
    print("explicit strong families, both endpoints, j=1..9")


if __name__ == "__main__":
    small_tables()
    all_short_boards()
    family_examples()
