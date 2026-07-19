#!/usr/bin/env python3
"""Exact referee for THM-1239 curvature-erasure guardrails.

Checks the two-blocker family from m=8, the sharper self-similar z=14a
one-blocker ladder from its exact threshold m=42, the deletion-quartet
non-BAD margin, and the explicit global lonely witness.  All calculations
use Fraction and explicit failures rather than ``assert``.
"""

from fractions import Fraction as F


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def gap(c: int, k: int) -> tuple[F, F]:
    return F(14 * k + 1, 14 * c), F(14 * k + 13, 14 * c)


def tooth(speed: int, address: int) -> tuple[F, F]:
    return F(14 * address - 1, 14 * speed), F(14 * address + 1, 14 * speed)


def base_cracks(m: int) -> tuple[tuple[str, F, F], ...]:
    a = 7 * m + 1
    glo, ghi = gap(a, m)
    intervals = {r: tooth(a + r, m + 1) for r in range(1, 7)}
    cracks: list[tuple[str, F, F]] = [("L", glo, intervals[6][0])]
    for r in range(5, 0, -1):
        cracks.append((str(r), intervals[r + 1][1], intervals[r][0]))
    cracks.append(("R", intervals[1][1], ghi))
    require(all(left < right for _, left, right in cracks), "positive cracks")
    return tuple(cracks)


def strictly_contains(container: tuple[F, F], target: tuple[F, F]) -> bool:
    return container[0] < target[0] and target[1] < container[1]


two_blocker_rows = 0
for m in range(8, 10001):
    a = 7 * m + 1
    c = a + 7
    d = 7 * a + 4
    cracks = base_cracks(m)
    left = cracks[0]
    right = cracks[-1]
    require(strictly_contains(tooth(c, m + 1), (left[1], left[2])), "c left blocker")
    require(strictly_contains(tooth(c, m + 2), (right[1], right[2])), "c right blocker")
    for label, lo, hi in cracks[1:-1]:
        r = int(label)
        require(
            strictly_contains(tooth(d, 7 * m + 7 - r), (lo, hi)),
            "d internal blocker",
        )
    two_blocker_rows += 1

# Exact u-coordinate self-similarity and the z=14a star.
one_blocker_rows = 0
for m in range(42, 10001):
    a = 7 * m + 1
    z = 14 * a
    cracks = base_cracks(m)
    addresses = tuple(range(14 * m + 1, 14 * m + 14, 2))
    require(len(addresses) == 7, "odd address ladder")
    for (_, lo, hi), address in zip(cracks, addresses):
        require(strictly_contains(tooth(z, address), (lo, hi)), "z star blocker")

    # The quartet {1,2,3,4}, translated by t, has last gap 1-3t.
    _, ghi = gap(a, m)
    require(1 - 3 * ghi > F(2, 7), "uniform non-BAD quartet margin")
    require(1 - 3 * ghi - F(2, 7) == F(28 * m - 29, 14 * a),
            "non-BAD exact margin")
    one_blocker_rows += 1

# m=42 is the exact one-blocker threshold: the r=2 crack controls it.
for m, expected in ((41, False), (42, True)):
    a = 7 * m + 1
    condition = F(21, a + 3) < F(1, 14)
    require(condition is expected, "sharp r=2 threshold")

# Check the compact u-coordinate formulas on a deterministic bank.
u_formula_rows = 0
for m in range(1, 501):
    a = 7 * m + 1
    glo, ghi = gap(a, m)
    require(14 * a * (glo - F(1, 7)) == -1, "u left endpoint")
    require(14 * a * (ghi - F(1, 7)) == 11, "u right endpoint")
    for r in range(1, 7):
        lo, hi = tooth(a + r, m + 1)
        require(14 * a * (lo - F(1, 7)) == F(a * (11 - 2 * r), a + r),
                "u tooth left")
        require(14 * a * (hi - F(1, 7)) == F(a * (13 - 2 * r), a + r),
                "u tooth right")
        u_formula_rows += 1

# Explicit full 13-speed sanity row: selected gap is locally erased, yet the
# packet has a different exact global lonely phase.
V = (1, 2, 3, 4, 295, 296, 297, 298, 299, 300, 301, 302, 4130)
p, q = 44, 199
distance_numerators = tuple(min((v * p) % q, q - ((v * p) % q)) for v in V)
expected = (44, 88, 67, 23, 45, 89, 66, 22, 22, 66, 89, 45, 33)
require(distance_numerators == expected, "explicit global witness residues")
require(14 * min(distance_numerators) > q, "explicit witness is lonely")

print("THM-1239 TWO/ONE-BLOCKER CURVATURE ERASURE EXACT AUDIT")
print(f"two-blocker rows checked (m=8..10000) = {two_blocker_rows}")
print(f"one-blocker z=14a rows checked (m=42..10000) = {one_blocker_rows}")
print(f"u-coordinate tooth formulas checked = {u_formula_rows}")
print("sharp one-blocker threshold = m=42 (a=295)")
print("one-blocker address word = 14m+1,14m+3,...,14m+13")
print("two-blocker incidence words = 1000001 and 0111110")
print("explicit 13-speed witness = t=44/199, minimum=22/199")
print("RESULT: PASS")
