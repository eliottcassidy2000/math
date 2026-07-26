#!/usr/bin/env python3
"""Exact physical sharpness companion for THM-2399.

The calculation uses Fraction arithmetic only.  It verifies a fully
centered, correctly valued (1,2,5) last-lane packet on one 49-point
orbit.  The local scalar cover holds at every orbit point, while the
THM-2396 clean set meets the orbit in exactly one point.
"""

from fractions import Fraction as F
from math import gcd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def valuation(n: int, prime: int) -> int:
    require(n > 0, "valuation input")
    exponent = 0
    while n % prime == 0:
        n //= prime
        exponent += 1
    return exponent


def circle_distance(value: F) -> F:
    residue = value % 1
    return min(residue, 1 - residue)


def danger(speed: int, point: F) -> bool:
    return circle_distance(speed * point) < F(1, 14)


def guard_danger(speed: int, point: F) -> bool:
    return circle_distance(speed * point) < F(1, 7)


def ordinary_word(speed: int, orbit: tuple[F, ...]) -> tuple[int, ...]:
    """One address in every normalized residue row."""

    addresses: list[list[int]] = [[] for _ in range(7)]
    for raw_index, point in enumerate(orbit):
        if danger(speed, point):
            index = (raw_index + 4) % 49
            addresses[index % 7].append(index // 7)
    require(all(len(row) == 1 for row in addresses), "ordinary word shape")
    return tuple(row[0] for row in addresses)


def guard_word(speed: int, orbit: tuple[F, ...]) -> tuple[tuple[int, int], ...]:
    """Two addresses in every normalized residue row."""

    addresses: list[list[int]] = [[] for _ in range(7)]
    for raw_index, point in enumerate(orbit):
        if guard_danger(speed, point):
            index = (raw_index + 4) % 49
            addresses[index % 7].append(index // 7)
    require(all(len(row) == 2 for row in addresses), "guard word shape")
    return tuple(tuple(sorted(row)) for row in addresses)


# A strict, rational base point.  Its denominator was chosen only to avoid
# every endpoint; the resulting word is stable on an open interval.
Z = F(8_104_501, 10_000_019)
ORBIT = tuple((Z + j) / 49 for j in range(49))

# Scalar row.  The blocker depths are exactly (1,2,5).
H = 1_273
Q_STAR = 77
LOWER_Q = (431, 643, 566, 720)
Q = (Q_STAR,) + LOWER_Q

C = (1, 13, 1_399_489)
c = (13, 169, 18_193_357)
V = C[2] // 49

require(C[2] == 49 * 13**4, "high quotient blocker")
require(c[2] == 13 * C[2] == 49 * 13**5, "high original blocker")
require(tuple(13 * value for value in C) == c, "blocker quotient law")
require(tuple(valuation(value, 13) for value in c) == (1, 2, 5), "profile")
require((C[0], C[1], c[0], c[1]) == (1, 13, 13, 169), "common core")
require(gcd(1, 91) == 1, "common-core unit")

require(valuation(Q_STAR, 7) == 1, "top septimal q depth")
require(all(valuation(value, 7) == 0 for value in (H,) + LOWER_Q), "lower layer")
require(valuation(C[2], 7) == valuation(c[2], 7) == 2, "high blocker depth")
require(all(value % 13 for value in (H,) + Q), "six thirteen-units")
require(len(set(Q)) == 5, "five distinct ordinary labels")
require(len(set((H,) + Q + c)) == 9, "nine distinct scalar speeds")
require(gcd(gcd(gcd(H, *Q), c[0]), gcd(c[1], c[2])) == 1, "primitive row")

# THM-2396's high-safe base condition.
require(not danger(V, Z), "V must be safe at the base")
require(not danger(13 * V, Z), "13V must be safe at the base")

# Exact normalized physical words.  These are not independently translated
# abstract masks: every one is evaluated at the same base Z.
A_WORD = ordinary_word(C[0], ORBIT)
B_WORD = ordinary_word(C[1], ORBIT)
C_WORD = ordinary_word(c[1], ORBIT)
G_WORD = guard_word(H, ORBIT)
LOWER_WORDS = tuple(ordinary_word(value, ORBIT) for value in LOWER_Q)

require(A_WORD == (0, 0, 0, 0, 0, 0, 0), "A word")
require(B_WORD == (1, 3, 5, 0, 2, 3, 5), "B word")
require(C_WORD == (0, 4, 1, 5, 1, 5, 2), "C word")
require(G_WORD == ((0, 1),) * 7, "guard word")
require(
    LOWER_WORDS
    == (
        (0, 4, 6, 3, 5, 2, 4),
        (1, 2, 3, 4, 4, 5, 6),
        (1, 5, 2, 6, 3, 6, 3),
        (1, 6, 4, 2, 6, 4, 2),
    ),
    "four lower-q words",
)

# The top q word is precisely normalized row zero.  Both high blockers are
# safe at all forty-nine physical points.
top_support = tuple(
    (raw_index + 4) % 49
    for raw_index, point in enumerate(ORBIT)
    if danger(Q_STAR, point)
)
require(len(top_support) == 7, "top word size")
require({index % 7 for index in top_support} == {0}, "top word row")
require(all(not danger(C[2], point) for point in ORBIT), "C3 orbit safety")
require(all(not danger(c[2], point) for point in ORBIT), "c3 orbit safety")

# On each of the six q_*-safe rows, guard plus four lower-q words form a
# six-address transversal.  There is one hole.  Five holes are B addresses;
# row three has the unique physical clean address, which is instead C.
holes: list[int] = []
for row in range(1, 7):
    occupied = set(G_WORD[row])
    occupied.update(word[row] for word in LOWER_WORDS)
    require(len(occupied) == 6, f"row {row} transversal")
    hole = next(iter(set(range(7)) - occupied))
    holes.append(hole)
require(tuple(holes) == (3, 5, 5, 2, 3, 5), "safe-row holes")
require(
    tuple(B_WORD[row] for row in range(1, 7)) == (3, 5, 0, 2, 3, 5),
    "B safe-row addresses",
)
require(holes[0:2] + holes[3:] == [3, 5, 2, 3, 5], "five B holes")
require(holes[2] == C_WORD[3] == 5, "unique C hole")

# Verify the actual local scalar cover directly, rather than inferring it
# only from the word table.  C_H is the guard-safe set.
local_cover_failures = []
clean_points = []
for raw_index, point in enumerate(ORBIT):
    guard_bit = guard_danger(H, point)
    q_bits = tuple(danger(value, point) for value in Q)
    blocker_bits = tuple(danger(value, point) for value in c)
    quotient_bits = tuple(danger(value, point) for value in C)
    K = int(guard_bit) + sum(danger(value, point) for value in LOWER_Q)

    if not guard_bit and not (any(q_bits) or any(blocker_bits)):
        local_cover_failures.append(raw_index)

    if (
        K == 0
        and not q_bits[0]
        and not blocker_bits[2]
        and not any(quotient_bits)
    ):
        index = (raw_index + 4) % 49
        clean_points.append((raw_index, index % 7, index // 7, point))

require(not local_cover_failures, "local scalar cover")
require(
    clean_points
    == [(34, 3, 5, F(348_105_147, 490_000_931))],
    "unique clean root",
)

# Every membership is strict.  The positive exact phase margin proves that
# the word and cover survive on an open base interval, so the witness is not
# an endpoint artefact.
phase_margins = []
for speed, radius in (
    ((H, F(1, 7)),)
    + tuple((value, F(1, 14)) for value in Q)
    + tuple((value, F(1, 14)) for value in C)
    + tuple((value, F(1, 14)) for value in c)
):
    for point in ORBIT:
        phase_margins.append(abs(circle_distance(speed * point) - radius))
strict_phase_margin = min(phase_margins)
require(strict_phase_margin > 0, "strict endpoint margin")

# This packet is deliberately not a global counterexample.  The displayed
# point is strictly safe for the guard and all eight ordinary/blocker speeds.
GLOBAL_SAFE_POINT = F(3, 14)
require(not guard_danger(H, GLOBAL_SAFE_POINT), "global hostile guard")
require(
    all(not danger(value, GLOBAL_SAFE_POINT) for value in Q + c),
    "global hostile ordinary factors",
)


print("theorem=THM-2399")
print("status=PROVED-CANDIDATE+VERIFIED-EXACT; independent-audit=PENDING")
print(f"profile=1,2,5; base={Z}; orbit_size={len(ORBIT)}")
print(f"H={H}; qstar={Q_STAR}; lower_q={','.join(map(str, LOWER_Q))}")
print(f"C={','.join(map(str, C))}; c={','.join(map(str, c))}")
print(f"top_rows=0; safe_row_holes={','.join(map(str, holes))}")
print("local_scalar_cover_failures=0")
print(
    "clean_count=1;"
    f" clean_raw_index={clean_points[0][0]};"
    f" clean_normalized_address=({clean_points[0][1]},{clean_points[0][2]})"
)
print(f"strict_phase_margin={strict_phase_margin}")
print("consequence=THM2396_orbitwise_N_S>=1_is_locally_sharp")
print(f"global_safe_point={GLOBAL_SAFE_POINT}; global_cover=FALSE")
print("row_decrement=0; LRC(14)=OPEN")
print("all_checks=PASS")
