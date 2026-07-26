#!/usr/bin/env python3
"""Exact finite atlas for THM-2448."""

from __future__ import annotations

from fractions import Fraction


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def first_unit_mismatch(left: int, right: int) -> int | None:
    """Cell 0 is all-safe; cells 1..5 fail first at that coordinate."""

    if left == right:
        return None
    if left == 0:
        return right
    if right == 0:
        return left
    return min(left, right)


def first_blocker_mismatch(
    left: tuple[int, int], right: tuple[int, int]
) -> int | None:
    for index, (a, b) in enumerate(zip(left, right), start=1):
        if a != b:
            return index
    return None


blocker_words = tuple(
    (source, target) for source in range(2) for target in range(2)
)
cells = tuple(
    (unit_cell, blocker_word)
    for unit_cell in range(6)
    for blocker_word in blocker_words
)
require(len(cells) == 24, "generic endpoint atlas changed size")

matched = 0
unit_transitions = [0] * 5
blocker_transitions = [0] * 2
for left in cells:
    for right in cells:
        if left == right:
            matched += 1
            continue
        unit = first_unit_mismatch(left[0], right[0])
        if unit is not None:
            unit_transitions[unit - 1] += 1
            continue
        blocker = first_blocker_mismatch(left[1], right[1])
        require(blocker is not None, "unequal labels lost their first mismatch")
        blocker_transitions[blocker - 1] += 1

require(matched == 24, "matched-cell census changed")
require(
    unit_transitions == [160, 128, 96, 64, 32],
    "guard/unit transition census changed",
)
require(
    blocker_transitions == [48, 24],
    "blocker transition census changed",
)
require(
    matched + sum(unit_transitions) + sum(blocker_transitions) == 24**2,
    "generic atlas did not partition",
)


# A matched first-failure cell i>0 leaves the later unit factors
# unspecified.  Split each omitted factor on both endpoint copies and
# stop as soon as an off-diagonal pair occurs.  If n factors remain,
# there are two off-diagonal terminal children and two diagonal children
# which continue:
#
#     T(0)=1,  T(n)=2+2T(n-1)=3*2^n-2.
#
# Thus the worst first-failure label is i=1, with four tail factors and
# 46 terminal cospan pieces.  Cell 0 is already all-safe.
def tail_terminal_count(number_of_tail_factors: int) -> int:
    require(number_of_tail_factors >= 0, "negative tail length")
    value = 1
    for _ in range(number_of_tail_factors):
        value = 2 + 2 * value
    require(
        value == 3 * 2**number_of_tail_factors - 2,
        "tail recurrence lost its closed form",
    )
    return value


tail_terminals_by_unit_cell = tuple(
    1 if unit_cell == 0 else tail_terminal_count(5 - unit_cell)
    for unit_cell in range(6)
)
require(
    tail_terminals_by_unit_cell == (1, 46, 22, 10, 4, 1),
    "tail-completion atlas changed",
)
# For a fixed left cell, 23 of the 24 initial right cells are already
# off-diagonal.  Only its single matched right cell is recursively
# refined.  Hence the complete decomposition has 23+T(n), not 24*T(n),
# terminal currents.
completed_denominators_by_unit_cell = tuple(
    23 + count for count in tail_terminals_by_unit_cell
)
require(
    completed_denominators_by_unit_cell == (24, 69, 45, 33, 27, 24),
    "completed fixed-left denominators changed",
)
generic_complete_uniform_denominator = max(completed_denominators_by_unit_cell)
require(
    generic_complete_uniform_denominator == 69,
    "generic completed-cospan denominator changed",
)


# In the source-completed ghost the source factor is already inserted on
# both endpoint copies.  Only five guard/unit roles and the other blocker
# remain: six first-failure labels times two blocker states.
ghost_cells = tuple(
    (unit_cell, target_blocker) for unit_cell in range(6) for target_blocker in range(2)
)
ghost_left = (0, 0)
ghost_matched = 0
ghost_unit_transitions = [0] * 5
ghost_blocker_transitions = 0
for right in ghost_cells:
    if right == ghost_left:
        ghost_matched += 1
    elif right[0] != 0:
        ghost_unit_transitions[right[0] - 1] += 1
    else:
        require(right == (0, 1), "ghost transition lost its blocker label")
        ghost_blocker_transitions += 1

require(len(ghost_cells) == 12, "ghost endpoint atlas changed size")
require(ghost_matched == 1, "ghost matched endpoint ceased to be unique")
require(
    ghost_unit_transitions == [2, 2, 2, 2, 2],
    "ghost unit transition census changed",
)
require(ghost_blocker_transitions == 1, "ghost blocker transition changed")


# Triangle bounds are sharp: all components may equal J/N.
generic_parts = (Fraction(1, 24),) * 24
ghost_parts = (Fraction(1, 12),) * 12
completed_parts = (Fraction(1, 69),) * 69
require(sum(generic_parts) == 1, "generic cospan conservation failed")
require(max(generic_parts) == Fraction(1, 24), "generic 1/24 sharpness failed")
require(sum(ghost_parts) == 1, "ghost cospan conservation failed")
require(max(ghost_parts) == Fraction(1, 12), "ghost 1/12 sharpness failed")
require(sum(completed_parts) == 1, "completed cospan conservation failed")
require(
    max(completed_parts) == Fraction(1, 69),
    "completed 1/69 sharpness failed",
)


# A matched term need not survive.  The whole coefficient can be carried
# by one off-diagonal transition.
two_cell_hostile = (Fraction(0), Fraction(1))
require(sum(two_cell_hostile) == 1, "off-diagonal hostile lost mass")
require(two_cell_hostile[0] == 0, "matched hostile unexpectedly survived")


print("THM-2448 exact companion")
print(f"generic_cells={len(cells)}")
print(f"generic_ordered_pairs={len(cells) ** 2}")
print(f"generic_matched={matched}")
print("generic_unit_transitions=" + ",".join(map(str, unit_transitions)))
print("generic_blocker_transitions=" + ",".join(map(str, blocker_transitions)))
print(f"generic_transition_total={sum(unit_transitions) + sum(blocker_transitions)}")
print(
    "tail_terminals_by_unit_cell="
    + ",".join(map(str, tail_terminals_by_unit_cell))
)
print(
    "completed_denominators_by_unit_cell="
    + ",".join(map(str, completed_denominators_by_unit_cell))
)
print(
    "generic_complete_uniform_denominator="
    + str(generic_complete_uniform_denominator)
)
print(f"ghost_cells={len(ghost_cells)}")
print(f"ghost_matched={ghost_matched}")
print("ghost_unit_transitions=" + ",".join(map(str, ghost_unit_transitions)))
print(f"ghost_blocker_transitions={ghost_blocker_transitions}")
print("generic_sharp_fraction=1/24")
print("completed_sharp_fraction=1/69")
print("ghost_sharp_fraction=1/12")
print("matched_zero_hostile=PASS")
print("ALL CHECKS PASSED")
