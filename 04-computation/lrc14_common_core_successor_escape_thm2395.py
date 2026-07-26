#!/usr/bin/env python3
"""Exact companion for THM-2395.

Checks the mixed-septimal shell lemma, the carry-labelled covariant-hole
automaton, the forced escape and role-jet constants, and an explicit strict
two-cycle hostile. All theorem arithmetic is exact.
"""

from collections import defaultdict
from fractions import Fraction as F


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def frac_part(x: F) -> F:
    return x - (x.numerator // x.denominator)


def circle_norm(x: F) -> F:
    r = frac_part(x)
    return min(r, 1 - r)


def danger(speed: int, x: F) -> bool:
    return circle_norm(speed * x) < F(1, 14)


def guard_danger(speed: int, x: F) -> bool:
    return circle_norm(speed * x) < F(1, 7)


def shell(speed: int, x: F) -> bool:
    return (not danger(speed, x)) and danger(13 * speed, x)


def boundary_grid(speeds: tuple[int, ...]) -> tuple[F, ...]:
    points = {F(0), F(1)}
    for speed in speeds:
        for tooth in range(speed):
            for sign in (-1, 1):
                points.add(frac_part(F(14 * tooth + sign, 14 * speed)))
    return tuple(sorted(points))


def exact_measure(speeds: tuple[int, ...], predicate) -> F:
    boundaries = boundary_grid(speeds)
    total = F()
    for left, right in zip(boundaries, boundaries[1:]):
        if left == right:
            continue
        midpoint = (left + right) / 2
        if predicate(midpoint):
            total += right - left
    return total


def danger_roots(speed: int, base: F) -> tuple[int, ...]:
    return tuple(
        root
        for root in range(7)
        if danger(speed, (base + root) / 7)
    )


def guard_roots(speed: int, base: F) -> tuple[int, ...]:
    return tuple(
        root
        for root in range(7)
        if guard_danger(speed, (base + root) / 7)
    )


def valuation(number: int, prime: int) -> int:
    out = 0
    while number % prime == 0:
        number //= prime
        out += 1
    return out


# ---------------------------------------------------------------------------
# 1. Mixed-septimal successor shells.
# ---------------------------------------------------------------------------

same_scale_overlap = exact_measure(
    (1, 13),
    lambda x: danger(1, x) and danger(13, x),
)
shell_mass = exact_measure((1, 13), lambda x: shell(1, x))
require(same_scale_overlap == F(1, 91), "wrong D_1/D_13 overlap")
require(shell_mass == F(12, 91), "wrong successor-shell mass")

address_equal_base_mass = 7 * same_scale_overlap
require(address_equal_base_mass == F(1, 13), "wrong address-equality mass")

mixed_shell_floor = (shell_mass - address_equal_base_mass) / 7
require(mixed_shell_floor == F(5, 637), "wrong mixed shell floor")

# Two direct physical controls, neither used as a theorem premise.
control_13_7 = exact_measure(
    (13, 169, 7, 91),
    lambda x: shell(13, x) and shell(7, x),
)
control_13_14 = exact_measure(
    (13, 169, 14, 182),
    lambda x: shell(13, x) and shell(14, x),
)
require(control_13_7 == F(72, 8281), "wrong shell control (13,7)")
require(control_13_14 == F(6, 637), "wrong shell control (13,14)")
require(control_13_7 >= mixed_shell_floor, "shell control below floor")
require(control_13_14 >= mixed_shell_floor, "shell control below floor")

shell_union_cap = 2 * shell_mass - mixed_shell_floor
require(shell_union_cap == F(163, 637), "wrong shell-union cap")
require(shell_union_cap == F(2119, 8281), "wrong converted shell cap")


# ---------------------------------------------------------------------------
# 2. Exact carrier and forced-escape arithmetic.
# ---------------------------------------------------------------------------

type_ii_mass = F(3335, 8281)
retained_type_ii = type_ii_mass - shell_union_cap
type_iii_cap = F(24, 169)
escape_base = retained_type_ii - type_iii_cap
escape_physical = escape_base / 7
fixed_role_physical = escape_physical / 5

require(retained_type_ii == F(1216, 8281), "wrong retained type-II mass")
require(type_iii_cap == F(1176, 8281), "wrong type-III cap")
require(escape_base == F(40, 8281), "wrong forced escape tax")
require(escape_physical == F(40, 57967), "wrong physical escape carrier")
require(fixed_role_physical == F(8, 57967), "wrong fixed-role carrier")


# ---------------------------------------------------------------------------
# 3. Carry-labelled covariant-hole automaton.
# ---------------------------------------------------------------------------

transition_types = defaultdict(set)

for current_type in ("I", "II", "III"):
    for a in range(7):
        for b in range(7):
            for c in range(7):
                if current_type == "I":
                    valid = b == c
                    hole = b
                elif current_type == "II":
                    valid = b != c
                    hole = b
                else:
                    valid = a == c and b != a
                    hole = a
                if not valid:
                    continue

                for carry in range(7):
                    image_hole = (carry - hole) % 7
                    next_a = (carry - b) % 7
                    next_b = (carry - c) % 7

                    # Under hole covariance, classify every possible next C
                    # address by the THM-2394 trichotomy.
                    for next_c in range(7):
                        next_type = None
                        if image_hole == next_b == next_c:
                            next_type = "I"
                        elif image_hole == next_b and next_c != next_b:
                            next_type = "II"
                        elif (
                            image_hole == next_a == next_c
                            and next_b != image_hole
                        ):
                            next_type = "III"
                        if next_type is not None:
                            transition_types[current_type].add(next_type)

require(transition_types["I"] == {"I", "II"}, "wrong I transitions")
require(transition_types["II"] == {"III"}, "wrong II transitions")
require(transition_types["III"] == {"I", "II"}, "wrong III transitions")


# ---------------------------------------------------------------------------
# 4. THM-2362 role-jet constants and fixed septimal-pair refinement.
# ---------------------------------------------------------------------------

rho = escape_physical
role_sum = F(11, 13) * rho
one_mode = F(11, 156) * rho
role_energy = F(121, 2028) * rho * rho

require(role_sum == F(440, 753571), "wrong role-mode sum")
require(one_mode == F(110, 2260713), "wrong one-mode floor")
require(role_energy == F(48400, 1703607756123), "wrong role energy")

fixed_pair_base = escape_base / 42
fixed_pair_physical = fixed_pair_base / 7
fixed_pair_mode = F(11, 156) * fixed_pair_physical
septimal_role_tensor = fixed_pair_mode / 7
septimal_charged_tensor = fixed_pair_mode / 49

require(fixed_pair_base == F(20, 173901), "wrong fixed-pair base floor")
require(
    fixed_pair_physical == F(20, 1217307),
    "wrong fixed-pair physical floor",
)
require(fixed_pair_mode == F(55, 47474973), "wrong fixed-pair mode")
require(
    septimal_role_tensor == F(55, 332324811),
    "wrong septimal-role tensor floor",
)
require(
    septimal_charged_tensor == F(55, 2326273677),
    "wrong charged tensor floor",
)


# ---------------------------------------------------------------------------
# 5. Strict local two-cycle hostile.
# ---------------------------------------------------------------------------

y0 = F(1, 24)
y1 = F(13, 24)
require(frac_part(13 * y0) == y1, "first base transition")
require(frac_part(13 * y1) == y0, "second base transition")

h = 1
guard_speed = 29
lower_q = (27, 55, 83, 71)
q_top = 14
C3 = 2 * 7**2 * 13**4
c3 = 13 * C3

require(valuation(q_top, 7) == 1, "wrong q* septimal depth")
require(valuation(C3, 7) == 2, "wrong C3 septimal depth")
require(valuation(c3, 13) == 5, "wrong c3 thirteen-adic depth")

expected_guard = (5, 6)
expected_q_roots = ((1,), (2,), (3,), (4,))
expected_core = {
    y0: ((0,), (1,), (0,)),
    y1: ((6,), (0,), (6,)),
}

for base in (y0, y1):
    require(guard_roots(guard_speed, base) == expected_guard, "guard hostile")
    require(
        tuple(danger_roots(speed, base) for speed in lower_q)
        == expected_q_roots,
        "lower-q hostile",
    )
    require(danger_roots(q_top, base) == (), "q* hostile not safe")
    require(danger_roots(C3, base) == (), "C3 hostile not safe")
    require(danger_roots(c3, base) == (), "c3 hostile not safe")
    require(
        tuple(danger_roots(speed, base) for speed in (h, 13 * h, 169 * h))
        == expected_core[base],
        "core hostile addresses",
    )

    k_word = [0] * 7
    for root in expected_guard:
        k_word[root] += 1
    for word in expected_q_roots:
        k_word[word[0]] += 1
    require(k_word == [0, 1, 1, 1, 1, 1, 1], "hostile K word")

# The selected physical hole is a genuine period-two orbit.
x0 = y0 / 7
x1 = y1 / 7
require(frac_part(13 * x0) == x1, "first physical hole transition")
require(frac_part(13 * x1) == x0, "second physical hole transition")

# Every displayed inequality is strict, hence persists on an open chamber.
for base in (y0, y1):
    for root in range(7):
        x = (base + root) / 7
        for speed in lower_q + (q_top, h, 13 * h, 169 * h, C3, c3):
            require(circle_norm(speed * x) != F(1, 14), "danger endpoint")
        require(
            circle_norm(guard_speed * x) != F(1, 7),
            "guard endpoint",
        )


print("theorem=THM-2395")
print("status=PROVED-CANDIDATE+VERIFIED-EXACT; independent-audit=PENDING")
print(
    f"shell_mass={shell_mass}; mixed_overlap>={mixed_shell_floor};"
    f" union<={shell_union_cap}"
)
print(
    f"type_II_retained>={retained_type_ii};"
    f" type_III_capacity<={type_iii_cap}"
)
print(
    f"forced_escape_base>={escape_base};"
    f" physical_R_1_to_0>={escape_physical}; fixed_role>={fixed_role_physical}"
)
print(
    "covariant_automaton="
    + ";".join(
        f"{key}->{','.join(sorted(transition_types[key]))}"
        for key in ("I", "II", "III")
    )
)
print(
    f"role_sum>={role_sum}; one_mode>={one_mode}; energy>={role_energy}"
)
print(
    f"fixed_pair_base>={fixed_pair_base};"
    f" physical>={fixed_pair_physical}; mode>={fixed_pair_mode}"
)
print(
    f"F7xF13_tensor>={septimal_role_tensor};"
    f" charged_tensor>={septimal_charged_tensor}"
)
print(
    "two_cycle=(1/24,typeIII,outer-owner)<->"
    "(13/24,typeII,middle-owner); strict_open=YES"
)
print("row_decrement=0; ledger=165; LRC(14)=OPEN")
print("all_checks=PASS")
