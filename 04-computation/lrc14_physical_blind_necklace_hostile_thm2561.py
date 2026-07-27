#!/usr/bin/env python3
"""Dependency-free exact referee for THM-2561.

The script verifies a positive interval on which one genuine guard comb and
five genuine ordinary danger combs realize THM-2558's blind mask {0,1,4}.
It enumerates every exact role wall, checks the isolated k_a failure and
q_*-safe sidecar, and independently recomputes all twelve lexicographic
selected heads.
"""

from fractions import Fraction
from math import gcd


P = 13
Z0 = Fraction(1, 97)
LEFT = Fraction(13, 1281)
RIGHT = Fraction(29, 2772)
OWNER_LEFT = Fraction(4117, 399854)
OWNER_RIGHT = Fraction(4129, 399854)

# The integer ``width`` is measured in fourteenths: width 1 is an ordinary
# D_w radius 1/14 and width 2 is the guard-failure radius 1/7.
ROLES = (
    ("guard_H", 183, 2),
    ("target_k_a", 95, 1),
    ("deep_q_star_k_b", 93, 1),
    ("neutral_u_1", 114, 1),
    ("neutral_u_2", 198, 1),
    ("neutral_u_3", 304, 1),
)

BLOCKERS = (
    ("owner_c_1", 13),
    ("safe_c_2", 13**2),
    ("safe_c_3", 13**5),
)

EXPECTED_FAILURES = {
    "guard_H": frozenset({10, 11, 12}),
    "target_k_a": frozenset({3}),
    "deep_q_star_k_b": frozenset({6}),
    "neutral_u_1": frozenset({5, 9}),
    "neutral_u_2": frozenset({8}),
    "neutral_u_3": frozenset({2, 7}),
}

EXPECTED_SAFE = frozenset({0, 1, 4})
EXPECTED_HEADS = (2, 2, 7, 8, 9, 7, 7, 9, 9, 11, 2, 12)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def floor_q(value):
    return value.numerator // value.denominator


def circle_distance(value):
    fractional = value - floor_q(value)
    return min(fractional, 1 - fractional)


def failure_mask(weight, width, z):
    radius = Fraction(width, 14)
    return frozenset(
        root
        for root in range(P)
        if circle_distance(Fraction(weight, P) * (z + root)) < radius
    )


def role_walls(weight, width):
    """All z in [0,1] at which one branch lies on a role boundary."""
    radius = Fraction(width, 14)
    walls = set()
    for root in range(P):
        # On z in [0,1], weight*(z+root)/13 lies in [0,weight].
        for integer in range(weight + 1):
            for sign in (-1, 1):
                z = Fraction(P, weight) * (integer + sign * radius) - root
                if 0 <= z <= 1:
                    walls.add(z)
    return walls


def blocker_danger(blocker, z):
    require(blocker % P == 0, "blocker is not depth-positive")
    return circle_distance(Fraction(blocker, P) * z) < Fraction(1, 14)


def blocker_walls(blocker):
    """Walls in the direct base; blocker status is root-independent."""
    reduced = blocker // P
    radius = Fraction(1, 14)
    walls = set()
    for integer in range(reduced + 1):
        for sign in (-1, 1):
            z = Fraction(integer + sign * radius, reduced)
            if 0 <= z <= 1:
                walls.add(z)
    return walls


def valuation_13(value):
    valuation = 0
    while value % P == 0:
        value //= P
        valuation += 1
    return valuation


def cyclic_word(mask, slope, start):
    return tuple(
        int((start + step * slope) % P in mask) for step in range(P)
    )


def selected_head(mask, slope):
    words = [(cyclic_word(mask, slope, start), start) for start in range(P)]
    maximal_word, alpha = max(words)
    require(sum(word == maximal_word for word, _ in words) == 1,
            "lexicographic maximum is not unique")
    run = 1
    while (alpha + run * slope) % P in mask:
        run += 1
    source = (alpha + (run - 1) * slope) % P
    head = (alpha + run * slope) % P
    require(source in mask and head not in mask,
            "selected edge is not occupied-to-empty")
    return head


def main():
    weights = [weight for _, weight, _ in ROLES]
    require(len(set(weights)) == 6, "role coefficients are not distinct")
    require(all(weight > 0 and weight % P for weight in weights),
            "role coefficient is not a positive thirteen-unit")
    require(weights[0] % 2 == 1, "guard is not odd")
    common_gcd = 0
    for weight in weights:
        common_gcd = gcd(common_gcd, weight)
    require(common_gcd == 1, "six-role tuple is not primitive")
    require(len({weight % P for weight in weights}) == 6,
            "role residues are not distinct modulo thirteen")
    require(tuple(valuation_13(c) for _, c in BLOCKERS) == (1, 2, 5),
            "blocker valuation profile is not (1,2,5)")

    all_walls = {Fraction(0), Fraction(1)}
    walls_by_role = {}
    for name, weight, width in ROLES:
        walls = role_walls(weight, width)
        walls_by_role[name] = walls
        all_walls.update(walls)

    previous_wall = max(wall for wall in all_walls if wall < Z0)
    next_wall = min(wall for wall in all_walls if wall > Z0)
    require(previous_wall == LEFT, "wrong nearest left role wall")
    require(next_wall == RIGHT, "wrong nearest right role wall")
    require(not any(LEFT < wall < RIGHT for wall in all_walls),
            "a role wall crosses the asserted constant-mask interval")
    require(LEFT in walls_by_role["guard_H"],
            "left endpoint is not the guard wall")
    require(RIGHT in walls_by_role["neutral_u_2"],
            "right endpoint is not the u_2 wall")

    blocker_walls_by_role = {
        name: blocker_walls(blocker) for name, blocker in BLOCKERS
    }
    all_blocker_walls = set().union(*blocker_walls_by_role.values())
    previous_blocker_wall = max(wall for wall in all_blocker_walls if wall < Z0)
    next_blocker_wall = min(wall for wall in all_blocker_walls if wall > Z0)
    require(previous_blocker_wall == OWNER_LEFT,
            "wrong nearest left blocker wall")
    require(next_blocker_wall == OWNER_RIGHT,
            "wrong nearest right blocker wall")
    require(LEFT < OWNER_LEFT < Z0 < OWNER_RIGHT < RIGHT,
            "exclusive-owner interval left the six-comb cell")
    blocker_status = {
        name: blocker_danger(blocker, Z0) for name, blocker in BLOCKERS
    }
    require(blocker_status == {
        "owner_c_1": True,
        "safe_c_2": False,
        "safe_c_3": False,
    }, "exclusive-owner blocker word changed")

    failures = {
        name: failure_mask(weight, width, Z0)
        for name, weight, width in ROLES
    }
    require(failures == EXPECTED_FAILURES, "physical failure masks changed")
    require(sum(len(mask) for mask in failures.values()) == 10,
            "failure multiplicity is not ten")
    require(
        all(
            failures[left].isdisjoint(failures[right])
            for index, left in enumerate(failures)
            for right in tuple(failures)[index + 1:]
        ),
        "the six failure masks do not form a partition",
    )
    failure_union = frozenset().union(*failures.values())
    safe_mask = frozenset(range(P)) - failure_union
    require(safe_mask == EXPECTED_SAFE, "wrong physical A_0 mask")

    target_root = 3
    require(failures["target_k_a"] == {target_root},
            "target-active failure is not isolated")
    require(target_root not in failures["deep_q_star_k_b"],
            "q_* is not safe at the target root")
    require(
        all(
            target_root not in failures[name]
            for name in ("guard_H", "neutral_u_1", "neutral_u_2", "neutral_u_3")
        ),
        "a neutral role cofails at the target root",
    )

    heads = tuple(selected_head(safe_mask, slope) for slope in range(1, P))
    require(heads == EXPECTED_HEADS, "all-slope selected-head word changed")
    head_image = frozenset(heads)
    require(head_image == {2, 7, 8, 9, 11, 12},
            "blind-necklace head image changed")
    require(target_root not in head_image,
            "a lexicographic selector found the target-active root")

    # The target-informed chord 0 -> 3 bypasses lexicographic blindness.
    source = 0
    slope = (target_root - source) % P
    require(source in safe_mask and target_root not in safe_mask,
            "target-informed chord is not occupied-to-empty")
    require(slope == 3, "wrong target-informed chord slope")
    require(int(target_root in safe_mask) - int(source in safe_mask) == -1,
            "Cayley chord gradient is not -1")

    interval_mass = RIGHT - LEFT
    branch_mass = interval_mass / P
    owner_interval_mass = OWNER_RIGHT - OWNER_LEFT
    owner_branch_mass = owner_interval_mass / P
    require(interval_mass == Fraction(53, 169092), "wrong base interval mass")
    require(branch_mass == Fraction(53, 2198196), "wrong branch mass")
    require(owner_interval_mass == Fraction(6, 199927),
            "wrong exclusive-owner base interval mass")
    require(owner_branch_mass == Fraction(6, 2599051),
            "wrong exclusive-owner branch mass")

    print("== THM-2561: physical blind-necklace interval ==")
    print("  coefficients H;q=" + str(weights[0]) + ";" + ",".join(map(str, weights[1:])))
    print("  residues_mod_13=" + ",".join(str(weight % P) for weight in weights))
    print(f"  base_interval=({LEFT},{RIGHT})")
    print(f"  base_interval_mass={interval_mass}")
    print(f"  one_root_branch_mass={branch_mass}")
    print(f"  exact_role_walls_enumerated={sum(len(walls) for walls in walls_by_role.values())}")
    for name, _, _ in ROLES:
        print(f"  {name}_failure=" + ",".join(map(str, sorted(failures[name]))))
    print("  A_0=0,1,4")
    print("  blocker_profile=1,2,5 coefficients=13,169,371293")
    print(f"  exclusive_owner_base_interval=({OWNER_LEFT},{OWNER_RIGHT})")
    print(f"  exclusive_owner_base_mass={owner_interval_mass}")
    print(f"  exclusive_owner_one_root_branch_mass={owner_branch_mass}")
    print("  blocker_word_on_owner_interval=danger,safe,safe")
    print("  k_a_root=3 isolated=YES q_star_safe=YES")
    print("  selected_heads_tau_1_to_12=" + ",".join(map(str, heads)))
    print("  selected_head_image=2,7,8,9,11,12")
    print("  target_root_3_seen_by_any_lex_selector=NO")
    print("  target_informed_chord=0->3 slope=3 gradient=-1")
    print("all exact checks=PASS")


if __name__ == "__main__":
    main()
