#!/usr/bin/env python3
"""Exact arithmetic referee for THM-2223.

The theorem's analytic input is the cited Vaaler upper polynomial from
THM-2085.  This companion checks all load-bearing rational arithmetic,
certificate minimality, balanced-base no-carry inequalities, and template
counts.  No Python ``assert`` is load bearing under ``python -O``.
"""

from fractions import Fraction


CHECKPOINT_BASE = 169
CHECKPOINTS = 4
BLOCKERS = 3
DIGIT_HEIGHT = 16
SCALAR_FLOOR = Fraction(961, 6930)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relation_free_bound(degree: int) -> Fraction:
    """The 3^4-term positive Selberg upper-product constant."""
    upper_constant = Fraction(1, 7) + Fraction(1, degree + 1)
    return BLOCKERS**CHECKPOINTS * upper_constant**CHECKPOINTS


def run() -> None:
    first_degree = next(
        degree
        for degree in range(1, 100)
        if relation_free_bound(degree) < SCALAR_FLOOR
    )
    require(first_degree == DIGIT_HEIGHT, "first degree changed")

    degree_15 = relation_free_bound(15)
    degree_16 = relation_free_bound(16)
    gap = SCALAR_FLOOR - degree_16
    require(
        degree_15 == Fraction(22667121, 157351936),
        "degree-fifteen bound changed",
    )
    require(
        degree_16 == Fraction(26873856, 200533921),
        "degree-sixteen bound changed",
    )
    require(
        gap == Fraction(925325143, 198528581790) and gap > 0,
        "scalar-floor gap changed",
    )

    place_sum = sum(
        CHECKPOINT_BASE**place for place in range(CHECKPOINTS)
    )
    coefficient_height = DIGIT_HEIGHT * place_sum
    require(place_sum == 4855540, "base-169 place sum changed")
    require(
        coefficient_height == 77688640,
        "coefficient-height invoice changed",
    )

    no_carry_margins = tuple(
        CHECKPOINT_BASE**highest
        - DIGIT_HEIGHT
        * sum(CHECKPOINT_BASE**place for place in range(highest))
        for highest in range(1, CHECKPOINTS)
    )
    require(
        no_carry_margins == (153, 25841, 4367113),
        "balanced-base no-carry margins changed",
    )
    require(all(value > 0 for value in no_carry_margins), "carry collision")

    options_per_checkpoint = 1 + BLOCKERS * 2 * DIGIT_HEIGHT
    directed_templates = options_per_checkpoint**CHECKPOINTS - 1
    require(options_per_checkpoint == 97, "checkpoint alphabet changed")
    require(directed_templates == 88529280, "template count changed")
    require(directed_templates % 2 == 0, "sign quotient not free")
    unoriented_templates = directed_templates // 2
    require(unoriented_templates == 44264640, "sign quotient changed")

    forbidden_profiles = tuple(
        (shallow, shallow, deepest)
        for shallow in range(6, 20)
        for deepest in range(shallow + 1, 20)
        if deepest - shallow >= 8
    )
    require(len(forbidden_profiles) == 21, "profile count changed")
    require(
        forbidden_profiles[0] == (6, 6, 14)
        and forbidden_profiles[-1] == (11, 11, 19),
        "profile boundary changed",
    )
    valuation_gap_spectrum = tuple(
        sorted(
            {
                abs(2 * (right - left) + right_epsilon - left_epsilon)
                for left in range(CHECKPOINTS)
                for right in range(CHECKPOINTS)
                if left != right
                for left_epsilon in range(2)
                for right_epsilon in range(2)
            }
        )
    )
    require(
        valuation_gap_spectrum == tuple(range(1, 8)),
        "valuation-gap spectrum changed",
    )

    print("FOUR-CHECKPOINT SELBERG CARRY GATE")
    print(f"scalar_floor={SCALAR_FLOOR}")
    print(f"first_successful_degree={first_degree}")
    print(f"degree15_relation_free_bound={degree_15}")
    print(f"degree16_relation_free_bound={degree_16}")
    print(f"degree16_floor_gap={gap}")
    print(f"checkpoint_base={CHECKPOINT_BASE}")
    print(f"place_sum={place_sum}")
    print(f"grouped_coefficient_height={coefficient_height}")
    print(f"no_carry_margins={no_carry_margins}")
    print(f"options_per_checkpoint={options_per_checkpoint}")
    print(f"directed_nonzero_templates={directed_templates}")
    print(f"unoriented_templates={unoriented_templates}")
    print(f"valuation_gap_spectrum={valuation_gap_spectrum}")
    print(f"excluded_equal-shallow_profiles={len(forbidden_profiles)}")
    print(f"profile_boundary={forbidden_profiles[0]}..{forbidden_profiles[-1]}")
    print("status=THM-2223_EXACT_ARITHMETIC_REFEREE")


if __name__ == "__main__":
    run()
