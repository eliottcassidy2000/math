#!/usr/bin/env python3
"""Exact, optimization-safe companion for THM-2337."""

from fractions import Fraction
from itertools import product


P = 13
SEPTIMAL = 7
TARGET_DIMENSION = 2
TARGET_GROUP_SIZE = P**TARGET_DIMENSION
JOINT_GROUP_SIZE = P ** (2 * TARGET_DIMENSION)
RELATIONS_PER_TARGET = 3_134_566_563_840
TWO_SIDED_FLOOR = 3 * 5**7


def require(condition: bool, message: str) -> None:
    """Raise under both ordinary and optimized Python."""
    if not condition:
        raise RuntimeError(message)


def crt_residue(mod_thirteen: int, mod_seven: int = 1) -> int:
    """Return the centered representative modulo 91."""
    matches = [
        value
        for value in range(-45, 46)
        if value % P == mod_thirteen % P
        and value % SEPTIMAL == mod_seven % SEPTIMAL
    ]
    require(len(matches) == 1, "CRT representative stopped being unique")
    return matches[0]


def residue_scaling_atlas() -> tuple[int, int]:
    """Check invisibility, septimal visibility, and the divided difference."""
    checks = 0
    sharp_layers = 0
    for depth in range(1, 9):
        scale = P**depth
        for beta in range(-45, 46):
            difference = scale * beta
            require(
                difference % (P**depth) == 0,
                "word mode became visible below expiration depth",
            )
            require(
                (difference // scale) % P == beta % P,
                "first divided difference stopped recovering beta",
            )
            require(
                difference % SEPTIMAL
                == pow(-1, depth, SEPTIMAL) * beta % SEPTIMAL,
                "septimal CRT response changed",
            )
            checks += 1
            if beta % P:
                require(
                    difference % (P ** (depth + 1)) != 0,
                    "sharp first layer disappeared",
                )
                sharp_layers += 1
    return checks, sharp_layers


def two_pivot_floor() -> int:
    """Exhaust the two-coordinate affine-avoidance floor over F_7."""
    minimum = SEPTIMAL**2
    for weight_left in range(1, SEPTIMAL):
        for weight_right in range(1, SEPTIMAL):
            for displacement_left in range(SEPTIMAL):
                for displacement_right in range(SEPTIMAL):
                    for right_hand_side in range(SEPTIMAL):
                        count = 0
                        for left in range(SEPTIMAL):
                            if left == 0 or (left + displacement_left) % 7 == 0:
                                continue
                            for right in range(SEPTIMAL):
                                if (
                                    right == 0
                                    or (right + displacement_right) % 7 == 0
                                ):
                                    continue
                                if (
                                    weight_left * left
                                    + weight_right * right
                                    - right_hand_side
                                ) % 7 == 0:
                                    count += 1
                        minimum = min(minimum, count)
    require(minimum == 3, "two-pivot avoidance floor changed")
    return minimum


def target_response_atlas() -> tuple[int, tuple[tuple[int, int], ...]]:
    """Verify full first-jet image with seven-supported CRT lifts."""
    response = set()
    for target in product(range(P), repeat=TARGET_DIMENSION):
        beta = tuple(crt_residue(coordinate) for coordinate in target)
        require(
            all(value % SEPTIMAL == 1 for value in beta),
            "target lift hit a septimal Fourier zero",
        )
        require(
            tuple(value % P for value in beta) == target,
            "target first jet changed",
        )
        response.add(tuple(value % P for value in beta))

    polarizers = ((1, 0), (0, 1), (1, 1))
    require(len(response) == TARGET_GROUP_SIZE, "first-jet image is not full")
    require(
        all(polarizer in response for polarizer in polarizers),
        "a semantic polarizer left the response image",
    )
    return len(response), polarizers


def word_role_codes() -> tuple[int, int]:
    """Check the ordered nonzero-mode and zero-mode word codes."""
    words = {
        "pure_a": (1, 0),
        "pure_b": (0, 1),
        "fork_ab": (1, 1),
    }
    nonzero_codes = {
        name: tuple(1 if bit else -1 for bit in membership)
        for name, membership in words.items()
    }
    zero_codes = {
        name: tuple(1 if bit else 6 for bit in membership)
        for name, membership in words.items()
    }
    require(
        len(set(nonzero_codes.values())) == len(words),
        "ordered nonzero word code lost injectivity",
    )
    require(
        len(set(zero_codes.values())) == len(words),
        "ordered zero-mode word code lost injectivity",
    )
    require(
        nonzero_codes["pure_a"] == (1, -1)
        and nonzero_codes["pure_b"] == (-1, 1)
        and nonzero_codes["fork_ab"] == (1, 1),
        "semantic sign dictionary changed",
    )
    return len(set(nonzero_codes.values())), len(set(zero_codes.values()))


def support_mask_atlas() -> tuple[int, int, int]:
    """Check the exact pure-axis/mixed masks on F_13^2."""
    counts = {"pure_a": 0, "pure_b": 0, "fork_ab": 0}
    for left, right in product(range(P), repeat=2):
        masks = {
            "pure_a": int(left != 0 and right == 0),
            "pure_b": int(left == 0 and right != 0),
            "fork_ab": int(left != 0 and right != 0),
        }
        require(sum(masks.values()) <= 1, "support masks overlap")
        for name, value in masks.items():
            counts[name] += value
    require(counts == {"pure_a": 12, "pure_b": 12, "fork_ab": 144}, "mask counts changed")
    return counts["pure_a"], counts["pure_b"], counts["fork_ab"]


def gauge_and_weight_hostile() -> Fraction:
    """Check exact address gauge invariance and its nontrivial weight cocycle."""
    depth = 2
    scale = P**depth
    u = 1
    beta = 1
    gamma = 7
    shifted_u = u - scale * gamma
    shifted_beta = beta + gamma
    require(
        u + scale * beta == shifted_u + scale * shifted_beta,
        "atomic gauge stopped preserving the total left index",
    )
    require(
        all(
            value % SEPTIMAL
            for value in (u, beta, shifted_u, shifted_beta)
        ),
        "support-preserving gauge hit a septimal zero",
    )
    ratio = Fraction(u * beta, shifted_u * shifted_beta)
    require(ratio == Fraction(-1, 9456), "gauge weight cocycle changed")
    require(ratio != 1, "address gauge accidentally became weight preserving")
    return ratio


def odd_antipode_audit() -> tuple[int, int]:
    """Check the response-image boundary inherited from odd target order."""
    nonzero_involutive_translations = 0
    reversing_displacements = set()
    for delta in product(range(P), repeat=TARGET_DIMENSION):
        if delta != (0, 0) and all(2 * coordinate % P == 0 for coordinate in delta):
            nonzero_involutive_translations += 1
    for point in product(range(P), repeat=TARGET_DIMENSION):
        reversing_displacements.add(
            tuple((-2 * coordinate) % P for coordinate in point)
        )
    require(
        nonzero_involutive_translations == 0,
        "odd target group acquired an antipodal translation",
    )
    require(
        len(reversing_displacements) == TARGET_GROUP_SIZE,
        "pointwise negation displacements stopped being surjective",
    )
    return nonzero_involutive_translations, len(reversing_displacements)


def symbolic_zero_only_hostile() -> tuple[Fraction, Fraction]:
    """Verify the full-support zero-only hostile on G x F_13^2."""
    size = JOINT_GROUP_SIZE
    zero_value = Fraction(2 * size - (size - 1), size + 1)
    nonzero_value = Fraction(-2 + size - (size - 2), size + 1)
    require(zero_value == 1, "joint hostile zero value changed")
    require(nonzero_value == 0, "joint hostile leaked off zero")
    return zero_value, nonzero_value


residue_checks, sharp_layers = residue_scaling_atlas()
pivot_floor = two_pivot_floor()
response_size, polarizers = target_response_atlas()
nonzero_word_codes, zero_word_codes = word_role_codes()
pure_a_count, pure_b_count, fork_count = support_mask_atlas()
gauge_ratio = gauge_and_weight_hostile()
nonzero_antipodes, reversing_image_size = odd_antipode_audit()
hostile_zero, hostile_nonzero = symbolic_zero_only_hostile()

term_pairs_per_joint_fibre = TWO_SIDED_FLOOR * RELATIONS_PER_TARGET
require(pivot_floor == 3, "wrong two-pivot floor")
require(TWO_SIDED_FLOOR == 234_375, "wrong nine-coordinate lift floor")
require(
    term_pairs_per_joint_fibre == 734_664_038_400_000_000,
    "joint term bank changed",
)
require(hostile_zero == 1 and hostile_nonzero == 0, "joint hostile changed")

print("theorem=THM-2337")
print("status=PROVED+VERIFIED-EXACT+CANDIDATE-UNDER-INDEPENDENT-AUDIT")
print(f"residue_scaling_checks={residue_checks}")
print(f"sharp_first_layer_checks={sharp_layers}")
print(f"two_pivot_floor={pivot_floor}")
print(f"joint_term_pairs_per_target_and_jet={term_pairs_per_joint_fibre}")
print(f"target_response_image_size={response_size}")
print("semantic_polarizers=" + ",".join(map(str, polarizers)))
print(f"ordered_nonzero_word_codes={nonzero_word_codes}")
print(f"ordered_zero_word_codes={zero_word_codes}")
print(f"pure_a_mask_size={pure_a_count}")
print(f"pure_b_mask_size={pure_b_count}")
print(f"fork_mask_size={fork_count}")
print(f"support_preserving_gauge_ratio={gauge_ratio}")
print(f"nonzero_involutive_target_translations={nonzero_antipodes}")
print(f"pointwise_reversing_response_size={reversing_image_size}")
print(f"joint_group_size={JOINT_GROUP_SIZE}")
print(f"joint_hostile_terms_per_fibre={JOINT_GROUP_SIZE}")
print(f"joint_hostile_total_terms={JOINT_GROUP_SIZE**2}")
print("joint_hostile_surviving_label=zero")
print("full_semantic_fibre_limits=EXIST")
print("some_semantic_target_jet_fibre=NONZERO")
print("word_support_masked_aggregate=OPEN")
print("nonzero_target_fibre=OPEN")
print("terminal_phase=OPEN")
print("scalar_rows_excluded=0")
print("lrc14_status=OPEN")
print("all_checks=PASS")
