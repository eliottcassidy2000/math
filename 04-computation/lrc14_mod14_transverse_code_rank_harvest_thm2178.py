#!/usr/bin/env python3
"""Exact arithmetic checks for THM-2178.

The universal theorem is proved symbolically in the markdown.  This
companion verifies its rational Fourier margin and two finite code controls.
"""

from fractions import Fraction


MODULUS = 14
ALPHABET = tuple(range(2, 13))
LENGTH = 13


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def codewords(matrix: tuple[tuple[int, ...], tuple[int, ...]]) -> set[tuple[int, ...]]:
    return {
        tuple(
            (y0 * matrix[0][index] + y1 * matrix[1][index]) % MODULUS
            for index in range(LENGTH)
        )
        for y0 in range(MODULUS)
        for y1 in range(MODULUS)
    }


def minimum_nonzero_weight(words: set[tuple[int, ...]]) -> int:
    return min(
        sum(entry != 0 for entry in word)
        for word in words
        if any(word)
    )


def allowed_kernel_count(
    matrix: tuple[tuple[int, ...], tuple[int, ...]],
) -> int:
    states: dict[tuple[int, int], int] = {(0, 0): 1}
    for index in range(LENGTH):
        next_states: dict[tuple[int, int], int] = {}
        for (first, second), count in states.items():
            for value in ALPHABET:
                target = (
                    (first + value * matrix[0][index]) % MODULUS,
                    (second + value * matrix[1][index]) % MODULUS,
                )
                next_states[target] = next_states.get(target, 0) + count
        states = next_states
    require(sum(states.values()) == len(ALPHABET) ** LENGTH, "DP lost words")
    return states.get((0, 0), 0)


def main() -> None:
    cosine_taylor = (
        Fraction(1)
        - Fraction(9, 98)
        + Fraction(81, 57624)
    )
    require(cosine_taylor == Fraction(17471, 19208), "Taylor value")
    require(cosine_taylor < Fraction(91, 100), "cosine rational cap")

    scaled_main = 11**4 * 50**4
    scaled_tail = 195 * 141**4
    fourier_margin = scaled_main - scaled_tail
    require(fourier_margin == 14431688605, "Fourier margin")
    require(fourier_margin > 0, "distance-four count did not close")

    count_lower_bound = (
        Fraction(11**9, 196)
        * Fraction(fourier_margin, 50**4)
    )
    require(
        count_lower_bound
        == Fraction(6805833364678152211, 245000000),
        "count lower bound",
    )
    require(count_lower_bound > 0, "count lower bound is not positive")

    # The terminal logarithm comparison needs only A_R >= 2.
    require(52 * 2 + 91 < 100 * 2, "Fejer terminal inequality")
    saturated_basis_cap = 26 * 105**2
    require(saturated_basis_cap == 286650, "saturated basis cap")

    affine = (
        tuple(1 for _ in range(LENGTH)),
        tuple(range(LENGTH)),
    )
    affine_words = codewords(affine)
    require(len(affine_words) == MODULUS**2, "affine code not injective")
    affine_distance = minimum_nonzero_weight(affine_words)
    require(affine_distance == 6, "affine code distance")
    affine_safe_count = allowed_kernel_count(affine)
    require(affine_safe_count == 176136387828, "affine safe count")

    low_distance = (
        (1, -1, 0) + (0,) * 10,
        (0, 1, -1) + (0,) * 10,
    )
    low_words = codewords(low_distance)
    require(len(low_words) == MODULUS**2, "low code not injective")
    low_code_distance = minimum_nonzero_weight(low_words)
    require(low_code_distance == 2, "low code distance")
    low_safe_count = allowed_kernel_count(low_distance)
    require(low_safe_count == 11**11, "low-distance safe count")

    print("THM-2178 exact mod-14 transverse-code checks")
    print(f"cosine_taylor={cosine_taylor}")
    print(f"cosine_cap_gap={Fraction(91, 100)-cosine_taylor}")
    print(f"scaled_fourier_margin={fourier_margin}")
    print(f"universal_count_lower_bound={count_lower_bound}")
    print(f"saturated_basis_cap={saturated_basis_cap}")
    print(f"affine_code_size={len(affine_words)}")
    print(f"affine_code_distance={affine_distance}")
    print(f"affine_allowed_kernel_count={affine_safe_count}")
    print(f"low_code_size={len(low_words)}")
    print(f"low_code_distance={low_code_distance}")
    print(f"low_allowed_kernel_count={low_safe_count}")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
