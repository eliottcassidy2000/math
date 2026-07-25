#!/usr/bin/env python3
"""Exact arithmetic referee for THM-2222's scalar parity tower.

This companion freezes the finite height, profile-count consequence, and the
canonical geometric-chain and adaptive cubic-moment controls.  It does not
enumerate the enormous S_4(B) universe and therefore does not claim the open
four-checkpoint inequality.
"""

from fractions import Fraction
from itertools import product
from math import comb


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def no_three_zeros_probability(length: int) -> Fraction:
    """Stationary two-state chain probability of no 000 in `length` bits."""
    require(length >= 1, "word length must be positive")
    transition = (
        (Fraction(145, 169), Fraction(24, 169)),
        (Fraction(144, 169), Fraction(25, 169)),
    )
    # State is (last bit, terminal run of zeros), with exact mass.
    state = {
        (0, 1): Fraction(6, 7),
        (1, 0): Fraction(1, 7),
    }
    for _ in range(1, length):
        next_state: dict[tuple[int, int], Fraction] = {}
        for (bit, zero_run), mass in state.items():
            for next_bit in (0, 1):
                next_run = zero_run + 1 if next_bit == 0 else 0
                if next_run >= 3:
                    continue
                key = (next_bit, next_run)
                next_state[key] = next_state.get(key, Fraction(0)) + (
                    mass * transition[bit][next_bit]
                )
        state = next_state
    return sum(state.values(), Fraction(0))


def four_checkpoint_histogram() -> tuple[Fraction, ...]:
    """Histogram of four 3-window OR bits in a stationary 6-bit word."""
    transition = (
        (Fraction(145, 169), Fraction(24, 169)),
        (Fraction(144, 169), Fraction(25, 169)),
    )
    stationary = (Fraction(6, 7), Fraction(1, 7))
    histogram = [Fraction(0) for _ in range(5)]
    for word in product((0, 1), repeat=6):
        mass = stationary[word[0]]
        for left, right in zip(word, word[1:]):
            mass *= transition[left][right]
        checkpoint_count = sum(any(word[r : r + 3]) for r in range(4))
        histogram[checkpoint_count] += mass
    require(sum(histogram, Fraction(0)) == 1, "histogram mass drift")
    return tuple(histogram)


def main() -> None:
    delta = Fraction(961, 6930)
    height_numerator = 91**12
    height_denominator = 224 * 13**5
    b6 = height_numerator // height_denominator
    require(b6 == 3_877_322_523_365_316, "B6 arithmetic drift")

    removed_profiles = sum(comb(c - 5, 2) for c in range(7, 20))
    require(removed_profiles == comb(15, 3) == 455, "profile count drift")
    require(1140 - 4 - removed_profiles == 681, "base conditional ledger drift")
    require(1140 - 6 - removed_profiles == 679, "current conditional ledger drift")

    transition = (
        (Fraction(145, 169), Fraction(24, 169)),
        (Fraction(144, 169), Fraction(25, 169)),
    )
    stationary = (Fraction(6, 7), Fraction(1, 7))
    pushed = tuple(
        sum(stationary[i] * transition[i][j] for i in (0, 1))
        for j in (0, 1)
    )
    require(pushed == stationary, "geometric-chain stationarity drift")

    # K_s for d=(1,169,169^2) is the event that a stationary Markov word of
    # length s+2 has no occurrence of 000.
    k3 = no_three_zeros_probability(5)
    k4 = no_three_zeros_probability(6)
    require(k3 == Fraction(916159, 4826809), "K3 fraction drift")
    require(k4 == Fraction(3385513, 33787663), "K4 fraction drift")
    require(
        k3 - delta == Fraction(1710418421, 33449786370),
        "K3 gap drift",
    )
    require(
        k4 - delta == -Fraction(1286905579, 33449786370),
        "K4 gap drift",
    )
    require(Fraction(1, 7) > delta, "three-checkpoint comb control drift")

    histogram = four_checkpoint_histogram()
    require(
        histogram[1] == Fraction(127310580000, 965009442943),
        "geometric p1 drift",
    )
    require(histogram[4] == k4, "geometric p4/K4 drift")
    moments = tuple(
        sum(Fraction(comb(k, r)) * histogram[k] for k in range(5))
        for r in range(4)
    )
    naive_upper = moments[3] / 4
    adaptive_upper = min(
        moments[3],
        moments[1] - 2 * moments[2] + 3 * moments[3],
    ) / 4
    require(
        naive_upper == Fraction(864315097, 5710115047),
        "naive cubic upper drift",
    )
    require(
        adaptive_upper == Fraction(128521281793, 965009442943),
        "adaptive cubic upper drift",
    )
    require(naive_upper > delta, "naive hostile control stopped failing")
    require(adaptive_upper < delta, "adaptive positive control stopped passing")
    require(
        delta - adaptive_upper
        == Fraction(5245941691819, 955359348513570),
        "adaptive upper margin drift",
    )

    print(f"delta5={delta}")
    print(f"B6={b6}")
    print(f"lambda1_ge_6_profile_count={removed_profiles}")
    print("conditional_post_depth3_ledger_before_depth4_closures=681")
    print("conditional_current_ledger_after_thm2213_thm2215=679")
    print("geometric_chain_transition=[[145,24],[144,25]]/169")
    print("geometric_chain_stationary=(6/7,1/7)")
    print(f"geometric_chain_K3={k3}")
    print(f"geometric_chain_K3_minus_delta={k3-delta}")
    print(f"geometric_chain_K4={k4}")
    print(f"geometric_chain_K4_minus_delta={k4-delta}")
    print(f"geometric_four_checkpoint_p1={histogram[1]}")
    print(f"geometric_naive_cubic_upper={naive_upper}")
    print(f"geometric_adaptive_cubic_upper={adaptive_upper}")
    print(f"delta_minus_adaptive_upper={delta-adaptive_upper}")
    print("status=THM-2222_EXACT_ARITHMETIC_REFEREE_OPEN_S4")


if __name__ == "__main__":
    main()
