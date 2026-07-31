#!/usr/bin/env python3
"""Exact controls for THM-2631's Wick-channel decoder no-go."""

from fractions import Fraction
from itertools import product
from math import factorial, gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def compositions(total, length, prefix=()):
    if length == 1:
        yield prefix + (total,)
        return
    for value in range(total + 1):
        yield from compositions(total - value, length - 1,
                                prefix + (value,))


def balanced_channels(exponents, moment):
    charges = tuple(a - b for a, b in exponents)
    channels = []
    for multiplicity in compositions(moment, len(exponents)):
        if sum(q * r for q, r in zip(charges, multiplicity)):
            continue
        radial = sum(a * r for (a, _b), r in zip(exponents, multiplicity))
        multinomial = factorial(moment)
        for r in multiplicity:
            multinomial //= factorial(r)
        channels.append((multiplicity, multinomial * factorial(radial)))
    return tuple(channels)


def private_rows(matrix):
    """Return one private row per column, or None."""
    rows = len(matrix)
    cols = len(matrix[0]) if rows else 0
    witnesses = []
    for col in range(cols):
        witness = next((row for row in range(rows)
                        if matrix[row][col] > 0
                        and all(matrix[row][other] == 0
                                for other in range(cols) if other != col)),
                       None)
        if witness is None:
            return None
        witnesses.append(witness)
    return tuple(witnesses)


def constructive_nonnegative_left_inverse(matrix, witnesses):
    rows = len(matrix)
    cols = len(matrix[0])
    left = [[Fraction(0) for _ in range(rows)] for _ in range(cols)]
    for col, row in enumerate(witnesses):
        left[col][row] = Fraction(1, matrix[row][col])
    product_matrix = [
        [sum(left[i][row] * matrix[row][j] for row in range(rows))
         for j in range(cols)]
        for i in range(cols)
    ]
    require(product_matrix == [
        [Fraction(int(i == j)) for j in range(cols)] for i in range(cols)
    ], "private-row decoder did not give the identity")


def binary_private_controls(rows, cols):
    total = 0
    decodable = 0
    for bits in product((0, 1), repeat=rows * cols):
        matrix = tuple(tuple(bits[row * cols + col] for col in range(cols))
                       for row in range(rows))
        total += 1
        witnesses = private_rows(matrix)
        if witnesses is not None:
            decodable += 1
            constructive_nonnegative_left_inverse(matrix, witnesses)
    return total, decodable


def main():
    # MISTAKE-211's exact first-return hostile:
    # P=a Z^6+b W^2+c W^18.
    exponents = ((6, 0), (0, 2), (0, 18))
    channels = {moment: balanced_channels(exponents, moment)
                for moment in range(1, 13)}
    counts = tuple(len(channels[moment]) for moment in range(1, 13))
    expected_four = (
        ((1, 3, 0), 4 * factorial(6)),
        ((3, 0, 1), 4 * factorial(18)),
    )
    require(channels[4] == expected_four,
            "MISTAKE-211 fourth-moment channels changed")
    require(all(not channels[moment] for moment in range(1, 4)),
            "the hostile acquired an earlier balanced return")

    # The obstruction persists at every return, not only the first one.
    # From 6 r1-2 r2-18 r3=0 and r1+r2+r3=m one obtains
    # m=4j and (r1,r2,r3)=(j+2t,3(j-t),t), 0<=t<=j.
    all_level_check = {}
    for moment in range(1, 49):
        actual = tuple(r for r, _weight
                       in balanced_channels(exponents, moment))
        if moment % 4:
            expected = ()
        else:
            j = moment // 4
            expected = tuple((j + 2 * t, 3 * (j - t), t)
                             for t in range(j + 1))
        require(actual == expected,
                f"all-level hostile parametrization failed at m={moment}")
        if expected:
            all_level_check[moment] = len(expected)

    a = Fraction(1)
    b = Fraction(1)
    c = -Fraction(factorial(6), factorial(18))
    fourth = (4 * factorial(6) * a * b**3
              + 4 * factorial(18) * a**3 * c)
    require(a and b and c and fourth == 0,
            "the two nonzero fourth-moment channels did not cancel")

    # A level-tagged moment incidence matrix has one row per degree and one
    # column per balanced channel.  Its row block at m is a single strictly
    # positive row of length |B_m|, so it has a private row exactly when that
    # level is singleton.
    nonempty_levels = tuple(moment for moment in range(1, 13)
                            if channels[moment])
    singleton_levels = tuple(moment for moment in nonempty_levels
                             if len(channels[moment]) == 1)
    collision_levels = tuple(moment for moment in nonempty_levels
                             if len(channels[moment]) > 1)
    require(4 in collision_levels and 4 not in singleton_levels,
            "the exact hostile lost its decoder obstruction")

    controls = tuple(binary_private_controls(rows, cols)
                     for rows, cols in ((2, 2), (2, 3), (3, 2), (3, 3)))
    require(controls == ((16, 2), (64, 0), (64, 18), (512, 6)),
            "binary private-row control census changed")

    # The two-channel coefficient gcd is useful only as a normalization
    # control; it does not separate the two monomials.
    coefficient_gcd = gcd(*(weight for _r, weight in channels[4]))
    require(coefficient_gcd == 4 * factorial(6),
            "fourth-moment coefficient content changed")

    print("THM-2631 exact homogeneous Wick-channel decoder controls")
    print(f"hostile_channel_counts_m1_to_m12={counts}")
    print("hostile_all_levels=m=4j has exactly j+1 channels; "
          f"checked_j1_to12={tuple(all_level_check.values())}")
    print(f"m4_channels={channels[4]} coefficient_gcd={coefficient_gcd}")
    print(f"nonempty_levels={nonempty_levels} singleton_levels={singleton_levels} collision_levels={collision_levels}")
    print(f"exact_cancellation_a_b_c={(a, b, c)} M4={fourth}")
    print(f"binary_private_decoder_controls_2x2_2x3_3x2_3x3={controls}")
    print("verdict=PASS: a colliding homogeneous moment row has no scalar linear channel decoder")
    print("SCOPE: nonlinear ideal/resultant/cumulant and whole-face Frobenius mechanisms are not excluded")


if __name__ == "__main__":
    main()
