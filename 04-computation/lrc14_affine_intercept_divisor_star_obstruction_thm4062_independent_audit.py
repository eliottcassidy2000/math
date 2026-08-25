#!/usr/bin/env python3
"""Independent exact hostile audit of THM-4062."""

from fractions import Fraction
from math import gcd, lcm


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def dist_to_integer(x):
    remainder = x % 1
    return min(remainder, 1 - remainder)


def mask(speed, y_representative, fold):
    return tuple(
        label
        for label in range(fold)
        if dist_to_integer(speed * (y_representative + label) / fold)
        < Fraction(1, 14)
    )


def depth_sign_by_euclid(a, q):
    """Compute the canonical depth sign without inverse parity."""
    require(0 < a < q and gcd(a, q) == 1, "reduced point required")
    numerator, denominator = a, q
    digits = []
    while numerator:
        digit, remainder = divmod(denominator, numerator)
        digits.append(digit)
        denominator, numerator = numerator, remainder
    require(digits[-1] >= 2, "canonical terminal digit")
    return -1 if (sum(digits) - 1) % 2 else 1


def least_cyclic_period(word):
    q = len(word)
    return min(
        period
        for period in range(1, q + 1)
        if q % period == 0
        and all(word[c] == word[(c + period) % q] for c in range(q))
    )


def owner_track_signed_word(P, q):
    """Build the signed word from THM-4042's tracks owner by owner."""
    inverse_P = pow(P, -1, q)
    plus_mass = Fraction(P - q, P * q)
    minus_mass = Fraction(P - q - 1, P * q)
    word = [Fraction(0) for _ in range(q)]
    for owner in range(1, q):
        if gcd(owner, q) != 1:
            continue
        sign = depth_sign_by_euclid(owner, q)
        inverse = pow(owner, -1, q)
        plus_track = -inverse % q
        minus_track = (1 - inverse_P) * inverse % q
        word[-plus_track % q] += plus_mass * sign
        word[-minus_track % q] += minus_mass * sign
    return tuple(word)


def packet_formula_signed_word(P, q):
    multiplier = (pow(P, -1, q) - 1) % q
    plus_mass = Fraction(P - q, P * q)
    minus_mass = Fraction(P - q - 1, P * q)
    word = []
    for c in range(q):
        total = Fraction(0)
        if gcd(c, q) == 1:
            total += plus_mass * depth_sign_by_euclid(c, q)
        for unit in range(1, q):
            if gcd(unit, q) == 1 and multiplier * unit % q == c:
                total += minus_mass * depth_sign_by_euclid(unit, q)
        word.append(total)
    return tuple(word)


def unsigned_word(P, q):
    multiplier = (pow(P, -1, q) - 1) % q
    plus_mass = Fraction(P - q, P * q)
    minus_mass = Fraction(P - q - 1, P * q)
    return tuple(
        plus_mass * (gcd(c, q) == 1)
        + minus_mass
        * sum(
            gcd(unit, q) == 1 and multiplier * unit % q == c
            for unit in range(1, q)
        )
        for c in range(q)
    )


def fourier_coefficient_histogram(word, k):
    """Coefficients of sum_c word[c] Z^(kc) in Q[Z]/(Z^q-1)."""
    q = len(word)
    result = [Fraction(0) for _ in range(q)]
    for c, coefficient in enumerate(word):
        result[k * c % q] += coefficient
    return tuple(result)


def claimed_fourier_histogram(P, q, k):
    multiplier = (pow(P, -1, q) - 1) % q
    plus_mass = Fraction(P - q, P * q)
    minus_mass = Fraction(P - q - 1, P * q)
    result = [Fraction(0) for _ in range(q)]
    for a in range(1, q):
        if gcd(a, q) != 1:
            continue
        sign = depth_sign_by_euclid(a, q)
        result[k * a % q] += plus_mass * sign
        result[k * multiplier * a % q] += minus_mass * sign
    return tuple(result)


def main():
    fold = 4
    y = Fraction(1, 11)
    body = (8, 12, 20, 24, 28, 32, 36, 40)
    cover = (2, 9, 11)
    leak = (2, 1, 3)
    pair = (4, 16)

    rows = (body + cover + pair, body + leak + pair)
    for row in rows:
        require(len(row) == 13 and len(set(row)) == 13, "distinct typed row")
        require(gcd(*row) == 1, "primitive typed row")
    require(
        all(
            tuple(
                sum(speed % divisor == 0 for speed in row[:11])
                for divisor in (1, 2, 4)
            )
            == (11, 9, 8)
            for row in rows
        ),
        "common divisor incidence",
    )
    pack = tuple(sorted({speed // 4 for speed in body} | {1, 4}))
    require(pack == tuple(range(1, 11)), "divided ten-pack")
    require(
        all(
            dist_to_integer(h * (y + label)) >= Fraction(1, 11)
            for h in pack
            for label in range(fold)
        ),
        "pack lifts safe",
    )

    cover_masks = tuple(mask(speed, y, fold) for speed in cover)
    leak_masks = tuple(mask(speed, y, fold) for speed in leak)
    require(cover_masks == ((0, 2), (3,), (1,)), "cover masks")
    require(leak_masks == ((0, 2), (0,), (0,)), "leak masks")
    require(set().union(*map(set, cover_masks)) == set(range(4)), "cover union")
    require(set().union(*map(set, leak_masks)) == {0, 2}, "leak union")
    require(
        all(
            dist_to_integer(Fraction(speed * 21, 22)) >= Fraction(1, 11)
            for speed in rows[0]
        ),
        "canonical hostile lonely at 21/22",
    )
    for speed in cover + leak:
        shifted = tuple(sorted((label - 1) % fold for label in mask(speed, y, fold)))
        require(mask(speed, y + 1, fold) == shifted, "label-gauge covariance")

    P = 11
    signed_periods = []
    unsigned_periods = []
    fourier_checks = 0
    for q in range(2, P):
        owner_word = owner_track_signed_word(P, q)
        packet_word = packet_formula_signed_word(P, q)
        require(owner_word == packet_word, f"owner/packet word q={q}")
        for k in range(q):
            require(
                fourier_coefficient_histogram(owner_word, k)
                == claimed_fourier_histogram(P, q, k),
                f"formal Fourier identity q={q}, k={k}",
            )
            fourier_checks += 1
        signed_periods.append(least_cyclic_period(owner_word))
        unsigned_periods.append(least_cyclic_period(unsigned_word(P, q)))

    signed_clock = lcm(*signed_periods)
    unsigned_clock = lcm(*unsigned_periods)
    require(tuple(signed_periods) == tuple(range(2, 11)), "full signed periods")
    require(
        tuple(unsigned_periods) == (2, 3, 4, 5, 6, 7, 4, 3, 10),
        "unsigned block profile",
    )
    require(signed_clock == 2520 and unsigned_clock == 420, "P=11 clocks")

    print("status=PASS_GAUGE_REPAIRED")
    print("typed_incidence=(11,9,8);pack=1..10;strict_masks=PASS")
    print("cover_masks", cover_masks)
    print("leak_masks", leak_masks)
    print("label_gauge_covariance=mask(y+1)=mask(y)-1_mod_4")
    print("signed_periods", tuple(signed_periods), "signed_clock", signed_clock)
    print("unsigned_periods", tuple(unsigned_periods), "unsigned_clock", unsigned_clock)
    print("owner_formula_and_formal_fourier_checks", fourier_checks)


if __name__ == "__main__":
    main()
