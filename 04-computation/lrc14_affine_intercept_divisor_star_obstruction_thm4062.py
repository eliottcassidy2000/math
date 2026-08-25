#!/usr/bin/env python3
"""Exact probes for the THM-4059 to LRC(14) complement bridge."""

from fractions import Fraction
from math import gcd, lcm


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def circle_distance(x):
    floor = x.numerator // x.denominator
    fraction = x - floor
    return min(fraction, 1 - fraction)


def danger_mask(delta, y_representative, d):
    return tuple(
        label
        for label in range(d)
        if circle_distance(
            Fraction(delta) * (y_representative + label) / d
        )
        < Fraction(1, 14)
    )


def packet(delta, d):
    common = gcd(delta, d)
    q = d // common
    return q, (delta // common) % q if q > 1 else 0


def stern_sign(a, q):
    require(q > 1 and gcd(a, q) == 1, "unit packet required")
    inverse = pow(a, -1, q)
    if q % 2:
        return 1 if (a + inverse) % 2 == 0 else -1
    carry = (a * inverse - 1) // q
    return 1 if (carry + 1) % 2 == 0 else -1


def signed_owner_word(P, q):
    """Depth-signed analogue of THM-4042 equation (3), q>=2."""
    require(2 <= q < P and all(P % p for p in range(2, P)), "prime P")
    multiplier = (pow(P, -1, q) - 1) % q
    plus = Fraction(P - q, P * q)
    minus = Fraction(P - q - 1, P * q)
    word = []
    for c in range(q):
        value = Fraction(0)
        if gcd(c, q) == 1:
            value += plus * stern_sign(c, q)
        for unit in range(1, q):
            if gcd(unit, q) == 1 and multiplier * unit % q == c:
                value += minus * stern_sign(unit, q)
        word.append(value)
    return tuple(word)


def bookkeeping_word(P):
    """THM-4059's declared +1 zero-phase weight, not a Stern depth."""
    return (Fraction(2 * P - 3, P),)


def least_period(word):
    q = len(word)
    return min(
        period
        for period in range(1, q + 1)
        if q % period == 0
        and all(word[i] == word[(i + period) % q] for i in range(q))
    )


def signed_phase_vector(P, phase):
    vector = [Fraction(0) for _ in range(P - 1)]
    vector[0] += bookkeeping_word(P)[0]
    for q in range(2, P):
        word = signed_owner_word(P, q)
        for c in range(q):
            vector[c] += word[(c - phase) % q]
    return tuple(vector)


def main():
    d = 4
    y = Fraction(1, 11)
    covering = (2, 9, 11)
    leaking = (2, 1, 3)
    packets_covering = tuple(packet(delta, d) for delta in covering)
    packets_leaking = tuple(packet(delta, d) for delta in leaking)
    masks_covering = tuple(danger_mask(delta, y, d) for delta in covering)
    masks_leaking = tuple(danger_mask(delta, y, d) for delta in leaking)
    union_covering = set().union(*map(set, masks_covering))
    union_leaking = set().union(*map(set, masks_leaking))

    require(packets_covering == packets_leaking, "packet aliases must agree")
    require(union_covering == set(range(d)), "first triple must cover C4")
    require(union_leaking != set(range(d)), "second triple must leak")
    require(
        tuple(stern_sign(a, q) for q, a in packets_covering) == (-1, -1, -1),
        "small equality packets have the expected constant depth sign",
    )
    for delta in covering + leaking:
        rotated = tuple(sorted((label - 1) % d for label in danger_mask(delta, y, d)))
        require(
            danger_mask(delta, y + 1, d) == rotated,
            "representative change must rotate the label gauge",
        )

    print("d4_packet_alias")
    print("phase_representative", y)
    print("covering_exceptions", covering)
    print("leaking_exceptions", leaking)
    print("common_packets", packets_covering)
    print("common_depth_signs", (-1, -1, -1))
    print("covering_masks", masks_covering, "union", tuple(sorted(union_covering)))
    print("leaking_masks", masks_leaking, "union", tuple(sorted(union_leaking)))
    print("missing_labels", tuple(sorted(set(range(d)) - union_leaking)))
    print("label_covariance", "D(y+1)=D(y)-1_mod_d")

    P = 11
    periods = {q: least_period(signed_owner_word(P, q)) for q in range(2, P)}
    require(periods[5] == 5 and periods[7] == 7, "odd-prime periods")
    require(periods[8] == 8 and periods[9] == 9, "restored prime-power depths")
    forced = lcm(periods[5], periods[7], periods[8], periods[9])
    require(forced == 2520, "signed P=11 lower bound")
    phase_vectors = tuple(signed_phase_vector(P, phase) for phase in range(forced))
    proper_divisors = tuple(divisor for divisor in range(1, forced) if forced % divisor == 0)
    require(
        all(
            any(
                phase_vectors[phase]
                != phase_vectors[(phase + divisor) % forced]
                for phase in range(forced)
            )
            for divisor in proper_divisors
        ),
        "every proper divisor shift must fail somewhere",
    )
    require(len(set(phase_vectors)) == forced, "all signed phase vectors distinct")
    print("p11_signed_owner_clock")
    print("bookkeeping_q1", bookkeeping_word(P), "not_a_stern_depth")
    print("block_periods", tuple(periods.items()))
    print("q8_word", signed_owner_word(P, 8))
    print("q9_word", signed_owner_word(P, 9))
    print("forced_period", forced)
    print("proper_divisor_shifts_rejected", len(proper_divisors))
    print("phase_vector_count", len(set(phase_vectors)))
    print("unsigned_period_from_THM4042", 420)


if __name__ == "__main__":
    main()
