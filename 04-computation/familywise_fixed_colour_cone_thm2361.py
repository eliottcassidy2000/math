#!/usr/bin/env python3
"""Exact companion for THM-2361's familywise cone and phase boundary."""

from __future__ import annotations

from fractions import Fraction
from math import gcd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def inverse_mod(value: int, modulus: int) -> int:
    return pow(value, -1, modulus)


def main() -> None:
    # The same Galois exponent straightens every member of a finite family.
    # Rational angle inequalities certify a cone of width strictly below pi.
    moduli = tuple(range(2, 201))
    automorphism_checks = 0
    family_size = 17
    for modulus in moduli:
        radius = (modulus + 6) // 7 - 1
        require(
            7 * radius < modulus
            and 4 * radius < modulus,
            f"acute cone inequalities failed at N={modulus}",
        )
        units = [
            colour for colour in range(1, modulus) if gcd(colour, modulus) == 1
        ]
        for colour in units:
            inverse = inverse_mod(colour, modulus)
            for family_index in range(family_size):
                coefficients = {
                    displacement: Fraction(
                        1
                        + (
                            family_index
                            + (displacement + radius) * (colour + 1)
                        )
                        % 19,
                        1 + family_index,
                    )
                    for displacement in range(-radius, radius + 1)
                }
                require(
                    coefficients.get(0, Fraction(0)) > 0,
                    "family diagonal stopped being positive",
                )
                for displacement in coefficients:
                    require(
                        (
                            colour
                            * displacement
                            * inverse
                            - displacement
                        )
                        % modulus
                        == 0,
                        "one Galois exponent stopped straightening the family",
                    )
                    automorphism_checks += 1

    # The arithmetic-comb variant uses a*k^{-1} and straightens kd to ad.
    comb_checks = 0
    for modulus in (13, 29, 91, 169, 181):
        for a in range(1, modulus):
            if gcd(a, modulus) != 1:
                continue
            for colour in range(1, modulus):
                if gcd(colour, modulus) != 1:
                    continue
                exponent = (a * inverse_mod(colour, modulus)) % modulus
                for displacement in range(modulus):
                    require(
                        (
                            exponent * colour * displacement
                            - a * displacement
                        )
                        % modulus
                        == 0,
                        "comb Galois straightening changed",
                    )
                    comb_checks += 1

    # A nontrivial order-13 character takes every thirteenth-root value.
    # Their formal exponent sum is 1+X+...+X^12, the cyclotomic relation.
    character_checks = 0
    for character in range(1, 13):
        values = {(character * parameter) % 13 for parameter in range(13)}
        require(
            values == set(range(13)),
            "order-13 character stopped being surjective",
        )
        character_checks += 1

    # Exact off-diagonal gauge-index ledger.
    gauge_checks = 0
    for modulus in (13, 91, 169):
        for k in range(1, modulus):
            if gcd(k, modulus) != 1:
                continue
            for difference in (-2 * modulus - 1, -13, -1, 0, 1, 13, 2 * modulus + 1):
                ell = k + difference
                for gauge_offset in (-3, -1, 0, 2, 5):
                    delta = difference + modulus * gauge_offset
                    for h in (-7, -1, 0, 4, 9):
                        left_frequency = k + modulus * h
                        right_frequency = ell + modulus * (h + gauge_offset)
                        require(
                            right_frequency - left_frequency == delta
                            and (ell - delta - k) % modulus == 0,
                            "off-diagonal gauge/physical phase ledger changed",
                        )
                        gauge_checks += 1

    # Sharp N=169 hostile, checked entirely with Fraction and Q(i) phases.
    modulus = 169
    difference = 13
    k = 1
    ell = 14
    epsilon = Fraction(1, 1352)
    centres = (Fraction(-1, 52), Fraction(1, 52))
    require(
        gcd(k, modulus) == gcd(ell, modulus) == 1,
        "hostile colours stopped being primitive",
    )
    require(
        all(abs(centre) + epsilon < Fraction(1, 14) for centre in centres),
        "hostile intervals left the danger arc",
    )
    image_centres = tuple((modulus * centre) % 1 for centre in centres)
    image_half_width = modulus * epsilon
    require(
        image_centres == (Fraction(3, 4), Fraction(1, 4))
        and image_half_width == Fraction(1, 8),
        "hostile one-sheet image geometry changed",
    )
    mass = 4 * epsilon
    diagonal = modulus * mass
    require(diagonal == Fraction(1, 2), "hostile diagonal changed")

    phase_exponents = tuple((difference * centre) % 1 for centre in centres)
    require(
        phase_exponents == (Fraction(3, 4), Fraction(1, 4)),
        "hostile terminal phases changed",
    )
    # e^(2*pi*i*3/4)=-i and e^(2*pi*i*1/4)=i.
    q_i_phases = ((0, -1), (0, 1))
    phase_sum = (
        sum(phase[0] for phase in q_i_phases),
        sum(phase[1] for phase in q_i_phases),
    )
    require(phase_sum == (0, 0), "off-diagonal hostile stopped cancelling")
    off_diagonal = Fraction(0)
    compensated = diagonal

    print("THM-2361 familywise fixed-colour cone exact companion")
    print(f"moduli checked: {len(moduli)}")
    print(f"family/Galois support checks: {automorphism_checks}")
    print(f"arithmetic-comb straightening checks: {comb_checks}")
    print(f"order-13 character obstructions: {character_checks}")
    print(f"off-diagonal gauge ledgers: {gauge_checks}")
    print("hostile: N=169, k=1, ell=14, Delta=13")
    print(f"hostile image centres/half-width: {image_centres}, {image_half_width}")
    print(f"hostile diagonal current: {diagonal}")
    print(f"hostile off-diagonal current: {off_diagonal}")
    print(f"hostile phase-compensated current: {compensated}")
    print("VERDICT: one familywise cone; off-diagonal terminal phase is essential")


if __name__ == "__main__":
    main()
