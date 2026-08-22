#!/usr/bin/env python3
"""Exact controls for THM-3208's positive sparse bispectrum line."""

from hashlib import sha256
from functools import lru_cache
from itertools import combinations
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


ROOT = Path(__file__).resolve().parents[1]
PINNED = {
    "04-computation/lrc14_sparse_root_bispectrum_current_thm2312.py":
        "6d80ff4460d720eff24e4339218c65bb3a21d2464dec6ffe73d5ea8bbadd8a4f",
    "04-computation/lrc14_normalized_unit_bispectrum_thm2802.py":
        "c220b84cd9574d882640fb3fa9e186a1121722256b37530349c4bb4449713341",
    "04-computation/lrc14_root_neutral_central_odd_bispectrum_thm3190.py":
        "4309bb71cbe5fd888f13a178f097826e9c0a8bf99edb566347d58b3c93840fa6",
}


def lf_hash(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


for relative, expected in PINNED.items():
    require(lf_hash(ROOT / relative) == expected,
            ("dependency hash drift", relative))


# Elements of Q(zeta_p) are represented in the basis
# 1,zeta,...,zeta^(p-2).  All controls use integer coefficients.
def monomial(prime, exponent, coefficient=1):
    exponent %= prime
    answer = [0] * (prime - 1)
    if exponent == prime - 1:
        for index in range(prime - 1):
            answer[index] -= coefficient
    else:
        answer[exponent] += coefficient
    return tuple(answer)


def add(left, right):
    return tuple(a + b for a, b in zip(left, right))


def scale(scalar, value):
    return tuple(scalar * entry for entry in value)


@lru_cache(maxsize=250000)
def multiply(prime, left, right):
    coefficients = [0] * prime
    for first, a in enumerate(left):
        for second, b in enumerate(right):
            coefficients[(first + second) % prime] += a * b
    tail = coefficients[-1]
    return tuple(coefficients[index] - tail for index in range(prime - 1))


def product(prime, *values):
    answer = monomial(prime, 0)
    for value in values:
        answer = multiply(prime, answer, value)
    return answer


@lru_cache(maxsize=None)
def conjugate(prime, value):
    answer = (0,) * (prime - 1)
    for exponent, coefficient in enumerate(value):
        answer = add(answer, monomial(prime, -exponent, coefficient))
    return answer


@lru_cache(maxsize=None)
def fourier(prime, profile, frequency):
    answer = (0,) * (prime - 1)
    for root, coefficient in enumerate(profile):
        answer = add(
            answer,
            monomial(prime, -frequency * root, coefficient),
        )
    return answer


@lru_cache(maxsize=100000)
def raw_bispectrum(prime, profile, first, second):
    return product(
        prime,
        fourier(prime, profile, first),
        fourier(prime, profile, second),
        conjugate(prime, fourier(prime, profile, first + second)),
    )


def normalized_bispectra_equal(prime, left, right):
    left_fourier = [fourier(prime, left, k) for k in range(prime)]
    right_fourier = [fourier(prime, right, k) for k in range(prime)]
    require(all(any(value) for value in left_fourier + right_fourier),
            "normalized bispectrum used outside the unit locus")
    for first in range(prime):
        for second in range(prime):
            total = (first + second) % prime
            lhs = product(
                prime,
                left_fourier[first], left_fourier[second],
                right_fourier[total], right_fourier[0],
            )
            rhs = product(
                prime,
                right_fourier[first], right_fourier[second],
                left_fourier[total], left_fourier[0],
            )
            if lhs != rhs:
                return False
    return True


def translate(profile, shift):
    size = len(profile)
    answer = [0] * size
    for root, coefficient in enumerate(profile):
        answer[(root + shift) % size] = coefficient
    return tuple(answer)


def positive_scale_translate(left, right):
    left_mass = sum(left)
    right_mass = sum(right)
    if left_mass <= 0 or right_mass <= 0:
        return False
    for shift in range(len(left)):
        moved = translate(left, shift)
        if all(right_mass * value == left_mass * target
               for value, target in zip(moved, right)):
            return True
    return False


def positive_sparse_profiles(prime, maximum_weight):
    answer = []
    for root in range(prime):
        for weight in range(1, maximum_weight + 1):
            profile = [0] * prime
            profile[root] = weight
            answer.append(tuple(profile))
    for first, second in combinations(range(prime), 2):
        for a in range(1, maximum_weight + 1):
            for b in range(1, maximum_weight + 1):
                profile = [0] * prime
                profile[first] = a
                profile[second] = b
                answer.append(tuple(profile))
    return tuple(answer)


unit_checks = 0
translation_checks = 0
whole_face_checks = 0
profile_banks = {
    5: positive_sparse_profiles(5, 3),
    7: positive_sparse_profiles(7, 2),
    13: positive_sparse_profiles(13, 2),
}

for prime, profiles in profile_banks.items():
    for profile in profiles:
        transforms = [fourier(prime, profile, k) for k in range(prime)]
        require(all(any(value) for value in transforms),
                ("positive sparse automatic-unit failure", prime, profile))
        unit_checks += prime

        # Check the exact Fourier covariance behind translation/scaling, then
        # the charged cubic law.  Normalized invariance follows by cancelling
        # these displayed factors and is also tested in the collision bank.
        for shift in range(prime):
            moved = tuple(2 * value for value in translate(profile, shift))
            for frequency in range(prime):
                expected = scale(2, multiply(
                    prime,
                    monomial(prime, -frequency * shift),
                    fourier(prime, profile, frequency),
                ))
                require(fourier(prime, moved, frequency) == expected,
                        "Fourier translation/scale covariance")
            test_pairs = ((0, 0), (1, 1), (1, 2),
                          (prime - 1, prime - 1))
            for first, second in test_pairs:
                require(
                    raw_bispectrum(prime, moved, first, second)
                    == scale(8, raw_bispectrum(
                        prime, profile, first, second)),
                    "raw charged-line scaling",
                )
            translation_checks += 1

        # Rebuild THM-2312's whole-face identity exactly in Q(zeta_p).
        whole = (0,) * (prime - 1)
        for first in range(1, prime):
            for second in range(1, prime):
                if (first + second) % prime:
                    whole = add(
                        whole,
                        raw_bispectrum(prime, profile, first, second),
                    )
        s1 = sum(profile)
        s2 = sum(value * value for value in profile)
        s3 = sum(value * value * value for value in profile)
        expected_scalar = prime * prime * s3 - 3 * prime * s1 * s2 + 2 * s1**3
        require(whole == monomial(prime, 0, expected_scalar),
                "whole-face identity")
        require(expected_scalar > 0, "positive sparse whole-face orientation")
        whole_face_checks += 1


# Independent finite collision controls for THM-2802's equality fibre.
collision_pair_checks = 0
collision_equalities = 0
for prime in [5, 7]:
    profiles = profile_banks[prime]
    for left_index, left in enumerate(profiles):
        for right in profiles[left_index:]:
            equal = normalized_bispectra_equal(prime, left, right)
            orbit = positive_scale_translate(left, right)
            require(equal == orbit,
                    ("normalized bispectrum collision fibre", prime))
            collision_pair_checks += 1
            collision_equalities += int(equal)


# Sharp boundaries.
even_uniform = (1, 1)
require(not any(fourier(2, even_uniform, 1)),
        "odd-radix hypothesis hostile")

uniform_thirteen = (1,) * 13
require(not any(fourier(13, uniform_thirteen, 1)),
        "arbitrary positive support need not be a unit")

signed = (1, -1) + (0,) * 11
require(not any(fourier(13, signed, 0)),
        "signed two-sheet zero-mode hostile")
signed_whole = (0,) * 12
for first in range(1, 13):
    for second in range(1, 13):
        if (first + second) % 13:
            signed_whole = add(
                signed_whole,
                raw_bispectrum(13, signed, first, second),
            )
require(not any(signed_whole), "signed whole-face zero hostile")

positive = (1,) + (0,) * 12
negative = (-1,) + (0,) * 12
require(normalized_bispectra_equal(13, positive, negative),
        "normalized bispectrum must forget scalar sign")
for first in range(13):
    for second in range(13):
        require(raw_bispectrum(13, negative, first, second)
                == scale(-1, raw_bispectrum(13, positive, first, second)),
                "raw bispectrum must retain scalar sign")


print("POSITIVE SPARSE UNIT-BISPECTRUM CLUTCH EXACT CONTROL")
print(f"automatic_unit_fourier_checks={unit_checks}")
print(f"positive_scale_translation_checks={translation_checks}")
print(f"whole_face_orientation_checks={whole_face_checks}")
print(f"finite_collision_pair_checks={collision_pair_checks}")
print(f"finite_collision_equalities={collision_equalities}")
print("normalized_fibre=positive_scale_times_cyclic_translation")
print("equal_mass_raw_line=positive_orientation")
print("signed_scalar_hostile=normalized_same,raw_opposite")
print("ALL EXACT CHECKS PASSED")
