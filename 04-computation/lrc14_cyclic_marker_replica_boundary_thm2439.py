#!/usr/bin/env python3
"""Exact companion for THM-2439.

This dependency-free audit computes the Boolean Moebius degree of the
translation-equivariant C13 lexicographic marker and verifies the first
non-dihedral homometric self-Gram hostile.
"""

from __future__ import annotations

from collections import defaultdict
from functools import lru_cache


P = 13
FULL = (1 << P) - 1


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


@lru_cache(maxsize=None)
def indicator(mask: int) -> tuple[int, ...]:
    return tuple((mask >> root) & 1 for root in range(P))


@lru_cache(maxsize=None)
def translate(mask: int, shift: int) -> int:
    result = 0
    for root in range(P):
        if (mask >> root) & 1:
            result |= 1 << ((root + shift) % P)
    return result


@lru_cache(maxsize=None)
def reflect(mask: int) -> int:
    result = 0
    for root in range(P):
        if (mask >> root) & 1:
            result |= 1 << ((-root) % P)
    return result


@lru_cache(maxsize=None)
def word(mask: int, start: int) -> tuple[int, ...]:
    bits = indicator(mask)
    return tuple(bits[(start + offset) % P] for offset in range(P))


@lru_cache(maxsize=None)
def marker(mask: int) -> int:
    require(mask not in (0, FULL), "marker called on a constant mask")
    words = tuple(word(mask, start) for start in range(P))
    maximum = max(words)
    winners = tuple(start for start, candidate in enumerate(words) if candidate == maximum)
    require(len(winners) == 1, "nonconstant prime-cyclic mask has a marker tie")
    require((mask >> winners[0]) & 1, "marker did not select a live root")
    return winners[0]


@lru_cache(maxsize=None)
def marker_zero(mask: int) -> int:
    if mask in (0, FULL):
        return 0
    return int(marker(mask) == 0)


def submasks(mask: int):
    submask = mask
    while True:
        yield submask
        if submask == 0:
            break
        submask = (submask - 1) & mask


# ---------------------------------------------------------------------------
# 1. Exact Boolean Moebius expansion.
# ---------------------------------------------------------------------------

coefficients: dict[int, int] = {}
for support in range(1 << P):
    coefficient = 0
    support_size = support.bit_count()
    for subset in submasks(support):
        sign = -1 if (support_size - subset.bit_count()) % 2 else 1
        coefficient += sign * marker_zero(subset)
    coefficients[support] = coefficient


def evaluate(mask: int, degree_cap: int) -> int:
    return sum(
        coefficients[support]
        for support in submasks(mask)
        if support.bit_count() <= degree_cap
    )


restricted_checks = 0
for mask in range(1 << P):
    if 1 <= mask.bit_count() <= 10:
        require(evaluate(mask, 10) == marker_zero(mask), "degree-ten formula failed")
        restricted_checks += 1

degree_ten_support = sum(1 << root for root in range(10))
degree_ten_zero_support = sum(1 << root for root in (0, 2, 3, 4, 5, 6, 7, 8, 9, 10))
require(coefficients[degree_ten_support] == 3, "degree-ten witness changed")
require(coefficients[degree_ten_zero_support] == 0, "degree-ten zero witness changed")
restricted_nonzero_coefficients = sum(
    coefficient != 0
    for support, coefficient in coefficients.items()
    if support.bit_count() <= 10
)
restricted_degree = max(
    support.bit_count()
    for support, coefficient in coefficients.items()
    if coefficient != 0 and support.bit_count() <= 10
)
require(restricted_nonzero_coefficients == 2753, "restricted coefficient count changed")
require(restricted_degree == 10, "restricted replica degree changed")

universal_checks = 0
for mask in range(1 << P):
    require(evaluate(mask, 13) == marker_zero(mask), "full Moebius formula failed")
    universal_checks += 1

universal_nonzero_coefficients = sum(value != 0 for value in coefficients.values())
universal_degree = max(
    support.bit_count() for support, value in coefficients.items() if value != 0
)
nonzero_coefficients_by_degree = tuple(
    sum(
        value != 0 and support.bit_count() == degree
        for support, value in coefficients.items()
    )
    for degree in range(1, 13)
)
marker_truths_by_size = tuple(
    sum(marker_zero(mask) for mask in range(1 << P) if mask.bit_count() == size)
    for size in range(1, 13)
)
require(universal_nonzero_coefficients == 2816, "universal coefficient count changed")
require(universal_degree == 12, "universal marker degree changed")
require(
    nonzero_coefficients_by_degree
    == (1, 6, 34, 119, 284, 493, 634, 594, 401, 187, 54, 9),
    "degreewise coefficient census changed",
)
require(
    marker_truths_by_size == (1, 6, 22, 55, 99, 132, 132, 99, 55, 22, 6, 1),
    "marker truth census changed",
)
require(coefficients[FULL] == 0, "degree-thirteen coefficient became nonzero")
require(coefficients[FULL ^ (1 << 6)] == -10, "first degree-twelve witness changed")
require(coefficients[FULL ^ (1 << 10)] == 2, "second degree-twelve witness changed")

equivariance_checks = 0
for mask in range(1, FULL):
    selected = marker(mask)
    for shift in range(P):
        require(
            marker(translate(mask, shift)) == (selected + shift) % P,
            "marker lost translation equivariance",
        )
        equivariance_checks += 1
require(equivariance_checks == 106470, "wrong marker-equivariance universe")


# ---------------------------------------------------------------------------
# 2. Homometric self-Gram hostile.
# ---------------------------------------------------------------------------

def mask_from_set(values: set[int]) -> int:
    return sum(1 << value for value in values)


@lru_cache(maxsize=None)
def autocorrelation(mask: int) -> tuple[int, ...]:
    values = {root for root in range(P) if (mask >> root) & 1}
    return tuple(
        sum(((root + shift) % P) in values for root in values)
        for shift in range(P)
    )


@lru_cache(maxsize=None)
def dihedral_orbit(mask: int) -> frozenset[int]:
    reflected = reflect(mask)
    return frozenset(
        [translate(mask, shift) for shift in range(P)]
        + [translate(reflected, shift) for shift in range(P)]
    )


A = mask_from_set({0, 1, 3, 9})
B = mask_from_set({1, 2, 5, 7})
expected_correlation = (4,) + (1,) * 12

require(autocorrelation(A) == expected_correlation, "A correlation changed")
require(autocorrelation(B) == expected_correlation, "B correlation changed")
require(B not in dihedral_orbit(A), "hostiles became dihedrally equivalent")
require(marker(A) == 0 and marker(B) == 1, "hostile markers changed")

A_marked_word = "".join(map(str, word(A, marker(A))))
B_marked_word = "".join(map(str, word(B, marker(B))))
require(A_marked_word == "1101000001000", "A marked word changed")
require(B_marked_word == "1100101000000", "B marked word changed")
require(A_marked_word != B_marked_word, "marked profiles became equal")
affine_equivalences = tuple(
    (unit, shift)
    for unit in range(1, P)
    for shift in range(P)
    if {
        (unit * root + shift) % P
        for root in range(P)
        if (A >> root) & 1
    }
    == {root for root in range(P) if (B >> root) & 1}
)
require(
    affine_equivalences == ((7, 7), (8, 7), (11, 7)),
    "affine-equivalence caveat changed",
)

# Correlation determines both the cyclic self-Gram distance and the
# exact Fourier power. Reduce its character polynomial modulo Phi_13
# without floating point.
translation_gram = tuple(
    2 * (A.bit_count() - overlap) for overlap in autocorrelation(A)
)
require(translation_gram == (0,) + (6,) * 12, "translation-Gram hostile changed")
computed_power_numerators = []
for character in range(P):
    if character == 0:
        computed_power_numerators.append(sum(autocorrelation(A)))
        continue
    permuted = [0] * P
    for shift, coefficient in enumerate(autocorrelation(A)):
        permuted[(character * shift) % P] = coefficient
    reduced = tuple(permuted[index] - permuted[12] for index in range(12))
    require(reduced[1:] == (0,) * 11, "power did not reduce to a rational")
    computed_power_numerators.append(reduced[0])
fourier_power_numerators = tuple(computed_power_numerators)
require(
    fourier_power_numerators == (16,) + (3,) * 12,
    "Fourier-power hostile changed",
)

# Verify that support sizes one through three have no autocorrelation
# class containing more than one dihedral orbit, while size four does.
first_homometric_support = None
homometric_signature_count_at_four = 0
for size in range(1, 5):
    signature_orbits: dict[tuple[int, ...], set[int]] = defaultdict(set)
    for mask in range(1 << P):
        if mask.bit_count() != size:
            continue
        orbit_key = min(dihedral_orbit(mask))
        signature_orbits[autocorrelation(mask)].add(orbit_key)
    ambiguous = {
        signature: orbits
        for signature, orbits in signature_orbits.items()
        if len(orbits) > 1
    }
    if ambiguous and first_homometric_support is None:
        first_homometric_support = size
    if size == 4:
        homometric_signature_count_at_four = len(ambiguous)
        require(
            autocorrelation(A) in ambiguous,
            "displayed hostile missing from size-four ambiguity bank",
        )
        require(
            len(ambiguous[autocorrelation(A)]) == 2,
            "displayed size-four signature no longer has exactly two orbits",
        )

require(first_homometric_support == 4, "first homometric support changed")
require(homometric_signature_count_at_four == 1, "size-four ambiguity count changed")


print("THM-2439 exact companion")
print(f"restricted_mask_checks={restricted_checks}")
print(f"restricted_marker_degree={restricted_degree}")
print(f"restricted_nonzero_coefficients={restricted_nonzero_coefficients}")
print(
    "degree_ten_witnesses="
    f"{coefficients[degree_ten_support]},{coefficients[degree_ten_zero_support]}"
)
print(f"universal_mask_checks={universal_checks}")
print(f"universal_marker_degree={universal_degree}")
print(f"universal_nonzero_coefficients={universal_nonzero_coefficients}")
print(
    "nonzero_coefficients_by_degree="
    + ",".join(map(str, nonzero_coefficients_by_degree))
)
print("marker_truths_by_size=" + ",".join(map(str, marker_truths_by_size)))
print(
    "degree_twelve_witnesses="
    f"{coefficients[FULL ^ (1 << 6)]},{coefficients[FULL ^ (1 << 10)]}"
)
print(f"marker_equivariance_checks={equivariance_checks}")
print(f"homometric_A_word={A_marked_word}")
print(f"homometric_B_word={B_marked_word}")
print(
    "affine_equivalences="
    + ";".join(f"{unit},{shift}" for unit, shift in affine_equivalences)
)
print("homometric_autocorrelation=4,1x12")
print("translation_gram=0,6x12")
print("normalized_fourier_powers=16/169,3/169x12")
print(f"first_nondihedral_homometric_support={first_homometric_support}")
print(f"size_four_ambiguous_signatures={homometric_signature_count_at_four}")
print("ALL CHECKS PASSED")
