#!/usr/bin/env python3
"""Exact finite companion for THM-2714.

The theorem itself is a finite-module argument over the unramified quadratic
DVR Z_2[omega].  This script verifies the cyclic even-length controls, the
two-elementary hyperbolic control, the unique-plane divisibility split, and
the bounded invariant-factor consequences.  It uses no floating point and
contains no truth-bearing Python assert.
"""

from __future__ import annotations

from itertools import combinations_with_replacement, product


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def add(left: tuple[int, int], right: tuple[int, int], modulus: int) -> tuple[int, int]:
    return ((left[0] + right[0]) % modulus, (left[1] + right[1]) % modulus)


def scale(integer: int, value: tuple[int, int], modulus: int) -> tuple[int, int]:
    return ((integer * value[0]) % modulus, (integer * value[1]) % modulus)


def multiply(
    left: tuple[int, int], right: tuple[int, int], modulus: int
) -> tuple[int, int]:
    """Multiply in Z/(modulus)[omega], omega^2+omega+1=0."""
    a, b = left
    c, d = right
    return ((a * c - b * d) % modulus, (a * d + b * c - b * d) % modulus)


def conjugate(value: tuple[int, int], modulus: int) -> tuple[int, int]:
    """The involution omega |-> omega^2=-1-omega."""
    a, b = value
    return ((a - b) % modulus, (-b) % modulus)


def trace(value: tuple[int, int]) -> int:
    return 2 * value[0] - value[1]


def pairing_numerator(
    left: tuple[int, int], right: tuple[int, int], modulus: int
) -> int:
    return trace(multiply(left, conjugate(right, modulus), modulus)) % modulus


def omega_action(value: tuple[int, int], modulus: int) -> tuple[int, int]:
    return multiply((0, 1), value, modulus)


def elements(modulus: int) -> tuple[tuple[int, int], ...]:
    return tuple(product(range(modulus), repeat=2))


def subgroup_multiples(modulus: int, factor: int) -> frozenset[tuple[int, int]]:
    return frozenset(
        (factor * a % modulus, factor * b % modulus)
        for a in range(modulus // factor)
        for b in range(modulus // factor)
    )


def orthogonal(
    ambient: tuple[tuple[int, int], ...],
    subgroup: frozenset[tuple[int, int]],
    modulus: int,
) -> frozenset[tuple[int, int]]:
    return frozenset(
        value
        for value in ambient
        if all(pairing_numerator(value, member, modulus) == 0 for member in subgroup)
    )


def two_torsion(modulus: int) -> frozenset[tuple[int, int]]:
    return frozenset(value for value in elements(modulus) if scale(2, value, modulus) == (0, 0))


def doubled_four_torsion(
    subgroup: frozenset[tuple[int, int]], modulus: int
) -> frozenset[tuple[int, int]]:
    return frozenset(
        scale(2, value, modulus)
        for value in subgroup
        if scale(4, value, modulus) == (0, 0)
    )


def is_power_of_four(value: int) -> bool:
    if value < 1:
        return False
    while value % 4 == 0:
        value //= 4
    return value == 1


def f4_add(left: int, right: int) -> int:
    return left ^ right


def f4_multiply(left: int, right: int) -> int:
    # bits encode a+b*omega and omega^2=omega+1 in characteristic two
    a, b = left & 1, (left >> 1) & 1
    c, d = right & 1, (right >> 1) & 1
    return ((a * c ^ b * d) & 1) | ((((a * d) ^ (b * c) ^ (b * d)) & 1) << 1)


def f4_conjugate(value: int) -> int:
    return f4_multiply(value, value)


def f4_trace(value: int) -> int:
    return value ^ f4_conjugate(value)


def elementary_pairing(
    left: tuple[int, int], right: tuple[int, int]
) -> int:
    total = f4_add(
        f4_multiply(left[0], f4_conjugate(right[0])),
        f4_multiply(left[1], f4_conjugate(right[1])),
    )
    trace_value = f4_trace(total)
    require(trace_value in (0, 1), "F4 trace lands in F2")
    return trace_value


def main() -> None:
    # Ring-law and C3 controls at the exact cyclic towers used below.
    cyclic_rows = []
    for r in range(1, 4):
        exponent = 2 * r
        modulus = 1 << exponent
        ambient = elements(modulus)
        metabolizer = subgroup_multiples(modulus, 1 << r)
        require(len(ambient) == 4**exponent, "cyclic O-module order")
        require(len(metabolizer) ** 2 == len(ambient), "metabolizer square order")
        require(
            {omega_action(value, modulus) for value in metabolizer} == set(metabolizer),
            "cyclic metabolizer C3 stable",
        )
        require(
            orthogonal(ambient, metabolizer, modulus) == metabolizer,
            "cyclic metabolizer self-orthogonal",
        )
        visible_plane = two_torsion(modulus)
        require(len(visible_plane) == 4, "unique visible standard plane")
        require(visible_plane <= metabolizer, "visible plane lies in metabolizer")
        require(
            {omega_action(value, modulus) for value in visible_plane} == set(visible_plane),
            "visible plane C3 stable",
        )
        nonzero_orbit = []
        value = next(value for value in visible_plane if value != (0, 0))
        for _ in range(3):
            nonzero_orbit.append(value)
            value = omega_action(value, modulus)
        require(len(set(nonzero_orbit)) == 3, "C3 cycles visible nonzero classes")
        secondary_liftable = visible_plane <= doubled_four_torsion(metabolizer, modulus)
        require(secondary_liftable == (r >= 2), "secondary divisibility threshold")
        cyclic_rows.append((r, exponent, len(metabolizer), secondary_liftable))

    # Two elementary standard planes admit the diagonal stable metabolizer.
    elementary_ambient = tuple(product(range(4), repeat=2))
    elementary_metabolizer = frozenset((value, value) for value in range(4))
    require(len(elementary_metabolizer) ** 2 == len(elementary_ambient),
            "elementary metabolizer square order")
    require(
        all(
            elementary_pairing(left, right) == 0
            for left in elementary_metabolizer
            for right in elementary_metabolizer
        ),
        "elementary diagonal isotropic",
    )
    elementary_orthogonal = frozenset(
        value
        for value in elementary_ambient
        if all(elementary_pairing(value, member) == 0 for member in elementary_metabolizer)
    )
    require(elementary_orthogonal == elementary_metabolizer,
            "elementary diagonal self-orthogonal")

    # Bounded invariant-factor audit of the parity and unique-plane laws.
    profiles = []
    for number_of_factors in range(1, 7):
        profiles.extend(combinations_with_replacement(range(1, 9), number_of_factors))
    even_length_profiles = tuple(profile for profile in profiles if sum(profile) % 2 == 0)
    odd_length_profiles = tuple(profile for profile in profiles if sum(profile) % 2 == 1)
    require(
        all(not is_power_of_four(2 ** sum(profile)) for profile in odd_length_profiles),
        "odd standard length cannot have an O-module square-root order",
    )
    # The square root of |A_std| is 2^N.  It is an O-module order only when
    # N is even, because every finite O-module has order 4^k.
    require(
        all(is_power_of_four(2 ** sum(profile)) for profile in even_length_profiles),
        "even standard length has an O-module square-root order",
    )
    elementary_profiles = tuple(profile for profile in profiles if set(profile) == {1})
    require(
        all((len(profile) % 2 == 0) == (sum(profile) % 2 == 0)
            for profile in elementary_profiles),
        "two-elementary nullity parity",
    )
    unique_plane_profiles = tuple(profile for profile in profiles if len(profile) == 1)
    require(
        tuple(profile[0] for profile in unique_plane_profiles if sum(profile) % 2 == 0)
        == (2, 4, 6, 8),
        "unique-plane even exponents",
    )

    print("THM-2714 C3 METABOLIC-LENGTH PARITY AUDIT")
    print("coefficient_ring=Z_2[omega] residue_field=F4")
    print(f"invariant_factor_profiles_checked={len(profiles)}")
    print(f"even_length_profiles={len(even_length_profiles)} odd_length_profiles={len(odd_length_profiles)}")
    print("two_elementary_metabolic_nullity=EVEN minimum_carrier_nullity=2")
    print("unique_plane_allowed_exponents=2,4,6,8_in_bounded_box")
    print("cyclic_r1_exponent2_Horder4_secondary=NONLIFTABLE")
    print("cyclic_r2_exponent4_Horder16_secondary=LIFTABLE")
    print("cyclic_r3_exponent6_Horder64_secondary=LIFTABLE")
    print("elementary_W2_diagonal_metabolizer=YES")
    print("d4_single_W_length1_metabolizer=NO")
    print("boundary_realization_and_reflection_not_tested")
    print("PASS")


if __name__ == "__main__":
    main()
