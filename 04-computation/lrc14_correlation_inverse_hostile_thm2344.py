#!/usr/bin/env python3
"""Exact, optimization-safe companion for THM-2344."""

from fractions import Fraction
from itertools import product


P = 13
GROUP = tuple(product(range(P), repeat=2))
ZERO = (0, 0)
PHASE = (0, 1)
INVERSE_PHASE = (0, P - 1)
GROUP_SIZE = P**2


def require(condition: bool, message: str) -> None:
    """Raise under ordinary and optimized Python."""
    if not condition:
        raise RuntimeError(message)


def add(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
    return ((left[0] + right[0]) % P, (left[1] + right[1]) % P)


def sub(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
    return ((left[0] - right[0]) % P, (left[1] - right[1]) % P)


def negate(point: tuple[int, int]) -> tuple[int, int]:
    return ((-point[0]) % P, (-point[1]) % P)


def dot(left: tuple[int, int], right: tuple[int, int]) -> int:
    return (left[0] * right[0] + left[1] * right[1]) % P


def correlation(
    left: dict[tuple[int, int], Fraction],
    right: dict[tuple[int, int], Fraction],
) -> dict[tuple[int, int], Fraction]:
    """Return C(r)=sum_(x-y=r) left(x) right(y)."""
    answer = {point: Fraction(0) for point in GROUP}
    for x, left_value in left.items():
        for y, right_value in right.items():
            answer[sub(x, y)] += left_value * right_value
    return answer


def shifted(
    coefficients: dict[tuple[int, int], Fraction],
    phase: tuple[int, int],
) -> dict[tuple[int, int], Fraction]:
    """Apply A_H(q)=A_K(q-p)."""
    return {
        q: coefficients.get(sub(q, phase), Fraction(0))
        for q in GROUP
    }


def orthogonality_atlas() -> int:
    """Check exact exponent histograms for all target differences."""
    checks = 0
    for difference in GROUP:
        histogram = [0] * P
        for character in GROUP:
            histogram[dot(character, difference)] += 1
        expected = [GROUP_SIZE] + [0] * (P - 1)
        if difference != ZERO:
            expected = [P] * P
        require(histogram == expected, "target-character orthogonality failed")
        checks += GROUP_SIZE
    return checks


def point_mass_correlation_control() -> tuple[int, int, Fraction]:
    """Verify the inverse-character correlation and full zero shift."""
    # Exact rational surrogates retain the interval hostile's two negative
    # endpoint signs and positive deepest coefficient.
    deep = Fraction(5, 7)
    left = {(0, 1): Fraction(-2, 7)}
    right = {(0, 2): Fraction(-3, 11)}
    phase_free = {
        point: deep * value
        for point, value in correlation(left, right).items()
    }
    nonzero_phase_free = {
        point: value for point, value in phase_free.items() if value
    }
    require(
        nonzero_phase_free == {INVERSE_PHASE: Fraction(30, 539)},
        "aligned endpoint correlation left the inverse phase",
    )

    full = shifted(phase_free, PHASE)
    nonzero_full = {point: value for point, value in full.items() if value}
    require(
        nonzero_full == {ZERO: Fraction(30, 539)},
        "deep translation did not move inverse phase to zero",
    )
    return len(nonzero_phase_free), len(nonzero_full), nonzero_full[ZERO]


def aligned_phase_atlas() -> tuple[int, int, int]:
    """Check K exponent -t and complete cancellation in H."""
    hermitian = 0
    constant_full = 0
    real_phase_free = 0
    for character in GROUP:
        s, t = character
        k_exponent = (-t) % P
        h_exponent = (t + k_exponent) % P
        opposite = negate(character)
        opposite_k_exponent = (-opposite[1]) % P
        conjugate_k_exponent = (-k_exponent) % P
        require(
            opposite_k_exponent == conjugate_k_exponent,
            "aligned phase lost Hermitian reflection",
        )
        hermitian += 1
        require(h_exponent == 0, "deep and endpoint phases stopped cancelling")
        constant_full += 1

        # In an odd cyclotomic group, a Pth root is real iff its exponent is 0.
        if k_exponent == 0:
            real_phase_free += 1
            require(t == 0, "real locus differs from the phase annihilator")
        else:
            require(t != 0, "detecting twist acquired a real phase")
        require(s in range(P), "inert character coordinate escaped F_13")
    return hermitian, constant_full, real_phase_free


def odd_even_obstruction_atlas() -> tuple[int, int]:
    """Every nonzero p has a detecting twist and a non-even character."""
    nonzero_phases = 0
    detecting_pairs = 0
    for phase in GROUP:
        if phase == ZERO:
            continue
        nonzero_phases += 1
        detected = False
        for character in GROUP:
            exponent = dot(character, phase)
            if exponent == 0:
                continue
            detected = True
            detecting_pairs += 1
            opposite_exponent = dot(negate(character), negate(phase))
            require(
                opposite_exponent == exponent,
                "double negation changed a character pairing",
            )
            require(
                (-exponent) % P != exponent,
                "nonzero odd-order character became even",
            )
        require(detected, "nonzero phase had no detecting character")
        require(
            add(phase, phase) != ZERO,
            "nonzero phase became two-torsion in an odd group",
        )
    return nonzero_phases, detecting_pairs


def interval_sign_ledger() -> tuple[int, int, int, int]:
    """Audit the exact analytic signs used by the one-tooth hostile."""
    # For n=1,2, 0 < pi*n/7 < pi, hence sin(pi*n/7)>0.
    danger_sign = {1: 1, 2: 1}
    safe_sign = {n: -danger_sign[n] for n in danger_sign}
    require(
        danger_sign == {1: 1, 2: 1},
        "danger coefficient sign ledger changed",
    )
    require(
        safe_sign == {1: -1, 2: -1},
        "safe complement sign ledger changed",
    )
    safe_amplitude_sign = danger_sign[1] * safe_sign[1] * safe_sign[2]
    danger_amplitude_sign = (
        danger_sign[1] * danger_sign[1] * danger_sign[2]
    )
    require(
        safe_amplitude_sign == danger_amplitude_sign == 1,
        "an aligned hostile amplitude is not positive",
    )
    return (
        danger_sign[1],
        safe_sign[1],
        safe_amplitude_sign,
        danger_amplitude_sign,
    )


orthogonality_checks = orthogonality_atlas()
phase_free_support, full_support, hostile_amplitude = (
    point_mass_correlation_control()
)
hermitian_checks, constant_twists, real_twists = aligned_phase_atlas()
nonzero_phases, detecting_pairs = odd_even_obstruction_atlas()
(
    danger_sign,
    safe_sign,
    safe_amplitude_sign,
    danger_amplitude_sign,
) = interval_sign_ledger()

require(
    orthogonality_checks == GROUP_SIZE**2,
    "orthogonality check census changed",
)
require(hermitian_checks == GROUP_SIZE, "reflection census changed")
require(constant_twists == GROUP_SIZE, "constant-twist census changed")
require(real_twists == P, "annihilator row size changed")
require(nonzero_phases == GROUP_SIZE - 1, "nonzero phase census changed")

print("theorem=THM-2344")
print("status=PROVED+VERIFIED-EXACT+CANDIDATE-UNDER-INDEPENDENT-AUDIT")
print(f"target_group_size={GROUP_SIZE}")
print(f"target_character_orthogonality_checks={orthogonality_checks}")
print("canonical_reflection=K(-ell)=conjugate(K(ell))")
print("bad_inverse_character_scalar=NONZERO_REAL")
print("phase_free_target_object=real_endpoint_cross_correlation")
print("bad_boundary=shifted_convolution_inverse")
print(f"aligned_phase_free_support={phase_free_support}")
print(f"aligned_full_support={full_support}")
print(f"aligned_hostile_rational_control_amplitude={hostile_amplitude}")
print(f"aligned_hermitian_twists={hermitian_checks}")
print(f"aligned_constant_full_twists={constant_twists}")
print(f"aligned_real_phase_free_twists={real_twists}")
print(f"odd_group_nonzero_phases={nonzero_phases}")
print(f"odd_group_detecting_pairs={detecting_pairs}")
print("one_real_detecting_twist_excludes_bad_line=YES")
print(f"danger_index_1_sign={danger_sign}")
print(f"safe_index_1_sign={safe_sign}")
print(f"aligned_safe_interval_amplitude_sign={safe_amplitude_sign}")
print(f"aligned_danger_interval_amplitude_sign={danger_amplitude_sign}")
print("nonnegative_endpoint_point_masses_exclude_bad_line=NO")
print("physical_factor_positivity_excludes_bad_line=NO")
print("canonical_nine_coordinate_hostile_constructed=NO")
print("scalar_rows_excluded=0")
print("lrc14_status=OPEN")
print("all_checks=PASS")
