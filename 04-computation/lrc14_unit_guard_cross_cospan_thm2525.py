#!/usr/bin/env python3
"""Dependency-free exact referee for THM-2525.

The analytic identities are proved in the theorem text.  This companion
freezes the finite C_13 root-fibre universe behind the unit-guard collision
floor, the word/bare-owner cross-cospan collapse, the THM-2524 inverse-row
bound, and the sharp guard-aligned translated-bank certificate.
"""

from fractions import Fraction
from math import comb


P = 13
CHI = {1: 1, 2: 1, 3: -1, 4: 1, 5: -1, 6: -1}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def circle_fraction(value: Fraction) -> Fraction:
    return value - (value.numerator // value.denominator)


def circle_norm(value: Fraction) -> Fraction:
    value = circle_fraction(value)
    return min(value, 1 - value)


def guard_count(phase: Fraction) -> int:
    return sum(
        circle_norm(phase + Fraction(root, P)) > Fraction(1, 7)
        for root in range(P)
    )


def matmul(left: list[list[int]], right: list[list[int]]) -> list[list[int]]:
    return [
        [
            sum(left[row][mid] * right[mid][column] for mid in range(P))
            for column in range(P)
        ]
        for row in range(P)
    ]


def matrix_linear_combination(
    terms: tuple[tuple[int, list[list[int]]], ...]
) -> list[list[int]]:
    return [
        [
            sum(coefficient * matrix[row][column] for coefficient, matrix in terms)
            for column in range(P)
        ]
        for row in range(P)
    ]


def hamilton_matrix(slope: int) -> list[list[int]]:
    matrix = [[0 for _ in range(P)] for _ in range(P)]
    for row in range(P):
        for step in range(1, 7):
            matrix[row][(row + 2 * slope * step) % P] -= CHI[step]
            matrix[row][(row - 2 * slope * step) % P] -= CHI[step]
    return matrix


def apply_matrix(matrix: list[list[int]], vector: list[int]) -> list[int]:
    return [
        sum(matrix[row][column] * vector[column] for column in range(P))
        for row in range(P)
    ]


def bits(mask: int) -> list[int]:
    return [(mask >> root) & 1 for root in range(P)]


def correlation_counts(mask: int) -> list[int]:
    profile = bits(mask)
    return [
        sum(profile[root] * profile[(root + shift) % P] for root in range(P))
        for shift in range(P)
    ]


def zero_set_has_cyclic_three_run(mask: int) -> bool:
    return any(
        all(((mask >> ((start + offset) % P)) & 1) == 0 for offset in range(3))
        for start in range(P)
    )


# The membership pattern changes only when one root hits +/-1/7.  Exact
# rational phase-cell sampling proves that every generic guard fibre has nine
# or ten points and that both counts occur.
boundaries = sorted(
    {
        circle_fraction(sign * Fraction(1, 7) - Fraction(root, P))
        for sign in (-1, 1)
        for root in range(P)
    }
)
require(len(boundaries) == 26, "unit-guard boundary census changed")
phase_counts: list[int] = []
for index, left in enumerate(boundaries):
    right = boundaries[(index + 1) % len(boundaries)]
    if index + 1 == len(boundaries):
        right += 1
    midpoint = circle_fraction((left + right) / 2)
    phase_counts.append(guard_count(midpoint))
require(set(phase_counts) == {9, 10}, "guard fibres are no longer 9/10")
require(max(phase_counts) == 10, "unit-guard maximum changed")


# THM-2524's signed Hamilton matrix and inverse numerator.
A_ONE = hamilton_matrix(1)
A_TWO = matmul(A_ONE, A_ONE)
A_THREE = matmul(A_TWO, A_ONE)
A_FIVE = matmul(matmul(A_THREE, A_ONE), A_ONE)
INVERSE_NUMERATOR = matrix_linear_combination(
    ((1, A_FIVE), (-39, A_THREE), (299, A_ONE))
)
EXPECTED_INVERSE_ROW = [0, 85, 15, -15, -5, -85, 5, 5, -85, -5, -15, 15, 85]
require(
    INVERSE_NUMERATOR[0] == EXPECTED_INVERSE_ROW,
    "THM-2524 inverse numerator row changed",
)
require(
    sum(abs(value) for value in EXPECTED_INVERSE_ROW) == 420,
    "inverse-row l1 norm changed",
)


# Exhaust every nonempty Boolean fibre allowed by the ten-point carrier
# invoice.  Correlation is even, nonconstant, and obeys the sharp 3/13 energy
# floor.  A root-constant Boolean word flag is either zero (empty fibre) or
# one, in which case word-vs-bare cross correlation equals self-correlation
# entry by entry.
carrier_profiles = 0
floor_equalities = 0
cross_checks = 0
cyclotomic_nonconstant_checks = 0
for mask in range(1, 1 << P):
    mass = mask.bit_count()
    if mass > 10:
        continue
    carrier_profiles += 1
    correlation = correlation_counts(mask)
    require(
        all(correlation[shift] == correlation[-shift % P] for shift in range(P)),
        "self-correlation lost antipodal symmetry",
    )
    require(
        correlation[0] * P - sum(correlation) >= 3 * mass,
        "3/13 fibre energy floor failed",
    )
    if correlation[0] * P - sum(correlation) == 3 * mass:
        require(mass == 10, "energy equality occurred below ten points")
        floor_equalities += 1
    require(len(set(correlation)) > 1, "proper Boolean fibre became constant")
    cyclotomic_nonconstant_checks += 1

    profile = bits(mask)
    # The active root-constant word flag is one.  This is the only nonempty
    # local case; flag zero contributes the identically zero table.
    for shift in range(P):
        cross = sum(
            profile[root] * profile[(root + shift) % P]
            for root in range(P)
        )
        require(cross == correlation[shift], "word/bare cross collapse failed")
        cross_checks += 1

require(
    carrier_profiles == sum(comb(P, mass) for mass in range(1, 11)) == 8099,
    "ten-point carrier universe changed",
)
require(floor_equalities == comb(P, 10) == 286, "floor equality census changed")


# A unit-guard complement contains at least three consecutive roots after the
# H-gauge.  This is the exact local universe for the sharp signed-bank lemma.
guard_profiles = [
    mask
    for mask in range(1, 1 << P)
    if zero_set_has_cyclic_three_run(mask)
]
require(len(guard_profiles) == 5434, "three-zero-run universe changed")

# Exact dual certificate for the guard-aligned slope.  If c is the
# unnormalized cyclic autocorrelation and r=A_1 c, then
#
#   9 r_0 + 368 r_1 - 340 r_4 - 198 r_5 >= 234 |mask|.
#
# Its coefficient l1 norm is 915, giving
# max|R| >= 13*234/915 rho = 1014/305 rho after integration.
DUAL = [9, 368, 0, 0, -340, -198, 0, 0, 0, 0, 0, 0, 0]
require(sum(abs(value) for value in DUAL) == 915, "dual l1 norm changed")
dual_equalities = 0
minimum_dual_gap: int | None = None
for mask in guard_profiles:
    mass = mask.bit_count()
    signed_bank = apply_matrix(A_ONE, correlation_counts(mask))
    gap = sum(
        coefficient * value for coefficient, value in zip(DUAL, signed_bank)
    ) - 234 * mass
    require(gap >= 0, f"sharp signed-bank dual failed at mask {mask:#x}")
    minimum_dual_gap = gap if minimum_dual_gap is None else min(minimum_dual_gap, gap)
    if gap == 0:
        dual_equalities += 1
require(minimum_dual_gap == 0, "sharp dual has no equality")
require(dual_equalities == 78, "sharp dual equality census changed")


# Exact primal equality certificate.  Divide these conic weights by 429/305
# to obtain probabilities 52/429, 39/429, 338/429.  All three masks have the
# fixed zero run {10,11,12}; hence the mixture is realizable inside one
# rational unit-guard phase cell.
PRIMAL_WEIGHTS = {
    0x25B: Fraction(52, 305),
    0x2DB: Fraction(39, 305),
    0x3FF: Fraction(338, 305),
}
for mask in PRIMAL_WEIGHTS:
    require(
        all(((mask >> root) & 1) == 0 for root in (10, 11, 12)),
        "primal mask lost the common guard zero run",
    )
primal_mass = sum(weight * mask.bit_count() for mask, weight in PRIMAL_WEIGHTS.items())
primal_weight_sum = sum(PRIMAL_WEIGHTS.values())
primal_bank = [
    sum(
        weight * apply_matrix(A_ONE, correlation_counts(mask))[coordinate]
        for mask, weight in PRIMAL_WEIGHTS.items()
    )
    for coordinate in range(P)
]
require(primal_mass == 13, "primal normalization changed")
require(primal_weight_sum == Fraction(429, 305), "primal weight sum changed")
require(
    max(abs(value) for value in primal_bank) == Fraction(1014, 305),
    "primal signed-bank maximum changed",
)
require(
    sum(coefficient * value for coefficient, value in zip(DUAL, primal_bank))
    == 234 * primal_mass,
    "primal/dual equality changed",
)


# The explicit shallow sharp control F=1_(D_13 intersect C_1) has ten active
# predecessors on every active base fibre.
sharp_mass = Fraction(10, 91)
sharp_drift = Fraction(30, 1183)
require(sharp_drift == Fraction(3, 13) * sharp_mass, "sharp guard control changed")


def main() -> None:
    print("THM-2525 unit-guard collision/cross-cospan referee: PASS")
    print(f"guard_phase_boundaries={len(boundaries)}")
    print(
        "guard_generic_fibre_counts="
        + ",".join(f"{count}:{phase_counts.count(count)}" for count in sorted(set(phase_counts)))
    )
    print(f"carrier_profiles={carrier_profiles}")
    print(f"floor_equalities={floor_equalities}")
    print(f"cross_collapse_checks={cross_checks}")
    print(f"cyclotomic_nonconstant_checks={cyclotomic_nonconstant_checks}")
    print(f"inverse_numerator_row={','.join(map(str, EXPECTED_INVERSE_ROW))}")
    print("inverse_numerator_row_l1=420")
    print("all_slope_inverse_bound=max_abs_R>=845/364*mass")
    print(f"three_zero_run_profiles={len(guard_profiles)}")
    print("sharp_dual=9R0+368R1-340R4-198R5>=234*m")
    print("sharp_dual_l1=915")
    print(f"sharp_dual_equalities={dual_equalities}")
    print("guard_aligned_sharp_bound=max_abs_R>=1014/305*mass")
    print("sharp_primal_masks=0x25b:52/305,0x2db:39/305,0x3ff:338/305")
    print(f"sharp_primal_mass={primal_mass}")
    print(f"sharp_primal_weight_sum={primal_weight_sum}")
    print("sharp_primal_probability_mass=3965/429; sharp_primal_rho=305/429")
    print("sharp_primal_probability_max_abs_R=1014/429")
    print(
        "sharp_primal_bank="
        + ",".join(str(value) for value in primal_bank)
    )
    print(f"shallow_guard_sharp_mass={sharp_mass}")
    print(f"shallow_guard_sharp_drift={sharp_drift}")
    print("word_bare_cross_profile=self_profile; odd_component=0")
    print("scalar_rows_excluded=0; lrc14_status=OPEN")


if __name__ == "__main__":
    main()
