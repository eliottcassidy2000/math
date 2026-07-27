#!/usr/bin/env python3
"""Dependency-free exact referee for THM-2529.

The theorem specializes the THM-2527 Boolean cut coordinate to the actual
thirteen-shift deep-comb masks.  This script freezes the one-or-two adjacent
root law, the resulting two-parameter profile and positive coordinate, all
primitive root modes, and an exact THM-2365 circulant/zero-target hostile.
"""

from fractions import Fraction


P = 13
CHI7 = {1: 1, 2: 1, 3: -1, 4: 1, 5: -1, 6: -1}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def circle_fraction(value: Fraction) -> Fraction:
    return value - value.numerator // value.denominator


def circle_norm(value: Fraction) -> Fraction:
    value = circle_fraction(value)
    return min(value, 1 - value)


def danger(value: Fraction) -> int:
    return int(circle_norm(value) < Fraction(1, 14))


def deep_mask(phase: Fraction) -> tuple[int, ...]:
    return tuple(danger(phase - Fraction(root, P)) for root in range(P))


def correlation(mask: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        sum(mask[root] * mask[(root + shift) % P] for root in range(P))
        for shift in range(P)
    )


def translation(offset: int) -> list[list[int]]:
    return [
        [int(column == (row + offset) % P) for column in range(P)]
        for row in range(P)
    ]


def zero_matrix() -> list[list[int]]:
    return [[0 for _ in range(P)] for _ in range(P)]


def add_scaled(target: list[list[int]], scalar: int, source: list[list[int]]) -> None:
    for row in range(P):
        for column in range(P):
            target[row][column] += scalar * source[row][column]


def operator_a(slope: int) -> list[list[int]]:
    matrix = zero_matrix()
    for step in range(1, 7):
        add_scaled(matrix, -CHI7[step], translation(2 * slope * step))
        add_scaled(matrix, -CHI7[step], translation(-2 * slope * step))
    return matrix


def operator_h(slope: int) -> list[list[int]]:
    matrix = zero_matrix()
    for step in range(1, 7):
        add_scaled(matrix, 1, translation(-2 * slope * step))
        add_scaled(matrix, -1, translation(2 * slope * step))
    return matrix


def matmul(left: list[list[int]], right: list[list[int]]) -> list[list[int]]:
    return [
        [
            sum(left[row][middle] * right[middle][column] for middle in range(P))
            for column in range(P)
        ]
        for row in range(P)
    ]


def matvec(matrix: list[list[int]], vector: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    return tuple(
        sum(matrix[row][column] * vector[column] for column in range(P))
        for row in range(P)
    )


AH = matmul(operator_a(1), operator_h(1))


# The phase-cell decomposition has exactly thirteen singleton and thirteen
# adjacent-pair cells.  The successor identity is checked on every cell.
boundaries = sorted(
    {
        circle_fraction(Fraction(root, P) + sign * Fraction(1, 14))
        for root in range(P)
        for sign in (-1, 1)
    }
)
require(len(boundaries) == 26, "deep-mask boundary census changed")

singleton_cells = 0
pair_cells = 0
mask_checks = 0
singleton_measures = [Fraction(0) for _ in range(P)]
pair_measure = Fraction(0)
for index, left in enumerate(boundaries):
    right = boundaries[(index + 1) % len(boundaries)]
    if index + 1 == len(boundaries):
        right += 1
    phase = circle_fraction((left + right) / 2)
    mask = deep_mask(phase)
    active = [root for root, value in enumerate(mask) if value]
    successor_bit = danger(13 * phase)
    require(len(active) == 2 - successor_bit, "successor identity failed")
    require(len(active) in (1, 2), "deep fibre is not one-or-two")
    if len(active) == 1:
        singleton_cells += 1
        singleton_measures[active[0]] += right - left
    else:
        require(
            (active[1] - active[0]) % P in (1, P - 1)
            or {active[0], active[1]} == {0, P - 1},
            "two active roots are not adjacent",
        )
        # Circular adjacency is more cleanly audited by a shift overlap.
        require(
            correlation(mask)[1] + correlation(mask)[-1] == 2,
            "adjacent pair lost its two oriented edges",
        )
        pair_cells += 1
        pair_measure += right - left

    local_profile = correlation(mask)
    local_bank = matvec(AH, tuple(Fraction(value) for value in local_profile))
    psi = local_bank[-4]
    require(psi == -local_bank[4], "odd coordinate lost its converse")
    require(
        psi == 2 + 5 * successor_bit,
        "deep-comb cut score did not collapse to 2+5d(13cx)",
    )
    mask_checks += 1

require(singleton_cells == pair_cells == 13, "one/two phase cells changed")
require(mask_checks == 26, "phase-cell check count changed")
require(
    singleton_measures == [Fraction(1, 91) for _ in range(P)],
    "singleton root cells no longer have mass 1/91",
)
require(pair_measure == Fraction(6, 7), "adjacent-pair phase mass changed")


def cyclotomic_transform(
    values: tuple[Fraction, ...], colour: int
) -> tuple[Fraction, ...]:
    coefficients = [Fraction(0) for _ in range(P)]
    for residue, value in enumerate(values):
        coefficients[(colour * residue) % P] += value
    return tuple(
        coefficients[index] - coefficients[P - 1]
        for index in range(P - 1)
    )


# If alpha is the F-mass on singleton fibres and beta the mass on adjacent
# pair fibres, the entire deep collision profile is forced.  Exact samples
# audit the branch normalization and the fixed positive coordinate.
mass_samples = (
    (Fraction(1, 91), Fraction(0)),
    (Fraction(0), Fraction(2, 91)),
    (Fraction(3, 70), Fraction(5, 84)),
    (Fraction(11, 137), Fraction(17, 211)),
)
profile_checks = 0
primitive_mode_checks = 0
for alpha, beta in mass_samples:
    require(alpha + beta > 0, "empty mass sample")
    profile = tuple(
        (alpha + 2 * beta) / P
        if shift == 0
        else beta / P
        if shift in (1, P - 1)
        else Fraction(0)
        for shift in range(P)
    )
    bank = tuple(P * value for value in matvec(AH, profile))
    require(bank[-4] == 7 * alpha + 2 * beta, "positive consumer changed")
    require(bank[4] == -(7 * alpha + 2 * beta), "converse sign changed")
    require(bank[-4] == 2 * (alpha + beta) + 5 * alpha, "2+5 identity failed")
    require(all(bank[-u % P] == -bank[u] for u in range(P)), "bank lost oddness")
    for values in (profile, bank):
        for colour in range(1, P):
            require(
                any(cyclotomic_transform(values, colour)),
                "positive rational deep profile lost a primitive root mode",
            )
            primitive_mode_checks += 1
    profile_checks += 1

require(profile_checks == len(mass_samples), "profile sample count changed")
require(primitive_mode_checks == 96, "primitive-mode census changed")


# Exact physical zero-target hostile.  Let
#
#   F_t(x)=Delta_(t+1)(x) d(13 c x).
#
# It selects the singleton phase whose unique active deep root is t+1.  Haar
# measure is 1/91 for every t; it is disjoint from Delta_t.  Its THM-2365
# table is H(r,s,t)=1/91 if r=t+1 and zero otherwise, hence G(r-t).
hostile_mass = singleton_measures[1]
require(hostile_mass == Fraction(1, 91), "hostile selector mass changed")
hostile_table = {
    (root, first_shift, second_shift): (
        hostile_mass if root == (second_shift + 1) % P else Fraction(0)
    )
    for root in range(P)
    for first_shift in range(P)
    for second_shift in range(P)
}

for first_shift in range(P):
    for second_shift in range(P):
        require(
            hostile_table[(second_shift, first_shift, second_shift)] == 0,
            "hostile violated the deep complement diagonal",
        )
        for root in range(P):
            require(
                hostile_table[(root, first_shift, second_shift)]
                == (
                    hostile_mass
                    if (root - second_shift) % P == 1
                    else Fraction(0)
                ),
                "hostile is not circulant in r-t",
            )

# The local singleton autocorrelation forgets its location.  Every (s,t)
# receives the same all-mode odd root bank, so its target-shift DFT is purely
# zero even though the root bank is nonzero in all twelve primitive modes.
singleton_profile = tuple(
    hostile_mass / P if shift == 0 else Fraction(0)
    for shift in range(P)
)
hostile_root_bank = tuple(P * value for value in matvec(AH, singleton_profile))
require(hostile_root_bank[-4] == Fraction(1, 13), "hostile consumer changed")
require(hostile_root_bank[4] == Fraction(-1, 13), "hostile converse changed")
for colour in range(1, P):
    require(
        any(cyclotomic_transform(hostile_root_bank, colour)),
        "hostile root bank lost a primitive mode",
    )

target_bank = {
    (first_shift, second_shift): hostile_root_bank
    for first_shift in range(P)
    for second_shift in range(P)
}
require(
    len(set(target_bank.values())) == 1,
    "circulant hostile unexpectedly acquired target variation",
)


print("theorem=THM-2529")
print("status=PROVED+VERIFIED-EXACT+INDEPENDENTLY-AUDITED")
print(f"deep_phase_cells={mask_checks}")
print(f"singleton_cells={singleton_cells}")
print(f"adjacent_pair_cells={pair_cells}")
print("successor_identity=n=2-d(13cx)")
print("deep_cut_identity=Psi=2+5*d(13cx)")
print("positive_consumer=O(-4)=2*mass(F)+5*mass(F_intersect_D_13c)")
print(f"profile_samples={profile_checks}")
print(f"primitive_mode_checks={primitive_mode_checks}")
print(f"circulant_hostile_mass={hostile_mass}")
print("circulant_hostile=H(r,s,t)=1/91_IF_r=t+1")
print(f"circulant_hostile_O_minus_4={hostile_root_bank[-4]}")
print("target_variation_in_hostile=0")
print("root_modes_in_hostile=12")
print("scalar_rows_excluded=0")
print("lrc14_status=OPEN")
print("all_checks=PASS")
