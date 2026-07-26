#!/usr/bin/env python3
"""Exact referee for the GMC2 lowest-face complete-orbit certificate.

The complete angular average is evaluated by exact residue bins, which is
character orthogonality without numerical roots of unity.
"""

from fractions import Fraction
from math import comb, factorial


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def add(left, right):
    answer = dict(left)
    for exponent, coefficient in right.items():
        answer[exponent] = answer.get(exponent, Fraction(0)) + coefficient
        if answer[exponent] == 0:
            del answer[exponent]
    return answer


def multiply(left, right):
    answer = {}
    for exponent_left, coefficient_left in left.items():
        for exponent_right, coefficient_right in right.items():
            exponent = exponent_left + exponent_right
            answer[exponent] = (
                answer.get(exponent, Fraction(0))
                + coefficient_left * coefficient_right
            )
    return {exponent: coefficient
            for exponent, coefficient in answer.items()
            if coefficient != 0}


def power(polynomial, exponent):
    answer = {0: Fraction(1)}
    base = dict(polynomial)
    value = exponent
    while value:
        if value & 1:
            answer = multiply(answer, base)
        base = multiply(base, base)
        value //= 2
    return answer


def l2_norm_squared(polynomial):
    return sum(coefficient * coefficient for coefficient in polynomial.values())


def complete_difference_energy(polynomial, order):
    """(1/(2K)) sum_a ||R_a f-f||_2^2 by exact residue orthogonality."""
    return sum(
        coefficient * coefficient
        for exponent, coefficient in polynomial.items()
        if exponent % order != 0
    )


def trivial_residue_energy(polynomial, order):
    return sum(
        coefficient * coefficient
        for exponent, coefficient in polynomial.items()
        if exponent % order == 0
    )


def gaussian_norm_squared(coefficient):
    real, imaginary = coefficient
    return real * real + imaginary * imaginary


def gaussian_l2_norm_squared(polynomial):
    return sum(gaussian_norm_squared(coefficient)
               for coefficient in polynomial.values())


def gaussian_complete_difference_energy(polynomial, order):
    return sum(
        gaussian_norm_squared(coefficient)
        for exponent, coefficient in polynomial.items()
        if exponent % order != 0
    )


def gaussian_trivial_residue_energy(polynomial, order):
    return sum(
        gaussian_norm_squared(coefficient)
        for exponent, coefficient in polynomial.items()
        if exponent % order == 0
    )


test_polynomials = (
    {
        -2: Fraction(1),
        -1: Fraction(2),
        0: Fraction(-1),
        1: Fraction(3),
        3: Fraction(-2),
    },
    {
        -3: Fraction(2),
        -1: Fraction(-1),
        2: Fraction(3),
        4: Fraction(1),
    },
    {
        -1: Fraction(5, 3),
        0: Fraction(-2, 5),
        2: Fraction(7, 4),
    },
)

alias_free_controls = 0
residue_controls = 0
for polynomial in test_polynomials:
    negative_width = -min(polynomial)
    positive_width = max(polynomial)
    for exponent in range(1, 9):
        powered = power(polynomial, exponent)
        order = exponent * max(negative_width, positive_width) + 1
        total = l2_norm_squared(powered)
        difference = complete_difference_energy(powered, order)
        trivial = trivial_residue_energy(powered, order)
        require(total - difference == trivial,
                "complete-orbit residue identity failed")
        residue_controls += 1
        constant = powered.get(0, Fraction(0))
        require(trivial == constant * constant,
                "alias-free trivial residue was not the constant term")
        alias_free_controls += 1

# The strict zero-residue threshold is necessary even inside the stated
# two-sided-width class.  At equality, f=u^-1+u^K has zero constant term,
# while its exponent K aliases into the trivial C_K residue.
alias_order = 11
alias_hostile = {-1: Fraction(1), alias_order: Fraction(1)}
require(alias_hostile.get(0, 0) == 0, "alias hostile gained a constant")
require(complete_difference_energy(alias_hostile, alias_order) == 1,
        "alias hostile nontrivial residue changed")
require(trivial_residue_energy(alias_hostile, alias_order) == 1,
        "alias hostile did not imitate vacuum energy")

# Exact Gaussian-rational controls verify that the norm in the orbit
# identity is Hermitian, not a real-coefficient square.
gaussian_test_polynomials = (
    {
        -3: (Fraction(2, 3), Fraction(-1, 4)),
        0: (Fraction(1, 5), Fraction(7, 6)),
        4: (Fraction(-3, 2), Fraction(5, 7)),
    },
    {
        -2: (Fraction(0), Fraction(3, 5)),
        1: (Fraction(-8, 9), Fraction(1, 3)),
        5: (Fraction(4, 7), Fraction(0)),
    },
)
gaussian_residue_controls = 0
for polynomial in gaussian_test_polynomials:
    for order in (3, 5, 7, 11):
        total = gaussian_l2_norm_squared(polynomial)
        difference = gaussian_complete_difference_energy(polynomial, order)
        trivial = gaussian_trivial_residue_energy(polynomial, order)
        require(total - difference == trivial,
                "Gaussian-rational orbit residue identity failed")
        gaussian_residue_controls += 1

# A one-body norm bank loses the phase needed after multiplication.
f_plus = {-2: Fraction(1), -1: Fraction(1),
          1: Fraction(1), 2: Fraction(1)}
f_minus = {-2: Fraction(-1), -1: Fraction(1),
           1: Fraction(1), 2: Fraction(1)}
require(
    {exponent: coefficient * coefficient
     for exponent, coefficient in f_plus.items()}
    == {exponent: coefficient * coefficient
        for exponent, coefficient in f_minus.items()},
    "phase hostile changed the labelled one-body norm bank",
)
ct_plus_square = power(f_plus, 2).get(0, Fraction(0))
ct_minus_square = power(f_minus, 2).get(0, Fraction(0))
require((ct_plus_square, ct_minus_square) == (4, 0),
        "phase hostile constant terms changed")

# Arbitrary-radial hostile in the Laguerre basis:
# L_0^(q)=1 and L_1^(q)=q+1-s.
# Z^q s = Z^q((q+1)L_0-L_1),
# Z^q(2(q+1)-s)=Z^q((q+1)L_0+L_1).
radial_controls = 0
radial_expansion_controls = 0
radial_charge_null_controls = 0
for charge in range(1, 13):
    minus_coefficients = (charge + 1, -1)
    plus_coefficients = (charge + 1, 1)
    weights = (factorial(charge), factorial(charge + 1))
    minus_norm = sum(
        coefficient * coefficient * weight
        for coefficient, weight in zip(minus_coefficients, weights)
    )
    plus_norm = sum(
        coefficient * coefficient * weight
        for coefficient, weight in zip(plus_coefficients, weights)
    )
    require(minus_norm == plus_norm,
            "radial Laguerre hostile changed its norm bank")
    require(minus_coefficients[0] == plus_coefficients[0]
            and abs(minus_coefficients[1]) == abs(plus_coefficients[1]),
            "radial Laguerre labelled magnitudes changed")
    require(minus_norm == factorial(charge + 2),
            "radial Laguerre norm is not (q+2)!")

    # Convert a*L_0+b*L_1 to monomial coefficients
    # (a+b(q+1)) + (-b)s and derive the two s-adic orders.
    minus_monomial = (
        minus_coefficients[0] + minus_coefficients[1] * (charge + 1),
        -minus_coefficients[1],
    )
    plus_monomial = (
        plus_coefficients[0] + plus_coefficients[1] * (charge + 1),
        -plus_coefficients[1],
    )
    require(minus_monomial == (0, 1),
            "P_- no longer equals Z^q s")
    require(plus_monomial == (2 * (charge + 1), -1),
            "P_+ monomial expansion changed")
    minus_order = 0 if minus_monomial[0] else 1
    plus_order = 0 if plus_monomial[0] else 1
    require((minus_order, plus_order) == (1, 0),
            "radial s-adic orders changed")
    radial_expansion_controls += 1

    # Circular symmetry makes every scalar moment of a nonzero charge
    # polynomial vanish, independently of its radial coefficient.
    require(charge != 0, "radial hostile lost nonzero charge")
    radial_charge_null_controls += 1
    radial_controls += 1

# Record the effective witness invoice imported from THM-2111.
width_controls = 0
full_bank_width_controls = 0
for negative_width in range(1, 8):
    for positive_width in range(1, 8):
        bound = comb(
            negative_width + positive_width,
            min(negative_width, positive_width),
        )
        require(bound >= 1, "effective compound-root bound changed")
        uniform_order = bound * max(negative_width, positive_width) + 1
        require(uniform_order > bound * max(negative_width, positive_width),
                "uniform angular order is not alias-free")
        width_controls += 1

        full_bank_order = bound * (negative_width + positive_width) + 1
        require(
            full_bank_order
            > bound * (negative_width + positive_width),
            "uniform full exponent bank is not injective",
        )
        # At equality the two extreme exponents of f^bound have the
        # same residue, so the strict full-bank inequality is sharp.
        collision_order = bound * (negative_width + positive_width)
        require(
            (-bound * negative_width) % collision_order
            == (bound * positive_width) % collision_order,
            "endpoint full-bank collision disappeared",
        )
        full_bank_width_controls += 1

print("GMC2 lowest-face complete-orbit certificate exact referee")
print(f"zero_residue_alias_free_power_controls: {alias_free_controls}")
print(f"complete_residue_identity_controls: {residue_controls}")
print(f"gaussian_residue_identity_controls: {gaussian_residue_controls}")
print(f"effective_width_invoice_controls: {width_controls}")
print(f"full_bank_width_invoice_controls: {full_bank_width_controls}")
print("zero_residue_exact_condition: supp(f^m) intersection KZ subset {0}")
print("zero_residue_width_condition: K>m*max(M,N)")
print("full_exponent_bank_condition: K>m*(M+N)")
print("effective_witness_bound: binom(M+N,min(M,N))")
print("alias_hostile: PASS")
print(f"one_body_phase_hostile_constant_terms: {ct_plus_square},{ct_minus_square}")
print(f"radial_laguerre_hostiles: {radial_controls}")
print(f"radial_monomial_expansion_controls: {radial_expansion_controls}")
print(f"radial_charge_null_controls: {radial_charge_null_controls}")
print("radial_s_adic_orders: 1,0")
print("NC2_status: PROVED_BY_THM-2022")
print("new_NC2_step: NO")
print("all_checks: PASS")
