#!/usr/bin/env python3
"""Exact companion for THM-3096's arbitrary-radial pair radical."""

from fractions import Fraction
from itertools import product
from math import comb, factorial

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def compositions(total, parts):
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in compositions(total - first, parts - 1):
            yield (first,) + tail


def moment_forms(lambdas, variables):
    forms = []
    slots = len(lambdas)
    for power in range(1, slots + 1):
        form = 0
        for alpha in compositions(power, slots):
            multinomial = factorial(power)
            radial_degree = 0
            denominator = 1
            monomial = 1
            for index, count in enumerate(alpha):
                multinomial //= factorial(count)
                radial_degree += count * lambdas[index]
                denominator *= factorial(lambdas[index]) ** count
                monomial *= variables[index] ** count
            form += sp.Rational(
                multinomial * factorial(radial_degree), denominator
            ) * monomial
        forms.append(sp.expand(form))
    return forms


def homogeneous_monomials(variables, degree):
    return [
        sp.prod(variable**exponent for variable, exponent in zip(variables, alpha))
        for alpha in compositions(degree, len(variables))
    ]


def coefficient_vector(poly, variables, monomials):
    expanded = sp.Poly(sp.expand(poly), *variables)
    return sp.Matrix([expanded.coeff_monomial(monomial) for monomial in monomials])


def ci_certificate(forms, variables, target, socle_cutoff):
    target_monomials = homogeneous_monomials(variables, socle_cutoff)
    columns = []
    for degree, form in enumerate(forms, start=1):
        for multiplier in homogeneous_monomials(variables, socle_cutoff - degree):
            columns.append(coefficient_vector(multiplier * form, variables, target_monomials))
    matrix = sp.Matrix.hstack(*columns)
    right = coefficient_vector(target, variables, target_monomials)
    solution_set = sp.linsolve((matrix, right))
    require(solution_set is not sp.EmptySet, "complete-intersection certificate")
    solution = next(iter(solution_set))
    free = set().union(*(entry.free_symbols for entry in solution))
    specialization = {symbol: 0 for symbol in free}
    vector = sp.Matrix([entry.subs(specialization) for entry in solution])
    require(matrix * vector == right, "literal certificate identity")
    return len(columns)


def evaluate_form(form, variables, values):
    return sp.Rational(form.subs(dict(zip(variables, values))))


def direct_gaussian_moment(lambdas, x_values, alpha_value, power):
    # alpha W and x_i/lambda_i! * Z (ZW)^(lambda_i-1).
    coefficients = [Fraction(alpha_value)] + [
        Fraction(x_values[index], factorial(lambdas[index]))
        for index in range(len(lambdas))
    ]
    z_degrees = [0] + list(lambdas)
    w_degrees = [1] + [value - 1 for value in lambdas]
    answer = Fraction(0)
    for counts in compositions(power, len(coefficients)):
        total_z = sum(count * degree for count, degree in zip(counts, z_degrees))
        total_w = sum(count * degree for count, degree in zip(counts, w_degrees))
        if total_z != total_w:
            continue
        multinomial = factorial(power)
        coefficient = Fraction(1)
        for count, value in zip(counts, coefficients):
            multinomial //= factorial(count)
            coefficient *= value**count
        answer += multinomial * factorial(total_z) * coefficient
    return answer


# Exact physical-resultant controls after F_1=sum x_i elimination.
support_bank = (
    (1,),
    (1, 3),
    (2, 7),
    (1, 2, 3),
    (1, 2, 4),
    (1, 3, 7),
    (2, 5, 9),
)
resultant_cells = 0
forms_bank = {}
for lambdas in support_bank:
    variables = sp.symbols(f"x0:{len(lambdas)}")
    forms = moment_forms(lambdas, variables)
    forms_bank[lambdas] = (variables, forms)
    require(forms[0] == sum(variables), "normalized first moment")
    if len(lambdas) == 1:
        resultant = 1
    elif len(lambdas) == 2:
        reduced = sp.expand(forms[1].subs(variables[0], -variables[1]))
        resultant = sp.Poly(reduced, variables[1]).coeff_monomial(variables[1] ** 2)
    else:
        reduced2 = sp.expand(forms[1].subs(variables[0], -variables[1] - variables[2]))
        reduced3 = sp.expand(forms[2].subs(variables[0], -variables[1] - variables[2]))
        affine2 = sp.expand(reduced2.subs(variables[2], 1))
        affine3 = sp.expand(reduced3.subs(variables[2], 1))
        require(sp.Poly(affine2, variables[1]).LC() != 0, "no quadratic infinity root")
        require(sp.Poly(affine3, variables[1]).LC() != 0, "no cubic infinity root")
        resultant = sp.resultant(affine2, affine3, variables[1])
    require(resultant != 0, "physical support resultant")
    resultant_cells += 1


# Hilbert series and socle degree for degrees 1,...,t.
hilbert_cells = 0
for slots in range(1, 9):
    socle_degree = slots * (slots - 1) // 2
    limit = socle_degree + 3
    numerator = [1] + [0] * limit
    for degree in range(1, slots + 1):
        updated = numerator[:]
        for index in range(degree, limit + 1):
            updated[index] -= numerator[index - degree]
        numerator = updated
    hilbert = []
    for degree in range(limit + 1):
        value = 0
        for index in range(degree + 1):
            value += numerator[index] * comb(degree - index + slots - 1, slots - 1)
        hilbert.append(value)
    require(all(value > 0 for value in hilbert[: socle_degree + 1]), "positive CI Hilbert range")
    require(all(value == 0 for value in hilbert[socle_degree + 1 :]), "socle cutoff")
    require(sum(hilbert) == factorial(slots), "complete-intersection length")
    hilbert_cells += 1


# Literal m^(D+1) subset (F_1,...,F_t) certificates on every checked support.
certificate_cells = 0
certificate_unknowns = 0
for lambdas in support_bank:
    if len(lambdas) > 3:
        continue
    variables, forms = forms_bank[lambdas]
    cutoff = len(lambdas) * (len(lambdas) - 1) // 2 + 1
    for variable in variables:
        certificate_unknowns += ci_certificate(forms, variables, variable**cutoff, cutoff)
        certificate_cells += 1


# Direct Wick expansion against the balanced-channel formula, including all
# odd moments and the full conjectural 2t cutoff.
moment_cells = 0
radical_grid_cells = 0
for lambdas in support_bank:
    variables, forms = forms_bank[lambdas]
    slots = len(lambdas)
    for x_values in product((-2, -1, 0, 1, 2), repeat=slots):
        if slots == 3 and sum(abs(value) for value in x_values) > 3:
            continue
        form_values = [evaluate_form(form, variables, x_values) for form in forms]
        for alpha_value in (-2, -1, 0, 1, 2):
            all_even_zero = True
            for power in range(1, 2 * slots + 1):
                direct = direct_gaussian_moment(lambdas, x_values, alpha_value, power)
                if power % 2:
                    predicted = 0
                else:
                    half = power // 2
                    predicted = (
                        comb(power, half)
                        * alpha_value**half
                        * form_values[half - 1]
                    )
                    if predicted != 0:
                        all_even_zero = False
                require(direct == predicted, "Wick radial moment identity")
                moment_cells += 1
            require(
                all_even_zero == (alpha_value == 0 or all(value == 0 for value in x_values)),
                "exact finite-grid radical",
            )
            radical_grid_cells += 1


# The neutral channel can cancel the visible pair channel at level two.
s = sp.symbols("s")
b = sp.I * sp.sqrt(2) * (s - 1)


def laplace(poly):
    expanded = sp.Poly(sp.expand(poly), s)
    return sp.expand(
        sum(coefficient * factorial(degree[0]) for degree, coefficient in expanded.terms())
    )


neutral_m1 = laplace(b)
neutral_m2 = laplace(b**2 + 2 * s)
neutral_m3 = laplace(b**3 + 6 * s * b)
require(neutral_m1 == 0 and neutral_m2 == 0, "neutral cancellation hostile")
require(sp.simplify(neutral_m3 - 2 * sp.I * sp.sqrt(2)) == 0, "neutral third moment")


print("THM-3096 PHYSICAL-SUPPORT COMPLETE INTERSECTION")
print(f"resultant_cells={resultant_cells} F1_elimination=PASS")
print(f"hilbert_cells={hilbert_cells} socle_D=t*(t-1)/2 length=t!=PASS")
print(f"certificate_cells={certificate_cells} unknown_columns={certificate_unknowns} literal_radicals=PASS")
print(f"moment_cells={moment_cells} odd_zero_even_factorization=PASS")
print(f"radical_grid_cells={radical_grid_cells} VJ=V(alpha)_union_V(c)=PASS")
print("cutoff=2t=(K-1)*R primitive_return_R=2")
print("certificate=(alpha*x_i)^(D+1)_in_<M2,...,M2t>")
print("neutral_hostile=M1=0;M2=0;M3=2*i*sqrt(2)")
print("scope=arbitrary_complex_radial_coefficients_on_resultant_good_support")
print("all_exact_checks=PASS")
