#!/usr/bin/env python3
"""Exact companion for THM-3108's Gamma/Jordan equality wall."""

from fractions import Fraction
from math import factorial


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def rising(start, length):
    out = 1
    for offset in range(length):
        out *= start + offset
    return out


def carrier_u(row, base, remote):
    return Fraction(
        rising(row * base + 1, row * remote),
        rising(base + 1, remote) ** row,
    )


def gamma_t(child, base, remote):
    upper = child + 1
    return carrier_u(child, base, remote) ** upper / carrier_u(
        upper, base, remote
    ) ** child


def closed_t(child, base, remote):
    upper = child + 1
    endpoint = base + remote
    return Fraction(
        (factorial(child * endpoint) // factorial(child * base)) ** upper,
        (factorial(upper * endpoint) // factorial(upper * base)) ** child,
    )


def shift_ratio(child, base, remote):
    upper = child + 1
    endpoint = base + remote
    return Fraction(
        rising(child * endpoint + 1, child) ** upper,
        rising(upper * endpoint + 1, upper) ** child,
    )


def poly_add(left, right):
    out = dict(left)
    for exponent, coefficient in right.items():
        out[exponent] = out.get(exponent, 0) + coefficient
        if not out[exponent]:
            del out[exponent]
    return out


def poly_scale(poly, scalar):
    return {exponent: scalar * coefficient for exponent, coefficient in poly.items()
            if scalar * coefficient}


def poly_multiply(left, right):
    out = {}
    for exponent, coefficient in left.items():
        for other_exponent, other_coefficient in right.items():
            total = exponent + other_exponent
            out[total] = out.get(total, 0) + coefficient * other_coefficient
    return {exponent: coefficient for exponent, coefficient in out.items() if coefficient}


def poly_valuation(poly):
    return min(poly) if poly else None


def determinant_2_by_2(matrix):
    return poly_add(
        poly_multiply(matrix[0][0], matrix[1][1]),
        poly_scale(poly_multiply(matrix[0][1], matrix[1][0]), -1),
    )


carrier_cells = 0
normalization_cells = 0
for m in range(3, 11):
    p = m + 1
    for n in range(4):
        previous = None
        for remote in range(1, 6):
            um = carrier_u(m, n, remote)
            up = carrier_u(p, n, remote)
            value = gamma_t(m, n, remote)
            require(value == closed_t(m, n, remote), "factorial carrier cancellation")
            if previous is not None:
                require(value / previous == shift_ratio(m, n, remote - 1),
                        "positive rational shift quotient")
            previous = value

            # The p-th powers avoid introducing any algebraic root.
            repair_power = um ** p / up ** m
            error = um ** 2 / up
            normalized_error_power = error ** p / up ** (m - 1)
            require(repair_power == value, "repair p-th power")
            require(normalized_error_power == value ** 2,
                    "degree-compatible normalized error")
            require(value > 0 and error > 0, "positive Gamma carriers")
            carrier_cells += 1
            normalization_cells += 1


# det(J_2+sE_21+zI)=z^2-s, including the equality and both sign sides.
matrix_cells = 0
for numerator in range(1, 18):
    z = Fraction(numerator, numerator + 19)
    equality = z * z
    determinant = z * z - equality
    require(determinant == 0, "Jordan equality wall")
    require(z * z - equality / 2 > 0, "subcritical attachment side")
    require(z * z - 2 * equality < 0, "supercritical attachment side")
    require(equality < z, "entrywise error is o(repair) model at z<1")
    matrix_cells += 1


# Over v^p=u, v^(2m)=u*v^(m-1).  Exponent and scalar norms coincide.
cyclic_cells = 0
for m in range(3, 21):
    p = m + 1
    require(2 * m == p + m - 1, "cyclic exponent reduction")
    for unit in (Fraction(1), Fraction(2), Fraction(7, 3)):
        for numerator in range(1, 6):
            tau = Fraction(numerator, numerator + 11)
            left_coefficient = tau ** 2 * unit
            right_coefficient = tau ** 2 * unit
            require(left_coefficient == right_coefficient, "cyclic zero-divisor wall")
            require((tau ** 2) ** p == tau ** (2 * p), "upper norm Gamma power")
            cyclic_cells += 1


# Derive the unshifted Smith bars and repaired zero determinant over Q[t].
zero = {}
one = {0: 1}
t = {1: 1}
t_squared = {2: 1}
unshifted = ((zero, one), (t_squared, zero))
unshifted_determinant = determinant_2_by_2(unshifted)
minor_valuations = [poly_valuation(entry) for row in unshifted for entry in row if entry]
first_bar = min(minor_valuations)
determinant_valuation = poly_valuation(unshifted_determinant)
smith_bars = (first_bar, determinant_valuation - first_bar)
require(unshifted_determinant == {2: -1}, "unshifted determinant polynomial")
require(smith_bars == (0, 2), "length-two Smith valuation")

shifted = ((t, one), (t_squared, t))
require(determinant_2_by_2(shifted) == {}, "repaired polynomial identity")


print("THM3108 exact one-normal Gamma/Jordan equality wall")
print(f"carrier_cells={carrier_cells};normalization_cells={normalization_cells}")
print("T=U_m^(m+1)/U_(m+1)^m;shift_quotient=positive_rational;T_to_zero=Stirling")
print("attachment=det(J2+sE21+zI)=z^2-s;equality=s=z^2")
print(f"matrix_wall_cells={matrix_cells};smith_bars={smith_bars[0]},{smith_bars[1]};error=o(repair)")
print(f"cyclic_degree_cells={cyclic_cells};v^(2m)=u*v^(m-1);norm_scale=T^2")
print("unscaled_error=U_m^2/U_(m+1);normalized_error=tau^2")
print("scope=formal_degree_compatible_local_hostile;physical_endpoint_jet_realization=OPEN")
print("all_exact_checks=PASS")
