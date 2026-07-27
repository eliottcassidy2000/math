#!/usr/bin/env python3
"""Exact nodal atlas for the quartic B--D ratio bank.

On the B=1, C=W=0 chart of THM-2262/2297, put x=y^2 and then

    z = x - 27*u + 120/7.

The spectral equation becomes a plane cubic H(u,z;D), with D occurring
only in the term 1959552*D*z.  This standard-library certificate proves:

* the singular-point elimination is z^2*S4(z);
* z=0 is precisely the central ratio D=25/126;
* the four nonzero singular points map bijectively to the four roots of
  THM-2311's quartic p_BD,3(D);
* both quartics are irreducible modulo 11 and hence over Q; and
* every quartic-ratio singularity has nondegenerate tangent cone; and
* restoring y through y^2=x gives a six-branch genus-two cover.

Thus every quartic-ratio quotient cubic is a geometrically irreducible
nodal cubic with rational normalization over its quartic ratio field,
while the original y-curve has genus two.  The genus/deck argument then
closes the four ratios in the inherited THM-2262 branch.
"""

from __future__ import annotations

from fractions import Fraction
from functools import reduce
import importlib.util
from math import gcd
from pathlib import Path


HERE = Path(__file__).resolve().parent
FACTOR_SCRIPT = HERE / "jc2_degree18_central_t_factor_exact_probe.py"
SPEC = importlib.util.spec_from_file_location("central_factor", FACTOR_SCRIPT)
FACTOR = importlib.util.module_from_spec(SPEC)
if SPEC.loader is None:
    raise RuntimeError("could not load exact polynomial helpers")
SPEC.loader.exec_module(FACTOR)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def evaluate(poly, value):
    out = Fraction(0)
    for coefficient in reversed(poly):
        out = out * value + coefficient
    return out


def monic(poly):
    poly = FACTOR.trim(poly[:])
    return FACTOR.scale(poly, Fraction(1, 1) / poly[-1])


def polynomial_gcd(first, second):
    first = FACTOR.trim(first[:])
    second = FACTOR.trim(second[:])
    while second != [0]:
        _, remainder = FACTOR.divmod_poly(first, second)
        first, second = second, remainder
    return monic(first)


def primitive_integer_coefficients(poly):
    denominator_lcm = 1
    for coefficient in poly:
        denominator_lcm = (
            denominator_lcm
            * coefficient.denominator
            // gcd(denominator_lcm, coefficient.denominator)
        )
    integers = [int(coefficient * denominator_lcm) for coefficient in poly]
    content = reduce(gcd, (abs(value) for value in integers if value))
    integers = [value // content for value in integers]
    if integers[-1] < 0:
        integers = [-value for value in integers]
    return integers


def reduce_mod(poly, modulus):
    return FACTOR.divmod_poly(poly, modulus)[1]


# ---------------------------------------------------------------------------
# Small multivariate expansion used only to check the coordinate change.


def mv_add(*polys):
    out = {}
    for poly in polys:
        for monomial, coefficient in poly.items():
            out[monomial] = out.get(monomial, Fraction(0)) + coefficient
    return {key: value for key, value in out.items() if value}


def mv_scale(poly, scalar):
    return {
        monomial: scalar * coefficient
        for monomial, coefficient in poly.items()
        if scalar * coefficient
    }


def mv_mul(first, second):
    out = {}
    for left_monomial, left in first.items():
        for right_monomial, right in second.items():
            monomial = tuple(
                left_monomial[index] + right_monomial[index]
                for index in range(3)
            )
            out[monomial] = out.get(monomial, Fraction(0)) + left * right
    return {key: value for key, value in out.items() if value}


def mv_power(poly, exponent):
    out = {(0, 0, 0): Fraction(1)}
    base = dict(poly)
    while exponent:
        if exponent & 1:
            out = mv_mul(out, base)
        base = mv_mul(base, base)
        exponent //= 2
    return out


def verify_coordinate_change() -> None:
    one = {(0, 0, 0): Fraction(1)}
    u = {(1, 0, 0): Fraction(1)}
    z = {(0, 1, 0): Fraction(1)}
    parameter = {(0, 0, 1): Fraction(1)}
    x = mv_add(z, mv_scale(u, 27), mv_scale(one, Fraction(-120, 7)))

    original = mv_add(
        mv_scale(mv_power(u, 3), -26040609),
        mv_scale(mv_power(u, 2), 49601160),
        mv_scale(mv_mul(x, mv_power(u, 2)), 1607445),
        mv_scale(u, -20995200),
        mv_scale(mv_mul(parameter, u), -52907904),
        mv_scale(mv_mul(x, u), -2857680),
        mv_scale(mv_mul(mv_power(x, 2), u), -138915),
        mv_scale(parameter, 33592320),
        mv_scale(x, 777600),
        mv_scale(mv_mul(parameter, x), 1959552),
        mv_scale(mv_power(x, 2), 78120),
        mv_scale(mv_power(x, 3), 1127),
    )

    expected = {
        (0, 0, 0): Fraction(27648000, 7),
        (0, 1, 0): Fraction(-907200),
        (1, 0, 0): Fraction(-37324800),
        (0, 1, 1): Fraction(1959552),
        (0, 2, 0): Fraction(20160),
        (1, 1, 0): Fraction(2993760),
        (2, 0, 0): Fraction(88179840),
        (0, 3, 0): Fraction(1127),
        (1, 2, 0): Fraction(-47628),
        (2, 1, 0): Fraction(-3429216),
        (3, 0, 0): Fraction(-61725888),
    }
    require(original == expected, "the x-to-z coordinate change drifted")


# ---------------------------------------------------------------------------
# Finite-field irreducibility certificate.


def mod_trim(poly, prime):
    poly = [value % prime for value in poly]
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def mod_sub(first, second, prime):
    out = [0] * max(len(first), len(second))
    for index, value in enumerate(first):
        out[index] = (out[index] + value) % prime
    for index, value in enumerate(second):
        out[index] = (out[index] - value) % prime
    return mod_trim(out, prime)


def mod_mul(first, second, prime):
    out = [0] * (len(first) + len(second) - 1)
    for i, left in enumerate(first):
        for j, right in enumerate(second):
            out[i + j] = (out[i + j] + left * right) % prime
    return mod_trim(out, prime)


def mod_divmod(numerator, denominator, prime):
    numerator = mod_trim(numerator[:], prime)
    denominator = mod_trim(denominator[:], prime)
    quotient = [0] * max(1, len(numerator) - len(denominator) + 1)
    inverse = pow(denominator[-1], -1, prime)
    while numerator != [0] and len(numerator) >= len(denominator):
        shift = len(numerator) - len(denominator)
        coefficient = numerator[-1] * inverse % prime
        quotient[shift] = coefficient
        for index, value in enumerate(denominator):
            numerator[index + shift] = (
                numerator[index + shift] - coefficient * value
            ) % prime
        numerator = mod_trim(numerator, prime)
    return mod_trim(quotient, prime), numerator


def mod_reduce(poly, modulus, prime):
    return mod_divmod(poly, modulus, prime)[1]


def mod_power(poly, exponent, modulus, prime):
    out = [1]
    base = mod_reduce(poly, modulus, prime)
    while exponent:
        if exponent & 1:
            out = mod_reduce(mod_mul(out, base, prime), modulus, prime)
        base = mod_reduce(mod_mul(base, base, prime), modulus, prime)
        exponent //= 2
    return out


def mod_gcd(first, second, prime):
    first = mod_trim(first[:], prime)
    second = mod_trim(second[:], prime)
    while second != [0]:
        _, remainder = mod_divmod(first, second, prime)
        first, second = second, remainder
    inverse = pow(first[-1], -1, prime)
    return [(value * inverse) % prime for value in first]


def quartic_irreducible_mod(poly, prime):
    poly = mod_trim(poly[:], prime)
    require(len(poly) == 5, "irreducibility test expects a quartic")
    variable = [0, 1]
    frobenius_two = mod_power(variable, prime**2, poly, prime)
    if len(mod_gcd(poly, mod_sub(frobenius_two, variable, prime), prime)) > 1:
        return False
    frobenius_four = mod_power(variable, prime**4, poly, prime)
    return mod_reduce(
        mod_sub(frobenius_four, variable, prime),
        poly,
        prime,
    ) == [0]


def main() -> None:
    verify_coordinate_change()

    # F=H-z*H_z.  At a singular point H_z=0, so H=0 iff F=0.
    f_u3 = Fraction(-61725888)
    f_u2 = Fraction(88179840)
    f_u1 = [Fraction(-37324800), 0, Fraction(47628)]
    f_u0 = [
        Fraction(27648000, 7),
        0,
        Fraction(-20160),
        Fraction(-2254),
    ]

    h_u2 = Fraction(-185177664)
    h_u1 = [Fraction(176359680), Fraction(-6858432)]
    h_u0 = [
        Fraction(-37324800),
        Fraction(2993760),
        Fraction(-47628),
    ]

    # Interpolate the exact u-resultant as a polynomial in z, then check
    # beyond the interpolation range.
    resultant_values = []
    for z_value in range(13):
        f_at_z = [
            evaluate(f_u0, z_value),
            evaluate(f_u1, z_value),
            f_u2,
            f_u3,
        ]
        hu_at_z = [
            evaluate(h_u0, z_value),
            evaluate(h_u1, z_value),
            h_u2,
        ]
        resultant_values.append(FACTOR.resultant(f_at_z, hu_at_z))
    singular_projection = FACTOR.interpolate_consecutive(resultant_values)
    for z_value in range(13, 21):
        direct = FACTOR.resultant(
            [
                evaluate(f_u0, z_value),
                evaluate(f_u1, z_value),
                f_u2,
                f_u3,
            ],
            [
                evaluate(h_u0, z_value),
                evaluate(h_u1, z_value),
                h_u2,
            ],
        )
        require(
            evaluate(singular_projection, z_value) == direct,
            "singular-projection interpolation failed a hostile value",
        )

    singular_projection_primitive = [
        0,
        0,
        51840000,
        21168000,
        0,
        123480,
        55223,
    ]
    require(
        primitive_integer_coefficients(singular_projection)
        == singular_projection_primitive,
        "singular-point projection factorization changed",
    )
    node_quartic = [Fraction(value) for value in singular_projection_primitive[2:]]
    require(
        polynomial_gcd(node_quartic, FACTOR.derivative(node_quartic))
        == [Fraction(1)],
        "node-coordinate quartic is not squarefree",
    )

    # A linear subresultant recovers the common u-coordinate.  Since the
    # leading u^2 coefficient of H_u is 3*f_u3, first remove u^3.
    k2 = FACTOR.add(
        [f_u2],
        FACTOR.scale(h_u1, Fraction(-1, 3)),
    )
    k1 = FACTOR.add(
        f_u1,
        FACTOR.scale(h_u0, Fraction(-1, 3)),
    )
    linear_coefficient = FACTOR.add(
        FACTOR.scale(k1, h_u2),
        FACTOR.scale(FACTOR.mul(k2, h_u1), -1),
    )
    constant_coefficient = FACTOR.add(
        FACTOR.scale(f_u0, h_u2),
        FACTOR.scale(FACTOR.mul(k2, h_u0), -1),
    )
    u_numerator = FACTOR.scale(constant_coefficient, -1)
    u_denominator = linear_coefficient
    require(
        polynomial_gcd(node_quartic, u_denominator) == [Fraction(1)],
        "linear subresultant vanishes on a node-coordinate root",
    )

    hu_on_node = FACTOR.add(
        FACTOR.add(
            FACTOR.scale(FACTOR.power(u_numerator, 2), h_u2),
            FACTOR.mul(
                FACTOR.mul(h_u1, u_numerator),
                u_denominator,
            ),
        ),
        FACTOR.mul(h_u0, FACTOR.power(u_denominator, 2)),
    )
    f_on_node = FACTOR.add(
        FACTOR.add(
            FACTOR.scale(FACTOR.power(u_numerator, 3), f_u3),
            FACTOR.scale(
                FACTOR.mul(
                    FACTOR.power(u_numerator, 2),
                    u_denominator,
                ),
                f_u2,
            ),
        ),
        FACTOR.add(
            FACTOR.mul(
                FACTOR.mul(f_u1, u_numerator),
                FACTOR.power(u_denominator, 2),
            ),
            FACTOR.mul(f_u0, FACTOR.power(u_denominator, 3)),
        ),
    )
    require(
        reduce_mod(hu_on_node, node_quartic) == [0]
        and reduce_mod(f_on_node, node_quartic) == [0],
        "subresultant node does not lie on F=H_u=0",
    )

    # H_z=0 recovers D as a rational function of z on the node quartic.
    base = [Fraction(907200), Fraction(-40320), Fraction(-3381)]
    linear_in_u = [Fraction(2993760), Fraction(-95256)]
    parameter_numerator = FACTOR.add(
        FACTOR.add(
            FACTOR.mul(base, FACTOR.power(u_denominator, 2)),
            FACTOR.scale(
                FACTOR.mul(
                    FACTOR.mul(linear_in_u, u_numerator),
                    u_denominator,
                ),
                -1,
            ),
        ),
        FACTOR.scale(FACTOR.power(u_numerator, 2), 3429216),
    )
    parameter_denominator = FACTOR.scale(
        FACTOR.power(u_denominator, 2),
        1959552,
    )
    require(
        polynomial_gcd(node_quartic, parameter_denominator)
        == [Fraction(1)],
        "node-to-parameter map has a pole",
    )

    ratio_values = [
        FACTOR.resultant(
            node_quartic,
            FACTOR.add(
                parameter_numerator,
                FACTOR.scale(parameter_denominator, -parameter),
            ),
        )
        for parameter in range(5)
    ]
    ratio_polynomial = FACTOR.interpolate_consecutive(ratio_values)
    for parameter in range(5, 9):
        direct = FACTOR.resultant(
            node_quartic,
            FACTOR.add(
                parameter_numerator,
                FACTOR.scale(parameter_denominator, -parameter),
            ),
        )
        require(
            evaluate(ratio_polynomial, parameter) == direct,
            "node-to-parameter interpolation failed a hostile value",
        )

    bd_quartic = [
        22656250,
        -772734375,
        7600635000,
        -30805790400,
        46376717184,
    ]
    require(
        primitive_integer_coefficients(ratio_polynomial) == bd_quartic,
        "node parameter bank is not p_BD,3",
    )
    bd_quartic_on_node = [Fraction(0)]
    for exponent, coefficient in enumerate(bd_quartic):
        bd_quartic_on_node = FACTOR.add(
            bd_quartic_on_node,
            FACTOR.scale(
                FACTOR.mul(
                    FACTOR.power(parameter_numerator, exponent),
                    FACTOR.power(
                        parameter_denominator,
                        4 - exponent,
                    ),
                ),
                coefficient,
            ),
        )
    require(
        reduce_mod(bd_quartic_on_node, node_quartic) == [0],
        "p_BD,3 does not vanish on the node-field parameter",
    )
    require(
        polynomial_gcd(
            [Fraction(value) for value in bd_quartic],
            FACTOR.derivative([Fraction(value) for value in bd_quartic]),
        )
        == [Fraction(1)],
        "p_BD,3 is not squarefree",
    )

    # The repeated z=0 solution maps to the central ratio.
    require(
        parameter_numerator[0] / parameter_denominator[0]
        == Fraction(25, 126),
        "z=0 no longer maps to the central ratio",
    )
    central_gcd = polynomial_gcd(
        [
            evaluate(f_u0, 0),
            evaluate(f_u1, 0),
            f_u2,
            f_u3,
        ],
        [
            evaluate(h_u0, 0),
            evaluate(h_u1, 0),
            h_u2,
        ],
    )
    require(
        central_gcd == [Fraction(-40, 63), Fraction(1)],
        "central singular point changed",
    )

    # A mod-11 Rabin certificate proves both quartics irreducible over Q.
    prime = 11
    require(
        quartic_irreducible_mod(
            [int(value) for value in node_quartic],
            prime,
        ),
        "node-coordinate quartic is reducible modulo 11",
    )
    require(
        quartic_irreducible_mod(bd_quartic, prime),
        "p_BD,3 is reducible modulo 11",
    )

    # Nonzero Hessian determinant gives two distinct tangent directions.
    huu_numerator = FACTOR.add(
        FACTOR.mul(
            [Fraction(176359680), Fraction(-6858432)],
            u_denominator,
        ),
        FACTOR.scale(u_numerator, -370355328),
    )
    huz_numerator = FACTOR.add(
        FACTOR.mul(
            [Fraction(2993760), Fraction(-95256)],
            u_denominator,
        ),
        FACTOR.scale(u_numerator, -6858432),
    )
    hzz_numerator = FACTOR.add(
        FACTOR.mul(
            [Fraction(40320), Fraction(6762)],
            u_denominator,
        ),
        FACTOR.scale(u_numerator, -95256),
    )
    hessian_numerator = FACTOR.add(
        FACTOR.mul(huu_numerator, hzz_numerator),
        FACTOR.scale(FACTOR.power(huz_numerator, 2), -1),
    )
    hessian_remainder = reduce_mod(hessian_numerator, node_quartic)
    require(
        primitive_integer_coefficients(hessian_remainder)
        == [-59428800, -20055000, -224273, 94325],
        "node Hessian remainder changed",
    )
    require(
        polynomial_gcd(node_quartic, hessian_remainder)
        == [Fraction(1)],
        "a quartic-ratio singularity lost its ordinary-node tangent cone",
    )

    # Restore the forgotten sign y through x=y^2.  On x=0 we have
    # z=-27u+120/7, so the derivative of the restricted cubic is the
    # directional derivative H_u-27H_z, not the changed-coordinate
    # partial H_u.  Its discriminant proves that this directional
    # derivative is nonzero at all three roots for every p_BD,3
    # parameter.
    zero_fibre_a = Fraction(-26040609)
    zero_fibre_b = [Fraction(49601160)]
    zero_fibre_c = [
        Fraction(-20995200),
        Fraction(-52907904),
    ]
    zero_fibre_d = [Fraction(0), Fraction(33592320)]
    zero_fibre_discriminant = FACTOR.add(
        FACTOR.add(
            FACTOR.mul(
                FACTOR.power(zero_fibre_b, 2),
                FACTOR.power(zero_fibre_c, 2),
            ),
            FACTOR.scale(
                FACTOR.power(zero_fibre_c, 3),
                -4 * zero_fibre_a,
            ),
        ),
        FACTOR.add(
            FACTOR.scale(
                FACTOR.mul(
                    FACTOR.power(zero_fibre_b, 3),
                    zero_fibre_d,
                ),
                -4,
            ),
            FACTOR.add(
                FACTOR.scale(
                    FACTOR.power(zero_fibre_d, 2),
                    -27 * zero_fibre_a**2,
                ),
                FACTOR.scale(
                    FACTOR.mul(
                        FACTOR.mul(zero_fibre_b, zero_fibre_c),
                        zero_fibre_d,
                    ),
                    18 * zero_fibre_a,
                ),
            ),
        ),
    )
    zero_fibre_discriminant_primitive = [
        -15625,
        236250,
        -1190700,
        2000376,
    ]
    require(
        primitive_integer_coefficients(zero_fibre_discriminant)
        == zero_fibre_discriminant_primitive,
        "x=0 fibre discriminant changed",
    )
    require(
        polynomial_gcd(
            [Fraction(value) for value in bd_quartic],
            zero_fibre_discriminant,
        )
        == [Fraction(1)],
        "a quartic ratio acquired a repeated x=0 point",
    )

    infinity_a = -26040609
    infinity_b = 1607445
    infinity_c = -138915
    infinity_d = 1127
    infinity_discriminant = (
        infinity_b**2 * infinity_c**2
        - 4 * infinity_a * infinity_c**3
        - 4 * infinity_b**3 * infinity_d
        - 27 * infinity_a**2 * infinity_d**2
        + 18 * infinity_a * infinity_b * infinity_c * infinity_d
    )
    require(
        infinity_discriminant == -153384762202971019112448,
        "infinity points collided",
    )
    branch_places = 3 + 3
    restored_genus = (branch_places - 2) // 2
    require(
        branch_places == 6 and restored_genus == 2,
        "restored y-cover genus count changed",
    )

    print("chart=B=1,C=0,W=0")
    print("coordinate_change=z=x-27u+120/7")
    print(
        "singular_projection_primitive_ascending="
        + ",".join(str(value) for value in singular_projection_primitive)
    )
    print(
        "node_coordinate_quartic_ascending="
        + ",".join(str(int(value)) for value in node_quartic)
    )
    print("node_coordinate_quartic=irreducible_mod_11")
    print("central_node=z=0,u=40/63,D=25/126")
    print(
        "node_parameter_quartic_ascending="
        + ",".join(str(value) for value in bd_quartic)
    )
    print("node_parameter_quartic=THM2311_p_BD3")
    print("node_parameter_quartic=irreducible_mod_11_and_squarefree")
    print("node_to_parameter_map=degree_four_field_isomorphism")
    print(
        "hessian_remainder_primitive_ascending="
        + ",".join(
            str(value)
            for value in [-59428800, -20055000, -224273, 94325]
        )
    )
    print("quartic_ratio_singularity=one_ordinary_node")
    print("projective_cubic=geometrically_irreducible_nodal")
    print("normalization=projective_line_over_quartic_ratio_field")
    print(
        "x_zero_fibre_discriminant_primitive_ascending="
        + ",".join(str(value) for value in zero_fibre_discriminant_primitive)
    )
    print("x_zero_transversality_derivative=H_u-27H_z")
    print("x_divisor_on_normalization=3_simple_zeros_plus_3_simple_poles")
    print("restored_y_curve=double_cover_branched_at_6_places")
    print("restored_y_curve_genus=2")
    print("status=BD_QUARTIC_GENUS_TWO_CLOSURE_CERTIFICATE")


if __name__ == "__main__":
    main()
