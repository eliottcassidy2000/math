#!/usr/bin/env python3
"""Exact singular atlas for the degree-18 central weight-30 factor.

This standard-library companion works on the B=1 chart of the central wall

    D = 25/126.

It imports the already audited coefficient list and rational-polynomial
helpers from ``jc2_degree18_central_t_factor_exact_probe.py`` and verifies:

* Disc_W(T) = constant * P8(C^2)^2 * P3(C^2)^3;
* finite-field Rabin certificates for P8(C^2) and P3(C^2);
* the complete affine singular locus of T=0;
* the node/cusp critical-section tests in the two quotient fields; and
* the residual degree-ten spectral branch multiplicity profiles.

The final irreducibility and Riemann--Hurwitz deductions are mathematical
consequences of the coefficient controls printed at the end; the script does
not claim to enumerate an infinite complex parameter curve.
"""

from fractions import Fraction
import importlib.util
from pathlib import Path


HERE = Path(__file__).resolve().parent
FACTOR_SCRIPT = HERE / "jc2_degree18_central_t_factor_exact_probe.py"
SPEC = importlib.util.spec_from_file_location("central_t_factor", FACTOR_SCRIPT)
FACTOR = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(FACTOR)


P8 = [
    -1087948019845507812500000000,
    92930586076313704687500000000,
    8954570361875562292640625000000,
    -313103637508400358291807187500000,
    13873912636573086466704217312500000,
    -588322473659000442440395871400000000,
    8764467368105916596083808036836642500,
    -68598223636291030733181729239483693550,
    57091662724422970717844096005443494577,
]

P3 = [
    72250000,
    26254935000,
    1530102783825,
    24450881901768,
]


def rational_monic(poly):
    poly = FACTOR.trim(poly[:])
    return FACTOR.scale(poly, 1 / poly[-1])


def rational_gcd(first, second):
    first = FACTOR.trim(first[:])
    second = FACTOR.trim(second[:])
    while second != [0]:
        _, remainder = FACTOR.divmod_poly(first, second)
        first, second = second, remainder
    return rational_monic(first)


def compose_square(poly):
    out = []
    for coefficient in poly:
        out.extend([coefficient, 0])
    return out[:-1]


def t_as_w_polynomial(c_value):
    out = [Fraction(0)] * 7
    for (i, j, k), coefficient in zip(
        FACTOR.T_MONOMIALS,
        FACTOR.T_COEFFICIENTS,
        strict=True,
    ):
        del i
        out[k] += Fraction(coefficient) * Fraction(c_value) ** j
    return FACTOR.trim(out)


def t_discriminant_value(c_value):
    poly = t_as_w_polynomial(c_value)
    degree = len(poly) - 1
    sign = -1 if degree * (degree - 1) // 2 & 1 else 1
    return (
        sign
        * FACTOR.resultant(poly, FACTOR.derivative(poly))
        / poly[-1]
    )


def reconstruct_t_discriminant():
    values = [t_discriminant_value(c_value) for c_value in range(51)]
    return FACTOR.interpolate_consecutive(values)


# ---------------------------------------------------------------------------
# Finite-field Rabin certificates


def mod_trim(poly, prime):
    poly = [value % prime for value in poly]
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def mod_add(first, second, prime):
    out = [0] * max(len(first), len(second))
    for index, value in enumerate(first):
        out[index] = (out[index] + value) % prime
    for index, value in enumerate(second):
        out[index] = (out[index] + value) % prime
    return mod_trim(out, prime)


def mod_sub(first, second, prime):
    return mod_add(first, [(-value) % prime for value in second], prime)


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


def prime_divisors(integer):
    out = []
    divisor = 2
    while divisor * divisor <= integer:
        if integer % divisor == 0:
            out.append(divisor)
        while integer % divisor == 0:
            integer //= divisor
        divisor += 1
    if integer > 1:
        out.append(integer)
    return out


def rabin_irreducible(poly, prime):
    poly = mod_trim(poly[:], prime)
    degree = len(poly) - 1
    variable = [0, 1]
    for divisor in prime_divisors(degree):
        frobenius = mod_power(
            variable,
            prime ** (degree // divisor),
            poly,
            prime,
        )
        if len(mod_gcd(poly, mod_sub(frobenius, variable, prime), prime)) > 1:
            return False
    final = mod_power(variable, prime**degree, poly, prime)
    return mod_reduce(mod_sub(final, variable, prime), poly, prime) == [0]


# ---------------------------------------------------------------------------
# Exact quotient field Q[C]/(P(C^2))


def q_trim(poly):
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def q_add(first, second):
    out = [Fraction(0)] * max(len(first), len(second))
    for index, value in enumerate(first):
        out[index] += value
    for index, value in enumerate(second):
        out[index] += value
    return q_trim(out)


def q_neg(poly):
    return [-value for value in poly]


def q_sub(first, second):
    return q_add(first, q_neg(second))


def q_mul(first, second):
    out = [Fraction(0)] * (len(first) + len(second) - 1)
    for i, left in enumerate(first):
        for j, right in enumerate(second):
            out[i + j] += left * right
    return q_trim(out)


def q_divmod(numerator, denominator):
    numerator = q_trim(numerator[:])
    denominator = q_trim(denominator[:])
    quotient = [Fraction(0)] * max(
        1, len(numerator) - len(denominator) + 1
    )
    while numerator != [0] and len(numerator) >= len(denominator):
        shift = len(numerator) - len(denominator)
        coefficient = numerator[-1] / denominator[-1]
        quotient[shift] = coefficient
        for index, value in enumerate(denominator):
            numerator[index + shift] -= coefficient * value
        q_trim(numerator)
    return q_trim(quotient), q_trim(numerator)


def q_extended_gcd(first, second):
    remainder_0, remainder_1 = first[:], second[:]
    bezout_0, bezout_1 = [Fraction(1)], [Fraction(0)]
    while remainder_1 != [0]:
        quotient, remainder = q_divmod(remainder_0, remainder_1)
        remainder_0, remainder_1 = remainder_1, remainder
        bezout_0, bezout_1 = (
            bezout_1,
            q_sub(bezout_0, q_mul(quotient, bezout_1)),
        )
    return remainder_0, bezout_0


class NumberField:
    """The field Q[C]/(P(C^2)), represented in the power basis."""

    def __init__(self, polynomial_in_c_squared):
        modulus = [Fraction(value) for value in polynomial_in_c_squared]
        self.modulus = [
            value / modulus[-1] for value in modulus
        ]

    def element(self, value):
        if isinstance(value, (int, Fraction)):
            value = [Fraction(value)]
        _, remainder = q_divmod(
            [Fraction(entry) for entry in value],
            self.modulus,
        )
        return tuple(remainder)

    def zero(self):
        return self.element(0)

    def one(self):
        return self.element(1)

    def add(self, first, second):
        return self.element(q_add(list(first), list(second)))

    def neg(self, value):
        return tuple(-entry for entry in value)

    def sub(self, first, second):
        return self.add(first, self.neg(second))

    def mul(self, first, second):
        return self.element(q_mul(list(first), list(second)))

    def inverse(self, value):
        gcd, bezout = q_extended_gcd(list(value), self.modulus)
        if len(gcd) != 1 or gcd[0] == 0:
            raise ZeroDivisionError("nonunit in advertised quotient field")
        return self.element([entry / gcd[0] for entry in bezout])

    def div(self, first, second):
        return self.mul(first, self.inverse(second))

    def scale(self, value, scalar):
        return self.mul(value, self.element(scalar))


# Polynomials over NumberField. They are used for W and for y.


def k_trim(poly, field):
    while len(poly) > 1 and poly[-1] == field.zero():
        poly.pop()
    return poly


def k_add(first, second, field):
    out = [field.zero()] * max(len(first), len(second))
    for index, value in enumerate(first):
        out[index] = field.add(out[index], value)
    for index, value in enumerate(second):
        out[index] = field.add(out[index], value)
    return k_trim(out, field)


def k_scale(poly, scalar, field):
    return k_trim([field.mul(value, scalar) for value in poly], field)


def k_mul(first, second, field):
    out = [field.zero()] * (len(first) + len(second) - 1)
    for i, left in enumerate(first):
        for j, right in enumerate(second):
            out[i + j] = field.add(
                out[i + j],
                field.mul(left, right),
            )
    return k_trim(out, field)


def k_power(poly, exponent, field):
    out = [field.one()]
    base = poly[:]
    while exponent:
        if exponent & 1:
            out = k_mul(out, base, field)
        base = k_mul(base, base, field)
        exponent //= 2
    return out


def k_divmod(numerator, denominator, field):
    numerator = k_trim(numerator[:], field)
    denominator = k_trim(denominator[:], field)
    quotient = [field.zero()] * max(
        1, len(numerator) - len(denominator) + 1
    )
    while (
        numerator != [field.zero()]
        and len(numerator) >= len(denominator)
    ):
        shift = len(numerator) - len(denominator)
        coefficient = field.div(numerator[-1], denominator[-1])
        quotient[shift] = coefficient
        for index, value in enumerate(denominator):
            numerator[index + shift] = field.sub(
                numerator[index + shift],
                field.mul(coefficient, value),
            )
        k_trim(numerator, field)
    return k_trim(quotient, field), k_trim(numerator, field)


def k_gcd(first, second, field):
    first = k_trim(first[:], field)
    second = k_trim(second[:], field)
    while second != [field.zero()]:
        _, remainder = k_divmod(first, second, field)
        first, second = second, remainder
    return k_scale(first, field.inverse(first[-1]), field)


def k_derivative(poly, field):
    return k_trim(
        [
            field.scale(poly[index], index)
            for index in range(1, len(poly))
        ],
        field,
    )


def k_evaluate(poly, value, field):
    out = field.zero()
    for coefficient in reversed(poly):
        out = field.add(field.mul(out, value), coefficient)
    return out


def t_partial(field, c_order=0, w_order=0):
    c_generator = field.element([0, 1])
    c_powers = [field.one()]
    for _ in range(10):
        c_powers.append(field.mul(c_powers[-1], c_generator))
    out = [field.zero()] * (7 - w_order)
    for (i, j, k), coefficient in zip(
        FACTOR.T_MONOMIALS,
        FACTOR.T_COEFFICIENTS,
        strict=True,
    ):
        del i
        if j < c_order or k < w_order:
            continue
        multiplier = coefficient
        for offset in range(c_order):
            multiplier *= j - offset
        for offset in range(w_order):
            multiplier *= k - offset
        out[k - w_order] = field.add(
            out[k - w_order],
            field.scale(c_powers[j - c_order], multiplier),
        )
    return k_trim(out, field)


def critical_section_data(field):
    t_poly = t_partial(field)
    t_w = t_partial(field, w_order=1)
    common = k_gcd(t_poly, t_w, field)
    if len(common) != 2 or common[-1] != field.one():
        raise RuntimeError("T/T_W gcd is not monic linear")
    w_value = field.neg(common[0])

    t_c = t_partial(field, c_order=1)
    _, remainder = k_divmod(t_c, common, field)
    if remainder != [field.zero()]:
        raise RuntimeError("vertical critical point is not singular")

    t_ww = k_evaluate(t_partial(field, w_order=2), w_value, field)
    t_cw = k_evaluate(
        t_partial(field, c_order=1, w_order=1),
        w_value,
        field,
    )
    t_cc = k_evaluate(t_partial(field, c_order=2), w_value, field)
    if t_ww == field.zero():
        raise RuntimeError("singular point lost the W-Morse direction")

    w_prime = field.neg(field.div(t_cw, t_ww))
    phi_2 = field.add(t_cc, field.mul(t_cw, w_prime))
    phi_3 = k_evaluate(t_partial(field, c_order=3), w_value, field)
    phi_3 = field.add(
        phi_3,
        field.scale(
            field.mul(
                k_evaluate(
                    t_partial(field, c_order=2, w_order=1),
                    w_value,
                    field,
                ),
                w_prime,
            ),
            3,
        ),
    )
    phi_3 = field.add(
        phi_3,
        field.scale(
            field.mul(
                k_evaluate(
                    t_partial(field, c_order=1, w_order=2),
                    w_value,
                    field,
                ),
                field.mul(w_prime, w_prime),
            ),
            3,
        ),
    )
    phi_3 = field.add(
        phi_3,
        field.mul(
            k_evaluate(t_partial(field, w_order=3), w_value, field),
            field.mul(field.mul(w_prime, w_prime), w_prime),
        ),
    )
    return common, w_value, t_ww, phi_2, phi_3


def residual_spectral_branch(field, w_value):
    c_generator = field.element([0, 1])
    central_d = Fraction(25, 126)
    cubic_a = field.element(-26040609)
    cubic_b = [
        field.element(49601160),
        field.zero(),
        field.element(1607445),
    ]
    cubic_c = [
        field.element(
            Fraction(-20995200) - Fraction(52907904) * central_d
        ),
        field.zero(),
        field.element(-2857680),
        field.zero(),
        field.element(-138915),
    ]
    cubic_d = [
        field.element(Fraction(33592320) * central_d),
        field.sub(
            field.scale(w_value, -5878656),
            field.scale(c_generator, 5598720),
        ),
        field.element(
            Fraction(777600) + Fraction(1959552) * central_d
        ),
        field.scale(c_generator, -435456),
        field.element(78120),
        field.zero(),
        field.element(1127),
    ]
    terms = [
        k_mul(
            k_power(cubic_b, 2, field),
            k_power(cubic_c, 2, field),
            field,
        ),
        k_scale(
            k_power(cubic_c, 3, field),
            field.scale(cubic_a, -4),
            field,
        ),
        k_scale(
            k_mul(k_power(cubic_b, 3, field), cubic_d, field),
            field.element(-4),
            field,
        ),
        k_scale(
            k_power(cubic_d, 2, field),
            field.scale(field.mul(cubic_a, cubic_a), -27),
            field,
        ),
        k_scale(
            k_mul(k_mul(cubic_b, cubic_c, field), cubic_d, field),
            field.scale(cubic_a, 18),
            field,
        ),
    ]
    branch = [field.zero()]
    for term in terms:
        branch = k_add(branch, term, field)
    if (
        len(branch) != 13
        or branch[0] != field.zero()
        or branch[1] != field.zero()
    ):
        raise RuntimeError("central spectral branch lost universal y^2")
    return branch[2:]


def squarefree_profile(poly, field):
    repeated = k_gcd(poly, k_derivative(poly, field), field)
    layer, remainder = k_divmod(poly, repeated, field)
    if remainder != [field.zero()]:
        raise RuntimeError("squarefree decomposition failed")
    profile = {}
    multiplicity = 1
    while layer != [field.one()]:
        overlap = k_gcd(layer, repeated, field)
        exact, remainder = k_divmod(layer, overlap, field)
        if remainder != [field.zero()]:
            raise RuntimeError("squarefree layer division failed")
        if exact != [field.one()]:
            profile[multiplicity] = len(exact) - 1
        layer = overlap
        repeated, remainder = k_divmod(repeated, overlap, field)
        if remainder != [field.zero()]:
            raise RuntimeError("repeated layer division failed")
        multiplicity += 1
    if repeated != [field.one()]:
        profile[multiplicity] = len(repeated) - 1
    return profile


def quotient_audit(name, primitive_factor):
    modulus = compose_square(primitive_factor)
    field = NumberField(modulus)
    common, w_value, t_ww, phi_2, phi_3 = critical_section_data(field)

    c_generator = field.element([0, 1])
    r_value = field.add(
        field.scale(c_generator, 20),
        field.scale(w_value, 21),
    )
    s_value = field.add(
        field.add(
            field.element(2888),
            field.scale(field.mul(c_generator, c_generator), 108864),
        ),
        field.add(
            field.scale(field.mul(c_generator, w_value), 571536),
            field.scale(field.mul(w_value, w_value), 750141),
        ),
    )
    if r_value == field.zero() or s_value == field.zero():
        raise RuntimeError(f"{name} singular stratum meets R*S=0 identically")

    residual = residual_spectral_branch(field, w_value)
    expected_central_residual = field.scale(
        field.mul(r_value, r_value),
        -27 * 26040609**2 * 279936**2,
    )
    if residual[0] != expected_central_residual:
        raise RuntimeError(
            f"{name} residual branch collides with the central fibre"
        )
    repeated = k_gcd(residual, k_derivative(residual, field), field)
    repeated_self = k_gcd(
        repeated,
        k_derivative(repeated, field),
        field,
    )
    profile = squarefree_profile(residual, field)

    return {
        "field": field,
        "common": common,
        "w": w_value,
        "t_ww": t_ww,
        "phi_2": phi_2,
        "phi_3": phi_3,
        "residual_degree": len(residual) - 1,
        "repeated_degree": len(repeated) - 1,
        "repeated_self_degree": len(repeated_self) - 1,
        "profile": profile,
    }


def main():
    # Exact discriminant reconstruction and squarefree decomposition.
    discriminant_c = reconstruct_t_discriminant()
    if len(discriminant_c) - 1 != 50:
        raise RuntimeError("Disc_W(T) lost degree 50 in C")
    if any(
        coefficient
        for exponent, coefficient in enumerate(discriminant_c)
        if exponent & 1
    ):
        raise RuntimeError("Disc_W(T) is not even in C")
    discriminant_x = [
        discriminant_c[2 * exponent]
        for exponent in range(26)
    ]
    expected = FACTOR.mul(
        FACTOR.power([Fraction(value) for value in P8], 2),
        FACTOR.power([Fraction(value) for value in P3], 3),
    )
    scalar = discriminant_x[-1] / expected[-1]
    if scalar == 0 or discriminant_x != FACTOR.scale(expected, scalar):
        raise RuntimeError("P8^2*P3^3 discriminant identity failed")

    p8_q = [Fraction(value) for value in P8]
    p3_q = [Fraction(value) for value in P3]
    if (
        len(rational_gcd(p8_q, FACTOR.derivative(p8_q))) != 1
        or len(rational_gcd(p3_q, FACTOR.derivative(p3_q))) != 1
        or len(rational_gcd(p8_q, p3_q)) != 1
    ):
        raise RuntimeError("P8/P3 squarefree-coprime check failed")

    p8_c = compose_square(P8)
    p3_c = compose_square(P3)
    if not rabin_irreducible(p8_c, 173):
        raise RuntimeError("P8(C^2) mod 173 is reducible")
    if not rabin_irreducible(p3_c, 11):
        raise RuntimeError("P3(C^2) mod 11 is reducible")

    p8_report = quotient_audit("P8", P8)
    p3_report = quotient_audit("P3", P3)

    if (
        p8_report["phi_2"] == p8_report["field"].zero()
        or p8_report["profile"] != {1: 6, 2: 2}
        or p8_report["repeated_degree"] != 2
        or p8_report["repeated_self_degree"] != 0
    ):
        raise RuntimeError("P8 node/residual profile failed")
    if (
        p3_report["phi_2"] != p3_report["field"].zero()
        or p3_report["phi_3"] == p3_report["field"].zero()
        or p3_report["profile"] != {1: 7, 3: 1}
        or p3_report["repeated_degree"] != 2
        or p3_report["repeated_self_degree"] != 1
    ):
        raise RuntimeError("P3 cusp/residual profile failed")

    expected_p3_selector = p3_report["field"].element(
        [
            0,
            Fraction(965, 714),
            0,
            Fraction(-104733, 3400),
            0,
            Fraction(-85562001, 53125),
        ]
    )
    if p3_report["common"][0] != expected_p3_selector:
        raise RuntimeError("P3 singular selector changed")

    if (
        p3_report["t_ww"][0]
        != Fraction(611393420358000000000000, 23)
        or p3_report["phi_3"][1]
        != Fraction(5660096667257733120000000000, 8993)
    ):
        raise RuntimeError("P3 local nonvanishing controls changed")

    if (
        p8_report["phi_2"][0]
        != Fraction(
            -22372418489981818832620603505383870089495590400000000000,
            232419884438325071097488244108192949,
        )
    ):
        raise RuntimeError("P8 Hessian control changed")

    # Analytic coefficient controls used after the exact atlas.
    central_fibre = [
        Fraction(46656000, 7),
        Fraction(-31492800),
        Fraction(49601160),
        Fraction(-26040609),
    ]
    central_expected = FACTOR.scale(
        FACTOR.power(
            [Fraction(-40, 63), Fraction(1)],
            3,
        ),
        Fraction(-26040609),
    )
    if central_fibre != central_expected:
        raise RuntimeError("central total-ramification identity failed")
    if (
        -5598720 != -279936 * 20
        or -5878656 != -279936 * 21
    ):
        raise RuntimeError("central G_y=-279936*R identity failed")

    infinity_cubic = [
        Fraction(1127),
        Fraction(-138915),
        Fraction(1607445),
        Fraction(-26040609),
    ]
    infinity_discriminant = (
        (-1) ** 3
        * FACTOR.resultant(
            infinity_cubic,
            FACTOR.derivative(infinity_cubic),
        )
        / infinity_cubic[-1]
    )
    if infinity_discriminant != -153384762202971019112448:
        raise RuntimeError("infinity cubic discriminant changed")

    # If u=a*y^2+b*y+c is a polynomial root, its y^5 coefficient is
    # b*L'(a). Once b=0, its odd coefficients are the displayed C/W terms.
    if FACTOR.T_COEFFICIENTS[-1] != 256000000000000000:
        raise RuntimeError("T(1,0,0) control changed")

    print("JC2 DEGREE-18 CENTRAL T SINGULAR ATLAS EXACT AUDIT")
    print("Disc_W_T_degree_C=50")
    print("factorization=constant*P8(C^2)^2*P3(C^2)^3")
    print("P8_degree_X=8; P3_degree_X=3")
    print("P8(C^2)_irreducible_mod_173_degree=16")
    print("P3(C^2)_irreducible_mod_11_degree=6")
    print(
        "P8_stratum=16 singular points; local_type=node; "
        "residual=6_simple+2_double"
    )
    print(
        "P3_stratum=6 singular points; local_type=cusp; "
        "residual=7_simple+1_triple"
    )
    print("R_and_S_nonzero_on_both_singular_strata=True")
    print(
        "smooth_T_locus_control="
        "universal_discriminant_smooth_locus=>8_simple+1_double"
    )
    print(
        "central_fibre="
        "-26040609*(u-40/63)^3; G_y=-279936*R; "
        "delta_bar(0)=-27*26040609^2*(279936*R)^2"
    )
    print(
        "infinity_cubic_discriminant="
        "-153384762202971019112448"
    )
    print(
        "spectral_irreducibility_control="
        "y5=b*L'(a); b=0=>y3=-435456*C and y1=-279936*R; "
        "T(1,0,0)!=0"
    )
    print("genus_lower_bounds=smooth>=3; P8>=2; P3>=3")
    print("PASS")


if __name__ == "__main__":
    main()
