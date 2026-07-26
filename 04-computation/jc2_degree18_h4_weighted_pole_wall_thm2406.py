#!/usr/bin/env python3
"""Exact splitting-field wall certificate for THM-2406."""

from __future__ import annotations

import sympy as sp
from sympy.polys.groebnertools import groebner as low_level_groebner
from sympy.polys.rings import ring


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> None:
    slope_symbol = sp.symbols("slope_symbol")
    slope_polynomial = sp.Poly(
        slope_symbol**3 + 2940 * slope_symbol + 30184,
        slope_symbol,
    )
    require(slope_polynomial.is_irreducible, "slope cubic irreducible")
    slope_discriminant = sp.discriminant(slope_polynomial)
    require(
        slope_discriminant == -126247730112,
        "slope cubic discriminant",
    )
    require(
        slope_discriminant == -23 * 74088**2,
        "slope discriminant square class",
    )

    epsilon_symbol = sp.symbols("epsilon_symbol")
    minimal = sp.Poly(
        epsilon_symbol**6
        + 10 * epsilon_symbol**4
        + 25 * epsilon_symbol**2
        + 23,
        epsilon_symbol,
    )
    require(minimal.is_irreducible, "splitting-field polynomial irreducible")
    require(
        sp.gcd(minimal, minimal.diff()).degree() == 0,
        "splitting-field polynomial separable",
    )
    field = sp.QQ.alg_field_from_poly(minimal, alias="epsilon")
    epsilon = field.unit
    slope_denominator = (
        field.convert(3) * epsilon**2 + field.convert(5)
    )
    require(slope_denominator != field.zero, "slope denominator nonzero")
    slope_denominator_inverse = field.one / slope_denominator
    require(
        slope_denominator * slope_denominator_inverse == field.one,
        "slope denominator inverse",
    )
    alpha0 = field.convert(154) * slope_denominator_inverse
    alpha1 = -alpha0 / field.convert(2) + field.convert(21) * epsilon
    alpha_infinity = -alpha0 / field.convert(2) - field.convert(21) * epsilon
    for slope in (alpha0, alpha1, alpha_infinity):
        require(
            slope**3 + field.convert(2940) * slope + field.convert(30184)
            == field.zero,
            "slope root",
        )
    require(alpha0 + alpha1 + alpha_infinity == field.zero, "slope sum")
    require(
        alpha1 - alpha_infinity == field.convert(42) * epsilon,
        "scaled slope difference",
    )
    require(
        len({str(alpha0), str(alpha1), str(alpha_infinity)}) == 3,
        "distinct slope expressions",
    )
    print("H4 WEIGHTED-POLE WALL AUDIT -- exact characteristic zero")
    print("slope cubic irreducible; discriminant = -23*74088^2")
    print("splitting field h(epsilon)=epsilon^6+10epsilon^4+25epsilon^2+23")
    print("slope denominator 3epsilon^2+5 is a unit in K")
    print("ordered slopes: alpha0=154/(3epsilon^2+5),")
    print("alpha1/alpha_inf=-alpha0/2 +/- 21epsilon")

    base, xr, yr, zr = ring("x,y,z", field)
    fractions = base.to_field()
    spectral, kappa, B, C, D, W, t = ring(
        "kappa,B,C,D,W,t",
        fractions,
    )

    def ground(value):
        return spectral.domain_new(fractions(value))

    x = ground(xr)
    y = ground(yr)
    z = ground(zr)
    s0 = ground(alpha0)
    s1 = ground(alpha1)
    sinf = ground(alpha_infinity)

    A = spectral.one + x * t + y * t**2 + z * t**3
    denominator = t * (t - 1)
    A1 = spectral.one + x + y + z
    A1p = x + 2 * y + 3 * z
    c0 = s0
    c1 = 2 * s0 * x
    c6 = sinf * z**2
    c5 = 2 * sinf * z * y
    sigma = s1 * A1**2 - (c0 + c1 + c5 + c6)
    tau = 2 * s1 * A1 * A1p - (c1 + 5 * c5 + 6 * c6)
    r0 = (
        c0
        + c1 * t
        + (3 * sigma - tau) * t**2
        + (tau - 2 * sigma) * t**3
        + c5 * t**5
        + c6 * t**6
    )
    r = r0 + kappa * denominator**2
    Ph = (
        245 * A**4
        + 1890 * B * A**2 * denominator**2
        + (-24300 * B**2 + 122472 * D) * denominator**4
    )
    Qh = (
        539 * A**6
        + 11340 * B * A**4 * denominator**2
        + 183708 * C * A**3 * denominator**3
        + (72900 * B**2 - 367416 * D) * A**2 * denominator**4
        + (2361960 * B * C + 2480058 * W) * A * denominator**5
    )
    numerator = r**3 + 12 * Ph * r + 56 * Qh
    require(numerator.degree(t) == 16, "Hermite-reduced numerator degree")
    print("homogenized numerator degree after pole jets:", numerator.degree(t))

    def jet(mode: str, order: int):
        if mode == "zero":
            return numerator.coeff_wrt(t, order)
        if mode == "one":
            return sum(
                (
                    numerator.coeff_wrt(t, degree)
                    * sp.binomial(degree, order)
                    for degree in range(order, 19)
                ),
                spectral.zero,
            )
        if mode == "infinity":
            return numerator.coeff_wrt(t, 18 - order)
        raise ValueError(mode)

    def constant(poly):
        return poly.coeff(spectral.one)

    zero2 = jet("zero", 2)
    one2 = jet("one", 2)
    z00 = constant(zero2)
    z0k = zero2.coeff(kappa)
    z0b = zero2.coeff(B)
    z10 = constant(one2)
    z1k = one2.coeff(kappa)
    z1b = one2.coeff(B)
    determinant = z0k * z1b - z0b * z1k
    k_value = (-z00 * z1b + z0b * z10) / determinant
    b_value = (-z0k * z10 + z00 * z1k) / determinant
    determinant_quotient = determinant / (fractions(xr + yr + zr + 1) ** 4)
    require(
        determinant_quotient.numer.degree() == 0
        and determinant_quotient.denom.degree() == 0,
        "quadratic determinant is not a constant times A1^4",
    )
    determinant_constant = (
        determinant_quotient.numer.coeff(base.one)
        / determinant_quotient.denom.coeff(base.one)
    )
    require(determinant_constant != field.zero, "quadratic determinant constant")
    require(
        determinant_constant * (field.one / determinant_constant)
        == field.one,
        "quadratic determinant constant inverse",
    )
    print("B-pivot determinant = nonzero constant * A3(1)^4")

    solved = {
        kappa: k_value,
        B: b_value,
    }

    def substitute(poly, mapping):
        out = poly
        for variable, value in mapping.items():
            out = out.subs(variable, spectral.domain_new(value))
        return out

    equations = [
        ("B_inf", constant(substitute(jet("infinity", 2), solved))),
    ]
    reconstruction_pivots = {}
    for variable, order, label in (
        (C, 3, "C"),
        (D, 4, "D"),
        (W, 5, "W"),
    ):
        zero_jet = substitute(jet("zero", order), solved)
        base_value = constant(zero_jet)
        coefficient = zero_jet.coeff(variable)
        require(coefficient != fractions.zero, f"{label} pivot")
        require(
            coefficient.numer.degree() == 0
            and coefficient.denom.degree() == 0,
            f"{label} pivot is not constant",
        )
        pivot_constant = (
            coefficient.numer.coeff(base.one)
            / coefficient.denom.coeff(base.one)
        )
        require(pivot_constant != field.zero, f"{label} pivot constant")
        require(
            pivot_constant * (field.one / pivot_constant) == field.one,
            f"{label} pivot inverse",
        )
        reconstruction_pivots[label] = pivot_constant
        solved[variable] = -base_value / coefficient
        for mode in ("one", "infinity"):
            value = constant(substitute(jet(mode, order), solved))
            equations.append((f"{label}_{mode}", value))
    equations.append(
        ("lock", constant(substitute(jet("zero", 6), solved)))
    )
    print("triangular spectral reconstruction: kappa,B,C,D,W")
    print(
        "constant reconstruction pivots are K-units:",
        ",".join(reconstruction_pivots),
    )

    # Lift the reduced numerators directly into K[eta,z,y,x], bypassing
    # expression-to-number-field conversion.
    gb_ring, eta, gz, gy, gx = ring(
        "eta,z,y,x",
        field,
        order="grevlex",
    )

    def lift_base(poly):
        return gb_ring.from_dict(
            {
                (0, monomial[2], monomial[1], monomial[0]): coefficient
                for monomial, coefficient in poly.items()
            }
        )

    # Four of the eight synchronization equations already form a wall-forcing
    # subsystem.  A full survivor satisfies this subsystem, so proving the wall
    # here loses nothing and avoids the much larger intermediate ideals made by
    # the unused matching equations.
    equation_by_name = dict(equations)
    forcing_names = ("lock", "B_inf", "C_one", "D_one")
    pole_unit = base.one + xr + yr + zr
    unit_factors = {
        "lock": base.one,
        "B_inf": zr**4,
        "C_one": pole_unit**3,
        "D_one": pole_unit**2,
    }

    def primitive_numerator(name):
        numerator = equation_by_name[name].numer
        quotients, remainder = numerator.div([unit_factors[name]])
        require(remainder == base.zero, f"{name} pole-unit division")
        return quotients[0]

    primitive_by_name = {
        name: primitive_numerator(name) for name in forcing_names
    }
    expected_profiles = {
        "lock": (6, 55),
        "B_inf": (2, 8),
        "C_one": (3, 14),
        "D_one": (4, 29),
    }
    for name in forcing_names:
        primitive = primitive_by_name[name]
        total_degree = max(sum(monomial) for monomial in primitive)
        print(
            "primitive profile",
            name,
            total_degree,
            len(primitive),
        )
        require(
            (total_degree, len(primitive)) == expected_profiles[name],
            f"{name} primitive profile",
        )
    print("localized unit factors:")
    print("  B_inf/a3^4 degree 2 terms 8")
    print("  C_one/A3(1)^3 degree 3 terms 14")
    print("  D_one/A3(1)^2 degree 4 terms 29")
    print("  sixth-order lock degree 6 terms 55")

    def pole_unit_denominator(poly):
        current = poly
        exponent = 0
        while current.degree() > 0:
            quotients, remainder = current.div([pole_unit])
            require(remainder == base.zero, "unexpected reconstruction denominator")
            current = quotients[0]
            exponent += 1
        constant_value = current.coeff(base.one)
        require(constant_value != field.zero, "zero reconstruction denominator")
        return exponent

    denominator_exponents = {
        str(variable): pole_unit_denominator(value.denom)
        for variable, value in solved.items()
    }
    print("reconstruction denominator A3(1)-exponents:", denominator_exponents)

    equation_polys = [
        lift_base(primitive_by_name[name]) for name in forcing_names
    ]
    saturation = eta * gz * (1 + gx + gy + gz) - 1
    basis = low_level_groebner([saturation, *equation_polys], gb_ring)
    print(
        "GB",
        ",".join(forcing_names),
        "size=", len(basis),
        flush=True,
    )
    require(len(basis) == 4, "triangular basis size")

    def contains_monomial(poly, monomial):
        return monomial in dict(poly.items())

    q_poly = next(
        poly
        for poly in basis
        if contains_monomial(poly, (0, 0, 0, 2))
    )
    eta_poly = next(
        poly
        for poly in basis
        if contains_monomial(poly, (1, 0, 0, 0))
    )
    z_poly = next(
        poly
        for poly in basis
        if contains_monomial(poly, (0, 1, 0, 0))
    )
    y_poly = next(
        poly
        for poly in basis
        if contains_monomial(poly, (0, 0, 1, 0))
    )
    require(
        len({str(poly) for poly in basis}) == 4,
        "triangular basis identification",
    )
    require(
        dict(eta_poly.items())[(1, 0, 0, 0)] == field.one,
        "eta basis coefficient is not one",
    )
    require(
        set(dict(q_poly.items())) == {
            (0, 0, 0, 2),
            (0, 0, 0, 1),
            (0, 0, 0, 0),
        },
        "quadratic quotient support",
    )
    for poly, leading in (
        (eta_poly, (1, 0, 0, 0)),
        (z_poly, (0, 1, 0, 0)),
        (y_poly, (0, 0, 1, 0)),
    ):
        require(
            set(dict(poly.items())) <= {
                leading,
                (0, 0, 0, 1),
                (0, 0, 0, 0),
            },
            "linear triangular support",
        )

    def field_text(value):
        return str(sp.factor(field.to_sympy(value)))

    print("triangular quotient over K:")
    print("  eta-leading coefficient=1, inverse=1 in K[x]/(q)")
    for label, poly in (
        ("q(x)", q_poly),
        ("eta", eta_poly),
        ("z", z_poly),
        ("y", y_poly),
    ):
        coefficient_map = dict(poly.items())
        linear = coefficient_map.get((0, 0, 0, 1), field.zero)
        constant_value = coefficient_map.get(
            (0, 0, 0, 0),
            field.zero,
        )
        if label == "q(x)":
            print(
                "  q=x^2+(" + field_text(linear) + ")*x+("
                + field_text(constant_value) + ")",
            )
        else:
            print(
                "  " + label + "+(" + field_text(linear) + ")*x+("
                + field_text(constant_value) + ")=0",
            )

    def reduce_fraction(value):
        numerator_remainder = lift_base(value.numer).div(basis)[1]
        denominator_remainder = lift_base(value.denom).div(basis)[1]
        # The final quotient is expected to be quadratic and a field.  For the
        # wall check it is enough to cross-multiply: a zero numerator remainder
        # proves vanishing wherever the original denominator is a unit.
        return numerator_remainder, denominator_remainder

    wall = (
        fractions(126) * solved[D]
        - fractions(25) * solved[B] ** 2
    )
    deep_wall = (
        fractions(20) * solved[B] * solved[C]
        + fractions(21) * solved[W]
    )
    wall_remainder, _ = reduce_fraction(wall)
    deep_remainder, _ = reduce_fraction(deep_wall)
    require(wall_remainder == 0, "common-root wall not forced")
    require(deep_remainder == 0, "deep wall not forced")
    print("normal form 126D-25B^2: 0")
    print("normal form 20BC+21W: 0")

    omitted = []
    for name, value in equations:
        if name in forcing_names:
            continue
        remainder, _ = reduce_fraction(value)
        omitted.append((name, remainder == 0))
    print(
        "omitted synchronization normal forms:",
        " ".join(
            f"{name}_vanishes={int(vanishes)}"
            for name, vanishes in omitted
        ),
    )

    def coefficient_mod_prime(value, prime, epsilon_residue):
        total = 0
        for coefficient in value.to_list():
            numerator_value = int(coefficient.numerator) % prime
            denominator_value = int(coefficient.denominator) % prime
            require(denominator_value != 0, "bad finite-place denominator")
            total = (
                total * epsilon_residue
                + numerator_value * pow(denominator_value, -1, prime)
            ) % prime
        return total

    def polynomial_mod_prime(poly, prime, epsilon_residue):
        return {
            monomial: coefficient_mod_prime(
                coefficient,
                prime,
                epsilon_residue,
            )
            for monomial, coefficient in poly.items()
            if coefficient_mod_prime(
                coefficient,
                prime,
                epsilon_residue,
            )
            != 0
        }

    finite_controls = (
        (
            59,
            1,
            (34, 4, 21),
            {
                "q": {
                    (0, 0, 0, 2): 1,
                    (0, 0, 0, 1): 18,
                    (0, 0, 0, 0): 24,
                },
                "eta": {
                    (1, 0, 0, 0): 1,
                    (0, 0, 0, 1): 9,
                    (0, 0, 0, 0): 19,
                },
                "z": {
                    (0, 1, 0, 0): 1,
                    (0, 0, 0, 1): 26,
                    (0, 0, 0, 0): 52,
                },
                "y": {
                    (0, 0, 1, 0): 1,
                    (0, 0, 0, 1): 7,
                    (0, 0, 0, 0): 40,
                },
            },
        ),
        (
            101,
            96,
            (60, 67, 75),
            {
                "q": {
                    (0, 0, 0, 2): 1,
                    (0, 0, 0, 1): 45,
                    (0, 0, 0, 0): 14,
                },
                "eta": {
                    (1, 0, 0, 0): 1,
                    (0, 0, 0, 1): 25,
                    (0, 0, 0, 0): 79,
                },
                "z": {
                    (0, 1, 0, 0): 1,
                    (0, 0, 0, 1): 64,
                    (0, 0, 0, 0): 27,
                },
                "y": {
                    (0, 0, 1, 0): 1,
                    (0, 0, 0, 1): 74,
                    (0, 0, 0, 0): 10,
                },
            },
        ),
    )
    for prime, epsilon_residue, slopes, expected in finite_controls:
        require(
            int(minimal.eval(epsilon_residue)) % prime == 0,
            "finite-place epsilon root",
        )
        computed_slopes = tuple(
            coefficient_mod_prime(slope, prime, epsilon_residue)
            for slope in (alpha0, alpha1, alpha_infinity)
        )
        require(computed_slopes == slopes, "finite-place slope ordering")
        actual = {
            "q": polynomial_mod_prime(q_poly, prime, epsilon_residue),
            "eta": polynomial_mod_prime(eta_poly, prime, epsilon_residue),
            "z": polynomial_mod_prime(z_poly, prime, epsilon_residue),
            "y": polynomial_mod_prime(y_poly, prime, epsilon_residue),
        }
        require(actual == expected, "finite-place triangular basis")
        q_coefficients = expected["q"]
        roots = [
            value
            for value in range(prime)
            if (
                value**2
                + q_coefficients[(0, 0, 0, 1)] * value
                + q_coefficients[(0, 0, 0, 0)]
            )
            % prime
            == 0
        ]
        print(
            f"finite-place control p={prime} epsilon={epsilon_residue} "
            f"slopes={slopes} q-roots={roots}",
        )

    def base_terms_mod_prime(poly, prime, epsilon_residue):
        return [
            (
                monomial,
                coefficient_mod_prime(
                    coefficient,
                    prime,
                    epsilon_residue,
                ),
            )
            for monomial, coefficient in poly.items()
        ]

    def evaluate_base_terms(terms, values, prime):
        return sum(
            coefficient
            * pow(values[0], monomial[0], prime)
            * pow(values[1], monomial[1], prime)
            * pow(values[2], monomial[2], prime)
            for monomial, coefficient in terms
        ) % prime

    prime = 59
    epsilon_residue = 1
    prelock_terms = [
        base_terms_mod_prime(
            primitive_by_name[name],
            prime,
            epsilon_residue,
        )
        for name in ("B_inf", "C_one", "D_one")
    ]
    lock_terms = base_terms_mod_prime(
        primitive_by_name["lock"],
        prime,
        epsilon_residue,
    )
    wall_numerator_terms = base_terms_mod_prime(
        wall.numer,
        prime,
        epsilon_residue,
    )
    wall_denominator_terms = base_terms_mod_prime(
        wall.denom,
        prime,
        epsilon_residue,
    )
    hostile = None
    for x_value in range(prime):
        if hostile is not None:
            break
        for y_value in range(prime):
            if hostile is not None:
                break
            for z_value in range(1, prime):
                values = (x_value, y_value, z_value)
                if (1 + x_value + y_value + z_value) % prime == 0:
                    continue
                if any(
                    evaluate_base_terms(terms, values, prime) != 0
                    for terms in prelock_terms
                ):
                    continue
                wall_denominator_value = evaluate_base_terms(
                    wall_denominator_terms,
                    values,
                    prime,
                )
                if wall_denominator_value == 0:
                    continue
                wall_numerator_value = evaluate_base_terms(
                    wall_numerator_terms,
                    values,
                    prime,
                )
                if wall_numerator_value == 0:
                    continue
                lock_value = evaluate_base_terms(
                    lock_terms,
                    values,
                    prime,
                )
                if lock_value == 0:
                    continue
                hostile = (
                    values,
                    wall_numerator_value
                    * pow(wall_denominator_value, -1, prime)
                    % prime,
                    lock_value,
                )
                break
    require(hostile is not None, "pre-lock off-wall hostile")
    print(
        "sharp lock hostile p=59 (x,y,z),J,lock:",
        hostile,
    )

    require(122472 * 25 == 24300 * 126, "P(0) wall arithmetic")
    require(367416 * 25 == 72900 * 126, "Q y2 wall arithmetic")
    require(2480058 * 20 == 2361960 * 21, "Q y wall arithmetic")
    print("wall handoff: P and Q share y=0; full wall closed by THM-2345")
    print("VERDICT: exact H4 pole subsystem lies on the closed deep wall")


if __name__ == "__main__":
    main()
