#!/usr/bin/env python3
"""Exact symbolic/finite companion for THM-3271."""

from itertools import product

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def centered_norm_polynomial(x, a3, a2, a1):
    return sp.cancel(
        -(
            20 * x**3
            + 15 * a3 * x**2
            + (18 * a2 - 3 * a3**2) * x
            + 27 * a1
            - 9 * a2 * a3
            + 2 * a3**3
        )
        / 27
    )


def quartic_centered_norm(f, x):
    coefficients = sp.Poly(f, x, domain=sp.QQ).all_coeffs()
    require(len(coefficients) == 5 and coefficients[0] == 1,
            "control quartic must be monic")
    _, a3, a2, a1, _ = coefficients
    return centered_norm_polynomial(x, a3, a2, a1)


def main():
    x, z = sp.symbols("x z")
    a3, a2, a1, a0 = sp.symbols("a3 a2 a1 a0")

    # Synthetic division and depression of the complementary cubic.
    b2 = a3 + x
    b1 = a2 + a3 * x + x**2
    b0 = a1 + a2 * x + a3 * x**2 + x**3
    mu = -b2 / 3
    synthetic_norm = sp.expand(-(mu**3 + b2 * mu**2 + b1 * mu + b0))
    universal_norm = centered_norm_polynomial(x, a3, a2, a1)
    require(sp.expand(synthetic_norm - universal_norm) == 0,
            "universal centered-norm formula failed")
    require(sp.Poly(universal_norm, x).degree() == 3,
            "universal packet lost its cubic leading term")
    require(sp.Poly(universal_norm, x).LC() == -sp.Rational(20, 27),
            "universal packet leading coefficient is wrong")

    # Rootwise proof and affine covariance.  The formula must equal the
    # centered product of the other three roots at every tautological sheet.
    roots = sp.symbols("r0:4")
    u, v = sp.symbols("u v", nonzero=True)
    root_quartic = sp.Poly(sp.prod(x - root for root in roots), x)
    _, ra3, ra2, ra1, _ = root_quartic.all_coeffs()
    root_formula = centered_norm_polynomial(x, ra3, ra2, ra1)
    affine_checks = 0
    for index, root in enumerate(roots):
        others = [roots[j] for j in range(4) if j != index]
        mean = sum(others) / 3
        product_norm = sp.expand(sp.prod(other - mean for other in others))
        require(sp.expand(root_formula.subs(x, root) - product_norm) == 0,
                "rootwise packet formula failed")

        transformed = [u * other + v for other in others]
        transformed_mean = sum(transformed) / 3
        transformed_norm = sp.expand(
            sp.prod(other - transformed_mean for other in transformed)
        )
        require(sp.expand(transformed_norm - u**3 * product_norm) == 0,
                "affine covariance failed")
        affine_checks += 1

    # A cyclic cubic 1+3 control.  The cubic has square discriminant 81,
    # and the packet characteristic polynomial has a unique simple rational
    # root.  Its spectral projector is the singleton idempotent.
    cubic = x**3 - 3 * x + 1
    cyclic_quartic = sp.expand((x - 2) * cubic)
    require(sp.discriminant(cubic, x) == 81,
            "cyclic cubic control has wrong discriminant")
    cyclic_norm = quartic_centered_norm(cyclic_quartic, x)
    characteristic = sp.Poly(
        sp.resultant(cyclic_quartic, z - cyclic_norm, x), z, domain=sp.QQ
    ).as_expr()
    fixed_norm = sp.cancel(cyclic_norm.subs(x, 2))
    require(fixed_norm == -1, "fixed component of cyclic control is wrong")
    moving_characteristic = sp.cancel(characteristic / (z - fixed_norm))
    require(sp.Poly(moving_characteristic, z, domain=sp.QQ).is_irreducible,
            "moving packet control should be irreducible cubic")
    spectral_denominator = sp.cancel(moving_characteristic.subs(z, fixed_norm))
    require(spectral_denominator == sp.Rational(424, 729),
            "spectral denominator is wrong")
    projector = sp.cancel(
        moving_characteristic.subs(z, cyclic_norm) / spectral_denominator
    )
    projector_remainder = sp.rem(
        sp.Poly(projector, x, domain=sp.QQ),
        sp.Poly(cyclic_quartic, x, domain=sp.QQ),
    ).as_expr()
    require(sp.expand(projector_remainder - (x**3 / 3 - x + sp.Rational(1, 3))) == 0,
            "spectral projector formula is wrong")
    require(sp.rem(sp.Poly(projector - 1, x), sp.Poly(x - 2, x)).as_expr() == 0,
            "projector is not one on the singleton factor")
    require(sp.rem(sp.Poly(projector, x), sp.Poly(cubic, x)).as_expr() == 0,
            "projector is not zero on the cubic factor")
    require(
        sp.rem(
            sp.Poly(projector**2 - projector, x, domain=sp.QQ),
            sp.Poly(cyclic_quartic, x, domain=sp.QQ),
        ).as_expr()
        == 0,
        "spectral projector is not idempotent",
    )

    # Exhaust a bounded hostile bank of irreducible cubic factors.  The
    # degree-three packet never identifies the rational singleton with a
    # conjugate moving value.
    separation_controls = 0
    for alpha, c2, c1, c0 in product(range(-2, 3), repeat=4):
        moving = sp.Poly(x**3 + c2 * x**2 + c1 * x + c0, x, domain=sp.QQ)
        if not moving.is_irreducible:
            continue
        f = sp.expand((x - alpha) * moving.as_expr())
        packet = sp.Poly(quartic_centered_norm(f, x), x, domain=sp.QQ)
        collision = sp.Poly(packet.as_expr() - packet.eval(alpha), x, domain=sp.QQ)
        require(sp.gcd(moving, collision).degree() == 0,
                "singleton/moving packet collision found")
        separation_controls += 1
    require(separation_controls == 360, "bounded 1+3 control count changed")

    # Sharp warning: without an irreducible 1+3 decomposition the scalar
    # characteristic polynomial can collide, even though the algebra-valued
    # packet still retains the sheets.
    collision_quartic = sp.expand((x**2 - 9) * (x**2 - 1))
    collision_norm = quartic_centered_norm(collision_quartic, x)
    collision_values = tuple(
        sp.cancel(collision_norm.subs(x, root)) for root in (3, -3, 1, -1)
    )
    require(collision_values == (0, 0, sp.Rational(160, 27), sp.Rational(-160, 27)),
            "collision hostile values changed")
    collision_characteristic = sp.factor(
        sp.resultant(collision_quartic, z - collision_norm, x)
    )
    require(
        sp.expand(
            collision_characteristic
            - z**2 * (27 * z - 160) * (27 * z + 160) / 729
        )
        == 0,
        "collision hostile characteristic polynomial changed",
    )

    print("THM3271 universal quartic centered-norm packet exact companion")
    print("formula_degree=3 leading=-20/27 synthetic_identity=True")
    print(f"rootwise_checks=4 affine_cube_covariance_checks={affine_checks}")
    print("cyclic_control: cubic_discriminant=81 fixed_norm=-1 unique_simple_root=True")
    print("spectral_denominator=424/729 projector=x^3/3-x+1/3 idempotent=True")
    print(f"irreducible_1plus3_separation_controls={separation_controls}")
    print("collision_hostile_values=(0,0,160/27,-160/27) scalar_packet_not_sheet_injective=True")


if __name__ == "__main__":
    main()
