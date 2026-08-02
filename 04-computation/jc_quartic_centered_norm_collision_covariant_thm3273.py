#!/usr/bin/env python3
"""Exact symbolic companion for THM-3273."""

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def main():
    x, z, p, q, r, u = sp.symbols("x z p q r u")
    quartic = x**4 + p * x**2 + q * x + r
    packet_map = -(20 * x**3 + 18 * p * x + 27 * q) / sp.Integer(27)
    packet_polynomial = sp.Poly(
        sp.resultant(quartic, z - packet_map, x), z
    ).as_expr()

    quartic_discriminant = sp.factor(sp.discriminant(quartic, x))
    packet_discriminant = sp.factor(sp.discriminant(packet_polynomial, z))
    collision_covariant = (
        9 * p**6
        - 1900 * p**4 * r
        + 1000 * p**3 * q**2
        + 110000 * p**2 * r**2
        - 1000000 * r**3
    )
    expected_discriminant = sp.factor(
        sp.Rational(2**12, 3**36)
        * collision_covariant**2
        * quartic_discriminant
    )
    require(sp.expand(packet_discriminant - expected_discriminant) == 0,
            "packet discriminant factorization failed")

    # Direct six-pair proof.  Impose sum roots=0, express p,q,r by the
    # elementary symmetric functions, and multiply the six packet/root
    # difference quotients.
    r0, r1, r2 = sp.symbols("r0 r1 r2")
    r3 = -r0 - r1 - r2
    roots = (r0, r1, r2, r3)
    root_quartic = sp.Poly(sp.prod(x - root for root in roots), x)
    _, root_a3, root_p, root_q, root_r = root_quartic.all_coeffs()
    require(sp.expand(root_a3) == 0, "root parametrization is not depressed")
    pair_product = sp.Integer(1)
    pair_checks = 0
    for i in range(4):
        for j in range(i + 1, 4):
            left = roots[i]
            right = roots[j]
            quotient = -(
                20 * (left**2 + left * right + right**2) + 18 * root_p
            ) / sp.Integer(27)
            direct = sp.cancel(
                (
                    packet_map.subs({x: left, p: root_p, q: root_q})
                    - packet_map.subs({x: right, p: root_p, q: root_q})
                )
                / (left - right)
            )
            require(sp.expand(direct - quotient) == 0,
                    "pair quotient formula failed")
            pair_product *= quotient
            pair_checks += 1
    root_covariant = collision_covariant.subs(
        {p: root_p, q: root_q, r: root_r}
    )
    require(
        sp.expand(pair_product + sp.Rational(2**6, 3**18) * root_covariant)
        == 0,
        "six pair quotients do not multiply to the collision covariant",
    )

    # Affine weight.  Depression removes translation; scaling roots by u
    # sends (p,q,r) to (u^2 p,u^3 q,u^4 r).
    scaled_covariant = sp.expand(
        collision_covariant.subs({p: u**2 * p, q: u**3 * q, r: u**4 * r}, simultaneous=True)
    )
    require(sp.expand(scaled_covariant - u**12 * collision_covariant) == 0,
            "collision covariant has wrong affine weight")

    # The tame special fibre (a,m,m,m), depressed after translating m to
    # zero, has a unit collision covariant whenever 2,3,a-m are units.
    d = sp.symbols("d")
    y = sp.symbols("y")
    special_quartic = sp.expand(x**3 * (x - d))
    depressed_special = sp.Poly(
        sp.expand(special_quartic.subs(x, y + d / 4)), y
    )
    _, special_a3, special_p, special_q, special_r = depressed_special.all_coeffs()
    require(sp.expand(special_a3) == 0, "special quartic was not depressed")
    special_covariant = sp.factor(
        collision_covariant.subs(
            {p: special_p, q: special_q, r: special_r}
        )
    )
    require(special_covariant == sp.Rational(27, 8) * d**12,
            "tame special-fibre covariant is wrong")

    # Sharp scalar collision hostile and a noncollision positive control.
    hostile = {p: -10, q: 0, r: 9}
    hostile_covariant = collision_covariant.subs(hostile)
    hostile_discriminant = quartic_discriminant.subs(hostile)
    require(hostile_covariant == 0 and hostile_discriminant != 0,
            "separable scalar-packet collision hostile failed")
    hostile_packet = sp.factor(packet_polynomial.subs(hostile))
    require(
        sp.expand(hostile_packet - z**2 * (27 * z - 160) * (27 * z + 160) / 729)
        == 0,
        "collision hostile packet polynomial changed",
    )

    positive = {p: -7, q: 6, r: 0}
    positive_covariant = collision_covariant.subs(positive)
    positive_discriminant = quartic_discriminant.subs(positive)
    positive_packet_discriminant = packet_discriminant.subs(positive)
    require(
        positive_covariant == -11289159
        and positive_discriminant != 0
        and positive_packet_discriminant != 0,
        "noncollision positive control failed",
    )

    print("THM3273 quartic centered-norm collision covariant exact companion")
    print("packet_degree=4 pair_quotient_checks=6")
    print("disc_packet=(2^12/3^36)*C^2*disc_quartic exact=True")
    print("C=9p^6-1900p^4r+1000p^3q^2+110000p^2r^2-1000000r^3")
    print("six_pair_quotient_product=-(2^6/3^18)*C")
    print("affine_weight_C=12 packet_discriminant_weight=36")
    print("special_fibre_C=27d^12/8")
    print("collision_hostile: p=-10 q=0 r=9 C=0 disc_f_nonzero=True")
    print("positive_control: p=-7 q=6 r=0 C=-11289159")


if __name__ == "__main__":
    main()
