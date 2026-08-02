#!/usr/bin/env python3
"""Exact symbolic/finite companion for THM-3272."""

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


def main():
    x, z, a, m, t = sp.symbols("x z a m t")
    d = a - m

    # Special fibre with roots (a,m,m,m).  The fixed-sheet norm is zero;
    # every moving-sheet norm is 2(a-m)^3/27.
    fixed_complement = (m, m, m)
    fixed_mean = sum(fixed_complement) / 3
    fixed_norm = sp.expand(
        sp.prod(root - fixed_mean for root in fixed_complement)
    )
    moving_complement = (a, m, m)
    moving_mean = sum(moving_complement) / 3
    moving_norm = sp.factor(
        sp.prod(root - moving_mean for root in moving_complement)
    )
    residue_value = sp.Rational(2, 27) * d**3
    require(fixed_norm == 0, "fixed special-fibre norm is not zero")
    require(sp.expand(moving_norm - residue_value) == 0,
            "moving special-fibre norm has wrong residue")
    delta_residue = sp.factor((-residue_value) ** 3)
    require(delta_residue == -sp.Rational(8, 19683) * d**9,
            "spectral denominator residue has wrong sign or exponent")

    # Universal tame Kummer control f=(X-a)((X-m)^3-t).  The centered packet
    # characteristic polynomial specializes exactly to the asserted split
    # residue, and its singleton value is t.
    cubic = (x - m) ** 3 - t
    quartic = sp.expand((x - a) * cubic)
    _, a3, a2, a1, _ = sp.Poly(quartic, x).all_coeffs()
    packet = centered_norm_polynomial(x, a3, a2, a1)
    characteristic = sp.Poly(
        sp.resultant(quartic, z - packet, x), z
    ).as_expr()
    singleton_value = sp.factor(packet.subs(x, a))
    require(singleton_value == t,
            "Kummer control singleton centered norm should equal t")
    special_characteristic = sp.factor(characteristic.subs(t, 0))
    require(
        sp.expand(special_characteristic - z * (z - residue_value) ** 3) == 0,
        "packet characteristic polynomial has wrong special fibre",
    )
    moving_characteristic = sp.cancel(characteristic / (z - singleton_value))
    spectral_denominator = sp.factor(
        moving_characteristic.subs(z, singleton_value)
    )
    require(
        sp.expand(spectral_denominator.subs(t, 0) - delta_residue) == 0,
        "Kummer spectral denominator has wrong residue",
    )

    # A concrete rational-function projector control (a,m)=(1,0).  The
    # packet spectral projector reduces to the ordinary factor idempotent.
    control_quartic = sp.expand((x - 1) * (x**3 - t))
    _, c3, c2, c1, _ = sp.Poly(control_quartic, x).all_coeffs()
    control_packet = centered_norm_polynomial(x, c3, c2, c1)
    control_characteristic = sp.Poly(
        sp.resultant(control_quartic, z - control_packet, x), z
    ).as_expr()
    control_singleton = sp.factor(control_packet.subs(x, 1))
    control_h = sp.cancel(control_characteristic / (z - control_singleton))
    control_delta = sp.factor(control_h.subs(z, control_singleton))
    control_projector = sp.cancel(
        control_h.subs(z, control_packet) / control_delta
    )
    field = sp.QQ.frac_field(t)
    projector_remainder = sp.factor(
        sp.rem(
            sp.Poly(control_projector, x, domain=field),
            sp.Poly(control_quartic, x, domain=field),
        ).as_expr()
    )
    expected_projector = sp.cancel((x**3 - t) / (1 - t))
    require(sp.expand(projector_remainder - expected_projector) == 0,
            "packet projector differs from the factor idempotent")
    require(sp.factor(control_delta.subs(t, 0)) == sp.Rational(-8, 19683),
            "concrete projector denominator has wrong unit residue")

    # Exact finite-residue controls.  With 2 and 3 invertible, every
    # separated pair has a nonzero packet gap.  In characteristic two all
    # such gaps collapse, while a=m collapses at every prime.
    good_controls = 0
    for prime in (5, 7, 11):
        inv27 = pow(27, -1, prime)
        for av in range(prime):
            for mv in range(prime):
                if av == mv:
                    continue
                c = (2 * pow(av - mv, 3, prime) * inv27) % prime
                delta = pow(-c, 3, prime)
                require(c != 0 and delta != 0,
                        "good residue control lost the unit gap")
                good_controls += 1
    require(good_controls == 172, "good residue control count changed")

    char2_collisions = 0
    for av, mv in ((0, 1), (1, 0)):
        c = (2 * pow(av - mv, 3, 2) * pow(27, -1, 2)) % 2
        require(c == 0, "characteristic-two hostile did not collide")
        char2_collisions += 1

    cross_resultant_collisions = 0
    for prime in (5, 7, 11):
        inv27 = pow(27, -1, prime)
        for av in range(prime):
            c = (2 * pow(av - av, 3, prime) * inv27) % prime
            require(c == 0, "nonseparated hostile did not collide")
            cross_resultant_collisions += 1
    require(cross_resultant_collisions == 23,
            "cross-resultant hostile count changed")

    print("THM3272 tame C3 centered-norm integral projector exact companion")
    print("special_fibre: fixed_norm=0 moving_norm=2(a-m)^3/27")
    print("delta_residue=-8(a-m)^9/19683 unit_when_2,3,d_are_units=True")
    print("kummer_family: singleton_norm=t Pbar=Z(Z-2d^3/27)^3")
    print("spectral_projector_equals_factor_idempotent=True")
    print(f"good_residue_controls={good_controls}")
    print(f"char2_collisions={char2_collisions} cross_resultant_collisions={cross_resultant_collisions}")
    print("char3_boundary=center_not_integral")


if __name__ == "__main__":
    main()
