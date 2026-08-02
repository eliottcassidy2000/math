#!/usr/bin/env python3
"""Exact symbolic companion for THM-3274."""

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def depressed_data(quartic, x):
    """Return the depressed coefficients and THM-3271 packet map."""
    poly = sp.Poly(sp.expand(quartic), x)
    coeffs = poly.all_coeffs()
    require(len(coeffs) == 5 and coeffs[0] == 1, "quartic must be monic")
    _, a3, a2, a1, a0 = coeffs
    require(sp.expand(a3) == 0, "quartic is not depressed")
    packet = -(20 * x**3 + 18 * a2 * x + 27 * a1) / sp.Integer(27)
    return a2, a1, a0, sp.expand(packet)


def collision_covariant(p, q, r):
    return sp.expand(
        9 * p**6
        - 1900 * p**4 * r
        + 1000 * p**3 * q**2
        + 110000 * p**2 * r**2
        - 1000000 * r**3
    )


def main():
    x, w, A = sp.symbols("x w A", nonzero=True)

    # THM-3201, lane 3 does not divide m.  After multiplying the roots by
    # pi^m, the four leading depressed roots are 0,A,A*zeta,A*zeta^2.
    # Their polynomial is X(X^3-A^3), which avoids adjoining zeta here.
    off_quartic = sp.expand(x * (x**3 - A**3))
    off_p, off_q, off_r, off_packet = depressed_data(off_quartic, x)
    off_characteristic = sp.factor(
        sp.resultant(off_quartic, w - off_packet, x)
    )
    off_expected = sp.expand(
        (w - A**3) * (w - sp.Rational(7, 27) * A**3) ** 3
    )
    require(sp.expand(off_characteristic - off_expected) == 0,
            "off-resonance scaled packet characteristic failed")
    off_fixed_value = sp.factor(off_packet.subs(x, 0))
    off_moving_value = sp.factor(off_packet.subs(x**3, A**3))
    require(off_fixed_value == A**3, "off-resonance fixed value failed")
    require(off_moving_value == sp.Rational(7, 27) * A**3,
            "off-resonance moving value failed")
    off_gap = sp.factor(off_fixed_value - off_moving_value)
    off_delta = sp.factor(off_gap**3)
    require(off_gap == sp.Rational(20, 27) * A**3,
            "off-resonance fixed/moving gap failed")
    require(
        collision_covariant(off_p, off_q, off_r) == 0,
        "off-resonance moving triple collision must lie on C=0",
    )

    # THM-3201, lane 3 divides m.  The scaled leading depressed roots are
    # -3A/4 and A/4 with the latter occurring three times.
    on_quartic = sp.expand((x + 3 * A / 4) * (x - A / 4) ** 3)
    on_p, on_q, on_r, on_packet = depressed_data(on_quartic, x)
    on_characteristic = sp.factor(
        sp.resultant(on_quartic, w - on_packet, x)
    )
    on_expected = sp.expand(
        w * (w + sp.Rational(2, 27) * A**3) ** 3
    )
    require(sp.expand(on_characteristic - on_expected) == 0,
            "resonant scaled packet characteristic failed")
    on_fixed_value = sp.factor(on_packet.subs(x, -3 * A / 4))
    on_moving_value = sp.factor(on_packet.subs(x, A / 4))
    require(on_fixed_value == 0, "resonant fixed value failed")
    require(on_moving_value == -sp.Rational(2, 27) * A**3,
            "resonant moving value failed")
    on_gap = sp.factor(on_fixed_value - on_moving_value)
    on_delta = sp.factor(on_gap**3)
    require(on_gap == sp.Rational(2, 27) * A**3,
            "resonant fixed/moving gap failed")
    on_covariant = sp.factor(collision_covariant(on_p, on_q, on_r))
    require(on_covariant == sp.Rational(27, 8) * A**12,
            "resonant collision covariant failed")

    # The cofactor is independent information even after the fixed packet
    # component has been decoded.  In L=K[y]/(y^3-t), f=(X-1)(X^3-t), take
    # c0=1/f'(1) and c=1/f'(y), so both physical Jacobians are one.  The
    # norm-one gauge (c0,c)->(lambda^-3 c0,lambda c) preserves the total
    # cofactor norm but changes the pointed physical-Jacobian ratio.
    lam, t = sp.symbols("lambda t", nonzero=True)
    d0 = 1 - t  # g(1)=1-t
    base_j0 = sp.cancel((1 / d0) * d0)
    base_jc = sp.Integer(1)
    twisted_j0 = sp.factor(lam ** -3 * base_j0)
    twisted_jc = sp.factor(lam * base_jc)
    twisted_ratio = sp.factor(twisted_jc / twisted_j0)
    norm_gauge = sp.factor(lam ** -3 * lam**3)
    require(base_j0 == base_jc == 1, "base Keller-compatible packet failed")
    require(norm_gauge == 1, "cofactor norm gauge is not preserved")
    require(twisted_ratio == lam**4, "twisted pointed ratio failed")

    lam_control = sp.Integer(2)
    ratio_control = twisted_ratio.subs(lam, lam_control)
    shifted_norm_control = (ratio_control - 1) ** 3
    require(ratio_control == 16 and shifted_norm_control == 3375,
            "lambda=2 cofactor hostile failed")

    # The root quartic and every packet scalar are unchanged by the cofactor
    # twist.  Record a nondegenerate depressed-covariant control for that
    # quartic so the hostile is not accidentally on the packet collision
    # divisor.
    y = sp.symbols("y")
    root_quartic = sp.expand((x - 1) * (x**3 - t))
    depressed = sp.Poly(sp.expand(root_quartic.subs(x, y + sp.Rational(1, 4))), y)
    _, depressed_a3, hostile_p, hostile_q, hostile_r = depressed.all_coeffs()
    require(depressed_a3 == 0, "cofactor hostile depression failed")
    hostile_covariant = sp.factor(
        collision_covariant(hostile_p, hostile_q, hostile_r)
    )
    hostile_covariant_expected = sp.factor(
        -sp.Rational(27, 64)
        * (125 * t - 1)
        * (8000 * t**2 - 475 * t + 8)
    )
    require(hostile_covariant == hostile_covariant_expected,
            "cofactor hostile packet covariant changed")
    require(hostile_covariant.subs(t, 0) == sp.Rational(27, 8),
            "cofactor hostile packet covariant is not a t-adic unit")

    print("THM3274 graph-quartic fixed decoder and cofactor hostile exact companion")
    print("off_resonance_scaled_roots_polynomial=X*(X^3-A^3)")
    print("off_resonance_packet=(W-A^3)*(W-7A^3/27)^3")
    print("off_resonance_fixed_moving_gap=20A^3/27")
    print(f"off_resonance_fixed_projector_denominator={off_delta}")
    print("off_resonance_C=0 moving_packet_triple_collision=True")
    print("resonant_scaled_roots=(-3A/4,A/4,A/4,A/4)")
    print("resonant_packet=W*(W+2A^3/27)^3")
    print("resonant_fixed_moving_gap=2A^3/27")
    print(f"resonant_fixed_projector_denominator={on_delta}")
    print("resonant_C=27A^12/8")
    print("cofactor_norm_gauge=(lambda^-3,lambda,lambda,lambda) product=1")
    print("twisted_pointed_Jacobian_ratio=lambda^4")
    print("lambda=2 ratio=16 shifted_norm=(16-1)^3=3375")
    print("cofactor_hostile_field=k((t)), y^3=t, k contains zeta_3")
    print("cofactor_hostile_C=-27(125t-1)(8000t^2-475t+8)/64 residue=27/8")


if __name__ == "__main__":
    main()
