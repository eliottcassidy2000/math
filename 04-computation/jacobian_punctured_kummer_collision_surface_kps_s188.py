#!/usr/bin/env python3
"""Exact audit of the punctured Kummer surface inside the THM-1300 map.

The distinguished curved collision component is C=0, where
C=2-3*x*y-x**2*z and F3=x*C.  Its coordinate ring is a Laurent polynomial
ring.  This script verifies that the restricted two-dimensional map is,
after explicit source and target coordinate changes, just identity times the
degree-two Kummer cover s |-> 4*s**2.
"""

from __future__ import annotations

import sympy as sp


def main() -> None:
    x, y, z = sp.symbols("x y z")
    s, v, b = sp.symbols("s v b", nonzero=True)

    u = 1 + x * y
    F1 = sp.expand(u**3 * z + y**2 * u * (4 + 3 * x * y))
    F2 = sp.expand(y + 3 * x * u**2 * z
                   + 3 * x * y**2 * (4 + 3 * x * y))
    F3 = sp.expand(2 * x - 3 * x**2 * y - x**3 * z)
    C = sp.expand(2 - 3 * x * y - x**2 * z)

    ambient_jac = sp.factor(sp.det(sp.Matrix([
        [sp.diff(F1, q) for q in (x, y, z)],
        [sp.diff(F2, q) for q in (x, y, z)],
        [sp.diff(F3, q) for q in (x, y, z)],
    ])))
    assert ambient_jac == -2
    assert sp.expand(F3 - x * C) == 0
    assert sp.expand(x * (3 * y + x * z) - (2 - C)) == 0

    # C=0 is parametrized by x=s^{-1}, y=v*s, z=s^2(2-3v).
    # Conversely s=x^{-1} is regular on C=0 because
    # x(3y+xz)=2 there, so x is a unit in characteristic zero.
    surface_sub = {
        x: 1 / s,
        y: v * s,
        z: s**2 * (2 - 3 * v),
    }
    assert sp.cancel(C.subs(surface_sub)) == 0
    alpha = sp.factor(sp.cancel(F1.subs(surface_sub)))
    beta = sp.factor(sp.cancel(F2.subs(surface_sub)))
    gamma = sp.factor(sp.cancel(F3.subs(surface_sub)))
    expected_alpha = s**2 * (v + 1) * (v + 2)
    expected_beta = 2 * s * (2 * v + 3)
    assert sp.expand(alpha - expected_alpha) == 0
    assert sp.expand(beta - expected_beta) == 0
    assert gamma == 0

    jac_sv = sp.factor(sp.det(sp.Matrix([
        [sp.diff(alpha, s), sp.diff(alpha, v)],
        [sp.diff(beta, s), sp.diff(beta, v)],
    ])))
    assert sp.expand(jac_sv + 2 * s**2) == 0

    surface_xy = {z: (2 - 3 * x * y) / x**2}
    alpha_xy = sp.factor(sp.cancel(F1.subs(surface_xy)))
    beta_xy = sp.factor(sp.cancel(F2.subs(surface_xy)))
    jac_xy = sp.factor(sp.cancel(
        sp.diff(alpha_xy, x) * sp.diff(beta_xy, y)
        - sp.diff(alpha_xy, y) * sp.diff(beta_xy, x)
    ))
    assert sp.cancel(jac_xy - 2 / x**3) == 0

    # Source Laurent automorphism: (s,v) -> (s,beta).
    source_b = beta
    source_change_jac = sp.factor(sp.det(sp.Matrix([
        [sp.diff(s, s), sp.diff(s, v)],
        [sp.diff(source_b, s), sp.diff(source_b, v)],
    ])))
    assert sp.expand(source_change_jac - 4 * s) == 0
    recovered_v = sp.cancel((b - 6 * s) / (4 * s))
    assert sp.cancel(source_b.subs(v, recovered_v) - b) == 0

    # Target polynomial automorphism: (alpha,beta) ->
    # (b=beta, delta=beta^2-16*alpha).
    aa, bb = sp.symbols("aa bb")
    target_delta = bb**2 - 16 * aa
    target_change_jac = sp.factor(sp.det(sp.Matrix([
        [sp.diff(bb, aa), sp.diff(bb, bb)],
        [sp.diff(target_delta, aa), sp.diff(target_delta, bb)],
    ])))
    assert target_change_jac == 16
    delta = sp.factor(beta**2 - 16 * alpha)
    assert sp.expand(delta - 4 * s**2) == 0
    assert sp.expand(alpha - (beta**2 - delta) / 16) == 0

    # In the new source coordinates (s,b), the map is (b,4s^2).
    normal_b = b
    normal_delta = 4 * s**2
    normal_jac = sp.factor(sp.det(sp.Matrix([
        [sp.diff(normal_b, s), sp.diff(normal_b, b)],
        [sp.diff(normal_delta, s), sp.diff(normal_delta, b)],
    ])))
    assert normal_jac == -8 * s

    # Deck involution in the original Laurent coordinates.
    deck = {s: -s, v: -v - 3}
    assert sp.expand(alpha.subs(deck, simultaneous=True) - alpha) == 0
    assert sp.expand(beta.subs(deck, simultaneous=True) - beta) == 0
    assert sp.expand((-s).subs(deck, simultaneous=True) - s) == 0

    p_plus = {s: 1, v: -sp.Rational(3, 2)}
    p_minus = {s: -1, v: -sp.Rational(3, 2)}
    image_plus = (sp.factor(alpha.subs(p_plus)),
                  sp.factor(beta.subs(p_plus)))
    image_minus = (sp.factor(alpha.subs(p_minus)),
                   sp.factor(beta.subs(p_minus)))
    assert image_plus == image_minus == (-sp.Rational(1, 4), 0)

    # The natural affine closure retains the formulas but ramifies at s=0.
    closure_jac = normal_jac
    assert closure_jac.subs(s, 0) == 0

    print("THM-1300 COLLISION SURFACE -- PUNCTURED KUMMER AUDIT")
    print(f"ambient det JF = {ambient_jac}; F3 = x*C with C = {C}")
    print("on C=0: x*(3*y+x*z)=2, hence x is a unit")
    print("Laurent coordinates: x=1/s, y=v*s, z=s^2*(2-3*v)")
    print(f"restricted (F1,F2) = ({alpha}, {beta})")
    print(f"det d(F1,F2)/d(s,v) = {jac_sv}")
    print(f"det d(F1,F2)/d(x,y) on C=0 = {jac_xy}")
    print(f"source change (s,v)->(s,b=F2): Jacobian = {source_change_jac}")
    print(f"target change (F1,F2)->(b,delta): Jacobian = {target_change_jac}")
    print(f"delta=F2^2-16*F1 = {delta}")
    print(f"normal form (s,b)->(b,4*s^2): Jacobian = {normal_jac}")
    print("deck involution: (s,v)->(-s,-v-3)")
    print(f"collision images at (s,v)=(+/-1,-3/2): {image_plus}")
    print("the natural A^2 closure ramifies on s=0; on s!=0 the map is etale")
    print("VERDICT: the curved collision component is exactly a degree-two punctured Kummer cover.")
    print("SCOPE: this is an etale near-counterexample on G_m x A^1, not a JC(2) counterexample.")


if __name__ == "__main__":
    main()
