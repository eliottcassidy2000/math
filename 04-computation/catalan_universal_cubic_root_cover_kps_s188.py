#!/usr/bin/env python3
"""Exact audit: THM-3545's Catalan branch is the universal cubic root cover.

Adjoin r=sqrt(1-3*kappa*w).  The formerly nonpolynomial Catalan map becomes
polynomial in (v,r), and affine source/target changes identify it with the
universal marked-root map (t,p)->(p,-t^3-p*t).  The script also verifies the
discriminant pullback, the three-point collision, and a first-jet obstruction
to correcting the ramification while fixing the ramification line.
"""

from __future__ import annotations

import sympy as sp


def main() -> None:
    v, r = sp.symbols("v r")
    kappa = sp.symbols("kappa", nonzero=True)

    # THM-3545 after adjoining r=sqrt(1-3*kappa*w).
    A = sp.Rational(2, 3) * (1 - r)
    B = 1 - r
    P = sp.expand(v**2 + A)
    Q = sp.expand(v**3 - v + v * B)
    jac_vr = sp.factor(sp.det(sp.Matrix([
        [sp.diff(P, v), sp.diff(P, r)],
        [sp.diff(Q, v), sp.diff(Q, r)],
    ])))
    assert sp.expand(jac_vr + sp.Rational(2, 3) * r) == 0

    # Since r^2=1-3*kappa*w, dr/dw=-3*kappa/(2r).  This is the
    # inverse-ramification factor that turns jac_vr into kappa.
    dr_dw = -3 * kappa / (2 * r)
    jac_vw = sp.cancel(jac_vr * dr_dw)
    assert jac_vw == kappa

    # Affine source and target changes.
    t = v
    p = sp.expand(2 * r - 3 * v**2)
    q = sp.expand(2 * Q)
    source_change_jac = sp.factor(sp.det(sp.Matrix([
        [sp.diff(t, v), sp.diff(t, r)],
        [sp.diff(p, v), sp.diff(p, r)],
    ])))
    assert source_change_jac == 2

    PP, QQ = sp.symbols("PP QQ")
    target_p = 2 - 3 * PP
    target_q = 2 * QQ
    target_change_jac = sp.factor(sp.det(sp.Matrix([
        [sp.diff(target_p, PP), sp.diff(target_p, QQ)],
        [sp.diff(target_q, PP), sp.diff(target_q, QQ)],
    ])))
    assert target_change_jac == -6
    assert sp.expand(p - (2 - 3 * P)) == 0
    assert sp.expand(q + t**3 + p * t) == 0

    # Universal marked-root cover and its branch data.
    tt, pp = sp.symbols("tt pp")
    universal_q = -tt**3 - pp * tt
    universal_jac = sp.factor(sp.det(sp.Matrix([
        [sp.diff(pp, tt), sp.diff(pp, pp)],
        [sp.diff(universal_q, tt), sp.diff(universal_q, pp)],
    ])))
    assert universal_jac == pp + 3 * tt**2
    assert sp.expand(universal_jac.subs({tt: t, pp: p}) - 2 * r) == 0

    X = sp.symbols("X")
    cubic = X**3 + pp * X + sp.symbols("qq")
    cubic_disc = sp.factor(sp.discriminant(cubic, X))
    assert cubic_disc == -4 * pp**3 - 27 * sp.symbols("qq")**2
    pulled_disc = sp.factor(cubic_disc.subs({pp: p, sp.symbols("qq"): q}))
    assert sp.expand(pulled_disc + 4 * r**2 * (8 * r - 9 * v**2)) == 0

    target_cusp = sp.factor((3 * P - 2)**3 - 27 * Q**2)
    assert sp.expand(4 * target_cusp - pulled_disc) == 0
    branch_P = sp.factor(P.subs(r, 0))
    branch_Q = sp.factor(Q.subs(r, 0))
    assert sp.expand(branch_P - (v**2 + sp.Rational(2, 3))) == 0
    assert branch_Q == v**3
    assert sp.expand((3 * branch_P - 2)**3 - 27 * branch_Q**2) == 0

    # The Catalan collision pair acquires the third marked root.
    triple = [
        (sp.Integer(0), -sp.Rational(1, 2)),
        (sp.Integer(1), sp.Integer(1)),
        (-sp.Integer(1), sp.Integer(1)),
    ]
    triple_images = [(sp.factor(P.subs({v: vv, r: rr})),
                      sp.factor(Q.subs({v: vv, r: rr})))
                     for vv, rr in triple]
    triple_jacs = [sp.factor(jac_vr.subs({v: vv, r: rr}))
                   for vv, rr in triple]
    assert triple_images == [(1, 0), (1, 0), (1, 0)]
    assert triple_jacs == [sp.Rational(1, 3),
                           -sp.Rational(2, 3),
                           -sp.Rational(2, 3)]

    # First-jet ramification-line gate.  Let R=p+3t^2.  If both corrected
    # components are congruent to the universal cover modulo R, write them
    # as p+R*a and q+R*b.  Their Jacobian on R=0 is forced to vanish at t=0.
    T, R, a0, b0 = sp.symbols("T R a0 b0")
    # Rows of the corrected differential after setting R=0.  Derivatives of
    # a,b disappear because they are multiplied by R.
    row_P = (6 * T * a0, 1 + a0)
    row_Q = (6 * T * b0, -T + b0)
    corrected_jac_on_R0 = sp.factor(
        row_P[0] * row_Q[1] - row_P[1] * row_Q[0]
    )
    assert sp.expand(corrected_jac_on_R0 + 6 * T * (T * a0 + b0)) == 0
    assert corrected_jac_on_R0.subs(T, 0) == 0

    print("CATALAN THICKENING -- UNIVERSAL CUBIC ROOT-COVER AUDIT")
    print(f"polynomialized map (P,Q) = ({P}, {Q})")
    print(f"det d(P,Q)/d(v,r) = {jac_vr}")
    print(f"dr/dw = {dr_dw}; det d(P,Q)/d(v,w) = {jac_vw}")
    print(f"source affine change (v,r)->(t,p): p={p}, Jacobian={source_change_jac}")
    print(f"target affine change (P,Q)->(p,q): Jacobian={target_change_jac}")
    print(f"universal normal form: (t,p)->(p,{universal_q})")
    print(f"universal Jacobian = {universal_jac}; pulled back = 2*r")
    print(f"cubic discriminant = {cubic_disc}")
    print(f"discriminant pullback = {pulled_disc}")
    print(f"ramification line r=0 maps as (P,Q)=({branch_P},{branch_Q})")
    print(f"triple fiber over (1,0): {triple}; Jacobians={triple_jacs}")
    print(f"if corrections fix R=0, corrected Jacobian there = {corrected_jac_on_R0}")
    print("VERDICT: the Catalan square root cancels the universal cubic cover's ramification factor.")
    print("SCOPE: the polynomialized map is ramified; it is a cubic-cover template, not JC(2).")


if __name__ == "__main__":
    main()
