#!/usr/bin/env python3
"""Exact audit of the torus quotient of the fixed THM-1300 Keller map.

The calculation is deliberately small.  It verifies the induced polynomial
map on the two categorical quotient planes, its Jacobian square, the two
quotient collision points, and the quadratic initial map at the ramified
collision point.  SymPy is used only as an exact polynomial arithmetic engine.
"""

from __future__ import annotations

import sympy as sp


def main() -> None:
    v, t = sp.symbols("v t")
    u = 1 + v

    # In the Laurent chart x != 0 the fixed map has
    # x^2 F1=A, x F2=B, F3=x C.  All quotient identities obtained below are
    # polynomial identities, so the chart derivation extends across x=0.
    A = sp.expand(t * u**3 + v**2 * u * (4 + 3 * v))
    B = sp.expand(v + 3 * t * u**2 + 3 * v**2 * (4 + 3 * v))
    C = sp.expand(2 - 3 * v - t)

    V = sp.expand(B * C)
    T = sp.expand(A * C**2)
    jac = sp.factor(sp.det(sp.Matrix([[sp.diff(V, v), sp.diff(V, t)],
                                      [sp.diff(T, v), sp.diff(T, t)]])))
    assert sp.expand(jac - 2 * C**2) == 0

    q0 = {v: sp.Integer(0), t: sp.Integer(0)}
    q1 = {v: -sp.Rational(3, 2), t: sp.Rational(13, 2)}
    q0_data = tuple(sp.factor(f.subs(q0)) for f in (A, B, C, V, T, jac))
    q1_data = tuple(sp.factor(f.subs(q1)) for f in (A, B, C, V, T, jac))
    assert q0_data == (0, 0, 2, 0, 0, 8)
    assert q1_data == (-sp.Rational(1, 4), 0, 0, 0, 0, 0)
    assert sp.factor(B.subs(t, 2 - 3 * v) - (4 * v + 6)) == 0
    assert sp.factor(A.subs(t, 2 - 3 * v) - (v + 1) * (v + 2)) == 0
    assert sp.expand(V.subs(t, 2 - 3 * v)) == 0
    assert sp.expand(T.subs(t, 2 - 3 * v)) == 0

    # Polynomiality is load-bearing.  The rational target recombination
    # (T/V^2,V) cancels C^2 and leaves a B-pole instead.
    rational_R = A / B**2
    rational_S = V
    rational_jac = sp.factor(
        sp.diff(rational_R, v) * sp.diff(rational_S, t)
        - sp.diff(rational_R, t) * sp.diff(rational_S, v)
    )
    assert sp.cancel(rational_jac + 2 / B**2) == 0

    # Local coordinates r=v+3/2 and c=C at q1.  The source coordinate change
    # has determinant -1, hence the local Jacobian is -2 c^2.
    r, c = sp.symbols("r c")
    local_sub = {v: r - sp.Rational(3, 2),
                 t: 2 - 3 * (r - sp.Rational(3, 2)) - c}
    A_local = sp.expand(A.subs(local_sub))
    B_local = sp.expand(B.subs(local_sub))
    V_local = sp.expand(c * B_local)
    T_local = sp.expand(c**2 * A_local)
    local_jac = sp.factor(sp.det(sp.Matrix([
        [sp.diff(V_local, r), sp.diff(V_local, c)],
        [sp.diff(T_local, r), sp.diff(T_local, c)],
    ])))
    assert sp.expand(local_jac + 2 * c**2) == 0

    def homogeneous_part(poly: sp.Expr, degree: int) -> sp.Expr:
        p = sp.Poly(poly, r, c, domain=sp.QQ)
        return sp.expand(sum(
            coeff * r**monom[0] * c**monom[1]
            for monom, coeff in p.terms()
            if sum(monom) == degree
        ))

    V2 = homogeneous_part(V_local, 2)
    T2 = homogeneous_part(T_local, 2)
    assert sp.expand(V2 - (4 * c * r - sp.Rational(3, 4) * c**2)) == 0
    assert sp.expand(T2 + sp.Rational(1, 4) * c**2) == 0
    initial_jac = sp.factor(sp.det(sp.Matrix([
        [sp.diff(V2, r), sp.diff(V2, c)],
        [sp.diff(T2, r), sp.diff(T2, c)],
    ])))
    assert sp.expand(initial_jac + 2 * c**2) == 0

    print("THM-1300 TORUS QUOTIENT -- EXACT RAMIFICATION AUDIT")
    print(f"A(v,t) = {sp.factor(A)}")
    print(f"B(v,t) = {sp.factor(B)}")
    print(f"C(v,t) = {C}")
    print(f"V(v,t) = B*C = {sp.factor(V)}")
    print(f"T(v,t) = A*C^2 = {sp.factor(T)}")
    print(f"det d(V,T)/d(v,t) = {jac}")
    print(f"q0=(0,0): (A,B,C,V,T,jac) = {q0_data}")
    print(f"q1=(-3/2,13/2): (A,B,C,V,T,jac) = {q1_data}")
    print(f"on C=0: B={sp.factor(B.subs(t, 2 - 3*v))}, "
          f"A={sp.factor(A.subs(t, 2 - 3*v))}")
    print("the whole line C=0 maps to (V,T)=(0,0)")
    print(f"rational escape (T/V^2,V): Jacobian = {rational_jac}")
    print(f"local A(r,c) = {A_local}")
    print(f"local B(r,c) = {B_local}")
    print(f"quadratic initial map = ({V2}, {T2})")
    print(f"local Jacobian = {local_jac}; initial Jacobian = {initial_jac}")
    print("VERDICT: the categorical quotient has a forced ramification square C^2.")
    print("SCOPE: this obstructs invariant polynomial compression only; it is not JC(2).")


if __name__ == "__main__":
    main()
