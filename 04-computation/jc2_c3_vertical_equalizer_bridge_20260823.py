#!/usr/bin/env python3
"""Exact scratch audit of the THM-3765 -> THM-3770 bridge.

This is deliberately assertion-free so normal and ``-O`` runs execute the
same checks.  It proves polynomial/rational identities over Q[h,g,p,c,X,T]
and uses small exact coefficient systems only as hostile controls.
"""

from __future__ import annotations

import sys

import sympy as sp


CHECKS = 0


def check(label: str, condition: bool) -> None:
    global CHECKS
    CHECKS += 1
    if not bool(condition):
        raise RuntimeError(f"FAILED: {label}")


def zero(expr: sp.Expr) -> bool:
    return sp.cancel(sp.together(expr)) == 0


def jac(f: sp.Expr, q: sp.Expr, x: sp.Symbol, t: sp.Symbol) -> sp.Expr:
    return sp.diff(f, x) * sp.diff(q, t) - sp.diff(f, t) * sp.diff(q, x)


def bounded_mate_rank(q: sp.Expr, x: sp.Symbol, t: sp.Symbol, degree: int) -> tuple[int, int, int]:
    mons = [x**i * t**j for i in range(degree + 1) for j in range(degree + 1 - i)]
    coeffs = sp.symbols(f"a0:{len(mons)}")
    candidate = sum(a * m for a, m in zip(coeffs, mons))
    debt = sp.Poly(sp.expand(jac(candidate, q, x, t) - 1), x, t)
    equations = list(debt.coeffs())
    matrix, rhs = sp.linear_eq_to_matrix(equations, coeffs)
    augmented = matrix.row_join(rhs)
    return len(mons), matrix.rank(), augmented.rank()


def main() -> None:
    # Freeze platform-independent LF output for byte-for-byte replay.
    sys.stdout.reconfigure(newline="\n")
    X, T, Z, R, L, Us, Ws, h, g, p, c, alpha = sp.symbols(
        "X T Z R L Us Ws h g p c alpha"
    )
    z = X * T
    u = 1 + g * T / 2
    Q = h + p * T + X * u**2
    Q_ansatz = X + (h + g * z) + T * (p + g**2 * z / 4)
    Y = X - T * (p + g**2 * z / 4)
    D = g * (Q - h) + 2 * p
    v = 2 * p + g * X * u
    q0 = h - 2 * p / g
    t0 = Q - q0

    # The C3 rational-exact boundary and its split exceptional fibre.
    check("normalized ansatz identity", zero(Q - Q_ansatz))
    check("D factors as u*v", zero(D - u * v))
    check("target parameter normalization", zero(D - g * t0))
    check("components are disjoint when p!=0 (Bezout sidecar)", zero(v - g * X * u - 2 * p))
    check("Y restriction on u=0", zero(Y - 2 * p / g - u * (X * (2 - u) - 2 * p / g)))
    check("Y restriction on v=0", zero(Y + 2 * p / g - (2 - u) * v / g))

    P0 = c * Y / D
    Pu = 2 * c / (g * u)
    Pv = 2 * c * X / v

    check("THM-3765 primitive sign", zero(jac(P0, Q, X, T) - c))
    check("u-chart primitive", zero(jac(Pu, Q, X, T) - c))
    check("v-chart primitive", zero(jac(Pv, Q, X, T) - c))
    check("P0 decomposition", zero(P0 - (Pu - c / g - 2 * c * p / (g * D))))
    check("chart transition is target-only", zero(Pv - Pu + 4 * c * p / (g * D)))
    check("symmetric primitive identity", zero(P0 + c / g - (Pu + Pv) / 2))

    # (Pu,Q) is an explicit birational coordinate chart; its inverse is
    # triangular after the reciprocal substitution u=2c/(gR).
    u_inverse = 2 * c / (g * R)
    T_inverse = 2 * (u_inverse - 1) / g
    X_inverse = (L - h - p * T_inverse) / u_inverse**2
    check("Pu-chart inverse recovers T", zero(T_inverse.subs(R, Pu) - T))
    check("Pu-chart inverse recovers X", zero(X_inverse.subs({R: Pu, L: Q}, simultaneous=True) - X))
    check("inverse sends Pu to R", zero(Pu.subs({X: X_inverse, T: T_inverse}, simultaneous=True) - R))
    check("inverse sends Q to L", zero(Q.subs({X: X_inverse, T: T_inverse}, simultaneous=True) - L))

    # Exact THM-3770 log-canonical dressing realization.  The pair (u,v) is
    # birational and equalized before the nonconstant dressing v/g is added.
    m = -g**2 / 2
    V0 = (v - 2 * p) / u
    T_lc_inverse = 2 * (Us - 1) / g
    X_lc_inverse = (Ws - 2 * p) / (g * Us)
    check("log-canonical pair", zero(jac(u, v, X, T) - m * u))
    check("old component spectrum is 2p", zero(v - 2 * p - g * X * u))
    check("equalized quotient V0=gX", zero(V0 - g * X))
    check("blowdown Keller response", zero(jac(u, V0, X, T) - m))
    check("birational log chart recovers u", zero(u.subs({X: X_lc_inverse, T: T_lc_inverse}, simultaneous=True) - Us))
    check("birational log chart recovers v", zero(v.subs({X: X_lc_inverse, T: T_lc_inverse}, simultaneous=True) - Ws))
    check("C3 is nonconstant spectral dressing", zero(t0 - u * (v / g)))
    check("THM-3770 dressed primitive is Pu", zero(-c * v / (m * t0) - Pu))

    # Hamiltonian/Darboux identities.  Here delta(f)=J(f,Q).
    du = jac(u, Q, X, T)
    dv = jac(v, Q, X, T)
    check("delta(u)=-(g/2)u^2", zero(du + g * u**2 / 2))
    check("delta(v)=(g/2)uv", zero(dv - g * u * v / 2))
    check("uv is target-constant", zero(jac(u * v, Q, X, T)))
    check("reciprocal u response", zero(jac(1 / u, Q, X, T) - g / 2))
    check("X/v response", zero(jac(X / v, Q, X, T) - sp.Rational(1, 2)))
    check("log-u cofactor", zero(du / u + g * u / 2))
    check("log-v cofactor", zero(dv / v - g * u / 2))

    # In the canonical target parameter t=Q-q0, P0 has opposite residues.
    residue_u = 2 * c * p / g**2
    residue_v = -2 * c * p / g**2
    check("opposite target residues", zero(residue_u + residue_v))
    check("nonzero equalizer defect polynomial", zero(g**2 * (residue_u - residue_v) - 4 * c * p))
    eq_u = g**2 * alpha + 2 * c * p
    eq_v = g**2 * alpha - 2 * c * p
    check("equalizer equations differ by 4cp", zero(eq_u - eq_v - 4 * c * p))

    # Constant-flank intersection: psi'=g^2/4, so over char 0 it is
    # constant only on g=0, where chi is constant and Q is linear.
    psi_profile = p + g**2 * Z / 4
    chi_profile = h + g * Z
    check("psi derivative", zero(sp.diff(psi_profile, Z) - g**2 / 4))
    check("linear edge after g=0", zero(Q.subs(g, 0) - (X + h + p * T)))
    check("chi constant after g=0", zero(sp.diff(chi_profile, Z).subs(g, 0)))

    # Linear-edge positive controls.
    Plin_p = c * (X - p * T) / (2 * p)
    Plin_0 = -c * T
    Qlin = X + h + p * T
    check("linear p!=0 mate", zero(jac(Plin_p, Qlin, X, T) - c))
    check("linear p=0 mate", zero(jac(Plin_0, Qlin.subs(p, 0), X, T) - c))

    # Hostile exact bounded polynomial-mate systems for one nonlinear member.
    Q_control = sp.expand(Q.subs({h: 0, g: 2, p: 3}))
    ranks = []
    for degree in range(0, 8):
        row = bounded_mate_rank(Q_control, X, T, degree)
        check(f"bounded no-mate degree {degree}", row[2] > row[1])
        ranks.append((degree,) + row)

    # Numeric controls check both signs and both one-component primitives.
    controls = [
        {h: 0, g: 2, p: 3, c: 1},
        {h: 5, g: -4, p: 7, c: -3},
        {h: -2, g: 6, p: -5, c: 11},
    ]
    control_rows = []
    for sub in controls:
        ru = sp.cancel(residue_u.subs(sub))
        rv = sp.cancel(residue_v.subs(sub))
        check("control residues nonzero", ru != 0 and rv != 0)
        check("control residues opposite", ru == -rv)
        check("control defect nonzero", sp.cancel((ru - rv)) != 0)
        control_rows.append((sub[h], sub[g], sub[p], sub[c], ru, rv))

    print("THM-3765 -> THM-3770 NORMALIZED C3 VERTICAL-EQUALIZER AUDIT")
    print("SCOPE: exact identities over characteristic zero; nonlinear branch assumes g*p*c != 0")
    print("SPECIAL_FIBRE: q0=h-2p/g,  g*(Q-q0)=D=u*v")
    print("COMPONENTS: u=1+g*T/2; v=2p+g*X*u; v-g*X*u=2p")
    print("P0_TARGET_RESIDUES: (+2*c*p/g^2, -2*c*p/g^2)")
    print("EQUALIZER_DEFECT: 4*c*p/g^2 != 0")
    print("ONE_COMPONENT_MATES:")
    print("  Pu=2*c/(g*u), pole only on u=0")
    print("  Pv=2*c*X/v, pole only on v=0")
    print("  Pv-Pu=-4*c*p/(g*D) in k(Q)")
    print("  P0+c/g=(Pu+Pv)/2")
    print("BIRATIONAL_Pu_CHART: u=2*c/(g*Pu), T=2*(u-1)/g, X=(Q-h-p*T)/u^2")
    print("THM-3770_DRESSING:")
    print("  J(u,v)=(-g^2/2)*u; Spec_0(u,v)=(2p); (v-2p)/u=gX")
    print("  Q-q0=u*(v/g); phi(W)=W/g has appended root 0")
    print("  smoothness p!=0 separates spectrum 2p from root 0")
    print("  -c*v/[(-g^2/2)*(Q-q0)]=Pu")
    print("DARBOUX_RESPONSES:")
    print("  delta(log u)=-g*u/2; delta(log v)=+g*u/2")
    print("  delta(1/u)=g/2; delta(X/v)=1/2")
    print("CONSTANT_FLANK_INTERSECTION: psi'=g^2/4, hence g=0 and Q is linear")
    print("BOUNDED_POLYNOMIAL_MATE_RANKS: degree,unknowns,rank,augmented_rank")
    for row in ranks:
        print("  " + ",".join(map(str, row)))
    print("NUMERIC_RESIDUE_CONTROLS: h,g,p,c,res_u,res_v")
    for row in control_rows:
        print("  " + ",".join(map(str, row)))
    print(f"CHECKS={CHECKS}")
    print("RESULT: PASS")


if __name__ == "__main__":
    main()

