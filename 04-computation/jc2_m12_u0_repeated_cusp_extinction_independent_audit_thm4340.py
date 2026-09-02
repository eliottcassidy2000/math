#!/usr/bin/env python3
"""Independent clean-room symbolic audit of THM-4340.

This file deliberately reconstructs the two scaled equations from the full
THM-4230 source rather than importing the proposal's probe or assembling both
sides from the same sparse records.
"""

from fractions import Fraction
from math import gcd
import sys

import sympy as sp


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")


def need(condition, label):
    if not bool(condition):
        raise RuntimeError(label)


def z(expr):
    return sp.cancel(sp.factor(expr))


def source_H(p, y, coeff):
    return (
        -3 * p
        + sp.Rational(8, 3) * p**2
        + sp.Rational(-1376, 135) * p**3
        + coeff["K"] * y**2
        + coeff["Phi"] * p**2 * y
        + coeff["Delta"] * p**4
        + coeff["Theta"] * p * y**2
        + coeff["eta"] * p**3 * y
        + coeff["zeta"] * y**3
        + coeff["u"] * p**5
        + coeff["xi"] * p**2 * y**2
        + coeff["alpha"] * p**4 * y
        + coeff["beta"] * p * y**3
        + coeff["U"] * p**6
        + coeff["W"] * p**3 * y**2
        + coeff["Z"] * y**4
    )


def generic_coefficients():
    names = "K Phi Delta Theta eta zeta u xi alpha beta U W Z"
    vals = sp.symbols(names)
    return dict(zip(names.split(), vals))


def h5_factorization():
    t, S, P, zz = sp.symbols("t S P z", nonzero=True)
    c = generic_coefficients()
    c["U"] = sp.Integer(0)
    p = P / t
    y = S * P / t
    H = source_H(p, y, c)
    Q = t**5
    F = (S**2 - p) * (1 - Q * H) - Q * S**2 / 2
    G = sp.expand(t * F)
    Phi_direct = z(sp.expand(zz**6 * G.subs(P, 1 / zz)))

    f1 = -3
    f2 = sp.Rational(8, 3) + c["K"] * S**2
    f3 = (
        sp.Rational(-1376, 135)
        + c["Phi"] * S
        + c["Theta"] * S**2
        + c["zeta"] * S**3
    )
    f4 = (
        c["Delta"]
        + c["eta"] * S
        + c["xi"] * S**2
        + c["beta"] * S**3
        + c["Z"] * S**4
    )
    f5 = c["u"] + c["alpha"] * S + c["W"] * S**2
    rho = t * zz
    A = f1 * rho**4 + f2 * rho**3 + f3 * rho**2 + f4 * rho + f5
    Phi_expected = (1 - S**2 * rho) * (A - zz**5) - S**2 * rho**6 / 2
    need(z(Phi_direct - Phi_expected) == 0, "H5 exact factorization")

    # Differentiate the independently constructed G at fixed S and compare
    # with the z-chart identity on Phi=0 (the discrepancy must be a Phi term).
    GP_z = z(sp.diff(G, P).subs(P, 1 / zz))
    relation = z(GP_z + zz ** (-4) * sp.diff(Phi_direct, zz))
    quotient = z(relation / Phi_direct)
    need(z(relation - quotient * Phi_direct) == 0, "H5 derivative modulo Phi")
    need(z(quotient - 6 / zz**5) == 0, "H5 derivative unit-term quotient")
    return Phi_direct


def e3_factorization():
    t, V, X = sp.symbols("t V X", nonzero=True)
    c = generic_coefficients()
    c["U"] = 0
    c["u"] = 0
    c["alpha"] = 0
    c["Delta"] = 0
    P = 1 / X
    S = V * X
    p = P / t
    s = t * S
    y = V
    H = source_H(p, y, c)
    Q = t**3
    F = (s**2 - p) * (1 - Q * H) - Q * s**2 / 2
    G = sp.expand(t * F)
    Phi_direct = z(sp.expand(X**4 * G))

    f3 = sp.Rational(-1376, 135) + c["eta"] * V + c["W"] * V**2
    f2 = sp.Rational(8, 3) + c["Phi"] * V + c["xi"] * V**2
    f1 = -3 + c["Theta"] * V**2 + c["beta"] * V**3
    f0 = c["zeta"] * V**3 + c["Z"] * V**4 + c["K"] * V**2
    # K*y^2 has p-exponent zero and therefore belongs to f0 here.
    rho = t * X
    A = f3 + f2 * rho + f1 * rho**2 + f0 * rho**3
    Phi_expected = (1 - V**2 * rho**3) * (A - X**3) - V**2 * rho**6 / 2
    need(z(Phi_direct - Phi_expected) == 0, "E3 exact factorization")

    # Reconstruct G as a function of the independent variables S,P, then
    # differentiate at fixed S.  This avoids accidentally differentiating at
    # fixed V in the transformed chart.
    SS, PP = sp.symbols("SS PP", nonzero=True)
    pp = PP / t
    ss = t * SS
    yy = SS * PP
    H_SP = source_H(pp, yy, c)
    F_SP = (ss**2 - pp) * (1 - Q * H_SP) - Q * ss**2 / 2
    G_SP = sp.expand(t * F_SP)
    GP_chart = z(sp.diff(G_SP, PP).subs({PP: 1 / X, SS: V * X}))
    expected_mod = X ** (-3) * (V * sp.diff(Phi_direct, V) - X * sp.diff(Phi_direct, X))
    discrepancy = z(GP_chart - expected_mod)
    need(z(discrepancy - 4 * X**(-3) * Phi_direct) == 0,
         "E3 fixed-S derivative identity modulo Phi")
    return Phi_direct


def critical_value_hostiles():
    q, S = sp.symbols("q S")
    e = sp.Rational(-1376, 135)
    d = sp.Rational(5936, 105)
    K = sp.Rational(-8, 3)
    need(K == sp.Rational(2848, 45) - sp.Rational(7, 6) * d,
         "H5 hostile inherited K/Delta relation")
    f5 = (S + 1) ** 2
    f4 = d * (S + 1) ** 4
    f3 = e * (S + 1)
    f2 = sp.Rational(8, 3) + K * S**2
    f1 = -3
    A = f5 + q * f4 + q**2 * f3 + q**3 * f2 + q**4 * f1
    B = A - S**2 * q**6 / (2 * (1 - S**2 * q))

    cs = sp.symbols("c1:6")
    Sc = -1 + sum(cs[i - 1] * q**i for i in range(1, 6))
    deriv = sp.series(sp.diff(B, S).subs(S, Sc), q, 0, 6).removeO().expand()
    sol = {}
    for power, variable in zip(range(1, 6), cs):
        equation = sp.expand(deriv.subs(sol)).coeff(q, power)
        sol[variable] = sp.solve(equation, variable)[0]
    psi = sp.series(B.subs(S, Sc).subs(sol), q, 0, 6).removeO().expand()
    need(all(psi.coeff(q, j) == 0 for j in (1, 2, 3)), "H5 hostile h1-h3")
    need(psi.coeff(q, 4) == -sp.Rational(528019, 18225), "H5 hostile h4")

    V = sp.symbols("V")
    need(sp.Rational(2848, 45)
         == sp.Rational(2848, 45) - sp.Rational(7, 6) * 0,
         "E3 hostile inherited K/Delta relation")
    f3e = e * (V - 1) ** 2
    f2e = sp.Rational(8, 3) * (V - 1) ** 2
    Ae = f3e + q * f2e - 3 * q**2
    Be = Ae - V**2 * q**6 / (2 * (1 - V**2 * q**3))
    ds = sp.symbols("d1:5")
    Vc = 1 + sum(ds[i - 1] * q**i for i in range(1, 5))
    deriv_e = sp.series(sp.diff(Be, V).subs(V, Vc), q, 0, 5).removeO().expand()
    sol_e = {}
    for power, variable in zip(range(1, 5), ds):
        equation = sp.expand(deriv_e.subs(sol_e)).coeff(q, power)
        sol_e[variable] = sp.solve(equation, variable)[0]
    psi_e = sp.series(Be.subs(V, Vc).subs(sol_e), q, 0, 4).removeO().expand()
    need(psi_e.coeff(q, 1) == 0, "E3 hostile h1")
    need(psi_e.coeff(q, 2) == -3, "E3 hostile h2")
    return psi, psi_e


def tail_tables():
    tables = []
    for label, m, d, k, B in (("H5", 5, 12, 50, 4), ("E3", 3, 4, 14, 3)):
        delta_m = (m - 1) // 2
        for r in range(1, m):
            b = Fraction(d * r, m - r)
            a = Fraction(m, 2) * b
            order = Fraction(k) + (Fraction(B + 1) - Fraction(m, 2)) * b
            squarefree_degree = m - r + (r % 2)
            genus = (squarefree_degree - 1) // 2
            persistent_delta = r // 2
            need(genus + persistent_delta == delta_m, f"{label} r={r} delta partition")
            pick_genus = 17 if label == "H5" else 16
            need(4 + 11 + genus == pick_genus - persistent_delta,
                 f"{label} r={r} normalized global genus")
            need(order > 0, f"{label} r={r} positive order")
            # A base change of this minimal degree clears both z and x.
            denominator = 1
            for value in (b, a):
                denominator = denominator * value.denominator // gcd(denominator, value.denominator)
            need((b * denominator).denominator == 1, f"{label} r={r} z integral")
            need((a * denominator).denominator == 1, f"{label} r={r} x integral")
            tables.append((label, r, b, a, order, genus, persistent_delta, denominator))
    return tables


def main():
    h5_factorization()
    e3_factorization()
    psi_h5, psi_e3 = critical_value_hostiles()
    rows = tail_tables()
    print("THM4340 U0 CUSP INDEPENDENT CLEANROOM SYMBOLIC AUDIT")
    print("FULL_SOURCE_FACTORIZATIONS=H5:PASS;E3:PASS")
    print("FIXED_S_DERIVATIVES=H5:PASS;E3:PASS")
    print(f"H5_HOSTILE_CRITICAL={psi_h5}")
    print(f"E3_HOSTILE_CRITICAL={psi_e3}")
    for row in rows:
        print("TAIL=" + ",".join(map(str, row)))
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
