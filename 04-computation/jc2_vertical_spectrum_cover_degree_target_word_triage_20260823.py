#!/usr/bin/env python3
"""Exact scratch probe for the THM-3770/71/72/73/74/75/76 junction.

The main calculation independently checks and four-parameter-scales the
``m=1`` row of THM-3774, including its quadratic target-field degree.  A
separate birational THM-3771 control shows why unequal vertical residues
alone cannot obstruct unrestricted rational target changes.  The script
contains no Python asserts so normal and optimized runs execute identical
checks.
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


def jac(f: sp.Expr, q: sp.Expr, x: sp.Symbol, y: sp.Symbol) -> sp.Expr:
    return sp.diff(f, x) * sp.diff(q, y) - sp.diff(f, y) * sp.diff(q, x)


def bounded_mate_rank(q: sp.Expr, x: sp.Symbol, y: sp.Symbol, degree: int) -> tuple[int, int, int]:
    mons = [x**i * y**j for i in range(degree + 1) for j in range(degree + 1 - i)]
    coeffs = sp.symbols(f"m0:{len(mons)}")
    candidate = sum(a * mon for a, mon in zip(coeffs, mons))
    debt = sp.Poly(sp.expand(jac(candidate, q, x, y) - 1), x, y)
    matrix, rhs = sp.linear_eq_to_matrix(list(debt.coeffs()), coeffs)
    return len(mons), matrix.rank(), matrix.row_join(rhs).rank()


def main() -> None:
    sys.stdout.reconfigure(newline="\n")

    x, y, L, R, Z, a, b, c, d = sp.symbols("x y L R Z a b c d")
    A = a + b * x * y
    B = c + d * x**3 * A
    U = x * A * B
    P = (2 * B - c) / x
    W = sp.cancel(U * P)
    kappa = -b * c**2

    # Uniform rational Keller/log-canonical identities.
    check("nonradial constant Jacobian", zero(jac(P, U, x, y) - kappa))
    check("blowdown polynomial", zero(W - A * B * (2 * B - c)))
    check("log-canonical blowdown", zero(jac(U, W, x, y) - b * c**2 * U))

    # Polynomialized Bezout identity: it forces a critical point onto x=0;
    # the remaining axis derivative is a*c.
    N = 2 * B - c
    bezout = (x * sp.diff(N, x) - N) * sp.diff(U, y) - x * sp.diff(N, y) * sp.diff(U, x)
    check("cleared Jacobian Bezout identity", zero(bezout - kappa * x**2))
    check("axis derivative U_x=a*c", zero(sp.diff(U, x).subs(x, 0) - a * c))
    check("axis derivative U_y=0", zero(sp.diff(U, y).subs(x, 0)))

    # Three pairwise-disjoint reduced components for a*b*c*d != 0.
    y_on_A = -a / (b * x)
    y_on_B = -(c + a * d * x**3) / (b * d * x**4)
    check("A on x=0", zero(A.subs(x, 0) - a))
    check("B on x=0", zero(B.subs(x, 0) - c))
    check("B on A=0", zero(B.subs(y, y_on_A) - c))
    check("A on B=0", zero(A.subs(y, y_on_B) + c / (d * x**3)))
    check("B_y on B component is nonzero monomial", zero(sp.diff(B, y) - b * d * x**4))

    # Componentwise U-principal coefficients of P are W|D.
    check("x-component coefficient", zero(W.subs(x, 0) - a * c**2))
    check("A-component coefficient", zero(W.subs(y, y_on_A)))
    check("B-component coefficient", zero(W.subs(y, y_on_B)))

    # Exact conic and quadratic target-field cover.
    S = c / x
    tau = S - P
    check("conic equation", zero(S**2 - P**2 + 4 * d * U))
    check("split conic product", zero(tau * (S + P) + 4 * d * U))
    check("conic parameter in source", zero(tau + 2 * d * x**2 * A))
    check("recover P from conic parameter", zero(P + tau / 2 + 2 * d * U / tau))
    check("quadratic x equation", zero(x**2 * (P**2 - 4 * d * U) - c**2))
    radicand = R**2 - 4 * d * L
    quadratic = radicand * Z**2 - c**2
    check("radicand is linear in target L", sp.Poly(radicand, L).degree() == 1)
    check("radicand is squarefree in target L", sp.gcd(sp.Poly(radicand, L), sp.Poly(sp.diff(radicand, L), L)).degree() == 0)
    check("quadratic discriminant", zero(sp.discriminant(quadratic, Z) - 4 * c**2 * radicand))

    # Recovery of the whole source field from (P,U,x).
    B_recovered = (c + x * P) / 2
    A_recovered = 2 * U / (x * (c + x * P))
    y_recovered = (A_recovered - a) / (b * x)
    check("recover B", zero(B_recovered - B))
    check("recover A", zero(A_recovered - A))
    check("recover y", zero(y_recovered - y))

    # Explicit quadratic deck involution in the birational (x,A) chart.
    x2 = -x
    A2 = -A * B / (c - B)
    B2 = c - B
    check("deck respects B definition", zero(B2 - (c + d * x2**3 * A2)))
    check("deck fixes U", zero(x2 * A2 * B2 - U))
    check("deck fixes P", zero((2 * B2 - c) / x2 - P))
    A3 = -A2 * B2 / (c - B2)
    check("deck is involutive", zero(A3 - A))

    # Exponent-three uniqueness in B_n=c+d*x^n*A, with the same P_n.
    exponent_rows = []
    for n in range(0, 8):
        Bn = c + d * x**n * A
        Un = x * A * Bn
        Pn = (2 * Bn - c) / x
        formula = -b * (c**2 + 2 * (3 - n) * Bn * (Bn - c))
        check(f"exponent formula n={n}", zero(jac(Pn, Un, x, y) - formula))
        exponent_rows.append((n, sp.factor(formula)))
    check("n=3 is constant", zero(exponent_rows[3][1] + b * c**2))

    # Sharp degenerations.
    check("d=0 becomes birational x=c/P", zero((x - c / P.subs(d, 0))))
    check("a=0 makes x-axis critical U_x", zero(sp.diff(U, x).subs({a: 0, x: 0})))
    check("a=0 makes x-axis critical U_y", zero(sp.diff(U, y).subs({a: 0, x: 0})))
    check("b=0 kills response", zero(kappa.subs(b, 0)))
    check("c=0 kills response", zero(kappa.subs(c, 0)))

    # Direct normalized smoothness and bounded polynomial-mate controls.
    U_control = sp.expand(U.subs({a: 1, b: 1, c: 1, d: 1}))
    critical_gb = sp.groebner([sp.diff(U_control, x), sp.diff(U_control, y)], x, y)
    check("normalized gradient ideal is unit", critical_gb.reduce(sp.Integer(1))[1] == 0)
    ranks = []
    for degree in range(0, 7):
        row = bounded_mate_rank(U_control, x, y, degree)
        check(f"bounded no polynomial mate degree {degree}", row[2] > row[1])
        ranks.append((degree,) + row)

    # A smallest nonlinear THM-3771 control.  Its unequal vertical addresses
    # do not survive arbitrary two-output rational target changes: the exact
    # inverse returns the polynomial source coordinates.
    xc, tc, pv, qv = sp.symbols("xc tc pv qv")
    Uc = xc
    Wc = xc + 3 * xc * tc + 1
    Qc = xc * (Wc - 2)
    Pc = sp.cancel(-Wc / (3 * Qc))
    check("cubic carrier constant Jacobian", zero(jac(Pc, Qc, xc, tc) - 1))
    check("cubic unequal address on U=0", zero((-Wc / 3).subs(xc, 0) + sp.Rational(1, 3)))
    check("cubic unequal address on W=2", zero((-Wc / 3).subs(tc, (1 - xc) / (3 * xc)) + sp.Rational(2, 3)))

    target_denominator = 3 * pv * qv + 2
    X_inverse = -qv / target_denominator
    T_inverse = (3 * pv * qv + 1) * target_denominator / (3 * qv) - sp.Rational(1, 3)
    check("cubic inverse target Jacobian", zero(jac(X_inverse, T_inverse, pv, qv) - 1))
    check("cubic inverse recovers x", zero(X_inverse.subs({pv: Pc, qv: Qc}, simultaneous=True) - xc))
    check("cubic inverse recovers t", zero(T_inverse.subs({pv: Pc, qv: Qc}, simultaneous=True) - tc))
    check("forward recovers p", zero(Pc.subs({xc: X_inverse, tc: T_inverse}, simultaneous=True) - pv))
    check("forward recovers q", zero(Qc.subs({xc: X_inverse, tc: T_inverse}, simultaneous=True) - qv))

    print("JC VERTICAL SPECTRUM / TARGET-WORD TRIAGE")
    print("SCOPE: exact characteristic-zero identities; nonradial three-component branch assumes a*b*c*d != 0")
    print("NONRADIAL_CONIC:")
    print("  A=a+b*x*y; B=c+d*x^3*A; U=x*A*B; P=(2*B-c)/x")
    print("  J(P,U)=-b*c^2; J(U,UP)=b*c^2*U")
    print("  vertical coefficient vector on (x,A,B)=(a*c^2,0,0)")
    print("  (P^2-4*d*U)*x^2=c^2; source recovery needs exactly quadratic x")
    print("  deck: x->-x, B->c-B, A->-A*B/(c-B)")
    print("  consequence: every birational symplectic target word preserves cover degree 2")
    print("EXPONENT_UNIQUENESS:")
    print("  J(P_n,U_n)=-b*[c^2+2*(3-n)*B_n*(B_n-c)]")
    print("  n=3 is the unique genuine d!=0 constant-response exponent")
    print("CUBIC_CARRIER_CONTROL:")
    print("  q=x*[x+3*x*t-1]; p=-[x+3*x*t+1]/(3*q); addresses (-1/3,-2/3)")
    print("  inverse: x=-q/(3*p*q+2), t=(3*p*q+1)*(3*p*q+2)/(3*q)-1/3")
    print("  consequence: unequal residues do not obstruct arbitrary rational target maps")
    print("BOUNDED_NONRADIAL_MATE_RANKS: degree,unknowns,rank,augmented_rank")
    for row in ranks:
        print("  " + ",".join(map(str, row)))
    print(f"CHECKS={CHECKS}")
    print("RESULT: PASS")


if __name__ == "__main__":
    main()
