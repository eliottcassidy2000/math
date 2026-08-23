#!/usr/bin/env python3
"""Exact scratch referee for Q=aX+bT+chi(XT).

This is a bounded control surface for an all-degree hand proof.  It uses
explicit exceptions rather than Python assertions so normal and optimized
runs exercise the same checks.
"""

from itertools import product

import sympy as sp


X, T, z = sp.symbols("X T z")


def check(condition, label):
    if not condition:
        raise RuntimeError(label)


def jacobian(P, Q):
    return sp.expand(sp.diff(P, X) * sp.diff(Q, T) - sp.diff(P, T) * sp.diff(Q, X))


def euler_inverse(poly, k):
    """Invert z*d/dz+k on Q[z]."""
    out = 0
    for (j,), coeff in sp.Poly(sp.expand(poly), z, domain=sp.QQ).terms():
        out += coeff * z**j / sp.Rational(j + k)
    return sp.expand(out)


def generic_identities():
    a, b = sp.symbols("a b", nonzero=True)
    c = sp.symbols("c0:5")
    chi = sum(c[j] * z**j for j in range(5))
    cp = sp.diff(chi, z)
    Q = a * X + b * T + chi.subs(z, X * T)
    check(sp.expand(sp.diff(Q, X) - (a + T * cp.subs(z, X * T))) == 0, "Q_X")
    check(sp.expand(sp.diff(Q, T) - (b + X * cp.subs(z, X * T))) == 0, "Q_T")

    # At a root z*cp^2=ab, X=-b/cp and T=-a/cp reconstruct a critical point.
    cp0 = sp.symbols("cp0", nonzero=True)
    x0, t0 = -b / cp0, -a / cp0
    check(sp.simplify(a + t0 * cp0) == 0, "critical Q_X")
    check(sp.simplify(b + x0 * cp0) == 0, "critical Q_T")
    check(sp.simplify(x0 * t0 - a * b / cp0**2) == 0, "critical product")

    # General three-sector torus eliminant.
    ph = 1 + 2 * z + 3 * z**2
    ps = 2 - z + z**2
    ch = z + 4 * z**2 + z**3
    w = sp.expand(z * ph * ps)
    E = sp.expand(sp.diff(w, z) ** 2 - w * sp.diff(ch, z) ** 2)
    xx = sp.symbols("xx", nonzero=True)
    # F(X,z)=X phi+chi+z psi/X; eliminate X from F_X=F_z=0.
    FX = sp.diff(xx * ph + ch + z * ps / xx, xx)
    Fz = sp.diff(xx * ph + ch + z * ps / xx, z)
    res = sp.factor(sp.resultant(sp.together(FX * xx**2), sp.together(Fz * xx), xx))
    # The resultant differs from E only by an explicit nonzero profile factor/sign.
    check(sp.simplify(res / E) != 0 and sp.rem(sp.Poly(res, z), sp.Poly(E, z)) == 0,
          "general torus eliminant divisor")
    return E, sp.factor(res / E)


def forced_weight_chain(chi, depth=8):
    """Build the unique negative-weight repair for Q=X+chi(XT)."""
    cp = sp.diff(chi, z)
    Fs = {1: sp.Integer(-1)}
    for k in range(2, depth + 1):
        rhs = sp.expand(-(k - 1) * Fs[k - 1] * cp)
        Fs[k] = euler_inverse(rhs, k)
        check(Fs[k] != 0, f"forced F_{k} vanished")
        check(sp.expand(z * sp.diff(Fs[k], z) + k * Fs[k] - rhs) == 0,
              f"Euler recurrence k={k}")

    Q = X + chi.subs(z, X * T)
    P = sum(T**k * Fs[k].subs(z, X * T) for k in range(1, depth + 1))
    debt = sp.expand(1 - depth * T**depth * Fs[depth].subs(z, X * T) * cp.subs(z, X * T))
    check(sp.expand(jacobian(P, Q) - debt) == 0, "truncated repair debt")
    return Fs, debt


def direct_mate_system(chi, degree):
    """Return whether a total-degree <=degree mate exists for Q=X+chi(XT)."""
    mons = [X**i * T**j for i in range(degree + 1) for j in range(degree + 1 - i)]
    coeffs = sp.symbols(f"u0:{len(mons)}")
    P = sum(c * m for c, m in zip(coeffs, mons))
    Q = X + chi.subs(z, X * T)
    defect = sp.Poly(sp.expand(jacobian(P, Q) - 1), X, T)
    equations = list(defect.coeffs())
    solution = sp.linsolve(equations, coeffs)
    return solution != sp.EmptySet


def smoothness_census():
    checked = smooth = singular = 0
    mismatches = []
    for a, b in product((-1, 0, 1), repeat=2):
        for c1, c2, c3 in product((-1, 0, 1), repeat=3):
            chi = c1 * z + c2 * z**2 + c3 * z**3
            Q = sp.expand(a * X + b * T + chi.subs(z, X * T))
            G = sp.groebner([sp.diff(Q, X), sp.diff(Q, T)], X, T, domain=sp.QQ)
            actual = G.reduce(sp.Integer(1))[1] == 0
            nonconstant = any((c1, c2, c3))
            if a != 0 and b != 0:
                predicted = not nonconstant
            elif (a != 0) ^ (b != 0):
                predicted = c1 == 0
            else:
                predicted = False
            checked += 1
            smooth += int(actual)
            singular += int(not actual)
            if actual != predicted:
                mismatches.append((a, b, c1, c2, c3, actual, predicted))
    check(not mismatches, f"smoothness mismatches: {mismatches[:3]}")
    return checked, smooth, singular


def main():
    E, quotient = generic_identities()
    print("GENERIC_TORUS_ELIMINANT_DEGREE", sp.degree(E, z))
    print("GENERIC_RESULTANT_QUOTIENT", quotient)

    profiles = [z**2, z**2 + z**3, z**2 + z**4, z**2 - z**4 + z**5]
    for chi in profiles:
        Fs, debt = forced_weight_chain(chi, depth=8)
        degrees = [sp.degree(Fs[k], z) for k in range(1, 9)]
        print("PROFILE", chi, "FORCED_DEGREES", degrees)
        print("FINAL_DEBT_TERMS", len(sp.Poly(debt, X, T).terms()))
        for degree in range(0, 7):
            check(not direct_mate_system(chi, degree),
                  f"unexpected mate for chi={chi}, degree={degree}")
        print("NO_MATE_TOTAL_DEGREE_LE_6", chi)

    checked, smooth, singular = smoothness_census()
    print("SMOOTHNESS_CENSUS", checked, "SMOOTH", smooth, "SINGULAR", singular)
    print("RESULT: PASS")


if __name__ == "__main__":
    main()

