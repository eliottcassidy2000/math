#!/usr/bin/env python3
"""Independent exact hostile controls for the proposed THM-3759.

The theorem is proved by hand in REPORT.md.  This script checks unrelated
symbolic consequences of that proof: the axis/interior boundary, the general
constant-Jacobian recurrence, the polynomial Poisson centralizer in bounded
degree, and the source-foliation degree distinction.
"""

import sympy as sp


X, T, z, u, v = sp.symbols("X T z u v")


def require(test, label):
    if not bool(test):
        raise RuntimeError(label)


def jac(P, Q):
    return sp.expand(sp.diff(P, X) * sp.diff(Q, T) - sp.diff(P, T) * sp.diff(Q, X))


def ideal_is_unit(polys):
    G = sp.groebner(polys, X, T, domain=sp.QQ)
    return G.reduce(sp.Integer(1))[1] == 0


def euler_inverse(poly, shift):
    ans = 0
    for (j,), coeff in sp.Poly(sp.expand(poly), z, domain=sp.QQ).terms():
        ans += coeff * z**j / sp.Rational(j + shift)
    return sp.expand(ans)


def forced_chain(a, c, chi, depth):
    cp = sp.diff(chi, z)
    F = {1: -sp.Rational(c, a)}
    for k in range(1, depth):
        rhs = -sp.Rational(k, a) * F[k] * cp
        F[k + 1] = euler_inverse(rhs, k + 1)
        require(F[k + 1] != 0, f"zero repair layer {k + 1}")
    Q = a * X + chi.subs(z, X * T)
    P = sum(T**k * F[k].subs(z, X * T) for k in range(1, depth + 1))
    terminal = -depth * T**depth * F[depth].subs(z, X * T) * cp.subs(z, X * T)
    require(sp.expand(jac(P, Q) - c - terminal) == 0, "forced-chain identity")
    return [sp.degree(F[k], z) for k in range(1, depth + 1)]


def linear_system_nullity(Q, degree, rhs):
    mons = [X**i * T**j for i in range(degree + 1) for j in range(degree + 1 - i)]
    coeff = sp.symbols(f"q0:{len(mons)}")
    R = sum(q * m for q, m in zip(coeff, mons))
    defect = sp.Poly(sp.expand(jac(R, Q) - rhs), X, T)
    equations = [sp.expand(e) for e in defect.coeffs()]
    A, bb = sp.linear_eq_to_matrix(equations, coeff)
    consistent = A.rank() == A.row_join(bb).rank()
    return consistent, len(coeff) - A.rank()


def main():
    aa, bb = sp.symbols("aa bb")
    chi_generic = z + 2 * z**2 + 3 * z**4
    Q_generic = aa * X + bb * T + chi_generic.subs(z, X * T)
    cp = sp.diff(chi_generic, z).subs(z, X * T)
    require(sp.expand(sp.diff(Q_generic, X) - aa - T * cp) == 0, "gradient X")
    require(sp.expand(sp.diff(Q_generic, T) - bb - X * cp) == 0, "gradient T")

    # Algebraically closed hostile controls: two flanks are singular; one
    # flank is smooth exactly when chi'(0)=0.
    interior = []
    for d in range(1, 7):
        chi = z**d + 2 * z ** (d + 1)
        Q = 2 * X + 3 * T + chi.subs(z, X * T)
        unit = ideal_is_unit([sp.diff(Q, X), sp.diff(Q, T)])
        require(not unit, f"interior unexpectedly smooth at d={d}")
        R = sp.expand(z * sp.diff(chi, z) ** 2 - 6)
        require(sp.degree(R, z) >= 1, "interior eliminant constant")
        require(sp.degree(sp.gcd(R, z * sp.diff(chi, z)), z) == 0,
                "root met excluded locus")
        interior.append(sp.degree(R, z))

    smooth_axes = []
    for d in range(2, 8):
        chi = z**d + z ** (d + 1)
        Q = 2 * X + chi.subs(z, X * T)
        unit = ideal_is_unit([sp.diff(Q, X), sp.diff(Q, T)])
        require(unit, f"axis unexpectedly singular at d={d}")
        smooth_axes.append(d)
    bad_axis = 2 * X + (z + z**4).subs(z, X * T)
    require(not ideal_is_unit([sp.diff(bad_axis, X), sp.diff(bad_axis, T)]),
            "chi'(0)!=0 boundary unexpectedly smooth")

    chains = []
    for chi in (z**2, z**2 + z**3, z**2 - z**4 + z**6):
        chains.append(forced_chain(a=2, c=3, chi=chi, depth=10))

    # Independent bounded centralizer/mate test for Q=X+(XT)^2.  The hand
    # descent predicts ker J(-,Q)=Q[Q].  In total degree N its dimension is
    # floor(N/4)+1, while J(R,Q)=1 is inconsistent at every N.
    Q0 = X + (X * T) ** 2
    centralizer = []
    for degree in range(0, 11):
        kernel_consistent, nullity = linear_system_nullity(Q0, degree, 0)
        mate_consistent, _ = linear_system_nullity(Q0, degree, 1)
        expected = degree // 4 + 1
        require(kernel_consistent and nullity == expected,
                f"centralizer mismatch at degree={degree}: {nullity} vs {expected}")
        require(not mate_consistent, f"bounded mate at degree={degree}")
        centralizer.append(nullity)

    # Foliation degree is degree in v over Q[u], not degree on every special
    # affine line.  Nonaxis and axis directions give 8 and 4 for d=4.
    Q4 = X + (X * T) ** 2 + (X * T) ** 4
    nonaxis = sp.expand(Q4.subs({X: u + v, T: 2 * u + 3 * v}, simultaneous=True))
    axis = sp.expand(Q4.subs({X: u + v, T: 2 * u + 1}, simultaneous=True))
    require(sp.degree(nonaxis, v) == 8, "nonaxis foliation degree")
    require(sp.degree(axis, v) == 4, "axis foliation degree")
    special_line = sp.expand(Q4.subs({X: v, T: 0}, simultaneous=True))
    require(sp.degree(special_line, v) == 1, "special-fibre degree-drop hostile")

    print("INTERIOR_ELIMINANT_DEGREES", interior)
    print("SMOOTH_AXIS_PROFILE_DEGREES", smooth_axes)
    print("FORCED_CHAIN_DEGREES", chains)
    print("CENTRALIZER_NULLITIES_DEG_0_TO_10", centralizer)
    print("FOLIATION_DEGREES", sp.degree(nonaxis, v), sp.degree(axis, v))
    print("SPECIAL_LINE_T_EQ_0_DEGREE", sp.degree(special_line, v))
    print("RESULT: PASS")


if __name__ == "__main__":
    main()

