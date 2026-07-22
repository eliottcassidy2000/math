#!/usr/bin/env python3
"""Exact symbolic controls for THM-2102's weighted-face descent."""
import sympy as sp


x, y = sp.symbols("x y")


def jac(f, g):
    return sp.expand(sp.diff(f, x) * sp.diff(g, y) - sp.diff(f, y) * sp.diff(g, x))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def weighted_component(f, weights, degree):
    a, b = weights
    poly = sp.Poly(sp.expand(f), x, y)
    answer = 0
    for (i, j), coefficient in poly.terms():
        if a * i + b * j == degree:
            answer += coefficient * x**i * y**j
    return sp.expand(answer)


def terminal_face_census():
    checks = 0
    for a in range(1, 5):
        for b in range(1, 5):
            W = a + b
            for D in range(1, W):
                E = W - D
                fmons = [x**i * y**j for i in range(W + 1) for j in range(W + 1)
                         if a * i + b * j == D]
                gmons = [x**i * y**j for i in range(W + 1) for j in range(W + 1)
                         if a * i + b * j == E]
                F = sum((i + 2) * monomial for i, monomial in enumerate(fmons))
                G = sum((i + 3) * monomial for i, monomial in enumerate(gmons))
                expected = sp.Poly(F, x, y).coeff_monomial(x) * sp.Poly(G, x, y).coeff_monomial(y)
                expected -= sp.Poly(F, x, y).coeff_monomial(y) * sp.Poly(G, x, y).coeff_monomial(x)
                require(jac(F, G) == expected, (a, b, D, E, F, G, jac(F, G), expected))
                checks += 1
    return checks


def main():
    # Power-free top face with lower terms: already triangular.
    f_triangular = y + x**3 + 2 * x**2 - 5
    require(jac(f_triangular, -x) == 1, "triangular power-free control")

    # Proper-power coordinate: terminal first defect exposes H and -x.
    H = y + x**2
    f = x + H**2
    g = H
    require(jac(f, g) == 1, "proper-power coordinate")
    A = B = sp.Integer(1)
    m, n = 2, 1
    P, Q = x, sp.Integer(0)
    L = sp.expand(A * m * H ** (m - 1) * Q - B * n * H ** (n - 1) * P)
    require(L == -x and jac(H, L) == 1, (L, jac(H, L)))

    # A genuine common approximate-root correction.
    H3 = y + x**3
    R = x**2
    K = H3 + R
    f3 = x + K**3
    g3 = K
    require(jac(f3, g3) == 1, "approximate-root coordinate")
    P3 = weighted_component(sp.expand(K**3), (1, 3), 8)
    Q3 = weighted_component(K, (1, 3), 2)
    require(sp.expand(P3 - 3 * H3**2 * R) == 0 and Q3 == R, (P3, Q3))
    L3 = sp.expand(3 * H3**2 * Q3 - P3)
    require(L3 == 0, L3)
    require(sp.cancel(P3 / (3 * H3**2)) == R and Q3 == R, "unique approximate root")

    # Synchronization without H-divisibility: the next bracket survives.
    P_bad = x**3
    Q_bad = sp.Rational(3, 2) * H * x**3
    L_bad = sp.expand(2 * H * Q_bad - 3 * H**2 * P_bad)
    bad_bracket = jac(H**2 + P_bad, H**3 + Q_bad)
    require(L_bad == 0, L_bad)
    require(sp.rem(sp.Poly(P_bad, y), sp.Poly(H, y)) != 0, "H should not divide P")
    require(bad_bracket == sp.Rational(9, 2) * x**5, bad_bracket)

    face_checks = terminal_face_census()
    print("triangular_power_free_jacobian=1")
    print("proper_power_terminal_L=-x terminal_jacobian=1")
    print("approximate_root_absorbed=True")
    print("synchronized_nonlift_bracket=9*x^5/2")
    print(f"terminal_face_checks={face_checks}")
    print("PASS")


if __name__ == "__main__":
    main()
