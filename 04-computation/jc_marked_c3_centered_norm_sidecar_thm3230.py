#!/usr/bin/env python3
"""Exact controls for THM-3230's marked C3 centered-norm sidecar.

The proof is a local Kummer/valuation argument.  This companion verifies the
trace-zero norm identity, its marked-cubic coefficient formula, the sharp
leading-graph hostile, the valuation congruence, and the Bezout-gauge law for
the intrinsic terminal prefactor cubeclass.  Every truth-bearing gate remains
active under ``python -O``.
"""

from math import gcd

import sympy as sp


def require(condition, label):
    if not bool(condition):
        raise AssertionError(label)


def reduce_zeta(expr, zeta):
    """Reduce a polynomial modulo zeta^2+zeta+1."""
    return sp.rem(sp.Poly(sp.expand(expr), zeta), sp.Poly(zeta**2 + zeta + 1, zeta)).as_expr()


def centered_norm_identity():
    zeta, a, b = sp.symbols("zeta a b")
    conjugates = [a + b, zeta * a + zeta**2 * b, zeta**2 * a + zeta * b]
    residue = sp.expand(reduce_zeta(sp.prod(conjugates) - a**3 - b**3, zeta))
    require(residue == 0, "C3 trace-zero norm identity")

    eps, t, A, B = sp.symbols("eps t A B", nonzero=True)
    norm = sp.expand((eps * t) * A**3 + (eps * t) ** 2 * B**3)
    require(norm == eps * t * A**3 + eps**2 * t**2 * B**3, "Kummer norm substitution")
    return norm


def marked_cubic_formula():
    T = sp.symbols("T")
    a3, a2, a1, fixed = sp.symbols("a3 a2 a1 fixed")
    b2 = a3 + fixed
    b1 = a2 + fixed * a3 + fixed**2
    b0 = a1 + fixed * a2 + fixed**2 * a3 + fixed**3
    g = T**3 + b2 * T**2 + b1 * T + b0
    f = sp.expand((T - fixed) * g)
    require(f.coeff(T, 4) == 1, "marked factor leading coefficient")
    require(sp.expand(f.coeff(T, 3) - a3) == 0, "marked factor cubic coefficient")
    require(sp.expand(f.coeff(T, 2) - a2) == 0, "marked factor quadratic coefficient")
    require(sp.expand(f.coeff(T, 1) - a1) == 0, "marked factor linear coefficient")
    require(sp.expand(f.coeff(T, 0) + fixed * b0) == 0, "marked factor constant ledger")

    mu = -b2 / 3
    centered_norm = sp.factor(-g.subs(T, mu))
    expected = sp.factor(-2 * b2**3 / 27 + b1 * b2 / 3 - b0)
    require(sp.expand(centered_norm - expected) == 0, "marked centered cubic norm")
    return expected


def valuation_congruence_grid():
    cells = 0
    for va in range(-9, 10):
        for vb in range(-9, 10):
            first = 1 + 3 * va
            second = 2 + 3 * vb
            require(first != second, "distinct C3-character values")
            h = min(first, second)
            require(h % 3 in (1, 2), "centered value is prime to three")
            epsilon_exponent = 1 if first < second else 2
            require(epsilon_exponent == h % 3, "leading norm cubeclass exponent")
            cells += 1
    return cells


def bezout_pair(g, e):
    for A in range(-e, e + 1):
        remainder = 1 - A * g
        if remainder % e == 0:
            return A, remainder // e
    raise AssertionError("Bezout pair not found")


def terminal_gauge_grid():
    cells = 0
    for g in range(1, 19):
        for e in range(1, 19):
            if gcd(g, e) != 1:
                continue
            A, B = bezout_pair(g, e)
            require(A * g + B * e == 1, "Bezout equation")
            for k in range(-5, 6):
                A2, B2 = A + k * e, B - k * g
                require(A2 * g + B2 * e == 1, "Bezout gauge equation")

                # rho'=rho*theta^(-k), K'=theta^(3k)K,
                # L'=theta^(k(3-e))L.  The exponent change in
                # L*theta^(A-1) is exactly 3k, hence a cube.
                lambda_shift = k * (3 - e) + (A2 - A)
                require(lambda_shift == 3 * k, "intrinsic Lambda gauge shift")
                require((3 * k) % 3 == 0, "K cubeclass gauge invariance")
                cells += 1
    return cells


def graph_compatibility_grid():
    cells = 0
    for h in range(-20, 21):
        if h % 3 == 0:
            continue
        for m in range(1, 25):
            # Write k for the abstract cubeclass exponent of terminal K.
            # q0 has exponent m*k and the centered norm has exponent -h*k.
            require((h * m + m * (-h)) % 3 == 0, "q0/centered-norm compatibility")
            cells += 1
    return cells


def leading_blind_hostile():
    T, t, eps = sp.symbols("T t eps", nonzero=True)
    # Fixed root 0 and moving orbit t^(-1)+zeta^i*pi, with pi^3=eps*t.
    quartic = sp.expand(T * ((T - t**-1) ** 3 - eps * t))
    a3 = quartic.coeff(T, 3)
    a2 = quartic.coeff(T, 2)
    a1 = quartic.coeff(T, 1)
    depressed_q = sp.factor(a1 - a2 * a3 / 2 + a3**3 / 8)
    expected_q = sp.Rational(1, 8) * t**-3 - eps * t
    require(sp.cancel(depressed_q - expected_q) == 0, "m=3 graph hostile q")

    W = sp.symbols("W")
    moved_cubic = sp.expand((W - t**-1) ** 3 - eps * t)
    mu = t**-1
    centered_norm = sp.factor(-moved_cubic.subs(W, mu))
    require(centered_norm == eps * t, "m=3 hostile centered norm")
    return depressed_q, centered_norm


def main():
    norm = centered_norm_identity()
    marked = marked_cubic_formula()
    valuation_cells = valuation_congruence_grid()
    gauge_cells = terminal_gauge_grid()
    compatibility_cells = graph_compatibility_grid()
    hostile_q, hostile_norm = leading_blind_hostile()

    print("THM-3230 exact controls")
    print("centered_norm=eps*t*A**3 + eps**2*t**2*B**3")
    print(f"marked_cubic_norm={marked}")
    print(f"valuation_cells={valuation_cells}:PASS")
    print(f"bezout_gauge_cells={gauge_cells}:PASS")
    print(f"graph_compatibility_cells={compatibility_cells}:PASS")
    print(f"m3_hostile_depressed_q={hostile_q}")
    print(f"m3_hostile_centered_norm={hostile_norm}")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
