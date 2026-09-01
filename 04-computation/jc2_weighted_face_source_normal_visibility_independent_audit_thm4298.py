#!/usr/bin/env python3
"""Independent SymPy audit for THM-4298.

This path starts from the generating-function substitution
H(z)=(1+z)^L F(z/(1+z)), reconstructs its coefficient matrix and inverse,
and separately checks the source-normal monomial expansion and M=12 walls.
It imports no code or data from the primary certificate.
"""

from __future__ import annotations

import sympy as sp


def require(condition: bool, label: str) -> None:
    if not condition:
        raise AssertionError(label)


def face(M: int) -> list[tuple[int, int]]:
    return [
        ((M - 3 * j) // 2, j)
        for j in range(M // 3 + 1)
        if M - 3 * j >= 0 and (M - 3 * j) % 2 == 0
    ]


def matrices(M: int) -> tuple[sp.Matrix, sp.Matrix]:
    monomials = face(M)
    size = len(monomials)
    forward = sp.zeros(size)
    inverse = sp.zeros(size)
    for s in range(size):
        for r, (i_r, j_r) in enumerate(monomials):
            q = s - r
            if 0 <= q <= i_r + j_r:
                forward[s, r] = sp.binomial(i_r + j_r, q)
                inverse[s, r] = (-1) ** q * sp.binomial(i_r + j_r, q)
    return forward, inverse


def audit_monomial_expansion() -> int:
    x, t = sp.symbols("x t")
    p = t * (1 + x**2 * t)
    y = x * t * p
    checks = 0
    for i in range(9):
        for j in range(9):
            actual = sp.Poly(sp.expand(p**i * y**j), x, t)
            expected = sp.Poly(
                sum(
                    sp.binomial(i + j, q) * x ** (j + 2 * q) * t ** (i + 2 * j + q)
                    for q in range(i + j + 1)
                ),
                x,
                t,
            )
            require(actual == expected, f"monomial expansion {(i, j)}")
            for (d, n), coefficient in actual.terms():
                require(coefficient != 0, "nonzero term")
                require(2 * n - d == 2 * i + 3 * j, "diagonal isolation")
                checks += 1
    return checks


def audit_generating_functions() -> tuple[int, int]:
    z, w = sp.symbols("z w")
    matrix_checks = 0
    symbolic_checks = 0
    for M in range(61):
        monomials = face(M)
        if not monomials:
            require(M == 1, "unique empty weight")
            continue
        forward, inverse = matrices(M)
        require(forward.det() == 1, f"determinant M={M}")
        require(forward * inverse == sp.eye(len(monomials)), f"right inverse M={M}")
        require(inverse * forward == sp.eye(len(monomials)), f"left inverse M={M}")
        rows = [(M + j) // 2 for _, j in monomials]
        require(rows[0] == (M + 1) // 2, "first row")
        require(rows[-1] == (2 * M) // 3, "last row")
        require(rows == list(range(rows[0], rows[-1] + 1)), "row interval")
        matrix_checks += len(monomials) ** 2

        if M <= 18:
            coefficients = sp.symbols(f"c0:{len(monomials)}")
            j0 = monomials[0][1]
            L = (M - j0) // 2
            F_poly = sum(coefficients[r] * w**r for r in range(len(monomials)))
            H_poly = sp.cancel((1 + z) ** L * F_poly.subs(w, z / (1 + z)))
            H_poly = sp.Poly(sp.expand(H_poly), z)
            h_vector = sp.Matrix(
                [H_poly.coeff_monomial(z**s) for s in range(len(monomials))]
            )
            require(
                h_vector == forward * sp.Matrix(coefficients),
                f"forward generating identity M={M}",
            )
            recovered = sp.cancel(
                (1 - w) ** L * H_poly.as_expr().subs(z, w / (1 - w))
            )
            require(
                sp.Poly(sp.expand(recovered - F_poly), w).is_zero,
                f"inverse generating identity M={M}",
            )
            symbolic_checks += len(monomials)
    return matrix_checks, symbolic_checks


def audit_m12() -> tuple[sp.Matrix, sp.Matrix]:
    U, W, Z = sp.symbols("U W Z")
    forward, inverse = matrices(12)
    require(face(12) == [(6, 0), (3, 2), (0, 4)], "M12 monomials")
    coefficients = sp.Matrix([U, W, Z])
    a, b, c = list(forward * coefficients)
    require(a == U and b == 6 * U + W and c == 15 * U + 5 * W + Z, "M12 channels")
    recovered = inverse * sp.Matrix([a, b, c])
    require(recovered == coefficients, "M12 recovery")
    require(sp.expand(U + W + Z - (c - 4 * b + 10 * a)) == 0, "Lambda wall")
    require(
        sp.expand(W**2 - 4 * U * Z - (b**2 + 8 * a * b - 24 * a**2 - 4 * a * c)) == 0,
        "D wall",
    )
    require(sp.Matrix([1, 6, 15]) != sp.Matrix([0, 1, 5]), "scalar-weight hostile")
    return forward, inverse


def main() -> None:
    monomial_checks = audit_monomial_expansion()
    matrix_checks, symbolic_checks = audit_generating_functions()
    forward, inverse = audit_m12()
    print("THM4298_WEIGHTED_FACE_VISIBILITY_INDEPENDENT_V1")
    print(f"MONOMIAL_TERMS={monomial_checks} DIAGONAL_IDENTITY=PASS")
    print(f"MATRIX_ENTRIES={matrix_checks} RANGE=0..60_except_1 DET=1 INVERSES=PASS")
    print(f"GENERATING_CHANNELS={symbolic_checks} RANGE=0..18 FORWARD_REVERSE=PASS")
    print(f"M12_TRANSFORM={forward.tolist()}")
    print(f"M12_INVERSE={inverse.tolist()}")
    print("M12_WALL_DICTIONARY=PASS")
    print("SCOPE lossless_face_flag_not_seam_entry_or_JC2")
    print("VERDICT PASS INDEPENDENT SYMBOLIC AUDIT")


if __name__ == "__main__":
    main()
