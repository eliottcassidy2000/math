#!/usr/bin/env python3
"""Exact controls for THM-2805."""

import ast
from itertools import combinations
from pathlib import Path

import sympy as sp


def require(condition, message):
    if not bool(condition):
        raise RuntimeError(message)


def has_asserts(path):
    return any(
        isinstance(node, ast.Assert)
        for node in ast.walk(ast.parse(Path(path).read_text(encoding="utf-8")))
    )


def reduce_mod(expr, variable, modulus):
    num, den = sp.cancel(expr).as_numer_denom()
    mod_poly = sp.Poly(modulus, variable, domain=sp.QQ)
    den_poly = sp.Poly(den, variable, domain=sp.QQ)
    require(sp.gcd(den_poly, mod_poly).degree() == 0, "nonunit denominator")
    inv = sp.invert(den_poly, mod_poly)
    return sp.rem(sp.Poly(num, variable, domain=sp.QQ) * inv, mod_poly).as_expr()


def reduce_x(expr, x, s, modulus):
    return all(
        reduce_mod(coef, s, modulus) == 0
        for coef in sp.Poly(sp.expand(expr), x).all_coeffs()
    )


def coprime_to(expr, s, modulus):
    num = sp.cancel(expr).as_numer_denom()[0]
    return sp.gcd(
        sp.Poly(num, s, domain=sp.QQ),
        sp.Poly(modulus, s, domain=sp.QQ),
    ).degree() == 0


def data_for(N, d, x, s):
    b = N - d
    p = sp.cancel((d * (N - 1) * s - d * (d - 1)) / (N * (N - 1)))
    h = [sp.Integer(1), s]
    for _ in range(2, N):
        h.append(sp.expand(s * h[-1] - p * h[-2]))

    coeffs = []
    H = sp.Integer(0)
    A_sum = sp.Integer(0)
    for j in range(b + 1):
        power = d + j - 1
        coeff = (d + j) * (-1) ** (b - j) * sp.binomial(b, j)
        coeffs.append((power, coeff))
        H += coeff * (h[power - 1] if power >= 1 else 0)
        A_sum += coeff * (h[power - 2] if power >= 2 else 0)
    H = sp.expand(H)
    A = sp.expand(p * A_sum)

    diagonal = sp.expand(
        N * (N - 1) * s**2
        - 4 * d * (N - 1) * s
        + 4 * d * (d - 1)
    )
    Q = sp.Poly(sp.cancel(H / diagonal), s, domain=sp.QQ).as_expr()
    z = sp.cancel(((N - 1) * s - (d - 1)) / N)

    E = sp.expand(x**2 - s * x + p)
    D = sp.expand(x**d * (x - 1) ** b)
    T = x * (x - 1)
    B = sp.expand(D + A * (x - z))
    C = sp.expand(-(N - 1) * A)
    S, remainder = sp.div(B, E**2, x)

    return {
        "b": b,
        "p": p,
        "h": h,
        "H": H,
        "diagonal": diagonal,
        "Q": Q,
        "z": z,
        "E": E,
        "D": D,
        "T": T,
        "A": A,
        "B": B,
        "C": C,
        "S": S,
        "remainder": remainder,
    }


def algebra_check(N, d, x, s):
    data = data_for(N, d, x, s)
    Q = data["Q"]
    require(sp.expand(data["H"] - data["diagonal"] * Q) == 0, "diagonal quotient")
    require(sp.degree(Q, s) == N - 4, "Maxwell degree")

    if N == 4:
        require(sp.degree(Q, s) == 0 and Q != 0, "empty N=4 middle chamber")
        return data

    require(sp.discriminant(Q, s) != 0, "bounded Q squarefree")
    require(reduce_x(data["remainder"], x, s, Q), "E^2 division")
    require(sp.degree(data["S"], x) == N - 4, "S degree")
    require(
        reduce_x(
            sp.diff(data["B"], x) * data["D"]
            - data["B"] * sp.diff(data["D"], x)
            - data["C"]
            * data["E"]
            * x ** (d - 1)
            * (x - 1) ** (N - d - 1),
            x,
            s,
            Q,
        ),
        "response derivative",
    )

    disc_e = sp.discriminant(data["E"], x)
    disc_s = sp.discriminant(data["S"], x) if N >= 6 else sp.Integer(1)
    res_se = sp.resultant(data["S"], data["E"] * x * (x - 1), x)
    res_e = sp.resultant(data["E"], x * (x - 1), x)
    require(
        all(
            coprime_to(gate, s, Q)
            for gate in (disc_e, data["A"], disc_s, res_se, res_e)
        ),
        "bounded converse gates",
    )
    if d == N - d:
        require(sp.expand(Q.subs(s, 1)) != 0, "equal-pole orientation fixed point")
    return data


def nielsen_count(N, d):
    b = N - d
    L = N - 1
    distances = {d, b}
    placements = [
        (left, right)
        for left, right in combinations(range(1, L), 2)
        if right - left in distances
    ]
    if d != b:
        require(len(placements) == N - 4, "unequal-pole Nielsen count")
        unmarked = len(placements)
        ordered = len(placements)
    else:
        require(len(placements) == N // 2 - 2, "equal-pole Nielsen count")
        unmarked = len(placements)
        ordered = 2 * len(placements)
        require(ordered == N - 4, "equal-pole ordered count")
    return unmarked, ordered


def main():
    require(not has_asserts(Path(__file__)), "truth-bearing Python assert found")
    x, s = sp.symbols("x s")
    print("THM-2805 GENERAL TWO-POLE E=2 MAXWELL AUDIT")
    print("N,d | degree Q | unmarked classes | ordered charts")

    cases = 0
    ordered_total = 0
    for N in range(4, 10):
        for d in range(2, N // 2 + 1):
            data = algebra_check(N, d, x, s)
            unmarked, ordered = nielsen_count(N, d)
            degree_q = int(sp.degree(data["Q"], s))
            require(degree_q == ordered, "algebra/Nielsen saturation")
            print(
                f"{N:2d},{d:1d} | {degree_q:8d} |"
                f" {unmarked:16d} | {ordered:14d}"
            )
            cases += 1
            ordered_total += ordered

            # Complementary pole charts are related by x -> 1-x, s -> 2-s.
            complement = data_for(N, N - d, x, s)
            transformed = sp.Poly(
                sp.expand(complement["Q"].subs(s, 2 - s)), s, domain=sp.QQ
            )
            original = sp.Poly(data["Q"], s, domain=sp.QQ)
            require(
                sp.cancel(transformed.LC() * original.as_expr() - original.LC() * transformed.as_expr())
                == 0,
                "pole-swap Maxwell symmetry",
            )

    # Combined with THM-2800's d=1 edge, the total h=2 count is binomial.
    for N in range(4, 31):
        total = N - 3
        for d in range(2, N // 2 + 1):
            total += nielsen_count(N, d)[0]
        require(total == (N - 2) * (N - 3) // 2, "total h=2 count")

    print(f"exact algebra/pole cases checked: {cases}")
    print(f"ordered Maxwell charts in checked middle corridors: {ordered_total}")
    print("all h=2 classes including d=1: binomial(N-2,2)")
    print("equal-pole charts pair under s -> 2-s with no fixed root")
    print("bitangent quotient and converse gates: PASS")
    print("assert_nodes=0")


if __name__ == "__main__":
    main()
