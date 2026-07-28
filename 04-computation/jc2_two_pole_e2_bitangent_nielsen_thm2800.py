#!/usr/bin/env python3
"""Exact controls for THM-2800.

The algebraic eliminant and the Nielsen count are checked independently.
All truth-bearing tests use ``require`` so optimized Python executes them.
"""

import ast
from itertools import combinations
from pathlib import Path

import sympy as sp


def require(condition, message):
    if not bool(condition):
        raise RuntimeError(message)


def has_asserts(path):
    tree = ast.parse(Path(path).read_text(encoding="utf-8"))
    return any(isinstance(node, ast.Assert) for node in ast.walk(tree))


def reduce_mod(expr, variable, modulus):
    num, den = sp.cancel(expr).as_numer_denom()
    mod_poly = sp.Poly(modulus, variable, domain=sp.QQ)
    den_poly = sp.Poly(den, variable, domain=sp.QQ)
    require(sp.gcd(den_poly, mod_poly).degree() == 0, "nonunit denominator")
    inv = sp.invert(den_poly, mod_poly)
    return sp.rem(sp.Poly(num, variable, domain=sp.QQ) * inv, mod_poly).as_expr()


def reduce_x_coefficients(expr, x, s, modulus):
    poly = sp.Poly(sp.expand(expr), x)
    return all(reduce_mod(coef, s, modulus) == 0 for coef in poly.all_coeffs())


def coprime_to(expr, variable, modulus):
    poly = sp.Poly(sp.cancel(expr).as_numer_denom()[0], variable, domain=sp.QQ)
    mod_poly = sp.Poly(modulus, variable, domain=sp.QQ)
    return sp.gcd(poly, mod_poly).degree() == 0


def recurrence_data(N, x, s):
    p = sp.Rational(N - 1, N) * s - sp.Rational(N - 2, N)
    h = [sp.Integer(1), s]
    for _ in range(2, N):
        h.append(sp.expand(s * h[-1] - p * h[-2]))
    H = sp.expand(N * h[N - 2] - (N - 1) * h[N - 3])
    diagonal = s - sp.Rational(2 * (N - 2), N)
    Q = sp.cancel(H / diagonal)
    require(sp.denom(Q) in sp.QQ, "Q denominator")
    Q = sp.Poly(Q, s, domain=sp.QQ).as_expr()

    E = sp.expand(x**2 - s * x + p)
    D = sp.expand(x ** (N - 1) * (x - 1))
    T = x * (x - 1)
    A = sp.expand(p * (N * h[N - 3] - (N - 1) * h[N - 4]))
    C = sp.expand(-(N - 1) * A)
    B = sp.expand(D + A * (x - p))
    quotient, remainder = sp.div(B, E**2, x)

    return {
        "p": p,
        "h": h,
        "H": H,
        "diagonal": diagonal,
        "Q": Q,
        "E": E,
        "D": D,
        "T": T,
        "A": A,
        "C": C,
        "B": B,
        "S": quotient,
        "remainder": remainder,
    }


def algebra_check(N, x, s):
    data = recurrence_data(N, x, s)
    Q = data["Q"]
    require(sp.degree(Q, s) == N - 3, "Q degree")
    require(sp.expand(data["H"] - data["diagonal"] * Q) == 0, "diagonal factor")
    require(sp.discriminant(Q, s) != 0, "Q squarefree in bounded control")
    require(reduce_x_coefficients(data["remainder"], x, s, Q), "E^2 division")
    require(sp.degree(data["S"], x) == N - 4, "S degree")
    require(
        reduce_x_coefficients(
            sp.diff(data["B"], x) * data["D"]
            - data["B"] * sp.diff(data["D"], x)
            - data["C"] * data["E"] * x ** (N - 2),
            x,
            s,
            Q,
        ),
        "response derivative identity",
    )
    require(
        reduce_x_coefficients(data["B"] - data["D"] - data["A"] * (x - data["p"]), x, s, Q),
        "linear third-fibre defect",
    )

    # Every gate is nonzero on every root of Q in the checked degrees.
    disc_e = sp.discriminant(data["E"], x)
    disc_s = sp.discriminant(data["S"], x) if N > 5 else sp.Integer(1)
    res_se = sp.resultant(data["S"], data["E"] * x * (x - 1), x)
    res_e = sp.resultant(data["E"], x * (x - 1), x)
    gates = [disc_e, data["A"], disc_s, res_se, res_e]
    require(all(coprime_to(gate, s, Q) for gate in gates), "converse gate")
    return data


def compose(left, right):
    return tuple(left[right[i]] for i in range(len(left)))


def inverse(perm):
    out = [0] * len(perm)
    for i, image in enumerate(perm):
        out[image] = i
    return tuple(out)


def cycle_type(perm):
    unseen = set(range(len(perm)))
    lengths = []
    while unseen:
        start = min(unseen)
        cur = start
        length = 0
        while cur in unseen:
            unseen.remove(cur)
            cur = perm[cur]
            length += 1
        lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


def generated_group(generators):
    identity = tuple(range(len(generators[0])))
    group = {identity}
    frontier = [identity]
    while frontier:
        current = frontier.pop()
        for generator in generators:
            nxt = compose(current, generator)
            if nxt not in group:
                group.add(nxt)
                frontier.append(nxt)
    return group


def is_transitive(generators):
    orbit = {0}
    frontier = [0]
    while frontier:
        point = frontier.pop()
        for generator in generators:
            image = generator[point]
            if image not in orbit:
                orbit.add(image)
                frontier.append(image)
    return len(orbit) == len(generators[0])


def normalized_nielsen(N):
    L = N - 1
    infinity = L
    sigma_one = tuple(list(range(1, L)) + [0, infinity])
    records = []
    all_pairs = 0
    for b, c in combinations(range(1, L), 2):
        all_pairs += 1
        sigma_zero = list(range(N))
        sigma_zero[infinity], sigma_zero[0] = 0, infinity
        sigma_zero[b], sigma_zero[c] = c, b
        sigma_zero = tuple(sigma_zero)
        pole = inverse(compose(sigma_zero, sigma_one))
        expected = tuple(sorted((c - b, L + 1 - (c - b)), reverse=True))
        require(cycle_type(pole) == expected, "cycle-length formula")
        if cycle_type(pole) == (L, 1):
            records.append((b, c, sigma_zero, pole))
    require(len(records) == N - 3, "Nielsen orbit count")
    require(all(c == b + 1 for b, c, _, _ in records), "adjacent-chord law")
    return records, all_pairs


def raw_n5_count():
    N = 5
    sigma_one = (1, 2, 3, 0, 4)
    candidates = []
    points = range(N)
    for four in combinations(points, 4):
        a, b, c, d = four
        pairings = [((a, b), (c, d)), ((a, c), (b, d)), ((a, d), (b, c))]
        for pairs in pairings:
            sigma_zero = list(range(N))
            for u, v in pairs:
                sigma_zero[u], sigma_zero[v] = v, u
            sigma_zero = tuple(sigma_zero)
            pole = inverse(compose(sigma_zero, sigma_one))
            if cycle_type(pole) == (4, 1) and is_transitive((sigma_zero, sigma_one)):
                candidates.append(sigma_zero)
    require(len(candidates) == 8, "raw N=5 candidate count")
    return candidates


def n5_algebra(x, a):
    relation = 5 * a**2 + 2 * a + 1
    u = (a - 1) / 2
    v = -(2 * a + 1) / 5
    C = 4 * (a - 2) / 25
    E = x**2 + u * x + v
    S = x - a
    D = x**4 * (x - 1)
    B = sp.expand(S * E**2)
    require(
        reduce_x_coefficients(B - D - (2 - a) * (x - v) / 25, x, a, relation),
        "N=5 third fibre",
    )
    require(
        reduce_x_coefficients(
            sp.diff(B, x) * D - B * sp.diff(D, x) - C * E * x**3,
            x,
            a,
            relation,
        ),
        "N=5 derivative",
    )
    values = [
        sp.discriminant(E, x),
        E.subs(x, 0),
        E.subs(x, 1),
        E.subs(x, a),
        E.subs(x, v),
        C,
    ]
    require(all(coprime_to(value, a, relation) for value in values), "N=5 gates")
    require(sp.discriminant(relation, a) == -16, "Q(i) discriminant")


def main():
    require(not has_asserts(Path(__file__)), "truth-bearing Python assert found")
    x, s, a = sp.symbols("x s a")

    print("THM-2800 TWO-POLE E=2 BITANGENT/NIELSEN AUDIT")
    print("N | Q_N(s) | Nielsen classes")
    algebra_rows = 0
    pair_rows = 0
    for N in range(4, 11):
        data = algebra_check(N, x, s)
        records, all_pairs = normalized_nielsen(N)
        algebra_rows += 1
        pair_rows += all_pairs
        if N <= 7:
            print(f"{N:2d} | {sp.factor(data['Q'])} | {len(records)}")
        else:
            print(f"{N:2d} | degree {sp.degree(data['Q'], s)} squarefree | {len(records)}")

    n4 = normalized_nielsen(4)[0][0]
    n5 = normalized_nielsen(5)[0]
    require(len(generated_group((n4[2], (1, 2, 0, 3)))) == 12, "N=4 A4 order")
    require(all(len(generated_group((row[2], (1, 2, 3, 0, 4)))) == 20 for row in n5), "N=5 group")
    raw_n5_count()
    n5_algebra(x, a)

    print(f"exact algebra degrees N=4..10: {algebra_rows}")
    print(f"normalized pair placements checked: {pair_rows}")
    print("N=4 monodromy order=12 (A4)")
    print("N=5 raw candidates=8, centralizer orbits=2, group order=20")
    print("N=5 field polynomial=5*a^2+2*a+1, discriminant=-16")
    print("bitangent eliminant and response converse gates: PASS")
    print("assert_nodes=0")


if __name__ == "__main__":
    main()
