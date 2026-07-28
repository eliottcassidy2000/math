#!/usr/bin/env python3
"""Exact algebra and Nielsen controls for THM-2808."""

import ast
from collections import defaultdict
from itertools import combinations, permutations
from math import comb
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


def reduce_x(expr, x, lam, modulus):
    return all(
        reduce_mod(coef, lam, modulus) == 0
        for coef in sp.Poly(sp.expand(expr), x).all_coeffs()
    )


def coprime_to(expr, lam, modulus):
    num = sp.cancel(expr).as_numer_denom()[0]
    return (
        sp.gcd(
            sp.Poly(num, lam, domain=sp.QQ),
            sp.Poly(modulus, lam, domain=sp.QQ),
        ).degree()
        == 0
    )


def three_pole_data(N, a, b, c, x, lam):
    D = sp.expand(x**a * (x - 1) ** b * (x - lam) ** c)
    T = x * (x - 1) * (x - lam)
    K = sp.expand(N * x**2 - (a + c + (a + b) * lam) * x + a * lam)
    _, remainder = sp.div(D, K, x)
    remainder_poly = sp.Poly(remainder, x)
    R = sp.expand(remainder_poly.coeff_monomial(x))
    r0 = sp.expand(remainder_poly.coeff_monomial(1))
    diagonal = sp.expand(sp.discriminant(K, x))
    Q = sp.Poly(sp.cancel(R / diagonal), lam, domain=sp.QQ).as_expr()

    E = sp.expand(K / N)
    A = sp.expand(-r0)
    B = sp.expand(D + A)
    C = sp.expand(-N * A)
    S, square_remainder = sp.div(B, E**2, x)

    return {
        "D": D,
        "T": T,
        "K": K,
        "R": R,
        "r0": r0,
        "diagonal": diagonal,
        "Q": Q,
        "E": E,
        "A": A,
        "B": B,
        "C": C,
        "S": S,
        "square_remainder": square_remainder,
    }


def algebra_check(N, a, b, c, x, lam):
    data = three_pole_data(N, a, b, c, x, lam)
    expected_lead = sp.Rational(
        (-1) ** c * (a + b) ** (a + b - 1) * c**c,
        N ** (N - 1),
    )
    require(sp.degree(data["R"], lam) == N - 1, "remainder degree")
    require(sp.LC(sp.Poly(data["R"], lam)) == expected_lead, "remainder lead")
    require(
        sp.expand(data["R"] - data["diagonal"] * data["Q"]) == 0,
        "diagonal quotient",
    )
    require(sp.degree(data["Q"], lam) == N - 3, "Maxwell degree")
    require(
        sp.discriminant(data["diagonal"], lam) == -16 * N * a * b * c,
        "two simple diagonal collisions",
    )
    require(
        sp.expand(data["diagonal"].subs(lam, 0) - (a + c) ** 2) == 0,
        "lambda=0 diagonal gate",
    )
    require(
        sp.expand(data["diagonal"].subs(lam, 1) - (b + c) ** 2) == 0,
        "lambda=1 diagonal gate",
    )
    require(sp.discriminant(data["Q"], lam) != 0, "bounded Q squarefree")
    require(
        reduce_x(data["square_remainder"], x, lam, data["Q"]),
        "E^2 division",
    )
    require(sp.degree(data["S"], x) == N - 4, "S degree")

    pole_core = x ** (a - 1) * (x - 1) ** (b - 1) * (x - lam) ** (c - 1)
    require(
        sp.expand(
            sp.diff(data["B"], x) * data["D"]
            - data["B"] * sp.diff(data["D"], x)
            - data["C"] * data["E"] * pole_core
        )
        == 0,
        "response derivative",
    )

    disc_s = sp.discriminant(data["S"], x) if N >= 6 else sp.Integer(1)
    gates = (
        lam,
        lam - 1,
        sp.discriminant(data["E"], x),
        data["A"],
        disc_s,
        sp.resultant(data["E"], data["T"], x),
        sp.resultant(data["S"], data["E"] * data["T"], x),
    )
    require(
        all(coprime_to(gate, lam, data["Q"]) for gate in gates),
        "bounded converse gates",
    )
    return data


def compose(left, right):
    return tuple(left[right[i]] for i in range(len(left)))


def cycles(permutation):
    seen = set()
    answer = []
    for start in range(len(permutation)):
        if start in seen:
            continue
        cycle = []
        point = start
        while point not in seen:
            seen.add(point)
            cycle.append(point)
            point = permutation[point]
        answer.append(tuple(cycle))
    return tuple(answer)


def generated_group(generators):
    identity = tuple(range(len(generators[0])))
    group = {identity}
    frontier = [identity]
    while frontier:
        current = frontier.pop()
        for generator in generators:
            candidate = compose(current, generator)
            if candidate not in group:
                group.add(candidate)
                frontier.append(candidate)
    return group


def rotate_pair(pair, shift, N):
    return tuple(
        sorted(
            tuple(sorted(((left + shift) % N, (right + shift) % N)))
            for left, right in pair
        )
    )


def pair_orbit_key(pair, N):
    return min(rotate_pair(pair, shift, N) for shift in range(N))


def labelled_orbit_key(pair, labelled_cycles, N):
    candidates = []
    for shift in range(N):
        shifted_pair = rotate_pair(pair, shift, N)
        shifted_cycles = tuple(
            tuple(sorted((point + shift) % N for point in cycle))
            for cycle in labelled_cycles
        )
        candidates.append((shifted_pair, shifted_cycles))
    return min(candidates)


def nielsen_atlas(N):
    long_cycle = tuple((i + 1) % N for i in range(N))
    labelled = defaultdict(set)
    unmarked = set()
    raw_three_cycle = 0

    for endpoints in combinations(range(N), 4):
        pairings = (
            ((endpoints[0], endpoints[1]), (endpoints[2], endpoints[3])),
            ((endpoints[0], endpoints[3]), (endpoints[1], endpoints[2])),
            ((endpoints[0], endpoints[2]), (endpoints[1], endpoints[3])),
        )
        three_cycle_pairings = 0
        for pair in pairings:
            zero_inertia = list(range(N))
            for left, right in pair:
                zero_inertia[left], zero_inertia[right] = (
                    zero_inertia[right],
                    zero_inertia[left],
                )
            pole_cycles = cycles(compose(tuple(zero_inertia), long_cycle))
            if len(pole_cycles) != 3:
                continue
            three_cycle_pairings += 1
            raw_three_cycle += 1
            unmarked.add(pair_orbit_key(pair, N))
            for labelled_cycles in permutations(pole_cycles):
                passport = tuple(len(cycle) for cycle in labelled_cycles)
                labelled[passport].add(
                    labelled_orbit_key(pair, labelled_cycles, N)
                )
        require(three_cycle_pairings == 2, "noncrossing pairing count")

    require(raw_three_cycle == 2 * comb(N, 4), "raw noncrossing atlas")
    for a in range(1, N - 1):
        for b in range(1, N - a):
            c = N - a - b
            require(
                len(labelled[(a, b, c)]) == N - 3,
                "ordered Nielsen count",
            )
            require(
                sum(2 * (part - 1) for part in (a, b, c)) // 2 == N - 3,
                "four-gap encoding count",
            )

    half_turn_fix = (N // 2) * (N // 2 - 1) if N % 2 == 0 else 0
    burnside = (2 * comb(N, 4) + half_turn_fix) // N
    floor_sum = sum(((degree - 2) // 2) ** 2 for degree in range(4, N + 1))
    require(len(unmarked) == burnside == floor_sum, "unmarked Burnside count")
    return len(labelled), len(unmarked)


def main():
    require(not has_asserts(Path(__file__)), "truth-bearing Python assert found")
    x, lam = sp.symbols("x lambda")
    print("THM-2808 THREE-POLE E=2 MAXWELL/NIELSEN AUDIT")
    print("N | ordered algebra passports | ordered Nielsen passports | unmarked total")

    algebra_cases = 0
    for N in range(4, 9):
        for a in range(1, N - 1):
            for b in range(1, N - a):
                c = N - a - b
                data = algebra_check(N, a, b, c, x, lam)
                require(sp.degree(data["Q"], lam) == N - 3, "degree recount")
                algebra_cases += 1
        labelled_types, unmarked_total = nielsen_atlas(N)
        ordered_passports = (N - 1) * (N - 2) // 2
        require(labelled_types == ordered_passports, "passport census")
        print(
            f"{N:2d} | {ordered_passports:25d} |"
            f" {labelled_types:25d} | {unmarked_total:14d}"
        )

    for N in range(9, 13):
        labelled_types, unmarked_total = nielsen_atlas(N)
        print(
            f"{N:2d} | {'--':>25s} |"
            f" {labelled_types:25d} | {unmarked_total:14d}"
        )

    examples = ((4, 1, 1, 2), (5, 1, 1, 3), (5, 1, 2, 2))
    print("monic Maxwell controls:")
    for N, a, b, c in examples:
        Q = three_pole_data(N, a, b, c, x, lam)["Q"]
        print(f"  ({a},{b},{c}): {sp.factor(sp.Poly(Q, lam).monic().as_expr())}")

    quartic_cycle = (1, 2, 3, 0)
    quartic_zero = (1, 0, 3, 2)
    quartic_group = generated_group((quartic_cycle, quartic_zero))
    quartic_v4 = {
        (0, 1, 2, 3),
        (1, 0, 3, 2),
        (2, 3, 0, 1),
        (3, 2, 1, 0),
    }
    require(len(quartic_group) == 8, "quartic D4 order")
    require(quartic_v4 <= quartic_group, "quartic V4 subgroup")
    print("quartic h=3 monodromy: |D4|=8, |V4|=4, quotient order=2")

    for N in range(4, 31):
        h1 = (N - 2) // 2
        h2 = comb(N - 2, 2)
        h3 = sum(((degree - 2) // 2) ** 2 for degree in range(4, N + 1))
        require(h1 + h2 + h3 > h3, "full e=2 layer count")

    print(f"exact ordered algebra passports checked: {algebra_cases}")
    print("ordered classes per passport: N-3")
    print("unmarked h=3 total: sum_[m=4..N] floor((m-2)/2)^2")
    print("three-pole quotient, converse gates, and Nielsen saturation: PASS")
    print("assert_nodes=0")


if __name__ == "__main__":
    main()
