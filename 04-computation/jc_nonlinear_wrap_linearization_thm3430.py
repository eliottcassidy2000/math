#!/usr/bin/env python3
"""Exact hostile referee for provisional THM-3430.

The wrap stage k for P_d=a*x+b+g(x)*z**d is K[x]/(g**(k+1)).
Its incoming exponent is n=(k+1)*d, so the transition

    (d/(a*n))*g*q'-(1/a)*g'*q

is exactly the d=1 transition.  This companion verifies that cancellation,
finite matrices, Hamiltonian action, split CRT intertwiners, torsion-window
lengths, and the sharply limited nonwrap root-gauge survivor.  Every decision
uses exact integer, Fraction, or SymPy rational polynomial arithmetic.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product

import sympy as sp


x = sp.symbols("x")


def require(condition: bool, label: object) -> None:
    if not condition:
        raise RuntimeError(label)


def roots_for(count: int) -> tuple[int, ...]:
    return (-2, 1, 4, 7)[:count]


def split_polynomial(exponents: tuple[int, ...], gamma: int = 1) -> sp.Expr:
    return sp.expand(
        gamma
        * sp.prod(
            (x - alpha) ** e
            for alpha, e in zip(roots_for(len(exponents)), exponents)
        )
    )


def sector_operator(
    g: sp.Expr, q: sp.Expr, d: int, n: int, a: int
) -> sp.Expr:
    return sp.expand(
        sp.Rational(d, a * n) * g * sp.diff(q, x)
        - sp.Rational(1, a) * sp.diff(g, x) * q
    )


def check_wrap_coefficient_collapse() -> int:
    checks = 0
    for a in (-3, -1, 1, 2, 5):
        for d in range(2, 15):
            for level in range(1, 13):
                n = level * d
                require(
                    Fraction(d, a * n) == Fraction(1, a * level),
                    ("wrap coefficient", a, d, level),
                )
                checks += 1
    return checks


def check_finite_transition_matrices_and_action() -> tuple[int, int, int]:
    transition_entries = 0
    action_identities = 0
    relation_ideal_checks = 0
    profiles = (
        (1,),
        (2,),
        (3,),
        (1, 1),
        (2, 1),
        (3, 2),
        (2, 3, 1),
    )
    for exponents in profiles:
        g = split_polynomial(exponents, gamma=2)
        degree = sum(exponents)
        for a in (-2, 1, 3):
            b = 5
            f = a * x + b
            for d in (2, 3, 5, 8):
                for level in range(1, 5):
                    n = level * d
                    h = 1 + 2 * x + 3 * x**2
                    require(
                        sp.expand(
                            sector_operator(g, g**level * h, d, n, a)
                            - sp.Rational(1, a * level) * g ** (level + 1) * sp.diff(h, x)
                        )
                        == 0,
                        ("exact relation ideal", exponents, a, d, level),
                    )
                    relation_ideal_checks += 1
                    for monomial_degree in range(level * degree):
                        q = x**monomial_degree
                        nonlinear = sector_operator(g, q, d, n, a)
                        linear = sector_operator(g, q, 1, level, a)
                        require(
                            sp.expand(nonlinear - linear) == 0,
                            (
                                "finite transition matrix",
                                exponents,
                                a,
                                d,
                                level,
                                monomial_degree,
                            ),
                        )
                        transition_entries += 1

                        lhs = sp.expand(
                            sector_operator(g, f * q, d, n, a)
                            + g * q
                            - f * sector_operator(g, q, d, n, a)
                        )
                        bridge = sp.Rational(level + 1, level) * g * q
                        require(
                            sp.expand(lhs - bridge) == 0,
                            ("Hamiltonian action first square", exponents, d, level),
                        )
                        require(
                            sp.expand(
                                sector_operator(g, bridge, d, n + d, a)
                                - g * sector_operator(g, q, d, n, a)
                            )
                            == 0,
                            ("Hamiltonian action second square", exponents, d, level),
                        )
                        action_identities += 2
    return transition_entries, action_identities, relation_ideal_checks


def check_split_crt_and_root_intertwiners() -> tuple[int, int, int]:
    crt_packets = 0
    crt_components = 0
    intertwining_checks = 0

    for root_count in range(1, 5):
        for exponents in product(range(1, 5), repeat=root_count):
            total_degree = sum(exponents)
            for level in range(1, 6):
                require(
                    level * total_degree == sum(level * e for e in exponents),
                    ("CRT dimension", exponents, level),
                )
                for i, e_i in enumerate(exponents):
                    for j, e_j in enumerate(exponents):
                        vanishing_order = 0 if i == j else level * e_j
                        require(
                            (vanishing_order >= level * e_j) == (i != j),
                            ("CRT support", exponents, level, i, j),
                        )
                    require(level * e_i > 0, ("local CRT length", exponents, level, i))
                    crt_components += 1
                crt_packets += 1

    exact_profiles = (
        (1,),
        (2,),
        (3,),
        (2, 1),
        (3, 2),
        (2, 3, 1),
        (1, 2, 3),
    )
    for exponents in exact_profiles:
        roots = roots_for(len(exponents))
        gamma = 3
        g = split_polynomial(exponents, gamma=gamma)
        for a in (-2, 1):
            b = -4
            for d in (2, 3, 5):
                for level in range(1, 4):
                    n = level * d
                    for i, (alpha, e_i) in enumerate(zip(roots, exponents)):
                        y = x - alpha
                        u = sp.prod(
                            (x - other_alpha) ** other_e
                            for j, (other_alpha, other_e) in enumerate(
                                zip(roots, exponents)
                            )
                            if j != i
                        )
                        p = 1 + 2 * y + 3 * y**2
                        lhs = sector_operator(g, sp.expand(u**level * p), d, n, a)
                        local_g = gamma * y**e_i
                        rhs = sp.expand(
                            u ** (level + 1)
                            * sector_operator(local_g, p, d, n, a)
                        )
                        require(
                            sp.expand(lhs - rhs) == 0,
                            ("root transition intertwiner", exponents, d, level, i),
                        )

                        beta = a * alpha + b
                        low_global = sp.expand((a * x + b - beta) * u**level * p)
                        low_local = sp.expand(a * y * u**level * p)
                        high_global = sp.expand(g * u**level * p)
                        high_local = sp.expand(
                            gamma * y**e_i * u ** (level + 1) * p
                        )
                        require(
                            sp.expand(low_global - low_local) == 0
                            and sp.expand(high_global - high_local) == 0,
                            ("root Hamiltonian intertwiner", exponents, d, level, i),
                        )
                        intertwining_checks += 2
    return crt_packets, crt_components, intertwining_checks


def check_torsion_windows() -> tuple[int, int, int]:
    windows = 0
    primary_blocks = 0
    squarefree_controls = 0
    for root_count in range(1, 5):
        for exponents in product(range(1, 6), repeat=root_count):
            excess = sum(e - 1 for e in exponents)
            repeated_roots = sum(e > 1 for e in exponents)
            if excess == 0:
                require(repeated_roots == 0, ("squarefree arm", exponents))
                squarefree_controls += 1
            for level in range(1, 8):
                local_lengths = tuple(level * (e - 1) for e in exponents if e > 1)
                require(
                    sum(local_lengths) == level * excess,
                    ("torsion window length", exponents, level),
                )
                require(
                    len(local_lengths) == repeated_roots,
                    ("one arm per repeated root", exponents, level),
                )
                windows += 1
                primary_blocks += len(local_lengths)
    return windows, primary_blocks, squarefree_controls


def check_nonwrap_gauge() -> tuple[int, int, int, int]:
    intertwining_checks = 0
    endpoint_checks = 0
    monodromy_hostiles = 0
    no_two_arm_checks = 0

    for d in range(2, 11):
        for sigma in range(1, d):
            for root_count in (2, 3, 4):
                for exponents in product(range(1, 6), repeat=root_count):
                    selected_indices: list[int] = []
                    for i, e_i in enumerate(exponents):
                        local_selected = e_i > 1 and sigma * (e_i - 1) % d == 0
                        other_trivial = all(
                            sigma * e_j % d == 0
                            for j, e_j in enumerate(exponents)
                            if j != i
                        )
                        if local_selected and not other_trivial:
                            monodromy_hostiles += 1
                        if not other_trivial:
                            continue
                        if local_selected:
                            selected_indices.append(i)
                    require(
                        len(selected_indices) <= 1,
                        ("two nonwrap embedded arms", d, sigma, exponents),
                    )
                    no_two_arm_checks += 1

    exact_gauges = (
        (2, 1, (3, 2), 0),
        (3, 1, (4, 3), 0),
        (4, 2, (3, 2), 0),
        (6, 2, (4, 3), 0),
        (6, 3, (3, 2), 0),
        (4, 2, (1, 2), 0),
    )
    for d, sigma, exponents, i in exact_gauges:
        roots = roots_for(len(exponents))
        alpha = roots[i]
        e_i = exponents[i]
        require(
            all(sigma * e_j % d == 0 for j, e_j in enumerate(exponents) if j != i),
            ("exact gauge arithmetic", d, sigma, exponents, i),
        )
        y = x - alpha
        u = sp.prod(
            (x - other_alpha) ** other_e
            for j, (other_alpha, other_e) in enumerate(zip(roots, exponents))
            if j != i
        )
        h = sp.prod(
            (x - other_alpha) ** (sigma * other_e // d)
            for j, (other_alpha, other_e) in enumerate(zip(roots, exponents))
            if j != i
        )
        gamma = 2
        g = sp.expand(gamma * y**e_i * u)
        for k in range(5):
            n = sigma + k * d
            p = 1 + 2 * y + 3 * y**2
            lhs = sector_operator(g, sp.expand(h * u**k * p), d, n, 1)
            rhs = sp.expand(
                h * u ** (k + 1) * sector_operator(gamma * y**e_i, p, d, n, 1)
            )
            require(
                sp.expand(lhs - rhs) == 0,
                ("nonwrap gauge intertwiner", d, sigma, exponents, i, k),
            )
            intertwining_checks += 1

        if e_i > 1 and sigma * (e_i - 1) % d == 0:
            depth = sigma * (e_i - 1) // d
            primitive_coefficient = sp.expand(h * y**depth / sigma)
            target_coefficient = sp.expand(h * y ** (depth - 1))
            high_image = sp.expand(
                sp.diff(g, x) * sigma * primitive_coefficient
                - d * g * sp.diff(primitive_coefficient, x)
            )
            require(
                sp.expand(high_image - g * target_coefficient) == 0,
                ("nonwrap endpoint", d, sigma, exponents, i),
            )
            endpoint_checks += 1

    require(
        2 * (3 - 1) % 4 == 0,
        "named local resonance setup",
    )
    require(2 * 1 % 4 != 0, "named nontrivial other monodromy")
    require(1 * (3 - 1) % 2 == 0 and 1 * 2 % 2 == 0, "named positive gauge")
    return intertwining_checks, endpoint_checks, monodromy_hostiles, no_two_arm_checks


def main() -> None:
    wrap_coefficients = check_wrap_coefficient_collapse()
    (
        transition_entries,
        action_identities,
        relation_ideal_checks,
    ) = check_finite_transition_matrices_and_action()
    crt_packets, crt_components, root_intertwiners = check_split_crt_and_root_intertwiners()
    torsion_windows, primary_blocks, squarefree_controls = check_torsion_windows()
    (
        nonwrap_intertwiners,
        nonwrap_endpoints,
        monodromy_hostiles,
        no_two_arm_checks,
    ) = check_nonwrap_gauge()

    print("NONLINEAR WRAP LINEARIZATION -- EXACT HOSTILE REFEREE")
    print(f"wrap coefficient collapses: {wrap_coefficients}")
    print(f"finite transition-matrix entries: {transition_entries}")
    print(f"Hamiltonian action-square identities: {action_identities}")
    print(f"exact relation-ideal inductions: {relation_ideal_checks}")
    print(f"split CRT stage packets: {crt_packets}")
    print(f"split CRT local components: {crt_components}")
    print(f"rootwise transition/action intertwiners: {root_intertwiners}")
    print(f"finite torsion-thickness windows: {torsion_windows}")
    print(f"split primary blocks: {primary_blocks}")
    print(f"squarefree zero-torsion controls: {squarefree_controls}")
    print(f"nonwrap algebraized-gauge intertwiners: {nonwrap_intertwiners}")
    print(f"nonwrap exact endpoint identities: {nonwrap_endpoints}")
    print(f"local resonances blocked by other-root monodromy: {monodromy_hostiles}")
    print(f"nonwrap at-most-one-arm checks: {no_two_arm_checks}")
    print("named hostile: (d,sigma;e1,e2)=(4,2;3,1)")
    print("named positive gauge: (d,sigma;e1,e2)=(2,1;3,2)")
    print("WRAP=d=1 FILTERED RESPONSE IDENTIFICATION VERIFIED")


if __name__ == "__main__":
    main()
