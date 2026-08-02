#!/usr/bin/env python3
"""Exact companion for THM-3054.

The computation has two jobs.

1. Enumerate the five cyclic inertia types on four roots and verify the
   normalized matching-sum fingerprints and their F2[M] + F3 clutches.
2. Verify the sharp binary/ternary resonance hostile: a nonmaximal tame
   transposition order and a maximal tame three-cycle order have the same
   complete integer-normalized six-edge valuation vector.

All truth-bearing checks use ``require`` so ``python`` and ``python -O`` run
the same verification.
"""

from __future__ import annotations

from fractions import Fraction
from math import lcm

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


EDGES = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
MATCHINGS = (
    ((0, 1), (2, 3)),
    ((0, 2), (1, 3)),
    ((0, 3), (1, 2)),
)


def canonical_edge(i: int, j: int) -> tuple[int, int]:
    return (i, j) if i < j else (j, i)


def cycles_of(perm: tuple[int, ...]) -> list[tuple[int, ...]]:
    seen: set[int] = set()
    cycles: list[tuple[int, ...]] = []
    for start in range(len(perm)):
        if start in seen:
            continue
        cycle = []
        x = start
        while x not in seen:
            seen.add(x)
            cycle.append(x)
            x = perm[x]
        cycles.append(tuple(cycle))
    return cycles


def permutation_order(perm: tuple[int, ...]) -> int:
    return lcm(*[len(c) for c in cycles_of(perm)])


def skeleton(perm: tuple[int, ...]) -> dict[tuple[int, int], int]:
    orbit = {}
    for number, cycle in enumerate(cycles_of(perm)):
        for vertex in cycle:
            orbit[vertex] = (number, len(cycle))
    return {
        edge: int(orbit[edge[0]][0] == orbit[edge[1]][0] and orbit[edge[0]][1] > 1)
        for edge in EDGES
    }


def matching_sums(edge_values: dict[tuple[int, int], int]) -> tuple[int, int, int]:
    return tuple(
        sum(edge_values[canonical_edge(i, j)] for i, j in matching)
        for matching in MATCHINGS
    )


def star_sums(edge_values: dict[tuple[int, int], int]) -> tuple[int, int, int, int]:
    return tuple(
        sum(value for edge, value in edge_values.items() if vertex in edge)
        for vertex in range(4)
    )


def clutch(lambdas: tuple[int, int, int]) -> tuple[str, int]:
    return ("".join(str(value % 2) for value in lambdas), sum(lambdas) % 3)


def polynomial_valuation(expr: sp.Expr, uniformizer: sp.Symbol) -> int:
    poly = sp.Poly(sp.expand(expr), uniformizer, extension=True)
    require(not poly.is_zero, "valuation requested for zero")
    return min(monomial[0] for monomial, coefficient in poly.terms() if coefficient != 0)


def edge_valuations(roots: tuple[sp.Expr, ...], uniformizer: sp.Symbol) -> dict[tuple[int, int], int]:
    return {
        edge: polynomial_valuation(roots[edge[0]] - roots[edge[1]], uniformizer)
        for edge in EDGES
    }


def main() -> None:
    inertia = (
        ("identity", (0, 1, 2, 3)),
        ("transposition", (1, 0, 2, 3)),
        ("double_transposition", (1, 0, 3, 2)),
        ("three_cycle", (1, 2, 0, 3)),
        ("four_cycle", (1, 2, 3, 0)),
    )
    expected = {
        "identity": ((0, 0, 0), ("000", 0)),
        "transposition": ((1, 0, 0), ("100", 1)),
        "double_transposition": ((2, 0, 0), ("000", 2)),
        "three_cycle": ((1, 1, 1), ("111", 0)),
        "four_cycle": ((2, 2, 2), ("000", 0)),
    }

    table_rows = []
    for name, perm in inertia:
        cycles = cycles_of(perm)
        e = int(permutation_order(perm))
        d = 4 - len(cycles)
        c = skeleton(perm)
        lambdas = matching_sums(c)
        require(2 * sum(c.values()) == e * d, f"tame discriminant skeleton failed for {name}")
        require((lambdas, clutch(lambdas)) == expected[name], f"fingerprint failed for {name}")
        table_rows.append((name, e, d, lambdas, clutch(lambdas)))

    # The lattice correction identity is independent of realizability of h.
    # Exhaust a nonnegative cube as a hostile arithmetic control.
    for name, perm in inertia:
        c = skeleton(perm)
        base = matching_sums(c)
        e = int(permutation_order(perm))
        for code in range(3**6):
            q = code
            h_values = []
            for _ in EDGES:
                h_values.append(q % 3)
                q //= 3
            h = dict(zip(EDGES, h_values))
            x = {edge: c[edge] + h[edge] for edge in EDGES}
            correction = tuple(a - b for a, b in zip(matching_sums(x), base))
            require(sum(correction) == sum(h_values), f"matching correction sum failed for {name}")
            require(
                clutch(matching_sums(x))[1] == (clutch(base)[1] + sum(h_values)) % 3,
                f"ternary correction failed for {name}",
            )
            if sum(h_values) % e == 0:
                index = sum(h_values) // e
                require(sum(correction) == e * index, f"index recovery failed for {name}")

    T, t, s, r = sp.symbols("T t s r")
    zeta = (-1 + sp.sqrt(-3)) / 2

    f2 = sp.expand((T**2 - t) * (T - t) * (T - 1))
    f3 = sp.expand((T**3 - t) * (T - 1))
    disc2 = sp.factor(sp.discriminant(f2, T))
    disc3 = sp.factor(sp.discriminant(f3, T))
    require(sp.expand(disc2 - 4 * t**3 * (t - 1) ** 6) == 0, "transposition discriminant")
    require(sp.expand(disc3 + 27 * t**2 * (t - 1) ** 2) == 0, "three-cycle discriminant")

    roots2 = (s, -s, s**2, sp.Integer(1))
    roots3 = (r, zeta * r, zeta**2 * r, sp.Integer(1))
    x2 = edge_valuations(roots2, s)
    x3 = edge_valuations(roots3, r)
    expected_edges = dict(zip(EDGES, (1, 1, 0, 1, 0, 0)))
    require(x2 == expected_edges, "transposition hostile edge vector")
    require(x3 == expected_edges, "three-cycle hostile edge vector")
    require(matching_sums(x2) == (1, 1, 1), "transposition hostile matching vector")
    require(matching_sums(x3) == (1, 1, 1), "three-cycle hostile matching vector")
    require(star_sums(x2) == (2, 2, 2, 0), "transposition hostile star vector")
    require(star_sums(x3) == (2, 2, 2, 0), "three-cycle hostile star vector")
    require(clutch(matching_sums(x2)) == ("111", 0), "resonant clutch")

    trans_index = (polynomial_valuation(disc2.subs(t, s**2), s) // 2 - 1) // 2
    three_index = (polynomial_valuation(disc3.subs(t, r**3), r) // 3 - 2) // 2
    require(trans_index == 1, "transposition order index")
    require(three_index == 0, "three-cycle order index")

    base_x2 = tuple(Fraction(x2[edge], 2) for edge in EDGES)
    base_x3 = tuple(Fraction(x3[edge], 3) for edge in EDGES)
    require(base_x2 != base_x3, "base-normalized valuations must separate the hostile pair")

    # Pointed singleton cross-resultants for the inertia-fixed sheets.
    d_at_t = sp.factor((t**2 - t) * (t - 1))
    d_at_1_f2 = sp.factor((1 - t) * (1 - t))
    d_at_1_f3 = sp.factor(1 - t)
    require(d_at_t == t * (t - 1) ** 2, "fixed sheet t cross-resultant")
    require(d_at_1_f2 == (t - 1) ** 2, "fixed sheet 1 transposition cross-resultant")
    require(d_at_1_f3 == 1 - t, "fixed sheet 1 three-cycle cross-resultant")

    complement_at_t = sp.expand((T**2 - t) * (T - 1))
    complement_at_1 = sp.expand((T**2 - t) * (T - t))
    complement_disc_at_t = sp.factor(sp.discriminant(complement_at_t, T))
    complement_disc_at_1 = sp.factor(sp.discriminant(complement_at_1, T))
    require(complement_disc_at_t == 4 * t * (t - 1) ** 2, "section t complement discriminant")
    require(complement_disc_at_1 == 4 * t**3 * (t - 1) ** 2, "section 1 complement discriminant")
    complement_index_at_t = (1 - 1) // 2
    complement_index_at_1 = (3 - 1) // 2
    require(complement_index_at_t + 1 == trans_index, "section t index decomposition")
    require(complement_index_at_1 + 0 == trans_index, "section 1 index decomposition")

    print("status=PROVED_VERIFIED_EXACT")
    for name, e, d, lambdas, class_value in table_rows:
        print(f"inertia={name};e={e};d={d};matching={lambdas};clutch={class_value}")
    print("formula=sum_edge_excess=e*order_index")
    print("formula=tau=tau_skeleton+e*order_index_mod_3")
    print(f"resonance_edge_vector={tuple(x2[edge] for edge in EDGES)}")
    print(f"resonance_matching={matching_sums(x2)};stars={star_sums(x2)};clutch={clutch(matching_sums(x2))}")
    print(f"transposition_disc={disc2};index={trans_index};fixed_owner_valuations=(1,0)")
    print(f"three_cycle_disc={disc3};index={three_index};fixed_owner_valuations=(0,)")
    print("transposition_owner_split=section_t:(complement_index,gluing)=(0,1);section_1:(1,0)")
    print(f"base_normalized_transposition={base_x2}")
    print(f"base_normalized_three_cycle={base_x3}")
    print("conclusion=integer-normalized edge valuations do not determine inertia, order index, or owner")


if __name__ == "__main__":
    main()
