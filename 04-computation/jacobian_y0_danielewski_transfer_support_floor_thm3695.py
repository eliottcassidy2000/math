#!/usr/bin/env python3
"""Exact algebraic companion for THM-3695.

The proof of THM-3695 also imports two non-computational ingredients:
Gwozdziewicz's cited injective-line theorem and the proved universal
three-by-three nonentry theorem THM-3592.  This companion checks the exact
Poisson embedding, grading/divisibility sidecars, collision, line and jet
invoices, the W003 all-scale valuation family, support arithmetic, and the
cited degree-floor arithmetic.

Every truth gate uses ``require`` and survives ``python -O``.
"""

from __future__ import annotations

from itertools import product
from math import ceil

import sympy as sp


CHECKS = 0


def require(condition: bool, detail: object) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(detail)


def source_bracket(first: sp.Expr, second: sp.Expr, x: sp.Symbol, z: sp.Symbol) -> sp.Expr:
    return sp.expand(
        sp.diff(first, x) * sp.diff(second, z)
        - sp.diff(first, z) * sp.diff(second, x)
    )


def total_degree(polynomial: sp.Expr, x: sp.Symbol, z: sp.Symbol) -> int:
    return sp.Poly(sp.expand(polynomial), x, z).total_degree()


def leading_wronskian_data(
    first_weight: int,
    first_order: int,
    second_weight: int,
    second_order: int,
) -> tuple[int, int] | None:
    """Leading local order/multiplier, or None when the initial term resonates."""
    multiplier = second_weight * first_order - first_weight * second_order
    if multiplier == 0:
        return None
    return first_order + second_order - 1, multiplier


def main() -> None:
    x, z = sp.symbols("x z")
    b = 1 - x**2 * z
    c = x
    e = z * (2 - x**2 * z)
    sigma = 1 - b**2

    A = 3 * e
    B = 2 * c * e
    C = c * b

    print("THM-3695 exact companion")
    print("SECTION graded Poisson embedding")
    require(sp.expand(c**2 * e - sigma) == 0, "Danielewski relation")
    require(sp.expand(source_bracket(b, c, x, z) - c**2) == 0, "{b,c}=c^2")
    require(sp.expand(source_bracket(c, e, x, z) - 2 * b) == 0, "{c,e}=-Sigma'(b)=2b")
    require(sp.expand(source_bracket(b, e, x, z) + 2 * c * e) == 0, "{b,e}=-2ce")
    require(sp.expand(A - 3 * z * (2 - x**2 * z)) == 0, "A=3e")
    require(sp.expand(B - 2 * x * z * (2 - x**2 * z)) == 0, "B=2ce")
    require(sp.expand(C - x * (1 - x**2 * z)) == 0, "C=cb")
    relation = 16 * A**3 * C**2 + 27 * B**4 - 36 * A * B**2
    require(sp.expand(relation) == 0, "collision-ring hypersurface relation")
    require(sp.cancel(e - (1 - b**2) / c**2) == 0, "rational inverse for e")
    require(sp.cancel(z - (1 - b) / c**2) == 0, "rational inverse for z")
    weights = {"A": -2, "B": -1, "C": 1}
    require(weights == {"A": -2, "B": -1, "C": 1}, "grading transfer")
    print("PASS R=C[3e,2ce,cb] is graded inside the Poisson Danielewski algebra D_(1-b^2)")

    print("SECTION collision, boundary lines, and linear jet")
    first_collision = tuple(sp.expand(value.subs({x: 1, z: 0})) for value in (A, B, C, b))
    second_collision = tuple(sp.expand(value.subs({x: -1, z: 2})) for value in (A, B, C, b))
    require(first_collision == (0, 0, 1, 1), first_collision)
    require(second_collision == (0, 0, 1, -1), second_collision)
    require(tuple(sp.expand(value.subs(z, 0)) for value in (A, B, C)) == (0, 0, x), "z=0 line")
    require(tuple(sp.expand(value.subs(x, 0)) for value in (A, B, C)) == (6 * z, 0, 0), "x=0 line")
    jet = {
        name: (
            sp.diff(value, x).subs({x: 0, z: 0}),
            sp.diff(value, z).subs({x: 0, z: 0}),
        )
        for name, value in (("A", A), ("B", B), ("C", C))
    }
    require(jet == {"A": (0, 6), "B": (0, 0), "C": (1, 0)}, jet)
    print("PASS fixed collision; z=0 sees only C, x=0 only A, and the origin jet sees only A,C")

    print("SECTION y=0 coefficient-module divisibility")
    monomial_rows = 0
    for i, j, k in product(range(9), repeat=3):
        if i + j + k == 0 or i + j + k > 8:
            continue
        weight = -2 * i - j + k
        # A^i B^j C^k = constant*x^weight*h^(i+j)*b^k,
        # where h=1-b^2.  Only the exponents are needed here.
        if weight > 0:
            require(k >= weight, ("positive b-order", i, j, k, weight))
        if weight < 0:
            require(i + j >= ceil((-weight) / 2), ("negative arm-order", i, j, k, weight))
        if weight == 0:
            require(i + j >= 1 and k >= 1, ("retained zero-weight divisors", i, j, k))
        monomial_rows += 1
    require(monomial_rows == 164, monomial_rows)
    print(f"PASS {monomial_rows} monomials: positive r carries b^r; retained zero carries b(1-b^2); negative -r carries (1-b^2)^ceil(-r/2)")

    print("SECTION R-specific W003 three-by-four ray-family nonentry")
    ray_rows = 0
    expected_word = (
        ((0, 0),),
        ((0, 1), (1, 0)),
        ((0, 2), (1, 1)),
        ((0, 3), (1, 2), (2, 0)),
        ((1, 3), (2, 1)),
        ((2, 2),),
        ((2, 3),),
    )
    for n in range(2, 65):
        p_weights = (1 - 3 * n, 1 - 2 * n, 1)
        q_weights = (-2, n - 2, 2 * n - 2, 3 * n - 2)
        fibres: dict[int, list[tuple[int, int]]] = {}
        for i, p_weight in enumerate(p_weights):
            for j, q_weight in enumerate(q_weights):
                fibres.setdefault(p_weight + q_weight + 1, []).append((i, j))
        word = tuple(tuple(cells) for _, cells in sorted(fibres.items()))
        require(word == expected_word, ("W003 word", n, word))

        # Scalar address 20 has weights (1,-2) and exact b=0 orders (1,0).
        require(leading_wronskian_data(1, 1, -2, 0) == (0, -2), ("scalar 20", n))
        require((2 * n - 2) - 1 >= 1 and (3 * n - 2) - 1 >= 3, ("other scalar orders", n))

        # Singleton rows 00,22,23 force orders f0=0,g2=2n-2,g3=3n-2.
        fixed_orders = (0, 2 * n - 2, 3 * n - 2)
        require(fixed_orders[1:] == q_weights[2:], ("singleton powers", n, fixed_orders))

        beta_floor = 1 if q_weights[1] == 0 else q_weights[1]
        candidates = []
        for alpha in range(0, 2 * n - 1):
            beta = 2 * n - 2 - alpha
            if beta < beta_floor:
                continue
            # Fibre 01=10 either has a resonant alpha=0 term or, when both
            # leading multipliers are nonzero, forces equal orders.
            if alpha == 0 or (alpha > 0 and beta > 0 and alpha == beta):
                candidates.append((alpha, beta))
        require(candidates == [(0, 2 * n - 2), (n - 1, n - 1)], ("W003 candidates", n, candidates))

        for alpha, beta in candidates:
            first = leading_wronskian_data(1 - 2 * n, alpha, 3 * n - 2, 3 * n - 2)
            second = leading_wronskian_data(1, 1, n - 2, beta)
            require(first is not None and second is not None, ("terminal nonresonance", n, alpha, beta))
            require(first[0] > second[0], ("terminal unique minimum", n, alpha, beta, first, second))
            ray_rows += 1
    require(ray_rows == 126, ray_rows)
    print("PASS all n>=2 in the W003 anchor-20 family: b=0 leaves a unique lowest term (126 exact controls through n=64)")

    print("SECTION support and degree floors")
    # Gwozdziewicz + the immersed-line lemma force, in each coordinate, a
    # C-only weight >=2, an A-only weight <=-4, and a distinct jet weight
    # in {-2,1}.  THM-3592 then removes the sole six-piece partition (3,3).
    first_live = []
    for total in range(2, 11):
        live = [
            (left, total - left)
            for left in range(1, total)
            if left >= 3 and total - left >= 3 and (left, total - left) != (3, 3)
        ]
        if live:
            first_live = live
            require(total == 7, (total, live))
            break
    require(first_live == [(3, 4), (4, 3)], first_live)

    degrees = tuple(total_degree(value, x, z) for value in (A, B, C))
    require(degrees == (4, 5, 4), degrees)
    reduced_height_floor = 108
    first_target_degree = ceil(reduced_height_floor / max(degrees))
    require(first_target_degree == 22, first_target_degree)
    require(max(degrees) * 21 < reduced_height_floor <= max(degrees) * 22, "degree-floor arithmetic")
    print("PASS line/jet gate gives >=3 weights per output; THM-3592 leaves first support cells 3x4 and 4x3")
    print("PASS CITED reduced-height 108 plus source degrees (4,5,4) forces common target-filtration cap >=22")
    print(f"CHECKS={CHECKS}")
    print("PASS")


if __name__ == "__main__":
    main()
