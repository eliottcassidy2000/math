#!/usr/bin/env python3
"""Exact companion for THM-3791's moving-root resonant-jet law."""

from __future__ import annotations

import hashlib
import itertools
import json

import sympy as sp


CHECKS = 0


def check(condition: object, label: str) -> None:
    global CHECKS
    if condition is not True and condition != sp.S.true:
        raise AssertionError(label)
    CHECKS += 1


def same(lhs: sp.Expr, rhs: sp.Expr, label: str) -> None:
    check(sp.cancel(sp.expand(lhs - rhs)) == 0, label)


c, b, e, aa = sp.symbols("c b e aa")


def jacobian_poisson(D: sp.Expr, lhs: sp.Expr, rhs: sp.Expr) -> sp.Expr:
    """Jacobian bracket in the oriented variable order (c,b,e)."""

    variables = (c, b, e)
    matrix = sp.Matrix(
        [
            [sp.diff(D, variable) for variable in variables],
            [sp.diff(lhs, variable) for variable in variables],
            [sp.diff(rhs, variable) for variable in variables],
        ]
    )
    return sp.expand(matrix.det())


def coefficient(poly: sp.Expr, exponent: int) -> sp.Expr:
    return sp.Poly(sp.expand(poly), c).coeff_monomial(c**exponent)


# Smooth Poisson packets, exact arm charts, and transition residues for five
# pole orders.  The controls contain lower, resonant, and higher root jets.
for n in range(1, 6):
    beta = [
        c ** (n - 1) + c ** (n + 1),
        5 + 2 * c ** (n - 1) + 2 * c ** (n + 2),
        11 - c ** (n - 1) + 3 * c ** (n + 2),
    ]
    gamma = sp.Integer(2)
    Sigma = sp.expand(gamma * sp.prod(b - root for root in beta))
    D = sp.expand(c**n * e - Sigma)

    same(jacobian_poisson(D, c, b), c**n, f"n={n} bracket c,b")
    same(
        jacobian_poisson(D, b, e),
        n * c ** (n - 1) * e - sp.diff(Sigma, c),
        f"n={n} bracket b,e",
    )
    same(
        jacobian_poisson(D, c, e),
        sp.diff(Sigma, b),
        f"n={n} bracket c,e",
    )
    same(jacobian_poisson(D, D, c), 0, f"n={n} D central c")
    same(jacobian_poisson(D, D, b), 0, f"n={n} D central b")
    same(jacobian_poisson(D, D, e), 0, f"n={n} D central e")

    constants = [sp.expand(root.subs(c, 0)) for root in beta]
    check(len(set(constants)) == 3, f"n={n} distinct constant roots")
    for i, root in enumerate(beta):
        other_product = sp.expand(gamma * sp.prod(b - beta[j] for j in range(3) if j != i))
        b_chart = sp.expand(root + c**n * aa)
        e_chart = sp.expand(aa * other_product.subs(b, b_chart))
        same(D.subs({b: b_chart, e: e_chart}), 0, f"n={n} chart relation i={i}")
        same(sp.diff(b_chart, aa) / c**n, 1, f"n={n} local symplectic i={i}")
        check(
            other_product.subs({c: 0, b: constants[i]}) != 0,
            f"n={n} retained arm coefficient i={i}",
        )
        check(
            sp.diff(D, b).subs({c: 0, b: constants[i]}) != 0,
            f"n={n} smooth arm i={i}",
        )

    resonant = [coefficient(root, n - 1) for root in beta]
    for i in range(3):
        for j in range(i + 1, 3):
            transition = sp.cancel((beta[i] - beta[j]) / c**n)
            same(
                (b - beta[j]) / c**n - (b - beta[i]) / c**n,
                transition,
                f"n={n} transition {i},{j}",
            )
            primitive_difference = sp.expand((beta[j] - beta[i]) / c**n)
            same(
                coefficient(sp.expand((beta[j] - beta[i])), n - 1),
                resonant[j] - resonant[i],
                f"n={n} residue coefficient {i},{j}",
            )
            check(
                sp.expand(primitive_difference).coeff(c, -1) == resonant[j] - resonant[i],
                f"n={n} Laurent residue {i},{j}",
            )

            # Every nonresonant monomial has an explicit rational primitive.
            difference_poly = sp.Poly(sp.expand(beta[j] - beta[i]), c)
            for (degree,), coeff_value in difference_poly.terms():
                if degree == n - 1:
                    continue
                primitive = coeff_value * c ** (degree - n + 1) / (degree - n + 1)
                same(
                    sp.diff(primitive, c),
                    coeff_value * c ** (degree - n),
                    f"n={n} nonresonant primitive {i},{j},deg={degree}",
                )


# Exact truncated-simplex Cech row.  The q=1 row begins at edges, so its
# edge-to-triangle kernel is not quotiented by the vertex image.
for h in range(2, 8):
    edges = list(itertools.combinations(range(h), 2))
    triangles = list(itertools.combinations(range(h), 3))
    edge_index = {edge: index for index, edge in enumerate(edges)}

    delta0 = sp.zeros(len(edges), h)
    for row, (i, j) in enumerate(edges):
        delta0[row, i] = -1
        delta0[row, j] = 1

    delta1 = sp.zeros(len(triangles), len(edges))
    for row, (i, j, k) in enumerate(triangles):
        delta1[row, edge_index[(j, k)]] = 1
        delta1[row, edge_index[(i, k)]] = -1
        delta1[row, edge_index[(i, j)]] = 1

    check(delta0.rank() == h - 1, f"h={h} vertex difference rank")
    check(delta1 * delta0 == sp.zeros(len(triangles), h), f"h={h} Cech square zero")
    check(len(edges) - delta1.rank() == h - 1, f"h={h} edge kernel dimension")

    rho = sp.Matrix([1] + [0] * (h - 1))
    edge_class = delta0 * rho
    check(edge_class != sp.zeros(len(edges), 1), f"h={h} nonconstant resonant class")
    check(delta1 * edge_class == sp.zeros(len(triangles), 1), f"h={h} resonant cocycle")
    check(delta0 * sp.ones(h, 1) == sp.zeros(len(edges), 1), f"h={h} constants vanish")


# Sharp specializations.
fixed_roots = [sp.Integer(0), sp.Integer(2), sp.Integer(5)]
fixed_n1 = [coefficient(root, 0) for root in fixed_roots]
fixed_n3 = [coefficient(root, 2) for root in fixed_roots]
check(len(set(fixed_n1)) == 3, "fixed-root exponent-one obstruction")
check(fixed_n3 == [0, 0, 0], "fixed-root higher-exponent exactness")

s = sp.symbols("s")
common_move = [root + s * c**3 for root in fixed_roots]
check(
    [coefficient(root, 3) for root in common_move] == [s, s, s],
    "common resonant translation invisible",
)

lambda_1, lambda_2 = sp.symbols("lambda_1 lambda_2")
thm3789_roots = [c**2, lambda_1, lambda_2]
check(
    [coefficient(root, 2) for root in thm3789_roots] == [1, 0, 0],
    "THM-3789 moving-axis class",
)


semantic = {
    "atlas": "Ui=A2_(c,ai); all multi-overlaps=D(c)=Gm_x_A1; ai=(b-beta_i(c))/c^n",
    "class": "[omega]=[([c^(n-1)]beta_i(c))_i] in k^h/k*1",
    "controls": "n1_fixed_nonexact;n>=2_fixed_exact;common_resonant_translation_exact;THM3789=(1,0,...,0)",
    "darboux": "nonconstant_resonant_vector_implies_no_polynomial_Darboux_pair",
    "hypercover": "q1_row_starts_at_edges;ker(edge_to_triangle)=vertex_differences=k^h/constants",
    "surface": "c^n e=gamma*product_i(b-beta_i(c)); beta_i(0)_pairwise_distinct",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

print("theorem=THM-3791-moving-root-danielewski-resonant-jet-de-rham-law")
print("surface=c^n*e=gamma*product_i(b-beta_i(c));beta_i(0)_distinct")
print("atlas=Ui=A2_(c,ai);ai=(b-beta_i(c))/c^n;multi_overlap=D(c)")
print("transition=aj-ai=(beta_i-beta_j)c^(-n)")
print("resonance=r_i=coefficient_[c^(n-1)]beta_i(c)")
print("derham=H2=k^h/k*1;[omega]=[r];exact_iff_r_constant")
print("darboux=nonconstant_r_implies_no_pair")
print("controls=n1_fixed_obstructed;n_ge_2_fixed_exact;THM3789_r=(1,0,...,0)")
print(f"CHECKS={CHECKS}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
