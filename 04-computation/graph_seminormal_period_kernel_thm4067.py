#!/usr/bin/env python3
"""Exact symbolic controls for THM-4067.

This companion checks the graph-seminormal dimension ledger, contact-thick
and ordinary-multiple-point defects, the THM-3696 missing derivative class,
and the fixed/mixed formal-contact obstructions in THM-4063's figure eight.
The arbitrary-degree proofs are in the theorem; finite tables are hostile
controls rather than extrapolations.
"""

from math import comb

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def a_valuation(polynomial, A, c):
    polynomial = sp.Poly(sp.expand(polynomial), A, c, domain=sp.QQ)
    if polynomial.is_zero:
        return None
    return min(monomial[0] for monomial, _ in polynomial.terms())


print("THM4067_PRIMARY_EXACT")

# Abstract edge-Poincare dimension ledger.  Functions have c-degree <=N+1,
# densities degree <=N.  Vertex equalities contribute 2E-V independent
# conditions and beta=E-V+1 cycle periods are independent.
print("abstract_graph_equalizer_ledger")
for name, vertices, edges in (
    ("tree", 4, 3),
    ("triangle", 3, 3),
    ("figure_eight", 5, 6),
):
    beta = edges - vertices + 1
    for cutoff in range(5):
        dim_functions = edges * (cutoff + 2)
        vertex_conditions = 2 * edges - vertices
        dim_equalizer = dim_functions - vertex_conditions
        derivative_rank = dim_equalizer - 1
        dim_densities = edges * (cutoff + 1)
        period_kernel = dim_densities - beta
        require(derivative_rank == period_kernel, (name, cutoff))
    print(
        f"{name}:V={vertices};E={edges};beta={beta};"
        "cutoffs_0_to_4=period_complete"
    )

# Contact-thick two-branch normalization and ordinary plane r-fold lengths.
print("conductor_defect_ledgers")
for contact in range(1, 9):
    delta = contact
    residue_gluing = 1
    defect = delta - residue_gluing
    require(defect == contact - 1, ("contact", contact))
    print(
        f"contact_m={contact}:delta={delta};residue_gluing=1;"
        f"seminormal_defect={defect}"
    )
for branches in range(2, 8):
    delta = comb(branches, 2)
    defect = delta - (branches - 1)
    require(defect == comb(branches - 1, 2), ("multiple point", branches))
    print(
        f"ordinary_plane_r={branches}:delta={delta};"
        f"seminormal_defect={defect}"
    )

# THM-3696: f=b^2 X is in the value equalizer and violates the sole jet law.
b = sp.symbols("b")
X = b * (1 - b**2)
triple_hostile = sp.expand(b**2 * X)
values = tuple(triple_hostile.subs(b, point) for point in (-1, 0, 1))
derivative = sp.diff(triple_hostile, b)
jet_law = (
    derivative.subs(b, -1)
    + 4 * derivative.subs(b, 0)
    + derivative.subs(b, 1)
)
require(values == (0, 0, 0), "THM-3696 value equalizer")
require(jet_law == -4, "THM-3696 missing derivative class")
print(f"thm3696_hostile=b^2X;values={values};jet_law={jet_law}")

# Reconstruct the six scaled supporting lines in U=u/epsilon, C=c/epsilon.
C, U = sp.symbols("C U")
origin = (sp.Rational(0), sp.Rational(0))
B1 = (sp.Rational(1), sp.Rational(0))
B2 = (sp.Rational(2), sp.Rational(1))
D1 = (sp.Rational(-1), sp.Rational(2))
D2 = (sp.Rational(-3), sp.Rational(1))
edge_vertices = (
    (origin, B1),
    (B1, B2),
    (B2, origin),
    (origin, D1),
    (D1, D2),
    (D2, origin),
)


def supporting_line(start, end):
    slope = sp.cancel((end[1] - start[1]) / (end[0] - start[0]))
    intercept = sp.cancel(start[1] - slope * start[0])
    return slope, intercept


lines = tuple(supporting_line(*edge) for edge in edge_vertices)
expected_lines = (
    (sp.Rational(0), sp.Rational(0)),
    (sp.Rational(1), sp.Rational(-1)),
    (sp.Rational(1, 2), sp.Rational(0)),
    (sp.Rational(-2), sp.Rational(0)),
    (sp.Rational(1, 2), sp.Rational(5, 2)),
    (sp.Rational(-1, 3), sp.Rational(0)),
)
require(lines == expected_lines, "six supporting lines")
print("figure_eight_scaled_lines=" + repr(lines))
require(lines[2][0] == lines[4][0], "parallel slopes")
require(lines[4][1] - lines[2][1] == sp.Rational(5, 2), "separation")
print("parallel_pair=e3,e5;separation=(5/2)*epsilon")

# Minimal graph-continuous primitive and its oriented cycle periods.
A, c, u = sp.symbols("A c u")
epsilon = sp.symbols("epsilon")
density = (0, 2, 1, 0, 0, 0)
c_increments = (epsilon, epsilon, -2 * epsilon, -epsilon, -2 * epsilon, 3 * epsilon)
period_1 = sum(density[index] * c_increments[index] for index in range(3))
period_2 = sum(density[index] * c_increments[index] for index in range(3, 6))
require(period_1 == 0 and period_2 == 0, "minimal witness periods")
primitive_endpoint_units = ((0, 0), (0, 2), (2, 0), (0, 0), (0, 0), (0, 0))
print(
    "minimal_witness=h=(0,2(c-epsilon),c,0,0,0);"
    f"density={density};periods=({period_1},{period_2});"
    f"endpoint_units={primitive_endpoint_units}"
)

# Common-target congruences.  On e3 and e5 the c derivative is the same
# tangential operator evaluated at u-values separated by (5/2)A^q.
fixed_checks = 0
mixed_checks = 0
surjection_checks = 0
for opening in range(1, 7):
    epsilon_q = A**opening
    for total_degree in range(7):
        for a_degree in range(total_degree + 1):
            for u_degree in range(total_degree - a_degree + 1):
                c_degree = total_degree - a_degree - u_degree
                target = A**a_degree * c**c_degree * u**u_degree
                restriction_3 = sp.expand(target.subs(u, c / 2))
                restriction_5 = sp.expand(
                    target.subs(u, c / 2 + sp.Rational(5, 2) * epsilon_q)
                )
                fixed_difference = sp.diff(restriction_3, c) - sp.diff(
                    restriction_5, c
                )
                fixed_valuation = a_valuation(fixed_difference, A, c)
                require(
                    fixed_valuation is None or fixed_valuation >= opening,
                    ("fixed contact valuation", opening, target),
                )
                fixed_checks += 1

                mixed_difference = sp.diff(restriction_3, A) - sp.diff(
                    restriction_5, A
                )
                mixed_valuation = a_valuation(mixed_difference, A, c)
                require(
                    mixed_valuation is None or mixed_valuation >= opening - 1,
                    ("mixed contact valuation", opening, target),
                )
                mixed_checks += 1

    # The fixed-kernel quotient is onto: phi=A^i c^n is realized by the
    # graph-continuous primitive H on e3 and a linear connector on e2.
    for a_degree in range(opening):
        for c_degree in range(6):
            phi = A**a_degree * c**c_degree
            H = sp.integrate(phi, (c, 0, c))
            value = sp.expand(H.subs(c, 2 * epsilon_q))
            connector = sp.cancel((value / epsilon_q) * (c - epsilon_q))
            require(not connector.has(sp.zoo), "connector regularity")
            require(connector.subs(c, epsilon_q) == 0, "connector first endpoint")
            require(
                sp.expand(connector.subs(c, 2 * epsilon_q) - value) == 0,
                "connector second endpoint",
            )
            require(sp.diff(H, c) == phi, "surjection primitive")
            surjection_checks += 1

# The A derivative loses exactly one contact order: target f=u realizes it.
for opening in range(1, 7):
    restriction_3 = c / 2
    restriction_5 = c / 2 + sp.Rational(5, 2) * A**opening
    sharp_difference = sp.diff(restriction_3, A) - sp.diff(restriction_5, A)
    require(
        sharp_difference == -sp.Rational(5, 2) * opening * A ** (opening - 1),
        ("sharp mixed loss", opening),
    )

print(
    f"target_contact_checks=fixed:{fixed_checks};mixed:{mixed_checks};"
    f"surjection:{surjection_checks}"
)
print("fixed_period_defect_surjects=R[[c]]/(A^q);all_q>=1;infinite_k_dimension")
print(
    "mixed_cokernel_surjects=R[[c]]/(A^(q-1));"
    "q>=2_nonzero_infinite_k_dimension;q=1_quotient_zero"
)
print("THM4063_displayed_figure_eight_period_completeness=REFUTED")
print("abstract_connection_and_ramification_theorems=SURVIVE")

