#!/usr/bin/env python3
"""Exact hostile referee for provisional THM-3427.

On the normalized generic fiber P=x+g(x)z**d=t, sector sigma has

    D(q z**sigma) = z**(sigma-1) L_sigma(q),
    L_sigma(q) = sigma*q + (t-x)*(sigma*(g'/g)*q-d*q').

For sigma<d this is the direct weight-sigma input.  For sigma=d the actual
input is weight zero; evaluation at x=t splits off its constant kernel and
shows that the same L_d image is valid.  The candidate theorem classifies
the constant target observer [z**(sigma-1)] and gives a polynomial sector
presentation.  At an accessible leading-degree resonance it also verifies
the normalized horizontal-section truncation and the complete root-value
interpolation formula for the low defect.

This companion uses SymPy only for exact Q(t)-linear systems and rational
identities.  Wrap residues and orbit-arrow checks use Fraction arithmetic.
No floating-point decision occurs.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product
from math import comb, factorial

import sympy as sp


x, t = sp.symbols("x t")


def require(condition: bool, label: object) -> None:
    if not condition:
        raise RuntimeError(label)


def multiplied_equation(g: sp.Expr, q: sp.Expr, d: int, sigma: int) -> sp.Expr:
    """Numerator of L_sigma(q)-1 after multiplication by g."""
    return sp.expand(
        sigma * g * q
        + (t - x) * (sigma * sp.diff(g, x) * q - d * g * sp.diff(q, x))
        - g
    )


def candidate_degree(d: int, sigma: int, exponents: tuple[int, ...]) -> int | None:
    numerator = sigma * (sum(exponents) - 1)
    if numerator % d:
        return None
    degree = numerator // d
    return degree if degree >= len(exponents) else None


def predicted_constant_exact(d: int, sigma: int, exponents: tuple[int, ...]) -> bool:
    return (
        len(exponents) == 1
        and exponents[0] > 1
        and sigma * (exponents[0] - 1) % d == 0
    )


def solve_nonwrap_profile(
    d: int, sigma: int, exponents: tuple[int, ...]
) -> tuple[bool, sp.Expr | None]:
    degree = candidate_degree(d, sigma, exponents)
    if degree is None:
        return False, None
    roots = tuple(2 * j - 1 for j in range(len(exponents)))
    g = sp.prod((x - alpha) ** e for alpha, e in zip(roots, exponents))
    coefficients = sp.symbols(f"c0:{degree + 1}")
    q = sum(coefficients[j] * x**j for j in range(degree + 1))
    equation = sp.Poly(multiplied_equation(g, q, d, sigma), x)
    solution_set = sp.linsolve(equation.all_coeffs(), coefficients)
    if solution_set == sp.EmptySet:
        return False, None
    solution = next(iter(solution_set), None)
    if solution is None:
        return False, None
    q_solution = sp.factor(q.subs(dict(zip(coefficients, solution))))
    require(
        sp.cancel(multiplied_equation(g, q_solution, d, sigma)) == 0,
        ("direct substitution", d, sigma, exponents),
    )
    return True, q_solution


def one_root_q(d: int, sigma: int, e: int, alpha: int) -> sp.Expr:
    """Closed coefficient q with L_sigma(q)=1 in a selected one-root sector."""
    require(e > 1 and sigma * (e - 1) % d == 0, ("not selected", d, sigma, e))
    depth = sigma * (e - 1) // d
    y = x - alpha
    tau = t - alpha
    out = 0
    for n in range(1, depth + 1):
        numerator = d ** (n - 1) * factorial(depth - 1) // factorial(depth - n)
        denominator = 1
        for j in range(depth - n, depth):
            denominator *= d * j + sigma
        out += sp.Rational(numerator, denominator) * (y / tau) ** n
    return sp.factor(out)


def check_closed_nonwrap_primitives() -> int:
    checks = 0
    for d in range(2, 10):
        for sigma in range(1, d):
            for e in range(2, 15):
                if sigma * (e - 1) % d:
                    continue
                for alpha in (-2, 3):
                    y = x - alpha
                    q = one_root_q(d, sigma, e, alpha)
                    reduced = sp.cancel(
                        sigma * q
                        + (t - x) * (sigma * e * q / y - d * sp.diff(q, x))
                        - 1
                    )
                    require(reduced == 0, ("closed nonwrap primitive", d, sigma, e, alpha))
                    checks += 1
    return checks


def check_nonwrap_direct_universe(
) -> tuple[int, int, int, list[tuple[int, int, tuple[int, ...], sp.Expr]]]:
    profiles = 0
    solved_systems = 0
    pruned = 0
    solutions: list[tuple[int, int, tuple[int, ...], sp.Expr]] = []
    for d in range(2, 7):
        for sigma in range(1, d):
            for root_count in range(1, 4):
                for exponents in product(range(1, 4), repeat=root_count):
                    profiles += 1
                    if candidate_degree(d, sigma, exponents) is None:
                        exists, q_solution = False, None
                        pruned += 1
                    else:
                        solved_systems += 1
                        exists, q_solution = solve_nonwrap_profile(d, sigma, exponents)
                    require(
                        exists == predicted_constant_exact(d, sigma, exponents),
                        ("nonwrap profile", d, sigma, exponents, q_solution),
                    )
                    if exists:
                        require(q_solution is not None, ("missing solution", d, sigma, exponents))
                        solutions.append((d, sigma, exponents, q_solution))
    return profiles, solved_systems, pruned, solutions


def truncated_product(left: list[Fraction], right: list[Fraction], top: int) -> list[Fraction]:
    out = [Fraction(0) for _ in range(top + 1)]
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            if i + j <= top:
                out[i + j] += a * b
    return out


def wrap_residues(roots: tuple[int, ...], exponents: tuple[int, ...]) -> tuple[Fraction, ...]:
    """Residues of dx/g by exact truncated local binomial products."""
    out = []
    for i, (alpha, e_i) in enumerate(zip(roots, exponents)):
        top = e_i - 1
        series = [Fraction(1)] + [Fraction(0) for _ in range(top)]
        for j, (beta, e_j) in enumerate(zip(roots, exponents)):
            if i == j:
                continue
            delta = alpha - beta
            factor = [
                Fraction((-1) ** k * comb(e_j + k - 1, k), delta ** (e_j + k))
                for k in range(top + 1)
            ]
            series = truncated_product(series, factor, top)
        out.append(series[top])
    return tuple(out)


def check_wrap_residue_universe() -> tuple[int, int]:
    profiles = 0
    zero_vectors = 0
    for root_count in range(1, 5):
        roots = tuple(3 * j - 2 for j in range(root_count))
        for exponents in product(range(1, 7), repeat=root_count):
            residues = wrap_residues(roots, exponents)
            exact = all(value == 0 for value in residues)
            predicted = root_count == 1 and exponents[0] > 1
            require(exact == predicted, ("wrap residue profile", roots, exponents, residues))
            profiles += 1
            zero_vectors += int(exact)
    return profiles, zero_vectors


def check_wrap_primitives_and_evaluation() -> tuple[int, int]:
    primitive_checks = 0
    evaluation_checks = 0

    for d in range(2, 10):
        for e in range(2, 14):
            for alpha, gamma in ((-3, sp.Rational(2, 3)), (2, sp.Rational(-4, 5))):
                    y = x - alpha
                    g = gamma * y**e
                    primitive = y ** (1 - e) / (d * gamma * (e - 1))
                    require(
                        sp.cancel(-d * g * sp.diff(primitive, x) - 1) == 0,
                        ("wrap primitive", d, e, alpha, gamma),
                    )
                    primitive_checks += 1

    # For a weight-zero input p in B, p-p(t) is divisible by x-t.  The
    # resulting q proves D(B)=D(z^d B), the load-bearing wrap sidecar.
    for k in range(1, 5):
        for m in range(0, 7):
            g = (x + 2) ** 2 * (x - 1) ** 3
            g_t = g.subs(x, t)
            numerator = sp.Poly(
                sp.expand((x**m + 2 * x + 3) * g_t**k - (t**m + 2 * t + 3) * g**k),
                x,
            )
            quotient, remainder = sp.div(numerator, sp.Poly(x - t, x))
            require(remainder.is_zero, ("wrap evaluation splitting", k, m, remainder))
            # Thus q=-quotient/[g^(k-1)g(t)^k] lies in B and reconstructs
            # p-p(t)=q(t-x)/g without a large rational simplification.
            require(
                sp.Poly(numerator.as_expr() - (x - t) * quotient.as_expr(), x).is_zero,
                ("wrap evaluation reconstruction", k, m),
            )
            evaluation_checks += 1
    return primitive_checks, evaluation_checks


def check_structural_coefficients() -> tuple[int, int, int]:
    local_checks = 0
    diagonal_checks = 0
    infinity_checks = 0
    for d in range(2, 14):
        for sigma in range(1, d + 1):
            for e in range(1, 14):
                for pole_order in range(0, 10):
                    require(
                        sigma * e + d * pole_order != 0,
                        ("local pole coefficient", d, sigma, e, pole_order),
                    )
                    local_checks += 1
            for r in range(1, 30):
                for degree in range(1, 20):
                    predicted = d * degree + sigma * (1 - r)
                    # Direct leading-term bookkeeping: sigma*q-x*(sigma*r-d*m)q/x.
                    direct = sigma - (sigma * r - d * degree)
                    require(predicted == direct, ("leading diagonal", d, sigma, r, degree))
                    diagonal_checks += 1

                    if sigma < d and d * degree == sigma * (r - 1):
                        for leading_degree in range(2, 20):
                            require(
                                d * leading_degree != sigma * r,
                                ("two infinity collision", d, sigma, r, degree, leading_degree),
                            )
                            infinity_checks += 1
    return local_checks, diagonal_checks, infinity_checks


def check_integral_arrow_exponents() -> int:
    checks = 0
    for d in range(2, 15):
        for sigma in range(1, d + 1):
            for e in range(2, 18):
                if sigma * (e - 1) % d:
                    continue
                depth = sigma * (e - 1) // d
                require(depth >= 1, ("positive depth", d, sigma, e))
                for ell in range(depth):
                    arrow = d * (ell + 1) - sigma * (e - 1)
                    require(
                        (arrow == 0) == (ell == depth - 1),
                        ("minimal annihilator arrow", d, sigma, e, ell, arrow),
                    )
                    checks += 1
    return checks


def resonant_defect_completion(d: int, sigma: int, exponents: tuple[int, ...]) -> int:
    """Verify the infinity truncation and complete low-defect interpolation."""
    root_count = len(exponents)
    degree = candidate_degree(d, sigma, exponents)
    require(degree is not None, ("missing resonance", d, sigma, exponents))
    roots = tuple(range(root_count))
    rad = sp.prod(x - alpha for alpha in roots)
    p_coefficients = sp.symbols(f"u0:{degree - root_count + 1}")
    p = sum(p_coefficients[j] * x**j for j in range(degree - root_count + 1))
    q = sp.expand(rad * p)
    logarithmic_q = sp.expand(
        sum(e * sp.cancel(q / (x - alpha)) for alpha, e in zip(roots, exponents))
    )
    image = sp.Poly(
        sp.expand(sigma * q + (t - x) * (sigma * logarithmic_q - d * sp.diff(q, x))),
        x,
    )
    equations = [image.coeff_monomial(x**j) for j in range(root_count, degree)]
    equations.append(p_coefficients[-1] - 1)
    solution_set = sp.linsolve(equations, p_coefficients)
    solution = next(iter(solution_set), None)
    require(solution is not None, ("defect solution", d, sigma, exponents))
    solution_map = dict(zip(p_coefficients, solution))
    p_solution = sp.expand(p.subs(solution_map))
    defect = sp.Poly(sp.expand(image.as_expr().subs(solution_map)), x)
    require(not defect.is_zero, ("zero resonant defect", d, sigma, exponents))

    excess = degree - root_count
    horizontal_coefficients: list[sp.Expr] = [sp.Integer(1)] + [
        sp.Integer(0) for _ in range(excess)
    ]
    horizontal_factors = [
        [
            sp.binomial(sp.Rational(sigma * e, d) - 1, n) * (-alpha) ** n
            for n in range(excess + 1)
        ]
        for alpha, e in zip(roots, exponents)
    ]
    horizontal_factors.append(
        [
            sp.binomial(-sp.Rational(sigma, d), n) * (-t) ** n
            for n in range(excess + 1)
        ]
    )
    for factor in horizontal_factors:
        horizontal_coefficients = [
            sp.expand(
                sum(
                    horizontal_coefficients[j] * factor[n - j]
                    for j in range(n + 1)
                )
            )
            for n in range(excess + 1)
        ]

    horizontal_truncation = sp.expand(
        sum(
            horizontal_coefficients[n] * x ** (excess - n)
            for n in range(excess + 1)
        )
    )
    require(
        sp.expand(p_solution - horizontal_truncation) == 0,
        ("horizontal truncation", d, sigma, exponents, p_solution),
    )

    interpolated_defect = sp.expand(
        sum(
            (t - alpha)
            * (sigma * e - d)
            * horizontal_truncation.subs(x, alpha)
            * sp.cancel(rad / (x - alpha))
            for alpha, e in zip(roots, exponents)
        )
    )
    require(
        sp.expand(defect.as_expr() - interpolated_defect) == 0,
        ("full defect interpolation", d, sigma, exponents, defect.as_expr()),
    )

    asymptotic_scalar = sp.prod(
        sp.Rational(d * (h - 1) + sigma, d * h) for h in range(1, excess + 1)
    )
    logarithmic_rad = sp.expand(
        sum(e * sp.cancel(rad / (x - alpha)) for alpha, e in zip(roots, exponents))
    )
    response_polynomial = sp.expand(sigma * logarithmic_rad - d * sp.diff(rad, x))
    top_t_coefficient = sp.Poly(defect.as_expr(), t).coeff_monomial(t ** (excess + 1))
    require(
        sp.expand(top_t_coefficient - asymptotic_scalar * response_polynomial) == 0,
        ("defect asymptotic", d, sigma, exponents, defect.as_expr()),
    )
    return defect.degree()


def check_resonant_defect_completions() -> int:
    checks = 0
    for d in range(2, 7):
        for sigma in range(1, d + 1):
            for root_count in range(1, 3):
                for exponents in product(range(1, 4), repeat=root_count):
                    if candidate_degree(d, sigma, exponents) is None:
                        continue
                    require(
                        resonant_defect_completion(d, sigma, exponents) == root_count - 1,
                        ("defect top low row", d, sigma, exponents),
                    )
                    checks += 1

    for d, sigma, exponents in (
        (2, 1, (3, 2, 2)),
        (3, 1, (4, 3, 3)),
        (4, 2, (3, 2, 2)),
        (5, 5, (2, 1, 1)),
    ):
        require(
            resonant_defect_completion(d, sigma, exponents) == len(exponents) - 1,
            ("three-root defect", d, sigma, exponents),
        )
        checks += 1
    return checks


def check_sharp_nonwrap_hostiles() -> int:
    checks = 0
    for d, sigma, exponents in (
        (2, 1, (3, 4)),
        (4, 2, (3, 4)),
        (6, 2, (3, 4)),
        (6, 3, (3, 4)),
    ):
        require(candidate_degree(d, sigma, exponents) is not None, ("hostile gate", d, sigma, exponents))
        exists, _ = solve_nonwrap_profile(d, sigma, exponents)
        require(not exists, ("multiroot hostile survived", d, sigma, exponents))
        checks += 1
    return checks


def main() -> None:
    primitive_checks = check_closed_nonwrap_primitives()
    profiles, systems, pruned, solutions = check_nonwrap_direct_universe()
    wrap_profiles, wrap_zero_vectors = check_wrap_residue_universe()
    wrap_primitives, wrap_evaluations = check_wrap_primitives_and_evaluation()
    local_checks, diagonal_checks, infinity_checks = check_structural_coefficients()
    arrow_checks = check_integral_arrow_exponents()
    defect_checks = check_resonant_defect_completions()
    hostile_checks = check_sharp_nonwrap_hostiles()

    print("ALL-SECTOR CONSTANT OBSERVER PACKET -- EXACT HOSTILE REFEREE")
    print(f"closed selected nonwrap primitives: {primitive_checks}")
    print(f"nonwrap ordered multiplicity/character profiles: {profiles}")
    print(f"nonwrap profiles removed by degree/root-count gate: {pruned}")
    print(f"exact nonwrap Q(t)-linear systems solved: {systems}")
    print(f"wrap residue profiles: {wrap_profiles}")
    print(f"wrap zero-residue vectors: {wrap_zero_vectors}")
    print(f"wrap rational primitive identities: {wrap_primitives}")
    print(f"wrap evaluation-at-t splittings: {wrap_evaluations}")
    print(f"local pole/intersection coefficients: {local_checks}")
    print(f"polynomial leading-diagonal identities: {diagonal_checks}")
    print(f"two-infinity incompatibility checks: {infinity_checks}")
    print(f"minimal integral arrow checks: {arrow_checks}")
    print(f"exact resonant defect completions: {defect_checks}")
    print(f"sharp multiroot nonwrap hostiles: {hostile_checks}")
    print("surviving nonwrap direct profiles and q(x,t):")
    for d, sigma, exponents, q_solution in solutions:
        print(f"  d={d}, sigma={sigma}, exponents={exponents}: q={q_solution}")
    print("multiroot constant-observer survivors (all characters): 0")
    print("CLASSIFICATION VERIFIED ON DECLARED FINITE UNIVERSES")


if __name__ == "__main__":
    main()
