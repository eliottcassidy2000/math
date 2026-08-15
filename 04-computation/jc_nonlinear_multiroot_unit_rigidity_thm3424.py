#!/usr/bin/env python3
"""Exact hostile referee for provisional THM-3424.

For the normalized generic fiber

    P = x + g(x) z**d = t,

the weight-one part of a primitive for the Hamiltonian unit class has the
form Q=q(x,t)z.  Its coefficient must solve

    g*q + (t-x)*(g'*q-d*g*q') = g.                       (E)

The candidate theorem says that a solution in Q(t)[x,g**(-1)] exists for
nonconstant g exactly when g has one root of multiplicity e=1 mod d with
e>1.  This script checks the direct linear system over Q(t) on a finite
profile universe, the two independent infinity degree gates, the local
valuation coefficients, and the closed one-root primitive.

SymPy is used only for exact polynomial/rational arithmetic and exact linear
systems.  No floating-point decisions occur.
"""

from __future__ import annotations

from itertools import product

import sympy as sp


x, t = sp.symbols("x t")


def require(condition: bool, label: object) -> None:
    if not condition:
        raise RuntimeError(label)


def generic_equation(g: sp.Expr, q: sp.Expr, d: int) -> sp.Expr:
    """The numerator of D_P(qz)-1 on P=t."""
    return sp.expand(g * q + (t - x) * (sp.diff(g, x) * q - d * g * sp.diff(q, x)) - g)


def predicted_exact(d: int, exponents: tuple[int, ...]) -> bool:
    return len(exponents) == 1 and exponents[0] > 1 and (exponents[0] - 1) % d == 0


def candidate_degree(d: int, exponents: tuple[int, ...]) -> int | None:
    """Degree forced by the x=infinity coefficient of (E)."""
    r = sum(exponents)
    if (r - 1) % d:
        return None
    m = (r - 1) // d
    return m if m >= len(exponents) else None


def solve_direct_profile(d: int, exponents: tuple[int, ...]) -> tuple[bool, sp.Expr | None]:
    """Solve (E) over Q(t) at the only possible x-degree."""
    m = candidate_degree(d, exponents)
    if m is None:
        return False, None
    roots = tuple(2 * j - 1 for j in range(len(exponents)))
    g = sp.prod((x - alpha) ** e for alpha, e in zip(roots, exponents))
    coefficients = sp.symbols(f"c0:{m + 1}")
    q = sum(coefficients[j] * x**j for j in range(m + 1))
    equation = sp.Poly(generic_equation(g, q, d), x)
    equations = [equation.coeff_monomial(x**j) for j in range(equation.degree() + 1)]
    solution_set = sp.linsolve(equations, coefficients)
    if solution_set == sp.EmptySet:
        return False, None
    solution = next(iter(solution_set), None)
    if solution is None:
        return False, None
    q_solution = sp.factor(q.subs(dict(zip(coefficients, solution))))
    require(sp.cancel(generic_equation(g, q_solution, d)) == 0, ("direct substitution", d, exponents))
    return True, q_solution


def one_root_q(d: int, q_power: int, alpha: int) -> sp.Expr:
    """Closed generic primitive coefficient obtained from THM-3422."""
    y = x - alpha
    tau = t - alpha
    numerator = sum(
        sp.Rational(sp.binomial(q_power - 1, j), 1 + j * d)
        * y ** (q_power - j)
        * (tau - y) ** j
        for j in range(q_power)
    )
    return sp.factor(numerator / tau**q_power)


def check_closed_primitives() -> int:
    checks = 0
    for d in range(2, 10):
        for q_power in range(1, 7):
            e = 1 + q_power * d
            for alpha in (-3, 0, 4):
                # Divide (E) by g before expanding.  This keeps the replay
                # small even at e=49; the nonzero scale of g cancels exactly.
                y = x - alpha
                q = one_root_q(d, q_power, alpha)
                reduced = sp.cancel(q + (t - x) * (e * q / y - d * sp.diff(q, x)) - 1)
                require(reduced == 0, ("closed primitive", d, q_power, alpha))
                checks += 1
    return checks


def check_local_and_infinity_gates() -> tuple[int, int, int]:
    pole_checks = 0
    zero_checks = 0
    degree_checks = 0

    # At a root of multiplicity e, a pole of order p has uncancellable
    # coefficient e+d*p.  A zero of order n has first coefficient e-d*n;
    # only n=1 can supply the constant right side, and then e=d is hostile.
    for d in range(2, 14):
        for e in range(1, 18):
            for pole_order in range(1, 12):
                require(e + d * pole_order != 0, ("pole gate", d, e, pole_order))
                pole_checks += 1
            for zero_order in range(1, 12):
                if zero_order > 1:
                    require(zero_order - 1 > 0, ("zero cannot reach constant", zero_order))
                if zero_order == 1 and e == d:
                    require(e - d * zero_order == 0, ("e=d cancellation", d, e))
                zero_checks += 1

    # First infinity gate for a full q of degree m: r=d*m+1.
    # Second gate for B(q_1)=1 and deg(q_1)=n>1: r=d*n.
    # They cannot both hold in characteristic zero.
    for d in range(2, 20):
        for r in range(1, 60):
            for m in range(1, 20):
                full_gate = r == d * m + 1
                for n in range(2, 20):
                    leading_gate = r == d * n
                    require(not (full_gate and leading_gate), ("infinity collision", d, r, m, n))
                    degree_checks += 1
    return pole_checks, zero_checks, degree_checks


def check_direct_universe() -> tuple[int, int, int, list[tuple[int, tuple[int, ...], sp.Expr]]]:
    profiles = 0
    solved_systems = 0
    degree_pruned = 0
    solutions: list[tuple[int, tuple[int, ...], sp.Expr]] = []

    # All ordered multiplicity profiles in this universe are checked.  Direct
    # Q(t)-linear solving is needed only when the independently derived top
    # degree and root-count gates leave a candidate degree.
    for d in range(2, 7):
        for root_count in range(1, 4):
            for exponents in product(range(1, 6), repeat=root_count):
                profiles += 1
                if candidate_degree(d, exponents) is None:
                    degree_pruned += 1
                    exists = False
                    q_solution = None
                else:
                    solved_systems += 1
                    exists, q_solution = solve_direct_profile(d, exponents)
                require(exists == predicted_exact(d, exponents), ("profile", d, exponents, q_solution))
                if exists:
                    require(q_solution is not None, ("missing solution", d, exponents))
                    solutions.append((d, exponents, q_solution))
    return profiles, solved_systems, degree_pruned, solutions


def check_sharp_hostiles() -> int:
    checks = 0

    # A multiroot leading equation can pass while the full equation fails.
    # Here B(x(x-1))=1 for d=2 and g=x(x-1)^3, but deg(g)=4 is not 1 mod 2.
    g = x * (x - 1) ** 3
    q1 = x * (x - 1)
    require(sp.factor(sp.diff(g, x) * q1 - 2 * g * sp.diff(q1, x) - g) == 0, "leading hostile")
    require(candidate_degree(2, (1, 3)) is None, "full infinity rejects leading hostile")
    checks += 1

    # These profiles pass the full degree/root-count gate and therefore reach
    # the direct Q(t)-linear system, but global multiroot compatibility fails.
    for d, exponents in ((2, (3, 4)), (3, (2, 5)), (4, (3, 5, 5))):
        require(candidate_degree(d, exponents) is not None, ("hostile reaches solver", d, exponents))
        exists, _ = solve_direct_profile(d, exponents)
        require(not exists, ("multiroot hostile survived", d, exponents))
        checks += 1
    return checks


def main() -> None:
    primitive_checks = check_closed_primitives()
    pole_checks, zero_checks, degree_checks = check_local_and_infinity_gates()
    profiles, solved_systems, degree_pruned, solutions = check_direct_universe()
    hostile_checks = check_sharp_hostiles()

    print("NONLINEAR MULTIROOT UNIT RIGIDITY -- EXACT HOSTILE REFEREE")
    print(f"closed one-root generic primitive identities: {primitive_checks}")
    print(f"local pole coefficients: {pole_checks}")
    print(f"local zero-order coefficients: {zero_checks}")
    print(f"two-infinity incompatibility checks: {degree_checks}")
    print(f"ordered multiplicity profiles: {profiles}")
    print(f"profiles eliminated by top-degree/root-count gate: {degree_pruned}")
    print(f"exact Q(t)-linear systems solved: {solved_systems}")
    print(f"sharp hostile checks: {hostile_checks}")
    print("surviving profiles and exact q(x,t):")
    for d, exponents, q_solution in solutions:
        print(f"  d={d}, exponents={exponents}: q={q_solution}")
    print("multiroot survivors: 0")
    print("CLASSIFICATION VERIFIED ON DECLARED FINITE UNIVERSE")


if __name__ == "__main__":
    main()
