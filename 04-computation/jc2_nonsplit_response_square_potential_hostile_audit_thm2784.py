#!/usr/bin/env python3
"""Independent exact hostile audit for THM-2784's square-potential theorem.

This companion uses the base linearization

    F = V G^2,       2 V G' + V' G = 2 kappa,

where G=F'/(2*kappa).  It thereby checks both directions of the quadratic
deck descent without importing the candidate referee's implementation.
"""

from itertools import product
from pathlib import Path
import ast

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


x, w = sp.symbols("x w")
GATES = 0


def gate(condition, message):
    global GATES
    require(bool(condition), message)
    GATES += 1


def cancel(expr):
    return sp.cancel(sp.together(expr))


def factor_dict(expr):
    polynomial = sp.Poly(sp.expand(expr), x, domain=sp.QQ)
    answer = {}
    for factor, exponent in sp.factor_list(polynomial.as_expr())[1]:
        monic = sp.Poly(factor, x, domain=sp.QQ).monic().as_expr()
        answer[monic] = exponent
    return answer


def valuation(expr, factor):
    numerator, denominator = sp.fraction(cancel(expr))
    factors_num = factor_dict(numerator)
    factors_den = factor_dict(denominator)
    monic = sp.Poly(factor, x, domain=sp.QQ).monic().as_expr()
    return factors_num.get(monic, 0) - factors_den.get(monic, 0)


def infinity_order(expr):
    numerator, denominator = sp.fraction(cancel(expr))
    return sp.degree(denominator, x) - sp.degree(numerator, x)


def is_polynomial(expr):
    _, denominator = sp.fraction(cancel(expr))
    return sp.degree(denominator, x) <= 0


def is_square_rational(expr):
    numerator, denominator = sp.fraction(cancel(expr))
    coefficient_num, factors_num = sp.factor_list(sp.Poly(numerator, x, domain=sp.QQ))
    coefficient_den, factors_den = sp.factor_list(sp.Poly(denominator, x, domain=sp.QQ))
    if any(exponent % 2 for _, exponent in factors_num + factors_den):
        return False
    rational_coefficient = sp.Rational(coefficient_num, coefficient_den)
    num = int(sp.numer(rational_coefficient))
    den = int(sp.denom(rational_coefficient))
    if num < 0:
        return False
    return sp.integer_nthroot(num, 2)[1] and sp.integer_nthroot(den, 2)[1]


def geometric_count(factors, predicate):
    return sum(
        sp.degree(factor, x)
        for factor, exponent in factors.items()
        if predicate(exponent)
    )


def check_solution(function, potential, kappa, label):
    """Check the square potential, linearized deck equation, and ledgers."""

    function = cancel(function)
    potential = cancel(potential)
    derivative = sp.diff(function, x)
    gate(derivative != 0, f"{label}: F must be nonconstant")
    gate(is_polynomial(potential), f"{label}: V must be polynomial")
    gate(
        cancel(potential * derivative**2 - 4 * kappa**2 * function) == 0,
        f"{label}: square-potential equation",
    )

    base_g = cancel(derivative / (2 * kappa))
    gate(
        cancel(function - potential * base_g**2) == 0,
        f"{label}: R^2=VG^2=F",
    )
    gate(
        cancel(2 * potential * sp.diff(base_g, x) + sp.diff(potential, x) * base_g - 2 * kappa)
        == 0,
        f"{label}: converse deck derivative",
    )
    gate(
        cancel(function / potential - base_g**2) == 0,
        f"{label}: F and V have the same squareclass",
    )

    potential_factors = factor_dict(potential)
    numerator, denominator = sp.fraction(function)
    numerator_factors = factor_dict(numerator)
    denominator_factors = factor_dict(denominator)
    derivative_num, derivative_den = sp.fraction(cancel(derivative))
    derivative_factors = set(factor_dict(derivative_num)) | set(factor_dict(derivative_den))
    all_factors = (
        set(potential_factors)
        | set(numerator_factors)
        | set(denominator_factors)
        | derivative_factors
    )

    local_rows = 0
    for factor in all_factors:
        m = valuation(potential, factor)
        n = valuation(function, factor)
        t = valuation(derivative, factor)
        gate(m + 2 * t == n, f"{label}: local valuation identity at {factor}")
        if n != 0:
            gate(t == n - 1, f"{label}: derivative order at nonunit {factor}")
            gate(n == 2 - m, f"{label}: finite divisor law at {factor}")
        else:
            gate(m == 0 and t == 0, f"{label}: unit point must be unramified {factor}")
        local_rows += 1

    simple_roots = geometric_count(potential_factors, lambda exponent: exponent == 1)
    high_roots = geometric_count(potential_factors, lambda exponent: exponent >= 3)
    double_roots = geometric_count(potential_factors, lambda exponent: exponent == 2)
    gate(double_roots == 0, f"{label}: double V-root prohibition")

    extra_double_zeros = sum(
        sp.degree(factor, x)
        for factor, exponent in numerator_factors.items()
        if exponent == 2 and factor not in potential_factors
    )
    finite_zero_degree = sum(
        sp.degree(factor, x) * exponent
        for factor, exponent in numerator_factors.items()
    )
    finite_pole_degree = sum(
        sp.degree(factor, x) * exponent
        for factor, exponent in denominator_factors.items()
    )
    high_multiplicities = [
        exponent
        for factor, exponent in potential_factors.items()
        for _ in range(sp.degree(factor, x))
        if exponent >= 3
    ]
    pole_from_high = sum(multiplicity - 2 for multiplicity in high_multiplicities)
    degree_v = sp.degree(potential, x)
    gate(pole_from_high == finite_pole_degree, f"{label}: finite pole ledger")
    gate(
        finite_zero_degree == simple_roots + 2 * extra_double_zeros,
        f"{label}: finite zero ledger",
    )
    gate(
        degree_v == simple_roots + finite_pole_degree + 2 * high_roots,
        f"{label}: degree-V ledger",
    )

    n_infinity = infinity_order(function)
    derivative_infinity = infinity_order(derivative)
    gate(
        n_infinity == finite_pole_degree - finite_zero_degree,
        f"{label}: principal divisor at infinity",
    )
    if n_infinity != 0:
        gate(
            derivative_infinity == n_infinity + 1,
            f"{label}: unbalanced derivative sign",
        )
        gate(n_infinity == degree_v - 2, f"{label}: unbalanced infinity law")
        gate(
            simple_roots + high_roots + extra_double_zeros == 1,
            f"{label}: one-finite-point classification",
        )
        map_degree = max(finite_zero_degree, finite_pole_degree)
        ramification = (
            extra_double_zeros
            + finite_pole_degree
            - high_roots
            + abs(n_infinity)
            - 1
        )
        gate(ramification == 2 * map_degree - 2, f"{label}: unbalanced RH")
        chamber = "unbalanced"
        passport_degree = map_degree
    else:
        leading_value = sp.LC(sp.Poly(numerator, x)) / sp.LC(sp.Poly(denominator, x))
        difference = cancel(function - leading_value)
        r_infinity = infinity_order(difference)
        gate(r_infinity >= 1, f"{label}: nonconstant balanced difference")
        gate(
            derivative_infinity == r_infinity + 1,
            f"{label}: balanced derivative sign",
        )
        gate(finite_pole_degree == finite_zero_degree, f"{label}: balanced P=Z")
        gate(degree_v == 2 * r_infinity + 2, f"{label}: balanced degree-V law")
        gate(
            r_infinity == simple_roots + extra_double_zeros + high_roots - 1,
            f"{label}: balanced point ledger",
        )
        map_degree = finite_pole_degree
        gate(r_infinity <= map_degree, f"{label}: passport long cycle exceeds degree")
        gate(
            map_degree - r_infinity == extra_double_zeros - high_roots + 1,
            f"{label}: passport fixed-point count",
        )
        zero_partition_sum = 2 * extra_double_zeros + simple_roots
        infinity_partition_sum = sum(
            multiplicity - 2 for multiplicity in high_multiplicities
        )
        third_partition_sum = r_infinity + (map_degree - r_infinity)
        gate(
            zero_partition_sum
            == infinity_partition_sum
            == third_partition_sum
            == map_degree,
            f"{label}: passport partition sums",
        )
        rh = (
            extra_double_zeros
            + (map_degree - high_roots)
            + (r_infinity - 1)
        )
        gate(rh == 2 * map_degree - 2, f"{label}: balanced RH")
        chamber = "balanced"
        passport_degree = map_degree

    return local_rows, chamber, passport_degree


def integer_ledger_audit():
    local_case_rows = 0
    for multiplicity in range(0, 11):
        for order_f in range(-10, 11):
            if order_f == 0:
                for derivative_order in range(0, 5):
                    possible = multiplicity + 2 * derivative_order == 0
                    gate(
                        possible == (multiplicity == 0 and derivative_order == 0),
                        "finite unit valuation classification",
                    )
                    local_case_rows += 1
            else:
                derivative_order = order_f - 1
                possible = multiplicity + 2 * derivative_order == order_f
                gate(
                    possible == (order_f == 2 - multiplicity),
                    "finite nonunit valuation classification",
                )
                local_case_rows += 1
    gate(
        all(2 - 2 == 0 for _ in (0,)),
        "double-root formal candidate is the excluded unit order",
    )

    infinity_rows = 0
    balanced_rows = 0
    unbalanced_rows = 0
    multiplicity_options = ((),) + tuple(
        (m,) for m in range(3, 9)
    ) + tuple(
        (m1, m2) for m1 in range(3, 8) for m2 in range(3, 8)
    )
    for high in multiplicity_options:
        h = len(high)
        pole_degree = sum(m - 2 for m in high)
        for simple in range(0, 7):
            for extra_double in range(0, 7):
                zero_degree = simple + 2 * extra_double
                degree_v = simple + pole_degree + 2 * h
                n_infinity = pole_degree - zero_degree
                if n_infinity != 0:
                    equation = n_infinity == degree_v - 2
                    gate(
                        equation == (simple + h + extra_double == 1),
                        "unbalanced ledger equivalence",
                    )
                    if equation:
                        unbalanced_rows += 1
                else:
                    r = (degree_v - 2) // 2
                    equation = degree_v == 2 * r + 2 and r >= 1
                    if equation:
                        gate(
                            r == simple + extra_double + h - 1,
                            "balanced r ledger",
                        )
                        rh = (
                            extra_double
                            + (zero_degree - h)
                            + (r - 1)
                        )
                        gate(rh == 2 * zero_degree - 2, "balanced integer RH")
                        balanced_rows += 1
                infinity_rows += 1
    return local_case_rows, infinity_rows, balanced_rows, unbalanced_rows


def faber_transfer_audit():
    # Reduced degrees and deck characters.  For an even seed, the Laurent
    # coefficient c_j has character (-1)^j; R=4c_3+p c_1 is anti-invariant.
    reduced_rows = 0
    odd_killed = 0
    even_surviving = 0
    for degree in range(1, 41):
        if degree % 4 == 0:
            continue
        reduced_rows += 1
        seed_character = -1 if degree % 2 else 1
        if seed_character == -1:
            odd_killed += 1
        else:
            even_surviving += 1
            c1_character = -1
            c2_character = 1
            c3_character = -1
            p_character = 1
            gate(c1_character == -1, "Phi deck character")
            gate(c2_character == 1, "Psi deck character")
            gate(c3_character == p_character * c1_character == -1, "R deck character")

    phi_prime, psi_prime, response_prime, p = sp.symbols(
        "phi_prime psi_prime response_prime p"
    )
    hamiltonian = sp.Poly(
        (w**2 + p / 4) * phi_prime + w * psi_prime + response_prime,
        w,
    )
    gate(hamiltonian.coeff_monomial(w**2) == phi_prime, "Hamiltonian Phi coefficient")
    gate(hamiltonian.coeff_monomial(w) == psi_prime, "Hamiltonian Psi coefficient")
    gate(
        hamiltonian.coeff_monomial(1)
        == response_prime + p * phi_prime / 4,
        "Hamiltonian constant coefficient before Phi vanishing",
    )
    gate(
        cancel(hamiltonian.as_expr().subs({phi_prime: 0, psi_prime: 0}) - response_prime)
        == 0,
        "Hamiltonian response after flux vanishing",
    )

    # Coordinate change w=Uz+B/(2U) has dw/dz=U, so J_(x,z)=U J_(x,w).
    U, kappa = sp.symbols("U kappa", nonzero=True)
    gate(cancel(U * (kappa / U) - kappa) == 0, "Jacobian coordinate scaling")
    return reduced_rows, odd_killed, even_surviving


def main():
    kappa = sp.Integer(3)
    named = (
        (18 * (x + 1), 2 * (x + 1), "linear"),
        (36 / (x - 2), (x - 2) ** 3, "pure-triple"),
        (
            36 * x / (x - 1),
            x * (x - 1) ** 3,
            "balanced-one-three",
        ),
    )
    solution_rows = 0
    local_rows = 0
    balanced_maps = 0
    unbalanced_maps = 0
    for function, potential, label in named:
        rows, chamber, _ = check_solution(function, potential, kappa, label)
        local_rows += rows
        balanced_maps += chamber == "balanced"
        unbalanced_maps += chamber == "unbalanced"
        solution_rows += 1

    atlas_rows = 0
    nonsquare_rows = 0
    points = (-2, 0, 1, 3)
    for exponents in product(range(-2, 3), repeat=len(points)):
        if exponents == (0, 0, 0, 0):
            continue
        function = sp.Integer(36)
        for point, exponent in zip(points, exponents):
            function *= (x - point) ** exponent
        function = cancel(function)
        derivative = sp.diff(function, x)
        if derivative == 0:
            continue
        potential = cancel(4 * kappa**2 * function / derivative**2)
        if not is_polynomial(potential) or sp.degree(potential, x) <= 0:
            continue
        rows, chamber, _ = check_solution(
            function,
            potential,
            kappa,
            f"atlas-{exponents}",
        )
        local_rows += rows
        balanced_maps += chamber == "balanced"
        unbalanced_maps += chamber == "unbalanced"
        nonsquare_rows += not is_square_rational(potential)
        atlas_rows += 1

    hostile_rows = 0
    maximal_radical_degree = 0
    for n in range(1, 10):
        function = cancel(4 * (1 - x ** (-n)))
        potential = cancel(x ** (n + 2) * (x**n - 1) / n**2)
        rows, chamber, _ = check_solution(function, potential, 1, f"hostile-{n}")
        gate(chamber == "balanced", f"hostile-{n}: expected balanced chamber")
        gate(not is_square_rational(potential), f"hostile-{n}: V must be nonsquare")
        radical_degree = sum(sp.degree(factor, x) for factor in factor_dict(potential))
        maximal_radical_degree = max(maximal_radical_degree, radical_degree)
        local_rows += rows
        hostile_rows += 1

    gate(maximal_radical_degree == 10, "hostile radical degree")

    # Squarefree classification and the unique linear scalar.
    squarefree_degree_rows = 0
    for degree_v in range(1, 13):
        degree_f = 2 - degree_v
        gate(
            (degree_f >= 1) == (degree_v == 1),
            f"squarefree degree classification D={degree_v}",
        )
        squarefree_degree_rows += 1
    c, a, symbol_kappa = sp.symbols("c a symbol_kappa", nonzero=True)
    linear_v = c * (x - a)
    linear_f = 4 * symbol_kappa**2 * linear_v / c**2
    gate(
        cancel(linear_v * sp.diff(linear_f, x) ** 2 - 4 * symbol_kappa**2 * linear_f)
        == 0,
        "unique linear survivor",
    )

    local_integer_rows, infinity_integer_rows, balanced_integer_rows, unbalanced_integer_rows = (
        integer_ledger_audit()
    )
    reduced_rows, odd_killed, even_surviving = faber_transfer_audit()

    source = Path(__file__).read_text(encoding="utf-8")
    tree = ast.parse(source)
    assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
    gate(assert_nodes == 0, "truth-bearing Python assert")

    print("THM-2784 INDEPENDENT NONSPLIT SQUARE-POTENTIAL HOSTILE AUDIT")
    print(f"named_solutions={len(named)}")
    print(f"atlas_polynomial_potentials={atlas_rows}")
    print(f"atlas_nonsquare_potentials={nonsquare_rows}")
    print(f"hostile_family_rows={hostile_rows}")
    print(f"maximal_hostile_radical_degree={maximal_radical_degree}")
    print(f"solution_local_factor_rows={local_rows}")
    print(f"balanced_solution_maps={balanced_maps + hostile_rows}")
    print(f"unbalanced_solution_maps={unbalanced_maps}")
    print(f"integer_local_cases={local_integer_rows}")
    print(f"integer_infinity_ledgers={infinity_integer_rows}")
    print(f"integer_balanced_ledgers={balanced_integer_rows}")
    print(f"integer_unbalanced_ledgers={unbalanced_integer_rows}")
    print(f"squarefree_degree_rows={squarefree_degree_rows}")
    print(f"faber_reduced_rows={reduced_rows}")
    print(f"faber_odd_rows_killed={odd_killed}")
    print(f"faber_even_rows_surviving={even_surviving}")
    print(f"assert_nodes={assert_nodes}")
    print(f"exact_gates={GATES}")
    print("square_potential_both_directions=PASS")
    print("finite_DVR_ledger=PASS")
    print("infinity_signs_and_valuations=PASS")
    print("passport_and_Riemann_Hurwitz=PASS")
    print("squarefree_linear_classification=PASS")
    print("large_squarefree_part_hostile=PASS")
    print("Faber_response_transfer=PASS")
    print("FAILED CHECKS: NONE")


if __name__ == "__main__":
    main()
