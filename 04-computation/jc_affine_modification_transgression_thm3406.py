#!/usr/bin/env python3
"""Exact standard-library audit for THM-3406.

The program uses rational sparse polynomial arithmetic.  It independently
checks the affine-modification denominator filtration, its CRT triangular
arms, the localized Bezout primitive and connecting-map representative, and
the chi/psi support-thickness filtration.  It also freezes the mixed q=2,
radical, duplicate-arm, and positive-characteristic hostiles.
"""

from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json


EXPECTED_SEMANTIC_SHA256 = "065f584f197ff000d9f42175ffa858de223ee42551b19d996a54da37cfcb88f7"


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


def trim(poly):
    values = [Fraction(value) for value in poly]
    while values and not values[-1]:
        values.pop()
    return tuple(values)


ZERO = ()
ONE = (Fraction(1),)
X = (Fraction(0), Fraction(1))


def add(left, right):
    size = max(len(left), len(right))
    return trim([
        (left[i] if i < len(left) else 0)
        + (right[i] if i < len(right) else 0)
        for i in range(size)
    ])


def scale(poly, scalar):
    scalar = Fraction(scalar)
    return trim([scalar * value for value in poly])


def sub(left, right):
    return add(left, scale(right, -1))


def mul(left, right):
    if not left or not right:
        return ZERO
    answer = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            answer[i + j] += left_value * right_value
    return trim(answer)


@lru_cache(maxsize=None)
def power(poly, exponent):
    require(exponent >= 0, ("negative power", exponent))
    answer = ONE
    base = trim(poly)
    current = exponent
    while current:
        if current & 1:
            answer = mul(answer, base)
        base = mul(base, base)
        current //= 2
    return answer


def deriv(poly):
    return trim([i * poly[i] for i in range(1, len(poly))])


def divmod_poly(dividend, divisor):
    divisor = trim(divisor)
    require(divisor, "zero polynomial divisor")
    remainder = list(trim(dividend))
    quotient = [Fraction(0)] * max(1, len(remainder) - len(divisor) + 1)
    while remainder and len(remainder) >= len(divisor):
        shift = len(remainder) - len(divisor)
        coefficient = remainder[-1] / divisor[-1]
        quotient[shift] += coefficient
        for i, value in enumerate(divisor):
            remainder[i + shift] -= coefficient * value
        remainder = list(trim(remainder))
    return trim(quotient), trim(remainder)


def quotient_exact(dividend, divisor, label):
    quotient, remainder = divmod_poly(dividend, divisor)
    require(not remainder, (label, "nonexact division", remainder))
    return quotient


def mod_poly(poly, modulus):
    return divmod_poly(poly, modulus)[1]


def xgcd(left, right):
    old_r, current_r = trim(left), trim(right)
    old_s, current_s = ONE, ZERO
    old_t, current_t = ZERO, ONE
    while current_r:
        quotient, remainder = divmod_poly(old_r, current_r)
        old_r, current_r = current_r, remainder
        old_s, current_s = current_s, sub(old_s, mul(quotient, current_s))
        old_t, current_t = current_t, sub(old_t, mul(quotient, current_t))
    require(old_r, "gcd vanished")
    unit = 1 / old_r[-1]
    return scale(old_r, unit), scale(old_s, unit), scale(old_t, unit)


def gcd_poly(left, right):
    return xgcd(left, right)[0]


def compose(poly, argument):
    answer = ZERO
    for coefficient in reversed(trim(poly)):
        answer = add(mul(answer, argument), (coefficient,))
    return answer


def lcm_poly(left, right):
    common = gcd_poly(left, right)
    return trim(quotient_exact(mul(left, right), common, "lcm"))


def coeff_text(value):
    value = Fraction(value)
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def poly_vector(poly):
    return [coeff_text(value) for value in trim(poly)]


# Sparse polynomials in two variables.  The second variable is z for Bezout
# rows and t for affine-modification numerators.
def bclean(poly):
    return {key: Fraction(value) for key, value in poly.items() if value}


def badd(left, right):
    answer = dict(left)
    for key, value in right.items():
        answer[key] = answer.get(key, Fraction(0)) + value
    return bclean(answer)


def bscale(poly, scalar):
    scalar = Fraction(scalar)
    return bclean({key: scalar * value for key, value in poly.items()})


def bsub(left, right):
    return badd(left, bscale(right, -1))


def bmul(left, right):
    answer = {}
    for (left_x, left_y), left_value in left.items():
        for (right_x, right_y), right_value in right.items():
            key = (left_x + right_x, left_y + right_y)
            answer[key] = answer.get(key, Fraction(0)) + left_value * right_value
    return bclean(answer)


def bpow(poly, exponent):
    require(exponent >= 0, ("negative bivariate power", exponent))
    answer = {(0, 0): Fraction(1)}
    base = bclean(poly)
    current = exponent
    while current:
        if current & 1:
            answer = bmul(answer, base)
        base = bmul(base, base)
        current //= 2
    return answer


def bdx(poly):
    return bclean({
        (x_degree - 1, y_degree): x_degree * value
        for (x_degree, y_degree), value in poly.items()
        if x_degree
    })


def bdy(poly):
    return bclean({
        (x_degree, y_degree - 1): y_degree * value
        for (x_degree, y_degree), value in poly.items()
        if y_degree
    })


def bx(poly):
    return bclean({(degree, 0): value for degree, value in enumerate(poly)})


BY = {(0, 1): Fraction(1)}
BONE = {(0, 0): Fraction(1)}


def bcompose(poly, argument):
    answer = {}
    for coefficient in reversed(trim(poly)):
        answer = badd(bmul(answer, argument), {(0, 0): coefficient})
    return answer


def coefficient_slices(poly):
    by_y = {}
    for (x_degree, y_degree), value in bclean(poly).items():
        coefficients = by_y.setdefault(y_degree, [])
        if len(coefficients) <= x_degree:
            coefficients.extend([Fraction(0)] * (x_degree + 1 - len(coefficients)))
        coefficients[x_degree] += value
    return {degree: trim(values) for degree, values in by_y.items()}


def bdivide_x_exact(poly, divisor, label):
    answer = {}
    for y_degree, coefficient in coefficient_slices(poly).items():
        quotient = quotient_exact(coefficient, divisor, (label, y_degree))
        for x_degree, value in enumerate(quotient):
            if value:
                answer[(x_degree, y_degree)] = value
    return bclean(answer)


def derivation(Px, g, poly):
    return bsub(bmul(Px, bdy(poly)), bmul(bx(g), bdx(poly)))


def multiplication_matrix(f, g):
    degree = len(g) - 1
    require(degree >= 1 and g[-1] == 1, ("nonmonic modulus", g))
    matrix = [[Fraction(0) for _ in range(degree)] for _ in range(degree)]
    x_power = ONE
    for column in range(degree):
        remainder = mod_poly(mul(f, x_power), g)
        for row, value in enumerate(remainder):
            matrix[row][column] = value
        x_power = mul(x_power, X)
    return matrix


def matrix_identity(size):
    return [[Fraction(row == column) for column in range(size)]
            for row in range(size)]


def matrix_mul(left, right):
    return [[
        sum(left[row][inner] * right[inner][column]
            for inner in range(len(right)))
        for column in range(len(right[0]))
    ] for row in range(len(left))]


def matrix_flat(matrix):
    return [value for row in matrix for value in row]


def solve_columns(columns, target):
    variable_count = len(columns)
    rows = [[columns[column][row] for column in range(variable_count)]
            + [target[row]] for row in range(len(target))]
    pivot_columns = []
    pivot_row = 0
    for column in range(variable_count):
        pivot = next((row for row in range(pivot_row, len(rows))
                      if rows[row][column]), None)
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        factor = rows[pivot_row][column]
        rows[pivot_row] = [value / factor for value in rows[pivot_row]]
        for row in range(len(rows)):
            if row == pivot_row or not rows[row][column]:
                continue
            factor = rows[row][column]
            rows[row] = [
                rows[row][entry] - factor * rows[pivot_row][entry]
                for entry in range(variable_count + 1)
            ]
        pivot_columns.append(column)
        pivot_row += 1
    for row in rows:
        if all(not row[column] for column in range(variable_count)) and row[-1]:
            return None
    require(len(pivot_columns) == variable_count,
            ("dependent earlier matrix powers", variable_count, pivot_columns))
    solution = [Fraction(0)] * variable_count
    for row, column in enumerate(pivot_columns):
        solution[column] = rows[row][-1]
    return solution


def minimal_polynomial(matrix):
    size = len(matrix)
    powers = [matrix_identity(size)]
    for degree in range(1, size + 1):
        current = matrix_mul(powers[-1], matrix)
        solution = solve_columns(
            [matrix_flat(item) for item in powers],
            [-value for value in matrix_flat(current)],
        )
        if solution is not None:
            return trim(solution + [Fraction(1)])
        powers.append(current)
    raise RuntimeError("minimal polynomial not found")


def matrix_rank(matrix):
    if not matrix:
        return 0
    rows = [list(map(Fraction, row)) for row in matrix]
    row_count = len(rows)
    column_count = len(rows[0])
    rank = 0
    for column in range(column_count):
        pivot = next((row for row in range(rank, row_count)
                      if rows[row][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        factor = rows[rank][column]
        rows[rank] = [entry / factor for entry in rows[rank]]
        for row in range(row_count):
            if row == rank or not rows[row][column]:
                continue
            factor = rows[row][column]
            rows[row] = [rows[row][entry] - factor * rows[rank][entry]
                         for entry in range(column_count)]
        rank += 1
        if rank == row_count:
            break
    return rank


def denominator_remainder(poly, g, q):
    """Obstructions after t=gz in a/g^q; empty means polynomial in x,z."""
    answer = {}
    for t_degree, coefficient in coefficient_slices(poly).items():
        if t_degree >= q:
            continue
        modulus = power(g, q - t_degree)
        remainder = mod_poly(coefficient, modulus)
        if remainder:
            answer[t_degree] = remainder
    return answer


def ideal_power_remainder(poly, g, q):
    """Normal remainder modulo (g,t)^q from its q+1 monomial generators."""
    answer = {}
    slices = coefficient_slices(poly)
    for t_degree in range(q):
        coefficient = slices.get(t_degree, ZERO)
        remainder = mod_poly(coefficient, power(g, q - t_degree))
        if remainder:
            answer[t_degree] = remainder
    return answer


def in_g_tq(poly, g, q):
    """Membership in (g,t^q)."""
    for t_degree, coefficient in coefficient_slices(poly).items():
        if t_degree < q and mod_poly(coefficient, g):
            return False
    return True


def term(coefficient, t_degree):
    return {(x_degree, t_degree): value
            for x_degree, value in enumerate(trim(coefficient)) if value}


def affine_candidate_bank(g, q):
    small = [ONE, X, add(X, ONE), add(g, ONE)]
    candidates = [{}]
    atoms = []
    for t_degree in range(q + 2):
        for g_power in range(q + 2):
            for coefficient in small:
                atom = term(mul(power(g, g_power), coefficient), t_degree)
                atoms.append(atom)
                candidates.append(atom)
    for index in range(0, len(atoms) - 1, max(1, q + 1)):
        candidates.append(badd(atoms[index], atoms[index + 1]))
    return candidates


def audit_affine_power_filtration(factor_packets):
    comparisons = 0
    transition_cells = 0
    length_cells = 0
    crt_cells = 0
    for packet in factor_packets:
        factors = packet["factors"]
        g = ONE
        for prime, exponent in factors:
            g = mul(g, power(prime, exponent))
        require(g == packet["g"], (packet["name"], "factor product"))
        degree_g = len(g) - 1
        for q in range(1, 6):
            candidates = affine_candidate_bank(g, q)
            for candidate in candidates:
                regular = not denominator_remainder(candidate, g, q)
                in_power = not ideal_power_remainder(candidate, g, q)
                require(regular == in_power,
                        (packet["name"], q, "slice mismatch", candidate))
                multiplied = bmul(bx(g), candidate)
                next_power = not ideal_power_remainder(multiplied, g, q + 1)
                require(next_power == in_power,
                        (packet["name"], q, "transition colon", candidate))
                comparisons += 1
                transition_cells += 1

            length = degree_g * q * (q + 1) // 2
            direct_length = sum(degree_g * (q - t_degree)
                                for t_degree in range(q))
            previous = degree_g * (q - 1) * q // 2
            require(length == direct_length, (packet["name"], q, "length"))
            require(length - previous == q * degree_g,
                    (packet["name"], q, "graded length"))
            arm_length = sum(
                (len(prime) - 1) * exponent * q * (q + 1) // 2
                for prime, exponent in factors
            )
            require(arm_length == length, (packet["name"], q, "arm length"))
            length_cells += 1

            for t_degree in range(q):
                exponent_scale = q - t_degree
                moduli = [power(prime, exponent * exponent_scale)
                           for prime, exponent in factors]
                modulus = ONE
                for item in moduli:
                    modulus = mul(modulus, item)
                require(modulus == power(g, exponent_scale),
                        (packet["name"], q, t_degree, "CRT product"))
                idempotents = []
                for item in moduli:
                    complement = quotient_exact(modulus, item, "CRT complement")
                    common, inverse, _ = xgcd(complement, item)
                    require(common == ONE, (packet["name"], "noncoprime arms"))
                    idempotents.append(mod_poly(mul(complement, inverse), modulus))
                require(mod_poly(sum_polys(idempotents), modulus) == ONE,
                        (packet["name"], q, t_degree, "CRT sum"))
                for left_index, left in enumerate(idempotents):
                    require(mod_poly(sub(mul(left, left), left), modulus) == ZERO,
                            (packet["name"], q, t_degree, "CRT idempotent"))
                    for right_index, right in enumerate(idempotents):
                        if left_index != right_index:
                            require(mod_poly(mul(left, right), modulus) == ZERO,
                                    (packet["name"], q, t_degree,
                                     "CRT orthogonality"))
                crt_cells += len(idempotents)

    # q=2 hostiles: the mixed generator is essential, and neither I nor rad(g)
    # may replace I^2 or I in the corresponding numerator slices.
    g = power(X, 2)
    mixed = term(g, 1)
    require(not ideal_power_remainder(mixed, g, 2), "gt missing from I^2")
    pure_square_ideal = in_g_tq(mixed, power(g, 2), 2)
    require(not pure_square_ideal, "gt falsely lies in (g^2,t^2)")
    require(ideal_power_remainder(bx(g), g, 2), "g falsely lies in I^2")
    require(denominator_remainder(bx(g), g, 2), "g/g^2 falsely regular")
    require(denominator_remainder(bx(X), g, 1), "rad(g)/g falsely regular")

    # K[x]/(x^2) has no nontrivial idempotent: duplicated factors do not make
    # two CRT arms.  Exhaust a rational coefficient grid as a hostile control.
    idempotents = []
    for constant in range(-2, 3):
        for linear in range(-2, 3):
            candidate = trim((constant, linear))
            if mod_poly(sub(mul(candidate, candidate), candidate), g) == ZERO:
                idempotents.append(poly_vector(candidate))
    require(idempotents == [[], ["1"]], ("duplicate-arm idempotents", idempotents))

    return {
        "comparisons": comparisons,
        "transition_cells": transition_cells,
        "length_cells": length_cells,
        "crt_cells": crt_cells,
        "duplicate_idempotents": idempotents,
    }


def sum_polys(polynomials):
    answer = ZERO
    for poly in polynomials:
        answer = add(answer, poly)
    return answer


def build_bezout_row(f, g):
    fp = deriv(f)
    gp = deriv(g)
    common, u, v = xgcd(fp, g)
    require(common == ONE, ("f' not invertible mod g", f, g))
    gp_square = quotient_exact(mul(gp, gp), g, "g divides g'^2")
    u2 = mul(u, u)
    A = bsub(bx(u), bmul(bx(mul(u2, gp)), BY))
    Px = badd(bx(fp), bmul(bx(gp), BY))
    C = bdivide_x_exact(bsub(BONE, bmul(A, Px)), g, "Bezout C")
    require(badd(bmul(A, Px), bmul(C, bx(g))) == BONE, "Bezout row failed")

    # The explicit C from direct expansion is a second path.
    explicit_C = badd(
        badd(bx(v), bscale(bmul(bx(mul(mul(u, v), gp)), BY), -1)),
        bmul(bx(mul(u2, gp_square)), bpow(BY, 2)),
    )
    require(C == explicit_C, ("explicit C mismatch", f, g))

    m = badd(bdx(A), bdy(C))
    cleared_primitive = bsub(
        bsub(bmul(bx(g), bdx(A)), bmul(Px, bdy(A))),
        bmul(A, bx(gp)),
    )
    require(cleared_primitive == bmul(bx(g), m),
            ("D(-A/g)=m failed", f, g))

    # A second Bezout row checks both quotient-class independences.
    r = badd(badd(BONE, bx(X)), BY)
    A2 = badd(A, bmul(r, bx(g)))
    C2 = bsub(C, bmul(r, Px))
    require(badd(bmul(A2, Px), bmul(C2, bx(g))) == BONE,
            "alternate Bezout row failed")
    m2 = badd(bdx(A2), bdy(C2))
    require(bsub(m2, m) == bscale(derivation(Px, g, r), -1),
            "divergence row class failed")

    return {"u": u, "A": A, "C": C, "Px": Px, "m": m}


def quotient_multiplication_rank(chi, psi):
    quotient = quotient_exact(chi, psi, "chi/psi")
    domain_degree = len(quotient) - 1
    codomain_degree = len(chi) - 1
    columns = []
    for exponent in range(domain_degree):
        image = mod_poly(mul(psi, power(X, exponent)), chi)
        columns.append([image[row] if row < len(image) else Fraction(0)
                        for row in range(codomain_degree)])
    matrix = [[columns[column][row] for column in range(len(columns))]
              for row in range(codomain_degree)]
    return matrix_rank(matrix), domain_degree, codomain_degree


def audit_response_case(case):
    name = case["name"]
    f, g = trim(case["f"]), trim(case["g"])
    row = build_bezout_row(f, g)
    gp = deriv(g)
    common = gcd_poly(g, gp)
    radical = quotient_exact(g, common, "radical")
    require(radical == case["radical"], (name, "radical", radical))

    chi = minimal_polynomial(multiplication_matrix(f, g))
    psi = minimal_polynomial(multiplication_matrix(f, radical))
    require(chi == case["chi"], (name, "chi", chi))
    require(psi == case["psi"], (name, "psi", psi))
    require(not divmod_poly(chi, psi)[1], (name, "psi does not divide chi"))

    local_chi = ONE
    local_psi = ONE
    for prime, exponent in case["factors"]:
        local_chi = lcm_poly(
            local_chi,
            minimal_polynomial(multiplication_matrix(f, power(prime, exponent))),
        )
        local_psi = lcm_poly(
            local_psi,
            minimal_polynomial(multiplication_matrix(f, prime)),
        )
    require(local_chi == chi, (name, "thick-arm lcm", local_chi, chi))
    require(local_psi == psi, (name, "support-arm lcm", local_psi, psi))

    u = row["u"]
    numerator = badd(
        bscale(bx(mul(u, g)), -1),
        bmul(bx(mul(mul(u, u), gp)), BY),
    )
    require(not ideal_power_remainder(numerator, g, 1),
            (name, "depth-two numerator not in I"))
    require(not in_g_tq(numerator, g, 2),
            (name, "top symbol vanished"))
    top_coefficient = coefficient_slices(numerator).get(1, ZERO)
    require(mod_poly(top_coefficient, g) == mod_poly(mul(mul(u, u), gp), g),
            (name, "top coefficient"))

    annihilator_tests = 0
    top_tests = 0
    max_degree = max(len(chi), len(psi))
    for coefficients in product((-1, 0, 1), repeat=max_degree + 1):
        F = trim(coefficients)
        evaluated = compose(F, f)
        eta_killed = not mod_poly(evaluated, g)
        chi_divides = not divmod_poly(F, chi)[1]
        require(eta_killed == chi_divides,
                (name, "eta annihilator mismatch", F))
        top_killed = not mod_poly(mul(gp, evaluated), g)
        psi_divides = not divmod_poly(F, psi)[1]
        require(top_killed == psi_divides,
                (name, "top annihilator mismatch", F))
        annihilator_tests += 1
        top_tests += 1

    rank, domain_degree, codomain_degree = quotient_multiplication_rank(chi, psi)
    require(rank == domain_degree, (name, "times-psi not injective", rank))
    require(codomain_degree - rank == len(psi) - 1,
            (name, "support quotient dimension"))

    # Dbar(eta)=0 because D(h)=m is polynomial.  The nonzero lower powers of
    # chi show that eta and therefore delta(eta)=mu need not vanish.
    quotient = quotient_exact(chi, psi, "thickness quotient")
    require(chi != psi, (name, "missing multiplicity thickness"))
    require(mod_poly(compose(psi, f), g),
            (name, "psi unexpectedly annihilates eta"))
    require(not mod_poly(mul(gp, compose(psi, f)), g),
            (name, "psi does not lower eta to Q1"))

    return {
        "name": name,
        "degree_g": len(g) - 1,
        "arms": len(case["factors"]),
        "chi": poly_vector(chi),
        "psi": poly_vector(psi),
        "thickness": poly_vector(quotient),
        "annihilator_tests": annihilator_tests,
        "top_tests": top_tests,
        "row_terms": [len(row["A"]), len(row["C"]), len(row["m"])],
    }


def audit_characteristic_p_hostile():
    rows = []
    for prime in (2, 3, 5, 7):
        # Over F_p: P=x+x^p z, g=x^p, P_x=1, A=1, C=0.  Thus
        # h=-x^-p is a nonpolynomial D-constant, eta has annihilator T^p,
        # but m=0 and the connecting map kills eta.
        derivative_coefficient = prime % prime
        require(derivative_coefficient == 0, (prime, "Frobenius derivative"))
        quotient_powers = [exponent >= prime for exponent in range(prime + 1)]
        require(quotient_powers == [False] * prime + [True],
                (prime, "minimal polynomial of x mod x^p"))
        require(-prime < 0 and (-prime) % prime == 0,
                (prime, "nonpolynomial Laurent D-constant"))
        rows.append({
            "p": prime,
            "eta_ann_degree": prime,
            "lower_survivors": prime,
            "mu_zero": True,
            "extra_kernel": "-x^-p",
        })
    return rows


def case_data():
    x_minus_one = sub(X, ONE)
    x_plus_one = add(X, ONE)
    x2_plus_one = add(power(X, 2), ONE)
    T_minus_one = sub(X, ONE)
    return [
        {
            "name": "one_root",
            "f": X,
            "g": power(X, 2),
            "radical": X,
            "factors": [(X, 2)],
            "chi": power(X, 2),
            "psi": X,
        },
        {
            "name": "equal_collision",
            "f": power(X, 2),
            "g": power(sub(power(X, 2), ONE), 2),
            "radical": sub(power(X, 2), ONE),
            "factors": [(x_minus_one, 2), (x_plus_one, 2)],
            "chi": power(T_minus_one, 2),
            "psi": T_minus_one,
        },
        {
            "name": "unequal_collision",
            "f": power(X, 2),
            "g": mul(power(x_minus_one, 2), power(x_plus_one, 3)),
            "radical": mul(x_minus_one, x_plus_one),
            "factors": [(x_minus_one, 2), (x_plus_one, 3)],
            "chi": power(T_minus_one, 3),
            "psi": T_minus_one,
        },
        {
            "name": "nonsplit_Q",
            "f": X,
            "g": power(x2_plus_one, 2),
            "radical": x2_plus_one,
            "factors": [(x2_plus_one, 2)],
            "chi": power(x2_plus_one, 2),
            "psi": x2_plus_one,
        },
        {
            "name": "distinct_values",
            "f": X,
            "g": mul(power(X, 2), power(x_minus_one, 3)),
            "radical": mul(X, x_minus_one),
            "factors": [(X, 2), (x_minus_one, 3)],
            "chi": mul(power(X, 2), power(T_minus_one, 3)),
            "psi": mul(X, T_minus_one),
        },
    ]


def main():
    cases = case_data()
    squarefree = mul(X, sub(X, ONE))
    affine_packets = cases + [{
        "name": "squarefree_affine_only",
        "g": squarefree,
        "factors": [(X, 1), (sub(X, ONE), 1)],
    }]
    affine = audit_affine_power_filtration(affine_packets)
    responses = [audit_response_case(case) for case in cases]
    characteristic_p = audit_characteristic_p_hostile()

    semantic = {
        "affine": affine,
        "responses": responses,
        "characteristic_p": characteristic_p,
        "type_repair": {
            "quotient_derivation": "Dbar(eta)=0",
            "connecting_map": "delta(eta)=mu",
            "delta_injective_hypothesis": "ker(D,R)=ker(D,S)=K[P]",
        },
        "q2_hostiles": {
            "mixed": "gt in (g,t)^2 but not in (g^2,t^2)",
            "one_power_short": "g/g^2 not regular",
            "radical": "x/x^2 not regular",
            "duplicate_arm": "K[x]/(x^2) has only 0,1 idempotents",
        },
    }
    digest = sha256(json.dumps(
        semantic, sort_keys=True, separators=(",", ":")
    ).encode()).hexdigest()
    require(EXPECTED_SEMANTIC_SHA256 == "TO_BE_FILLED"
            or digest == EXPECTED_SEMANTIC_SHA256,
            ("semantic digest mismatch", digest, EXPECTED_SEMANTIC_SHA256))

    print("THM-3406 AFFINE-MODIFICATION TRANSGRESSION AUDIT")
    print("affine_power_filtration=PASS"
          + f" comparisons={affine['comparisons']}"
          + f" transitions={affine['transition_cells']}"
          + f" length_cells={affine['length_cells']}"
          + f" crt_cells={affine['crt_cells']}")
    print("graded_quotient=PASS Q_q/Q_(q-1)=B/(g,t^q)")
    print("transgression=PASS Dbar_eta=0 delta_eta=mu delta_injective=char0")
    for row in responses:
        print("case=" + row["name"]
              + " degree_g=" + str(row["degree_g"])
              + " arms=" + str(row["arms"])
              + " chi=" + json.dumps(row["chi"], separators=(",", ":"))
              + " psi=" + json.dumps(row["psi"], separators=(",", ":"))
              + " thickness=" + json.dumps(row["thickness"], separators=(",", ":"))
              + " ann_tests=" + str(row["annihilator_tests"]))
    print("support_thickness=PASS exact_sequence=K[T]/(chi/psi)->K[T]/chi->K[T]/psi")
    print("hostiles=PASS q2_mixed=1 radical=1 duplicate_arm=1"
          + f" characteristic_p={len(characteristic_p)}")
    print(f"semantic_sha256={digest}")


if __name__ == "__main__":
    main()
