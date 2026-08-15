#!/usr/bin/env python3
"""Exact hostile audit for THM-3386's linear-z divergence torsion law."""

from fractions import Fraction
from hashlib import sha256
from itertools import permutations
import json


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


def trim(poly):
    values = [Fraction(value) for value in poly]
    while values and not values[-1]:
        values.pop()
    return tuple(values)


UZERO = ()
UONE = (Fraction(1),)
UX = (Fraction(0), Fraction(1))


def uadd(left, right):
    size = max(len(left), len(right))
    return trim([
        (left[i] if i < len(left) else 0)
        + (right[i] if i < len(right) else 0)
        for i in range(size)
    ])


def uscale(poly, scalar):
    scalar = Fraction(scalar)
    return trim([scalar * value for value in poly])


def usub(left, right):
    return uadd(left, uscale(right, -1))


def umul(left, right):
    if not left or not right:
        return UZERO
    answer = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            answer[i + j] += left_value * right_value
    return trim(answer)


def upow(poly, exponent):
    require(exponent >= 0, ("negative univariate power", exponent))
    answer = UONE
    base = poly
    current = exponent
    while current:
        if current & 1:
            answer = umul(answer, base)
        base = umul(base, base)
        current //= 2
    return answer


def uderiv(poly):
    return trim([i * poly[i] for i in range(1, len(poly))])


def udivmod(dividend, divisor):
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


def umod(poly, modulus):
    return udivmod(poly, modulus)[1]


def uxgcd(left, right):
    old_r, current_r = trim(left), trim(right)
    old_s, current_s = UONE, UZERO
    old_t, current_t = UZERO, UONE
    while current_r:
        quotient, remainder = udivmod(old_r, current_r)
        old_r, current_r = current_r, remainder
        old_s, current_s = current_s, usub(old_s, umul(quotient, current_s))
        old_t, current_t = current_t, usub(old_t, umul(quotient, current_t))
    require(old_r, "gcd vanished")
    factor = 1 / old_r[-1]
    return uscale(old_r, factor), uscale(old_s, factor), uscale(old_t, factor)


def ucompose(poly, argument):
    answer = UZERO
    for coefficient in reversed(poly):
        answer = uadd(umul(answer, argument), (coefficient,))
    return answer


def coeff_text(value):
    value = Fraction(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def poly_vector(poly):
    return [coeff_text(value) for value in trim(poly)]


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


def bmul(left, right):
    answer = {}
    for (left_x, left_z), left_value in left.items():
        for (right_x, right_z), right_value in right.items():
            key = (left_x + right_x, left_z + right_z)
            answer[key] = answer.get(key, Fraction(0)) + left_value * right_value
    return bclean(answer)


def bpow(poly, exponent):
    require(exponent >= 0, ("negative bivariate power", exponent))
    answer = {(0, 0): Fraction(1)}
    base = poly
    current = exponent
    while current:
        if current & 1:
            answer = bmul(answer, base)
        base = bmul(base, base)
        current //= 2
    return answer


def bdx(poly):
    return bclean({(ix - 1, iz): ix * value
                   for (ix, iz), value in poly.items() if ix})


def bdz(poly):
    return bclean({(ix, iz - 1): iz * value
                   for (ix, iz), value in poly.items() if iz})


def bx(poly):
    return bclean({(i, 0): value for i, value in enumerate(poly)})


BONE = {(0, 0): Fraction(1)}
BZ = {(0, 1): Fraction(1)}


def bcompose(poly, argument):
    answer = {}
    for coefficient in reversed(poly):
        answer = badd(bmul(answer, argument), {(0, 0): coefficient})
    return answer


def bdivide_x_exact(poly, divisor):
    by_z = {}
    for (ix, iz), value in poly.items():
        coefficients = by_z.setdefault(iz, [])
        if len(coefficients) <= ix:
            coefficients.extend([Fraction(0)] * (ix + 1 - len(coefficients)))
        coefficients[ix] += value
    answer = {}
    for iz, coefficients in by_z.items():
        quotient, remainder = udivmod(trim(coefficients), divisor)
        require(not remainder, ("nonexact x division", iz, remainder))
        for ix, value in enumerate(quotient):
            if value:
                answer[(ix, iz)] = value
    return bclean(answer)


def derivation(Px, g, poly):
    return badd(bmul(Px, bdz(poly)), bscale(bmul(bx(g), bdx(poly)), -1))


def multiplication_matrix(f, g):
    degree = len(g) - 1
    require(degree >= 1 and g[-1] == 1, ("modulus not monic", g))
    matrix = [[Fraction(0) for _ in range(degree)] for _ in range(degree)]
    x_power = UONE
    for column in range(degree):
        remainder = umod(umul(f, x_power), g)
        for row, value in enumerate(remainder):
            matrix[row][column] = value
        x_power = umul(x_power, UX)
    return matrix


def matrix_identity(size):
    return [[Fraction(i == j) for j in range(size)] for i in range(size)]


def matrix_mul(left, right):
    rows = len(left)
    inner = len(right)
    columns = len(right[0])
    return [[sum(left[i][k] * right[k][j] for k in range(inner))
             for j in range(columns)] for i in range(rows)]


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
        scale = rows[pivot_row][column]
        rows[pivot_row] = [value / scale for value in rows[pivot_row]]
        for row in range(len(rows)):
            if row == pivot_row or not rows[row][column]:
                continue
            factor = rows[row][column]
            rows[row] = [rows[row][j] - factor * rows[pivot_row][j]
                         for j in range(variable_count + 1)]
        pivot_columns.append(column)
        pivot_row += 1
    for row in rows:
        if all(not row[column] for column in range(variable_count)) and row[-1]:
            return None
    require(len(pivot_columns) == variable_count,
            ("dependent prior powers", len(pivot_columns), variable_count))
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
            [matrix_flat(power) for power in powers],
            [-value for value in matrix_flat(current)],
        )
        if solution is not None:
            return trim(solution + [Fraction(1)])
        powers.append(current)
    raise RuntimeError("minimal polynomial not found")


def permutation_sign(permutation):
    inversions = sum(permutation[i] > permutation[j]
                     for i in range(len(permutation))
                     for j in range(i + 1, len(permutation)))
    return -1 if inversions % 2 else 1


def characteristic_polynomial(matrix):
    size = len(matrix)
    answer = UZERO
    for permutation in permutations(range(size)):
        term = UONE
        for row, column in enumerate(permutation):
            entry = (-matrix[row][column], Fraction(row == column))
            term = umul(term, trim(entry))
        answer = uadd(answer, uscale(term, permutation_sign(permutation)))
    return answer


def build_row(f, g):
    fp = uderiv(f)
    gp = uderiv(g)
    common, U, V = uxgcd(fp, g)
    require(common == UONE, ("gradient gcd", common))
    gp_square_quotient, remainder = udivmod(umul(gp, gp), g)
    require(not remainder, ("root multiplicity gate", g, remainder))

    U2 = umul(U, U)
    A = badd(bx(U), bscale(bmul(bx(umul(U2, gp)), BZ), -1))
    B = badd(
        badd(bx(V), bscale(bmul(bx(umul(umul(U, V), gp)), BZ), -1)),
        bmul(bx(umul(U2, gp_square_quotient)), bpow(BZ, 2)),
    )
    Px = badd(bx(fp), bmul(bx(gp), BZ))
    require(badd(bmul(A, Px), bmul(B, bx(g))) == BONE,
            ("Bezout row", f, g))
    m = badd(bdx(A), bdz(B))
    cleared = badd(
        badd(bscale(bmul(Px, bdz(A)), -1), bmul(bx(g), bdx(A))),
        bscale(bmul(A, bx(gp)), -1),
    )
    require(cleared == bmul(bx(g), m), ("localized primitive", f, g))
    return A, B, Px, m


def case_data():
    T_minus_one = usub(UX, UONE)
    T2_plus_one = uadd(upow(UX, 2), UONE)
    x_minus_one = usub(UX, UONE)
    x_plus_one = uadd(UX, UONE)
    return [
        {"name": "constant", "f": uadd(upow(UX, 2), UONE), "g": UONE,
         "chi": UONE, "char": UONE, "lower": []},
        {"name": "one_root", "f": UX, "g": upow(UX, 2),
         "chi": upow(UX, 2), "char": upow(UX, 2), "lower": [UX]},
        {"name": "equal_collision", "f": upow(UX, 2),
         "g": upow(usub(upow(UX, 2), UONE), 2),
         "chi": upow(T_minus_one, 2), "char": upow(T_minus_one, 4),
         "lower": [T_minus_one]},
        {"name": "unequal_collision", "f": upow(UX, 2),
         "g": umul(upow(x_minus_one, 2), upow(x_plus_one, 3)),
         "chi": upow(T_minus_one, 3), "char": upow(T_minus_one, 5),
         "lower": [T_minus_one, upow(T_minus_one, 2)]},
        {"name": "nonsplit_Q", "f": UX, "g": upow(T2_plus_one, 2),
         "chi": upow(T2_plus_one, 2), "char": upow(T2_plus_one, 2),
         "lower": [T2_plus_one]},
        {"name": "distinct_values", "f": UX,
         "g": umul(upow(UX, 2), upow(x_minus_one, 3)),
         "chi": umul(upow(UX, 2), upow(T_minus_one, 3)),
         "char": umul(upow(UX, 2), upow(T_minus_one, 3)),
         "lower": [umul(UX, T_minus_one)]},
    ]


def audit_case(case):
    name, f, g = case["name"], trim(case["f"]), trim(case["g"])
    require(g and g[-1] == 1, (name, "nonmonic g", g))
    if len(g) == 1:
        require(case["chi"] == UONE, (name, "constant chi"))
        return {"name": name, "degree_g": 0, "chi": ["1"],
                "char": ["1"], "lower_hostiles": 0}

    A, B, Px, m = build_row(f, g)
    P = badd(bx(f), bmul(bx(g), BZ))
    matrix = multiplication_matrix(f, g)
    chi = minimal_polynomial(matrix)
    char = characteristic_polynomial(matrix)
    require(chi == trim(case["chi"]), (name, "minimal polynomial", chi))
    require(char == trim(case["char"]), (name, "characteristic", char))

    evaluated = ucompose(chi, f)
    _, remainder = udivmod(evaluated, g)
    require(not remainder, (name, "evaluation kernel", remainder))
    FP = bcompose(chi, P)
    numerator = bscale(bmul(A, FP), -1)
    q = bdivide_x_exact(numerator, g)
    require(derivation(Px, g, q) == bmul(FP, m),
            (name, "killing primitive"))

    for lower in case["lower"]:
        _, lower_remainder = udivmod(ucompose(lower, f), g)
        require(lower_remainder, (name, "lower power unexpectedly kills", lower))

    return {
        "name": name,
        "degree_g": len(g) - 1,
        "chi": poly_vector(chi),
        "char": poly_vector(char),
        "lower_hostiles": len(case["lower"]),
        "row_terms": [len(A), len(B), len(m)],
    }


def main():
    rows = [audit_case(case) for case in case_data()]
    semantic = {"cases": rows, "case_count": len(rows)}
    digest = sha256(json.dumps(
        semantic, sort_keys=True, separators=(",", ":")
    ).encode()).hexdigest()

    print("THM-3386 LINEAR-Z CANONICAL DIVERGENCE AUDIT")
    print(f"cases=PASS count={len(rows)}")
    for row in rows:
        print("case=" + row["name"]
              + " degree_g=" + str(row["degree_g"])
              + " chi=" + json.dumps(row["chi"], separators=(",", ":"))
              + " char=" + json.dumps(row["char"], separators=(",", ":"))
              + " lower_hostiles=" + str(row["lower_hostiles"]))
    print("annihilator=PASS evaluation_kernel=minpoly_of_f_mod_g")
    print("collision_law=PASS coincident_values_use_max_not_sum")
    print(f"semantic_sha256={digest}")


if __name__ == "__main__":
    main()
