#!/usr/bin/env python3
"""Exact hostile controls for THM-3348's linear-z puncture response.

Standard-library ``Fraction`` arithmetic only.  These finite checks audit
signs, residue bases, and sharp valuations; the uniform theorem is proved in
THM-3348 and is not inferred from this bank.
"""

from fractions import Fraction as F
from hashlib import sha256
from math import comb
from pathlib import Path


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def trim(poly):
    poly = [F(value) for value in poly]
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def multiply(first, second, cap=None):
    size = len(first) + len(second) - 1
    if cap is not None:
        size = min(size, cap + 1)
    result = [F(0)] * size
    for i, left in enumerate(first):
        for j, right in enumerate(second):
            if i + j < size:
                result[i + j] += left * right
    return trim(result)


def power(poly, exponent):
    result = [F(1)]
    while exponent:
        if exponent & 1:
            result = multiply(result, poly)
        poly = multiply(poly, poly)
        exponent //= 2
    return result


def derivative(poly):
    return trim([i * poly[i] for i in range(1, len(poly))] or [0])


def evaluate(poly, value):
    result = F(0)
    for coefficient in reversed(poly):
        result = result * value + coefficient
    return result


def exact_divide(dividend, divisor):
    dividend, divisor = trim(dividend), trim(divisor)
    quotient = [F(0)] * max(1, len(dividend) - len(divisor) + 1)
    remainder = dividend[:]
    while len(remainder) >= len(divisor) and remainder != [0]:
        shift = len(remainder) - len(divisor)
        coefficient = remainder[-1] / divisor[-1]
        quotient[shift] += coefficient
        for j, value in enumerate(divisor):
            remainder[shift + j] -= coefficient * value
        remainder = trim(remainder)
    require(remainder == [0], ("nonexact division", dividend, divisor, remainder))
    return trim(quotient)


def factor(root):
    return [-F(root), F(1)]


def determinant(matrix):
    matrix = [row[:] for row in matrix]
    answer = F(1)
    for column in range(len(matrix)):
        pivot = next(
            (row for row in range(column, len(matrix)) if matrix[row][column]),
            None,
        )
        require(pivot is not None, ("singular matrix", matrix))
        if pivot != column:
            matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
            answer = -answer
        value = matrix[column][column]
        answer *= value
        for j in range(column, len(matrix)):
            matrix[column][j] /= value
        for row in range(column + 1, len(matrix)):
            value = matrix[row][column]
            for j in range(column, len(matrix)):
                matrix[row][j] -= value * matrix[column][j]
    return answer


def qtext(value):
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def vector_text(vector):
    return "(" + ",".join(qtext(value) for value in vector) + ")"


def matrix_text(matrix):
    return "[" + ",".join(vector_text(row) for row in matrix) + "]"


def theta_residue(roots, multiplicities, leading, index):
    """Residue of ``-dx/g`` at one root by reciprocal power series."""
    root, order = roots[index], multiplicities[index]
    cofactor = [F(leading)]
    for other, (candidate, multiplicity) in enumerate(zip(roots, multiplicities)):
        if other == index:
            continue
        delta = root - candidate
        shifted = [
            F(comb(multiplicity, j)) * delta ** (multiplicity - j)
            for j in range(multiplicity + 1)
        ]
        cofactor = multiply(cofactor, shifted, order - 1)
    inverse = [F(1) / cofactor[0]]
    for n in range(1, order):
        inverse.append(
            -sum(
                (cofactor[k] if k < len(cofactor) else 0) * inverse[n - k]
                for k in range(1, n + 1)
            )
            / cofactor[0]
        )
    return -inverse[order - 1]


def puncture_control(name, roots, multiplicities, leading=F(1)):
    roots = tuple(map(F, roots))
    count = len(roots)
    squarefree, polynomial = [F(1)], [F(leading)]
    for root, multiplicity in zip(roots, multiplicities):
        squarefree = multiply(squarefree, factor(root))
        polynomial = multiply(polynomial, power(factor(root), multiplicity))
    quotient = exact_divide(polynomial, squarefree)
    squarefree_derivative = derivative(squarefree)
    matrix = [
        [root**j / evaluate(squarefree_derivative, root) for j in range(count)]
        for root in roots
    ]
    det = determinant(matrix)
    require(det != 0, (name, "zero residue determinant"))
    for j in range(count):
        monomial = [F(0)] * j + [F(1)]
        eta = multiply(quotient, monomial)
        require(
            multiply(eta, squarefree) == multiply(polynomial, monomial),
            (name, "eta identity", j),
        )
    residues = [
        theta_residue(roots, multiplicities, leading, index)
        for index in range(count)
    ]
    if len(polynomial) - 1 >= 2:
        require(sum(residues, F(0)) == 0, (name, "residue sum", residues))
    all_zero = all(value == 0 for value in residues)
    predicted_zero = count == 1 and multiplicities[0] >= 2
    require(all_zero == predicted_zero, (name, "theta criterion", residues))
    print(
        f"PUNCTURE {name}: r={count} det={qtext(det)} "
        f"theta_res={vector_text(residues)} M={matrix_text(matrix)}"
    )
    return name, tuple(residues), det, tuple(tuple(row) for row in matrix)


# Laurent bivariate polynomials in (u,z), keyed by exponent pairs.
def clean(poly):
    return {key: F(value) for key, value in poly.items() if value}


def bivariate_add(first, second, scale=F(1)):
    result = dict(first)
    for key, value in second.items():
        result[key] = result.get(key, F(0)) + scale * value
    return clean(result)


def bivariate_scale(poly, scalar):
    return clean({key: F(scalar) * value for key, value in poly.items()})


def bivariate_multiply(first, second):
    result = {}
    for (i, j), left in first.items():
        for (k, ell), right in second.items():
            key = i + k, j + ell
            result[key] = result.get(key, F(0)) + left * right
    return clean(result)


def bivariate_power(poly, exponent):
    result = {(0, 0): F(1)}
    while exponent:
        if exponent & 1:
            result = bivariate_multiply(result, poly)
        poly = bivariate_multiply(poly, poly)
        exponent //= 2
    return result


def bivariate_derivative(poly, coordinate):
    result = {}
    for (i, j), coefficient in poly.items():
        exponent = i if coordinate == 0 else j
        if exponent:
            key = (i - 1, j) if coordinate == 0 else (i, j - 1)
            result[key] = coefficient * exponent
    return clean(result)


def hamiltonian(p, leading, degree, q):
    p_z = {(degree, 0): F(leading)}
    return bivariate_add(
        bivariate_multiply(bivariate_derivative(p, 0), bivariate_derivative(q, 1)),
        bivariate_multiply(p_z, bivariate_derivative(q, 0)),
        F(-1),
    )


def one_root_control(degree, order, leading, initial):
    phi = {} if order is None else {
        (order, 0): F(initial),
        (order + 1, 0): F(2, 3),
        (order + 3, 0): F(-1, 5),
    }
    p = bivariate_add(phi, {(degree, 1): F(leading)})
    nu = degree if order is None else min(order, degree)
    annihilator_power = (degree - 2 + nu) // nu
    q0 = {(1 - degree, 0): F(1) / (F(leading) * (degree - 1))}
    unit = bivariate_add({(0, 0): F(1)}, bivariate_scale(p, F(2, 3)))
    unit = bivariate_add(unit, bivariate_scale(bivariate_power(p, 2), F(-1, 5)))
    for exponent in (annihilator_power - 1, annihilator_power, annihilator_power + 1):
        response = bivariate_multiply(bivariate_power(p, exponent), unit)
        primitive = bivariate_multiply(response, q0)
        require(
            hamiltonian(p, leading, degree, primitive) == response,
            ("Hamiltonian identity", degree, order, exponent),
        )
        valuation = min(i for i, _ in primitive)
        require(
            valuation == exponent * nu - (degree - 1),
            ("valuation", degree, order, exponent, valuation),
        )
        require(
            all(i >= 0 for i, _ in primitive) == (exponent >= annihilator_power),
            ("sharp polynomial threshold", degree, order, exponent),
        )
    return nu, annihilator_power


def main():
    print("EXACT CONTROL ONLY -- finite checks do not prove the uniform theorem")
    rows = (
        puncture_control("double", [0], [2]),
        puncture_control("two-double", [0, 1], [2, 2]),
        puncture_control("three-double", [0, 1, 2], [2, 2, 2]),
        puncture_control("simple-boundary", [F(2, 3)], [1], F(-5, 2)),
        puncture_control(
            "mixed", [F(-3, 2), F(1, 3), F(7, 2)], [2, 1, 4], F(3, 5)
        ),
        puncture_control("irregular-repeated", [-3, 1, 4, 7], [2, 3, 4, 2], F(-2, 7)),
    )

    cases = 0
    samples = []
    leadings = [F(1), F(-3, 2), F(5, 7)]
    initials = [F(1), F(-4, 3), F(7, 5)]
    for degree in range(2, 11):
        orders = sorted({1, 2, max(1, degree - 1), degree, degree + 2}) + [None]
        for index, order in enumerate(orders):
            nu, annihilator_power = one_root_control(
                degree,
                order,
                leadings[(degree + index) % len(leadings)],
                initials[(2 * degree + index) % len(initials)],
            )
            cases += 1
            if (degree, order) in ((2, 1), (5, 2), (7, 7), (10, None)):
                samples.append(
                    (degree, "inf" if order is None else order, nu, annihilator_power)
                )
    require(cases == 51, cases)
    print(f"ONE_ROOT samples (d,e,nu,N): {samples}")
    print(
        f"ONE_ROOT hostile bank passed: {cases} cases, d=2..10; "
        "sharp N-1 failure, N success, N+1 success; exact D(Q)=H(P)"
    )
    print(f"semantic_sha256={sha256(repr((rows, samples, cases)).encode()).hexdigest()}")
    source = Path(__file__)
    print(f"source_sha256={sha256(source.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()}")
    print("PASS")


if __name__ == "__main__":
    main()
