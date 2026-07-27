#!/usr/bin/env python3
"""Finite-field scout for the remaining degree-18 central T=0 wall.

This is an exact diagnostic, not a proof of a complex factorization.  It
specializes B=C=1, D=25/126 in the THM-2297 spectral cubic, constructs the
degree-12 branch discriminant, interpolates its repeated-root resultant as
a polynomial in W, divides the known R and S factors, and reports the
squarefree multiplicity profile of the residual factor.
"""


WEIGHT_30_MONOMIALS = [
    (i, j, k)
    for i in range(16)
    for j in range(11)
    for k in range(7)
    if 2 * i + 3 * j + 5 * k == 30
]


def trim(poly):
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def add(first, second, prime):
    out = [0] * max(len(first), len(second))
    for index, value in enumerate(first):
        out[index] = (out[index] + value) % prime
    for index, value in enumerate(second):
        out[index] = (out[index] + value) % prime
    return trim(out)


def scale(poly, scalar, prime):
    return trim([(scalar * value) % prime for value in poly])


def mul(first, second, prime):
    out = [0] * (len(first) + len(second) - 1)
    for i, left in enumerate(first):
        for j, right in enumerate(second):
            out[i + j] = (out[i + j] + left * right) % prime
    return trim(out)


def power(poly, exponent, prime):
    out = [1]
    base = poly[:]
    while exponent:
        if exponent & 1:
            out = mul(out, base, prime)
        base = mul(base, base, prime)
        exponent //= 2
    return out


def derivative(poly, prime):
    if len(poly) == 1:
        return [0]
    return trim([(index * poly[index]) % prime for index in range(1, len(poly))])


def divmod_poly(numerator, denominator, prime):
    numerator = trim(numerator[:])
    denominator = trim(denominator[:])
    if denominator == [0]:
        raise ZeroDivisionError
    if len(numerator) < len(denominator):
        return [0], numerator
    quotient = [0] * (len(numerator) - len(denominator) + 1)
    inverse = pow(denominator[-1], prime - 2, prime)
    while numerator != [0] and len(numerator) >= len(denominator):
        shift = len(numerator) - len(denominator)
        coefficient = numerator[-1] * inverse % prime
        quotient[shift] = coefficient
        for index, value in enumerate(denominator):
            numerator[index + shift] = (
                numerator[index + shift] - coefficient * value
            ) % prime
        trim(numerator)
    return trim(quotient), trim(numerator)


def exact_quotient(numerator, denominator, prime):
    quotient, remainder = divmod_poly(numerator, denominator, prime)
    if remainder != [0]:
        raise RuntimeError("polynomial division was not exact")
    return quotient


def gcd_poly(first, second, prime):
    first = trim(first[:])
    second = trim(second[:])
    while second != [0]:
        _, remainder = divmod_poly(first, second, prime)
        first, second = second, remainder
    inverse = pow(first[-1], prime - 2, prime)
    return scale(first, inverse, prime)


def resultant(first, second, prime):
    first = trim(first[:])
    second = trim(second[:])
    m = len(first) - 1
    n = len(second) - 1
    if second == [0]:
        return 0
    if n == 0:
        return pow(second[0], m, prime)
    _, remainder = divmod_poly(first, second, prime)
    if remainder == [0]:
        return 0
    r = len(remainder) - 1
    sign = prime - 1 if (m * n) & 1 else 1
    head = sign * pow(second[-1], m - r, prime) % prime
    return head * resultant(second, remainder, prime) % prime


def interpolate_consecutive(values, prime):
    """Lagrange interpolation at x=0,...,len(values)-1."""
    out = [0]
    count = len(values)
    for i, value in enumerate(values):
        basis = [1]
        denominator = 1
        for j in range(count):
            if i == j:
                continue
            basis = mul(basis, [(-j) % prime, 1], prime)
            denominator = denominator * (i - j) % prime
        coefficient = value * pow(denominator, prime - 2, prime) % prime
        out = add(out, scale(basis, coefficient, prime), prime)
    return trim(out)


def solve_square_linear(rows, values, prime):
    size = len(rows)
    augmented = [
        [entry % prime for entry in row] + [value % prime]
        for row, value in zip(rows, values, strict=True)
    ]
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if augmented[row][column]),
            None,
        )
        if pivot is None:
            return None
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        inverse = pow(augmented[column][column], prime - 2, prime)
        augmented[column] = [
            entry * inverse % prime for entry in augmented[column]
        ]
        for row in range(size):
            if row == column:
                continue
            factor = augmented[row][column]
            if factor:
                augmented[row] = [
                    (left - factor * right) % prime
                    for left, right in zip(
                        augmented[row], augmented[column], strict=True
                    )
                ]
    return [augmented[index][-1] for index in range(size)]


def recover_normalized_t(prime):
    size = len(WEIGHT_30_MONOMIALS)
    for attempt in range(50):
        rows = []
        values = []
        for offset in range(size):
            index = attempt * size + offset
            b_value = 1 + index % 11
            c_value = 1 + (3 * index + index // 11) % 13
            w_value = 1 + (5 * index * index + 7 * index + 3) % 29
            r_value = (20 * b_value * c_value + 21 * w_value) % prime
            s_value = (
                2888 * b_value**5
                + 108864 * b_value**2 * c_value**2
                + 571536 * b_value * c_value * w_value
                + 750141 * w_value**2
            ) % prime
            if not r_value or not s_value:
                break
            repeated = residual_resultant_value(
                prime, b_value, c_value, w_value
            )
            denominator = (
                pow(r_value, 6, prime) * pow(s_value, 3, prime)
            ) % prime
            t_value = repeated * pow(denominator, prime - 2, prime) % prime
            rows.append(
                [
                    pow(b_value, i, prime)
                    * pow(c_value, j, prime)
                    * pow(w_value, k, prime)
                    % prime
                    for i, j, k in WEIGHT_30_MONOMIALS
                ]
            )
            values.append(t_value)
        if len(rows) != size:
            continue
        solution = solve_square_linear(rows, values, prime)
        if solution is None or not solution[0]:
            continue
        inverse = pow(solution[0], prime - 2, prime)
        normalized = [coefficient * inverse % prime for coefficient in solution]

        # A fresh point is an exact interpolation control.
        b_value, c_value, w_value = 17, 19, 23
        r_value = (20 * b_value * c_value + 21 * w_value) % prime
        s_value = (
            2888 * b_value**5
            + 108864 * b_value**2 * c_value**2
            + 571536 * b_value * c_value * w_value
            + 750141 * w_value**2
        ) % prime
        repeated = residual_resultant_value(prime, b_value, c_value, w_value)
        denominator = pow(r_value, 6, prime) * pow(s_value, 3, prime) % prime
        raw_t = repeated * pow(denominator, prime - 2, prime) % prime
        predicted_normalized = sum(
            coefficient
            * pow(b_value, i, prime)
            * pow(c_value, j, prime)
            * pow(w_value, k, prime)
            for coefficient, (i, j, k) in zip(
                normalized, WEIGHT_30_MONOMIALS, strict=True
            )
        ) % prime
        raw_scale = solution[0]
        if predicted_normalized * raw_scale % prime != raw_t:
            raise RuntimeError("weight-thirty interpolation control failed")
        return normalized
    raise RuntimeError("failed to find a nonsingular interpolation grid")


def branch_discriminant(prime, w_value, b_value=1, c_value=1):
    d_value = 25 * b_value * b_value * pow(126, prime - 2, prime) % prime
    a = -26040609 % prime
    bpoly = [49601160 * b_value % prime, 0, 1607445 % prime]
    cpoly = [
        (-20995200 * b_value**2 - 52907904 * d_value) % prime,
        0,
        -2857680 * b_value % prime,
        0,
        -138915 % prime,
    ]
    dpoly = [
        33592320 * b_value * d_value % prime,
        (-5878656 * w_value - 5598720 * b_value * c_value) % prime,
        (777600 * b_value**2 + 1959552 * d_value) % prime,
        -435456 * c_value % prime,
        78120 * b_value % prime,
        0,
        1127 % prime,
    ]
    terms = [
        mul(power(bpoly, 2, prime), power(cpoly, 2, prime), prime),
        scale(power(cpoly, 3, prime), -4 * a, prime),
        scale(mul(power(bpoly, 3, prime), dpoly, prime), -4, prime),
        scale(power(dpoly, 2, prime), -27 * a * a, prime),
        scale(mul(mul(bpoly, cpoly, prime), dpoly, prime), 18 * a, prime),
    ]
    out = [0]
    for term in terms:
        out = add(out, term, prime)
    if len(out) - 1 != 12:
        raise RuntimeError("branch discriminant lost degree twelve")
    if out[0] != 0 or out[1] != 0:
        raise RuntimeError("central branch lost its universal y^2 factor")
    residual = trim(out[2:])
    if len(residual) - 1 != 10:
        raise RuntimeError("central residual discriminant lost degree ten")
    return residual


def residual_resultant_value(prime, b_value, c_value, w_value):
    branch = branch_discriminant(prime, w_value, b_value, c_value)
    return resultant(branch, derivative(branch, prime), prime)


def squarefree_profile(poly, prime):
    derivative_poly = derivative(poly, prime)
    repeated = gcd_poly(poly, derivative_poly, prime)
    layer = exact_quotient(poly, repeated, prime)
    profile = []
    multiplicity = 1
    while layer != [1]:
        overlap = gcd_poly(layer, repeated, prime)
        exact = exact_quotient(layer, overlap, prime)
        if exact != [1]:
            profile.append((len(exact) - 1, multiplicity))
        layer = overlap
        repeated = exact_quotient(repeated, overlap, prime)
        multiplicity += 1
    if repeated != [1]:
        profile.append((len(repeated) - 1, f">={multiplicity} or pth-power"))
    return profile


def audit_prime(prime):
    # After removing the universal y^2, the degree-ten residual resultant
    # has weight 10*9=90 and W has weight five.
    degree_bound = 18
    values = []
    for w_value in range(degree_bound + 1):
        branch = branch_discriminant(prime, w_value)
        values.append(resultant(branch, derivative(branch, prime), prime))
    repeated_resultant = interpolate_consecutive(values, prime)

    r_factor = [20 % prime, 21 % prime]
    s_factor = [
        (2888 + 108864) % prime,
        571536 % prime,
        750141 % prime,
    ]
    residual = repeated_resultant
    r_order = 0
    s_order = 0
    while True:
        quotient, remainder = divmod_poly(residual, r_factor, prime)
        if remainder != [0]:
            break
        residual = quotient
        r_order += 1
    while True:
        quotient, remainder = divmod_poly(residual, s_factor, prime)
        if remainder != [0]:
            break
        residual = quotient
        s_order += 1
    return {
        "prime": prime,
        "resultant_degree_W": len(repeated_resultant) - 1,
        "R_order": r_order,
        "S_order": s_order,
        "residual_degree_W": len(residual) - 1,
        "residual_squarefree_profile": squarefree_profile(residual, prime),
    }


def main():
    reports = [audit_prime(prime) for prime in (1000003, 1000033, 1000037)]
    if not all(report == {**reports[0], "prime": report["prime"]} for report in reports):
        raise RuntimeError("finite-field profiles disagree")
    print("JC2 DEGREE-18 CENTRAL T=0 MODULAR SCOUT")
    for report in reports:
        print(report)
    print("Exact finite-field diagnostic only; no C-factorization claim.")
    print("PASS")


if __name__ == "__main__":
    main()
