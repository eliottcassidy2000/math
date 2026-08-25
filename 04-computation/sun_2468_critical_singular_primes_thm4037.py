from __future__ import annotations

from fractions import Fraction
from itertools import permutations, product
from math import comb, factorial


TARGET = 896_315_812_331_399


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


# Polynomials are coefficient lists in ascending order.
def trim(poly):
    answer = list(poly)
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return answer


def add(left, right):
    size = max(len(left), len(right))
    answer = [left[0] * 0 for _ in range(size)]
    for i in range(size):
        if i < len(left):
            answer[i] += left[i]
        if i < len(right):
            answer[i] += right[i]
    return trim(answer)


def neg(poly):
    return [-coefficient for coefficient in poly]


def sub(left, right):
    return add(left, neg(right))


def mul(left, right):
    answer = [left[0] * right[0] * 0 for _ in range(len(left) + len(right) - 1)]
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            answer[i + j] += a * b
    return trim(answer)


def scale(poly, scalar):
    return trim([scalar * coefficient for coefficient in poly])


def power(poly, exponent: int):
    answer = [poly[0] * 0 + 1]
    base = list(poly)
    while exponent:
        if exponent & 1:
            answer = mul(answer, base)
        base = mul(base, base)
        exponent >>= 1
    return answer


def derivative(poly):
    if len(poly) == 1:
        return [poly[0] * 0]
    return trim([i * poly[i] for i in range(1, len(poly))])


def evaluate(poly, value):
    answer = poly[-1] * 0
    for coefficient in reversed(poly):
        answer = answer * value + coefficient
    return answer


def compose(outer, inner):
    answer = [inner[0] * 0]
    for coefficient in reversed(outer):
        answer = add(mul(answer, inner), [coefficient])
    return trim(answer)


def divmod_field(numerator, denominator):
    numerator = trim(numerator)
    denominator = trim(denominator)
    require(denominator != [0], "polynomial division by zero")
    if len(numerator) < len(denominator):
        return [numerator[0] * 0], numerator
    quotient = [numerator[0] * 0 for _ in range(len(numerator) - len(denominator) + 1)]
    while len(numerator) >= len(denominator) and numerator != [0]:
        shift = len(numerator) - len(denominator)
        coefficient = numerator[-1] / denominator[-1]
        quotient[shift] += coefficient
        subtractor = [numerator[0] * 0] * shift + scale(denominator, coefficient)
        numerator = sub(numerator, subtractor)
    return trim(quotient), trim(numerator)


def monic(poly):
    poly = trim(poly)
    return scale(poly, 1 / poly[-1])


def q_numerator(r: int):
    answer = [1]
    for odd in range(1, 2 * r, 2):
        answer = mul(answer, [-odd * odd, 1])
    return answer


def q_denominator(r: int) -> int:
    return 2 ** (2 * r) * factorial(2 * r)


CRITICAL_VALUE_POLYNOMIALS = {
    1: [Fraction(1, 8), Fraction(1)],
    2: mul(
        [Fraction(1, 24), Fraction(1)],
        [Fraction(-3, 128), Fraction(1)],
    ),
    3: mul(
        [Fraction(5, 1024), Fraction(1)],
        [Fraction(-1, 6075), Fraction(4, 243), Fraction(1)],
    ),
    4: mul(
        [Fraction(-35, 32768), Fraction(1)],
        [
            Fraction(-1, 14049280),
            Fraction(-3, 89600),
            Fraction(9, 640),
            Fraction(1),
        ],
    ),
}


def check_critical_value_polynomials() -> int:
    checks = 0
    for r in range(1, 5):
        q = [Fraction(value, q_denominator(r)) for value in q_numerator(r)]
        critical_squared = mul([Fraction(0), Fraction(1)], derivative(q))
        substituted = compose(CRITICAL_VALUE_POLYNOMIALS[r], q)
        _, remainder = divmod_field(substituted, critical_squared)
        require(remainder == [Fraction(0)], f"V_{r}(Q_{r}) critical remainder")
        require(
            len(CRITICAL_VALUE_POLYNOMIALS[r]) - 1 == r,
            f"V_{r} degree",
        )
        require(len(critical_squared) - 1 == r, f"critical polynomial {r} degree")
        checks += 3
    return checks


# Determinants over Q[T], used only for Sylvester matrices of size at most 5.
def permutation_sign(permutation) -> int:
    inversions = 0
    for i in range(len(permutation)):
        for j in range(i + 1, len(permutation)):
            inversions += permutation[i] > permutation[j]
    return -1 if inversions % 2 else 1


def determinant_poly(matrix):
    size = len(matrix)
    answer = [Fraction(0)]
    for permutation in permutations(range(size)):
        term = [Fraction(permutation_sign(permutation))]
        for row, column in enumerate(permutation):
            term = mul(term, matrix[row][column])
        answer = add(answer, term)
    return trim(answer)


def resultant_poly_t(left_y, right_y_with_t):
    m = len(left_y) - 1
    n = len(right_y_with_t) - 1
    size = m + n
    zero = [Fraction(0)]
    matrix = [[zero for _ in range(size)] for _ in range(size)]
    left_descending = [[coefficient] for coefficient in reversed(left_y)]
    right_descending = list(reversed(right_y_with_t))
    for row in range(n):
        for j, coefficient in enumerate(left_descending):
            matrix[row][row + j] = coefficient
    for row in range(m):
        for j, coefficient in enumerate(right_descending):
            matrix[n + row][row + j] = coefficient
    return determinant_poly(matrix)


def shifted_right_polynomial(right_z, shift):
    # Substitute z=T-shift-y and return coefficients in y, each in Q[T].
    t_minus_shift = [Fraction(-shift), Fraction(1)]
    answer = [[Fraction(0)] for _ in range(len(right_z))]
    for j, coefficient in enumerate(right_z):
        for y_degree in range(j + 1):
            contribution = scale(
                power(t_minus_shift, j - y_degree),
                coefficient * comb(j, y_degree) * ((-1) ** y_degree),
            )
            answer[y_degree] = add(answer[y_degree], contribution)
    return [trim(coefficient) for coefficient in answer]


Y3_BRANCHES = {
    "y0": [Fraction(5, 1024), Fraction(1)],
    "yalg": [Fraction(-1, 6075), Fraction(4, 243), Fraction(1)],
}

Z4_BRANCHES = {
    "z0": [Fraction(-35, 32768), Fraction(1)],
    "zalg": [
        Fraction(-1, 14049280),
        Fraction(-3, 89600),
        Fraction(9, 640),
        Fraction(1),
    ],
}

X4_VALUES = {
    "x5": Fraction(-1, 24),
    "x0": Fraction(3, 128),
}

A1 = Fraction(-1, 8)


BRANCH_POLYNOMIALS = {
    ("x5", "y0", "z0"): [Fraction(16759, 98304), Fraction(1)],
    ("x5", "y0", "zalg"): [
        Fraction(1356500621657, 248598075801600),
        Fraction(51240223, 550502400),
        Fraction(2707, 5120),
        Fraction(1),
    ],
    ("x5", "yalg", "z0"): [
        Fraction(586758629753, 19568944742400),
        Fraction(1384135, 3981312),
        Fraction(1),
    ],
    ("x5", "yalg", "zalg"): [
        Fraction(1559088324932338366889, 47793709368329895936000000),
        Fraction(745140272406905311, 677348488780185600000),
        Fraction(4976475632582171, 322546899419136000),
        Fraction(2322164773226827, 20159181213696000),
        Fraction(27250459141, 56435097600),
        Fraction(27929, 25920),
        Fraction(1),
    ],
    ("x0", "y0", "z0"): [Fraction(3453, 32768), Fraction(1)],
    ("x0", "y0", "zalg"): [
        Fraction(12538469091, 9207336140800),
        Fraction(6780741, 183500800),
        Fraction(1707, 5120),
        Fraction(1),
    ],
    ("x0", "yalg", "z0"): [
        Fraction(8399238139, 724775731200),
        Fraction(865735, 3981312),
        Fraction(1),
    ],
    ("x0", "yalg", "zalg"): [
        Fraction(20304796494162528520321, 9667310299885395247104000000),
        Fraction(34731461719314266609, 308268823338182246400000),
        Fraction(137969412357579659, 55048004167532544000),
        Fraction(9546566984049607, 322546899419136000),
        Fraction(5521913633, 28217548800),
        Fraction(4451, 6480),
        Fraction(1),
    ],
}


PARTIAL_CANDIDATE_FACTORS = {
    ("x5", "y0", "z0"): {5: 1, 499: 1, 1753: 1, 21617: 1, 931932289: 1},
    ("x5", "y0", "zalg"): {73: 1, 34667: 1},
    ("x5", "yalg", "z0"): {307: 1, 607: 1},
    ("x5", "yalg", "zalg"): {31: 2, 5503241: 1},
    ("x0", "y0", "z0"): {5: 1, 7: 1, 42667: 1, 199447: 1, 98610539: 1},
    ("x0", "y0", "zalg"): {13: 1, 3221: 1, 35837: 1, 508471: 1},
    ("x0", "yalg", "z0"): {7: 2, 109: 1},
    ("x0", "yalg", "zalg"): {31: 1, 401: 1, 479: 1},
}


def valuation(number: int, prime: int) -> int:
    exponent = 0
    while number % prime == 0:
        number //= prime
        exponent += 1
    return exponent


def is_prime_64(number: int) -> bool:
    if number < 2:
        return False
    small_primes = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for prime in small_primes:
        if number % prime == 0:
            return number == prime
    d = number - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for base in (2, 325, 9375, 28178, 450775, 9780504, 1795265022):
        if base % number == 0:
            continue
        value = pow(base, d, number)
        if value in (1, number - 1):
            continue
        for _ in range(s - 1):
            value = value * value % number
            if value == number - 1:
                break
        else:
            return False
    return True


def check_additive_resultant_branches():
    checks = 0
    numerators = {}
    for key, expected in BRANCH_POLYNOMIALS.items():
        x_name, y_name, z_name = key
        shift = A1 + X4_VALUES[x_name]
        right = shifted_right_polynomial(Z4_BRANCHES[z_name], shift)
        computed = monic(resultant_poly_t(Y3_BRANCHES[y_name], right))
        require(computed == expected, f"additive resultant branch {key}")
        value = evaluate(expected, Fraction(TARGET))
        numerator = abs(value.numerator)
        numerators[key] = numerator
        for prime, expected_exponent in PARTIAL_CANDIDATE_FACTORS[key].items():
            require(is_prime_64(prime), f"candidate factor {prime} is prime")
            require(
                valuation(numerator, prime) == expected_exponent,
                f"candidate valuation {prime} in {key}",
            )
            checks += 2
        checks += 1
    return checks, numerators


# Finite-field polynomial arithmetic.
def mod_poly(poly, prime):
    return trim([int(coefficient) % prime for coefficient in poly])


def divmod_mod(numerator, denominator, prime):
    numerator = mod_poly(numerator, prime)
    denominator = mod_poly(denominator, prime)
    require(denominator != [0], "modular polynomial division by zero")
    if len(numerator) < len(denominator):
        return [0], numerator
    quotient = [0] * (len(numerator) - len(denominator) + 1)
    inverse = pow(denominator[-1], -1, prime)
    while len(numerator) >= len(denominator) and numerator != [0]:
        shift = len(numerator) - len(denominator)
        coefficient = numerator[-1] * inverse % prime
        quotient[shift] = (quotient[shift] + coefficient) % prime
        for i, value in enumerate(denominator):
            numerator[shift + i] = (numerator[shift + i] - coefficient * value) % prime
        numerator = trim(numerator)
    return trim(quotient), trim(numerator)


def mul_mod(left, right, modulus, prime):
    _, remainder = divmod_mod(mul(left, right), modulus, prime)
    return remainder


def pow_mod_poly(base, exponent, modulus, prime):
    answer = [1]
    base = divmod_mod(base, modulus, prime)[1]
    while exponent:
        if exponent & 1:
            answer = mul_mod(answer, base, modulus, prime)
        base = mul_mod(base, base, modulus, prime)
        exponent >>= 1
    return answer


def gcd_mod(left, right, prime):
    left = mod_poly(left, prime)
    right = mod_poly(right, prime)
    while right != [0]:
        left, right = right, divmod_mod(left, right, prime)[1]
    if left == [0]:
        return [0]
    return scale(left, pow(int(left[-1]), -1, prime))


def split_linear_polynomial(poly, prime):
    poly = mod_poly(poly, prime)
    if len(poly) == 2:
        return [(-poly[0] * pow(poly[1], -1, prime)) % prime]
    for shift in range(128):
        base = [shift % prime, 1]
        half_power = pow_mod_poly(base, (prime - 1) // 2, poly, prime)
        factor = gcd_mod(poly, sub(half_power, [1]), prime)
        degree = len(factor) - 1
        if 0 < degree < len(poly) - 1:
            quotient, remainder = divmod_mod(poly, factor, prime)
            require(remainder == [0], "finite-field split remainder")
            return sorted(
                split_linear_polynomial(factor, prime)
                + split_linear_polynomial(quotient, prime)
            )
    raise RuntimeError(f"failed deterministic linear split mod {prime}: {poly}")


def roots_in_prime_field(poly, prime):
    # gcd(poly, X^p-X) removes non-rational and repeated roots.
    x = [0, 1]
    xp_minus_x = sub(pow_mod_poly(x, prime, poly, prime), x)
    rational_part = gcd_mod(poly, xp_minus_x, prime)
    if rational_part == [1]:
        return []
    roots = split_linear_polynomial(rational_part, prime)
    require(
        all(evaluate(poly, root) % prime == 0 for root in roots),
        f"root substitution mod {prime}",
    )
    require(len(roots) == len(rational_part) - 1, f"root completeness mod {prime}")
    return sorted(set(roots))


def square_root_mod_prime(value: int, prime: int):
    value %= prime
    if value == 0:
        return 0
    require(pow(value, (prime - 1) // 2, prime) == 1, "Tonelli input nonsquare")
    if prime % 4 == 3:
        root = pow(value, (prime + 1) // 4, prime)
    else:
        q = prime - 1
        s = 0
        while q % 2 == 0:
            q //= 2
            s += 1
        nonresidue = 2
        while pow(nonresidue, (prime - 1) // 2, prime) != prime - 1:
            nonresidue += 1
        c = pow(nonresidue, q, prime)
        root = pow(value, (q + 1) // 2, prime)
        t = pow(value, q, prime)
        m = s
        while t != 1:
            i = 1
            t2 = t * t % prime
            while t2 != 1:
                t2 = t2 * t2 % prime
                i += 1
                require(i < m, "Tonelli-Shanks loop")
            b = pow(c, 1 << (m - i - 1), prime)
            root = root * b % prime
            t = t * b * b % prime
            c = b * b % prime
            m = i
    require(root * root % prime == value, "Tonelli-Shanks output")
    return root


def critical_squared_polynomial(r: int):
    return mul([0, 1], derivative(q_numerator(r)))


def critical_role_data(prime: int, r: int):
    q = q_numerator(r)
    roots = roots_in_prime_field(critical_squared_polynomial(r), prime)
    denominator_inverse = pow(q_denominator(r), -1, prime)
    data = []
    for squared in roots:
        if squared == 0:
            u_roots = (0,)
        elif pow(squared, (prime - 1) // 2, prime) == 1:
            root = square_root_mod_prime(squared, prime)
            u_roots = tuple(sorted({root, (-root) % prime}))
        else:
            continue
        value = evaluate(q, squared) * denominator_inverse % prime
        data.append((squared, value, u_roots))
    return data


CANDIDATE_PRIMES = sorted(
    {
        prime
        for factors in PARTIAL_CANDIDATE_FACTORS.values()
        for prime in factors
        if prime > 8
    }
)

EXPECTED_SINGULAR_COUNTS = {
    31: 4,
    499: 2,
    35837: 2,
    42667: 1,
    199447: 1,
    508471: 2,
    931932289: 2,
    98610539: 1,
}

EXPECTED_FIRST_RESIDUALS = {
    31: [7, 7, 7, 7],
    499: [356, 356],
    35837: [27183, 27183],
    42667: [25517],
    199447: [69875],
    508471: [151075, 151075],
    931932289: [529232077, 529232077],
    98610539: [70179963],
}


def binomial_derivative_mod(top: int, k: int, prime: int) -> int:
    # p>k makes the rational polynomial derivative denominator invertible.
    numerator_poly = [1]
    for j in range(k):
        numerator_poly = mul(numerator_poly, [-j, 1])
    derivative_value = evaluate(derivative(numerator_poly), top) % prime
    return derivative_value * pow(factorial(k), -1, prime) % prime


def singular_points_and_nonlifts(prime: int):
    roles = [critical_role_data(prime, r) for r in range(1, 5)]
    points = []
    for squared_combo in product(*roles):
        if sum(item[1] for item in squared_combo) % prime != TARGET % prime:
            continue
        for u_combo in product(*(item[2] for item in squared_combo)):
            top_indices = tuple(
                ((u + (2 * r - 1)) * pow(2, -1, prime)) % prime
                for r, u in enumerate(u_combo, 1)
            )
            exact_sum = sum(comb(top, 2 * r) for r, top in enumerate(top_indices, 1))
            require((exact_sum - TARGET) % prime == 0, f"singular target mod {prime}")
            for r, top in enumerate(top_indices, 1):
                require(
                    binomial_derivative_mod(top, 2 * r, prime) == 0,
                    f"gradient coordinate r={r} mod {prime}",
                )
            first_residual = ((exact_sum - TARGET) // prime) % prime
            points.append((top_indices, first_residual))
    return sorted(points)


def check_candidate_primes_and_nonlifts() -> int:
    checks = 0
    actual = {}
    for prime in CANDIDATE_PRIMES:
        require(is_prime_64(prime), f"candidate primality {prime}")
        points = singular_points_and_nonlifts(prime)
        count = len(points)
        expected = EXPECTED_SINGULAR_COUNTS.get(prime, 0)
        require(count == expected, f"singular point count mod {prime}: {count} != {expected}")
        if points:
            residuals = sorted(residual for _, residual in points)
            require(
                residuals == EXPECTED_FIRST_RESIDUALS[prime],
                f"first-lift residuals mod {prime}",
            )
            require(all(residual != 0 for residual in residuals), f"unexpected lift mod {prime}^2")
            actual[prime] = count
            checks += len(points) * 6
        checks += 2
    require(actual == EXPECTED_SINGULAR_COUNTS, "actual singular-prime classification")
    return checks


def check_centered_lift_law() -> int:
    centered_integer = 32768 * TARGET + 3453
    centered_factors = (5, 7, 42667, 199447, 98610539)
    product_value = 1
    for prime in centered_factors:
        product_value *= prime
    require(product_value == centered_integer, "centered branch factorization")
    require(all(valuation(centered_integer, prime) == 1 for prime in centered_factors), "squarefree centered branch")
    # For p>8, p^2 divisibility is exactly the all-centered first-lift condition.
    require(
        all(centered_integer % (prime * prime) != 0 for prime in centered_factors if prime > 8),
        "centered p^2 nonlift",
    )
    return 3 + len(centered_factors)


def main() -> None:
    critical_checks = check_critical_value_polynomials()
    resultant_checks, numerators = check_additive_resultant_branches()
    singular_checks = check_candidate_primes_and_nonlifts()
    centered_checks = check_centered_lift_law()
    displayed = ",".join(str(prime) for prime in sorted(EXPECTED_SINGULAR_COUNTS))
    rejected = ",".join(
        str(prime) for prime in CANDIDATE_PRIMES if prime not in EXPECTED_SINGULAR_COUNTS
    )
    print(f"critical_value_checks={critical_checks}")
    print(f"additive_resultant_checks={resultant_checks}")
    print(f"branch_numerators={len(numerators)}")
    print(f"candidate_primes={len(CANDIDATE_PRIMES)}")
    print(f"singular_nonlift_checks={singular_checks}")
    print(f"centered_lift_checks={centered_checks}")
    print(f"displayed_singular_primes={displayed}")
    print(f"square_or_field_rejected_candidates={rejected}")
    print("all_checks=PASS")


if __name__ == "__main__":
    main()
