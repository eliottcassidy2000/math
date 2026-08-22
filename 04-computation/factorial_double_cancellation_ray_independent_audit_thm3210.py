#!/usr/bin/env python3
"""Independent exact audit for THM-3210.

The implementation uses hand-written sparse polynomial arithmetic, exact
Sylvester determinants, Euclidean gcds, and finite-field irreducibility tests.
It imports neither SymPy nor the primary THM-3210 companion.
"""

from fractions import Fraction
from itertools import permutations


def require(condition, data):
    if not condition:
        raise RuntimeError(data)


# ---------------------------------------------------------------------------
# Sparse Q[N,D,V] polynomials.  Exponent triples are (N,D,V).

def mclean(poly):
    return {exponents: coefficient for exponents, coefficient in poly.items()
            if coefficient}


def mconst(value):
    value = Fraction(value)
    return {} if not value else {(0, 0, 0): value}


def mvar(slot):
    exponents = [0, 0, 0]
    exponents[slot] = 1
    return {tuple(exponents): Fraction(1)}


def madd(left, right):
    answer = dict(left)
    for exponents, coefficient in right.items():
        answer[exponents] = answer.get(exponents, Fraction(0)) + coefficient
    return mclean(answer)


def mneg(poly):
    return {exponents: -coefficient
            for exponents, coefficient in poly.items()}


def msub(left, right):
    return madd(left, mneg(right))


def mmul(left, right):
    answer = {}
    for powers_left, coefficient_left in left.items():
        for powers_right, coefficient_right in right.items():
            powers = tuple(a + b for a, b in zip(powers_left, powers_right))
            answer[powers] = answer.get(powers, Fraction(0)) \
                + coefficient_left * coefficient_right
    return mclean(answer)


def mpow(poly, exponent):
    answer = MONE
    for _ in range(exponent):
        answer = mmul(answer, poly)
    return answer


def mproduct(values):
    answer = MONE
    for value in values:
        answer = mmul(answer, value)
    return answer


MZERO = mconst(0)
MONE = mconst(1)
N = mvar(0)
D = mvar(1)
V = mvar(2)
DELTA = msub(MONE, mmul(mconst(4), mmul(D, V)))


def shifted(poly, offset):
    return madd(poly, mconst(offset))


def weights(index):
    index_plus_one = shifted(index, 1)
    b_value = mmul(mmul(index, index_plus_one), DELTA)
    u_value = mneg(b_value)
    c_value = msub(msub(D, index), MONE)
    alpha_value = mproduct((mconst(2), D, index_plus_one,
                            madd(mmul(mconst(2), index), MONE), V))
    beta_value = mmul(D, b_value)
    return u_value, c_value, alpha_value, beta_value


def amplitude(length, index):
    continuant = [MONE]
    if length >= 2:
        continuant.append(weights(shifted(index, 1))[2])
    for rank in range(2, length):
        _, _, alpha_value, beta_value = weights(shifted(index, rank))
        continuant.append(madd(
            mmul(alpha_value, continuant[rank - 1]),
            mmul(mmul(D, beta_value), continuant[rank - 2])))
    total = MZERO
    for exit_time in range(1, length):
        _, c_value, _, _ = weights(shifted(index, exit_time))
        tail = mproduct(weights(shifted(index, tail_index))[0]
                        for tail_index in range(exit_time + 1, length))
        total = madd(total,
                     mproduct((c_value, continuant[exit_time - 1], tail)))
    return total


# Exact general length-three formula and its rational root graph.
E3 = amplitude(3, N)
expected_e3 = mmul(
    shifted(N, 2),
    madd(
        mneg(mproduct((msub(D, shifted(N, 2)), shifted(N, 3), DELTA))),
        mproduct((mconst(2), madd(mmul(mconst(2), N), mconst(3)),
                  V, D, msub(D, shifted(N, 3))))))
require(E3 == expected_e3, "general E3 formula")


def coefficient_in_v(poly, degree):
    return {(power_n, power_d, 0): coefficient
            for (power_n, power_d, power_v), coefficient in poly.items()
            if power_v == degree}


e3_constant = coefficient_in_v(E3, 0)
e3_linear = coefficient_in_v(E3, 1)
root_numerator = mmul(shifted(N, 3), msub(D, shifted(N, 2)))
root_denominator = mproduct((mconst(2), D,
    msub(mmul(madd(mmul(mconst(4), N), mconst(9)), D),
         mmul(shifted(N, 3), madd(mmul(mconst(4), N), mconst(7))))))
require(madd(mmul(e3_constant, root_denominator),
             mmul(e3_linear, root_numerator)) == MZERO,
        "length-three rational cancellation graph")


# ---------------------------------------------------------------------------
# Sparse Q[N] arithmetic for the symbolic ray d=N+4,
# v=(N+3)/(3(N+4)(2N+5)).

def uclean(poly):
    return {degree: coefficient for degree, coefficient in poly.items()
            if coefficient}


def uconst(value):
    value = Fraction(value)
    return {} if not value else {0: value}


def uadd(left, right):
    answer = dict(left)
    for degree, coefficient in right.items():
        answer[degree] = answer.get(degree, Fraction(0)) + coefficient
    return uclean(answer)


def uneg(poly):
    return {degree: -coefficient for degree, coefficient in poly.items()}


def usub(left, right):
    return uadd(left, uneg(right))


def umul(left, right):
    answer = {}
    for degree_left, coefficient_left in left.items():
        for degree_right, coefficient_right in right.items():
            degree = degree_left + degree_right
            answer[degree] = answer.get(degree, Fraction(0)) \
                + coefficient_left * coefficient_right
    return uclean(answer)


def upow(poly, exponent):
    answer = UONE
    for _ in range(exponent):
        answer = umul(answer, poly)
    return answer


def uproduct(values):
    answer = UONE
    for value in values:
        answer = umul(answer, value)
    return answer


UZERO = uconst(0)
UONE = uconst(1)
UX = {1: Fraction(1)}


def ushift(offset):
    return uadd(UX, uconst(offset))


RAY_D = ushift(4)
RAY_V_NUMERATOR = ushift(3)
RAY_V_DENOMINATOR = uproduct(
    (uconst(3), ushift(4), uadd(umul(uconst(2), UX), uconst(5))))


def clear_ray_denominator(poly):
    max_v_degree = max((powers[2] for powers in poly), default=0)
    numerator = UZERO
    for (power_n, power_d, power_v), coefficient in poly.items():
        term = uproduct((uconst(coefficient), upow(UX, power_n),
                         upow(RAY_D, power_d),
                         upow(RAY_V_NUMERATOR, power_v),
                         upow(RAY_V_DENOMINATOR,
                              max_v_degree - power_v)))
        numerator = uadd(numerator, term)
    return numerator, upow(RAY_V_DENOMINATOR, max_v_degree)


ray_values = {}
for length in range(2, 6):
    ray_values[length] = clear_ray_denominator(amplitude(length, N))

require(ray_values[2][0] == umul(uconst(2), ray_values[2][1]), "ray E2")
require(ray_values[3][0] == UZERO, "ray E3")
require(ray_values[4][0] == UZERO, "ray E4")

e5_expected_numerator = uneg(uproduct((
    uconst(4), ushift(2), upow(ushift(3), 2), ushift(4),
    uadd(umul(uconst(2), UX), uconst(3)),
    uadd(uadd(umul(uconst(10), upow(UX, 3)),
              umul(uconst(101), upow(UX, 2))),
         uadd(umul(uconst(336), UX), uconst(366))))))
e5_expected_denominator = umul(
    uconst(27), upow(uadd(umul(uconst(2), UX), uconst(5)), 2))
require(umul(ray_values[5][0], e5_expected_denominator)
        == umul(e5_expected_numerator, ray_values[5][1]), "ray E5")


# ---------------------------------------------------------------------------
# Fixed-n resultants in Q[D], computed by exact Sylvester determinants.

def specialize_n(poly, n_value):
    answer = {}
    for (power_n, power_d, power_v), coefficient in poly.items():
        key = (power_d, power_v)
        answer[key] = answer.get(key, Fraction(0)) \
            + coefficient * Fraction(n_value) ** power_n
    return {key: coefficient for key, coefficient in answer.items()
            if coefficient}


def coefficients_in_v(poly):
    degree_v = max((power_v for _, power_v in poly), default=-1)
    answer = []
    for power_v in range(degree_v + 1):
        answer.append({power_d: coefficient
                       for (power_d, current_v), coefficient in poly.items()
                       if current_v == power_v})
    return answer


def permutation_sign(permutation):
    inversions = sum(permutation[left] > permutation[right]
                     for left in range(len(permutation))
                     for right in range(left + 1, len(permutation)))
    return -1 if inversions % 2 else 1


def determinant(matrix):
    size = len(matrix)
    answer = UZERO
    for permutation in permutations(range(size)):
        term = uconst(permutation_sign(permutation))
        for row, column in enumerate(permutation):
            term = umul(term, matrix[row][column])
        answer = uadd(answer, term)
    return answer


def resultant_v(left, right):
    left_coefficients = list(reversed(coefficients_in_v(left)))
    right_coefficients = list(reversed(coefficients_in_v(right)))
    left_degree = len(left_coefficients) - 1
    right_degree = len(right_coefficients) - 1
    size = left_degree + right_degree
    matrix = [[UZERO for _ in range(size)] for _ in range(size)]
    for row in range(right_degree):
        for offset, coefficient in enumerate(left_coefficients):
            matrix[row][row + offset] = coefficient
    for local_row in range(left_degree):
        row = right_degree + local_row
        for offset, coefficient in enumerate(right_coefficients):
            matrix[row][local_row + offset] = coefficient
    return determinant(matrix)


def udegree(poly):
    return max(poly, default=-1)


def udivmod(numerator, denominator):
    require(denominator != UZERO, "polynomial division by zero")
    quotient = UZERO
    remainder = dict(numerator)
    denominator_degree = udegree(denominator)
    denominator_lead = denominator[denominator_degree]
    while remainder and udegree(remainder) >= denominator_degree:
        degree = udegree(remainder) - denominator_degree
        coefficient = remainder[udegree(remainder)] / denominator_lead
        term = {degree: coefficient}
        quotient = uadd(quotient, term)
        remainder = usub(remainder, umul(term, denominator))
    return quotient, remainder


def umonic(poly):
    require(poly != UZERO, "monic zero")
    lead = poly[udegree(poly)]
    return {degree: coefficient / lead for degree, coefficient in poly.items()}


def ugcd(left, right):
    while right:
        _, remainder = udivmod(left, right)
        left, right = right, remainder
    return umonic(left)


quartics = {
    1: [1008, -672, 288, -96, 13],
    2: [3600, -1800, 600, -160, 17],
    3: [3300, -1320, 360, -80, 7],
}
resultant_constants = {1: 2160, 2: 8960, 3: 81000}
irreducibility_primes = {1: 17, 2: 19, 3: 13}
resultant_rows = []


# Elementary F_p[x] routines for the quartic irreducibility certificates.
def fp_trim(poly, prime):
    answer = [coefficient % prime for coefficient in poly]
    while answer and answer[-1] == 0:
        answer.pop()
    return answer


def fp_divmod(numerator, denominator, prime):
    numerator = fp_trim(numerator, prime)
    denominator = fp_trim(denominator, prime)
    require(denominator, "finite-field division by zero")
    quotient = [0] * max(1, len(numerator) - len(denominator) + 1)
    inverse_lead = pow(denominator[-1], -1, prime)
    while len(numerator) >= len(denominator):
        degree = len(numerator) - len(denominator)
        coefficient = numerator[-1] * inverse_lead % prime
        quotient[degree] = coefficient
        for offset, value in enumerate(denominator):
            numerator[degree + offset] = (
                numerator[degree + offset] - coefficient * value) % prime
        numerator = fp_trim(numerator, prime)
    return fp_trim(quotient, prime), numerator


def fp_gcd(left, right, prime):
    while right:
        _, remainder = fp_divmod(left, right, prime)
        left, right = right, remainder
    inverse_lead = pow(left[-1], -1, prime)
    return fp_trim([value * inverse_lead for value in left], prime)


def fp_mul_mod(left, right, modulus, prime):
    product_value = [0] * (len(left) + len(right) - 1)
    for left_degree, left_value in enumerate(left):
        for right_degree, right_value in enumerate(right):
            product_value[left_degree + right_degree] = (
                product_value[left_degree + right_degree]
                + left_value * right_value) % prime
    return fp_divmod(product_value, modulus, prime)[1]


def fp_pow_mod(base, exponent, modulus, prime):
    answer = [1]
    while exponent:
        if exponent & 1:
            answer = fp_mul_mod(answer, base, modulus, prime)
        base = fp_mul_mod(base, base, modulus, prime)
        exponent //= 2
    return answer


for index in (1, 2, 3):
    poly3 = specialize_n(amplitude(3, mconst(index)), index)
    poly4 = specialize_n(amplitude(4, mconst(index)), index)
    poly5 = specialize_n(amplitude(5, mconst(index)), index)
    # amplitude(..., constant index) contains no N, so specialize_n's first
    # argument is a no-op on its mathematical content.
    res34 = resultant_v(poly3, poly4)
    res35 = resultant_v(poly3, poly5)
    quartic = {degree: Fraction(coefficient)
               for degree, coefficient in enumerate(quartics[index])}
    expected_res34 = uproduct((
        uconst(resultant_constants[index]), upow(UX, 2),
        usub(UX, uconst(index + 4)), quartic))
    require(res34 == expected_res34, ("resultant E3/E4", index))
    require(ugcd(res34, res35) == upow(UX, 2),
            ("no triple cancellation", index))

    prime = irreducibility_primes[index]
    modulus = fp_trim(quartics[index], prime)
    x_poly = [0, 1]
    x_p2_minus_x = fp_trim(
        [(value - (1 if degree == 1 else 0)) % prime
         for degree, value in enumerate(
             fp_pow_mod(x_poly, prime ** 2, modulus, prime))], prime)
    require(fp_gcd(modulus, x_p2_minus_x, prime) == [1],
            ("quartic degree-1/2 factor", index, prime))
    x_p4_minus_x = fp_trim(
        [(value - (1 if degree == 1 else 0)) % prime
         for degree, value in enumerate(
             fp_pow_mod(x_poly, prime ** 4, modulus, prime))], prime)
    require(x_p4_minus_x == [],
            ("quartic Frobenius closure", index, prime))
    resultant_rows.append((index, resultant_constants[index], prime))


def evaluate(poly, n_value, d_value, v_value):
    return sum(coefficient * Fraction(n_value) ** power_n
               * Fraction(d_value) ** power_d
               * Fraction(v_value) ** power_v
               for (power_n, power_d, power_v), coefficient in poly.items())


require(evaluate(E3, 1, 5, Fraction(4, 105)) == 0,
        "THM-3186 hostile point")
require(evaluate(amplitude(4, mconst(1)), 1, 6, Fraction(4, 105)) != 0,
        "off-ray positive control")


print("THM-3210 INDEPENDENT DOUBLE-CANCELLATION AUDIT")
print("implementation=custom Q[N,D,V] + Sylvester + F_p arithmetic")
print("general_length3_formula_and_rational_locus=PASS")
print("ray_identities=(E2=2,E3=0,E4=0,E5=closed_nonzero)")
print("visibility_profile=(visible,invisible,invisible,visible)")
print("resultant_irreducibility_rows=" + repr(resultant_rows))
print("gcd_resultants_E34_E35=d^2_for_n=1,2,3")
print("quartics_irreducible_over_Q_by_mod_primes=(17,19,13)")
print("hostile_and_off_ray_positive_control=PASS")
print("scope=Q(n,d)[v]_degree_and_scalar_visibility_not_PRS_depth")
print("INDEPENDENT EXACT AUDIT PASSED")
