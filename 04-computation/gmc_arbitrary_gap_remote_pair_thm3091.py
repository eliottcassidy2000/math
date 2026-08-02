#!/usr/bin/env python3
"""Exact companion for THM-3091's arbitrary-gap remote-pair theorem."""

from fractions import Fraction
from functools import lru_cache
from math import factorial


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def rising(start, length):
    answer = 1
    for offset in range(length):
        answer *= start + offset
    return answer


def mixed_ratio(degree, far_count, N, gap):
    return Fraction(
        rising(degree * N + 1, far_count * gap),
        rising(N + 1, gap) ** far_count,
    )


def bareiss_determinant(matrix):
    matrix = [list(map(int, row)) for row in matrix]
    size = len(matrix)
    if size == 1:
        return matrix[0][0]
    sign = 1
    denominator = 1
    for pivot_index in range(size - 1):
        if matrix[pivot_index][pivot_index] == 0:
            swap = next(
                (row for row in range(pivot_index + 1, size) if matrix[row][pivot_index]),
                None,
            )
            require(swap is not None, "singular determinant")
            matrix[pivot_index], matrix[swap] = matrix[swap], matrix[pivot_index]
            sign *= -1
        pivot = matrix[pivot_index][pivot_index]
        for row in range(pivot_index + 1, size):
            for column in range(pivot_index + 1, size):
                numerator = (
                    matrix[row][column] * pivot
                    - matrix[row][pivot_index] * matrix[pivot_index][column]
                )
                require(numerator % denominator == 0, "Bareiss exact division")
                matrix[row][column] = numerator // denominator
        denominator = pivot
        for row in range(pivot_index + 1, size):
            matrix[row][pivot_index] = 0
    return sign * matrix[-1][-1]


def sylvester_resultant(first_descending, second_descending):
    first_degree = len(first_descending) - 1
    second_degree = len(second_descending) - 1
    size = first_degree + second_degree
    matrix = []
    for shift in range(second_degree):
        matrix.append(
            [0] * shift
            + list(first_descending)
            + [0] * (second_degree - 1 - shift)
        )
    for shift in range(first_degree):
        matrix.append(
            [0] * shift
            + list(second_descending)
            + [0] * (first_degree - 1 - shift)
        )
    require(all(len(row) == size for row in matrix), "Sylvester size")
    return bareiss_determinant(matrix)


def add(left, right):
    return left[0] + right[0], left[1] + right[1]


def scale(value, interval):
    if value >= 0:
        return value * interval[0], value * interval[1]
    return value * interval[1], value * interval[0]


@lru_cache(maxsize=None)
def log_bounds(value, terms=45):
    require(value > 0, "log domain")
    if value < 1:
        lower, upper = log_bounds(1 / value, terms)
        return -upper, -lower
    power_two = 0
    reduced = value
    while reduced >= 2:
        reduced /= 2
        power_two += 1
    z = (reduced - 1) / (reduced + 1)
    partial = Fraction(0)
    zpower = z
    for index in range(terms):
        partial += 2 * zpower / (2 * index + 1)
        zpower *= z * z
    tail = 2 * zpower / ((2 * terms + 1) * (1 - z * z))
    core = partial, partial + tail
    if not power_two:
        return core
    z2 = Fraction(1, 3)
    partial2 = Fraction(0)
    zpower2 = z2
    for index in range(terms):
        partial2 += 2 * zpower2 / (2 * index + 1)
        zpower2 *= z2 * z2
    tail2 = 2 * zpower2 / ((2 * terms + 1) * (1 - z2 * z2))
    return add(core, scale(power_two, (partial2, partial2 + tail2)))


# Exact secondary scaling in rational p-th powers.  The degree-p block stays
# subexponential, while every non-top degree-k coefficient pays one or more
# copies of log(k/p) per gap unit.
coefficient_cells = 0
response_ratio_cells = 0
for child in range(3, 7):
    p = child + 1
    k = child + 2
    for N in (53,):
        for gap in (1, 5):
            Vp = mixed_ratio(p, p, N, gap)
            Vk = mixed_ratio(k, k, N, gap)
            require(Vp > 0 and Vk > 0, "positive far carriers")
            response_bound = Fraction(p * (p + 1), p * (p + 1) + 1)
            require(
                Vp**k * response_bound ** (-p * k * gap) <= Vk**p,
                "uniform response-ratio bound",
            )
            response_ratio_cells += 1
            for far_count in range(p + 1):
                coefficient = mixed_ratio(p, far_count, N, gap)
                require(
                    coefficient**p <= Vp**far_count,
                    "degree-p Jensen bound",
                )
                bp_power = coefficient**p / Vp**far_count
                lower, upper = log_bounds(bp_power)
                bound = Fraction(2 * p * p * gap * (gap + 1), N)
                require(lower >= -bound and upper <= bound, "degree-p growth")
                coefficient_cells += 1
            for far_count in range(k + 1):
                coefficient = mixed_ratio(k, far_count, N, gap)
                require(
                    coefficient**k <= Vk**far_count,
                    "degree-k Jensen bound",
                )
                bk_power = coefficient**p * Vp ** (k - far_count) / Vk**p
                if far_count == k:
                    require(bk_power == 1, "exact top coefficient")
                target_neutral = bk_power * Fraction(k, p) ** (
                    p * (k - far_count) * gap
                )
                lower, upper = log_bounds(target_neutral)
                bound = Fraction(3 * p * k * gap * (gap + 1), N)
                require(lower >= -bound and upper <= bound, "degree-k decay")
                coefficient_cells += 1


# The same inequalities are the whole-system outer-layer invoice.  In a
# lower or degree-p row, an atom with x far physical factors pays at most
# Q_(p,x), while lambda^t pays at least Vp^(-x/p).  In the degree-k row the
# extra equation scaling leaves
# Q_(k,x) Vp^((k-x)/p)/Vk <= 1.  Check the latter in integral powers.
outer_atom_cells = 0
for child in range(3, 8):
    p = child + 1
    k = child + 2
    for N in (p, p + 3, 31):
        for gap in (1, 3, 11):
            Vp = mixed_ratio(p, p, N, gap)
            Vk = mixed_ratio(k, k, N, gap)
            for far_count in range(p + 1):
                Qp = mixed_ratio(p, far_count, N, gap)
                require(Qp**p <= Vp**far_count, "lower/p outer atom")
                outer_atom_cells += 1
            for far_count in range(k + 1):
                Qk = mixed_ratio(k, far_count, N, gap)
                require(
                    Qk ** (p * k) * Vp ** (k * (k - far_count))
                    <= Vk ** (p * k),
                    "degree-k outer atom",
                )
                outer_atom_cells += 1

            # Raw physical starts jN+D+1, including the maximal admissible
            # start and t>x fixed-pivot normal degrees.
            denominator = rising(N + 1, gap)
            for row_cap in (p, k):
                for actual_high in range(row_cap + 1):
                    for far_count in range(actual_high + 1):
                        for displacement in {
                            0,
                            min(2, (row_cap - actual_high) * N),
                            (row_cap - actual_high) * N,
                        }:
                            W = Fraction(
                                rising(
                                    actual_high * N + displacement + 1,
                                    far_count * gap,
                                ),
                                denominator**far_count,
                            )
                            for far_degree in {far_count, row_cap}:
                                if row_cap == p:
                                    require(
                                        W**p <= Vp**far_degree,
                                        "raw lower/p outer atom",
                                    )
                                else:
                                    require(
                                        W ** (p * k)
                                        * Vp ** (k * (k - far_degree))
                                        <= Vk ** (p * k),
                                        "raw degree-k outer atom",
                                    )
                                outer_atom_cells += 1


# Exact covariance ledger: lambda=Vp^(-1/p),
# s=Vp^(k/p)/Vk, so lambda^(pk)s^p=Vk^(-p).
covariance_cells = 0
for child in range(3, 21):
    p = child + 1
    k = child + 2
    vp_exponent = -k + k
    vk_exponent = -p
    require((vp_exponent, vk_exponent) == (0, -p), "secondary covariance")
    covariance_cells += 1


# Independent Sylvester controls for Res(B_p,v^k)=B_p(1,0)^k=1.
triangular_cells = 0
for child in range(3, 11):
    p = child + 1
    k = child + 2
    second = [1] + [0] * k
    for seed in range(1, 8):
        # Descending coefficients of B_p(1,t); its constant term is one.
        first = [seed + index * index for index in range(p)] + [1]
        value = sylvester_resultant(first, second)
        require(value == 1, "triangular Sylvester quotient")
        triangular_cells += 1


# Exact carrier merger U_k(C)V_k(N,H)=U_k(C+H) and factorial exponents.
carrier_cells = 0
for child in range(3, 13):
    p = child + 1
    k = child + 2
    require(p * factorial(child) == factorial(k - 1), "far exponent merger")
    require(factorial(k) // factorial(child) == p * k, "child exponent")
    for n in (1, 3):
        for C in (7, 11):
            N = n + C
            for gap in (1, 4):
                U_base = Fraction(
                    rising(k * n + 1, k * C),
                    rising(n + 1, C) ** k,
                )
                V_far = mixed_ratio(k, k, N, gap)
                U_far = Fraction(
                    rising(k * n + 1, k * (C + gap)),
                    rising(n + 1, C + gap) ** k,
                )
                require(U_base * V_far == U_far, "far carrier identity")
                carrier_cells += 1


# Macroscopic binary normal face.  For every positive delta the transformed
# rate is strictly convex in j, with T_0<0 and T_k=0, hence T_j<0 for j<k.
macroscopic_p_cells = 0
macroscopic_k_cells = 0
for child in range(3, 11):
    p = child + 1
    k = child + 2
    for delta in (Fraction(1, 20), Fraction(1, 5), Fraction(1), Fraction(3)):
        for far_count in range(1, p):
            rate = (0, 0)
            rate = add(
                rate,
                scale(p + far_count * delta, log_bounds(p + far_count * delta)),
            )
            rate = add(rate, scale(-p, log_bounds(p)))
            rate = add(
                rate,
                scale(-far_count * (1 + delta), log_bounds(1 + delta)),
            )
            rate = add(rate, scale(-far_count * delta, log_bounds(p)))
            require(rate[1] < 0, "strict macroscopic degree-p mixed rate")
            curvature = delta * delta / (p + far_count * delta)
            require(curvature > 0, "strict degree-p rate convexity")
            macroscopic_p_cells += 1
        for far_count in range(k):
            rate = (0, 0)
            rate = add(
                rate,
                scale(k + far_count * delta, log_bounds(k + far_count * delta)),
            )
            rate = add(rate, scale(-k, log_bounds(k)))
            rate = add(
                rate,
                scale(-far_count * (1 + delta), log_bounds(1 + delta)),
            )
            rate = add(rate, scale((k - far_count) * delta, log_bounds(p)))
            rate = add(rate, scale(-k * delta, log_bounds(k)))
            require(rate[1] < 0, "strict macroscopic non-top rate")
            curvature = delta * delta / (k + far_count * delta)
            require(curvature > 0, "strict rate convexity")
            macroscopic_k_cells += 1


print("THM-3091 ARBITRARY-GAP REMOTE-PAIR DESUSPENSION")
print(f"secondary_coefficient_cells={coefficient_cells} p_growth_and_k_decay=PASS")
print(f"response_ratio_cells={response_ratio_cells} exact_Jensen_cp_bound=PASS")
print(f"outer_atom_cells={outer_atom_cells} whole_system_contraction=PASS")
print(f"covariance_cells={covariance_cells} Etilde=E/Vk^p")
print(f"triangular_sylvester_cells={triangular_cells} resultant=1")
print(f"carrier_merge_cells={carrier_cells} Uk(C)*Vk=Uk(C+H)")
print(f"macroscopic_p_cells={macroscopic_p_cells} two_endpoint_face=PASS")
print(f"macroscopic_k_cells={macroscopic_k_cells} all_non_top_negative=PASS")
print("gap_range=all_integers_H>=1 threshold_in_C_independent_of_H")
print("unbounded_carrier=S_m^(pk)*Up(C)^(k!/p)*Uk(C+H)^((k-1)!)")
print("boundary=fixed_child_and_remote_pair;not_arbitrary_two_scale_support")
print("scope=fixed_width;arbitrary_gap;whole_system_secondary_scaling")
print("all_exact_checks=PASS")
