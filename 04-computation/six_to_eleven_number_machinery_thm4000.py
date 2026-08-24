#!/usr/bin/env python3
"""Exact companion for THM-4000.

Universe and controls
---------------------
* integer polynomials on consecutive sample nodes;
* bases 2..30 and moduli 2..512 for centered/tripotent checks;
* Farey--Pell matrices through depth 15;
* power-of-two cyclotomic folds through order 32;
* the first four theta coefficients of A_1,A_2,A_3 stabilized by E8;
* odd-cycle packing partitions at covered sizes 6..11.

All checks use Python integers.  Positive controls are the base-7 and base-10
packets.  Hostiles are the factor-two loss from residues at base 7, the extra
tripotents modulo 8 and 10, the first OCF scalarization kernel at covered size
9, and the false Z^8=E8 theta identification corrected in MISTAKE-471.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations, product
from math import comb, factorial, gcd, isqrt, lcm
from pathlib import Path
import hashlib
import random


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def poly_eval(coeffs: tuple[int, ...] | list[int], x: int) -> int:
    value = 0
    for coefficient in reversed(coeffs):
        value = value * x + coefficient
    return value


def integer_binomial(n: int, k: int) -> int:
    """The polynomial binomial n(n-1)...(n-k+1)/k! for every integer n."""
    numerator = 1
    for offset in range(k):
        numerator *= n - offset
    return numerator // factorial(k)


def forward_heads(values: list[int]) -> list[int]:
    row = values[:]
    heads = []
    while row:
        heads.append(row[0])
        row = [row[i + 1] - row[i] for i in range(len(row) - 1)]
    return heads


def consecutive_remainder_value(
    coeffs: tuple[int, ...] | list[int], a: int, m: int, at: int
) -> int:
    values = [poly_eval(coeffs, a + j) for j in range(m + 1)]
    heads = forward_heads(values)
    return sum(integer_binomial(at - a, j) * heads[j] for j in range(m + 1))


def observer_weights(a: int, m: int, at: int) -> list[int]:
    weights = []
    for k in range(m + 1):
        values = [int(j == k) for j in range(m + 1)]
        heads = forward_heads(values)
        weights.append(
            sum(integer_binomial(at - a, j) * heads[j] for j in range(m + 1))
        )
    return weights


def factorint(n: int) -> dict[int, int]:
    n = abs(n)
    factors: dict[int, int] = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            factors[p] = factors.get(p, 0) + 1
            n //= p
        p = 3 if p == 2 else p + 2
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors


def multiplicative_order(a: int, modulus: int) -> int:
    gate(gcd(a, modulus) == 1, "multiplicative order needs a unit")
    x = 1
    for order in range(1, modulus + 1):
        x = (x * a) % modulus
        if x == 1:
            return order
    raise RuntimeError("order search failed")


def tripotents(n: int) -> list[int]:
    return [x for x in range(n) if (x * x * x - x) % n == 0]


def predicted_tripotent_count(n: int) -> int:
    factors = factorint(n)
    two_exp = factors.pop(2, 0)
    two_factor = 1 if two_exp == 0 else (2 if two_exp == 1 else (3 if two_exp == 2 else 5))
    return two_factor * 3 ** len(factors)


def is_odd_prime_power(n: int) -> bool:
    factors = factorint(n)
    return len(factors) == 1 and 2 not in factors


def matrix_mul(a: list[list[int]], b: list[list[int]]) -> list[list[int]]:
    return [
        [sum(a[i][k] * b[k][j] for k in range(len(b))) for j in range(len(b[0]))]
        for i in range(len(a))
    ]


def matrix_pow(a: list[list[int]], n: int) -> list[list[int]]:
    result = [[1, 0], [0, 1]]
    base = a
    while n:
        if n & 1:
            result = matrix_mul(result, base)
        base = matrix_mul(base, base)
        n //= 2
    return result


def chebyshev_pair(b: int, n: int) -> tuple[int, int]:
    """Return T_n(b), U_(n-1)(b), with U_-1=0."""
    if n == 0:
        return 1, 0
    t_prev, t = 1, b
    u_prev, u = 0, 1
    for _ in range(1, n):
        t_prev, t = t, 2 * b * t - t_prev
        u_prev, u = u, 2 * b * u - u_prev
    return t, u


def determinant_bareiss(matrix: list[list[int]]) -> int:
    a = [row[:] for row in matrix]
    n = len(a)
    if n == 0:
        return 1
    sign = 1
    previous = 1
    for k in range(n - 1):
        if a[k][k] == 0:
            pivot = next((i for i in range(k + 1, n) if a[i][k] != 0), None)
            if pivot is None:
                return 0
            a[k], a[pivot] = a[pivot], a[k]
            sign *= -1
        pivot_value = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                a[i][j] = (a[i][j] * pivot_value - a[i][k] * a[k][j]) // previous
        previous = pivot_value
    return sign * a[-1][-1]


def smith_diagonal_from_minors(matrix: list[list[int]]) -> list[int]:
    """Recover the Smith diagonal of a full-rank square integer matrix."""
    size = len(matrix)
    previous_divisor = 1
    diagonal = []
    for minor_size in range(1, size + 1):
        determinantal_divisor = 0
        for rows in combinations(range(size), minor_size):
            for columns in combinations(range(size), minor_size):
                minor = [[matrix[i][j] for j in columns] for i in rows]
                determinantal_divisor = gcd(
                    determinantal_divisor, abs(determinant_bareiss(minor))
                )
        if determinantal_divisor == 0 or determinantal_divisor % previous_divisor:
            raise RuntimeError("invalid determinantal-divisor chain")
        diagonal.append(determinantal_divisor // previous_divisor)
        previous_divisor = determinantal_divisor
    return diagonal


def sigma3(n: int) -> int:
    return sum(d**3 for d in range(1, n + 1) if n % d == 0)


def a_cartan(rank: int) -> list[list[int]]:
    return [
        [2 if i == j else (-1 if abs(i - j) == 1 else 0) for j in range(rank)]
        for i in range(rank)
    ]


def theta_a(rank: int, max_exponent: int) -> list[int]:
    counts = [0] * (max_exponent + 1)
    bound = isqrt(2 * max_exponent)
    for prefix in product(range(-bound, bound + 1), repeat=rank):
        coordinates = prefix + (-sum(prefix),)
        norm = sum(x * x for x in coordinates)
        gate(norm % 2 == 0, "A_r must be even")
        exponent = norm // 2
        if exponent <= max_exponent:
            counts[exponent] += 1
    return counts


def odd_partitions(total: int, minimum: int = 3) -> list[tuple[int, ...]]:
    answers: list[tuple[int, ...]] = []

    def rec(remaining: int, least: int, parts: list[int]) -> None:
        if remaining == 0:
            answers.append(tuple(parts))
            return
        for part in range(least, remaining + 1, 2):
            rec(remaining - part, part, parts + [part])

    rec(total, minimum, [])
    return answers


print("THM-4000 exact companion: six-to-eleven number machinery")
print()

# ---------------------------------------------------------------------------
# 1. Consecutive-node sampler and its factorial Smith lattice
# ---------------------------------------------------------------------------

sampler_cases = 0
for coeffs in product((-1, 0, 1), repeat=7):
    for a in (-2, -1, 0):
        for m in (2, 5):
            for at in (a - 5, a - 1, a + m + 1, a + m + 3, a + m + 7):
                modulus = 1
                for j in range(m + 1):
                    modulus *= at - (a + j)
                remainder = consecutive_remainder_value(coeffs, a, m, at)
                gate((poly_eval(coeffs, at) - remainder) % modulus == 0,
                     "consecutive sampler congruence")
                sampler_cases += 1

for m in range(0, 8):
    falling_eval = [
        [factorial(k) * comb(j, k) if j >= k else 0 for k in range(m + 1)]
        for j in range(m + 1)
    ]
    gate(determinant_bareiss(falling_eval) ==
         __import__("functools").reduce(lambda x, y: x * y,
                                         (factorial(k) for k in range(m + 1)), 1),
         "factorial evaluation determinant")

six_falling_eval = [
    [factorial(k) * comb(j, k) if j >= k else 0 for k in range(6)]
    for j in range(6)
]
smith_six = smith_diagonal_from_minors(six_falling_eval)
gate(smith_six == [1, 1, 2, 6, 24, 120], "six-node Smith diagonal")
print("consecutive_sampler_cases=", sampler_cases, sep="")
print("evaluation_snf_six_nodes=", smith_six, sep="")
print("optimality_control=P=product(X-a-j) attains the full modulus")

# The requested six-number packet.
nodes_6_11 = list(range(-1, 5))
weights_6_11 = observer_weights(-1, 5, 10)
gate(weights_6_11 == [-252, 1386, -3080, 3465, -1980, 462],
     "six-observer weights")
full_modulus = 1
for modulus in range(6, 12):
    full_modulus *= modulus
naive_modulus = lcm(*range(6, 12))
gate(full_modulus == 332640 and naive_modulus == 27720 and
     full_modulus // naive_modulus == 12, "six-observer modulus ledger")
for degree in range(7):
    compiled_moment = sum(
        weight * node**degree for weight, node in zip(weights_6_11, nodes_6_11)
    )
    expected_error = 0 if degree <= 5 else full_modulus
    gate(10**degree - compiled_moment == expected_error,
         "six-observer sharp monomial control")

rng = random.Random(4000)
digit_cases = 0
for _ in range(10000):
    digits = tuple(rng.randrange(10) for _ in range(rng.randrange(1, 25)))
    samples = [poly_eval(digits, a) for a in nodes_6_11]
    compiled = sum(w * v for w, v in zip(weights_6_11, samples))
    gate((poly_eval(digits, 10) - compiled) % full_modulus == 0,
         "decimal six-observer compiler")
    digit_cases += 1

print("six_observer_nodes=", nodes_6_11, sep="")
print("six_observer_moduli=", list(reversed(range(6, 12))), sep="")
print("six_observer_weights=", weights_6_11, sep="")
print("six_observer_modulus=332640; residue_only_lcm=27720; lift_factor=12")
print("decimal_digit_controls=", digit_cases, sep="")

# ---------------------------------------------------------------------------
# 2. The split cubic, centered bases, and tripotents
# ---------------------------------------------------------------------------

split_eval_matrix = [[1, -1, 1], [1, 0, 0], [1, 1, 1]]
gate(abs(determinant_bareiss(split_eval_matrix)) == 2,
     "split-cubic evaluation index")
gate(smith_diagonal_from_minors(split_eval_matrix) == [1, 1, 2],
     "split-cubic Smith diagonal")
print("split_cubic_eval_snf=[1,1,2]; image_condition=P(-1)=P(1) mod 2")

for coeffs in product(range(-2, 3), repeat=6):
    p_minus = poly_eval(coeffs, -1)
    p_zero = poly_eval(coeffs, 0)
    p_plus = poly_eval(coeffs, 1)
    gate((p_plus - p_minus) % 2 == 0, "odd divided difference integral")
    gate((p_plus + p_minus - 2 * p_zero) % 2 == 0,
         "even divided difference integral")
    for b in (7, 10):
        reduced = (
            b * (b + 1) // 2 * p_plus
            + (1 - b * b) * p_zero
            + b * (b - 1) // 2 * p_minus
        )
        gate((poly_eval(coeffs, b) - reduced) % (b**3 - b) == 0,
             "split-cubic base compiler")

for b in range(2, 31):
    gate(gcd(b - 1, b) == gcd(b, b + 1) == 1, "neighbor gcd")
    gate(gcd(b - 1, b + 1) == gcd(b - 1, 2), "endpoint resultant gcd")
    if b % 2 == 0:
        e_plus = b * b * (b + 1) // 2
        e_zero = 1 - b * b
        e_minus = b * (b - 1) * (b + 2) // 2
        modulus = b**3 - b
        gate(e_plus % (b - 1) == 1 % (b - 1) and e_plus % b == 0 and
             e_plus % (b + 1) == 0, "CRT plus idempotent")
        gate(e_zero % (b - 1) == 0 and e_zero % b == 1 and
             e_zero % (b + 1) == 0, "CRT zero idempotent")
        gate(e_minus % (b - 1) == 0 and e_minus % b == 0 and
             e_minus % (b + 1) == 1, "CRT minus idempotent")
        gate((e_plus - e_minus - b) % modulus == 0, "CRT base decomposition")

print("base7_compiler=P(7)=28P(1)-48P(0)+21P(-1) mod 336")
print("base7_residue_only_modulus=lcm(6,7,8)=168; missing_sidecar_factor=2")
print("base10_compiler=P(10)=55P(1)-99P(0)+45P(-1) mod 990")
print("base10_crt=P(10)=550P(1)-99P(0)+540P(-1) mod 990")

for n in range(2, 513):
    roots = tripotents(n)
    gate(len(roots) == predicted_tripotent_count(n), "tripotent CRT count")
    for x in roots:
        e = x * x % n
        gate(e * e % n == e and x * e % n == x,
             "tripotent has idempotent support and is a corner involution")
    if n >= 3:
        centered_exhausts = set(roots) == {0, 1, n - 1}
        gate(centered_exhausts == (n == 4 or is_odd_prime_power(n)),
             "centered tripotent exhaustiveness classification")

root_table = {n: tripotents(n) for n in range(6, 12)}
gate(root_table == {
    6: [0, 1, 2, 3, 4, 5],
    7: [0, 1, 6],
    8: [0, 1, 3, 5, 7],
    9: [0, 1, 8],
    10: [0, 1, 4, 5, 6, 9],
    11: [0, 1, 10],
}, "six-to-eleven tripotent table")

fixed_divisor = 0
for x in range(-30, 31):
    fixed_divisor = gcd(fixed_divisor, x**3 - x)
gate(fixed_divisor == 6 and 2**3 - 2 == 6, "fixed divisor of X^3-X")
print("tripotent_roots_6_to_11=", root_table, sep="")
print("fixed_divisor(X^3-X)=6; universal_residue_moduli_are_exactly_divisors_of_6")
print("centered_roots_exhaust_for_n_ge_3_iff=odd_prime_power_or_4")

# ---------------------------------------------------------------------------
# 3. Farey cusp ray and the Pell/Chebyshev depth tower
# ---------------------------------------------------------------------------

cusp_step = [[0, 1], [-1, 2]]
farey_checks = 0
pell_checks = 0
for b in range(2, 31):
    a_b = [[b, b - 1], [b + 1, b]]
    a_next = [[b + 1, b], [b + 2, b + 1]]
    gate(determinant_bareiss(a_b) == 1, "Farey matrix determinant")
    gate(matrix_mul(cusp_step, a_b) == a_next, "parabolic cusp transition")
    gate(Fraction(b, b + 1) - Fraction(b - 1, b) == Fraction(1, b * (b + 1)),
         "Farey interval width")
    farey_checks += 1
    for n in range(1, 16):
        t_n, u_n = chebyshev_pair(b, n)
        predicted = [[t_n, (b - 1) * u_n], [(b + 1) * u_n, t_n]]
        gate(matrix_pow(a_b, n) == predicted, "Chebyshev matrix power")
        gate(t_n * t_n - (b * b - 1) * u_n * u_n == 1,
             "Pell norm identity")
        lower = Fraction((b - 1) * u_n, t_n)
        upper = Fraction(t_n, (b + 1) * u_n)
        gate(upper - lower == Fraction(1, (b + 1) * t_n * u_n),
             "Pell Farey bracket width")
        gate(lower * lower < Fraction(b - 1, b + 1) < upper * upper,
             "Pell Farey bracket")
        if n % 2 == 1:
            m = (n - 1) // 2
            _, u_m = chebyshev_pair(b, m + 1)
            u_m_minus = 0 if m == 0 else chebyshev_pair(b, m)[1]
            gate(t_n - 1 == (b - 1) * (u_m + u_m_minus) ** 2,
                 "odd Chebyshev minus square")
            gate(t_n + 1 == (b + 1) * (u_m - u_m_minus) ** 2,
                 "odd Chebyshev plus square")
        pell_checks += 1

def odd_pell_centers(b: int, count: int) -> list[int]:
    return [chebyshev_pair(b, 2 * m + 1)[0] for m in range(count)]

print("farey_edges: 6/7<7/8 and 9/10<10/11 lie on one cusp-1 parabolic orbit")
print("farey_base_checks=", farey_checks, "; pell_depth_checks=", pell_checks, sep="")
print("base7_odd_pell_centers=", odd_pell_centers(7, 3), sep="")
print("base10_odd_pell_centers=", odd_pell_centers(10, 3), sep="")
print("base10_identity=X-1=9*A^2 and X+1=11*B^2 at every odd depth")

# ---------------------------------------------------------------------------
# 4. The genuine arithmetic period-eight observer
# ---------------------------------------------------------------------------

cyclotomic_cases = 0
for b in range(2, 31):
    for k in range(1, 6):
        half_order = 2 ** (k - 1)
        modulus = b**half_order + 1
        gate(pow(b, half_order, modulus) == modulus - 1, "cyclotomic half-turn")
        gate(pow(b, 2 * half_order, modulus) == 1, "cyclotomic full turn")
        coeffs = tuple(rng.randrange(-5, 6) for _ in range(4 * half_order + 3))
        folded = [0] * half_order
        for i, coefficient in enumerate(coeffs):
            block, residue = divmod(i, half_order)
            folded[residue] += coefficient if block % 2 == 0 else -coefficient
        gate((poly_eval(coeffs, b) - poly_eval(folded, b)) % modulus == 0,
             "power-of-two cyclotomic fold")
        cyclotomic_cases += 1

phi8_data = {}
for b in (7, 10):
    modulus = b**4 + 1
    factors = factorint(modulus)
    orders = {q: multiplicative_order(b % q, q) for q in factors if q % 2 == 1}
    gate(all(order == 8 and q % 8 == 1 for q, order in orders.items()),
         "Phi8 prime splitting")
    phi8_data[b] = (modulus, factors, orders)

print("power_two_cyclotomic_checks=", cyclotomic_cases, sep="")
print("phi8_data=", phi8_data, sep="")
print("phi8_fold=P(b) mod (b^4+1) is the alternating sum of four-digit blocks")

# ---------------------------------------------------------------------------
# 5. Arithmetic rank-eight stabilization, with the Z8 hostile
# ---------------------------------------------------------------------------

max_q = 4
e8_theta = [1] + [240 * sigma3(n) for n in range(1, max_q + 1)]
stable_rows = []
for rank in (1, 2, 3):
    theta = theta_a(rank, max_q)
    stabilized = [
        sum(theta[j] * e8_theta[n - j] for j in range(n + 1))
        for n in range(max_q + 1)
    ]
    discriminant = determinant_bareiss(a_cartan(rank))
    gate(discriminant == rank + 1, "A_r discriminant")
    gate(stabilized[1] == 240 + rank * (rank + 1),
         "E8-stabilized root shell")
    stable_rows.append((rank, rank + 8, discriminant, theta[1], stabilized[1]))

z8_norm2 = 16 * sum((-1) ** (2 + d) * d**3 for d in (1, 2))
z8_norm4 = 16 * sum((-1) ** (4 + d) * d**3 for d in (1, 2, 4))
gate((z8_norm2, e8_theta[1], z8_norm4, e8_theta[2]) ==
     (112, 240, 1136, 2160),
     "MISTAKE-471 theta hostile")
print("E8_stabilization_rows=(rank,rank+8,disc,roots_before,roots_after)")
print(stable_rows)
print("theta_hostile_common_shells: norm2 Z8/E8=112/240; norm4=1136/2160")

# ---------------------------------------------------------------------------
# 6. Unrelated repo connection: the OCF fugacity scalarization kernel
# ---------------------------------------------------------------------------

packing_rows = {}
for covered in range(6, 12):
    partitions = odd_partitions(covered)
    lengths = sorted({len(parts) for parts in partitions})
    kernels = []
    if lengths:
        least = lengths[0]
        for length in lengths[1:]:
            vector = {least: 2 ** (length - least), length: -1}
            gate(sum(coefficient * 2**ell for ell, coefficient in vector.items()) == 0,
                 "OCF scalarization kernel")
            kernels.append(vector)
    packing_rows[covered] = (partitions, lengths, kernels)

gate(packing_rows[9][1] == [1, 3] and packing_rows[11][1] == [1, 3] and
     packing_rows[10][1] == [2], "OCF 9-10-11 boundary")
print("ocf_packing_rows_6_to_11=", packing_rows, sep="")
print("ocf_kernel_at_9_and_11=4*e_length1-e_length3; no length kernel at 10")

# ---------------------------------------------------------------------------
# 7. Unrelated repo connection: Moon--Busch as a min-plus semigroup machine
# ---------------------------------------------------------------------------

strong_floor = {}
strong_floor_optimizers = {}
for order in range(6, 12):
    weight = order - 1
    feasible = [
        (3**a * 5**b, a, b)
        for b in range(weight // 3 + 1)
        for a in range(weight // 2 + 1)
        if 2 * a + 3 * b == weight
    ]
    optimum, a_opt, b_opt = min(feasible)
    largest_parity_b = max(
        b for b in range(weight // 3 + 1) if b % 2 == weight % 2
    )
    gate(b_opt == largest_parity_b and 2 * a_opt + 3 * b_opt == weight,
         "Moon--Busch min-plus optimizer")
    strong_floor[order] = optimum
    strong_floor_optimizers[order] = (a_opt, b_opt)

gate(strong_floor == {6: 15, 7: 25, 8: 45, 9: 75, 10: 125, 11: 225},
     "Moon--Busch six-to-eleven floor values")
print("moon_busch_minplus_floors_6_to_11=", strong_floor, sep="")
print("moon_busch_optimizers_(a,b)=", strong_floor_optimizers, sep="")

source_sha = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
print()
print("checks=", CHECKS, sep="")
print("source_sha256=", source_sha, sep="")
print("status=PASS (integer exact; theorem proofs and scope are in THM-4000)")
