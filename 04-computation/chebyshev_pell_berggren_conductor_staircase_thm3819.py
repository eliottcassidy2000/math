#!/usr/bin/env python3
"""Assertion-free exact companion for THM-3819.

The global statements in THM-3819 follow from the displayed Chebyshev,
Pell, Berggren, and module algebra.  This deterministic companion audits the
indexing, parity, operation conventions, consequence objects, finite-clock
cocycle, scalar observer, and sharp hostile boundaries.  ``gate`` is used
instead of ``assert`` so ordinary and optimized Python execute identical
checks.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from math import gcd


CHEBYSHEV_DEPTH = 18
PELL_DEPTH = 240
LRC_DEPTH = 60
MODULUS = 13

CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def trim(poly: list[int]) -> list[int]:
    result = list(poly)
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def poly_add(left: list[int], right: list[int]) -> list[int]:
    size = max(len(left), len(right))
    result = [0] * size
    for index in range(size):
        if index < len(left):
            result[index] += left[index]
        if index < len(right):
            result[index] += right[index]
    return trim(result)


def poly_neg(poly: list[int]) -> list[int]:
    return [-coefficient for coefficient in poly]


def poly_sub(left: list[int], right: list[int]) -> list[int]:
    return poly_add(left, poly_neg(right))


def poly_scale(poly: list[int], scalar: int) -> list[int]:
    return trim([scalar * coefficient for coefficient in poly])


def poly_mul(left: list[int], right: list[int]) -> list[int]:
    result = [0] * (len(left) + len(right) - 1)
    for left_index, left_coefficient in enumerate(left):
        for right_index, right_coefficient in enumerate(right):
            result[left_index + right_index] += (
                left_coefficient * right_coefficient
            )
    return trim(result)


def poly_eval(poly: list[int], value: int) -> int:
    result = 0
    for coefficient in reversed(poly):
        result = result * value + coefficient
    return result


def chebyshev_tables(depth: int) -> tuple[list[list[int]], list[list[int]]]:
    """Return C_0..C_depth and U_0..U_depth at tau=2z-1."""
    tau = [-1, 2]
    first = [[1], tau]
    second = [[1], poly_scale(tau, 2)]
    for _index in range(1, depth):
        first.append(
            poly_sub(poly_scale(poly_mul(tau, first[-1]), 2), first[-2])
        )
        second.append(
            poly_sub(poly_scale(poly_mul(tau, second[-1]), 2), second[-2])
        )
    return first[: depth + 1], second[: depth + 1]


def pell_numbers(limit: int) -> list[int]:
    values = [0, 1]
    while len(values) <= limit:
        values.append(2 * values[-1] + values[-2])
    return values


def triangular(index: int) -> int:
    return index * (index + 1) // 2


def pell_update(x_value: int, q_value: int) -> tuple[int, int]:
    return 3 * x_value + 8 * q_value, x_value + 3 * q_value


def evaluated_update(g_value: int, h_value: int) -> tuple[int, int]:
    return 3 * g_value + 2 * h_value, 4 * g_value + 3 * h_value


def distance(speed: int, phase: Fraction) -> Fraction:
    residue = (speed * phase.numerator) % phase.denominator
    return Fraction(min(residue, phase.denominator - residue), phase.denominator)


def mobius(value: Fraction) -> Fraction:
    return Fraction(1, 1) / (4 - 2 * value)


def factor_integer(value: int) -> dict[int, int]:
    remaining = value
    prime = 2
    factors: dict[int, int] = {}
    while prime * prime <= remaining:
        while remaining % prime == 0:
            factors[prime] = factors.get(prime, 0) + 1
            remaining //= prime
        prime = 3 if prime == 2 else prime + 2
    if remaining > 1:
        factors[remaining] = factors.get(remaining, 0) + 1
    return factors


def inert_cube_free_pair(pair: tuple[int, int]) -> bool:
    left, right = pair
    common = gcd(left, right)
    primitive_sum = (left + right) // common
    total_factors = factor_integer(left + right)
    primitive_factors = factor_integer(primitive_sum)
    return all(prime % 3 == 2 for prime in total_factors) and all(
        exponent <= 2 for exponent in primitive_factors.values()
    )


maximum_pell_index = max(2 * PELL_DEPTH + 2, 4 * LRC_DEPTH + 5)
PELL = pell_numbers(maximum_pell_index)

semantic_rows: list[str] = []

# Polynomial Chebyshev/Pell identity and evaluation.
FIRST, SECOND = chebyshev_tables(CHEBYSHEV_DEPTH)
tau_minus_one = [-2, 2]
tau_plus_one = [0, 2]
z_poly = [0, 1]
z_minus_one = [-1, 1]
polynomial_rows: list[tuple[int, int, int]] = []
for depth in range(1, CHEBYSHEV_DEPTH + 1):
    cheb_first = FIRST[depth]
    cheb_second = SECOND[depth - 1]
    g_poly = poly_add(cheb_first, poly_mul(tau_minus_one, cheb_second))
    h_poly = poly_add(cheb_first, poly_mul(tau_plus_one, cheb_second))
    pell_identity = poly_sub(
        poly_mul(z_poly, poly_mul(g_poly, g_poly)),
        poly_mul(z_minus_one, poly_mul(h_poly, h_poly)),
    )
    gate(pell_identity == [1], f"polynomial Pell identity at depth {depth}")
    g_value = poly_eval(g_poly, 2)
    h_value = poly_eval(h_poly, 2)
    gate(g_value == PELL[2 * depth + 1], "G(2) Pell indexing")
    gate(h_value == PELL[2 * depth + 1] + PELL[2 * depth],
         "H(2) Pell indexing")
    gate(h_value * h_value - 2 * g_value * g_value == -1,
         "evaluated negative Pell norm")
    polynomial_rows.append((depth, g_value, h_value))

# The perturbation keeps the same evaluation and destroys the polynomial law.
g_one = [-3, 4]
h_one = [-1, 4]
g_perturbed = poly_add(g_one, [-2, 1])
perturbed_identity = poly_sub(
    poly_mul(z_poly, poly_mul(g_perturbed, g_perturbed)),
    poly_mul(z_minus_one, poly_mul(h_one, h_one)),
)
perturbed_target = poly_mul(z_minus_one, [-1, -17, 9])
gate(poly_eval(g_one, 2) == poly_eval(g_perturbed, 2) == 5,
     "evaluation hostile keeps G(2)")
gate(poly_eval(h_one, 2) == 7, "evaluation hostile keeps H(2)")
gate(perturbed_identity == perturbed_target,
     "evaluation hostile polynomial invoice")
gate(perturbed_identity != [1], "evaluation hostile destroys Pell law")

# All integer state, forest-return, Pythagorean, and defect-shape gates.
first_rows: list[tuple[int, int, int, int, int, int]] = []
return_times: list[int] = []
for depth in range(1, PELL_DEPTH + 1):
    x_value = PELL[2 * depth] + PELL[2 * depth - 1]
    gate(PELL[2 * depth] % 2 == 0, "even Pell coordinate parity")
    q_value = PELL[2 * depth] // 2
    g_value = PELL[2 * depth + 1]
    h_value = PELL[2 * depth + 1] + PELL[2 * depth]
    a_value = PELL[2 * depth - 1]
    b_value = g_value
    gate(g_value % 2 == 1 and h_value % 2 == 1,
         "evaluated coordinates are odd")
    gate(x_value == 2 * g_value - h_value, "recover x from evaluation")
    gate(q_value == (h_value - g_value) // 2,
         "recover q from evaluation")
    gate(2 * q_value == h_value - g_value, "q recovery parity")
    gate(a_value == 3 * g_value - 2 * h_value,
         "recover left Markov factor")
    gate(b_value == g_value, "recover right Markov factor")
    gate(h_value * h_value - 2 * g_value * g_value == -1,
         "negative Pell norm")
    gate(x_value * x_value - 8 * q_value * q_value == 1,
         "positive Pell norm")
    rank = (x_value + 1) // 2
    gate(2 * rank - 1 == x_value, "rank integrality")
    delta = triangular(rank - 1)
    gate(delta == q_value * q_value, "square triangular defect")
    triple = (x_value, 4 * q_value * q_value, 4 * q_value * q_value + 1)
    gate(triple[0] * triple[0] + triple[1] * triple[1]
         == triple[2] * triple[2], "Pythagorean image")
    gate(gcd(gcd(triple[0], triple[1]), triple[2]) == 1,
         "primitive Pythagorean image")
    gate(triple[1] + triple[2] == x_value * x_value,
         "outer odd-square root")
    gate(triple[2] - triple[1] == 1, "unit inner odd-square root")
    tree_depth = rank - 2
    conductor_exponent = rank - 1
    gate(conductor_exponent == tree_depth + 1,
         "conductor exponent versus forest depth")
    gate(delta == conductor_exponent * rank // 2,
         "defect staircase area")
    if depth < PELL_DEPTH:
        next_x, next_q = pell_update(x_value, q_value)
        next_g, next_h = evaluated_update(g_value, h_value)
        next_rank = (next_x + 1) // 2
        gate(next_x == 2 * next_g - next_h, "updated evaluation inverse")
        gate(next_q * 2 == next_h - next_g, "updated q inverse")
        gate(next_rank == rank + h_value, "forest return exponent")
        gate(x_value + 2 * h_value == next_x,
             "odd-root L-power return")
        gate(a_value == 3 * g_value - 2 * h_value,
             "Markov source coordinate")
        gate((b_value, 6 * b_value - a_value)
             == (PELL[2 * depth + 1], PELL[2 * depth + 3]),
             "Markov mutation")
        appended_count = next_rank - rank
        appended_length = appended_count * (rank + next_rank - 1) // 2
        next_delta = triangular(next_rank - 1)
        gate(appended_count == h_value,
             "defect summand count equals return time")
        gate(appended_length == next_delta - delta,
             "defect summand length invoice")
    if depth <= 5:
        first_rows.append((depth, rank, q_value, x_value, g_value, h_value))
        return_times.append(h_value)

# LRC scalar observer and its alternating Mobius semiconjugacy.
lrc_rows: list[tuple[int, Fraction, Fraction]] = []
for depth in range(1, LRC_DEPTH + 1):
    x_value = PELL[2 * depth] + PELL[2 * depth - 1]
    q_value = PELL[2 * depth] // 2
    a_value = x_value - 2 * q_value
    b_value = x_value + 2 * q_value
    minus_index = 4 * depth - 1
    plus_index = 4 * depth + 1
    minus_direct = Fraction(
        PELL[minus_index] - PELL[minus_index - 1] + 1,
        2 * (PELL[minus_index] + 1),
    )
    plus_direct = Fraction(
        PELL[plus_index] - PELL[plus_index - 1] + 1,
        2 * (PELL[plus_index] + 1),
    )
    minus_formula = Fraction(a_value, x_value)
    plus_formula = Fraction(x_value, 2 * b_value)
    gate(minus_direct == minus_formula, "minus LRC factor formula")
    gate(plus_direct == plus_formula, "plus LRC factor formula")
    gate(mobius(minus_formula) == plus_formula,
         "first Mobius half-step")
    next_x, next_q = pell_update(x_value, q_value)
    next_minus = Fraction(next_x - 2 * next_q, next_x)
    gate(mobius(plus_formula) == next_minus,
         "second Mobius half-step")
    if depth <= 3:
        lrc_rows.append((depth, minus_formula, plus_formula))

# Same optimum, different labelled carry packet.
phase_two = Fraction(1, 3)
packet_two = tuple(distance(speed, phase_two) for speed in (1, 2))
packet_three = tuple(distance(speed, phase_two) for speed in (1, 2, 5))
gate(min(packet_two) == min(packet_three) == Fraction(1, 3),
     "scalar observer hostile maximum")
gate(len(packet_two) != len(packet_three),
     "scalar observer hostile packet cardinality")

# Mod-thirteen clock conjugacy and state-dependent forest cocycle.
norm_one_states = [
    (x_value, q_value)
    for x_value in range(MODULUS)
    for q_value in range(MODULUS)
    if (x_value * x_value - 8 * q_value * q_value) % MODULUS == 1
]
gate(len(norm_one_states) == 14, "mod-thirteen norm-one fibre size")
for initial_x, initial_q in norm_one_states:
    x_value, q_value = initial_x, initial_q
    visited: list[tuple[int, int]] = []
    half_sum = 0
    full_sum = 0
    for step in range(14):
        visited.append((x_value, q_value))
        g_value = (x_value + 2 * q_value) % MODULUS
        h_value = (x_value + 4 * q_value) % MODULUS
        next_x = (3 * x_value + 8 * q_value) % MODULUS
        next_q = (x_value + 3 * q_value) % MODULUS
        next_g = (3 * g_value + 2 * h_value) % MODULUS
        next_h = (4 * g_value + 3 * h_value) % MODULUS
        gate(next_g == (next_x + 2 * next_q) % MODULUS,
             "mod-thirteen evaluated conjugacy G")
        gate(next_h == (next_x + 4 * next_q) % MODULUS,
             "mod-thirteen evaluated conjugacy H")
        gate(next_x == (x_value + 2 * h_value) % MODULUS,
             "mod-thirteen forest cocycle")
        if step < 7:
            half_sum = (half_sum + h_value) % MODULUS
        full_sum = (full_sum + h_value) % MODULUS
        x_value, q_value = next_x, next_q
        if step == 6:
            gate((x_value, q_value)
                 == ((-initial_x) % MODULUS, (-initial_q) % MODULUS),
                 "mod-thirteen central-sign half-period")
    gate(len(set(visited)) == 14, "mod-thirteen sharp orbit")
    gate((x_value, q_value) == (initial_x, initial_q),
         "mod-thirteen full period")
    gate(half_sum == (-initial_x) % MODULUS,
         "mod-thirteen half cocycle sum")
    gate(full_sum == 0, "mod-thirteen full cocycle sum")

# Native-operation and coordinate-type hostiles.
gate(triangular(1) == 1 and triangular(2) == 3,
     "single L-edge square-selector hostile")
gate(return_times[:2] == [7, 41], "nonconstant return-time hostile")
gate(pell_update(3, 1) == (17, 6), "Pell sidecar coordinate hostile")
gate((3 + 2 * 7, 1) == (17, 1), "forest inner-root hostile")
gate((17, 6) != (17, 1), "second-coordinate semantics diverge")

# Actual monomial conductor rings at consecutive marked returns are
# incomparable inside their common normalization k[t].
def in_two_generator_semigroup(value: int, first: int, second: int) -> bool:
    return any((value - count * first) % second == 0
               for count in range(value // first + 1))


gate(in_two_generator_semigroup(50, 9, 10), "t^50 belongs to A_9")
gate(not in_two_generator_semigroup(51, 9, 10),
     "t^51 does not belong to A_9")
gate(not in_two_generator_semigroup(9, 50, 51),
     "t^9 does not belong to A_50")

# Exact intersection with THM-3793's l1<=356 support-two ratios.
pell_cap_pairs: list[tuple[int, int]] = []
depth = 1
while True:
    pair = (PELL[2 * depth - 1], PELL[2 * depth + 1])
    if sum(pair) > 356:
        break
    pell_cap_pairs.append(pair)
    depth += 1
gate(pell_cap_pairs == [(1, 5), (5, 29), (29, 169)],
     "complete Pell intersection with support-two cap")
admissible_pairs = [pair for pair in pell_cap_pairs if inert_cube_free_pair(pair)]
gate(admissible_pairs == [(5, 29)], "inert cube-free Pell intersection")
cube_address = 5**3 + 29**3
representations: list[tuple[int, int]] = []
for right in range(2, 100):
    for left in range(1, right):
        if left**3 + right**3 == cube_address:
            representations.append((left, right))
gate(cube_address == 24_514, "two-cube address value")
gate(representations == [(5, 29)], "two-cube singleton fibre")
gate(78 * len(admissible_pairs) == 78, "labelled placement count")

semantic_rows.extend(
    [
        f"poly={polynomial_rows}",
        f"first={first_rows}",
        f"lrc={lrc_rows}",
        f"mod13={norm_one_states}",
        f"cap={pell_cap_pairs}",
        f"admissible={admissible_pairs}",
        f"cube={cube_address}:{representations}",
        f"perturb={perturbed_identity}",
    ]
)
semantic_digest = sha256("\n".join(semantic_rows).encode("ascii")).hexdigest()

print("THM-3819 EXACT COMPANION")
print(
    "UNIVERSE "
    f"chebyshev_depth={CHEBYSHEV_DEPTH} "
    f"pell_depth={PELL_DEPTH} "
    f"lrc_depth={LRC_DEPTH} "
    f"mod13_norm_one={len(norm_one_states)}"
)
print("FIRST_ROWS n:m:q:x:g:h")
for row in first_rows:
    print(":".join(str(entry) for entry in row))
print("RETURN_TIMES " + ",".join(str(value) for value in return_times))
print(
    "LRC_ROWS "
    + ";".join(
        f"{depth}:{minus.numerator}/{minus.denominator},"
        f"{plus.numerator}/{plus.denominator}"
        for depth, minus, plus in lrc_rows
    )
)
print("MOD13 half_sum=-x full_sum=0 states=14 PASS")
print(f"SUPPORT2 pairs={pell_cap_pairs} inert={admissible_pairs}")
print(f"CUBE_ADDRESS {cube_address} representations={representations}")
print(
    "HOSTILES selector_nonclosure,fixed_word,coordinate_type,"
    "evaluation_loss,ring_incomparability,scalar_packet PASS"
)
print(f"CHECKS={CHECKS}")
print(f"SEMANTIC_SHA256={semantic_digest}")
print("STATUS=PASS")
