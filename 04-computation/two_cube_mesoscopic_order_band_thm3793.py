#!/usr/bin/env python3
"""Exact and asymptotic controls for the two-cube mesoscopic order band."""

from __future__ import annotations

import hashlib
import itertools
import math
import sys
from fractions import Fraction


sys.stdout.reconfigure(encoding="utf-8", newline="\n")
GATES = 0


def gate(condition: bool, label: str) -> None:
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(label)


def primes_through(bound: int) -> list[int]:
    result = []
    for value in range(2, bound + 1):
        if all(value % divisor for divisor in range(2, math.isqrt(value) + 1)):
            result.append(value)
    return result


def inert_bank(bound: int) -> list[int]:
    return [prime for prime in primes_through(bound) if prime >= 5 and prime % 3 == 2]


def elementary(values: list[Fraction], order: int) -> Fraction:
    coefficients = [Fraction(1)] + [Fraction(0)] * order
    for value in values:
        for degree in range(order, 0, -1):
            coefficients[degree] += value * coefficients[degree - 1]
    return coefficients[order]


def direct_elementary(values: list[Fraction], order: int) -> Fraction:
    return sum(
        (math.prod(choice) for choice in itertools.combinations(values, order)),
        Fraction(0),
    )


def log_poisson_band(mean: float, lower: int, upper: int) -> float:
    logs = [index * math.log(mean) - math.lgamma(index + 1) - mean for index in range(lower, upper + 1)]
    maximum = max(logs)
    return maximum + math.log(sum(math.exp(value - maximum) for value in logs))


# Complete positive-pair hostile control through a three-prime bank.  Every
# constructed value is checked against the entire x+y<=max_d universe, not
# merely against the construction family.
small_bank = inert_bank(17)
gate(small_bank == [5, 11, 17], "small inert bank")
products: list[tuple[int, int]] = []
for order in range(1, len(small_bank) + 1):
    for choice in itertools.combinations(small_bank, order):
        products.append((order, math.prod(choice)))
max_d = max(product for _, product in products)
pair_fibres: dict[int, list[tuple[int, int]]] = {}
for x in range(1, max_d):
    for y in range(x + 1, max_d + 1):
        value = x**3 + y**3
        pair_fibres.setdefault(value, []).append((x, y))

constructed: dict[int, tuple[int, int, int]] = {}
row_floor = Fraction(0)
elementary_floor = Fraction(0)
for order, pair_sum in products:
    row_floor += Fraction(2, 5 * pair_sum)
    elementary_floor += Fraction(2, 5 * pair_sum)
    row_surrogate = Fraction(0)
    for x in range(1, (pair_sum - 1) // 2 + 1):
        y = pair_sum - x
        value = x**3 + y**3
        gate(pair_fibres[value] == [(x, y)], f"complete singleton fibre d={pair_sum},x={x}")
        gate(value not in constructed, f"cross-layer disjointness m={value}")
        constructed[value] = (order, x, y)
        row_surrogate += Fraction(1, pair_sum**2)
    gate(row_surrogate == Fraction(pair_sum - 1, 2 * pair_sum**2), f"row count d={pair_sum}")
    gate(row_surrogate >= Fraction(2, 5 * pair_sum), f"row mass floor d={pair_sum}")

reciprocals = [Fraction(1, prime) for prime in small_bank]
sum_layers = sum(
    (Fraction(2, 5) * elementary(reciprocals, order) for order in range(1, 4)),
    Fraction(0),
)
gate(sum_layers == elementary_floor, "all-order elementary decomposition")
gate(len({order for order, _, _ in constructed.values()}) == 3, "three disjoint layers")

# Independent coefficient recursion/direct-combination agreement and the
# ordered-tuple collision lower bound on a larger bank.
large_bank = inert_bank(101)
large_reciprocals = [Fraction(1, prime) for prime in large_bank]
A = sum(large_reciprocals, Fraction(0))
B = sum((value * value for value in large_reciprocals), Fraction(0))
for order in range(0, min(8, len(large_bank)) + 1):
    exact = elementary(large_reciprocals, order)
    gate(exact == direct_elementary(large_reciprocals, order), f"elementary routes j={order}")
    if order >= 2:
        collision_lower = (
            A**order - math.comb(order, 2) * B * A ** (order - 2)
        ) / math.factorial(order)
        gate(exact >= collision_lower, f"collision union bound j={order}")

# The mesoscopic specialization.  The first cutoff reproduces the narrower
# audited window.  The second lies above the proxy mean by L^(2/3), so a
# Chebyshev bound captures asymptotically the entire artificial Poisson lower
# tail and gives the sharper canonical constant.
C_P = 14 / math.e + 0.5 + (math.log(3) + 2 / 9) / 2
normal_band = 0.5 * (math.erf(0 / math.sqrt(2)) - math.erf(-1 / math.sqrt(2)))
predicted = (7 / 20) * math.exp(-C_P) * math.sqrt(2 / 3) * normal_band
predicted_sharp = (7 / 20) * math.exp(-C_P) * math.sqrt(2 / 3)
samples = []
sharp_samples = []
for L in (1_000_000.0, 10_000_000.0, 100_000_000.0):
    J = math.floor(L / 2 - math.log(L))
    a = 0.5 * (L - math.log(3 * J)) - C_P
    lower = J - math.floor(math.sqrt(J))
    gate(2 <= lower <= J <= a, f"mesoscopic ordering L={L}")
    collision_factor = 1 - J * J / (8 * a * a)
    gate(collision_factor >= 7 / 8, f"uniform collision factor L={L}")
    probability = math.exp(log_poisson_band(a, lower, J))
    scaled = (2 / 5) * collision_factor * probability * math.exp(
        a - L / 2 + 0.5 * math.log(L)
    )
    samples.append((int(L), J, a, probability, scaled))

    J_plus = math.floor(L / 2 - 0.5 * math.log(L) + L ** (2 / 3))
    a_plus = 0.5 * (L - math.log(3 * J_plus)) - C_P
    theta = L / 2 - 0.5 * math.log(L) + L ** (2 / 3) - J_plus
    gap = J_plus - a_plus
    exact_gap = L ** (2 / 3) + 0.5 * math.log(3 * J_plus / L) + C_P - theta
    delta = 1 - J_plus * J_plus / (8 * a_plus * a_plus)
    poisson_tail_budget = a_plus / (gap * gap) + math.exp(-a_plus) * (1 + a_plus)
    probability_floor = 1 - poisson_tail_budget
    sharp_scaled = (
        (2 / 5)
        * delta
        * probability_floor
        * math.exp(-C_P)
        * math.sqrt(L / (3 * J_plus))
    )
    gate(L / 2 < J_plus <= 0.51 * L, f"sharp cutoff range L={L}")
    gate(0.49 * L <= a_plus < L / 2, f"sharp proxy range L={L}")
    gate(abs(gap - exact_gap) < 1e-8, f"sharp floor identity L={L}")
    gate(gap >= L ** (2 / 3), f"sharp upper-tail gap L={L}")
    gate(0 < delta < 7 / 8, f"sharp collision factor L={L}")
    gate(J_plus / a_plus < math.sqrt(8), f"sharp collision positivity L={L}")
    gate(0 < probability_floor < 1, f"sharp Chebyshev floor L={L}")
    sharp_samples.append((int(L), J_plus, a_plus, delta, poisson_tail_budget, sharp_scaled))
gate(samples[-1][-1] > predicted * 0.99, "scaled lower approaches predicted constant")
gate(abs(samples[-1][3] - normal_band) < 0.002, "Poisson band approaches normal interval")
gate(sharp_samples[-1][-1] > predicted_sharp * 0.98, "sharp lower approaches full-tail constant")
gate(predicted_sharp > 2.9 * predicted, "sharp constant strictly improves narrow band")

semantic = (
    "all squarefree inert prime-product orders give disjoint singleton layers",
    "H(Z^(3J)) is at least two-fifths times the sum of e_j for j<=J",
    "the collision lower bound is uniform for every j<=a_X",
    "J=L/2-log L lies a-o(sqrt(a)) below the Poisson mean a",
    "J_plus=L/2-(log L)/2+L^(2/3) captures the full artificial Poisson lower tail",
    "the mesoscopic band gives the explicit sharp audited constant times sqrt(log X/log log X)",
)
semantic_hash = hashlib.sha256(repr(semantic).encode()).hexdigest()
print("TWO_CUBE_MESOSCOPIC_BAND", f"small_bank={small_bank};constructed={len(constructed)};layers=3")
print("TWO_CUBE_EXACT_LAYER_FLOOR", str(sum_layers))
for L, J, a, probability, scaled in samples:
    print(
        "TWO_CUBE_ASYMPTOTIC_SAMPLE",
        f"L={L};J={J};a={a:.12f};band_probability={probability:.12f};scaled={scaled:.12e}",
    )
for L, J, a, delta, tail_budget, scaled in sharp_samples:
    print(
        "TWO_CUBE_SHARP_CUTOFF",
        f"L={L};J={J};a={a:.12f};delta={delta:.12f};chebyshev_tail={tail_budget:.12e};scaled={scaled:.12e}",
    )
print("TWO_CUBE_PREDICTED_LIMINF_CONSTANT", f"{predicted:.12e}")
print("TWO_CUBE_SHARP_LIMINF_CONSTANT", f"{predicted_sharp:.15e}")
print("SEMANTIC_SHA256", semantic_hash)
print("GATES", GATES)
print("RESULT", "PASS")
