#!/usr/bin/env python3
"""Exact referee for THM-3053 using only integer and rational arithmetic."""

from fractions import Fraction
from itertools import product
from math import factorial


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def rising(x, length):
    out = Fraction(1)
    for step in range(length):
        out *= x + step
    return out


def determinant(matrix):
    a = [[Fraction(entry) for entry in row] for row in matrix]
    size = len(a)
    det = Fraction(1)
    for column in range(size):
        pivot = next((row for row in range(column, size) if a[row][column]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            a[column], a[pivot] = a[pivot], a[column]
            det = -det
        value = a[column][column]
        det *= value
        for j in range(column, size):
            a[column][j] /= value
        for row in range(column + 1, size):
            scale = a[row][column]
            if scale:
                for j in range(column, size):
                    a[row][j] -= scale * a[column][j]
    return det


def prefix_sums(inventory):
    total = 0
    out = []
    for entry in inventory:
        total += entry
        out.append(total)
    return out


def direct_moment(a, inventory, m):
    out = Fraction(1)
    for j, exponent in enumerate(inventory):
        out *= rising(a + j, m) ** exponent
    return out


def adjacent_beta_gamma_moment(a, inventory, m):
    """Canonical factorization using prefix copies of Beta(a+j,1)."""
    prefixes = prefix_sums(inventory)
    require(all(value >= 0 for value in prefixes), "prefix-infeasible inventory")
    last = len(inventory) - 1
    out = rising(a + last, m) ** prefixes[last]
    for j in range(last):
        out *= (rising(a + j, m) / rising(a + j + 1, m)) ** prefixes[j]
    return out


def greedy_long_edges(inventory):
    """Pair each negative unit with an earlier positive unit."""
    supply = []
    beta_edges = []
    for j, entry in enumerate(inventory):
        if entry > 0:
            supply.extend([j] * entry)
        else:
            for _ in range(-entry):
                if not supply:
                    return None
                beta_edges.append((supply.pop(), j))
    return beta_edges, tuple(supply)


def ledger_from_transport(length, transport):
    edges, gamma_levels = transport
    ledger = [0] * length
    for left, right in edges:
        ledger[left] += 1
        ledger[right] -= 1
    for level in gamma_levels:
        ledger[level] += 1
    return tuple(ledger)


print("THM-3053 BETA-GAMMA PREFIX TRANSPORT AND MULTIPLICATIVE HOLOTOPY")

# Exhaust the small signed-inventory box.  Prefix feasibility, greedy ordered
# transport, and the canonical adjacent factorization agree exactly.
inventory_cells = 0
feasible_cells = 0
for length in range(1, 6):
    for inventory in product(range(-3, 4), repeat=length):
        prefixes = prefix_sums(inventory)
        feasible = all(value >= 0 for value in prefixes)
        transport = greedy_long_edges(inventory)
        require((transport is not None) == feasible, "prefix/transport equivalence failed")
        if feasible:
            require(ledger_from_transport(length, transport) == inventory,
                    "transport divergence failed")
            feasible_cells += 1
        inventory_cells += 1
print(f"inventory_cells={inventory_cells} prefix_feasible={feasible_cells}")

# Materialized moment identities and strict generalized-Hankel controls.
moment_cells = 0
hankel_cells = 0
representatives = (
    (1,), (0, 1), (2, -1), (1, 0, 2), (3, -2, 1),
    (2, -1, 0, 2), (1, 2, -2, 1), (3, -1, -1, 2),
)
for a in (Fraction(1, 2), Fraction(1), Fraction(5, 3)):
    for inventory in representatives:
        require(all(value >= 0 for value in prefix_sums(inventory)), "bad representative")
        sequence = []
        for m in range(13):
            direct = direct_moment(a, inventory, m)
            factored = adjacent_beta_gamma_moment(a, inventory, m)
            require(direct == factored, "canonical Beta-Gamma moment identity failed")
            sequence.append(direct)
            moment_cells += 1
        for size in (1, 2, 3):
            for row_shift in range(3):
                for column_shift in range(3):
                    block = [[sequence[row_shift + column_shift + i + j]
                              for j in range(size)] for i in range(size)]
                    require(determinant(block) > 0, "strict Hankel control failed")
                    hankel_cells += 1
print(f"moment_identity_cells={moment_cells} strict_hankel_cells={hankel_cells}")

# Beta edge subdivision/coalescence: B_(i,j) B_(j,k) and B_(i,k) have
# identical Hausdorff moments.
clutch_cells = 0
for a in (Fraction(2, 3), Fraction(1), Fraction(7, 4)):
    for i in range(4):
        for j in range(i + 1, 5):
            for k in range(j + 1, 6):
                for m in range(9):
                    left = rising(a + i, m) / rising(a + j, m)
                    left *= rising(a + j, m) / rising(a + k, m)
                    right = rising(a + i, m) / rising(a + k, m)
                    require(left == right, "Beta clutch identity failed")
                    clutch_cells += 1
print(f"beta_clutch_cells={clutch_cells}")

# Specialize to THM-3047's adjacent transfer exponents.  Prefix feasibility
# telescopes to d_0<=A and d_m<=I for m>=1.
transfer_cells = 0
for slot_count in range(2, 6):
    harmonic = sum(Fraction(1, j) for j in range(1, slot_count + 1))
    a_count = int(factorial(slot_count) * (harmonic - 1))
    b_count = int(factorial(slot_count) * (slot_count + 1 - 2 * harmonic))
    interior = a_count + b_count
    for d in product(range(-2, 5), repeat=4):
        e = [a_count, b_count, 0, 0, 0]
        extended_d = list(d) + [0]
        inventory = []
        previous = 0
        for j in range(5):
            inventory.append(e[j] + previous - extended_d[j])
            previous = extended_d[j]
        feasible = all(value >= 0 for value in prefix_sums(inventory))
        bounds = d[0] <= a_count and all(value <= interior for value in d[1:])
        require(feasible == bounds, "THM-3047 transfer-prefix criterion failed")
        transfer_cells += 1
print(f"formal_width_transfer_cells={transfer_cells}")

# Prefix nonnegativity is not necessary for Stieltjes representability once
# convex mixtures, rather than products alone, are admitted.
mixture_escape_cells = 0
hausdorff_cells = 0
for a in (Fraction(1, 2), Fraction(1), Fraction(5, 3), Fraction(7, 2)):
    sequence = []
    for m in range(11):
        direct = direct_moment(a, (1, -2, 1), m)
        mixture = a / (a + 1) + Fraction(1, 1) / (a + 1) * a / (a + m)
        require(direct == mixture, "convex-mixture prefix escape failed")
        sequence.append(direct)
        mixture_escape_cells += 1
    for size in (2, 3):
        require(determinant([[sequence[i + j] for j in range(size)] for i in range(size)]) > 0,
                "convex-mixture escape lost strict Hankel positivity")
    for start in range(6):
        for order in range(1, 6):
            difference = sum((-1) ** (order - j) * factorial(order) // (factorial(j) * factorial(order - j))
                             * sequence[start + j] for j in range(order + 1))
            require((-1) ** order * difference > 0,
                    "convex-mixture escape lost Hausdorff complete monotonicity")
            hausdorff_cells += 1
print(f"mixture_escape_cells={mixture_escape_cells} hausdorff_cells={hausdorff_cells}")

# Sharp first-negative-prefix controls.  The first fails H2 immediately; the
# second remains strictly log-convex through H2 but fails H3.
a = Fraction(1)
affine_hostile = [direct_moment(a, (-1, 1), m) for m in range(3)]
affine_h2 = determinant([[affine_hostile[i + j] for j in range(2)] for i in range(2)])
require(affine_h2 == -1, "first-prefix H2 hostile changed")

hidden_hostile = [direct_moment(a, (-1, 2), m) for m in range(5)]
hidden_h2 = determinant([[hidden_hostile[i + j] for j in range(2)] for i in range(2)])
hidden_h3 = determinant([[hidden_hostile[i + j] for j in range(3)] for i in range(3)])
require(hidden_h2 == 2 and hidden_h3 == -24, "first-prefix H3 hostile changed")
for m in range(1, 4):
    require(hidden_hostile[m - 1] * hidden_hostile[m + 1] > hidden_hostile[m] ** 2,
            "hidden hostile lost strict log-convexity")
print(f"negative_prefix_hostiles=H2:{affine_h2},hidden_H2:{hidden_h2},hidden_H3:{hidden_h3}")
print("all_exact_checks=PASS")
