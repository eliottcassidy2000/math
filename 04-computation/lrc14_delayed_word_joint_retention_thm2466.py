#!/usr/bin/env python3
"""Exact finite companion for THM-2466.

The mathematical proof is the BV/Fourier argument in the theorem.  This
companion uses only integers and Fraction arithmetic to replay its constants,
representative dilation integrals, root-clock identity, and sharp hostiles.
Every truth-bearing check uses ``require`` and remains active under
``python -O``.
"""

from fractions import Fraction
from math import lcm


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


P = 13


def cyclic_variation(values):
    """Total variation of a function constant on equal cyclic cells."""
    return sum(
        abs(values[(i + 1) % len(values)] - values[i])
        for i in range(len(values))
    )


def mean(values):
    return sum(values, Fraction(0)) / len(values)


def dilated_product_integral(left, word, dilation):
    """Integral of left(x) word(dilation*x) for equal-cell step data."""
    require(dilation >= 1, "dilation must be positive")
    n_left = len(left)
    n_word = len(word)
    mesh = lcm(n_left, n_word * dilation)
    total = Fraction(0)
    for cell in range(mesh):
        left_index = (cell * n_left) // mesh
        # Both functions are constant on this common half-open mesh cell.
        word_index = ((cell * n_word * dilation) // mesh) % n_word
        total += left[left_index] * word[word_index]
    return total / mesh


print("THM-2466 DELAYED-WORD JOINT RETENTION -- exact finite audit")

# One rational Boolean word with positive mean and finite cyclic variation.
word = tuple(Fraction(x) for x in (1, 0, 1, 0, 0, 1, 0))
mu = mean(word)
v_word = cyclic_variation(word)
require(mu == Fraction(3, 7), "word mean")
require(v_word == 6, "word variation")
print(f"word: mean={mu}, variation={v_word}")

# Fixed signed drift and nonnegative service densities.  Denominator ten is
# coprime to the word denominator and to thirteen, so the common meshes are
# genuine rather than aligned toy grids.
drift = tuple(
    Fraction(x) for x in (2, -1, 3, 0, -2, 4, 1, -1, 2, 1)
)
service = tuple(
    Fraction(x) for x in (0, 2, 1, 4, 0, 3, 1, 0, 2, 5)
)
q = mean(drift)
mass = mean(service)
v_drift = cyclic_variation(drift)
v_service = cyclic_variation(service)
require(q == Fraction(9, 10), "drift mean")
require(mass == Fraction(9, 5), "service mean")
require(q != 0 and mass > 0, "positive marginals")

rows = []
for k in range(1, 5):
    root_N = P ** (k - 1)
    q_k = dilated_product_integral(drift, word, root_N)
    m_k = dilated_product_integral(service, word, root_N)
    drift_error = abs(q_k - mu * q)
    service_error = abs(m_k - mu * mass)
    drift_bound = Fraction(v_word * v_drift, 12 * root_N)
    service_bound = Fraction(v_word * v_service, 12 * root_N)
    drift_discrepancy_bound = mu * (1 - mu) * v_drift / root_N
    service_discrepancy_bound = mu * (1 - mu) * v_service / root_N
    require(drift_error <= drift_bound, "drift BV inequality")
    require(service_error <= service_bound, "service BV inequality")
    require(drift_error <= drift_discrepancy_bound,
            "drift Boolean-discrepancy inequality")
    require(service_error <= service_discrepancy_bound,
            "service Boolean-discrepancy inequality")
    rows.append((k, q_k, m_k, drift_error, service_error))

print(
    "fixed densities:",
    f"drift mean={q}, variation={v_drift};",
    f"service mean={mass}, variation={v_service}",
)
for k, q_k, m_k, drift_error, service_error in rows:
    print(
        f"k={k}: q_k={q_k}, M_k={m_k},",
        f"errors=({drift_error},{service_error})",
    )
require(all(q_k != 0 and m_k > 0 for _, q_k, m_k, _, _ in rows[1:]),
        "one delayed word did not retain both observables")

# Exact root constancy: Q(13^k (y+u)/13)=Q(13^(k-1)y).
root_checks = 0
for k in range(1, 7):
    for y_num in (1, 8, 17, 41):
        y = Fraction(y_num, 113)
        for u in range(P):
            left = P**k * (y + u) / P
            right = P ** (k - 1) * y
            left_frac = left - left.numerator // left.denominator
            right_frac = right - right.numerator // right.denominator
            require(left_frac == right_frac, "root-constant clock identity")
            root_checks += 1
print(f"root-constant one-power-offset checks: {root_checks}")

# A fixed-clock hostile which later dilation repairs.  At N=1 the density is
# exactly the complement of the word.  Freeze it, then change only the clock.
complement = tuple(1 - x for x in word)
require(mean(complement) == Fraction(4, 7), "hostile positive mean")
require(dilated_product_integral(complement, word, 1) == 0,
        "fixed-clock hostile should vanish")
recovered = []
for n in (P, P**2, P**3):
    value = dilated_product_integral(complement, word, n)
    require(value > 0, "frozen hostile failed to recover after delay")
    recovered.append(value)
print("fixed-clock hostile: zero at N=1; frozen-density recovery:", *recovered)

# If the density may change with the clock, positivity never forces a hit.
for k in range(1, 5):
    n = P ** (k - 1)
    # On the full 7n-cell mesh this is literally Q(ny).
    moving_word = tuple(word[i % len(word)] for i in range(len(word) * n))
    moving_complement = tuple(1 - x for x in moving_word)
    require(mean(moving_complement) == 1 - mu, "moving density mean")
    require(dilated_product_integral(moving_complement, word, n) == 0,
            "full-mesh clock-dependent hostile")
    pointwise_product = [
        moving_word[i] * moving_complement[i]
        for i in range(len(moving_word))
    ]
    require(sum(pointwise_product) == 0, "clock-dependent hostile")
print("clock-dependent hostile: positive mean and zero hit at four clocks")

# THM-2459 -> THM-2466 -> THM-2457 invoices.
drift_den = 4 * 63001
service_den = 2 * 16384
root_energy_den = service_den**2 * 342732
root_coefficient_den = service_den * 2028
require(drift_den == 252004, "joint drift denominator")
require(service_den == 32768, "joint service denominator")
require(root_energy_den == 368005682823168, "root energy denominator")
require(root_coefficient_den == 66453504, "root coefficient denominator")
print(
    "THM-2459/2457 invoice:",
    f"D >= mu^2 D0/{drift_den},",
    f"M >= mu M0/{service_den},",
    f"root energy denominator={root_energy_den},",
    f"root coefficient denominator={root_coefficient_den}",
)

# Target neutrality and both cofinal septimal parity classes.
parities = set()
for k in range(1, 9):
    R = P**k
    require(R % P == 0, "target-neutral delayed clock")
    residue = R % 7
    require(residue == (6 if k % 2 else 1), "septimal parity")
    parities.add(residue)
require(parities == {1, 6}, "both parity classes")
print("clock gauge: target-neutral; septimal residues alternate 6,1 cofinally")

print("fixed delay is not prescribed first expiration")
print("common oriented root section remains a required sidecar")
print("all exact checks passed")
