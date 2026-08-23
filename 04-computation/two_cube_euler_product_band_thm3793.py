#!/usr/bin/env python3
"""Exact controls for the two-cube Euler-product mesoscopic band."""

from __future__ import annotations

import hashlib
import math
import sys
from fractions import Fraction


sys.stdout.reconfigure(newline="\n")
GATES = 0


def require(condition: bool, label: str) -> None:
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(label)


def primes_to(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for p in range(2, math.isqrt(limit) + 1):
        if sieve[p]:
            sieve[p * p : limit + 1 : p] = b"\x00" * (((limit - p * p) // p) + 1)
    return [n for n in range(2, limit + 1) if sieve[n]]


def elementary_coefficients(primes: list[int]) -> list[Fraction]:
    coefficients = [Fraction(1)]
    for p in primes:
        coefficients.append(Fraction(0))
        for j in range(len(coefficients) - 1, 0, -1):
            coefficients[j] += coefficients[j - 1] / p
    return coefficients


def bernoulli_order_probabilities(primes: list[int]) -> list[Fraction]:
    probabilities = [Fraction(1)]
    for p in primes:
        q = Fraction(1, p + 1)
        old = probabilities
        probabilities = [Fraction(0)] * (len(old) + 1)
        for j, value in enumerate(old):
            probabilities[j] += value * (1 - q)
            probabilities[j + 1] += value * q
    return probabilities


# Exact Euler-product/Bernoulli normalization.  This is the native law for
# the complete sum of elementary symmetric layers, unlike the artificial
# Poisson law used by the collision lower bound.
bank = [p for p in primes_to(101) if p >= 5 and p % 3 == 2]
elementary = elementary_coefficients(bank)
probabilities = bernoulli_order_probabilities(bank)
euler = Fraction(1)
for p in bank:
    euler *= Fraction(p + 1, p)
require(sum(elementary, Fraction(0)) == euler, "elementary Euler product")
require(sum(probabilities, Fraction(0)) == 1, "Bernoulli probability mass")
for j in range(len(elementary)):
    require(euler * probabilities[j] == elementary[j], f"Bernoulli layer j={j}")

mean = sum((Fraction(1, p + 1) for p in bank), Fraction(0))
variance = sum((Fraction(p, (p + 1) ** 2) for p in bank), Fraction(0))
mean_from_distribution = sum((j * probabilities[j] for j in range(len(probabilities))), Fraction(0))
second = sum((j * j * probabilities[j] for j in range(len(probabilities))), Fraction(0))
require(mean_from_distribution == mean, "Bernoulli mean")
require(second - mean * mean == variance, "Bernoulli variance")
require(variance <= mean, "variance bounded by mean")

for J in range(1, len(bank) + 1):
    left = sum(elementary[1 : J + 1], Fraction(0))
    right = euler * sum(probabilities[1 : J + 1], Fraction(0))
    require(left == right, f"truncated Euler identity J={J}")

# Floor-sensitive asymptotic ledger.  Williams' Mertens theorem gives
# mu(Z)=0.5*log log Z+O(1); these rows verify the deterministic displacement
# J_+ - 0.5*log log Z and its Chebyshev scale.
cutoff_rows = []
for L in (10**6, 10**7, 10**8, 10**9, 10**10):
    J = math.floor(L / 2 - 0.5 * math.log(L) + L ** (2 / 3))
    half_loglog_z = 0.5 * (L - math.log(3 * J))
    displacement = J - half_loglog_z
    require(J > L / 2, f"cutoff above half L={L}")
    require(displacement >= L ** (2 / 3) - 1, f"Mertens-mean displacement L={L}")
    require(displacement / math.sqrt(L) > 10, f"many standard deviations L={L}")
    require(L / (displacement * displacement) < 0.011, f"Chebyshev decay L={L}")
    cutoff_rows.append((L, J, half_loglog_z, displacement, L / displacement**2))

# Independent finite approximants to the Williams constant.  Q includes p=2;
# E excludes p=2, exactly as in THM-3793.  The theorem uses the limit and the
# cited Mertens product, not these bounded decimals.
EULER_GAMMA = 0.577215664901532860606512090082402431
L_CHI_MINUS3 = math.pi / (3 * math.sqrt(3))
numeric_rows = []
for limit in (10_000, 100_000, 1_000_000):
    q_product = 1.0
    e_product = 1.0
    for p in primes_to(limit):
        if p % 3 != 2:
            continue
        q_product *= 1 - 1 / (p * p)
        if p >= 5:
            e_product *= 1 + 1 / p
    constant_from_q = (2 / 3) * math.sqrt(
        2 * math.exp(EULER_GAMMA) * q_product / (3 * L_CHI_MINUS3)
    )
    constant_direct = e_product / math.sqrt(math.log(limit))
    critical_constant = (2 / 5) * math.sqrt(2 / 3) * constant_from_q
    require(abs(constant_direct - constant_from_q) < 0.0002, f"constant routes limit={limit}")
    numeric_rows.append((limit, q_product, constant_from_q, constant_direct, critical_constant))

require(0.7854 < numeric_rows[-1][2] < 0.7859, "Euler constant numerical control")
require(0.2564 < numeric_rows[-1][4] < 0.2568, "critical constant numerical control")

semantic = (
    "sum_j e_j is an Euler product",
    "normalized elementary layers are independent Bernoulli orders with q_p=1/(p+1)",
    "mean and variance are O(log log Z)",
    "J_plus exceeds the mean by L^(2/3)+O(1)",
    "Williams Mertens constant replaces the collision-union-bound constant",
)
semantic_sha256 = hashlib.sha256(repr(semantic).encode()).hexdigest()
print("TWO_CUBE_EULER_PRODUCT_BAND_PRIMARY_20260823")
print("status=PROVED_GIVEN_CITED_WILLIAMS_MERTENS_AP;no_support_asymptotic")
print(f"finite_bank_size={len(bank)};euler={euler};mean={mean};variance={variance}")
print("cutoff_rows=" + repr(tuple((row[0], row[1], round(row[3], 6), f"{row[4]:.6e}") for row in cutoff_rows)).replace(" ", ""))
print("numeric_rows=" + repr(tuple((row[0], f"{row[2]:.12f}", f"{row[3]:.12f}", f"{row[4]:.12f}") for row in numeric_rows)).replace(" ", ""))
print(f"semantic_sha256={semantic_sha256}")
print(f"gates={GATES}")
print("ALL CHECKS PASSED")
