#!/usr/bin/env python3
"""Import-independent hostile audit of the two-cube Euler-product band.

This companion deliberately does not import the candidate checker.  It uses
brute subset enumeration for the finite Bernoulli law and an independently
written incremental prime generator for the numerical controls.
"""

from __future__ import annotations

import hashlib
import math
import sys
from fractions import Fraction
from pathlib import Path


sys.stdout.reconfigure(newline="\n")
GATES = 0


def require(condition: bool, label: str) -> None:
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(label)


def incremental_primes(limit: int) -> list[int]:
    """Return primes through limit by trial division against earlier primes."""
    if limit < 2:
        return []
    answer = [2]
    for candidate in range(3, limit + 1, 2):
        prime = True
        for p in answer[1:]:
            if p * p > candidate:
                break
            if candidate % p == 0:
                prime = False
                break
        if prime:
            answer.append(candidate)
    return answer


def brute_layers(bank: list[int]) -> tuple[list[Fraction], list[Fraction]]:
    """Enumerate every subset, separately accumulating weights/probabilities."""
    weights = [Fraction(0) for _ in range(len(bank) + 1)]
    probabilities = [Fraction(0) for _ in range(len(bank) + 1)]
    for mask in range(1 << len(bank)):
        order = mask.bit_count()
        weight = Fraction(1)
        probability = Fraction(1)
        for i, p in enumerate(bank):
            if mask & (1 << i):
                weight *= Fraction(1, p)
                probability *= Fraction(1, p + 1)
            else:
                probability *= Fraction(p, p + 1)
        weights[order] += weight
        probabilities[order] += probability
    return weights, probabilities


# Preserve the exact provisional-draft fingerprint as provenance, and pin the
# promoted primary without importing or executing it.  The canonical audit is
# deliberately independent of the disposable `.scratch` tree.
repo = Path(__file__).resolve().parents[1]
draft_sha256 = "4aecf426fc577154dfb6765fa454dc0dffd0b3215bc94d7c6c5c3189d2bb95d5"
primary_path = repo / "04-computation/two_cube_euler_product_band_thm3793.py"
primary_sha256 = hashlib.sha256(primary_path.read_bytes()).hexdigest()
require(
    len(draft_sha256) == 64 and all(c in "0123456789abcdef" for c in draft_sha256),
    "audited draft fingerprint",
)
require(
    primary_sha256 == "a4bd1a16d175ae54b4fa37d0bca01ad629b5c0be69dd08bde0fea4d79e6722d1",
    "audited primary bytes",
)


# Exact Bernoulli normalization, including empty/singleton/full layers and
# several nontrivial truncations.  Brute enumeration is algorithmically
# independent of the candidate's coefficient recurrences.
all_primes = incremental_primes(101)
full_bank = [p for p in all_primes if p >= 5 and p % 3 == 2]
bank_summaries: list[tuple[int, Fraction, Fraction, Fraction]] = []
for bank_size in (0, 1, 4, 8, len(full_bank)):
    bank = full_bank[:bank_size]
    weights, probabilities = brute_layers(bank)
    euler = math.prod((Fraction(p + 1, p) for p in bank), start=Fraction(1))
    require(sum(weights, Fraction(0)) == euler, f"Euler weight bank={bank_size}")
    require(sum(probabilities, Fraction(0)) == 1, f"probability mass bank={bank_size}")
    for j, (weight, probability) in enumerate(zip(weights, probabilities)):
        require(euler * probability == weight, f"native layer bank={bank_size},j={j}")
    for cutoff in range(bank_size + 1):
        require(
            sum(weights[1 : cutoff + 1], Fraction(0))
            == euler * sum(probabilities[1 : cutoff + 1], Fraction(0)),
            f"truncated native law bank={bank_size},J={cutoff}",
        )
    mean_distribution = sum(
        (j * probabilities[j] for j in range(bank_size + 1)), Fraction(0)
    )
    second_distribution = sum(
        (j * j * probabilities[j] for j in range(bank_size + 1)), Fraction(0)
    )
    mean_terms = sum((Fraction(1, p + 1) for p in bank), Fraction(0))
    variance_terms = sum((Fraction(p, (p + 1) ** 2) for p in bank), Fraction(0))
    require(mean_distribution == mean_terms, f"mean bank={bank_size}")
    require(
        second_distribution - mean_distribution * mean_distribution == variance_terms,
        f"variance bank={bank_size}",
    )
    require(variance_terms <= mean_terms, f"variance upper bound bank={bank_size}")
    require(probabilities[0] == 1 / euler, f"zero layer bank={bank_size}")
    bank_summaries.append((bank_size, euler, mean_terms, variance_terms))

# A hostile alternative q_p=1/p does not give the claimed native law.
hostile_bank = full_bank[:4]
hostile_euler = math.prod(
    (Fraction(p + 1, p) for p in hostile_bank), start=Fraction(1)
)
wrong_zero_probability = math.prod(
    (Fraction(p - 1, p) for p in hostile_bank), start=Fraction(1)
)
require(hostile_euler * wrong_zero_probability != 1, "reject q_p=1/p normalization")


# Exact finite Euler-product bookkeeping.  Q and P_2 include p=2, whereas
# E excludes it.  The factor 2/3 is therefore forced at every finite cutoff.
factor_rows: list[tuple[int, Fraction, Fraction]] = []
for bound in (5, 17, 101):
    primes = incremental_primes(bound)
    class_one = [p for p in primes if p % 3 == 1]
    class_two = [p for p in primes if p % 3 == 2]
    p_one = math.prod((Fraction(p - 1, p) for p in class_one), start=Fraction(1))
    p_two = math.prod((Fraction(p - 1, p) for p in class_two), start=Fraction(1))
    q_two = math.prod(
        (Fraction(p * p - 1, p * p) for p in class_two), start=Fraction(1)
    )
    e_excluding_two = math.prod(
        (Fraction(p + 1, p) for p in class_two if p >= 5), start=Fraction(1)
    )
    require(e_excluding_two == Fraction(2, 3) * q_two / p_two, f"p=2 invoice {bound}")
    require(Fraction(3, 2) * e_excluding_two == q_two / p_two, f"included E {bound}")
    full_mertens = math.prod(
        (Fraction(p - 1, p) for p in primes), start=Fraction(1)
    )
    require(full_mertens == Fraction(2, 3) * p_one * p_two, f"prime 3 invoice {bound}")
    # The symmetrically truncated chi_-3 Euler product is P_2/(P_1 Q_32).
    l_truncated = math.prod(
        (
            Fraction(p, p - 1)
            if p % 3 == 1
            else Fraction(p, p + 1)
            for p in primes
            if p != 3
        ),
        start=Fraction(1),
    )
    require(l_truncated == p_two / (p_one * q_two), f"L-product quotient {bound}")
    factor_rows.append((bound, e_excluding_two, q_two))


# Independent numerical routes: an Euler product, the absolutely convergent
# Q product inserted into the derived exact constant, and the elementary
# Dirichlet series evaluation L(1,chi_-3)=pi/(3 sqrt(3)).
EULER_GAMMA = 0.577215664901532860606512090082402431
L_EXACT = math.pi / (3 * math.sqrt(3))
l_terms = 300_000
l_series = math.fsum(1 / (3 * n + 1) - 1 / (3 * n + 2) for n in range(l_terms))
require(abs(l_series - L_EXACT) < 4e-7, "L(1,chi_-3) series value")

numeric_rows: list[tuple[int, float, float, float, float]] = []
for bound in (10_000, 50_000, 200_000):
    primes = incremental_primes(bound)
    q_value = math.prod(1 - 1 / (p * p) for p in primes if p % 3 == 2)
    e_value = math.prod(1 + 1 / p for p in primes if p >= 5 and p % 3 == 2)
    c_from_q = Fraction(2, 3) * math.sqrt(
        2 * math.exp(EULER_GAMMA) * q_value / (3 * L_EXACT)
    )
    c_direct = e_value / math.sqrt(math.log(bound))
    c_critical = Fraction(2, 5) * math.sqrt(Fraction(2, 3)) * c_from_q
    require(abs(float(c_from_q) - c_direct) < 2e-4, f"constant route bound={bound}")
    # Omitting the 2/3 invoice produces exactly a 3/2 multiplicative error.
    c_wrong = math.sqrt(2 * math.exp(EULER_GAMMA) * q_value / (3 * L_EXACT))
    require(abs(c_wrong / float(c_from_q) - 1.5) < 2e-15, f"hostile p=2 bound={bound}")
    numeric_rows.append((bound, q_value, float(c_from_q), c_direct, float(c_critical)))
require(0.7854 < numeric_rows[-1][2] < 0.7859, "E_32 numerical range")
require(0.2564 < numeric_rows[-1][4] < 0.2568, "critical constant numerical range")

# Algebraic reconstruction from the two limiting relations.  If
# M=C_1 C_2=(3/2)e^-gamma and C_1/C_2=1/(LQ), then
# C_2=sqrt(M L Q), and (2/3)Q/C_2 is exactly the displayed E_32.
q_last = numeric_rows[-1][1]
mertens_pair = 1.5 * math.exp(-EULER_GAMMA)
c_two = math.sqrt(mertens_pair * L_EXACT * q_last)
c_solved = Fraction(2, 3) * q_last / c_two
c_displayed = Fraction(2, 3) * math.sqrt(
    2 * math.exp(EULER_GAMMA) * q_last / (3 * L_EXACT)
)
require(abs(float(c_solved) - float(c_displayed)) < 2e-15, "constant algebra")


# Floor-sensitive all-X ledger.  The construction is made separately for
# every sufficiently large real X; it is not a sparse subsequence argument.
cutoff_rows: list[tuple[int, int, float, float, float]] = []
for L in (10**6, 10**7, 10**8, 10**9, 10**10):
    raw = L / 2 - 0.5 * math.log(L) + L ** (2 / 3)
    J = math.floor(raw)
    theta = raw - J
    half_loglog_z = 0.5 * (L - math.log(3 * J))
    displacement = J - half_loglog_z
    ledger = L ** (2 / 3) + 0.5 * math.log(3 * J / L) - theta
    require(abs(displacement - ledger) < 2e-6, f"floor identity L={L}")
    require(0 <= theta < 1, f"floor range L={L}")
    require(displacement > 0, f"positive upper-tail gap L={L}")
    require(displacement / math.sqrt(L) > 9, f"gap exceeds sd L={L}")
    chebyshev_proxy = L / displacement**2
    require(chebyshev_proxy < 0.011, f"Chebyshev proxy L={L}")
    ratio = math.sqrt(L / (3 * J))
    require(abs(ratio - math.sqrt(Fraction(2, 3))) < 0.01, f"height ratio L={L}")
    cutoff_rows.append((L, J, displacement, chebyshev_proxy, ratio))

# Explicit moderate all-real-X controls in logarithmic coordinates.  Setting
# log Z=log X/(3J) gives Z^(3J)=X identically and Z>=11 in these samples.
all_x_rows: list[tuple[int, int, float]] = []
for L in (10, 20, 50, 100, 200):
    raw = L / 2 - 0.5 * math.log(L) + L ** (2 / 3)
    J = math.floor(raw)
    log_x = math.exp(L)
    log_z = log_x / (3 * J)
    require(J >= 1, f"positive J L={L}")
    require(log_z >= math.log(11), f"Z>=11 L={L}")
    require(abs(3 * J * log_z / log_x - 1) < 3e-16, f"exact log height L={L}")
    all_x_rows.append((L, J, math.log(log_z)))


semantic = (
    "finite normalized elementary layers are exactly Bernoulli orders with q_p=1/(p+1)",
    "Williams fixed-modulus Mertens products plus ordinary Mertens and the ordered chi_-3 Euler product give E_32",
    "the excluded prime 2 contributes the exact factor 2/3 while Q_32 includes p=2",
    "the X-dependent cutoff applies THM-3793 at every sufficiently large real X and captures Euler mass 1-o(1)",
    "the conclusion is only a lower bound for deduplicated weighted H and gives no support-count asymptotic",
)
semantic_sha256 = hashlib.sha256(repr(semantic).encode()).hexdigest()

print("TWO_CUBE_EULER_PRODUCT_BAND_INDEPENDENT_AUDIT_20260823")
print("verdict=PASS_GIVEN_WILLIAMS_FIXED_MODULUS_MERTENS_AP_IMPORT")
print(f"draft_sha256={draft_sha256}")
print(f"primary_sha256={primary_sha256}")
print(
    "full_bank="
    + repr(
        (
            len(full_bank),
            str(bank_summaries[-1][1]),
            str(bank_summaries[-1][2]),
            str(bank_summaries[-1][3]),
        )
    ).replace(" ", "")
)
print(
    "factor_rows="
    + repr(tuple((bound, str(euler), str(q_value)) for bound, euler, q_value in factor_rows)).replace(
        " ", ""
    )
)
print(f"L_series_terms={l_terms};L_series={l_series:.12f};L_exact={L_EXACT:.12f}")
print(
    "numeric_rows="
    + repr(
        tuple(
            (bound, f"{q_value:.12f}", f"{c_q:.12f}", f"{c_direct:.12f}", f"{c_h:.12f}")
            for bound, q_value, c_q, c_direct, c_h in numeric_rows
        )
    ).replace(" ", "")
)
print(
    "cutoff_rows="
    + repr(
        tuple(
            (L, J, f"{gap:.6f}", f"{tail:.6e}", f"{ratio:.12f}")
            for L, J, gap, tail, ratio in cutoff_rows
        )
    ).replace(" ", "")
)
print("all_X_mode=pointwise_every_sufficiently_large_real_X;Z_power_3J_equals_X")
print("scope=deduplicated_weighted_H_lower_bound_only;no_support_count_asymptotic")
print(f"semantic_sha256={semantic_sha256}")
print(f"gates={GATES}")
print("ALL CHECKS PASSED")
