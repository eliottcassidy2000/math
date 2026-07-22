#!/usr/bin/env python3
"""Exact arithmetic referee for the strengthened part of THM-2051.

The analytic Fejer/BV estimate is proved in THM-2051.  This referee checks all
finite rational arithmetic used by its two versions, and the finite core of
the sharp THM-965 pair-covariance floor.

Tournament Analysis is not forced here.  The faithful combinatorial object is
a height-filtered signed relation hypergraph: a hyperedge retains its support
and every integer coefficient.  Turning it into a binary tournament destroys
the coefficient vector that character orthogonality tests, and there is no
canonical antisymmetric switch or Hamiltonian path preserving the theorem.
"""

from fractions import Fraction as F
from math import comb, gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def fold14(r):
    r %= 14
    return r * (14 - r)


def pair_delta_reduced(a, b):
    """Centered pair moment from THM-965 after gcd reduction."""
    require(0 < a < b and gcd(a, b) == 1, "pair must be positive, ordered, coprime")
    return F(fold14(a + b) - fold14(b - a), 196 * a * b)


def centered_coefficient(s):
    return (-1) ** s * sum(
        (F(-1, 7) ** j) * comb(13 - s, j) for j in range(6 - s)
    )


expected_coefficients = {
    2: F(24, 343),
    3: F(-24, 49),
    4: F(-2, 7),
    5: F(-1),
}
actual_coefficients = {s: centered_coefficient(s) for s in range(2, 6)}
require(actual_coefficients == expected_coefficients, "centered B5 coefficients failed")

base = sum(F((-1) ** k * comb(13, k), 7**k) for k in range(6))
require(base == F(2052, 16807), "equilibrium B5 value failed")

# If ab>=27, THM-965 and 0<=fold14<=49 give delta>=-1/(4ab)>=-1/108,
# which is strictly above -6/637.  The remaining proof core is finite.
pair_floor = F(-6, 637)
require(F(-1, 108) > pair_floor, "large-product pair tail comparison failed")

finite_core = [
    (a, b, pair_delta_reduced(a, b))
    for a in range(1, 27)
    for b in range(a + 1, 27)
    if a * b <= 26 and gcd(a, b) == 1
]
require(len(finite_core) == 36, "unexpected finite pair core size")
require(min(delta for _, _, delta in finite_core) == pair_floor, "pair floor failed")
require(
    [(a, b) for a, b, delta in finite_core if delta == pair_floor] == [(1, 13)],
    "pair-floor equality classifier failed",
)

group_minima = []
for a in sorted({a for a, _, _ in finite_core}):
    delta, b = min((d, b) for aa, b, d in finite_core if aa == a)
    group_minima.append((a, b, delta))
require(
    group_minima
    == [(1, 13, F(-6, 637)), (2, 11, F(-4, 539)),
        (3, 8, F(-1, 392)), (4, 5, F(2, 245))],
    "finite-core grouped minima failed",
)

def telescope_cost(s):
    return comb(13, s) * abs(actual_coefficients[s]) * s * F(6, 7) ** (s - 1)


costs = {s: telescope_cost(s) for s in range(2, 6)}
K_all = sum(costs.values(), F(0))
require(K_all == F(1477008, 343), "all-support telescope constant failed")

N_all = 2**20 + 1
error_all = K_all * F(43, 2 * N_all)
require(error_all == F(31755672, 359661911), "all-support error failed")
margin_all = base - error_all
require(margin_all == F(595652076, 17623433639) and margin_all > 0,
        "all-support margin failed")

# Stronger theorem: pay every pair exactly, then Fejer-annihilate only support 3..5.
pair_contribution_floor = expected_coefficients[2] * comb(13, 2) * pair_floor
require(pair_contribution_floor == F(-864, 16807), "aggregate pair floor failed")
remaining_base = base + pair_contribution_floor
require(remaining_base == F(1188, 16807), "post-pair base failed")

K_high = sum((costs[s] for s in range(3, 6)), F(0))
require(K_high == F(10316592, 2401), "higher-support telescope constant failed")
N_high = 2**21 + 1
error_high = K_high * F(45, 2 * N_high)
require(error_high == F(25791480, 559473817), "higher-support error failed")
margin_high = remaining_base - error_high
require(margin_high == F(96283836, 3916316719) and margin_high > 0,
        "higher-support margin failed")

# Diagnostic beyond the proof cutoff: the same sharp floor persists in a broad box.
box_min = min(
    (pair_delta_reduced(a, b), a, b)
    for a in range(1, 501)
    for b in range(a + 1, 501)
    if gcd(a, b) == 1
)
require(box_min == (pair_floor, 1, 13), "broad pair-floor diagnostic failed")

print("THM-2051 FEJER--BV HIGHER-CIRCUIT REFEREE")
print(f"centered coefficients={actual_coefficients}")
print(f"equilibrium B5={base}")
print(f"pair finite core size={len(finite_core)} grouped minima={group_minima}")
print(f"pair floor={pair_floor} equality reduced ratio=1:13")
print(f"all-support K={K_all} error<{error_all} margin>{margin_all}")
print(f"pair contribution floor={pair_contribution_floor} remaining base={remaining_base}")
print(f"higher-support K={K_high} error<{error_high} margin>{margin_high}")
print("TOURNAMENT ANALYSIS=not applicable: signed coefficient hyperedges are indispensable")
print("RESULT=PASS")
