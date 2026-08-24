#!/usr/bin/env python3
"""Exact controls for THM-3995's body-maximum parity-hole tariff."""

from fractions import Fraction
from math import gcd

import sympy as sp


def require(label, condition):
    if not condition:
        raise RuntimeError(f"FAIL  {label}")
    print(f"PASS  {label}")


# THM-3878's exact obstruction for the scale-two (p,q)=(1,9) row.
boundaries = (
    Fraction(2, 21),
    Fraction(8, 63),
    Fraction(55, 63),
    Fraction(19, 21),
)
component_lengths = (
    boundaries[1] - boundaries[0],
    boundaries[3] - boundaries[2],
)
require(
    "two obstruction components each have length 2/63",
    component_lengths == (Fraction(2, 63),) * 2,
)
require("obstruction measure is 4/63", sum(component_lengths) == Fraction(4, 63))


# Entering a safe owner interval uses 14a+1; leaving uses 14a-1.  After
# clearing denominators at the four oriented target boundaries, the event
# numerators are the following.  For odd t they are always odd, independently
# of a,u,M, so equality with the target boundary is impossible.
def parity_rows(t, a, u, M):
    return (
        3 * t * (14 * a + 1) - 4 * u - 42 * u * M,
        9 * t * (14 * a - 1) - 16 * u - 126 * u * M,
        9 * t * (14 * a + 1) - 110 * u - 126 * u * M,
        3 * t * (14 * a - 1) - 38 * u - 42 * u * M,
    )


t_sym, a_sym, u_sym, M_sym = sp.symbols("t a u M", integer=True)
enter = t_sym * (14 * a_sym + 1) / (14 * u_sym)
leave = t_sym * (14 * a_sym - 1) / (14 * u_sym)
symbolic_rows = (
    sp.expand(42 * u_sym * (enter - sp.Rational(2, 21) - M_sym)),
    sp.expand(126 * u_sym * (leave - sp.Rational(8, 63) - M_sym)),
    sp.expand(126 * u_sym * (enter - sp.Rational(55, 63) - M_sym)),
    sp.expand(42 * u_sym * (leave - sp.Rational(19, 21) - M_sym)),
)
require(
    "all four oriented wall differences reconstruct the stated numerators",
    all(
        sp.expand(lhs - rhs) == 0
        for lhs, rhs in zip(
            symbolic_rows, parity_rows(t_sym, a_sym, u_sym, M_sym)
        )
    ),
)

parity_checks = 0
for t_residue in (1, 3):
    for a_residue in (0, 1):
        for u_residue in (0, 1):
            for M_residue in (0, 1):
                values = parity_rows(t_residue, a_residue, u_residue, M_residue)
                if not all(value % 2 == 1 for value in values):
                    raise RuntimeError(
                        (t_residue, a_residue, u_residue, M_residue, values)
                    )
                parity_checks += 4
print(f"PASS  all four cleared event numerators are odd ({parity_checks} residue checks)")

# Odd t, not odd U, is load-bearing.  If even t were admitted, oriented walls
# could coincide with target boundaries and a parity hole could disappear.
require(
    "even-t hostile hits the first and fourth target boundaries",
    parity_rows(12, 0, 9, 0)[0] == 0
    and parity_rows(12, 0, 9, -1)[3] == 0,
)


# The exact U-aware clockwise sweep.  A modulo residue is the directed wall
# distance numerator after choosing the correct circular lift M.  The four
# directions are enter-after, exit-before, enter-after, exit-before.
event_checks = 0
eventwise_equalities = [False, False, False, False]
for U_value in range(11, 51):
    for t_value in range(U_value, 4 * U_value // 3 + 1):
        if t_value % 2 == 0 or 3 * t_value >= 4 * U_value:
            continue

        holes = (
            Fraction(1, 42 * U_value),
            Fraction(1, 126 * U_value),
            Fraction(1, 126 * U_value),
            Fraction(1, 42 * U_value),
        )
        if not (
            boundaries[0] + holes[0] < boundaries[1] - holes[1]
            and boundaries[2] + holes[2] < boundaries[3] - holes[3]
        ):
            raise RuntimeError((U_value, t_value, "overlapping trims"))

        for u_value in range(1, U_value + 1):
            for a_value in range(u_value):
                raw = parity_rows(t_value, a_value, u_value, 0)
                moduli = (42 * u_value, 126 * u_value, 126 * u_value, 42 * u_value)
                directed = (
                    raw[0] % moduli[0],
                    (-raw[1]) % moduli[1],
                    raw[2] % moduli[2],
                    (-raw[3]) % moduli[3],
                )
                for j, residue in enumerate(directed):
                    if residue == 0 or residue % 2 == 0:
                        raise RuntimeError(
                            (U_value, t_value, u_value, a_value, j, residue)
                        )
                    # residue/(d*u) >= 1/(d*U) iff residue*U >= u.
                    if residue * U_value < u_value:
                        raise RuntimeError(
                            (U_value, t_value, u_value, a_value, j, residue)
                        )
                    if residue == 1 and u_value == U_value:
                        eventwise_equalities[j] = True
                    event_checks += 1
require(
    f"U-aware clockwise wall sweep ({event_checks} oriented events)",
    eventwise_equalities == [True, True, True, True],
)

# Explicit equality controls for the four individual wall constants.  These
# do not assert simultaneous equality for one body or an actual LRC failure.
sharp_rows = (
    parity_rows(45, 1, 44, 1)[0],
    parity_rows(41, 8, 40, 8)[1],
    parity_rows(41, 32, 40, 32)[2],
    parity_rows(45, 43, 44, 43)[3],
)
require("all four wall constants are individually attained", sharp_rows == (1, -1, 1, -1))

# Simultaneous masking is real: owner 1 enters and owner 13 exits at y=1/14.
# It cannot create an earlier positive cell; a first positive jump still has
# at least one entering event.
masked_enter = Fraction(14 * 0 + 1, 14 * 1)
masked_exit = Fraction(14 * 1 - 1, 14 * 13)
require("oppositely oriented owner walls can mask at one point", masked_enter == masked_exit)


# All four equalities cannot be owned by u=U: rows 1/4 would require
# 3t-4U=1 mod 42, while rows 2/3 require 9t+16U=1 mod 126.  Multiplying the
# first by 3 forces 28U=-2 mod 126, impossible modulo gcd(28,126)=14.
require("same-maximum four-wall equality congruence is impossible", gcd(28, 126) == 14 and (-2) % 14 != 0)


def target_event(row_index, t_value, u_value, target):
    """Return one (a,M) attaining a target numerator, or None."""

    modulus = (42, 126, 126, 42)[row_index] * u_value
    for a_value in range(u_value):
        raw = parity_rows(t_value, a_value, u_value, 0)[row_index]
        if (raw - target) % modulus == 0:
            return a_value, (raw - target) // modulus
    return None


# Nevertheless the aggregate 4/(63U) coefficient is asymptotically sharp at
# the oriented-event arithmetic layer.  This still does not realize equality
# of the support cap or produce an LRC failure.
family_checks = 0
for n_value in range(21):
    if n_value % 13 == 11:
        continue
    U_value = 126 * n_value + 122
    t_value = U_value + 13
    owners = (U_value, U_value - 1, U_value - 1, U_value)
    targets = (1, -1, 1, -1)
    if not (t_value % 2 == 1 and U_value <= t_value and 3 * t_value < 4 * U_value):
        raise RuntimeError((n_value, U_value, t_value))
    if gcd(t_value, U_value) != 1 or gcd(t_value, U_value - 1) != 1:
        raise RuntimeError((n_value, U_value, t_value, "gcd"))
    if any(
        target_event(j, t_value, owners[j], targets[j]) is None
        for j in range(4)
    ):
        raise RuntimeError((n_value, U_value, t_value, "target event"))
    ratio = Fraction(U_value, 8) * (
        Fraction(3, owners[0])
        + Fraction(1, owners[1])
        + Fraction(1, owners[2])
        + Fraction(3, owners[3])
    )
    if ratio != 1 + Fraction(1, 4 * (U_value - 1)):
        raise RuntimeError((n_value, ratio))
    family_checks += 1
require(
    f"asymptotically sharp U/U-1 event family ({family_checks} exact controls)",
    family_checks == 20,
)


U = sp.symbols("U", integer=True, positive=True)
t = sp.symbols("t", integer=True, positive=True)
s_U = sp.Rational(4, 63) - sp.Rational(4, 63) / U
s_t = sp.Rational(4, 63) - sp.Rational(4, 63) / t
require("body-maximum support cap simplifies", sp.factor(s_U - 4 * (U - 1) / (63 * U)) == 0)
require(
    "each closed trimmed core has length 2(U-1)/(63U)",
    sp.factor(
        (sp.Rational(8, 63) - 1 / (126 * U))
        - (sp.Rational(2, 21) + 1 / (42 * U))
        - 2 * (U - 1) / (63 * U)
    )
    == 0,
)
require(
    "old-minus-new support cap is 4(t-U)/(63tU)",
    sp.factor(s_t - s_U - 4 * (t - U) / (63 * t * U)) == 0,
)
c_U = sp.factor((1 - s_U) / s_U)
require(
    "body-maximum variance coefficient simplifies",
    sp.factor(c_U - (59 * U + 4) / (4 * (U - 1))) == 0,
)


# Independent finite dynamic-programming controls for the integer optimization:
# among q equally weighted cells, total mass M, and at most B positive cells,
# the minimum sum of squares is obtained at adjacent integers.  This includes
# k=0, where only M of the B allowed cells need be positive.
integer_checks = 0
for q in range(2, 10):
    max_total = 4 * q
    dp = {(0, 0): 0}
    for _ in range(q):
        nxt = {}
        for (positive, total), square_sum in dp.items():
            for value in range(max_total - total + 1):
                key = (positive + (value > 0), total + value)
                candidate = square_sum + value * value
                if candidate < nxt.get(key, 10**9):
                    nxt[key] = candidate
        dp = nxt
    for B in range(1, q + 1):
        for M in range(0, 4 * B + 1):
            actual = min(
                square_sum
                for (positive, total), square_sum in dp.items()
                if total == M and positive <= B
            )
            k, remainder = divmod(M, B)
            predicted = (B - remainder) * k * k + remainder * (k + 1) ** 2
            if actual != predicted:
                raise RuntimeError((q, B, M, actual, predicted))
            integer_checks += 1
print(f"PASS  integer support-envelope dynamic controls ({integer_checks} cases)")

plateau_checks = 0
for smaller_cap in range(1, 8):
    for larger_cap in range(smaller_cap + 1, 10):
        for total_mass in range(smaller_cap + 1):
            small_k, small_r = divmod(total_mass, smaller_cap)
            large_k, large_r = divmod(total_mass, larger_cap)
            small_min = (smaller_cap - small_r) * small_k**2 + small_r * (small_k + 1) ** 2
            large_min = (larger_cap - large_r) * large_k**2 + large_r * (large_k + 1) ** 2
            if small_min != total_mass or large_min != total_mass:
                raise RuntimeError((smaller_cap, larger_cap, total_mass))
            plateau_checks += 1
require(f"k=0 support-cap plateau ({plateau_checks} controls)", plateau_checks == 147)

small_cap, large_cap, total_mass = 5, 7, 6
small_k, small_r = divmod(total_mass, small_cap)
large_k, large_r = divmod(total_mass, large_cap)
small_min = (small_cap - small_r) * small_k**2 + small_r * (small_k + 1) ** 2
large_min = (large_cap - large_r) * large_k**2 + large_r * (large_k + 1) ** 2
require("support shrink is strict immediately beyond the smaller-cap plateau", small_min > large_min)

# A hostile in which the integer correction is strictly positive.
B_hostile, M_hostile = 5, 7
k_hostile, rem_hostile = divmod(M_hostile, B_hostile)
require("integer correction hostile is positive", 0 < rem_hostile < B_hostile and k_hostile == 1)


# Exact discrepancy and BV algebra.  Here m=t*mu is the mean of N_t and
# disc_t(G)=Var(N_t)/t^2.
mu, theta = sp.symbols("mu theta", nonnegative=True)
variance_tariff = (t * mu) ** 2 * c_U + s_U * theta * (1 - theta)
disc_tariff = mu**2 * c_U + s_U * theta * (1 - theta) / t**2
require("variance tariff divides to the discrepancy tariff", sp.factor(variance_tariff / t**2 - disc_tariff) == 0)
require(
    "BV gate coefficient simplifies",
    sp.factor(1 / (3 * c_U) - 4 * (U - 1) / (3 * (59 * U + 4))) == 0,
)

print("RESULT obstruction_measure=4/63")
print("RESULT parity_hole_measure=4/(63U)")
print("RESULT support_cap=4(U-1)/(63U)")
print("RESULT event_spacing_sharpness=asymptotic_not_global_failure_equality")
print("RESULT k0_plateau=exact_envelope_unchanged_for_m<=s_U")
print(
    "RESULT sufficient_disc_gate=disc_t(G)<mu(G)^2*(59U+4)/(4(U-1))"
    "+(s_U/t^2)*theta*(1-theta)"
)
print("RESULT sufficient_BV_gate=t*mu(G)/r>sqrt(4(U-1)/(3(59U+4)))")
print("ALL EXACT CHECKS PASSED")
