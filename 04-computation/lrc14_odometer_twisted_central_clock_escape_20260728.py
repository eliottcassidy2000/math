#!/usr/bin/env python3
"""Exact central-tooth scout for odometer-twisted dilation handoffs.

THM-2684 proves three-event nilpotence for the fixed handoff D(x)={13x} on
the complete THM-2584 three-tooth rail envelope.  THM-2657 classifies the
physical carry/root translations k/13^6 with k a 13-unit.  This scout tests
the affine composition

    T_k(x) = {13*x + k/13^6}.

It proves two sharply scoped facts at the envelope/nearest-clock/state level.
A fixed nonzero twist creates one positive clock flip but cannot create two
consecutive flips.  Alternating k=-14,+14 instead has the exact central
two-cycle 1/2+1/13^6 <-> 1/2-1/13^6 and positive finite-horizon
neighbourhoods of every prescribed length.  The translation quotient labels
are delta=2k mod 13, namely 11 and 2.

No labelled rail, present factor, delayed word, unit, semantic endpoint, or
row consequence is tested here.
"""

from fractions import Fraction


P = 13
Q = 7
R = P**6
CENTRE = Fraction(1, 2)
TOOTH_RADIUS = Fraction(1, 28)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def circle(value):
    return value % 1


def handoff(value, lift):
    return circle(P * value + Fraction(lift, R))


def clock(value):
    return int((Q * circle(value) + Fraction(1, 2)) // 1) % Q


def carry(value):
    return int((R * circle(value)) // 1) % P


def future_digit(value):
    fractional = circle(R * circle(value))
    return int((P * fractional) // 1)


def in_central_tooth(value):
    return CENTRE - TOOTH_RADIUS <= value < CENTRE + TOOTH_RADIUS


def fixed_twist_threshold(lift, iterate):
    """Initial displacement u for which the iterate displacement is zero."""
    tau = Fraction(lift, R)
    return -tau * Fraction(P**iterate - 1, P**iterate * (P - 1))


PLUS = CENTRE + Fraction(1, R)
MINUS = CENTRE - Fraction(1, R)
LIFTS = (-14, 14)

require(handoff(PLUS, LIFTS[0]) == MINUS,
        "negative odometer lift lost the central two-cycle")
require(handoff(MINUS, LIFTS[1]) == PLUS,
        "positive odometer lift lost the central two-cycle")
require(all(in_central_tooth(value) for value in (PLUS, MINUS)),
        "central two-cycle left the THM-2684 middle tooth")
require((clock(PLUS), clock(MINUS)) == (4, 3),
        "central two-cycle lost its nonconstant nearest clocks")
require((carry(PLUS), carry(MINUS)) == (7, 5),
        "central two-cycle predecessor carries changed")
require((future_digit(PLUS), future_digit(MINUS)) == (6, 6),
        "central two-cycle future digits changed")
require(
    carry(handoff(PLUS, -14)) == (future_digit(PLUS) - 14) % P
    and carry(handoff(MINUS, 14)) == (future_digit(MINUS) + 14) % P,
    "affine handoff digit covariance failed",
)
require(tuple((2 * lift) % P for lift in LIFTS) == (11, 2),
        "THM-2657 quotient labels changed")

# A fixed twist has iterate displacements
# u_n=P^n(u+tau/(P-1))-tau/(P-1).  Its zero thresholds are monotone.
# Hence an initial u between the first two thresholds makes u_1,u_2 have
# opposite signs, but u_2,u_3 have the same sign.  The first flip interval
# has exact length |tau|/P^2.
fixed_rows = []
for lift in LIFTS:
    z1 = fixed_twist_threshold(lift, 1)
    z2 = fixed_twist_threshold(lift, 2)
    z3 = fixed_twist_threshold(lift, 3)
    if lift > 0:
        require(z3 < z2 < z1,
                "positive fixed-twist zero thresholds lost monotonicity")
    else:
        require(z1 < z2 < z3,
                "negative fixed-twist zero thresholds lost monotonicity")
    ordered = sorted((z1, z2))
    flip_length = ordered[1] - ordered[0]
    require(flip_length == Fraction(abs(lift), R * P**2),
            "fixed-twist first-flip length changed")
    midpoint = sum(ordered, Fraction(0)) / 2
    u1 = P * midpoint + Fraction(lift, R)
    u2 = P * u1 + Fraction(lift, R)
    u3 = P * u2 + Fraction(lift, R)
    require(u1 * u2 < 0 and u2 * u3 > 0,
            "fixed twist unexpectedly sustained two clock flips")
    require(all(abs(value) < TOOTH_RADIUS for value in (midpoint, u1, u2, u3)),
            "fixed-twist hostile left the central tooth")
    fixed_rows.append((lift, z1, z2, z3, flip_length))

# For the alternating lifts, perturbing PLUS by e perturbs the nth state by
# P^n e.  Radius 1/(3 R P^H) therefore keeps every state through H at least
# 2/(3R) from the clock boundary and inside the central tooth.  This supplies
# a positive finite-horizon interval for every H, although the exact two-cycle
# itself is repelling and there is no common positive interval for all H.
HORIZONS = (3, 7, 12)
horizon_rows = []
for horizon in HORIZONS:
    radius = Fraction(1, 3 * R * P**horizon)
    for sign in (-1, 1):
        value = PLUS + sign * radius
        clocks = [clock(value)]
        require(in_central_tooth(value), "horizon endpoint left central tooth")
        for step in range(horizon):
            value = handoff(value, LIFTS[step % 2])
            require(in_central_tooth(value),
                    "alternating finite-horizon orbit left central tooth")
            clocks.append(clock(value))
        require(all(left != right for left, right in zip(clocks, clocks[1:])),
                "alternating finite-horizon clocks stopped alternating")
    interval_length = 2 * radius
    horizon_rows.append((horizon, radius, interval_length))

print("LRC14 odometer-twisted central-clock escape")
print("status=VERIFIED-EXACT SCOUT; envelope/clock/state level only")
print(f"p={P} q={Q} R={R} central_tooth=[13/28,15/28)")
print(f"central_cycle=({PLUS},{MINUS})")
print(f"alternating_lifts={LIFTS} quotient_root_steps={(11, 2)}")
print(f"cycle_clocks={(clock(PLUS), clock(MINUS))} carries={(carry(PLUS), carry(MINUS))} future_digits={(future_digit(PLUS), future_digit(MINUS))}")
for lift, z1, z2, z3, length in fixed_rows:
    print(
        f"fixed_lift={lift}:zero_thresholds=({z1},{z2},{z3}):"
        f"first_flip_length={length}:second_consecutive_flip=EMPTY"
    )
for horizon, radius, length in horizon_rows:
    print(
        f"alternating_horizon={horizon}:initial_radius={radius}:"
        f"positive_interval_length={length}:all_clocks_nonconstant=True"
    )
print("mechanism=fixed twist crosses the central clock wall once; alternating odometer digits create a repelling two-cycle")
print("not_tested=labelled rails; present factors; delayed words; primitive units; semantic transition; row exclusion; LRC14")
print("ALL CHECKS PASSED")
