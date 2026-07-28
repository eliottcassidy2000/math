#!/usr/bin/env python3
"""Exact central-tooth scout for odometer-twisted dilation handoffs.

THM-2684 proves three-event nilpotence for the fixed handoff D(x)={13x} on
the complete THM-2584 three-tooth rail envelope.  THM-2657 classifies the
physical carry/root translations k/13^6 with k a 13-unit.  This scout tests
the affine composition

    T_k(x) = {13*x + k/13^6}.

It first records why the minimal alternating lifts k=-14,+14 are only a
cosmetic orbit-clock cycle: their intrinsic stored edges are diagonal and
their owner/shallow interface does not glue.  The larger lawful lifts

    k=+/-(13^5+1)

instead exchange 1/2 +/- k/(14*13^6).  Their intrinsic stored clock edges are
4->3 and 3->4, and each owner is exactly the following shallow clock.  Positive
finite-horizon neighbourhoods exist at every prescribed length.  The
translation quotient labels are again 11 and 2.

No labelled rail, present factor, delayed word, unit, semantic endpoint, or
row consequence is tested here.  In particular, this is a lawful clock/state
carrier, not yet a positive full-packet transition.
"""

from fractions import Fraction


P = 13
Q = 7
R = P**6
S = P**5
CENTRE = Fraction(1, 2)
TOOTH_RADIUS = Fraction(1, 28)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def circle(value):
    return value % 1


def handoff(value, lift):
    return circle(P * value + Fraction(lift, R))


def dilate(value):
    return circle(P * value)


def clock(value):
    return int((Q * circle(value) + Fraction(1, 2)) // 1) % Q


def shallow(value):
    return clock(dilate(value))


def owner(value):
    return clock(dilate(dilate(value)))


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


SMALL_PLUS = CENTRE + Fraction(1, R)
SMALL_MINUS = CENTRE - Fraction(1, R)
SMALL_LIFTS = (-14, 14)
require(handoff(SMALL_PLUS, -14) == SMALL_MINUS
        and handoff(SMALL_MINUS, 14) == SMALL_PLUS,
        "minimal lifts lost their orbit-clock two-cycle")
require((shallow(SMALL_PLUS), owner(SMALL_PLUS)) == (4, 4)
        and (shallow(SMALL_MINUS), owner(SMALL_MINUS)) == (3, 3),
        "minimal-lift stored-edge hostile changed")
require(owner(SMALL_PLUS) != shallow(SMALL_MINUS)
        and owner(SMALL_MINUS) != shallow(SMALL_PLUS),
        "minimal lifts unexpectedly acquired the existing clock interface")

LIFT_MAGNITUDE = S + 1
LIFTS = (-LIFT_MAGNITUDE, LIFT_MAGNITUDE)
AMPLITUDE = Fraction(LIFT_MAGNITUDE, (P + 1) * R)
PLUS = CENTRE + AMPLITUDE
MINUS = CENTRE - AMPLITUDE
require(handoff(PLUS, LIFTS[0]) == MINUS,
        "negative odometer lift lost the lawful central two-cycle")
require(handoff(MINUS, LIFTS[1]) == PLUS,
        "positive odometer lift lost the lawful central two-cycle")
require(all(in_central_tooth(value) for value in (PLUS, MINUS)),
        "lawful central two-cycle left the THM-2684 middle tooth")
require((shallow(PLUS), owner(PLUS)) == (4, 3)
        and (shallow(MINUS), owner(MINUS)) == (3, 4),
        "lawful central two-cycle lost its nonconstant stored edges")
require(owner(PLUS) == shallow(MINUS)
        and owner(MINUS) == shallow(PLUS),
        "lawful central two-cycle lost owner-to-next-shallow gluing")
require((carry(PLUS), carry(MINUS)) == (7, 5),
        "lawful central two-cycle predecessor carries changed")
require((future_digit(PLUS), future_digit(MINUS)) == (6, 6),
        "lawful central two-cycle future digits changed")
require(
    carry(handoff(PLUS, LIFTS[0]))
    == (future_digit(PLUS) + LIFTS[0]) % P
    and carry(handoff(MINUS, LIFTS[1]))
    == (future_digit(MINUS) + LIFTS[1]) % P,
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
for lift in SMALL_LIFTS:
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
# P^n e.  Radius 1/(3 R P^H) keeps every state and its two intrinsic clock
# images away from their boundaries and inside the central tooth.  This gives
# a positive finite-horizon interval for every H on which stored edges are
# nonconstant and owner-to-next-shallow gluing holds.  The exact two-cycle is
# repelling, so there is no common positive interval for all H.
HORIZONS = (3, 7, 12)
horizon_rows = []
for horizon in HORIZONS:
    radius = Fraction(1, 3 * R * P**horizon)
    for sign in (-1, 1):
        value = PLUS + sign * radius
        require(in_central_tooth(value), "horizon endpoint left central tooth")
        for step in range(horizon):
            next_value = handoff(value, LIFTS[step % 2])
            require(in_central_tooth(next_value),
                    "alternating finite-horizon orbit left central tooth")
            require(shallow(value) != owner(value),
                    "alternating finite-horizon stored edge became constant")
            require(owner(value) == shallow(next_value),
                    "alternating finite-horizon interface stopped gluing")
            value = next_value
    interval_length = 2 * radius
    horizon_rows.append((horizon, radius, interval_length))

print("LRC14 odometer-twisted central-clock escape")
print("status=VERIFIED-EXACT SCOUT; envelope/intrinsic-clock/state level only")
print(f"p={P} q={Q} R={R} central_tooth=[13/28,15/28)")
print(f"minimal_cycle_control=({SMALL_PLUS},{SMALL_MINUS}):stored_edges=((4,4),(3,3)):existing_interface=False")
print(f"lawful_central_cycle=({PLUS},{MINUS}) amplitude={AMPLITUDE}")
print(f"alternating_lifts={LIFTS} quotient_root_steps={(11, 2)}")
print(f"cycle_stored_edges={((shallow(PLUS), owner(PLUS)), (shallow(MINUS), owner(MINUS)))} interface_glues=True")
print(f"cycle_carries={(carry(PLUS), carry(MINUS))} future_digits={(future_digit(PLUS), future_digit(MINUS))}")
for lift, z1, z2, z3, length in fixed_rows:
    print(
        f"fixed_lift={lift}:zero_thresholds=({z1},{z2},{z3}):"
        f"first_flip_length={length}:second_consecutive_flip=EMPTY"
    )
for horizon, radius, length in horizon_rows:
    print(
        f"alternating_horizon={horizon}:initial_radius={radius}:"
        f"positive_interval_length={length}:stored_edges_nonconstant=True:"
        f"interfaces_glue=True"
    )
print("mechanism=minimal lifts are cosmetic; lifts +/-(13^5+1) align intrinsic owner/shallow clocks on a repelling two-cycle")
print("not_tested=labelled rail fibre; present factors; delayed words; primitive units; semantic endpoint; row exclusion; LRC14")
print("ALL CHECKS PASSED")
