#!/usr/bin/env python3
"""Exact phase-lift and slope-14 controls for the THM-2693 delayed word.

This is a bounded sidecar, not an LRC(14) row exclusion.  It checks three
distinct objects which must not be conflated:

1. the minimal carrier ``D_(13^3) intersect D_14^c`` under the one-level
   deeper lift ``F_s(y)={13y+s/13}``;
2. the full raw guard-sector union ``W`` under phase-labelled ``F_s`` words;
3. the unchanged full ``W`` under the heterogeneous slope ``y -> {14y}``.

The first object has twelve exact interior fixed points.  The second remains
uniformly nilpotent at depth four because the unchanged high-speed safety is
``338`` times the target phase.  The third crosses the sharp no-wrap boundary
and survives through fourteen states, but still dies at state fifteen.
"""

from fractions import Fraction

import lrc14_cross_time_target_future_diagonal_thm2616 as cross
import lrc14_odometer_alternating_lift_labelled_tail_scout_20260728 as scout
import lrc14_successor_halfcell_carry_no_go_thm2623 as prior


P = 13
R = P**6
C2 = P**3
C3 = 2 * P**5
RHO = Fraction(1, 14)


EXPECTED_PHASE_DEPTH3 = {
    "minimum": ((3, 0), 6_732, Fraction(18_513, 439_239_619)),
    "maximum": ((1, 8), 8_294, Fraction(3_509, 67_575_326)),
}

EXPECTED_SLOPE_14 = (
    (1, 47_484, Fraction(604_725_613_249, 11_455_265_301_480)),
    (2, 35_220, Fraction(7_074_988_411, 2_545_614_511_440)),
    (3, 21_206, Fraction(1_419_717_323, 11_879_534_386_720)),
    (4, 13_148, Fraction(1_006_391, 188_992_592_516)),
    (5, 8_400, Fraction(2_567_587, 10_583_585_180_896)),
    (6, 5_552, Fraction(32_003, 2_795_664_010_048)),
    (7, 3_596, Fraction(10_337, 19_569_648_070_336)),
    (8, 2_602, Fraction(15_007, 547_950_145_969_408)),
    (9, 1_794, Fraction(81, 59_932_047_215_404)),
    (10, 1_080, Fraction(3_119, 53_699_114_305_001_984)),
    (11, 696, Fraction(4_033, 1_503_575_200_540_055_552)),
    (12, 324, Fraction(1_867, 21_050_052_807_560_777_728)),
    (13, 106, Fraction(603, 294_700_739_305_850_888_192)),
    (14, 20, Fraction(109, 4_125_810_350_281_912_434_688)),
    (15, 0, Fraction(0)),
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def circle(x):
    return x % 1


def centered(x):
    value = circle(x)
    return value - 1 if value >= Fraction(1, 2) else value


def danger(x):
    value = centered(x)
    return -RHO <= value < RHO


def phase_lift(s, y):
    return circle(P * y + Fraction(s, P))


def shifted_base(intervals, s, period):
    """Return ``{u : u+s/13 belongs to intervals}`` on an integer grid."""
    require(period % P == 0, "base grid does not resolve the phase shift")
    offset = s * (period // P)
    out = []
    for left, right in intervals:
        start = (left - offset) % period
        end = start + (right - left)
        if end <= period:
            out.append((start, end))
        else:
            out.extend(((start, period), (0, end - period)))
    return scout.merge_intervals(out)


def phase_depth_three_atlas(base, period):
    """Enumerate all 13^2 two-edge words on one common exact grid."""
    targets = [shifted_base(base, s, period) for s in range(P)]
    scale = P**2
    state_zero = [(scale * left, scale * right) for left, right in base]
    rows = {}
    for first in range(P):
        state_one_target = [
            (P * left, P * right) for left, right in targets[first]
        ]
        level_two = scout.periodic_integer_restrict(
            state_zero, state_one_target, P * period
        )
        for second in range(P):
            level_three = scout.periodic_integer_restrict(
                level_two, targets[second], period
            )
            mass = Fraction(
                sum(right - left for left, right in level_three),
                scale * period,
            )
            rows[first, second] = (len(level_three), mass)
    return rows, targets


def phase_word_row(base, targets, period, word):
    sets = [base] + [targets[s] for s in word]
    intervals, grid = scout.integer_itinerary(sets)
    mass = Fraction(sum(right - left for left, right in intervals), grid)
    return len(intervals), mass


def slope_rows(full_word, slope, depth):
    current = full_word
    rows = []
    for state_count in range(1, depth + 1):
        if state_count > 1:
            current = scout.fraction_pullback_intersection(
                current, full_word, slope ** (state_count - 1)
            )
        rows.append(
            (state_count, len(current), scout.interval_mass(current))
        )
        if not current:
            break
    return tuple(rows)


def main():
    # The one-level-deeper lift on x induces the asserted affine tail map.
    # For each nonzero phase it has an exact physical fixed lift whose tail is
    # y=s/13.  That point is interior to the minimal carrier and is rejected
    # by the full carrier only at the high-speed safety factor.
    inverse_twelve = pow(12, -1, R)
    fixed_lifts = []
    margins = []
    for s in range(1, P):
        y = Fraction(s, P)
        require(phase_lift(s, y) == y, "tail fixed point failed")
        require(danger(C2 * y), "fixed point left the target tooth")
        require(not danger(14 * y), "fixed point lost speed-14 safety")
        require(danger(C3 * y), "full-word high-speed hostile disappeared")

        odometer = (-s * inverse_twelve) % R
        x = (Fraction(odometer) + y) / R
        lifted = circle(P * x + Fraction(s, P * R))
        require(lifted == x, "physical C_13^7 fixed lift failed")
        require(circle(R * x) == y, "fixed lift has the wrong delayed tail")
        fixed_lifts.append((s, x, y))
        margins.append(abs(centered(14 * y)) - RHO)

    require(min(margins) == Fraction(1, 182),
            "nearest nonzero phase margin changed")

    # Build the exact raw guard-sector union from the canonical carrier.
    module, *_ = cross.build_carrier_data()
    sector_words = prior.sector_words(module)
    base = scout.merge_intervals(sector_words[0] + sector_words[1])
    period = cross.T
    full_word = [
        (Fraction(left, period), Fraction(right, period))
        for left, right in base
    ]
    require((C2, C3, C3 // C2) == (P**3, 2 * P**5, 338),
            "canonical target/high-speed ratio changed")
    require(C3 % C2 == 0 and C3 // C2 < P**3,
            "uniform four-target high-speed obstruction changed")

    # All two-edge phase policies survive to state three.  Three exact
    # hostile controls are replayed at state four; the proof for every one of
    # the 13^3 words is the phase-independent 338*z contraction above.
    depth_three, targets = phase_depth_three_atlas(base, period)
    require(len(depth_three) == P**2 and
            all(count > 0 for count, _ in depth_three.values()),
            "some phase word died before state three")
    minimum = min(depth_three.items(), key=lambda item: item[1][1])
    maximum = max(depth_three.items(), key=lambda item: item[1][1])
    require((minimum[0], minimum[1][0], minimum[1][1])
            == EXPECTED_PHASE_DEPTH3["minimum"],
            "depth-three minimum changed")
    require((maximum[0], maximum[1][0], maximum[1][1])
            == EXPECTED_PHASE_DEPTH3["maximum"],
            "depth-three maximum changed")
    depth_four_controls = {
        word: phase_word_row(base, targets, period, word)
        for word in ((0, 0, 0), (1, 2, 3), (12, 7, 4))
    }
    require(all(row == (0, Fraction(0))
                for row in depth_four_controls.values()),
            "a depth-four hostile phase word became positive")

    # The slope-14 map is the first integer slope above the no-wrap threshold.
    # It greatly delays the full-word zero but does not create recurrence.
    slope_14 = slope_rows(full_word, 14, 15)
    require(slope_14 == EXPECTED_SLOPE_14,
            "slope-14 exact tail transcript changed")
    require(not danger(Fraction(1, 15) * 338),
            "target-only slope-14 two-cycle lost its high-speed escape")
    require(danger(Fraction(1, 15)) and danger(Fraction(-1, 15)),
            "target-only slope-14 two-cycle left the target tooth")

    print("LRC14 C13^7 phase lift and slope-14 boundary")
    print("status=VERIFIED-EXACT SIDECAR; not row exclusion or LRC14")
    print(f"parameters=(p={P},R={R},c2={C2},c3={C3},c3/c2={C3//C2})")
    print("physical_lift=T_s^(7)(x)={13x+s/13^7}; tail F_s(y)={13y+s/13}")
    print(f"minimal_fixed_points={len(fixed_lifts)}/12; min_speed14_safe_margin={min(margins)}")
    print(f"s=1_fixed_lift=(x={fixed_lifts[0][1]},y={fixed_lifts[0][2]})")
    print("minimal_carrier=12 recurrent interior fixed points for s=1,...,12")
    print("full_W_uniform_stop=four target teeth force z in I/13^3; c3-safety fails because c3*y=338z and 338<13^3")
    print("phase_language=(depth1 1/1, depth2 13/13, depth3 169/169, depth4 0/2197 by uniform proof)")
    print(f"depth3_min={EXPECTED_PHASE_DEPTH3['minimum']}")
    print(f"depth3_max={EXPECTED_PHASE_DEPTH3['maximum']}")
    print(f"depth4_exact_controls={depth_four_controls}")
    print("slope14_full_W_rows=")
    for row in slope_14:
        print(row)
    print("slope14_target_control=(1/15,-1/15) is a target-danger 2-cycle with 338/15 high-speed-safe, but full W dies at state15")
    print("scope=the C13^7 lift is a genuine circle map but lacks inherited THM-2657 carry/root and present/clock typing")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
