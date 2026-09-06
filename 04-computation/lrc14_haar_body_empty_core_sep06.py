#!/usr/bin/env python3
"""Exact, bounded ten-body Haar probes at weak clearance 1/14.

Primary: intersections of closed safe teeth, including isolated points.
Independent: union of danger teeth for measure; wall/midpoint scan for the
complete closed geometry. No old computation module is imported.
Run: python3 -B 04-computation/lrc14_haar_body_empty_core_sep06.py
"""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations
from math import gcd

DELTA = Q(1, 14)
TARGET = Q(6, 77)
TAIL = (1, 5, 11)
CHECKS = 0


def need(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def merge(intervals):
    result = []
    for left, right in sorted(intervals):
        if result and left <= result[-1][1]:
            result[-1] = result[-1][0], max(right, result[-1][1])
        else:
            result.append((left, right))
    return tuple(result)


def intersect(first, second):
    i = j = 0
    result = []
    while i < len(first) and j < len(second):
        left, right = max(first[i][0], second[j][0]), min(first[i][1], second[j][1])
        if left <= right:
            result.append((left, right))
        x, y = first[i][1], second[j][1]
        i += x <= y
        j += y <= x
    return merge(result)


@lru_cache(None)
def safe_teeth(speed):
    return tuple((Q(14*k + 1, 14*speed), Q(14*k + 13, 14*speed)) for k in range(speed))


@lru_cache(None)
def safe(body):
    if not body:
        return ((Q(0), Q(1)),)
    return intersect(safe(body[:-1]), safe_teeth(body[-1]))


def measure(intervals):
    return sum((right - left for left, right in intervals), Q(0))


def clearance(body, phase):
    return min(min((speed*phase) % 1, 1 - (speed*phase) % 1) for speed in body)


def danger_union_measure(body):
    # Open/closed endpoint choices have no effect on this independent measure.
    teeth = []
    for speed in body:
        for owner in range(speed + 1):
            left, right = Q(14*owner - 1, 14*speed), Q(14*owner + 1, 14*speed)
            left, right = max(Q(0), left), min(Q(1), right)
            if left < right:
                teeth.append((left, right))
    return 1 - measure(merge(teeth))


def wall_safe(body):
    walls = {Q(0), Q(1)}
    for speed in body:
        for owner in range(speed + 1):
            for sign in (-1, 1):
                point = Q(14*owner + sign, 14*speed)
                if 0 <= point <= 1:
                    walls.add(point)
    walls = sorted(walls)
    pieces = [(point, point) for point in walls if clearance(body, point) >= DELTA]
    for left, right in zip(walls, walls[1:]):
        if clearance(body, (left + right) / 2) >= DELTA:
            pieces.append((left, right))
    return merge(pieces)


def completed(body):
    row = tuple(sorted(tuple(3*x for x in body) + TAIL))
    need(len(set(row)) == 13, "ten divided speeds and three distinct exceptional tails")
    return row


def small_clock_escapes(row):
    return tuple(q for q in range(2, 15) if all(speed % q for speed in row))


def first_rational_clock(row, maximum=100):
    for denominator in range(2, maximum + 1):
        for numerator in range(1, denominator // 2 + 1):
            if gcd(numerator, denominator) != 1:
                continue
            phase = Q(numerator, denominator)
            gap = clearance(row, phase)
            if gap >= DELTA:
                return phase, gap
    return None


def show_body(label, body, geometry=False):
    intervals = safe(body)
    mass = measure(intervals)
    need(mass == danger_union_measure(body), "independent danger-union measure " + label)
    need(intervals == wall_safe(body), "independent complete closed-wall geometry " + label)
    points = tuple(left for left, right in intervals if left == right)
    positive = sum(left < right for left, right in intervals)
    print("BODY", label, body, "mu", mass, "mu-minus-6/77", mass-TARGET,
          "positive_components", positive, "isolated", tuple(map(str, points)))
    if geometry:
        print("CLOSED_COMPONENTS", label, tuple((str(left), str(right)) for left, right in intervals))
    return mass


def main():
    print("FINITE-EXACT / recovered obstruction: ten-body universal 6/77 Haar floor is false")
    print("threshold weak ||cy||>=1/14; phase y in R/Z; normalized Haar measure")
    canonical = tuple(range(1, 11))
    small = (1, 2, 3, 4, 5, 7, 8, 9, 11, 13)
    minimum13 = (1, 2, 3, 5, 7, 8, 9, 11, 12, 13)
    clock_hostile = (1, 3, 4, 10, 11, 13, 14, 16, 17, 18)
    need(show_body("consecutive_1_to_10", canonical) == Q(1217, 8820), "inherited AP measure")
    for scale in (2, 7):
        need(show_body("dilation_" + str(scale), tuple(scale*x for x in canonical)) == Q(1217, 8820),
             "Haar dilation preserves full measure")
    show_body("shifted_AP_2_to_11", tuple(range(2, 12)))
    show_body("odd_AP_1_to_19", tuple(range(1, 20, 2)))
    for far in (11, 13, 21, 101):
        show_body("one_far_" + str(far), tuple(range(1, 10)) + (far,))
    need(show_body("first_height13_hostile", small, True) == Q(21514, 315315) < TARGET,
         "small exact refutation")
    need(show_body("minimum_height13", minimum13, True) == Q(14249, 252252) < TARGET,
         "stronger inherited height13 obstruction")
    need(show_body("clock_filtered_hostile", clock_hostile, True) == Q(534689, 7796880) < TARGET,
         "necessary small-clock sieve does not imply Haar floor")

    # Explain the small near-AP mechanism by exact mass removed at each step.
    growing = (1, 2, 3, 5, 7, 8, 9)
    previous = measure(safe(growing))
    print("REMOVAL_CHAIN base", growing, "mu", previous)
    for added in (11, 12, 13):
        growing += (added,)
        current = measure(safe(growing))
        print("REMOVAL_CHAIN add", added, "mu", current, "removed", previous-current)
        previous = current
    need(previous == Q(14249, 252252), "small near-AP removal mechanism")

    # Complete, very small first-height audit: choose ten labels from 1..13.
    census = [(measure(safe(body)), body) for body in combinations(range(1, 14), 10)]
    below12 = [(m, body) for m, body in census if body[-1] <= 12]
    need(len(census) == 286 and len(below12) == 66, "exact small universes")
    need(all(m >= TARGET for m, _ in below12), "no smaller-height ten-body hostile")
    need(min(census) == (Q(14249, 252252), minimum13), "complete height13 minimum")
    need(sum(m < TARGET for m, _ in census) == 12, "height13 hostile count")
    for m, body in census:
        need(m == danger_union_measure(body), "independent complete height13 measure census")
    print("COMPLETE_HEIGHT13", len(census), "bodies; max<=12", len(below12),
          "all pass; belowfloor", 12, "minimum", min(census))

    # Necessary entry sieve for the fixed sharp tail: every denominator<=14
    # has a divisible speed in 3C union T. This is not sufficient entry.
    first_clock_hostile = None
    for height in range(14, 19):
        count = 0
        leader = (Q(1), None)
        for prefix in combinations(range(1, height), 9):
            body = prefix + (height,)
            row = tuple(3*x for x in body) + TAIL
            if small_clock_escapes(row):
                continue
            count += 1
            m = measure(safe(body))
            if m < leader[0]:
                leader = m, body
            if m < TARGET:
                first_clock_hostile = body
                break
        print("CLOCK_FILTERED_HEIGHT", height, "examined", count,
              "minimum_examined", leader, "stopped_at_first_hostile", first_clock_hostile is not None)
        if first_clock_hostile is not None:
            break
    need(first_clock_hostile == clock_hostile, "first filtered hostile and smallest filtered height")

    for label, body, expected_phase, expected_gap, expected_mass in [
        ("small", small, Q(1, 10), Q(1, 10), Q(23818, 945945)),
        ("clock_filtered", clock_hostile, Q(9, 19), Q(2, 19), Q(8131, 194480)),
    ]:
        row = completed(body)
        intervals = safe(row)
        need(intervals == wall_safe(row), "full-row independent closed geometry")
        need(measure(intervals) == danger_union_measure(row) == expected_mass, "full-row mass")
        phase, gap = first_rational_clock(row)
        need((phase, gap) == (expected_phase, expected_gap), "actual positive completion phase")
        print("COMPLETION", label, row, "mu_safe", expected_mass,
              "small_unit_clock_escapes", small_clock_escapes(row),
              "first_rational_clock", str(phase), "clearance", str(gap),
              "body_phase", str((3*phase) % 1))
    need(not small_clock_escapes(completed(clock_hostile)), "fixed-clock inheritance exclusion")

    # Inherited THM-4032 selector hostile, with the original phase retained.
    phase = Q(2, 11)
    need(clearance(canonical, phase) >= DELTA, "inherited safe divided-pack phase")
    bad = tuple(tuple(j for j in range(3) if clearance((tail,), (phase+j)/3) < DELTA) for tail in TAIL)
    need(bad == ((0,), (1,), (2,)), "three original lifts spoiled by separate owners")
    need(clearance(completed(canonical), Q(1, 14)) == DELTA, "inherited alternative full-row phase")
    print("SELECTOR_HOSTILE canonical body y=2/11; tail bad labels", bad,
          "alternative full-row x=1/14 clearance1/14")
    print("CHECKS", CHECKS)
    print("PASS: universal Haar floor refuted; necessary clock-sieved floor also refuted; completions remain safe")
    print("SCOPE: recovered exact bodies, new threshold application; no new LRC family or entry theorem")


if __name__ == "__main__":
    main()
