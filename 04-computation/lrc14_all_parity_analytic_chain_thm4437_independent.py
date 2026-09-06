#!/usr/bin/env python3
"""Exact arithmetic checks for the parity-free analytic reduction."""

from fractions import Fraction as Q
from itertools import product


CHECKS = 0


def need(condition, payload):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(payload)


def dot(x, y):
    return sum(a * b for a, b in zip(x, y)) % 3


def cross(x, y):
    return (
        (x[1] * y[2] - x[2] * y[1]) % 3,
        (x[2] * y[0] - x[0] * y[2]) % 3,
        (x[0] * y[1] - x[1] * y[0]) % 3,
    )


def mod3_audit():
    unit_speeds = list(product((1, 2), repeat=3))
    relations = 0
    cases = {0: 0, 1: 0}
    observations = 0
    for w in unit_speeds:
        live = [C for C in product((1, 2), repeat=3) if dot(C, w) == 0]
        need(len(live) == 2, ("live plane", w, live))
        for v in product(range(3), repeat=3):
            if v == (0, 0, 0) or dot(v, w):
                continue
            zeros = sum(x == 0 for x in v)
            need(zeros in (0, 1), ("two-zero relation", w, v))
            cases[zeros] += 1
            relations += 1
            for C in live:
                vc = cross(v, C)
                deltas = {(vc[i] * pow(w[i], -1, 3)) % 3 for i in range(3)}
                need(len(deltas) == 1, ("cross not parallel", w, v, C, vc))
                delta = next(iter(deltas))
                longitudinal = sum(
                    all((C[i] + t * v[i]) % 3 for i in range(3))
                    for t in range(3)
                )
                if zeros == 0:
                    need(delta == 0 and longitudinal == 2,
                         ("unit dichotomy", w, v, C, delta, longitudinal))
                else:
                    need(delta != 0 and longitudinal == 1,
                         ("one-zero dichotomy", w, v, C, delta, longitudinal))
                observations += 1
    return relations, cases, observations


def defect_count_audit():
    max_s = 10000
    empty_onezero = []
    for S in range(1, max_s + 1):
        bound = (3 * S - 1) // 14
        unit = sum(d % 3 == 0 for d in range(-bound, bound + 1))
        onezero = 2 * bound + 1 - unit
        need(Q(4, 3) * unit < Q(4, 21) * S + Q(4, 3),
             ("unit intercept", S, unit))
        need(onezero < Q(2, 7) * S + Q(4, 3),
             ("one-zero intercept", S, onezero))
        if onezero == 0:
            empty_onezero.append(S)
        if S >= 5:
            need(onezero >= 2, ("missing +/-1 defects", S, onezero))
    need(empty_onezero == [1, 2, 3, 4], ("empty defect range", empty_onezero))
    return max_s, empty_onezero


def cutoff_audit():
    threshold = lambda S: Q(308, 31) * S + Q(4312, 93)
    quadratic = lambda S: Q(3, 16) * S * S
    need(threshold(57) == Q(56980, 93) < 613,
         ("small-S cutoff", threshold(57)))
    need(quadratic(58) > threshold(58),
         ("large-S base", quadratic(58), threshold(58)))
    # q(S)-threshold(S) has positive discrete increments for S>=58.
    for S in range(58, 10000):
        need(quadratic(S) > threshold(S), ("large-S cutoff", S))
    need(612 % 3 == 0 and 611 % 3 and 613 % 3,
         "eligible-head boundary")
    return threshold(57), quadratic(58), threshold(58)


def main():
    relations, cases, observations = mod3_audit()
    max_s, empty = defect_count_audit()
    t57, q58, t58 = cutoff_audit()
    large_coefficient = Q(6, 49) + Q(4, 7 * 19)
    need(large_coefficient == Q(142, 931) < Q(15, 98),
         ("M>=19", large_coefficient))
    print("CLEANROOM PARITY-FREE ANALYTIC CHECKS")
    print("MOD3 relations", relations, "zero_count_cases", cases,
          "live_observations", observations)
    print("DEFECT_COUNTS checked_S_through", max_s,
          "onezero_empty_exactly", empty)
    print("M_GE_19", large_coefficient, "TARGET", Q(15, 98))
    print("CUTOFF threshold57", t57, "quadratic58", q58,
          "threshold58", t58, "tail_c>=613 head_c<=611")
    print("CHECKS", CHECKS)
    print("PASS")


if __name__ == "__main__":
    main()
