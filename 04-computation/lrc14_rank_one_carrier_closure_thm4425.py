#!/usr/bin/env python3
"""Exact finite and arithmetic verifier for THM-4425.

The all-height step is elementary.  This companion independently rebuilds
the complete c<35 base from the defining integer carrier lattice and checks
the parity/ternary-unit coefficient classification used by the proof.
"""

from fractions import Fraction as Q
from itertools import combinations, product
from math import gcd
import sys


sys.stdout.reconfigure(newline="\n")


TARGET = Q(6, 77)
CHECKS = 0


class VerificationError(RuntimeError):
    pass


def need(condition, payload):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise VerificationError(payload)


def strict_floor(numerator, denominator):
    return (numerator - 1) // denominator


def carriers(w):
    a, b, c = w
    bounds = (
        strict_floor(3 * (b + c), 14),
        strict_floor(3 * (a + c), 14),
        strict_floor(3 * (a + b), 14),
    )
    result = []
    for first in range(-bounds[0], bounds[0] + 1):
        for second in range(-bounds[1], bounds[1] + 1):
            numerator = -(a * first + b * second)
            if numerator % c:
                continue
            third = numerator // c
            row = (first, second, third)
            if abs(third) > bounds[2] or not all(value % 3 for value in row):
                continue
            need(a * first + b * second + c * third == 0,
                 ("kernel", w, row))
            result.append(row)
    return tuple(result)


def primitive_direction(row):
    divisor = gcd(gcd(abs(row[0]), abs(row[1])), abs(row[2]))
    direction = tuple(value // divisor for value in row)
    for value in direction:
        if value:
            if value < 0:
                direction = tuple(-entry for entry in direction)
            break
    return direction


def rank_one(rows):
    return not rows or len({primitive_direction(row) for row in rows}) == 1


def complete_line(w, direction):
    bounds = (
        strict_floor(3 * (w[1] + w[2]), 14 * abs(direction[0])),
        strict_floor(3 * (w[0] + w[2]), 14 * abs(direction[1])),
        strict_floor(3 * (w[0] + w[1]), 14 * abs(direction[2])),
    )
    last = min(bounds)
    result = set()
    for address in range(1, last + 1):
        if address % 3:
            row = tuple(address * value for value in direction)
            result.add(row)
            result.add(tuple(-value for value in row))
    return result


def projections(w, rows):
    a, b, c = w
    denominator = 14 * a * b * c
    sums = [0, 0, 0]
    for first, second, third in rows:
        margins = (
            3 * (a + b) - 14 * abs(third),
            3 * (a + c) - 14 * abs(second),
            3 * (b + c) - 14 * abs(first),
        )
        need(all(value > 0 for value in margins),
             ("roof", w, (first, second, third), margins))
        cap = 6 * a * b
        terms = (
            min(c * margins[0], cap),
            min(b * margins[1], cap),
            min(a * margins[2], cap),
        )
        for index in range(3):
            sums[index] += terms[index]
    return tuple(Q(value, denominator) for value in sums)


def coefficient_gate():
    patterns = set()
    signed_relations = 0
    for magnitudes in product((1, 2), repeat=3):
        if gcd(gcd(*magnitudes[:2]), magnitudes[2]) != 1:
            continue
        if sum(magnitudes) % 2:
            continue
        need(all(value % 3 for value in magnitudes),
             ("ternary-unit", magnitudes))
        patterns.add(tuple(sorted(magnitudes)))
        for signs in product((-1, 1), repeat=3):
            row = tuple(signs[i] * magnitudes[i] for i in range(3))
            need(sum(row) % 2 == 0, ("signed-parity", row))
            signed_relations += 1
    need(patterns == {(1, 1, 2)}, ("small-patterns", patterns))
    return patterns, signed_relations


def finite_base():
    values = tuple(value for value in range(1, 35, 2) if value % 3)
    triples = 0
    empty = 0
    one_direction = 0
    non_norm_four = 0
    equality = []
    leader = None
    semantic = []
    for w in combinations(values, 3):
        if gcd(gcd(w[0], w[1]), w[2]) != 1:
            continue
        triples += 1
        rows = carriers(w)
        if not rows:
            empty += 1
        if not rank_one(rows):
            continue
        one_direction += 1
        sums = projections(w, rows)
        best = min(sums)
        need(best <= TARGET, ("base-target", w, rows, sums))
        if best == TARGET:
            equality.append(w)
        if rows:
            direction = primitive_direction(rows[0])
            expected = complete_line(w, direction)
            need(set(rows) == expected,
                 ("complete-natural-address", w, direction,
                  sorted(set(rows).symmetric_difference(expected))))
            norm_four = tuple(sorted(abs(value) for value in direction)) == (1, 1, 2)
            if not norm_four:
                non_norm_four += 1
        else:
            direction = None
            norm_four = False
        row = (best, w, len(rows), direction, sums)
        if leader is None or row > leader:
            leader = row
        semantic.append((w, len(rows), direction, sums))
    need(equality == [(1, 5, 11)], ("base-equality", equality))
    return values, triples, empty, one_direction, non_norm_four, leader, semantic


def arithmetic_gate():
    threshold = Q(308, 9)
    need(Q(35) >= threshold, ("height-threshold", threshold))
    need(Q(33) < threshold, ("sharp-integer-predecessor", threshold))
    # If max |u_i| >= 4, B < 3c/28 and the block count gives
    # N < c/7+4/3 <= 2c/11 at c >= 308/9.
    need(Q(1, 7) + Q(4, 3 * 35) <= Q(2, 11),
         ("normalized-count-bound", Q(1, 7) + Q(4, 105)))
    return threshold


def main():
    patterns, signed = coefficient_gate()
    threshold = arithmetic_gate()
    values, triples, empty, one_direction, non_norm_four, leader, semantic = finite_base()
    print("LRC14 RANK-ONE CARRIER CLOSURE THM-4425 VERIFIER")
    print("status=PASS PROVED_RELATIVE+FINITE_EXACT_BASE; universal_projection=OPEN; LRC14=OPEN")
    print("small_coefficient_patterns=%s signed_parity_controls=%d" %
          (sorted(patterns), signed))
    print("analytic_non_norm4_gate=max_abs_u>=4 B<3c/28 N<c/7+4/3")
    print("automatic_height_threshold=%s first_eligible_height=35" % threshold)
    print("finite_base_values=%s triples=%d empty=%d rank_at_most_one=%d non_norm4_nonempty=%d" %
          (values, triples, empty, one_direction, non_norm_four))
    print("finite_base_leader=w:%s,best:%s,carriers:%d,direction:%s,sums:%s" %
          (leader[1], leader[0], leader[2], leader[3],
           tuple(str(value) for value in leader[4])))
    print("finite_base_equality=[(1,5,11)] semantic_rows=%d" % len(semantic))
    print("checks=%d" % CHECKS)
    print("verdict=PASS")


if __name__ == "__main__":
    main()
