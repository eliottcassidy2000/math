#!/usr/bin/env python3
"""Standalone exact referee for THM-4425's rank-one carrier closure.

No theorem implementation is imported.  Every check uses integers or
fractions and remains active under ``python -O``.
"""

from fractions import Fraction as Q
from itertools import combinations, product
from math import gcd
import sys


sys.stdout.reconfigure(newline="\n")
TARGET = Q(6, 77)
CHECKS = 0


class RefereeFailure(RuntimeError):
    pass


def need(condition, payload):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RefereeFailure(payload)


def strict_floor(num, den):
    return (num - 1) // den


def eligible_values(height):
    return tuple(x for x in range(1, height + 1, 2) if x % 3)


def carriers(w):
    a, b, c = w
    box = (
        strict_floor(3 * (b + c), 14),
        strict_floor(3 * (a + c), 14),
        strict_floor(3 * (a + b), 14),
    )
    out = set()
    for x in range(-box[0], box[0] + 1):
        for y in range(-box[1], box[1] + 1):
            nz = -a * x - b * y
            if nz % c:
                continue
            z = nz // c
            if abs(z) <= box[2] and x % 3 and y % 3 and z % 3:
                need(a * x + b * y + c * z == 0, ("kernel", w, x, y, z))
                out.add((x, y, z))
    return out


def primitive_direction(row):
    d = gcd(gcd(abs(row[0]), abs(row[1])), abs(row[2]))
    out = tuple(x // d for x in row)
    if next(x for x in out if x) < 0:
        out = tuple(-x for x in out)
    return out


def directions(rows):
    return {primitive_direction(row) for row in rows}


def exact_line_dictionary(w, u):
    bounds = tuple(Q(3 * (sum(w) - w[i]), 14 * abs(u[i])) for i in range(3))
    B = min(bounds)
    K = strict_floor(B.numerator, B.denominator)
    predicted = {
        tuple(k * x for x in u)
        for k in range(-K, K + 1)
        if k and k % 3
    }
    return B, predicted


def projection_sums(w, rows):
    a, b, c = w
    den = 14 * a * b * c
    nums = [0, 0, 0]
    for C in rows:
        margins = (
            3 * (a + b) - 14 * abs(C[2]),
            3 * (a + c) - 14 * abs(C[1]),
            3 * (b + c) - 14 * abs(C[0]),
        )
        need(all(x > 0 for x in margins), ("roof", w, C, margins))
        terms = (
            min(c * margins[0], 6 * a * b),
            min(b * margins[1], 6 * a * b),
            min(a * margins[2], 6 * a * b),
        )
        for i in range(3):
            nums[i] += terms[i]
    return tuple(Q(x, den) for x in nums)


def coefficient_audit():
    survivors = set()
    for magnitudes in product((1, 2), repeat=3):
        if gcd(gcd(magnitudes[0], magnitudes[1]), magnitudes[2]) != 1:
            continue
        if sum(magnitudes) % 2 == 0:
            survivors.add(tuple(sorted(magnitudes)))
    need(survivors == {(1, 1, 2)}, ("small-directions", survivors))

    # Three is forbidden by the live residue gate.  Thus every other
    # primitive full-support direction has a coordinate of magnitude >=4.
    for row in product(range(-3, 4), repeat=3):
        if not all(x and x % 3 for x in row):
            continue
        if gcd(gcd(abs(row[0]), abs(row[1])), abs(row[2])) != 1:
            continue
        if sum(row) % 2:
            continue
        need(tuple(sorted(abs(x) for x in row)) == (1, 1, 2),
             ("missed-small-direction", row))
    return survivors


def scan(height):
    triples = empty = rank_at_most_one = nonnorm4 = 0
    equality = []
    exact_dictionary_rows = 0
    analytic_tail_controls = 0
    best_nonnorm4 = None
    for w in combinations(eligible_values(height), 3):
        if gcd(gcd(w[0], w[1]), w[2]) != 1:
            continue
        triples += 1
        rows = carriers(w)
        ds = directions(rows)
        if not rows:
            empty += 1
        if len(ds) > 1:
            continue
        rank_at_most_one += 1
        if rows:
            u = next(iter(ds))
            need(all(x % 3 for x in u), ("unit-direction", w, u))
            B, predicted = exact_line_dictionary(w, u)
            need(rows == predicted, ("exact-line-dictionary", w, rows ^ predicted))
            exact_dictionary_rows += 1
            is_norm4 = tuple(sorted(abs(x) for x in u)) == (1, 1, 2)
            if not is_norm4:
                nonnorm4 += 1
                max_abs = max(abs(x) for x in u)
                need(max_abs >= 4, ("coefficient-gap", w, u))
                need(B < Q(3 * w[2], 28), ("B-bound", w, u, B))
                row = (Q(len(rows), w[2]), w, len(rows), u, B)
                if best_nonnorm4 is None or row > best_nonnorm4:
                    best_nonnorm4 = row
                if w[2] >= 35:
                    need(Q(len(rows)) <= Q(2 * w[2], 11),
                         ("tail-count", w, len(rows)))
                    analytic_tail_controls += 1
        sums = projection_sums(w, rows)
        need(min(sums) <= TARGET, ("bounded-projection", w, sums))
        if min(sums) == TARGET:
            equality.append(w)
    return {
        "triples": triples,
        "empty": empty,
        "rank": rank_at_most_one,
        "nonnorm4": nonnorm4,
        "dictionary": exact_dictionary_rows,
        "tail": analytic_tail_controls,
        "equality": equality,
        "best_nonnorm4": best_nonnorm4,
    }


def main():
    patterns = coefficient_audit()
    threshold = Q(308, 9)
    need(Q(1, 7) + Q(4, 3 * 35) <= Q(2, 11), "threshold-35")
    need(Q(1, 7) + Q(4, 3 * 33) > Q(2, 11), "predecessor-33")
    low = scan(31)
    need(low["triples"] == 165 and low["rank"] == 158,
         ("low-counts", low))
    need(low["equality"] == [(1, 5, 11)], ("low-equality", low["equality"]))
    wide = scan(79)
    need(wide["triples"] == 2910 and wide["rank"] == 1304,
         ("wide-counts", wide))
    need(wide["equality"] == [(1, 5, 11)], ("wide-equality", wide["equality"]))

    print("THM-4425 INDEPENDENT RANK-ONE REFEREE")
    print("status=PASS exact-line-dictionary+arithmetic-tail+finite-base; LRC14=OPEN")
    print("small_live_primitive_patterns=%s next_possible_max_abs=4" % sorted(patterns))
    print("exact_identity=Lambda(w)={k*u:abs(k)<B and k_mod_3_nonzero}")
    print("analytic_chain=B<3c/28; N=2R_<(B)<c/7+4/3<=2c/11")
    print("real_threshold=%s first_eligible_height=35" % threshold)
    print("H31 triples=%d empty=%d rank_at_most_one=%d nonnorm4=%d exact_dictionaries=%d" %
          (low["triples"], low["empty"], low["rank"], low["nonnorm4"], low["dictionary"]))
    print("H79 triples=%d empty=%d rank_at_most_one=%d nonnorm4=%d exact_dictionaries=%d tail_controls=%d" %
          (wide["triples"], wide["empty"], wide["rank"], wide["nonnorm4"],
           wide["dictionary"], wide["tail"]))
    best = wide["best_nonnorm4"]
    print("H79 best_nonnorm4_density=%s w=%s N=%d u=%s B=%s" %
          (best[0], best[1], best[2], best[3], best[4]))
    print("unique_equality=(1,5,11) checks=%d" % CHECKS)
    print("verdict=PASS")


if __name__ == "__main__":
    main()
