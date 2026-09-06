#!/usr/bin/env python3
"""Sharp all-height mixed-parity sum-ray certificate; complete c<33 head.

The head is all primitive 0<a<b, c=a+b<33 with 3 not dividing abc.
Independent paths: closed residue-count sum, raw lattice carriers, literal
six-sheet interval networks, and exact periodic-antiderivative identities.
The imported earlier script supplies the last two paths, not this proof.
"""

from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from math import gcd
from pathlib import Path


SPEC = spec_from_file_location(
    "literal_raw_referee", Path(__file__).with_name("lrc14_one_ray_overnight_hexagon_sep05.py")
)
BASE = module_from_spec(SPEC)
SPEC.loader.exec_module(BASE)
CHECKS = 0
SHARP = Q(6, 55)
ODD_TARGET = Q(6, 77)


def need(value, payload):
    global CHECKS
    CHECKS += 1
    if not value:
        raise AssertionError(payload)


def count(k):
    return k - k // 3


def total(k):
    return k * (k + 1) // 2 - 3 * (k // 3) * (k // 3 + 1) // 2


def discrepancy(t):
    return t.numerator // t.denominator - (t / 3).numerator // (t / 3).denominator - Q(2, 3) * t


def primitive(t):
    r = t - 3 * (t // 3)
    if r <= 1:
        return -r * r / 3
    if r <= 2:
        return r - 1 - r * r / 3
    return -(3 - r) ** 2 / 3


def closed_data(a, b):
    c = a + b
    k = (3 * c - 1) // 14
    n, s = count(k), total(k)
    q = Q(3, 7 * c)
    transitions = (Q(3 * a, 14), Q(3 * b, 14), Q(3 * (a * a + b * b), 14 * c))
    projections = []
    for pair, threshold in zip(((b, c), (a, c), (a, b)), transitions):
        j, ell = pair
        flat = min(k, threshold.numerator // threshold.denominator)
        nf, sf = count(flat), total(flat)
        projections.append(2 * (q * nf + Q(3 * (j + ell) * (n - nf) - 14 * (s - sf), 14 * j * ell)))
    B = Q(3 * c, 14)
    ratio = Q(a, c)
    continuous = ((9 + 3 * ratio) / 98, (12 - 3 * ratio) / 98,
                  Q(6, 49) * (1 - ratio + ratio * ratio))
    for index, pair in enumerate(((b, c), (a, c), (a, b))):
        end_value = q / 2 if index < 2 else Q(0)
        exact = continuous[index] + 2 * end_value * discrepancy(B) + Q(2, pair[0] * pair[1]) * (primitive(B) - primitive(transitions[index]))
        need(exact == projections[index], ("full boundary-discrepancy identity", a, b, index))
    physical = Q(9, 98) + Q(2, b * c) * (primitive(transitions[1]) - primitive(transitions[0])) + Q(2, a * b) * (primitive(B) - primitive(transitions[1]))
    need(abs(physical - Q(9, 98)) <= Q(2, 3 * a * b), ("physical endpoint error", a, b))
    return tuple(projections), physical


def audit(a, b, literal=True):
    c = a + b
    w = (a, b, c)
    live = BASE.carriers(w)
    bound = (3 * c - 1) // 14
    expected = {(sign * k, sign * k, -sign * k)
                for k in range(1, bound + 1) if k % 3 for sign in (-1, 1)}
    need(live == expected, ("complete sum-ray dictionary", w, live ^ expected))
    raw = BASE.projection_data(w, live)
    closed = closed_data(a, b)
    need(raw == closed, ("raw versus closed count", w, raw, closed))
    if literal:
        need(BASE.literal_projection_data(w) == raw, ("literal six-sheet paths", w))
    need(min(raw[0]) == min(raw[0][0], raw[0][2]), ("middle selector dominated", w))
    return raw


def main():
    # Every strict raw cutoff which is integral is a deleted multiple of 3.
    for c in range(1, 200):
        B = Q(3 * c, 14)
        need(B.denominator != 1 or B.numerator % 3 == 0, ("strict cutoff", c))
    for denominator in range(1, 15):
        for numerator in range(42 * denominator + 1):
            t = Q(numerator, denominator)
            need(Q(-2, 3) <= discrepancy(t) <= Q(2, 3), ("residue discrepancy range", t))
            need(Q(-1, 3) <= primitive(t) <= 0, ("periodic primitive range", t))
            need(primitive(t + 3) == primitive(t), ("periodic primitive", t))
    tail_a = Q(39, 392) + Q(2, 7 * 33) + Q(8, 9 * 33**2)
    tail_c = Q(39, 392) + Q(32, 9 * 33**2)
    need(tail_a < SHARP and tail_c < SHARP, "all-height tail at33")
    need(Q(9, 98) - Q(2, 3 * 48) > ODD_TARGET, "unbounded actual-mass hostile at49")
    rows = []
    digest = sha256()
    for c in range(2, 33):
        for a in range(1, (c + 1) // 2):
            b = c - a
            if a % 3 and b % 3 and c % 3 and gcd(a, b) == 1:
                projections, physical = audit(a, b)
                need(min(projections) <= SHARP, ("complete head bound", a, b))
                row = (a, b, c, projections, physical)
                rows.append(row)
                digest.update((repr(row) + "\n").encode())
    eq_network = [row[:3] for row in rows if min(row[3]) == SHARP]
    eq_physical = [row[:3] for row in rows if row[4] == SHARP]
    need(len(rows) == 42, ("complete head count", len(rows)))
    need(eq_network == eq_physical == [(1, 10, 11)], ("sharp equality", eq_network, eq_physical))
    print("COMPLETE PRIMITIVE SUM-RAY HEAD c<33 ROWS", len(rows))
    print("SEMANTIC SHA256", digest.hexdigest())
    print("SHARP min-projection AND actual-mass", SHARP, "EQUALITY", eq_network)
    print("TAIL AT33", tail_a, tail_c, "<", SHARP)
    for w in ((2, 5, 7), (1, 10, 11), (2, 47, 49), (1, 49, 50), (7, 43, 50), (31, 67, 98)):
        projections, physical = audit(w[0], w[1])
        if w[2] >= 49:
            need(physical > ODD_TARGET, ("actual parity obstruction", w))
        print("CONTROL", w, "E", projections, "ACTUAL", physical)
    # The discarded boundary term really matters at the sharp hostile.
    continuous_a = Q(51, 539)
    need(SHARP - continuous_a > Q(2, 3 * 10 * 11), "pure quadratic projection error is false")
    print("CHECKS", CHECKS, "INDEPENDENT RAW/LITERAL CHECKS", BASE.CHECKS)
    print("PROVED sharp sum-ray6/55; every primitive sum ray c>=49 actually exceeds6/77")
    print("OPEN arbitrary-body entry/synchronization; no inherited6/55 Haar floor")


if __name__ == "__main__":
    main()
