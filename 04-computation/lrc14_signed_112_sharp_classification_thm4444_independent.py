#!/usr/bin/env python3
"""Clean-room literal referee for the signed (1,1,2) classification.

This program imports only the pre-existing mixed-parity literal sheet engine,
not the candidate producer.  It constructs signed directions independently,
enumerates every typed row through height 611, checks raw carriers against the
predicted primitive ray, and recomputes all sharp and 6/77 boundary rows.
"""

from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations, permutations, product
import importlib.util
from math import gcd
from pathlib import Path


PINNED = "b99af309ff6f8643dedf923f5ee8d86d67b32ff2b0d6510209c565820894f399"
SOURCE = Path("04-computation/lrc14_parity_empty_core_sep06.py")
R = Q(3, 14)
TARGET = Q(6, 77)
CHECKS = 0


def need(test, detail):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(detail)


def load_engine():
    need(sha256(SOURCE.read_bytes()).hexdigest() == PINNED, "engine hash")
    spec = importlib.util.spec_from_file_location("literal_parity_engine", SOURCE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def directions():
    answer = set()
    for magnitudes in set(permutations((1, 1, 2))):
        for signs in product((-1, 1), repeat=3):
            u = tuple(magnitudes[i] * signs[i] for i in range(3))
            if u[0] < 0:
                u = tuple(-x for x in u)
            answer.add(u)
    return answer


ALL_DIRECTIONS = directions()
SORTED_DIRECTIONS = {
    (2, 1, -1): 1,
    (1, -2, 1): 2,
    (1, 2, -1): 3,
}


def typed(w):
    return (0 < w[0] < w[1] < w[2]
            and gcd(gcd(w[0], w[1]), w[2]) == 1
            and all(value % 3 for value in w))


def signed_relations(w):
    return tuple(u for u in ALL_DIRECTIONS
                 if sum(u[i] * w[i] for i in range(3)) == 0)


def rows(height):
    answer = []
    for c in range(1, height + 1):
        if c % 3 == 0:
            continue
        for a in range(1, c):
            if a % 3 == 0:
                continue
            for b in {
                c - 2 * a,
                (a + c) // 2 if (a + c) % 2 == 0 else 0,
                (c - a) // 2 if (c - a) % 2 == 0 else 0,
            }:
                w = (a, b, c)
                if not typed(w):
                    continue
                relations = signed_relations(w)
                need(len(relations) == 1, ("unique signed cone", w, relations))
                u = relations[0]
                need(u in SORTED_DIRECTIONS, ("unexpected direction", w, u))
                answer.append((w, u, SORTED_DIRECTIONS[u]))
    need(len(answer) == len({w for w, _, _ in answer}), "duplicate rows")
    return sorted(answer)


def predicted_carriers(w, u):
    bounds = tuple((3 * (sum(w) - w[i]) - 1) // 14 for i in range(3))
    largest = min(bounds[i] // abs(u[i]) for i in range(3))
    return {
        tuple(sign * k * u[i] for i in range(3))
        for k in range(1, largest + 1) if k % 3
        for sign in (-1, 1)
    }


def continuum(w, u):
    """Generic breakpoint integration, independent of cone formulas."""
    z = tuple(Q(value, w[2]) for value in w)
    lines = []
    support = None
    points = {Q(0)}
    for i in range(3):
        j, k = tuple(index for index in range(3) if index != i)
        alpha = R * (z[j] + z[k]) / (z[j] * z[k])
        beta = Q(abs(u[i]), 1) / (z[j] * z[k])
        lines.append((alpha, beta))
        zero = alpha / beta
        support = zero if support is None else min(support, zero)
        cap = (alpha - 2 * R) / beta
        if cap > 0:
            points.add(cap)
    points.add(support)
    for i in range(3):
        for j in range(i):
            ai, bi = lines[i]
            aj, bj = lines[j]
            if bi != bj:
                crossing = (ai - aj) / (bi - bj)
                if 0 < crossing < support:
                    points.add(crossing)
    points = sorted(point for point in points if 0 <= point <= support)
    integrals = [Q(0), Q(0), Q(0)]
    envelope = Q(0)

    def value(i, x):
        alpha, beta = lines[i]
        return max(Q(0), min(2 * R, alpha - beta * x))

    for left, right in zip(points, points[1:]):
        midpoint = (left + right) / 2
        owner = min(range(3), key=lambda i: value(i, midpoint))
        envelope += (right - left) * (value(owner, left) + value(owner, right)) / 2
        for i in range(3):
            integrals[i] += (right - left) * (value(i, left) + value(i, right)) / 2
    return tuple(integrals), envelope


def update(table, key, value, payload):
    old = table.get(key)
    if old is None or value > old[0]:
        table[key] = (value, [payload])
    elif value == old[0]:
        old[1].append(payload)


def main():
    engine = load_engine()
    need(len(ALL_DIRECTIONS) == 12, ("direction count", ALL_DIRECTIONS))
    sample_viable = set()
    for w in combinations((x for x in range(1, 101) if x % 3), 3):
        if gcd(gcd(*w[:2]), w[2]) != 1:
            continue
        sample_viable.update(signed_relations(w))
    need(sample_viable == set(SORTED_DIRECTIONS), ("viable cones", sample_viable))

    data = rows(611)
    counts = {1: 0, 2: 0, 3: 0}
    leaders = {}
    above = []
    equal = []
    for w, u, family in data:
        E, physical = engine.native(w)
        Er, pr, carriers = engine.raw(w)
        need((E, physical) == (Er, pr), ("literal/raw", w))
        need(carriers == predicted_carriers(w, u), ("complete ray", w, u))
        selected = next(i for i in range(3) if abs(u[i]) == 2)
        need(physical == E[selected] == min(E), ("selector", w, u, E, physical))
        integrals, envelope = continuum(w, u)
        need(integrals[selected] == envelope == R * R,
             ("continuum integral", w, u, integrals, envelope))
        error = Q(4, 7 * w[2])
        need(Q(3, 49) - error <= physical < Q(3, 49) + error,
             ("quadrature band", w, physical, error))
        counts[family] += 1
        update(leaders, family, physical, (w, E))
        if physical > TARGET:
            above.append((w, physical))
        elif physical == TARGET:
            equal.append((w, physical))

    expected = {
        1: (Q(58, 833), [((5, 7, 17), (Q(58, 833), Q(12, 119), Q(346, 4165)))]),
        2: (Q(11, 140), [((2, 11, 20), (Q(131, 1540), Q(11, 140), Q(3, 35)))]),
        3: (TARGET, [((1, 5, 11), (TARGET, TARGET, TARGET))]),
    }
    need(counts == {1: 9482, 2: 14214, 3: 4742}, ("counts", counts))
    need(leaders == expected, ("leaders", leaders))
    need(above == [((2, 11, 20), Q(11, 140))], ("above", above))
    need(equal == [((1, 5, 11), TARGET)], ("equal", equal))
    need(Q(3, 49) + Q(4, 7 * 35) < TARGET, "target tail cutoff")
    print("INDEPENDENT_LITERAL_SIGNED_112_REFEREE=PASS")
    print("pinned_engine_sha256", PINNED)
    print("directions", len(ALL_DIRECTIONS), "viable_sorted", sorted(sample_viable))
    print("H611_rows", len(data), "family_counts", counts)
    print("family_leaders", leaders)
    print("above_6_77", above)
    print("equal_6_77", equal)
    print("continuum_selected_integral", R * R, "universal_target_tail_height", 35)
    print("new_checks", CHECKS, "engine_checks", engine.CHECKS)


if __name__ == "__main__":
    main()
