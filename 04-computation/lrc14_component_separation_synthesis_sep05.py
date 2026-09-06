#!/usr/bin/env python3
"""Exact audit of all-height separated-component identities and the 44/13 tail.

No repository mathematics is imported.  All checks remain live under -O.
The finite census is supporting evidence, not the proof of the infinite tail.
"""

from argparse import ArgumentParser
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations
from math import gcd


CHECKS = 0
TARGET = Q(6, 77)


def check(value, detail):
    global CHECKS
    CHECKS += 1
    if not value:
        raise RuntimeError(detail)


def bernoulli(x):
    x %= 1
    return x * x - x + Q(1, 6)


def correction(p, q):
    return 3 * sum((sign * bernoulli(Q(1, 3) + shift)
                    for shift, sign in ((Q(p-q, 14), 1), (Q(q-p, 14), 1),
                                        (Q(p+q, 14), -1), (Q(-p-q, 14), -1))), Q(0))


def pair_formula(p, q):
    radius = Q(p+q, 14)
    number = (radius + Q(1, 3)).__floor__() + (radius - Q(1, 3)).__floor__() + 1
    mass = Q(6, 49) + correction(p, q) / (p*q)
    roof = Q(6, p*q) * sum((max(Q(0), min(Q(p, 7), radius-abs(Q(z)+Q(1, 3))))
                            for z in range(-int(radius)-2, int(radius)+3)), Q(0))
    check(mass == roof, ("Bernoulli/roof", p, q, mass, roof))
    return mass, 6*number


@lru_cache(None)
def base_sheet(speed, sheet):
    """Integer endpoints on the circle cut at zero, denominator 42*speed."""
    denominator = 42*speed
    pieces = []
    for owner in range(speed):
        center = 14*((3*owner-speed*sheet) % (3*speed))
        left, right = center-3, center+3
        if left < 0:
            pieces.extend(((0, right), (denominator+left, denominator)))
        elif right > denominator:
            pieces.extend(((left, denominator), (0, right-denominator)))
        else:
            pieces.append((left, right))
    return tuple(sorted(pieces))


def intersect(a, b):
    i = j = 0
    result = []
    while i < len(a) and j < len(b):
        left, right = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if left < right:
            result.append((left, right))
        ai, bj = a[i][1], b[j][1]
        i += ai <= bj
        j += bj <= ai
    return result


def contacts(a, b):
    i = j = 0
    result = []
    while i < len(a) and j < len(b):
        left, right = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if left < right:
            result.append((i, j, right-left))
        ai, bj = a[i][1], b[j][1]
        i += ai <= bj
        j += bj <= ai
    return result


def triple(w):
    denominator = 42*w[0]*w[1]*w[2]
    sheets = {(i, s): tuple((left*(denominator//(42*w[i])), right*(denominator//(42*w[i])))
                            for left, right in base_sheet(w[i], s))
              for i in range(3) for s in range(3)}
    bounds = []
    physical_controls = []
    for i, j in combinations(range(3), 2):
        k = 3-i-j
        capacity = physical = pair_mass = pair_count = 0
        for assignment in permutations(range(3)):
            a = intersect(sheets[i, assignment[i]], sheets[j, assignment[j]])
            b = sheets[k, assignment[k]]
            pair_count += len(a)
            pair_mass += sum(right-left for left, right in a)
            edges = contacts(a, b)
            loads_a, loads_b = [0]*len(a), [0]*len(b)
            degrees_a, degrees_b = [0]*len(a), [0]*len(b)
            for x, y, overlap in edges:
                load = min(a[x][1]-a[x][0], b[y][1]-b[y][0])
                loads_a[x] += load
                loads_b[y] += load
                degrees_a[x] += 1
                degrees_b[y] += 1
                capacity += load
                physical += overlap
            # These loads attain every individual edge upper bound.  Their
            # feasibility is an exact primal + matching upper certificate.
            check(all(load <= right-left for load, (left, right) in zip(loads_a, a)),
                  ("pair capacities", w, (i, j), assignment))
            check(all(load <= right-left for load, (left, right) in zip(loads_b, b)),
                  ("third capacities", w, (i, j), assignment))
            check(all(degrees_a[x] == 1 or degrees_b[y] == 1 for x, y, _ in edges),
                  ("star forest", w, (i, j), assignment))
            if w[k] >= max(w[i], w[j]) and w[k] < 6*max(w[i], w[j]):
                check(all(d <= 1 for d in degrees_a), ("single contact", w, (i, j)))
        common = gcd(w[i], w[j])
        pm, pn = pair_formula(w[i]//common, w[j]//common)
        check(Q(pair_mass, denominator) == pm, ("physical pair mass", w, (i, j)))
        check(pair_count == common*pn, ("physical pair count", w, (i, j)))
        check(7*capacity*w[k] <= pair_mass*w[k] + Q(8, 7)*pair_count*denominator,
              ("all-height counting envelope", w, (i, j)))
        bounds.append(Q(capacity, denominator))
        physical_controls.append(Q(physical, denominator))
    check(len(set(physical_controls)) == 1, ("independent pair grouping", w))
    check(physical_controls[0] <= min(bounds), ("upper bound", w))
    if 13*w[2] >= 44*w[1]:
        check(bounds[0] < TARGET, ("44/13 tail certificate", w, bounds[0]))
    return tuple(bounds), physical_controls[0]


def main():
    parser = ArgumentParser()
    parser.add_argument("--height", type=int, default=79)
    args = parser.parse_args()
    # Finite residue proof obligations in the analytic argument.
    corrections = [correction(p, q) for p in range(1, 14, 2) for q in range(1, 14, 2)]
    check(max(corrections) == Q(26, 49), "49-residue correction maximum")
    for t in range(7):
        r = Q(t, 7)
        n = (r+Q(1, 3)).__floor__() + (r-Q(1, 3)).__floor__() + 1
        check(n - 2*r <= Q(4, 7), ("seven count residues", t))
    short_pairs = [(p, q) for p in range(1, 16, 2) for q in range(p+2, 16, 2)
                   if p % 3 and q % 3 and gcd(p, q) == 1 and p*q < 16]
    check(short_pairs == [(1, 5), (1, 7), (1, 11), (1, 13)], "small product universe")
    for p, q in short_pairs:
        mass, _ = pair_formula(p, q)
        check(mass <= Q(12, 77), ("small-product mass", p, q))
    for q in (5, 7, 11):
        for p in range(1, q, 2):
            if p % 3 and gcd(p, q) == 1:
                _, number = pair_formula(p, q)
                check(Q(number, q) < Q(24, 13), ("small count denominator", p, q))
    check(Q(6, 49)+Q(26, 49*16) < Q(12, 77), "large product mass bound")
    check(Q(12, 7)+Q(12, 7*13) == Q(24, 13), "large denominator count bound")
    check(Q(12, 539)+Q(32, 637) < TARGET, "c>=6b bound")
    check(Q(24, 91)/Q(44, 13) == TARGET, "single-contact threshold")
    print("ANALYTIC finite obligations: 49 Bernoulli residues; 7 count residues; 4 small mass pairs; q=5,7,11 count bases")
    print("PAIR SHARP BOUNDS: T<=12/77 at p:q=1:11; N/b<=24/13 at p:q=11:13")
    print("ALL-HEIGHT: star-forest capacity=edge-min sum; c/b>=44/13 implies U_(a,b)<6/77 for odd ternary units")

    controls = ((1, 5, 11), (1, 19, 79), (11, 13, 17), (11, 13, 43),
                (11, 13, 47), (11, 13, 79), (1, 11, 331), (7, 11, 997))
    for w in controls:
        bound, mass = triple(w)
        print("CONTROL", w, "pair bounds", tuple(map(str, bound)), "physical", str(mass))
    check(triple((1, 5, 11))[1] == TARGET, "equality hostile")
    check(triple((1, 19, 79))[1] < min(triple((1, 19, 79))[0]), "strict information-loss hostile")

    speeds = [v for v in range(1, args.height+1, 2) if v % 3]
    rows = tails = strict = equal = 0
    digest = sha256()
    best = (Q(0), None)
    for w in combinations(speeds, 3):
        if gcd(gcd(w[0], w[1]), w[2]) != 1:
            continue
        bounds, mass = triple(w)
        selected = min(bounds)
        rows += 1
        tails += 13*w[2] >= 44*w[1]
        strict += selected < TARGET
        equal += selected == TARGET
        if selected > best[0]:
            best = selected, w
        digest.update(repr((w, bounds, mass)).encode())
    print("FINITE-EXACT UNIVERSE: distinct positive odd ternary units, primitive, max<=", args.height)
    print("ROWS", rows, "TAIL ROWS", tails, "NETWORK BELOW/EQUAL/ABOVE", strict, equal, rows-strict-equal)
    print("MAX SELECTED", str(best[0]), best[1])
    print("SEMANTIC SHA256", digest.hexdigest())
    print("CHECKS", CHECKS)
    print("OPEN: c/b<44/13 at unbounded height; entry; synchronization; LRC(14)")


if __name__ == "__main__":
    main()
