#!/usr/bin/env python3
"""Exact producer audit for the two THM-4391 nonunit l1=16 LRC14 sectors.

The script uses only the Python standard library.  Every check stays live under
optimized Python.  It proves the finite part of the sharp atlas and compares
the quotient roof with a definition-level circle-cell implementation.
"""

from fractions import Fraction as F
from itertools import combinations_with_replacement
from math import gcd


R = F(3, 14)
PATTERNS = ((7, 7, 2), (5, 7, 4))  # h,m,l in m*b=s*h*p+t*l*q
EXPECTED = {
    (7, 7, 2): {
        "maximum": F(304, 12397),
        "winner": (1, 23, 77),
        "parts": {-3: F(31, 12397), 0: F(22, 1127), 3: F(31, 12397)},
        "cutoff": 366,
        "presentations": 1754,
        "triples": 877,
    },
    (5, 7, 4): {
        "maximum": F(2178, 91945),
        "winner": (5, 37, 71),
        "parts": {-3: F(2, 2485), 0: F(58, 2627), 3: F(2, 2485)},
        "cutoff": 416,
        "presentations": 2389,
        "triples": 2388,
    },
}


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def egcd(a, b):
    if b == 0:
        return a, 1, 0
    g, x, y = egcd(b, a % b)
    return g, y, x - (a // b) * y


def determinant_section(p, b, k):
    """Return n_p,n_b with b*n_p-p*n_b=k; gcd(p,b)=1."""
    g, u, v = egcd(b, p)
    check(g == 1, ("noncyclic-pb", p, b, g))
    n_p, n_b = u * k, -v * k
    check(b * n_p - p * n_b == k, ("bad-section", p, b, k))
    return n_p, n_b


def component_length(p, b, q, h, m, ell, s, t, delta, k):
    numer_pq = m * k + p * delta
    numer_bq = s * h * k + b * delta
    check(numer_pq % ell == 0, ("pq-nonintegral", numer_pq, ell))
    check(numer_bq % ell == 0, ("bq-nonintegral", numer_bq, ell))
    terms = (
        2 * R / p,
        2 * R / b,
        2 * R / q,
        R / p + R / b - F(abs(k), p * b),
        R / p + R / q - F(abs(numer_pq), ell * p * q),
        R / b + R / q - F(abs(numer_bq), ell * b * q),
    )
    return max(F(0), min(terms))


def quotient_measure(row):
    p, b, q, h, m, ell, s, t = row
    bound = R * (p + b)
    limit = bound.numerator // bound.denominator + 1
    parts = {}
    terms = {}
    for delta in (-3, 0, 3):
        subtotal = F(0)
        listed = []
        residue = None
        for k in range(-limit, limit + 1):
            if k % 3 == 0:
                continue
            n_p, n_b = determinant_section(p, b, k)
            numerator = m * n_b - s * h * n_p - delta
            if numerator % (t * ell):
                continue
            n_q = numerator // (t * ell)
            check(m * n_b - s * h * n_p - t * ell * n_q == delta,
                  ("bad-defect", row, delta, k))
            check((m * k + p * delta) % ell == 0,
                  ("bad-pq-integrality", row, delta, k))
            check((s * h * k + b * delta) % ell == 0,
                  ("bad-bq-integrality", row, delta, k))
            if residue is None:
                residue = k % ell
            check(k % ell == residue, ("multiple-ell-residues", row, delta, k))
            value = component_length(p, b, q, h, m, ell, s, t, delta, k)
            if value:
                listed.append((k, value))
                subtotal += value
        parts[delta] = subtotal
        terms[delta] = tuple(listed)
        check(subtotal == sum((v for _, v in listed), F(0)), "subtotal mismatch")
    check(parts[-3] == parts[3], ("sign-symmetry", row, parts))
    return sum(parts.values(), F(0)), parts, terms


def nearest_integer(x):
    # Called only at open-cell midpoints, never at half-integers.
    n = x.numerator // x.denominator
    return n + (2 * (x - n) > 1)


def physical_measure(speeds):
    """Definition-level y-circle comb measure, independent of quotient formula."""
    walls = {F(0), F(1)}
    for w in speeds:
        for n in range(0, w + 1):
            for sign in (-1, 1):
                wall = F(n, w) + sign * R / w
                if 0 < wall < 1:
                    walls.add(wall)
    walls = sorted(walls)
    total = F(0)
    live_cells = 0
    for lo, hi in zip(walls, walls[1:]):
        y = (lo + hi) / 2
        owners = []
        eligible = True
        for w in speeds:
            n = nearest_integer(w * y)
            error = w * y - n
            if not abs(error) < R:
                eligible = False
                break
            owners.append((-pow(w, -1, 3) * n) % 3)
        if eligible and len(set(owners)) == 3:
            total += hi - lo
            live_cells += 1
    return total, live_cells


def presentations(h, m, ell, cutoff):
    units = tuple(x for x in range(1, cutoff, 2) if x % 3)
    rows = []
    for p in units:
        for b in units:
            if p == b:
                continue
            for s in (-1, 1):
                for t in (-1, 1):
                    numerator = m * b - s * h * p
                    if numerator % (t * ell):
                        continue
                    q = numerator // (t * ell)
                    if not (0 < q < cutoff and q % 2 and q % 3):
                        continue
                    if q in (p, b) or gcd(gcd(p, b), q) != 1:
                        continue
                    check(gcd(p, b) == 1,
                          ("even-ell-did-not-force-coprime-pair", p, b, q, ell))
                    rows.append((p, b, q, h, m, ell, s, t))
    return tuple(rows)


def integral_affine(lo, hi, slope, intercept):
    return slope * (hi * hi - lo * lo) / 2 + intercept * (hi - lo)


def strip_area(h, m, ell, delta):
    """Area in the (e_p,e_b) square where |delta-sh e_p+m e_b|<ell*r."""
    check(h <= m, ("area-order", h, m))
    plateau = (m - h) * R
    support = (m + h) * R
    lo = delta - ell * R
    hi = delta + ell * R
    cuts = sorted({lo, hi, -support, -plateau, F(0), plateau, support})
    area_uv = F(0)
    for left, right in zip(cuts, cuts[1:]):
        a, b = max(left, lo), min(right, hi)
        if not a < b:
            continue
        mid = (a + b) / 2
        if abs(mid) >= support:
            continue
        if abs(mid) <= plateau:
            area_uv += 2 * h * R * (b - a)
        elif mid > 0:
            area_uv += integral_affine(a, b, F(-1), support)
        else:
            area_uv += integral_affine(a, b, F(1), support)
    return area_uv / (h * m)


def width_cutoff(target, bulk):
    x = F(18, 7) / (target - bulk)
    return x.numerator // x.denominator + 1


def coefficient_shell():
    rows = []
    for a, b, c in combinations_with_replacement(range(1, 15), 3):
        if a + b + c != 16 or any(x % 3 == 0 for x in (a, b, c)):
            continue
        if gcd(gcd(a, b), c) != 1:
            continue
        rows.append((a, b, c))
    return tuple(rows)


def short_owner_degenerate_relations(winner):
    rows = []
    for c0 in range(-12, 13):
        for c1 in range(-12, 13):
            for c2 in range(-12, 13):
                c = (c0, c1, c2)
                if 0 in c or sum(abs(x) for x in c) > 12:
                    continue
                if sum(x * y for x, y in zip(c, winner)):
                    continue
                if gcd(gcd(abs(c0), abs(c1)), abs(c2)) != 1:
                    continue
                first = next(x for x in c if x)
                if first < 0:
                    continue
                rows.append(c)
    return tuple(sorted(set(rows), key=lambda c: (sum(map(abs, c)), c)))


def main():
    expected_shell = (
        (1, 1, 14), (1, 2, 13), (1, 4, 11), (1, 5, 10),
        (1, 7, 8), (2, 7, 7), (4, 5, 7),
    )
    check(coefficient_shell() == expected_shell,
          ("coefficient-shell", coefficient_shell()))
    print(f"L1_16_PRIMITIVE_FULL_SUPPORT_PATTERNS={coefficient_shell()}")
    print("COEFFICIENT_ONE_FIVE plus EXTRA=(2,7,7),(4,5,7)")

    expected_areas = {
        (7, 7, 2): {-3: F(9, 4802), 0: F(117, 2401), 3: F(9, 4802)},
        (5, 7, 4): {-3: F(9, 3430), 0: F(171, 1715), 3: F(9, 3430)},
    }
    check_count = 2
    for pattern in PATTERNS:
        h, m, ell = pattern
        areas = {delta: strip_area(h, m, ell, delta) for delta in (-3, 0, 3)}
        check(areas == expected_areas[pattern], ("strip-areas", pattern, areas))
        bulk = F(2, 3 * ell) * sum(areas.values(), F(0))
        check(bulk == F(6, 343), ("bulk", pattern, bulk))
        expected = EXPECTED[pattern]
        cutoff = width_cutoff(expected["maximum"], bulk)
        check(cutoff == expected["cutoff"], ("cutoff", pattern, cutoff))
        check(bulk + F(18, 7 * cutoff) < expected["maximum"],
              ("tail-does-not-close", pattern))
        check(bulk + F(18, 7 * (cutoff - 1)) >= expected["maximum"],
              ("cutoff-not-minimal", pattern))

        rows = presentations(h, m, ell, cutoff)
        triples = {tuple(sorted(row[:3])) for row in rows}
        check(len(rows) == expected["presentations"],
              ("presentation-count", pattern, len(rows)))
        check(len(triples) == expected["triples"],
              ("triple-count", pattern, len(triples)))
        physical_cache = {}
        by_triple = {}
        maxima = []
        for row in rows:
            value, parts, terms = quotient_measure(row)
            triple = tuple(sorted(row[:3]))
            if triple not in physical_cache:
                physical_cache[triple] = physical_measure(triple)[0]
            check(value == physical_cache[triple],
                  ("physical-disagreement", row, value, physical_cache[triple]))
            if triple in by_triple:
                check(by_triple[triple] == value,
                      ("chart-disagreement", triple, by_triple[triple], value))
            by_triple[triple] = value
            if value == expected["maximum"]:
                maxima.append((triple, row, parts, terms))
            check(value <= expected["maximum"],
                  ("larger-finite-row", pattern, triple, value))
            check_count += 10 + sum(len(v) for v in terms.values())
        winner_rows = [x for x in maxima if x[0] == expected["winner"]]
        check(maxima and len({x[0] for x in maxima}) == 1,
              ("nonunique-winner", pattern, maxima))
        check(winner_rows[0][2] == expected["parts"],
              ("winner-parts", pattern, winner_rows[0][2]))
        print(
            f"PATTERN={pattern} maximum={expected['maximum']} "
            f"winner={expected['winner']} parts={expected['parts']} "
            f"bulk={bulk} cutoff={cutoff} "
            f"presentations={len(rows)} triples={len(triples)} "
            f"maximizing_presentations={len(maxima)}"
        )

    hostile_relations = {
        (1, 23, 77): short_owner_degenerate_relations((1, 23, 77)),
        (5, 37, 71): short_owner_degenerate_relations((5, 37, 71)),
    }
    check((-8, -3, 1) in tuple(tuple(-x for x in c) for c in hostile_relations[(1, 23, 77)])
          or (8, 3, -1) in hostile_relations[(1, 23, 77)],
          ("missing-277-short-relation", hostile_relations[(1, 23, 77)]))
    check((8, -3, 1) in hostile_relations[(5, 37, 71)],
          ("missing-457-short-relation", hostile_relations[(5, 37, 71)]))
    for winner, relations in hostile_relations.items():
        check(all(any(x % 3 == 0 for x in c) for c in relations),
              ("unexpected-owner-admissible-short-relation", winner, relations))
        print(f"WINNER_SHORT_OWNER_DEGENERATE_RELATIONS {winner} {relations}")

    print(f"CHECKS_AT_LEAST={check_count}")
    print("TAIL_ENVELOPE=6/343+18/(7W)")
    print("STATUS=PASS; these two sectors only; universal comb bound and LRC14 OPEN")


if __name__ == "__main__":
    main()
