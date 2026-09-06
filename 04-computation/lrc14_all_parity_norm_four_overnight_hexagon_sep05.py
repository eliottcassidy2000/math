#!/usr/bin/env python3
"""Complete parity-expanded signed112 endpoint head, pq<81.

Full triples are primitive; coefficient-one endpoints need not be coprime.
The ten even-endpoint cases are essential. No primary norm-four producer
is imported: raw/literal engines and the inherited periodic primitive are
used as separate checks of the endpoint formula.
"""

from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from math import gcd
from pathlib import Path


SPEC = spec_from_file_location("boundary", Path(__file__).with_name(
    "lrc14_sum_ray_boundary_discrepancy_overnight_hexagon_sep05.py"))
BOUNDARY = module_from_spec(SPEC)
SPEC.loader.exec_module(BOUNDARY)
BASE = BOUNDARY.BASE
CHECKS = 0
TARGET = Q(6, 77)
SHARP = Q(11, 140)


def need(value, payload):
    global CHECKS
    CHECKS += 1
    if not value:
        raise AssertionError(payload)


def audit(p, q, sign):
    need(p < q and (q + sign * p) % 2 == 0, "integer endpoint branch")
    r = (q + sign * p) // 2
    coefficients = (1, -2, 1) if sign == 1 else (1, 2, -1)
    ordered = sorted(zip((p, r, q), coefficients))
    w = tuple(t[0] for t in ordered)
    d = tuple(t[1] for t in ordered)
    need(len(set(w)) == 3 and all(t % 3 for t in w)
         and gcd(gcd(w[0], w[1]), w[2]) == 1, ("full primitive triple", w))
    need(gcd(p, q) in (1, 2), ("endpoint gcd", w, p, q))
    if gcd(p, q) == 2:
        need(r % 2 == 1, ("even endpoint sidecar", w))
    A, B = Q(3 * (q - p), 28), Q(3 * (q + p), 28)
    K = (B.numerator - 1) // B.denominator
    expected = {tuple(s * k * t for t in d)
                for k in range(1, K + 1) if k % 3 for s in (-1, 1)}
    live = BASE.carriers(w)
    need(live == expected, ("complete full-lattice ray", w, live ^ expected))
    projections, physical = BASE.projection_data(w, live)
    need(BASE.literal_projection_data(w) == (projections, physical),
         ("literal six-sheet network", w))
    endpoint_formula = Q(3, 49) + Q(4, p * q) * (BOUNDARY.primitive(B) - BOUNDARY.primitive(A))
    need(physical == endpoint_formula, ("parity-free endpoint quadrature", w))
    need(projections[w.index(r)] == physical == min(projections),
         ("coefficient-two selector is physical", w, projections, physical))
    need(physical <= Q(3, 49) + Q(4, 3 * p * q), ("endpoint product error", w))
    return w, projections, physical


def main():
    need(Q(3, 49) + Q(4, 3 * 81) == Q(925, 11907) < TARGET,
         "infinite product tail")
    rows = []
    endpoint_counts, parity_counts = Counter(), Counter()
    digest = sha256()
    for p in range(1, 9):
        for q in range(p + 1, 81 // p + 1):
            if p * q >= 81 or (q - p) % 2:
                continue
            for sign in (-1, 1):
                r = (q + sign * p) // 2
                w = tuple(sorted((p, r, q)))
                if (len(set(w)) != 3 or not all(t % 3 for t in w)
                        or gcd(gcd(w[0], w[1]), w[2]) != 1):
                    continue
                w, projections, physical = audit(p, q, sign)
                need(physical <= SHARP, ("sharp complete head", w))
                if w != (2, 11, 20):
                    need(physical <= TARGET, ("only old-target exception", w))
                record = (p, q, r, sign, w, projections, physical)
                rows.append(record)
                endpoint_counts[gcd(p, q)] += 1
                parity_counts[sum(t % 2 == 0 for t in w)] += 1
                digest.update((repr(record) + "\n").encode())
    need(len(rows) == len({row[4] for row in rows}) == 40, "complete distinct endpoint head")
    need(endpoint_counts == {1: 30, 2: 10}, ("both endpoint gcd types", endpoint_counts))
    need(parity_counts == {0: 13, 1: 17, 2: 10}, ("complete parity head", parity_counts))
    sharp = [row[4] for row in rows if row[6] == SHARP]
    old_equal = [row[4] for row in rows if row[6] == TARGET]
    need(sharp == [(2, 11, 20)] and old_equal == [(1, 5, 11)],
         ("exact equality loci", sharp, old_equal))
    print("COMPLETE ALL-PARITY ENDPOINT HEAD pq<81", len(rows))
    print("ENDPOINT GCD", dict(sorted(endpoint_counts.items())),
          "NUMBER EVEN COORDINATES", dict(sorted(parity_counts.items())))
    print("SEMANTIC SHA256", digest.hexdigest())
    print("SHARP min-projection AND actual", SHARP, "EQUALITY", sharp)
    print("OLD6/77 EQUALITY", old_equal, "SOLE ABOVE-TARGET EXCEPTION", sharp)
    print("TAIL AT81", Q(925, 11907), "MARGIN TO6/77", TARGET - Q(925, 11907))
    for p, q, sign in ((1, 11, -1), (2, 20, 1), (1, 121, 1),
                       (2, 124, -1), (4, 202, 1), (1, 101, -1)):
        w, projections, physical = audit(p, q, sign)
        if p * q >= 81:
            need(physical < TARGET, ("large product control", w))
        print("CONTROL", (p, q, sign), w, "E", projections, "ACTUAL", physical)
    original = BASE.literal_projection_data((2, 11, 20))
    for scale in (2, 5):
        need(BASE.literal_projection_data(tuple(scale * t for t in (2, 11, 20))) == original,
             ("common unit dilation", scale))
    print("CHECKS", CHECKS, "INDEPENDENT RAW/LITERAL CHECKS", BASE.CHECKS)
    print("PROVED all-parity signed112; min E=actual; unique11/140 maximum and one6/77 exception")
    print("OPEN actual arbitrary-body entry/synchronization; no body Haar floor supplied")


if __name__ == "__main__":
    main()
