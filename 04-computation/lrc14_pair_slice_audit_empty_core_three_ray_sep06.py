#!/usr/bin/env python3
"""Independent bounded audit of actual-zero-coordinate norm20 addendum.

Uses the prior independent auditor's rational polygon clipper and direct
integer boxes, not the pair producer's rectangle or speed parameterization.
Enumerates the entire eligible speed universe through height79 before filtering
by gcd-reduced pair ratios. No change to either audited producer is required.
"""

import ast
import re
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations, permutations
from math import gcd
from pathlib import Path

from lrc14_relation_slice_audit_empty_core_three_ray_sep06 import sliced_polygon, full_box, projections

CHECKS = 0


def need(test, detail):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(detail)


def polygon_slope(p, q, defects):
    best = Q(0)
    for raw in set(permutations((p, -q, 0))):
        # Rotate to a nonzero elimination coordinate; the cube and width
        # consumer are invariant under this simultaneous permutation.
        pivot = next(i for i in range(3) if raw[i])
        v = raw[pivot:]+raw[:pivot]
        width_coordinate = next(i for i in range(1, 3) if v[i])
        polygons = [sliced_polygon(v, Q(d), -Q(3, 14), Q(3, 14)) for d in defects]
        for w in sliced_polygon(v, Q(0), Q(0), Q(1)):
            total = Q(0)
            for polygon in polygons:
                scalars = []
                for e in polygon:
                    cross = (w[1]*e[2]-w[2]*e[1], w[2]*e[0]-w[0]*e[2],
                             w[0]*e[1]-w[1]*e[0])
                    scalars.append(cross[width_coordinate]/v[width_coordinate])
                total += (max(scalars)-min(scalars))/3
            best = max(best, total)
    return best


def main():
    path = Path(__file__).parents[1]/'05-knowledge/results/lrc14_pair_relation_empty_core_certificate_sep06.out'
    text = path.read_text()
    records = {}
    for line in text.splitlines():
        found = re.match(r'pair (\(.*?\)) defects (\(.*?\)) error_lengths .* slope (\S+) intercept (\S+) cutoff (\S+) head_rows (\d+)', line)
        if found:
            p, ds, slope, intercept, cutoff, count = found.groups()
            records[ast.literal_eval(p)] = (ast.literal_eval(ds), Q(slope), int(intercept), Q(cutoff), int(count))
    patterns = {(p, q) for p in range(1, 20) for q in range(p+1, 20)
                if p+q <= 20 and gcd(p, q) == 1 and p % 2 and q % 2 and p*q % 3}
    need(len(patterns) == 11 and patterns == set(records), ('complete coefficient pattern universe', patterns))
    R = Q(3, 14)
    for p, q in sorted(patterns):
        defects = tuple(d for d in range(-5, 6) if 14*abs(d) < 3*(p+q) and d % 3)
        closed_slope = (2*R*p*len(defects)+sum(min(2*p*R, R*(p+q)-abs(d)) for d in defects))/(3*p*q)
        poly_slope = polygon_slope(p, q, defects)
        need(closed_slope == poly_slope, ('closed formula and clipped polygon agree', p, q))
        cutoff = len(defects)/(Q(2, 11)-closed_slope)
        need((defects, closed_slope, len(defects), cutoff) == records[p, q][:4],
             ('complete exact compiler equality', p, q, closed_slope, cutoff))
    need(max(record[3] for record in records.values()) == Q(21021, 263), 'largest cutoff')
    rows = {}
    members = Counter()
    universe = 0
    for w in combinations([x for x in range(1, 80, 2) if x % 3], 3):
        if gcd(gcd(w[0], w[1]), w[2]) != 1:
            continue
        universe += 1
        found = set()
        for a, b in combinations(w, 2):
            d = gcd(a, b)
            pattern = a//d, b//d
            if pattern in records and w[2] < records[pattern][3]:
                found.add(pattern)
        if found:
            rows[w] = found
            members.update(found)
    need(universe == 2910, ('complete H79 speed universe', universe))
    need(all(members[p] == record[4] for p, record in records.items()),
         ('direct gcd-ratio finite heads', dict(members)))
    digest = sha256()
    equalities = []
    for w in sorted(rows, key=lambda w: (w[2], w[0], w[1])):
        live = full_box(w)
        values, physical = projections(w, live)
        need(min(values) <= Q(6, 77), ('actual finite projection consumer', w, values))
        if min(values) == Q(6, 77):
            equalities.append(w)
        digest.update(repr((w, tuple(sorted(live)), values, physical)).encode())
    expected = re.search(r'complete_head_digest (\w+)', text).group(1)
    need(digest.hexdigest() == expected, ('complete consumer digest', digest.hexdigest()))
    need(equalities == [(1, 5, 11)], ('unique equality', equalities))
    print('PASS independent actual-zero-coordinate relation addendum audit')
    print('PATTERNS 11; all rectangle slopes independently verified by clipped polygons and closed formula')
    print('COMPLETE SPEED UNIVERSE', universe, 'height<=79')
    print('HEAD MEMBERSHIPS', sum(members.values()), 'UNIQUE HEAD ROWS', len(rows))
    print('EXACT CONSUMERS: all raw integer-box carriers, all three projections, physical mass')
    print('UNIQUE EQUALITY', equalities)
    print('COMPLETE HEAD SEMANTIC SHA256', digest.hexdigest())
    print('CHECKS', CHECKS)


if __name__ == '__main__':
    main()
