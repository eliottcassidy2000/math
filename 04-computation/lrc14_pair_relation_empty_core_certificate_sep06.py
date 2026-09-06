#!/usr/bin/env python3
"""Actual-zero-coordinate addendum to the independently audited norm20 proof.

Imports the frozen core's exact carrier/projection routines as declared
dependencies; every head also has an independent literal-sheet comparison.
The pair interval width is computed directly, without polytope enumeration.
"""

from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from math import gcd

from lrc14_empty_core_certificate_sep06 import raw_carriers, projection_data, primitive, cross
from lrc14_one_ray_overnight_hexagon_sep05 import literal_projection_data

CHECKS = 0
R = Q(3, 14)
TARGET = Q(6, 77)


def need(test, detail):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(detail)


def pair_record(p, q):
    bound = (3*(p+q)-1)//14
    defects = tuple(d for d in range(-bound, bound+1) if d % 3)
    lengths = tuple(min(R, (p*R-d)/q)-max(-R, (-p*R-d)/q)
                    for d in defects)
    need(all(length > 0 for length in lengths), ('strict nonempty slices', p, q))
    slope = sum(2*R/q+length/p for length in lengths)/3
    intercept = len(defects)
    need(slope < Q(2, 11), ('strict count slope', p, q, slope))
    cutoff = intercept/(Q(2, 11)-slope)
    return (p, q), defects, lengths, slope, intercept, cutoff


def pair_head(record):
    (p, q), defects, lengths, slope, intercept, cutoff = record
    limit = (cutoff.numerator-1)//cutoff.denominator
    rows = {}
    for s in range(1, limit//q+1, 2):
        if s % 3 == 0:
            continue
        for t in range(1, limit+1, 2):
            if t % 3 == 0 or t in (p*s, q*s) or gcd(s, t) != 1:
                continue
            w = tuple(sorted((p*s, q*s, t)))
            v = tuple(p if wi == q*s else -q if wi == p*s else 0 for wi in w)
            need(sum(vi*wi for vi, wi in zip(v, w)) == 0, ('pair relation', w, v))
            rows.setdefault(w, []).append((v, s, t))
    return rows


def main():
    patterns = tuple((p, q) for p in range(1, 20, 2)
                     for q in range(p+2, 21-p, 2)
                     if p*q % 3 and gcd(p, q) == 1)
    need(len(patterns) == 11, ('complete support-two pattern list', patterns))
    records = [pair_record(p, q) for p, q in patterns]
    need(max(r[-1] for r in records) == Q(21021, 263), 'maximum finite cutoff')
    memberships = {}
    for record in records:
        rows = pair_head(record)
        for w, presentations in rows.items():
            memberships.setdefault(w, []).append((record, presentations))
        print('pair', record[0], 'defects', record[1], 'error_lengths', record[2],
              'slope', record[3], 'intercept', record[4], 'cutoff', record[5],
              'head_rows', len(rows))
    kinds = Counter()
    equalities = set()
    digest = sha256()
    memberships_count = sum(len(rows) for rows in memberships.values())
    for w in sorted(memberships, key=lambda w:(w[2], w[0], w[1])):
        live = raw_carriers(w)
        projected, physical = projection_data(w, live)
        need((projected, physical) == literal_projection_data(w), ('literal comparison', w))
        need(min(projected) <= TARGET, ('selected finite target', w, projected))
        if min(projected) == TARGET:
            equalities.add(w)
        directions = {primitive(C) for C in live}
        norm4 = any(sum(map(abs, v)) == 4 for v in directions)
        kinds['empty' if not directions else 'norm4' if norm4 else
              'other_one' if len(directions) == 1 else 'multi'] += 1
        if not norm4:
            need(max(projected) < TARGET, ('every projection outside norm4', w, projected))
        for record, presentations in memberships[w]:
            (p, q), defects, lengths, slope, intercept, cutoff = record
            need(len(live) < slope*w[2]+intercept, ('compiled pair count', w, p, q))
            for v, s, t in presentations:
                counts = Counter()
                for C in live:
                    z = cross(v, C)
                    delta = Q(z[0], w[0])
                    need(delta.denominator == 1 and z == tuple(delta*wi for wi in w),
                         ('integer pair defect', w, v, C))
                    need(int(delta) in defects, ('allowed pair defect', w, v, C, delta))
                    counts[int(delta)] += 1
                for delta, length in zip(defects, lengths):
                    width = 2*R*s+t*length/p
                    need(counts[delta] < width/3+1,
                         ('strict pair affine interval count', w, v, delta, width))
        digest.update(repr((w, tuple(sorted(live)), projected, physical)).encode())
    need(equalities == {(1, 5, 11)}, ('sole pair-head target equality', equalities))
    print('patterns', len(patterns), 'head_memberships', memberships_count)
    print('unique_head_rows', len(memberships), 'support_classes', dict(sorted(kinds.items())))
    print('head_equalities', tuple(sorted(equalities)))
    print('complete_head_digest', digest.hexdigest())
    print('explicit_checks', CHECKS)
    print('PASS: all support-two primitive relation shells through l1=20 satisfy the network target')


if __name__ == '__main__':
    main()
