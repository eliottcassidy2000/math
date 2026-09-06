#!/usr/bin/env python3
"""Independent audit of the concurrent short-relation slice certificate.

Uses rational 2D polygon clipping instead of cube-edge vertex enumeration,
the third cross-product coordinate instead of the first, direct height-first
triple enumeration instead of solving speed relations, and full integer-box
carriers. Producer output is read only to compare frozen exact records.
"""

import ast
import re
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations, permutations, product
from math import gcd
from pathlib import Path

CHECKS = 0


def need(test, detail):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(detail)


def clip(poly, A, B, C):
    """Keep A*x+B*y<=C, retaining exact lower-dimensional boundaries."""
    if not poly:
        return []
    out = []
    for p, q in zip(poly, poly[1:]+poly[:1]):
        fp = A*p[0]+B*p[1]-C
        fq = A*q[0]+B*q[1]-C
        if fp <= 0:
            out.append(p)
        if (fp < 0 < fq) or (fq < 0 < fp):
            t = fp/(fp-fq)
            out.append((p[0]+t*(q[0]-p[0]), p[1]+t*(q[1]-p[1])))
    return list(dict.fromkeys(out))


def sliced_polygon(v, delta, lo, hi):
    poly = [(lo, lo), (hi, lo), (hi, hi), (lo, hi)]
    lower, upper = sorted((v[0]*lo, v[0]*hi))
    poly = clip(poly, -v[1], -v[2], upper-delta)
    poly = clip(poly, v[1], v[2], delta-lower)
    return tuple(((delta-v[1]*y-v[2]*z)/v[0], y, z) for y, z in poly)


def sectors(pattern):
    return {tuple(-a[i] if i == j else a[i] for i in range(3))
            for a in permutations(pattern) for j in range(3)}


def compile_independently(pattern):
    unit = all(x % 3 for x in pattern)
    deltas = tuple(d for d in range(-5, 6)
                   if 14*abs(d) < 3*sum(pattern) and (d % 3 == 0) == unit)
    factor = Q(2, 3) if unit else Q(1, 3)
    intercept = (Q(4, 3) if unit else Q(1))*len(deltas)
    best = Q(0)
    for v in sectors(pattern):
        errors = [sliced_polygon(v, Q(d), -Q(3, 14), Q(3, 14)) for d in deltas]
        for w in sliced_polygon(v, Q(0), Q(0), Q(1)):
            value = Q(0)
            for points in errors:
                coordinates = [(w[0]*e[1]-w[1]*e[0])/v[2] for e in points]
                value += factor*(max(coordinates)-min(coordinates))
            best = max(best, value)
    cutoff = None if best >= Q(2, 11) else intercept/(Q(2, 11)-best)
    return deltas, best, intercept, cutoff


def cross(u, v):
    return (u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0])


def residue_audit():
    for v in product(range(3), repeat=3):
        if sum(x == 0 for x in v) > 1:
            continue
        unit = all(v)
        for w in product((1, 2), repeat=3):
            if sum(v[i]*w[i] for i in range(3)) % 3:
                continue
            for delta in range(3):
                solutions = [C for C in product(range(3), repeat=3)
                             if sum(C[i]*w[i] for i in range(3)) % 3 == 0
                             and all((cross(v, C)[i]-delta*w[i]) % 3 == 0 for i in range(3))]
                need(len(solutions) == 3, ('affine F3 fiber cardinality', v, w, delta))
                live = sum(all(C) for C in solutions)
                expected = (2 if delta == 0 else 0) if unit else (0 if delta == 0 else 1)
                need(live == expected, ('exact owner residue fiber', v, w, delta, live))


def full_box(w):
    a, b, c = w
    bounds = [(3*(sum(w)-v)-1)//14 for v in w]
    out = set()
    for x in range(-bounds[0], bounds[0]+1):
        for y in range(-bounds[1], bounds[1]+1):
            if not x % 3 or not y % 3 or (a*x+b*y) % c:
                continue
            z = -(a*x+b*y)//c
            if z % 3 and abs(z) <= bounds[2]:
                out.add((x, y, z))
    return out


def projections(w, live):
    denominator = 14*w[0]*w[1]*w[2]
    cap = 6*w[0]*w[1]
    totals, physical = [0, 0, 0], 0
    for C in live:
        terms = [min(cap, w[i]*(3*(sum(w)-w[i])-14*abs(C[i]))) for i in range(3)]
        need(min(terms) > 0, ('open raw carrier roof', w, C))
        physical += min(terms)
        totals = [totals[i]+terms[i] for i in range(3)]
    return tuple(Q(t, denominator) for t in totals), Q(physical, denominator)


def main():
    output = Path(__file__).parents[1]/'05-knowledge/results/lrc14_empty_core_certificate_sep06.out'
    text = output.read_text()
    records = {}
    for line in text.splitlines():
        found = re.match(r'sector (\(.*?\)) rho \S+ defects (\(.*?\)) slope (\S+) intercept (\S+) cutoff (\S+) head (\d+) ', line)
        if found:
            p, ds, s, intercept, cutoff, n = found.groups()
            records[ast.literal_eval(p)] = (ast.literal_eval(ds), Q(s), Q(intercept), Q(cutoff), int(n))
    patterns = {tuple(sorted(p)) for p in product(range(1, 19), repeat=3)
                if 4 <= sum(p) <= 20 and sum(p) % 2 == 0
                and gcd(gcd(p[0], p[1]), p[2]) == 1
                and sum(x % 3 == 0 for x in p) <= 1}
    need(len(patterns) == 73 and set(records) == patterns-{(1, 1, 2)},
         ('independent complete pattern universe', len(patterns), len(records)))
    for pattern in sorted(patterns):
        actual = compile_independently(pattern)
        if pattern == (1, 1, 2):
            need(actual[-1] is None, ('norm-four count exception', actual))
        else:
            need(actual == records[pattern][:4], ('polygon-clipping compiler equality', pattern, actual))
    residue_audit()
    max_height = max((r[3].numerator-1)//r[3].denominator for r in records.values())
    speed_bank = [v for v in range(1, max_height+1, 2) if v % 3]
    signed = {p: sectors(p) for p in records}
    members = Counter()
    proof_rows = []
    complete_universe = 0
    for w in combinations(speed_bank, 3):
        if gcd(gcd(w[0], w[1]), w[2]) != 1:
            continue
        complete_universe += 1
        found = []
        for p, record in records.items():
            if w[2] < record[3] and any(sum(v[i]*w[i] for i in range(3)) == 0 for v in signed[p]):
                members[p] += 1
                found.append(p)
        if found:
            proof_rows.append(w)
    need(all(members[p] == record[4] for p, record in records.items()),
         ('independent complete per-pattern head counts', dict(members)))
    digest = sha256()
    equalities = []
    for w in sorted(proof_rows, key=lambda w: (w[2], w[0], w[1])):
        live = full_box(w)
        values, physical = projections(w, live)
        need(min(values) <= Q(6, 77), ('direct exact head consumer', w, values))
        if min(values) == Q(6, 77):
            equalities.append(w)
        digest.update(repr((w, tuple(sorted(live)), values, physical)).encode())
    claimed_digest = re.search(r'complete_head_digest (\w+)', text).group(1)
    need(digest.hexdigest() == claimed_digest, ('complete head semantic equality', digest.hexdigest()))
    need(equalities == [(1, 5, 11)], ('equality universe', equalities))
    print('PASS independent short-relation slice certificate audit')
    print('POLYTOPE COMPILER: all 73 patterns, every signed permutation, exact 2D clipping and third cross coordinate')
    print('RESIDUES: complete F3 kernel and affine defect fibers, both owner classes')
    print('DIRECT SPEED UNIVERSE', complete_universe, 'maximum height', max_height)
    print('HEAD MEMBERSHIPS', sum(members.values()), 'UNIQUE HEAD ROWS', len(proof_rows))
    print('EXACT CONSUMER: full integer boxes; three projections; physical mass; unique equality', equalities)
    print('COMPLETE HEAD SEMANTIC SHA256', digest.hexdigest())
    print('CHECKS', CHECKS)


if __name__ == '__main__':
    main()
