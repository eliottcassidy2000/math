#!/usr/bin/env python3
"""Complete hereditary gcd relaxation with the clock32/96 shadow exclusions.

The analytic shadow proofs are in the companion note. This independent
source reimplements the inherited CRT/Hunter arithmetic and retains all
labelled-signature outputs. It is not an actual-phase cover enumeration.
"""
from argparse import ArgumentParser
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, combinations_with_replacement
from math import gcd, lcm
from pathlib import Path
import json

GATES = 0


def check(ok, label):
    global GATES
    if not ok:
        raise RuntimeError(label)
    GATES += 1


def k(q):
    return (q+6)//7


def meet(c, q, r):
    b = gcd(q, r)
    A, u = divmod(k(q), b)
    B, v = divmod(k(r), b)
    return c//lcm(q, r)*(b*A*B+A*v+B*u+max(0, u+v-b))


@lru_cache(None)
def cap(c, gs):
    qs = [c//g for g in gs]
    # Prim implementation, unlike the original producer's Kruskal pass.
    used = {0}
    credit = 0
    while len(used) < len(qs):
        value, j = max((meet(c, qs[i], qs[j]), j)
                       for i in used for j in range(len(qs)) if j not in used)
        used.add(j)
        credit += value
    return sum(g*k(q) for g, q in zip(gs, qs))-credit


def admissible(c, gs):
    if gcd(*gs) != 1:
        return False
    for a in range(2, c+1):
        if c % a == 0:
            residual = tuple(sorted(gcd(a, g) for g in gs if g % a))
            if cap(a, residual) < a:
                return False
    return True


def compile_levels():
    previous = (1,)
    levels = {}
    excluded = {}
    for d in range(1, 7):
        orders = [q for q in range(1, 6*d//(7-d)+1) if d*k(q) >= q]
        clocks = sorted({g*q for g in previous for q in orders})
        rows, removed = [], []
        proposals = 0
        for c in clocks:
            options = [g for g in previous if c % g == 0]
            for gs in combinations_with_replacement(options, d):
                proposals += 1
                if not admissible(c, gs):
                    continue
                if (d == 5 and c == 32) or (d == 6 and c == 96):
                    removed.append([c, list(gs)])
                else:
                    rows.append([c, list(gs)])
        previous = tuple(sorted({c for c, gs in rows}))
        levels[str(d)] = dict(body_size=13-d, gcds=list(previous), profiles=rows,
                             candidate_clocks=len(clocks), proposals=proposals)
        excluded[str(d)] = removed
    check(excluded['5'] == [[32, [1,1,2,4,4]]], 'complete clock32 candidates')
    check(excluded['6'] == [[96, [1,1,4,4,6,12]], [96, [1,3,4,4,6,12]]],
          'complete clock96 candidates after inheritance')
    check([max(x['gcds']) for x in levels.values()] == [1,2,4,9,30,90],
          'refined gcd ceilings')
    check([len(x['profiles']) for x in levels.values()] == [1,2,5,19,109,1213],
          'refined full profile counts')
    return levels, excluded


def block(c, q, start, step=1):
    shadow = {(start+step*j) % q for j in range(k(q))}
    return {j for j in range(c) if j % q in shadow}


def positive_controls():
    definitions = {
        30: [(15,2), (15,7), (15,12), (10,0), (10,5)],
        90: [(45,8), (45,23), (45,38), (30,3), (30,18), (15,0)],
    }
    for c, words in definitions.items():
        sets = [block(c, q, a) for q, a in words]
        check(set.union(*sets) == set(range(c)), 'positive full padded cover')
        check(sum(map(len, sets)) == c, 'positive disjoint padded partition')
        check(gcd(*(c//q for q, a in words)) == 1, 'primitive signature')
        # Realize each signature with an honestly safe thirteen-speed row.
        gs = tuple(c//q for q, a in words)
        body = [c*i for i in range(1, 14-len(gs))]
        tails = [g*(1+c*(1+14*(i+1))) for i,g in enumerate(gs)]
        row = body+tails
        check(len(row) == len(set(row)) == 13 and gcd(*row) == 1, 'row type')
        check(all(min(v % (14*c), (-v) % (14*c)) >= c for v in row),
              'safe x=1/(14c), independent of padded-cover realization')
    actual = {
        30: (8, 112, (25082,24992,85712,11073,123)),
        90: (7, 126, (542,55082,25292,211773,30513,51126)),
    }
    den = 1009
    expected_maxima = {30: (1,2,3,3,30,30), 90: (1,2,3,3,6,90)}
    for c, (m, p, tails) in actual.items():
        body = [c*i for i in range(1, m+1)]
        check(all(14*min(i*p % den, (-i*p) % den) > den for i in range(1,m+1)),
              'strictly body-safe rational y')
        masks = []
        for w, (q, start) in zip(tails, definitions[c]):
            literal = {j for j in range(c)
                       if 14*min(w*(p+j*den) % (c*den),
                                 (-w*(p+j*den)) % (c*den)) < c*den}
            check(literal == block(c, q, start), 'literal actual mask equals target')
            check(gcd(c, w) == c//q, 'actual exact effective order')
            masks.append(literal)
        check(len(set.union(*masks)) == sum(map(len,masks)) == c,
              'actual strict-danger partition of all lifts')
        check(len(set(body+list(tails))) == 13 and gcd(*(body+list(tails))) == 1,
              'actual full row distinct and primitive')
        check(all(min(v % (14*c), (-v) % (14*c)) >= c for v in body+list(tails)),
              'same actual row is safe at x=1/(14c)')
        maxima = tuple(max(gcd(*s) for s in combinations(body+list(tails),n))
                       for n in range(12,6,-1))
        check(maxima == expected_maxima[c], 'all global subset gcd maxima')
    return definitions, actual, expected_maxima


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--write-profiles', type=Path)
    args = parser.parse_args()
    levels, excluded = compile_levels()
    covers, actual, subset_maxima = positive_controls()
    payload = dict(threshold='1/14', predicate='necessary relaxed signatures only',
                   levels=levels, shadow_exclusions=excluded)
    raw = (json.dumps(payload, sort_keys=True, separators=(',', ':'))+'\n').encode()
    if args.write_profiles:
        args.write_profiles.write_bytes(raw)
    print('SCOPE: proved shadow exclusions plus finite-exact hereditary relaxation; LRC14 OPEN')
    for d, row in levels.items():
        print('tails', d, 'body', row['body_size'], 'gcds', row['gcds'],
              'profiles', len(row['profiles']), 'candidate_clocks', row['candidate_clocks'],
              'proposals', row['proposals'])
    print('shadow_exclusions', json.dumps(excluded, sort_keys=True))
    print('positive_padded_partitions', json.dumps(covers, sort_keys=True))
    print('actual_body_phase_partitions', json.dumps(actual, sort_keys=True), 'denominator',1009)
    print('actual_row_subset_gcd_maxima_sizes_12_to_7',json.dumps(subset_maxima,sort_keys=True))
    print('explicit_gates', GATES)
    print('profile_json_sha256', sha256(raw).hexdigest())


if __name__ == '__main__':
    main()
