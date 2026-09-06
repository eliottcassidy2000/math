#!/usr/bin/env python3
"""Exact hereditary gcd sieve with padded CRT overlaps at threshold 1/14.

Universe: all abstract gcd signatures with d=1,...,6 tails, after the
proved pivot reduction; these are necessary arithmetic filters, not
counterexamples. No repository mathematical producer is imported.
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


def teeth(q):
    return (q + 6) // 7


@lru_cache(None)
def divisors(c):
    return tuple(a for a in range(2, c + 1) if c % a == 0)


def pair_lower(c, q, r):
    """Minimum CRT meet of two padded consecutive-unit orbit blocks."""
    a = gcd(q, r)
    A, u = divmod(teeth(q), a)
    B, v = divmod(teeth(r), a)
    return (c // lcm(q, r)) * (a*A*B + A*v + B*u + max(0, u+v-a))


@lru_cache(None)
def union_upper(c, gs):
    """Hunter upper bound; gs are exact gcds with this chosen clock."""
    qs = tuple(c // g for g in gs)
    total = sum(g * teeth(q) for g, q in zip(gs, qs))
    edges = sorted(((pair_lower(c, qs[i], qs[j]), i, j)
                    for i, j in combinations(range(len(gs)), 2)), reverse=True)
    parent = list(range(len(gs)))
    def root(i):
        while parent[i] != i:
            i = parent[i]
        return i
    for weight, i, j in edges:
        x, y = root(i), root(j)
        if x != y:
            parent[x] = y
            total -= weight
    return total


def admissible(c, gs):
    if gcd(*gs) != 1:
        return False
    for a in divisors(c):
        residual = tuple(sorted(gcd(a, g) for g in gs if g % a != 0))
        if union_upper(a, residual) < a:
            return False
    return True


EXPECTED_GCDS = (
    (1,),
    (1, 2),
    (1, 2, 3, 4),
    (1, 2, 3, 4, 6, 8, 9),
    (1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 15, 16, 18, 24, 30, 32),
    (1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 15, 16, 17, 18, 20,
     22, 23, 24, 25, 27, 29, 30, 32, 33, 34, 36, 40, 44, 45,
     46, 48, 50, 51, 54, 58, 60, 64, 66, 72, 88, 90, 96),
)
EXPECTED_COUNTS = (1, 2, 5, 19, 110, 1217)


def compile_profiles():
    previous = (1,)
    levels = {}
    for d in range(1, 7):
        # beta(q)>=1/d implies q <= floor(6d/(7-d)); no empirical cutoff.
        orders = tuple(q for q in range(1, 6*d//(7-d) + 1)
                       if d*teeth(q) >= q)
        candidates = sorted({g*q for g in previous for q in orders})
        profiles = []
        proposed = 0
        for c in candidates:
            options = tuple(g for g in previous if c % g == 0)
            for gs in combinations_with_replacement(options, d):
                proposed += 1
                if admissible(c, gs):
                    profiles.append([c, list(gs)])
        current = tuple(sorted({c for c, _ in profiles}))
        check(current == EXPECTED_GCDS[d-1], f'complete gcd set d={d}')
        check(len(profiles) == EXPECTED_COUNTS[d-1], f'profile count d={d}')
        levels[str(d)] = dict(body_size=13-d, pivot_orders=list(orders),
                             candidate_clocks=len(candidates), proposals=proposed,
                             gcds=list(current), profiles=profiles)
        previous = current
    return levels


def literal_pair_controls():
    # Every order pair1..28, both orientation gauges, every relative shift
    # modulo their common order. Literal masks on the common labelled grid.
    cases = 0
    for q in range(1, 29):
        for r in range(1, 29):
            c, a = lcm(q, r), gcd(q, r)
            lower = pair_lower(c, q, r)
            attained = c + 1
            for sign in (-1, 1):
                first = {j for j in range(c) if j % q < teeth(q)}
                for shift in range(a):
                    second = {j for j in range(c)
                              if (sign*j-shift) % r < teeth(r)}
                    meet = len(first & second)
                    check(meet >= lower, 'literal CRT lower bound')
                    attained = min(attained, meet)
                    cases += 1
            check(attained == lower, 'pair lower bound is attained in relaxation')
    return cases


def safe_at_fourteen_clock(c, gs, body_size):
    # Realizes the gcd signature, with a rational witness x=1/(14c).
    body = tuple(c*i for i in range(1, body_size+1))
    tails = tuple(g*(1+c*(1+14*(i+1))) for i, g in enumerate(gs))
    speeds = body + tails
    check(len(set(speeds)) == 13 and gcd(*speeds) == 1, 'safe realization type')
    check(tuple(gcd(c, w) for w in tails) == gs, 'realization signature')
    check(all(14*min(v % (14*c), (-v) % (14*c)) >= 14*c for v in speeds),
          'explicit weak-safe witness')
    return list(tails)


def controls(levels):
    check(pair_lower(288, 16, 9) == 12, 'coprime-grid forced overlap')
    check(sum(g*teeth(36//g) for g in (1,4,4,6,9)) == 37,
          'scalar hostile passes before overlap')
    check(union_upper(36, (1,4,4,6,9)) < 36, 'overlap excludes old c36')
    check(union_upper(288, (1,1,4,4,18,32)) < 288, 'overlap excludes old c288')
    check(not admissible(4, (1,1,1,4)), 'absorption closes inherited-only q1')
    check(not admissible(8, (1,1,1,4)), 'proper-divisor absorption matters')
    # The maximum surviving clocks are achieved by honest safe rows.
    named = {}
    for d in (4, 5, 6):
        rows = levels[str(d)]['profiles']
        c = max(row[0] for row in rows)
        for clock, gs in rows:
            if clock == c:
                named[str((c, tuple(gs)))] = safe_at_fourteen_clock(c, tuple(gs), 13-d)
    # d=7: primes retain arbitrary clock size in this arithmetic relaxation.
    # The proof for all primes is in the note; this finite list is a control.
    for p in (2,3,5,7,11,13,17,29,31,101,1009):
        check(admissible(p, (1,)*7), 'seven-tail prime boundary')
        safe_at_fourteen_clock(p, (1,)*7, 6)
    return named


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--write-profiles', type=Path)
    args = parser.parse_args()
    cases = literal_pair_controls()
    levels = compile_profiles()
    named = controls(levels)
    payload = dict(threshold='1/14', predicate='necessary arithmetic filters only', levels=levels)
    raw = (json.dumps(payload, sort_keys=True, separators=(',', ':'))+'\n').encode()
    if args.write_profiles:
        args.write_profiles.write_bytes(raw)
    print('SCOPE: analytic pivot + finite-exact complete arithmetic compiler; LRC14 OPEN')
    for d, row in levels.items():
        print('tails', d, 'body', row['body_size'], 'gcds', row['gcds'],
              'profiles', len(row['profiles']), 'candidate_clocks', row['candidate_clocks'],
              'proposals', row['proposals'])
    print('literal_pair_cases', cases, 'explicit_gates', GATES)
    print('named_safe_maximum_clock_realizations', json.dumps(named, sort_keys=True))
    print('profile_json_sha256', sha256(raw).hexdigest())


if __name__ == '__main__':
    main()
