#!/usr/bin/env python3
"""Standalone exact controls for odd-braid support and H4 single cycles.
No mathematical implementation is imported. All gates remain active with -O.
"""
from itertools import combinations, permutations
from math import gcd
from hashlib import sha256
import json

GATES = 0

def need(value, message):
    global GATES
    GATES += 1
    if not value:
        raise RuntimeError(message)

def mul(a, b):
    return tuple(a[j] for j in b)

def alt(a, b, n):
    out = tuple(range(len(a)))
    for j in range(n):
        out = mul(out, a if j % 2 == 0 else b)
    return out

def braid(a, b, n):
    return alt(a, b, n) == alt(b, a, n)

def supp(a):
    return frozenset(i for i, j in enumerate(a) if i != j)

def cycle(d, values):
    out = list(range(d))
    values = tuple(values)
    for x, y in zip(values, values[1:] + values[:1]):
        out[x] = y
    return tuple(out)

def cycles(d, m):
    for values in combinations(range(d), m):
        for tail in permutations(values[1:]):
            yield cycle(d, (values[0],) + tail)

def orbit(gens):
    seen = {0}
    todo = [0]
    while todo:
        x = todo.pop()
        for g in gens:
            y = g[x]
            if y not in seen:
                seen.add(y)
                todo.append(y)
    return seen

report = {}
# Unfiltered ordered pairs in every S_d, 2 <= d <= 5, with odd labels3,5,7.
# Every pair satisfying the literal relation is retained, including identity
# and multi-cycle permutations. No equal-cycle-type filter is imposed.
bank = []
for d in range(2, 6):
    perms = list(permutations(range(d)))
    row = {'degree': d, 'ordered_pairs': len(perms) ** 2, 'solutions': {}}
    for n in (3, 5, 7):
        r = (n - 1) // 2
        accepted = 0
        for a in perms:
            A = supp(a)
            for b in perms:
                if not braid(a, b, n):
                    continue
                accepted += 1
                B = supp(b)
                need((r + 1) * len(A & B) >= len(B), 'odd support bound')
                need(len(A) == len(B), 'odd braid conjugacy')
                for x in B:
                    y = x
                    met = False
                    for _ in range(r + 1):
                        met |= y in A
                        y = b[y]
                    need(met, 'forbidden consecutive outside run')
        row['solutions'][n] = accepted
    bank.append(row)
report['unfiltered_pair_bank'] = bank

# All relative positions of a single m-cycle against a fixed standard one:
# j shared labels chosen arbitrarily; all m-j fresh labels are present;
# all (m-1)! cycles on each support are retained, for2 <= m <=7.
# Simultaneous relabeling fixes the first standard cycle and can place the
# fresh letters at the chosen labels, so every union degree2m-j is represented.
pair_rows = []
for m in range(2, 8):
    row = {'length': m, 'examined': 0, 'braid3': {}, 'braid5': {}}
    for j in range(1, m + 1):
        d = 2 * m - j
        a = cycle(d, range(m))
        A = supp(a)
        counts = {3: 0, 5: 0}
        for shared in combinations(range(m), j):
            values = shared + tuple(range(m, d))
            for tail in permutations(values[1:]):
                b = cycle(d, (values[0],) + tail)
                row['examined'] += 1
                for n in (3, 5):
                    if braid(a, b, n):
                        counts[n] += 1
                        need(((n + 1) // 2) * j >= m, 'single-cycle support bound')
                        if n == 3 and 2 * j == m:
                            B = supp(b)
                            need(all((x in B) != (a[x] in B) for x in A), 'half-overlap alternation')
        for n in (3, 5):
            if counts[n]:
                row['braid' + str(n)][j] = counts[n]
    pair_rows.append(row)
report['single_cycle_pair_bank'] = pair_rows

# Independent full tuple census. The first m-cycle is fixed by simultaneous
# conjugacy. Enumerate all other m-cycles on each declared label universe;
# the filters are exactly the six H4 relations. In particular transitivity
# is checked after, not assumed from a union statistic or a quotient.
quad_rows = []
for m in range(2, 6):
    for d in range(m + 1, min(3 * m, 10) + 1):
        a = cycle(d, range(m))
        A = supp(a)
        pool = list(cycles(d, m))
        near = [b for b in pool if braid(a, b, 3)]
        comm = [c for c in pool if mul(a, c) == mul(c, a)]
        accepted = transitive = 0
        for b in near:
            for c in comm:
                if not braid(b, c, 3):
                    continue
                for e in comm:
                    if mul(b, e) != mul(e, b) or not braid(c, e, 5):
                        continue
                    accepted += 1
                    need(A == supp(b) == supp(c) == supp(e), 'full H4 tuple support')
                    transitive += len(orbit((a, b, c, e))) == d
        need(transitive == 0, 'proper single-cycle transitive H4 action')
        quad_rows.append([m, d, len(pool), len(near), len(comm), accepted, transitive])
report['full_H4_tuple_bank'] = quad_rows

# Sharpness and missing-relation controls.
a = cycle(5, (0, 1, 2)); b = cycle(5, (2, 3, 4))
need(braid(a, b, 5), 'five-letter positive braid5 hostile')
need(len(supp(a) & supp(b)) == 1, 'sharp third-support equality')
need(2 * len(supp(a) & supp(b)) < len(supp(a)), 'ordinary half-bound must fail')
need(not braid(a, b, 3), 'five-cusp cannot inherit ordinary relation')
chain = [cycle(5, (i, i + 1)) for i in range(4)]
a, b, c, e = chain
need(braid(a, b, 3) and braid(b, c, 3), 'A4 chain ordinary edges')
need(all(mul(x, y) == mul(y, x) for x, y in ((a, c), (a, e), (b, e))), 'A4 chain commuting edges')
need(braid(c, e, 3) and not braid(c, e, 5), 'terminal five-relation is necessary')
need(len(orbit(chain)) == 5 and len(supp(a)) == 2, 'A4 chain transitive fixed-letter hostile')
for m in range(2, 10):
    a = cycle(m, range(m))
    need(braid(a, a, 3) and braid(a, a, 5), 'same-support positive relation')
    need(len(orbit((a, a, a, a))) == m and not (set(range(m)) - supp(a)), 'fixed-free positive action')

# Direct conflict bank at the nontrivial equality boundary: c alternates I
# and its complement. Test every single cycle d on I plus m/2 fresh labels.
# Every such candidate has all required support positions; none may braid5.
conflicts = []
for m in (2, 4, 6, 8):
    d = 3 * m // 2
    c = cycle(d, range(m))
    I = tuple(range(0, m, 2))
    values = I + tuple(range(m, d))
    tried = 0
    for tail in permutations(values[1:]):
        e = cycle(d, (values[0],) + tail)
        tried += 1
        need(not braid(c, e, 5), 'alternating half-intersection conflict')
    conflicts.append([m, tried])
report['alternating_conflict_bank'] = conflicts
semantic = json.dumps(report, sort_keys=True, separators=(',', ':')).encode()
print(json.dumps(report, sort_keys=True, indent=2))
print('semantic_sha256=' + sha256(semantic).hexdigest())
print('always_active_exact_gates=' + str(GATES))
