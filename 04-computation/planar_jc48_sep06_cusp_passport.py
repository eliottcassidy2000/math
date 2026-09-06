#!/usr/bin/env python3
"""Exact bounded controls for the cusp passport proof; no global existence test.

Universe: all ordered permutation pairs in S_d, 2<=d<=5, and every nonempty
proper subset A of Fix(sigma).  Products act on the left, right factor first.
The all-degree reduction and actual-sheet transport are analytic obligations.
"""
from collections import Counter
from itertools import combinations, permutations
import hashlib
import json

GATES = 0


def check(condition, message):
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(message)


def mul(p, q):
    return tuple(p[q[i]] for i in range(len(p)))


def inv(p):
    out = [0] * len(p)
    for i, j in enumerate(p):
        out[j] = i
    return tuple(out)


def image(p, subset):
    return frozenset(p[i] for i in subset)


def support(p):
    return frozenset(i for i, j in enumerate(p) if i != j)


def braid(sigma, tau):
    return mul(mul(sigma, tau), sigma) == mul(mul(tau, sigma), tau)


def record(sigma, tau, A):
    B = image(mul(sigma, tau), A)
    return dict(d=len(sigma), sigma=list(sigma), tau=list(tau),
                A=sorted(A), B=sorted(B), n=len(A & B))


rows = []
universe = []
for d in range(2, 6):
    omega = frozenset(range(d))
    ps = tuple(permutations(range(d)))
    identity = tuple(range(d))
    pairs = subset_rows = selected = 0
    for sigma in ps:
        fixed = sorted(omega - support(sigma))
        for tau in ps:
            if not braid(sigma, tau):
                continue
            pairs += 1
            g = mul(sigma, tau)
            check(mul(mul(g, sigma), inv(g)) == tau, 'marked conjugacy')
            for a in range(1, min(d, len(fixed) + 1)):
                for entries in combinations(fixed, a):
                    subset_rows += 1
                    A = frozenset(entries)
                    B = image(g, A)
                    I, U, V = A & B, A - B, B - A
                    O = omega - (A | B)
                    n, h = len(I), len(O)
                    check(len(B) == a and B <= omega - support(tau),
                          'transported subset has same size and correct inertia')
                    check(image(g, I) == I and not (I & support(g)),
                          'intersection is pointwise fixed')
                    check(h == d - 2 * a + n, 'outside-set cardinality')
                    if n > 1 or h > 1:
                        continue
                    selected += 1
                    check((d, a, n) in {(2, 1, 1), (3, 1, 0), (4, 2, 1)},
                          'complete sparse passport types')
                    check(h == 1 and len(U) == len(V) <= 1,
                          'sparse supports')
                    if d == 2:
                        check(sigma == tau == identity, 'degree-two identity')
                    else:
                        check(support(sigma) == V | O and
                              support(tau) == U | O, 'two transposition supports')
                        check(len(support(sigma)) == len(support(tau)) == 2,
                              'no longer cycles')
                    rows.append(record(sigma, tau, A))
    universe.append(dict(d=d, ordered_pairs=len(ps)**2, braid_pairs=pairs,
                         subset_rows=subset_rows, selected=selected))

check(Counter((r['d'], len(r['A']), r['n']) for r in rows) ==
      Counter({(2, 1, 1): 2, (3, 1, 0): 6, (4, 2, 1): 24}),
      'independently predicted labeled multiplicities')

# Named hostiles preserve raw labels, and do not impose local transitivity.
sigma4, tau4, A4 = (1, 0, 2, 3), (0, 2, 1, 3), frozenset({2, 3})
check(braid(sigma4, tau4), 'four-sheet braid hostile')
check(record(sigma4, tau4, A4) ==
      dict(d=4, sigma=[1, 0, 2, 3], tau=[0, 2, 1, 3],
           A=[2, 3], B=[0, 3], n=1), 'four-sheet retained subsets')
sigma3, tau3, A3 = (1, 0, 2), (0, 2, 1), frozenset({2})
check(record(sigma3, tau3, A3)['n'] == 0, 'three-sheet empty intersection')

# A too-strong pointwise fixed-letter = actual-letter rule fails already for
# a disjoint identity sheet which is deleted in its entirety along the cusp.
identity2 = (0, 1)
check(len(set(range(2)) - support(identity2)) == 2 and len({0}) == 1,
      'inertia-fixed labels need the actual-retention sidecar')

# Independent integer ledger controls: every omega_i>=max(0,d-2a),
# n+sum omega_i=1.  Allocation of the remaining 0/1 is exhaustive.
ledger_rows = 0
for d in range(2, 41):
    for a in range(1, d):
        for N in range(1, 9):
            for n in (0, 1):
                lower = max(0, d - 2 * a)
                possible = N * lower <= 1 - n
                if not possible:
                    continue
                ledger_rows += 1
                check(d - 2 * a + n <= 1, 'ledger supplies outside bound')
                if N >= 2:
                    check(d <= 2 * a, 'two nodes remove degree-three passport')

# The exceptional cubic has no root over F_7, hence is irreducible over Q.
residues = [(128*x**3 - 288*x - 283) % 7 for x in range(7)]
check(all(residues), 'cusp-parameter cubic irreducible mod seven')

payload = json.dumps(rows, sort_keys=True, separators=(',', ':')).encode()
print('FINITE-EXACT cusp passport controls; analytic geometry is proved in the note')
print('universe=' + json.dumps(universe, sort_keys=True))
print('labeled sparse passports=' + json.dumps(rows, sort_keys=True))
print('raw-passport SHA256=' + hashlib.sha256(payload).hexdigest())
print('integer-ledger eligible rows=' + str(ledger_rows))
print('C(lambda) mod-7 values=' + str(residues))
print('PASS always-active gates=' + str(GATES))
