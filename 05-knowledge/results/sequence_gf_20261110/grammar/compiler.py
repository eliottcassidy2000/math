#!/usr/bin/env python3
"""Exact structural examples, not sequence fitting. Standard library only."""
from collections import Counter, deque
from itertools import permutations, product
from math import comb, factorial
import json


def tournament(n, backwards=()):
    a = [[int(i < j) for j in range(n)] for i in range(n)]
    for i, j in backwards:
        a[i][j], a[j][i] = 1, 0
    return a


def profile(a):
    """Ordered covers via literal permutations and legal cuts."""
    n = len(a)
    f = [0] * (n + 1)
    for p in permutations(range(n)):
        b = sum(not a[p[i]][p[i + 1]] for i in range(n - 1))
        for c in range(b + 1, n + 1):
            f[c] += comb(n - 1 - b, c - 1 - b)
    return f


def qvalue(f, j):
    return sum(f[c] * comb(j, c) for c in range(1, min(j, len(f) - 1) + 1))


def kernel(q, g):
    """First g vertices are growing modules, remaining b are actual fixed vertices.

    This finite compiler requires q[:g] transitive. Enumerate its quotient
    words, using fixed vertices once and growing labels at least once.
    """
    n, b = len(q), len(q) - g
    if any(q[i][j] and q[j][k] and q[k][i]
           for i, j, k in permutations(range(g), 3)):
        raise ValueError('Growing quotient is cyclic: no bounded-run compiler.')
    current = {(0, (0,) * g, -1): 1}
    out = Counter()
    # At most b+1 increasing growing runs and b fixed vertices.
    for _ in range(g * (b + 1) + b):
        following = Counter()
        for (mask, counts, last), ways in current.items():
            for nxt in range(n):
                if last >= 0 and not q[last][nxt]:
                    continue
                if nxt < g:
                    if counts[nxt] >= b + 1:
                        continue
                    c = list(counts)
                    c[nxt] += 1
                    state = (mask, tuple(c), nxt)
                else:
                    bit = 1 << (nxt - g)
                    if mask & bit:
                        continue
                    state = (mask | bit, counts, nxt)
                following[state] += ways
                if state[0] == (1 << b) - 1 and min(state[1]) > 0:
                    out[state[1]] += ways
        current = following
    return out


def compile_bases(q, seeds):
    weights = kernel(q, len(seeds))
    profiles = [profile(a) for a in seeds]
    bases = Counter()
    for c, w in weights.items():
        for js in product(*(range(1, x + 1) for x in c)):
            coefficient, base = w, 1
            for i, (ci, j) in enumerate(zip(c, js)):
                coefficient *= (-1) ** (ci - j) * comb(ci, j)
                base *= qvalue(profiles[i], j)
            bases[base] += coefficient
    return weights, {k: v for k, v in sorted(bases.items()) if v}


def join_copies(a, r):
    n = len(a)
    return [[a[i % n][j % n] if i // n == j // n else int(i // n < j // n)
             for j in range(n * r)] for i in range(n * r)]


def substitute(q, blocks):
    labels = [(i, v) for i, a in enumerate(blocks) for v in range(len(a))]
    return [[blocks[i][u][v] if i == j else q[i][j]
             for j, v in labels] for i, u in labels]


def hamilton(a):
    n = len(a)
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not dp[mask][v]:
                continue
            for w in range(n):
                if not mask >> w & 1 and a[v][w]:
                    dp[mask | (1 << w)][w] += dp[mask][v]
    return sum(dp[-1])


def polynomial(roots):
    """Descending coefficients of product (X-root)."""
    p = [1]
    for root in roots:
        q = [0] * (len(p) + 1)
        for i, c in enumerate(p):
            q[i] += c
            q[i + 1] -= root * c
        p = q
    return p


def one_apex(word):
    h, s, p = 1, 0, 0
    for bit in word:
        if bit == 0:
            h += s + p
        s, p = 2 * s + bit, bit
    return h


def one_apex_tournament(word):
    n = len(word)
    a = tournament(n + 1)
    for i, bit in enumerate(word):
        a[n][i], a[i][n] = bit, 1 - bit
    return a


def capped_fiber(h):
    """Reachable exact finite recognizer for the fiber H=h."""
    cap = h + 1
    states = [(1, 0, 0)]
    indices = {states[0]: 0}
    transitions = []
    for v in states:
        row = []
        for bit in (0, 1):
            a, s, p = v
            a = min(cap, a + (s + p if not bit else 0))
            nxt = (a, min(cap, 2 * s + bit), bit)
            if nxt not in indices:
                indices[nxt] = len(states)
                states.append(nxt)
            row.append(indices[nxt])
        transitions.append(row)
    counts, values = [1] + [0] * (len(states) - 1), []
    for _ in range(18):
        values.append(sum(c for c, v in zip(counts, states) if v[0] == h))
        nxt = [0] * len(states)
        for c, row in zip(counts, transitions):
            for j in row:
                nxt[j] += c
        counts = nxt
    return len(states), values


def main():
    k1 = tournament(1)
    c3 = tournament(3, [(2, 0)])
    cases = [('one_growing_C3_seed', c3, [c3]),
             ('two_growing_singletons', c3, [k1, k1]),
             ('two_growing_mixed', c3, [c3, k1])]
    report = {'convention': 'r>=1; counted objects are literal directed paths', 'contexts': []}
    for name, q, seeds in cases:
        weights, bases = compile_bases(q, seeds)
        terms = [sum(c * lam ** r for lam, c in bases.items()) for r in range(1, 9)]
        controls = []
        for r in range(1, 4):
            blocks = [join_copies(a, r) for a in seeds] + [k1] * (len(q) - len(seeds))
            literal = hamilton(substitute(q, blocks))
            if literal != terms[r - 1]:
                raise ArithmeticError((name, r, literal, terms[r - 1]))
            controls.append({'r': r, 'vertices': sum(map(len, blocks)), 'H': literal})
        recurrence = polynomial(bases)
        for i in range(len(recurrence) - 1, len(terms)):
            if sum(c * terms[i - j] for j, c in enumerate(recurrence)):
                raise ArithmeticError('Recurrence failure')
        report['contexts'].append({'name': name,
            'kernel': {','.join(map(str, k)): v for k, v in sorted(weights.items())},
            'exponential_coefficients': bases, 'terms_r_1_to_8': terms,
            'minimal_characteristic_polynomial_descending': recurrence,
            'vertex_DP_controls': controls})
    cycle_controls = []
    for r in range(1, 5):
        value = hamilton(substitute(c3, [join_copies(k1, r)] * 3))
        floor = 3 * factorial(r) ** 3
        if value < floor:
            raise ArithmeticError('Cyclic lower bound failure')
        cycle_controls.append({'r': r, 'H': value, 'factorial_floor': floor})
    report['cyclic_hostile'] = cycle_controls
    for n in range(7):
        for word in product((0, 1), repeat=n):
            if one_apex(word) != hamilton(one_apex_tournament(word)):
                raise ArithmeticError(('One-apex control', word))
    ceiling = 511
    cores = {h: Counter() for h in range(3, ceiling + 1, 2)}
    for n in range(2, (ceiling - 1).bit_length() + 2):
        for inner in product((0, 1), repeat=n - 2):
            word = (1,) + inner + (0,)
            h = one_apex(word)
            if h <= ceiling:
                cores[h][n] += 1
    fibers = []
    for h in [1, 3, 5, 9, 11, 13, 15, 21, 25]:
        size, values = capped_fiber(h)
        numerator = {0: 1} if h == 1 else dict(sorted(cores[h].items()))
        expected = [sum(c * (n - k + 1) for k, c in numerator.items() if k <= n)
                    for n in range(len(values))]
        if values != expected:
            raise ArithmeticError(('Fiber mismatch', h, values, expected))
        fibers.append({'H': h, 'capped_states': size, 'fiber_numerator': numerator,
                       'fiber_denominator': '(1-z)^2', 'counts_n_0_to_10': values[:11]})
    periodic = [one_apex((1, 0) * r) for r in range(10)]
    if periodic != [(4 ** (r + 1) + 6 * r + 5) // 9 for r in range(10)]:
        raise ArithmeticError('Periodic formula failure')
    report['one_apex'] = {'literal_controls': 127, 'fibers': fibers,
        'periodic_10_terms_r_0_to_9': periodic,
        'periodic_10_GF': '(1-3z)/((1-4z)(1-z)^2)',
        'direct_support_through_511_size': 1 + sum(bool(v) for v in cores.values()),
        'direct_support_through_511_first_misses': [h for h, v in cores.items() if not v][:20]}
    report['status'] = 'All new-formula controls passed; no all-odd support claim.'
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
