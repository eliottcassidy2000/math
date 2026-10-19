#!/usr/bin/env python3
"""Independent whole-poset extension counting for the all-c two-slave theorem.

No imports from the main verifier; standard library only. Finite checks are
implementation audits, not substitutes for the symbolic all-c proof.
"""
from functools import lru_cache
from random import Random
from fractions import Fraction
from pathlib import Path
import json

def close(p):
    p = list(p)
    for k in range(len(p)):
        for j in range(len(p)):
            if p[j] >> k & 1:
                p[j] |= p[k]
    if not all((not p[j] >> j & 1 for j in range(len(p)))):
        raise AssertionError('exact certificate check failed')
    return p

def exact_score(p, x, y, common):
    n = len(p)
    full = (1 << n) - 1
    cm = sum((1 << t for t in common))

    @lru_cache(None)
    def f(mask):
        if mask == full:
            return (1, 0)
        E = N = 0
        for v in range(n):
            if mask >> v & 1 or p[v] & ~mask:
                continue
            e, s = f(mask | 1 << v)
            added = int(bool(mask >> x & 1)) + int(not mask >> y & 1) if cm >> v & 1 else 0
            E += e
            N += s + added * e
        return (E, N)
    ans = f(0)
    return (ans, f.cache_info().currsize)

def make_actual(c, trial, rng):
    na = rng.randrange(4)
    nb = rng.randrange(4)
    common = list(range(4, 4 + c))
    A = list(range(4 + c, 4 + c + na))
    B = list(range(4 + c + na, 4 + c + na + nb))
    n = 4 + c + na + nb
    p = [0] * n

    def edge(a, b):
        p[b] |= 1 << a
    edge(0, 3)
    edge(2, 1)
    density = (0, 0.2, 0.5, 1)[trial % 4]
    for i, t in enumerate(common):
        for u in common[i + 1:]:
            if rng.random() < density:
                edge(t, u)
        if rng.random() < 0.4:
            edge(2, t)
        if rng.random() < 0.4:
            edge(t, 3)
    if rng.random() < 0.5:
        edge(2, 3)
    for side in (A, B):
        for i, a in enumerate(side):
            for b in side[i + 1:]:
                if rng.random() < 0.4:
                    edge(a, b)
    for a in A:
        for t in (0, 1, 3):
            edge(a, t)
        for t in [2] + common:
            if rng.random() < 0.5:
                edge(a, t)
        for b in B:
            edge(a, b)
    for b in B:
        for t in (0, 1, 2):
            edge(t, b)
        for t in common + [3]:
            if rng.random() < 0.5:
                edge(t, b)
    return (close(p), common, na, nb)

def main():
    rng = Random(202610190642)
    rows = []
    maxn = maxstates = 0
    minimum = None
    for c in range(4, 8):
        best = None
        smallest = None
        for trial in range(100):
            p, common, na, nb = make_actual(c, trial, rng)
            inc = lambda t: {u for u in range(len(p)) if u != t and (not (p[t] >> u & 1 or p[u] >> t & 1))}
            if not inc(0) == {1, 2, *common}:
                raise AssertionError('exact certificate check failed')
            if not inc(1) == {0, 3, *common}:
                raise AssertionError('exact certificate check failed')
            (E, N), states = exact_score(p, 0, 1, common)
            delta = 4 * c * E - 3 * N
            if not delta >= 0:
                raise AssertionError((c, p, E, N))
            ratio = Fraction(N, E)
            gap = Fraction(delta, E)
            best = ratio if best is None else max(best, ratio)
            smallest = gap if smallest is None else min(smallest, gap)
            maxn = max(maxn, len(p))
            maxstates = max(maxstates, states)
        rows.append({'common_count': c, 'actual_posets': 100, 'maximum_observed_sum': str(best), 'minimum_normalized_deficit': str(smallest)})
    missing = 0
    for c in range(9):
        for L in range(c + 1):
            p = [0, 0, 1 | sum((1 << 3 + j for j in range(L)))] + [sum((1 << 3 + i for i in range(j))) for j in range(c)]
            (E, N), _ = exact_score(p, 0, 1, range(3, c + 3))
            if not E == (c + 3) * (c - L + 1) * (c + L + 2) // 2:
                raise AssertionError('exact certificate check failed')
            if not 7 * c * E - 6 * N == c * L * (L + 1) * (c - L + 1):
                raise AssertionError('exact certificate check failed')
            missing += 1
    result = {'status': 'PASS', 'actual_exterior_results': rows, 'actual_posets_total': 400, 'maximum_order': maxn, 'largest_subset_dp_state_count': maxstates, 'missing_slave_formula_checks': missing}
    out = Path(__file__).with_name('two_slave_allc_exteriors_independent_audit_20261019.json')
    print(json.dumps(result, indent=2))
if __name__ == '__main__':
    main()
