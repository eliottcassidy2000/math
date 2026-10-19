#!/usr/bin/env python3
"""Direct exact-polynomial proof for zero lower and at most three upper slaves.

Standard library only. The polynomial is derived from shuffle counts, not fitted.
The independent recurrences are implementation audits, not the all-size proof.
"""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations_with_replacement, combinations, permutations
from random import Random
from pathlib import Path
import json

class Algebra:

    def __init__(self, n):
        self.n = n
        self.zero = (0,) * n

    def C(self, q):
        return {self.zero: F(q)} if q else {}

    def V(self, i):
        e = list(self.zero)
        e[i] = 1
        return {tuple(e): F(1)}

    def add(self, *ps):
        out = defaultdict(F)
        for p in ps:
            for e, q in p.items():
                out[e] += q
        return {e: q for e, q in out.items() if q}

    def scale(self, p, q):
        return {e: v * q for e, v in p.items() if v * q}

    def sub(self, p, q):
        return self.add(p, self.scale(q, -1))

    def mul(self, *ps):
        out = self.C(1)
        for p in ps:
            new = defaultdict(F)
            for e, q in out.items():
                for f, r in p.items():
                    new[tuple((a + b for a, b in zip(e, f)))] += q * r
            out = {e: q for e, q in new.items() if q}
        return out

    def power(self, p, k):
        return self.mul(*[p] * k)

    def choose(self, p, k):
        out = self.C(1)
        for j in range(k):
            out = self.scale(self.mul(out, self.sub(p, self.C(j))), F(1, j + 1))
        return out

    def evaluate(self, p, vs):
        return sum((q * prod((v ** a for v, a in zip(vs, e))) for e, q in p.items()))

def prod(xs):
    ans = 1
    for x in xs:
        ans *= x
    return ans

def compositions(n, k):
    if k == 1:
        yield (n,)
        return
    for j in range(n + 1):
        for tail in compositions(n - j, k - 1):
            yield ((j,) + tail)

def derive(m):
    a = Algebra(m + 1)
    vs = [a.V(i) for i in range(m + 1)]
    c = a.add(*vs)
    E0 = {}
    N = {}
    D = {}
    cuts = 0
    for hs in compositions(m + 1, m + 1):
        blocks = [i for i, h in enumerate(hs) for _ in range(h)]
        if any((blocks[j] < j for j in range(1, m + 1))):
            continue
        cuts += 1
        weight = a.C(1)
        ranks = []
        offset = {}
        for i, (h, v) in enumerate(zip(hs, vs)):
            width = v if i == m else a.sub(v, a.C(1))
            weight = a.mul(weight, a.choose(a.add(width, a.C(h)), h))
            ranks += [a.add(offset, a.scale(width, F(j, h + 1))) for j in range(1, h + 1)]
            offset = a.add(offset, v)
        E0 = a.add(E0, weight)
        nword = a.add(a.scale(a.mul(c, a.add(a.scale(c, 3), a.C(2 * m + 5))), F(1, 2)), a.scale(a.mul(a.add(c, a.C(m + 1)), ranks[0]), -1), *ranks[1:])
        dword = a.add(a.scale(a.mul(c, a.sub(a.C(2 * m + 1), c)), F(1, 2)), a.scale(a.mul(a.add(c, a.C(m + 1)), ranks[0]), 3), a.scale(a.add(*ranks[1:]), -3))
        N = a.add(N, a.mul(weight, nword))
        D = a.add(D, a.mul(weight, dword))
    E = a.mul(a.add(c, a.C(m + 2)), E0)
    if not D == a.sub(a.scale(a.mul(c, E), 4), a.scale(N, 3)):
        raise AssertionError('exact certificate check failed')
    return (a, vs, c, E, N, D, cuts)

def positivity(m, a, vs, c, E, D):
    if m == 0:
        if not D == a.mul(c, E):
            raise AssertionError('exact certificate check failed')
        return {'certificate': 'D=cE; root-exchange symmetry', 'positive_terms': len(D)}
    if m <= 2:
        if not all((q > 0 for q in D.values())):
            raise AssertionError('exact certificate check failed')
        return {'certificate': 'All monomial coefficients are positive', 'positive_terms': len(D)}
    if not m == 3:
        raise AssertionError('exact certificate check failed')
    q, t = (vs[1], vs[2])
    H = a.add(a.scale(a.power(q, 3), 6), a.scale(a.mul(a.power(q, 2), t), 6), a.scale(a.power(q, 2), 18), a.scale(a.mul(q, a.power(t, 2)), 4), a.scale(a.mul(q, t), 14), a.scale(q, 12), a.power(t, 3), a.scale(a.power(t, 2), 5), a.scale(t, 6))
    boundary = a.scale(a.mul(a.add(t, a.C(1)), a.sub(a.add(q, t), a.C(1)), H), F(1, 24))
    actual_boundary = {e: v for e, v in D.items() if e[0] == 0 and e[3] == 0}
    if not boundary == actual_boundary:
        raise AssertionError('exact certificate check failed')
    remainder = a.sub(D, boundary)
    if not all((v > 0 for v in remainder.values())):
        raise AssertionError('exact certificate check failed')
    if not all((e[0] > 0 or e[3] > 0 for e in remainder)):
        raise AssertionError('exact certificate check failed')
    if not (all((v > 0 for v in H.values())) and a.zero not in H):
        raise AssertionError('exact certificate check failed')
    return {'certificate': 'Positive p/s remainder plus integer-positive boundary factor', 'all_deficit_terms': len(D), 'positive_remainder_terms': len(remainder), 'boundary_terms': len(boundary), 'boundary_factor': '(t+1)(q+t-1)H(q,t)/24', 'H_terms': [{'powers': list(e), 'coefficient': str(v)} for e, v in sorted(H.items())]}

def threshold_count(c, L):
    m = len(L)

    @lru_cache(None)
    def f(i, j, mask):
        if i == c and j == m and (mask == 3):
            return (1, 0)
        E = N = 0

        def add(state, score=0):
            nonlocal E, N
            e, n = f(*state)
            E += e
            N += n + score * e
        if i < c:
            add((i + 1, j, mask), int(bool(mask & 1)) + int(not mask & 2))
        if j < m and mask & 1 and (i >= L[j]):
            add((i, j + 1, mask))
        if not mask & 1:
            add((i, j, mask | 1))
        if not mask & 2:
            add((i, j, mask | 2))
        return (E, N)
    return f(0, 0, 0)

def threshold_audit(derived, maxc=20):
    rows = []
    for m, item in enumerate(derived):
        a, vs, cvar, Epol, Npol, Dpol, cuts = item
        checked = 0
        for c in range(maxc + 1):
            for L in combinations_with_replacement(range(c + 1), m):
                if m:
                    gap = (L[0],) + tuple((L[i] - L[i - 1] for i in range(1, m))) + (c - L[-1],)
                else:
                    gap = (c,)
                E, N = threshold_count(c, L)
                if not a.evaluate(Epol, gap) == E:
                    raise AssertionError((m, c, L, 'E'))
                if not a.evaluate(Npol, gap) == N:
                    raise AssertionError((m, c, L, 'N'))
                if not a.evaluate(Dpol, gap) == 4 * c * E - 3 * N >= 0:
                    raise AssertionError((m, c, L, 'D'))
                checked += 1
        rows.append({'upper_slaves': m, 'threshold_cases': checked, 'largest_common_count': maxc})
    return rows

def close(p):
    p = list(p)
    for k in range(len(p)):
        for j in range(len(p)):
            if p[j] >> k & 1:
                p[j] |= p[k]
    if not all((not p[j] >> j & 1 for j in range(len(p)))):
        raise AssertionError('exact certificate check failed')
    return p

def actual_count(p, common):
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
            added = int(bool(mask & 1)) + int(not mask & 2) if cm >> v & 1 else 0
            E += e
            N += s + added * e
        return (E, N)
    return f(0)

def actual_exteriors(trials=100):
    rng = Random(202610190652)
    rows = []
    maxn = 0
    for m in range(4):
        bestgap = None
        for trial in range(trials):
            c = rng.randrange(9)
            na = rng.randrange(4)
            nb = rng.randrange(4)
            U = list(range(2, 2 + m))
            Z = list(range(2 + m, 2 + m + c))
            A = list(range(2 + m + c, 2 + m + c + na))
            B = list(range(2 + m + c + na, 2 + m + c + na + nb))
            n = 2 + m + c + na + nb
            p = [0] * n
            maxn = max(maxn, n)

            def edge(a, b):
                p[b] |= 1 << a
            for a in U:
                edge(0, a)
            for side in (U, Z, A, B):
                for i, a in enumerate(side):
                    for b in side[i + 1:]:
                        if rng.random() < 0.4:
                            edge(a, b)
            for z in Z:
                for a in U:
                    if rng.random() < 0.4:
                        edge(z, a)
            for a in A:
                for v in [0, 1] + U:
                    edge(a, v)
                for z in Z:
                    if rng.random() < 0.5:
                        edge(a, z)
                for b in B:
                    edge(a, b)
            for b in B:
                for v in (0, 1):
                    edge(v, b)
                for v in Z + U:
                    if rng.random() < 0.5:
                        edge(v, b)
            p = close(p)
            inc = lambda t: {u for u in range(n) if u != t and (not (p[t] >> u & 1 or p[u] >> t & 1))}
            if not inc(0) == {1, *Z}:
                raise AssertionError('exact certificate check failed')
            if not inc(1) == {0, *Z, *U}:
                raise AssertionError('exact certificate check failed')
            E, N = actual_count(p, Z)
            delta = 4 * c * E - 3 * N
            if not delta >= 0:
                raise AssertionError((m, c, p, E, N))
            gap = F(delta, E)
            bestgap = gap if bestgap is None else min(bestgap, gap)
        rows.append({'upper_slaves': m, 'actual_posets': trials, 'minimum_normalized_deficit': str(bestgap)})
    return {'total_actual_posets': 4 * trials, 'maximum_order': maxn, 'families': rows}

def zero_common_suffix_certificates():
    families = []
    for m in range(4):
        n = m + 2
        possible = list(combinations(range(2, n), 2))
        seen = set()
        cases = []
        for mask in range(1 << len(possible)):
            p = [0, 0] + [1] * m
            for j, (a, b) in enumerate(possible):
                if mask >> j & 1:
                    p[b] |= 1 << a
            p = tuple(close(p))
            if p in seen:
                continue
            seen.add(p)
            orders = []
            for o in permutations(range(n)):
                at = {v: i for i, v in enumerate(o)}
                if all((not p[v] >> u & 1 or at[u] < at[v] for u in range(n) for v in range(n))):
                    orders.append(o)
            inc = [(a, b) for a, b in combinations(range(n), 2) if not (p[a] >> b & 1 or p[b] >> a & 1)]
            suffixes = {()}
            for o in orders:
                last = max(o.index(0), o.index(1))
                suffixes.update((o[j:] for j in range(last + 1, n)))
            suffixes = sorted(suffixes, key=lambda v: (len(v), v))
            generators = []
            for v in suffixes:
                compatible = [o for o in orders if not v or o[-len(v):] == v]
                generators.append((v, compatible))
            for M in orders:
                at = {v: i for i, v in enumerate(M)}
                if at[0] > at[1]:
                    continue
                span = at[1] - at[0]
                if any((abs(at[a] - at[b]) > span for a, b in inc)):
                    continue
                certificate = None
                for a, b in inc:
                    if at[a] > at[b]:
                        a, b = (b, a)
                    columns = []
                    for v, compatible in generators:
                        E = len(compatible)
                        N = sum((o.index(a) < o.index(b) for o in compatible))
                        columns.append({'suffix': list(v), 'extensions': E, 'forward_comparison_count': N, 'margin_3N_minus_2E': 3 * N - 2 * E})
                    if all((col['margin_3N_minus_2E'] <= 0 for col in columns)):
                        certificate = {'oriented_pair': [a, b], 'generator_columns': columns}
                        break
                if not certificate is not None:
                    raise AssertionError((m, p, M))
                cases.append({'predecessor_masks': list(p), 'candidate_majority_order': list(M), 'root_span': span, **certificate})
        families.append({'upper_slaves': m, 'naturally_labelled_upper_posets': len(seen), 'admissible_maxspan_orders': len(cases), 'single_pair_certificates': cases})
    if not [f['naturally_labelled_upper_posets'] for f in families] == [1, 1, 2, 7]:
        raise AssertionError('exact certificate check failed')
    if not [f['admissible_maxspan_orders'] for f in families] == [1, 2, 6, 48]:
        raise AssertionError('exact certificate check failed')
    return families

def encode(p):
    return [{'powers': list(e), 'coefficient': str(v)} for e, v in sorted(p.items())]

def main():
    derived = [derive(m) for m in range(4)]
    certs = []
    for m, (a, vs, c, E, N, D, cuts) in enumerate(derived):
        cert = positivity(m, a, vs, c, E, D)
        certs.append({'upper_slaves': m, 'valid_shuffle_cuts': cuts, 'deficit_degree': max(map(sum, D), default=0), 'positivity': cert, 'extension_polynomial': encode(E), 'aggregate_numerator_polynomial': encode(N), 'deficit_polynomial': encode(D)})
    if not [d[-1] for d in derived] == [1, 2, 5, 14]:
        raise AssertionError('exact certificate check failed')
    thresholds = threshold_audit(derived)
    exterior = actual_exteriors()
    zero_common = zero_common_suffix_certificates()
    E, N = threshold_count(1, (0, 0, 0, 1))
    if not (E, N) == (35, 48):
        raise AssertionError('exact certificate check failed')
    result = {'status': 'PASS', 'claim': 'Universal aggregate 4c/3 for no lower slaves and at most three upper slaves, or its order dual.', 'direct_symbolic_certificates': certs, 'independent_threshold_dp': thresholds, 'actual_exterior_audit': exterior, 'zero_common_maxspan_suffix_certificates': zero_common, 'four_upper_slave_obstruction': {'common_count': 1, 'thresholds': [0, 0, 0, 1], 'extensions': E, 'aggregate_numerator': N, 'deficit': -4}}
    out = Path(__file__).with_name('one_sided_three_slave_aggregate_audit_20261019.json')
    print(json.dumps(result, indent=2))
if __name__ == '__main__':
    main()
