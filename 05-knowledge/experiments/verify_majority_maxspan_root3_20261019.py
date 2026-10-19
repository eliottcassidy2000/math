#!/usr/bin/env python3
"""Independent standard-library audit of the max-span/root-three reduction."""
from itertools import permutations, combinations
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from random import Random
from pathlib import Path
import json
SHAPES = [(0, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 0), (1, 1, 0), (1, 0, 1), (1, 1, 1)]
VARS = ('C', 'Z', 'U', 'W', 'D', 'F', 'V', 'T')
LOW = {(): None, ('b',): 0, ('z',): 1, ('b', 'z'): 2, ('z', 'b'): 3}
UP = {(): None, ('a',): 4, ('z',): 5, ('z', 'a'): 6, ('a', 'z'): 7}

def relations(shape):
    ba, bz, za = shape
    return [('x', 'a'), ('b', 'y')] + ([('b', 'a')] if ba else []) + ([('b', 'z')] if bz else []) + ([('z', 'a')] if za else [])

def core_words(shape):
    for word in permutations('bxzya'):
        idx = {a: word.index(a) for a in word}
        if all((idx[a] < idx[b] for a, b in relations(shape))):
            i, j = sorted((idx['x'], idx['y']))
            yield (word, idx, word[:i], word[j + 1:])

def mon(poly, coeff, *letters):
    powers = [0] * 8
    for letter in letters:
        powers[VARS.index(letter)] += 1
    poly[tuple(powers)] += coeff

def expected(shape):
    p = defaultdict(int)
    for a, b in [('C', 'D'), ('C', 'V'), ('D', 'U')]:
        mon(p, 2, a, b)
    for a in ('C', 'D'):
        mon(p, 2, a)
    for a in ('U', 'V'):
        mon(p, 3, a)

    def X(r):
        mon(p, 1, 'C', 'F')
        mon(p, 2, 'C', 'T')
        mon(p, 1, 'C')
        mon(p, r, 'F')
        mon(p, 3, 'T')

    def Y(r):
        mon(p, 1, 'D', 'Z')
        mon(p, 2, 'D', 'W')
        mon(p, 1, 'D')
        mon(p, r, 'Z')
        mon(p, 3, 'W')
    if shape == (0, 0, 0):
        X(3)
        Y(3)
    elif shape == (0, 1, 0):
        X(3)
        mon(p, 2)
    elif shape == (0, 0, 1):
        Y(3)
        mon(p, 2)
    elif shape == (1, 0, 0):
        X(2)
        Y(2)
        mon(p, 4)
    elif shape == (1, 1, 0):
        X(2)
        mon(p, 3)
    elif shape == (1, 0, 1):
        Y(2)
        mon(p, 3)
    elif shape == (1, 1, 1):
        mon(p, 2)
    return {k: v for k, v in p.items() if v}

def symbolic_audit():
    result = []
    for shape in SHAPES:
        p = defaultdict(int)
        count = 0
        for word, idx, pre, suf in core_words(shape):
            count += 1
            coef = 4 - 3 * ((idx['x'] < idx['z']) + (idx['z'] < idx['y']))
            features = [q for q in (LOW[pre], UP[suf]) if q is not None]
            for mask in range(1 << len(features)):
                power = [0] * 8
                for j, q in enumerate(features):
                    if mask >> j & 1:
                        power[q] += 1
                p[tuple(power)] += coef
        p = {k: v for k, v in p.items() if v}
        if not p == expected(shape):
            raise AssertionError((shape, p, expected(shape)))
        if not all((v > 0 for v in p.values())):
            raise AssertionError('exact certificate check failed')
        result.append({'shape': list(shape), 'core_word_count': count, 'positive_monomials': len(p), 'minimum_coefficient': min(p.values())})
    return result

def close(pred):
    p = list(pred)
    for k in range(len(p)):
        for j in range(len(p)):
            if p[j] >> k & 1:
                p[j] |= p[k]
    if not all((not p[i] >> i & 1 for i in range(len(p)))):
        raise AssertionError('exact certificate check failed')
    return p

def count(pred):
    pred = tuple(pred)
    n = len(pred)
    full = (1 << n) - 1

    @lru_cache(None)
    def f(mask):
        if mask == full:
            return 1
        ans = 0
        for i in range(n):
            if not mask >> i & 1 and pred[i] & ~mask == 0:
                ans += f(mask | 1 << i)
        return ans
    return f(0)

def induced_count(pred, vertices, ordered=()):
    vs = list(vertices)
    at = {v: i for i, v in enumerate(vs)}
    p = [sum((1 << at[u] for u in vs if pred[v] >> u & 1)) for v in vs]
    for a, b in zip(ordered, ordered[1:]):
        p[at[b]] |= 1 << at[a]
    return count(p)

def compare_count(pred, a, b):
    p = pred[:]
    p[b] |= 1 << a
    return count(p)

def random_poset(shape, rng):
    na = rng.randrange(5)
    nb = rng.randrange(5)
    n = na + 5 + nb
    A = list(range(na))
    core = {a: na + i for i, a in enumerate('bxzya')}
    B = list(range(na + 5, n))
    p = [0] * n

    def edge(a, b):
        p[b] |= 1 << a
    for a, b in relations(shape):
        edge(core[a], core[b])
    for side in (A, B):
        for a, b in combinations(side, 2):
            if rng.random() < 0.35:
                edge(a, b)
    for a in A:
        for letter in ('x', 'y', 'a'):
            edge(a, core[letter])
        for letter in ('b', 'z'):
            if rng.random() < 0.5:
                edge(a, core[letter])
        for b in B:
            edge(a, b)
    for b in B:
        for letter in ('b', 'x', 'y'):
            edge(core[letter], b)
        for letter in ('z', 'a'):
            if rng.random() < 0.5:
                edge(core[letter], b)
    return (close(p), A, core, B)

def exterior_audit(trials=60):
    rng = Random(202610190625)
    checked = 0
    maxn = 0
    smallest = None
    for shape in SHAPES:
        for _ in range(trials):
            p, A, c, B = random_poset(shape, rng)
            n = len(p)
            maxn = max(maxn, n)
            x, z, y = (c[a] for a in 'xzy')
            incomps = lambda a: {b for b in range(n) if a != b and (not p[a] >> b & 1) and (not p[b] >> a & 1)}
            if not incomps(x) == {c[a] for a in 'byz'}:
                raise AssertionError('exact certificate check failed')
            if not incomps(y) == {c[a] for a in 'xza'}:
                raise AssertionError('exact certificate check failed')
            e0 = count(p)
            nxz = compare_count(p, x, z)
            nzy = compare_count(p, z, y)
            factored = [0, 0, 0]
            cachelo = {}
            cacheup = {}
            for word, idx, pre, suf in core_words(shape):
                if pre not in cachelo:
                    ordered = tuple((c[a] for a in pre))
                    cachelo[pre] = induced_count(p, A + list(ordered), ordered)
                if suf not in cacheup:
                    ordered = tuple((c[a] for a in suf))
                    cacheup[suf] = induced_count(p, list(ordered) + B, ordered)
                wt = cachelo[pre] * cacheup[suf]
                factored[0] += wt
                factored[1] += wt * (idx['x'] < idx['z'])
                factored[2] += wt * (idx['z'] < idx['y'])
            if not factored == [e0, nxz, nzy]:
                raise AssertionError((shape, p, factored, [e0, nxz, nzy]))
            delta = 4 * e0 - 3 * (nxz + nzy)
            if not delta >= 0:
                raise AssertionError('exact certificate check failed')
            gap = Fraction(delta, e0)
            smallest = gap if smallest is None else min(smallest, gap)
            checked += 1
    return {'posets_checked': checked, 'maximum_order': maxn, 'all_exact_factorizations_match': True, 'minimum_normalized_deficit': str(smallest)}

def selection_audit():
    counts = {}
    for n in range(2, 6):
        pairs = list(combinations(range(n), 2))
        seen = set()
        checked = 0
        for mask in range(1 << len(pairs)):
            p = [0] * n
            for j, (a, b) in enumerate(pairs):
                if mask >> j & 1:
                    p[b] |= 1 << a
            p = tuple(close(p))
            if p in seen:
                continue
            seen.add(p)
            inc = [(a, b) for a, b in pairs if not p[b] >> a & 1]
            if not inc:
                continue
            span = max((b - a for a, b in inc))
            for x, y in inc:
                if y - x != span:
                    continue
                if not p[x] & ~p[y] == 0:
                    raise AssertionError('exact certificate check failed')
                ux = {j for j in range(n) if p[j] >> x & 1}
                uy = {j for j in range(n) if p[j] >> y & 1}
                if not uy <= ux:
                    raise AssertionError('exact certificate check failed')
                common = [z for z in range(n) if z not in (x, y) and (not (p[x] >> z & 1 or p[z] >> x & 1 or p[y] >> z & 1 or p[z] >> y & 1))]
                if not all((x < z < y for z in common)):
                    raise AssertionError('exact certificate check failed')
                checked += 1
        counts[str(n)] = {'naturally_labelled_posets': len(seen), 'maxspan_pairs_checked': checked}
    return counts

def adjacency_obstruction_audit():
    p = [0, 1, 3, 0, 8]
    order = [0, 3, 1, 4, 2]
    N = count(p)
    if not N == 10:
        raise AssertionError('exact certificate check failed')
    data = []
    for i, a in enumerate(order):
        for b in order[i + 1:]:
            q = compare_count(p, a, b)
            if not 2 * q > N:
                raise AssertionError('exact certificate check failed')
            data.append([a, b, q])
    for a, b in zip(order, order[1:]):
        ua = {j for j in range(5) if p[j] >> a & 1}
        ub = {j for j in range(5) if p[j] >> b & 1}
        if not (p[a] & ~p[b] or not ub <= ua):
            raise AssertionError('exact certificate check failed')
    return {'extension_count': N, 'majority_order': order, 'forward_counts': data, 'no_adjacent_forward_critical_pair': True}
if __name__ == '__main__':
    result = {'status': 'PASS', 'seven_exact_certificates': symbolic_audit(), 'actual_exterior_factorization': exterior_audit(), 'maximum_span_selection': selection_audit(), 'adjacency_obstruction': adjacency_obstruction_audit()}
    target = Path(__file__).with_name('majority_maxspan_root3_audit_20261019.json')
    print(json.dumps(result, indent=2))
