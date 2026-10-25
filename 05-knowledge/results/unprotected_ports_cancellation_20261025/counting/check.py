"""Exact finite controls for the cycle-count moment transfer."""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations
from math import comb, factorial, gcd, prod
import json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def two(e):
    return Fraction(1 << e) if e >= 0 else Fraction(1, 1 << -e)


def beta(k, e):
    return (-1) ** e * comb(k + e - 1, e) if k else int(e == 0)


def moments_target(n):
    if type(n) is not int or n < 0:
        raise ValueError('target must be a nonnegative integer')
    if n == 0:
        return Fraction(1), {'keys': 1, 'scalar_cells': 1, 'rows': []}
    top = n if n % 2 else n - 1
    states = {(n, (0,) * ((top + 1) // 2)): [Fraction(0)] + [Fraction(1)] * n}
    rows = []
    for ell in range(top, 0, -2):
        nxt = {}
        for (remaining, defect), vector in states.items():
            selected = defect[(ell - 1) // 2] if ell <= remaining else 0
            for k in range(remaining // ell + 1):
                for e in (range((remaining - ell * k) // ell + 1) if k else (0,)):
                    rest = remaining - ell * (k + e)
                    if ell == 1 and rest:
                        continue
                    future = range(1, min(ell - 2, rest) + 1, 2)
                    dnext = tuple(defect[(d - 1) // 2] + k * (d - gcd(ell, d)) for d in future)
                    key = (rest, dnext)
                    out = nxt.setdefault(key, [Fraction(0)] * (rest + 1))
                    common = beta(k, e) / Fraction(ell ** k * factorial(k))
                    exponent = k * (ell - 1) // 2 + ell * k * (k - 1) // 2 - k * selected
                    for r in (range(1, rest + 1) if rest else (0,)):
                        source_index = ell * k + r
                        out[r] += common * two(exponent + r * k) * vector[source_index]
        nxt = {key: vec for key, vec in nxt.items() if any(vec)}
        rows.append({'length': ell, 'keys': len(states),
                     'scalar_cells': sum(max(1, key[0]) for key in states)})
        states = nxt
    need(set(states) <= {(0, ())}, ('terminal moment key', n))
    value = states.get((0, ()), [Fraction(0)])[0]
    return value, {'keys': sum(row['keys'] for row in rows),
                   'scalar_cells': sum(row['scalar_cells'] for row in rows), 'rows': rows}


def old_target(n):
    if n == 0:
        return Fraction(1), {'keys': 1}
    top = n if n % 2 else n - 1
    states = {(n, (0,) * ((top + 1) // 2)): Fraction(1)}
    keys = 0
    rows = []
    for ell in range(top, 0, -2):
        nxt = defaultdict(Fraction)
        keys += len(states)
        rows.append({'length': ell, 'keys': len(states)})
        for (remaining, response), w in states.items():
            selected = response[(ell - 1) // 2] if ell <= remaining else 0
            for k in range(remaining // ell + 1):
                exponent = k * (ell - 1) // 2 + ell * k * (k - 1) // 2 + k * selected
                common = w * two(exponent) / (ell ** k * factorial(k))
                for e in (range((remaining - ell * k) // ell + 1) if k else (0,)):
                    rest = remaining - ell * (k + e)
                    if ell == 1 and rest:
                        continue
                    response_next = tuple(response[(d - 1) // 2] + k * gcd(ell, d)
                                          for d in range(1, min(ell - 2, rest) + 1, 2))
                    nxt[(rest, response_next)] += common * beta(k, e)
        states = {key: w for key, w in nxt.items() if w}
    return states.get((0, ()), Fraction(0)), {'keys': keys, 'rows': rows}


def adaptive_target(n):
    """Sparse m weights until their support exceeds the remaining horizon."""
    if n == 0:
        return Fraction(1), {'keys': 1, 'scalar_cells': 1, 'rows': []}
    top = n if n % 2 else n - 1
    states = {(n, (0,) * ((top + 1) // 2)): ('sparse', {0: Fraction(1)})}
    rows = []
    for ell in range(top, 0, -2):
        nxt = {}

        def dense_from_sparse(remaining, support):
            if remaining == 0:
                return [sum(support.values(), Fraction(0))]
            return [Fraction(0)] + [sum((w * two(r*m) for m, w in support.items()), Fraction(0))
                                    for r in range(1, remaining + 1)]

        def add_sparse(key, m, w):
            remaining = key[0]
            if remaining:
                need((w > 0) == ((n-remaining-m) % 2 == 0), 'prefix sign coherence')
            if key not in nxt:
                nxt[key] = ('sparse', {}) if remaining else ('dense', [Fraction(0)])
            kind, data = nxt[key]
            if kind == 'sparse':
                if m in data:
                    need((data[m] > 0) == (w > 0), 'same-m terms never cancel')
                data[m] = data.get(m, Fraction(0)) + w
                if len(data) > remaining:
                    nxt[key] = ('dense', dense_from_sparse(remaining, data))
            elif remaining:
                for r in range(1, remaining + 1):
                    data[r] += w * two(r*m)
            else:
                data[0] += w

        def add_dense(key, vector):
            remaining = key[0]
            if key not in nxt:
                nxt[key] = ('dense', [Fraction(0)] * (remaining + 1))
            kind, data = nxt[key]
            if kind == 'sparse':
                data = dense_from_sparse(remaining, data)
                nxt[key] = ('dense', data)
            for r in (range(1, remaining+1) if remaining else (0,)):
                data[r] += vector[r]

        rows.append({'length': ell, 'keys': len(states),
                     'scalar_cells': sum(len(data) if kind == 'sparse' else max(1,key[0])
                                         for key,(kind,data) in states.items())})
        for (remaining, defect), (kind, data) in states.items():
            selected = defect[(ell-1)//2] if ell <= remaining else 0
            for k in range(remaining//ell+1):
                for e in (range((remaining-ell*k)//ell+1) if k else (0,)):
                    rest = remaining-ell*(k+e)
                    if ell == 1 and rest:
                        continue
                    dnext = tuple(defect[(d-1)//2]+k*(d-gcd(ell,d))
                                  for d in range(1,min(ell-2,rest)+1,2))
                    key = (rest,dnext)
                    factor = beta(k,e)/Fraction(ell**k*factorial(k))
                    exponent = k*(ell-1)//2+ell*k*(k-1)//2-k*selected
                    if kind == 'sparse':
                        for m,w in data.items():
                            add_sparse(key,m+k,w*factor*two(exponent+ell*k*m))
                    else:
                        vector = [Fraction(0)]*(rest+1)
                        for r in (range(1,rest+1) if rest else (0,)):
                            vector[r] = factor*two(exponent+r*k)*data[ell*k+r]
                        # Retain zero nonterminal moment vectors: their semantic
                        # m support still exists and is used in the storage proof.
                        add_dense(key,vector)
        states = {key: entry for key,entry in nxt.items()
                  if key[0] or entry[1][0]}
    need(set(states) <= {(0,())}, 'adaptive terminal key')
    value = states.get((0,()),('dense',[Fraction(0)]))[1][0]
    return value, {'keys':sum(row['keys'] for row in rows),
                   'scalar_cells':sum(row['scalar_cells'] for row in rows),'rows':rows}


def partitions(n, low=1):
    if n == 0:
        yield ()
    for first in range(low, n + 1, 2):
        for tail in partitions(n - first, first):
            yield (first,) + tail


def direct(n):
    total = Fraction(0)
    for actual in range(n + 1):
        excess = n - actual
        for p in partitions(actual):
            mult = Counter(p)
            coeff = [1] + [0] * excess
            for d, k in mult.items():
                out = [0] * (excess + 1)
                for i, value in enumerate(coeff):
                    for e in range((excess - i) // d + 1):
                        out[i + d * e] += value * (-1) ** e * comb(k + e - 1, e)
                coeff = out
            exponent = sum((d - 1) // 2 for d in p) + sum(gcd(a, b) for a, b in combinations(p, 2))
            denominator = prod(d ** k * factorial(k) for d, k in mult.items())
            total += Fraction((1 << exponent) * coeff[excess], denominator)
    return total


def rank(matrix):
    rows = [[Fraction(v) for v in row] for row in matrix]
    pivots = 0
    for col in range(len(rows[0]) if rows else 0):
        pivot = next((i for i in range(pivots, len(rows)) if rows[i][col]), None)
        if pivot is None:
            continue
        rows[pivot], rows[pivots] = rows[pivots], rows[pivot]
        value = rows[pivots][col]
        rows[pivots] = [v / value for v in rows[pivots]]
        for i in range(pivots + 1, len(rows)):
            value = rows[i][col]
            rows[i] = [a - value * b for a, b in zip(rows[i], rows[pivots])]
        pivots += 1
    return pivots


def main():
    samples = []
    known = [1, 1, 0, 1, 2, 5, 26, 255, 4602, 147095, 8228360]
    for n in range(19):
        value, new = moments_target(n)
        old, before = old_target(n)
        adaptive, compressed = adaptive_target(n)
        need(value == old == adaptive == direct(n), ('four coefficient paths', n))
        need(value.denominator == 1, ('integer core count', n))
        need(new['keys'] <= before['keys'], ('coarser key graph', n))
        if n:
            need(all(a['length'] == b['length'] and a['scalar_cells'] <= b['keys']
                     for a,b in zip(compressed['rows'],before['rows'])), ('stagewise storage domination',n))
        if n < len(known):
            need(value == known[n], ('inherited core prefix', n))
        if n in (6, 10, 14, 18):
            samples.append({'n': n, 'old_scalar_keys': before['keys'],
                            'new_vector_keys': new['keys'], 'new_scalar_cells': new['scalar_cells'],
                            'adaptive_scalar_cells':compressed['scalar_cells'],
                            'adaptive_keys':compressed['keys']})
    comparisons = 0
    for s in range(1, 8):
        for h in range(1, 10):
            matrix = [[2 ** (r * m) for m in range(1, h + 1)] for r in range(1, s + 1)]
            need(rank(matrix) == min(s, h), ('exact moment rank', s, h))
            comparisons += 1
    alpha = (8, -6, 1)
    moment = lambda r: sum(a * 2 ** (r * m) for m, a in enumerate(alpha, 1))
    need(moment(1) == moment(2) == 0 and moment(3) == 192, 'horizon hostile')
    for ell, length in ((3, 9), (5, 15), (7, 105)):
        for m in range(1, 10):
            need(all(d * m - m * gcd(d, length) == 0 for d in range(1, ell + 1, 2)),
                 'realizable same-defect multiplicities')
    for bad in (-1, True, 1.5, '3'):
        try:
            moments_target(bad)
        except ValueError:
            pass
        else:
            raise RuntimeError('bad target accepted')
    print(json.dumps({'status': 'PASS', 'coefficient_range': [0, 18], 'vandermonde_checks': comparisons,
                      'strict_linear_kernel': {'cycle_counts': [1, 2, 3], 'coefficients': list(alpha),
                                               'moments_1_2_3': [moment(1), moment(2), moment(3)]},
                      'work_counts_not_timings': samples,
                      'scope': 'adaptive scalar-cell domination; no bit-cost or runtime improvement claim'}, indent=2))


if __name__ == '__main__':
    main()
