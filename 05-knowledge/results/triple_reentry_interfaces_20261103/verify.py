#!/usr/bin/env python3
"""Exact fixed-source unit-interface certificates. Standard library only.

Counts use the OLD-permutation boundary identity inherited from
mrb_square_probe_spectrum_20261103/verify.py. independent_audit.py instead
counts NEW permutations. Local assertions use explicit exceptions.
"""
import argparse
from functools import lru_cache
import hashlib
import itertools
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def require(ok, why):
    if not ok:
        raise RuntimeError(why)


def rows_from_code(n, code):
    rows = [0] * n
    for bit, (i, j) in enumerate(itertools.combinations(range(n), 2)):
        a, b = (i, j) if code >> bit & 1 else (j, i)
        rows[a] |= 1 << b
    return rows


def strong(rows):
    reach = [r | 1 << i for i, r in enumerate(rows)]
    for k in range(len(rows)):
        for i in range(len(rows)):
            if reach[i] >> k & 1:
                reach[i] |= reach[k]
    return all(r == (1 << len(rows)) - 1 for r in reach)


def response(rows):
    n = len(rows)
    start, end = [0] * n, [0] * n
    q = [[0] * n for _ in rows]
    h = 0
    for word in itertools.permutations(range(n)):
        bad = [i for i in range(n - 1)
               if not rows[word[i]] >> word[i + 1] & 1]
        if not bad:
            h += 1
            start[word[0]] += 1
            end[word[-1]] += 1
            for a, b in zip(word, word[1:]):
                q[a][b] += 1
        elif len(bad) == 1:
            i = bad[0]
            q[word[i]][word[i + 1]] += 1
    linear = [end[i] - start[i] + sum(q[i]) for i in range(n)]
    pair = {(i, j): q[i][j] + q[j][i]
            for i, j in itertools.combinations(range(n), 2)}
    values = [h + sum(linear[i] for i in range(n) if s >> i & 1)
              - sum(v for (i, j), v in pair.items() if s >> i & 1 and s >> j & 1)
              for s in range(1 << n)]
    require(values[0] == values[-1] == h, 'constant cuts')
    require(all(v > 0 and v % 2 for v in values), 'positive odd raw counts')
    return values


@lru_cache(None)
def support(x):
    out = set()
    p = 2
    while p * p <= x:
        if x % p == 0:
            out.add(p)
            while x % p == 0:
                x //= p
        p += 1
    if x > 1:
        out.add(x)
    return frozenset(out)


def reduced_edges(edges):
    edges = sorted(set(edges), key=lambda e: (len(e), tuple(sorted(e))))
    return [e for i, e in enumerate(edges) if not any(d <= e for d in edges[:i])]


def min_product(edges, bound=None):
    """Complete finite branch-and-bound, returns bound if no smaller solution."""
    edges = reduced_edges(edges)
    if any(not e for e in edges):
        return bound
    if bound is None:
        bound = math.prod(set().union(*edges)) + 1
    def search(es, prod):
        nonlocal bound
        if prod >= bound:
            return
        if not es:
            bound = prod
            return
        for p in sorted(es[0]):
            search([e for e in es if p not in e], prod * p)
    search(edges, 1)
    return bound


def minimal_transversals(edges):
    """Independent incremental algorithm, keeping the inclusion-minimal sets."""
    candidates = {frozenset()}
    for edge in reduced_edges(edges):
        candidates = set(reduced_edges(s | {p} for s in candidates for p in edge))
    return sorted((sorted(s) for s in candidates), key=lambda s: (math.prod(s), s))


def cost(mask, n, two_sided):
    w = mask.bit_count()
    return min(w, n - w) if two_sided else w


def min_obstruction(values, width, two_sided, bound=None):
    n = (len(values) - 1).bit_length()
    low = {v for s, v in enumerate(values) if cost(s, n, two_sided) <= width}
    result = None
    for target in sorted(set(values) - low):
        edges = [support(v) - support(target) for v in low]
        if any(not e for e in edges):
            continue
        m = min_product(edges, bound)
        if bound is None or m < bound:
            bound = m
            result = {'modulus': m, 'target': target,
                      'mask': values.index(target)}
    return result


def describe(rows, width, two_sided):
    values = response(rows)
    witness = min_obstruction(values, width, two_sided)
    require(witness is not None, 'expected obstruction')
    m = witness['modulus']
    n = len(rows)
    units = [s for s, v in enumerate(values) if math.gcd(v, m) == 1]
    require(min(cost(s, n, two_sided) for s in units) == width + 1,
            'exact minimum interface cost')
    low = {v for s, v in enumerate(values) if cost(s, n, two_sided) <= width}
    full_edges = [support(v) for v in set(values)]
    low_edges = [support(v) for v in low]
    witness.update(rows=rows, H=values[0], strong=strong(rows), values=values,
                   width=width, two_sided=two_sided, unit_masks=units,
                   low_minimal_blockers=minimal_transversals(low_edges),
                   full_minimal_blockers=minimal_transversals(full_edges))
    restricted = [support(v) - support(witness['target']) for v in low]
    require(min(map(math.prod, minimal_transversals(restricted))) == m,
            'two optimization algorithms disagree')
    return witness


def small_controls():
    path_values = []
    previous, current = 1, 3
    for m in range(13):
        brute = sum(2 ** s.bit_count() for s in range(1 << m) if not s & (s << 1))
        binomial = sum(math.comb(m + 1 - k, k) * 2 ** k
                       for k in range((m + 1) // 2 + 1))
        closed = (2 ** (m + 2) - (-1) ** (m + 2)) // 3
        require(brute == binomial == closed, 'path index identity')
        if m == 0:
            require(brute == previous, 'path empty base')
        elif m == 1:
            require(brute == current, 'path singleton base')
        else:
            previous, current = current, current + 2 * previous
            require(brute == current, 'path recurrence')
        path_values.append(brute)
    divisor_examples = []
    for x in (1335, 4005, 6675, 118815):
        divisors = [d for d in range(2, x) if x % d == 0]
        f = len(divisors)
        s = sum(all(d % (p * p) for p in support(d)) for d in divisors)
        u = sum(len(support(d)) == 1 and next(iter(support(d))) == d for d in divisors)
        divisor_examples.append({'N': x, 'F': f, 'S': s, 'U': u})
        require(math.prod(support(x)) == 1335, 'radical promotion')
        require((f == s + u) == (x != 1335), 'balance witness')
    # Same two-source triangle transfer sees a two-path-cover sidecar.
    return {'path_values': path_values, 'divisor_examples': divisor_examples}


def census():
    out = []
    for n in range(1, 7):
        total = 0
        strong_total = 0
        # At order six only strong parents are in the declared universe.
        minima = {'fixed': None, 'either': None, 'strong_fixed': None}
        digest = hashlib.sha256()
        for code in range(1 << (n * (n - 1) // 2)):
            rows = rows_from_code(n, code)
            is_strong = strong(rows)
            if n == 6 and not is_strong:
                continue
            total += 1
            strong_total += is_strong
            values = response(rows)
            digest.update(json.dumps([code, values], separators=(',', ':')).encode() + b'\n')
            for key, two_sided in (('fixed', False), ('either', True), ('strong_fixed', False)):
                if key == 'strong_fixed' and not is_strong:
                    continue
                old = minima[key]
                found = min_obstruction(values, 2, two_sided,
                                        old['modulus'] if old else None)
                if found:
                    minima[key] = dict(found, code=code, rows=rows, H=values[0])
        out.append({'n': n, 'universe': 'strong labelled' if n == 6 else 'all labelled',
                    'total': total, 'strong_total': strong_total, 'minima': minima,
                    'response_sha256': digest.hexdigest()})
    require(out[4]['minima']['fixed']['modulus'] == 255, 'order5 minimum')
    require(out[4]['minima']['strong_fixed'] is None, 'order5 strong boundary')
    require(out[5]['minima']['fixed']['modulus'] == 1335, 'order6 fixed minimum')
    require(out[5]['minima']['either']['modulus'] == 72320955, 'order6 either minimum')
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--census', action='store_true')
    ap.add_argument('--write', action='store_true')
    args = ap.parse_args()
    if args.census:
        payload, filename = census(), 'census.json'
    else:
        payload = {'controls': small_controls(), 'examples': [
            describe([4, 1, 2, 7, 15], 2, False),
            describe([32, 21, 1, 7, 13, 30], 2, False),
            describe([32, 9, 3, 5, 15, 30], 2, True),
            describe([244, 61, 232, 1, 236, 8, 170, 42], 3, True)]}
        require(set(payload['examples'][1]['values']) == set(payload['examples'][2]['values']),
                'same full image, different cost profiles')
        filename = 'receipt.json'
    serialized = json.dumps(payload, indent=2, sort_keys=True) + '\n'
    if args.write:
        (HERE / filename).write_text(serialized)
    else:
        require((HERE / filename).read_text() == serialized, 'receipt changed: ' + filename)
    print(json.dumps({'status': 'PASS', 'file': filename,
                      'sha256': hashlib.sha256(serialized.encode()).hexdigest()}))


if __name__ == '__main__':
    main()
