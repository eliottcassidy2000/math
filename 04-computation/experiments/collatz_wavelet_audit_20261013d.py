#!/usr/bin/env python3
"""Independent strict-orbit controls and source-preserving wavelet packets.

No numpy, parity-bit input, randomness model, or universal convergence premise.
The finite orbit controls terminate explicitly. All checks survive python -O.
"""
import argparse
from collections import defaultdict
from functools import lru_cache
import hashlib
import json
from pathlib import Path


def require(test, message):
    if not test:
        raise RuntimeError(message)


@lru_cache(None)
def parity(n):
    value, steps = n, 0
    while value != 1:
        value = (3 * value + 1) // 2 if value % 2 else value // 2
        steps += 1
        require(steps < 100000, 'strict orbit cap')
    return 1 if steps % 2 == 0 else -1


def addresses(q, r):
    m = 3 ** q
    return ((m + 1) // 2 * r) % m, ((3 * m + 1) // 2 * (3 * r + 1)) % (3 * m)


def psi(q, r, n):
    return 3 * (n % (3 ** q) == r) - (n % (3 ** (q - 1)) == r % (3 ** (q - 1)))


@lru_cache(None)
def sparse_h(q, r, n):
    p = 3 ** (q - 1)
    start = r % p or p
    return sum(psi(q, r, j) * parity(j) for j in range(start, n + 1, p))


def packet_step(row, reduce=True):
    nxt = defaultdict(int)
    for (k, r), c in row.items():
        e, o = addresses(k + 1, r)
        nxt[k, e] += c
        nxt[k + 1, o] += c
    if not reduce:
        return dict(nxt)
    groups = {}
    for (k, r), c in nxt.items():
        p = 3 ** k
        groups.setdefault((k, r % p), [0, 0, 0])[r // p] = c
    result = {}
    for (k, a), cs in groups.items():
        median = sorted(cs)[1]
        p = 3 ** k
        for j, c in enumerate(cs):
            if c != median:
                result[k, a + j * p] = c - median
    return result


def gram_audit(census):
    # Entirely separate path: direct strict orbits and literal prefix arrays.
    limit = 3 * (1 << 11)
    ts = [0] + [parity(n) for n in range(1, limit + 1)]
    hs = {}
    for q in range(1, 5):
        for r in range(3 ** q):
            values = [0]
            for n in range(1, limit + 1):
                values.append(values[-1] + psi(q, r, n) * ts[n])
            hs[q, r] = values
    checks = 0
    for q in range(1, 4):
        for r in range(3 ** q):
            e, o = addresses(q, r)
            for n in range(0, (1 << 12) + 1):
                require(hs[q, r][n] == -hs[q, e][n // 2] - hs[q + 1, o][3 * ((n + 1) // 2)], 'all-N identity')
                checks += 1
        p = 3 ** (q - 1)
        for a in range(p):
            for n in range(0, limit + 1):
                require(sum(hs[q, a + j * p][n] for j in range(3)) == 0, 'sibling zero')
    rows = 0
    for item in census['rows']:
        if item['k'] > 10:
            continue
        q, r, x = item['q'], item['r'], item['X']
        e, o = addresses(q, r)
        sums = dict.fromkeys(('xx', 'yy', 'xy', 'even', 'carry_cross', 'carry_square', 'child_energy'), 0)
        for m in range(x, 2 * x):
            a, b, even = hs[q, e][m], hs[q + 1, o][3 * m], hs[q, r][2 * m]
            odd = hs[q, r][2 * m + 1]
            z = odd - even
            values = (a*a, b*b, a*b, even*even, even*z, z*z, even*even + odd*odd)
            for key, value in zip(sums, values):
                sums[key] += value
        require(all(sums[key] == int(item[key]) for key in sums), 'independent Gram row')
        rows += 1
    # Same second moments do not specify their cross-term.
    require(1*1 + 1*1 + 2*(1*1) == 4, 'positive Gram hostile')
    require(1*1 + (-1)*(-1) + 2*(1*(-1)) == 0, 'negative Gram hostile')
    return {'strict_orbit_universe': limit, 'all_N_identities': checks,
            'independently_rebuilt_rows': rows, 'row_k_range': [0, 10]}


def packet_audit(max_depth):
    row = {(0, 2): 1}
    raw = dict(row)
    records = []
    old_norm = 1
    for depth in range(1, max_depth + 1):
        row = packet_step(row)
        norm = sum(abs(c) for c in row.values())
        require(norm <= 2 * old_norm, 'normalized norm monotonicity')
        old_norm = norm
        if depth <= 10:
            raw = packet_step(raw, reduce=False)
            for m in range(1, 5):
                expected = sparse_h(1, 2, (1 << depth) * m)
                reduced = (-1)**depth * sum(c * sparse_h(k+1, r, 3**k*m) for (k, r), c in row.items())
                direct = (-1)**depth * sum(c * sparse_h(k+1, r, 3**k*m) for (k, r), c in raw.items())
                require(expected == reduced == direct, 'packet preserves aligned source')
        by_k = defaultdict(lambda: [0, 0, 0, 0])
        for (k, _), c in row.items():
            a = by_k[k]
            a[0] += 1
            a[1] += abs(c)
            a[2] += c*c
            a[3] = max(a[3], abs(c))
        records.append({'depth': depth, 'states': len(row), 'l1': norm,
                        'raw_words': 1 << depth,
                        'by_odd_count': [{'k': k, 'states': v[0], 'l1': v[1], 'l2_square': v[2], 'max_coefficient': v[3]}
                                         for k, v in sorted(by_k.items())]})
    # Every residual and every ordered word through depth 8, no sampled words.
    carry_checks = 0
    for depth in range(1, 9):
        for s in range(1 << depth):
            states = [(0, 2, s)]
            for _ in range(depth):
                nxt = []
                for k, r, value in states:
                    e, o = addresses(k + 1, r)
                    nxt.extend(((k, e, value // 2), (k + 1, o, 3 * ((value + 1) // 2))))
                states = nxt
            for k, r, residual in states:
                require(0 <= residual < 2 * 3**k, 'uniform residual carry')
                for m in (0, 1):
                    change = sparse_h(k+1, r, 3**k*m + residual) - sparse_h(k+1, r, 3**k*m)
                    require(abs(change) <= 4, 'four-unit terminal carry')
                carry_checks += 1
            actual = (-1)**depth * sum(sparse_h(k+1, r, 3**k + residual) for k, r, residual in states)
            require(actual == sparse_h(1, 2, (1 << depth) + s), 'unaligned word source')
    e, o = addresses(1, 2)
    eo = addresses(2, o)[0]
    oe = addresses(1, e)[1]
    require((eo, oe) == (4, 2), 'ordered address hostile')
    depth_hostiles = []
    for q in range(1, 7):
        norms = []
        for r in range(3**q):
            probe = {(q-1, r): 1}
            for _ in range(3):
                probe = packet_step(probe)
            norms.append(sum(abs(c) for c in probe.values()))
        require(set(norms) == ({7} if q == 1 else {8}), 'depth-blind contraction hostile')
        depth_hostiles.append({'q': q, 'all_residues': 3**q, 'three_step_norm': norms[0]})
    # Carry-polynomial formula, including leading even letters and k=0.
    polynomial_checks = 0
    for depth in range(1, 11):
        for code in range(1 << depth):
            q, r, b, k = 3, 17, 0, 0
            for j in range(depth):
                e, o = addresses(q, r)
                if (code >> j) & 1:
                    b, k, q, r = 3*b+(1 << j), k+1, q+1, o
                else:
                    r = e
            mod = 3**q
            require(r == (pow(2, -depth, mod)*(3**k*17+b)) % mod, 'carry polynomial')
            require(b <= 2**(depth-k)*(3**k-2**k), 'fixed-cost carry upper bound')
            polynomial_checks += 1
    return {'rows': records, 'unaligned_word_carry_checks': carry_checks,
            'ordered_address_hostile': {'E_after_O': eo, 'O_after_E': oe},
            'depth_hostiles': depth_hostiles, 'carry_polynomial_checks': polynomial_checks}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--census', type=Path, required=True)
    p.add_argument('--depth', type=int, default=24)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    require(1 <= args.depth <= 24, 'finite packet scope')
    data = args.census.read_bytes()
    census = json.loads(data)
    result = {'status': 'FINITE-EXACT controls; no extrapolation of packet norm',
              'census_sha256': hashlib.sha256(data).hexdigest(),
              'independent_gram': gram_audit(census), 'packets': packet_audit(args.depth)}
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'gram': result['independent_gram'], 'last_packet': result['packets']['rows'][-1],
                      'carry_checks': result['packets']['unaligned_word_carry_checks']}))


if __name__ == '__main__':
    main()
