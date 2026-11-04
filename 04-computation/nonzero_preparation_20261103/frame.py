#!/usr/bin/env python3
"""Literal controls for the fixed-core terminal resource theorem.

All labelled tournaments through five and every nonempty port set are used.
Small path tables come from literal path permutations. The separate growing
zero-channel controls use forward subset DP, without importing prior producers.
"""
from itertools import combinations, permutations
from pathlib import Path
import argparse
import hashlib
import json

CHECKS = 0


def need(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)


def literal_tables(rows):
    n = len(rows)
    h = [0] * (1 << n)
    e = [[0] * n for _ in h]
    s = [[0] * n for _ in h]
    h[0] = 1
    for length in range(1, n + 1):
        for word in permutations(range(n), length):
            if all(rows[a] >> b & 1 for a, b in zip(word, word[1:])):
                mask = sum(1 << a for a in word)
                h[mask] += 1
                e[mask][word[-1]] += 1
                s[mask][word[0]] += 1
    return h, e, s


def dp_tables(rows):
    n = len(rows)
    e = [[0] * n for _ in range(1 << n)]
    s = [[0] * n for _ in e]
    for v in range(n):
        e[1 << v][v] = s[1 << v][v] = 1
    for mask in range(1, 1 << n):
        if mask & (mask - 1) == 0:
            continue
        for v in range(n):
            if mask >> v & 1:
                old = mask ^ (1 << v)
                e[mask][v] = sum(e[old][u] for u in range(n)
                                 if rows[u] >> v & 1)
                s[mask][v] = sum(s[old][u] for u in range(n)
                                 if rows[v] >> u & 1)
    h = [sum(x) for x in e]
    h[0] = 1
    return h, e, s


def signature(tables, ports):
    h, e, s = tables
    full = len(h) - 1
    names = [v for v in range(len(e[0])) if ports >> v & 1]
    ee = [sum(row[v] for v in names) for row in e]
    ss = [sum(row[v] for v in names) for row in s]
    return (h[full], sum(h[a] * h[full ^ a] for a in range(full + 1)),
            sum(ee[a] * h[full ^ a] for a in range(full + 1)),
            sum(h[a] * ss[full ^ a] for a in range(full + 1)),
            sum(ee[a] * ss[full ^ a] for a in range(full + 1)),
            ee[full], ss[full])


def birth(rows, outgoing):
    n = len(rows)
    return tuple(row | ((1 << n) if not (outgoing >> v & 1) else 0)
                 for v, row in enumerate(rows)) + (outgoing,)


def c2(rows, a):
    x = len(rows)
    return birth(birth(rows, 1 << a), 1 << x), x


def join(left, right):
    a, b = len(left), len(right)
    return tuple(row | (((1 << b) - 1) << a) for row in left) + tuple(row << a for row in right)


def attach(rows, ports, k):
    prefix = tuple(sum(1 << j for j in range(i + 1, k)) for i in range(k))
    joined = join(prefix, rows)
    root = len(joined)
    marked = ports << k
    return tuple(row | ((1 << root) if marked >> v & 1 else 0)
                 for v, row in enumerate(joined)) + (((1 << root) - 1) ^ marked,)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=Path(__file__).with_suffix('.json'))
    args = parser.parse_args()
    cases = port_cases = roots = 0
    small = []
    for n in range(1, 6):
        edges = list(combinations(range(n), 2))
        for bits in range(1 << len(edges)):
            rows = [0] * n
            for q, (a, b) in enumerate(edges):
                if bits >> q & 1:
                    rows[a] |= 1 << b
                else:
                    rows[b] |= 1 << a
            tab = literal_tables(rows)
            h, ends, starts = tab
            full = (1 << n) - 1
            need(h[full] % 2 == 1, 'odd Hamilton count')
            for ports in range(1, 1 << n):
                sig = signature(tab, ports)
                need(sig[2] >= ports.bit_count() * sig[0], 'charged cut inequality')
                port_cases += 1
            for v in range(n):
                sig = signature(tab, 1 << v)
                q = h[full ^ (1 << v)]
                delta = sig[2] - sig[0] - q + sig[6]
                degree = sum(row >> v & 1 for row in rows)
                need(q % 2 == 1 and q >= sig[6], 'deletion parity and cut injection')
                need(delta >= degree * (degree - 1) and delta % 2 == 0,
                     'root defect decomposition')
                need(sig[4] == 0, 'singleton joint vanishes')
                if sig[2] == sig[0]:
                    need(sig[6] == q and q % 2 == 1, 'zero-gap start parity')
                roots += 1
            cases += 1
        small.append({'order': n, 'labelled': 1 << len(edges)})
    # The actual terminal chain preserves u=0, while its auxiliary migrates.
    rows, auxiliary = (2, 4, 1), 2
    zero = []
    D = (2, 4, 57, 3, 11, 27)
    need(signature(dp_tables(D), 1 << 4) == (17, 224, 28, 44, 0, 2, 4), 'core')
    for m in range(1, 6):
        sig = signature(dp_tables(rows), 1)
        N = 3 ** m
        need((sig[0], sig[2], sig[4], sig[6]) == (N, N, 0, N // 3), 'zero terminal')
        if m <= 2:
            for k in (1, 2):
                suffix = join(D, rows)
                full = attach(suffix, (1 << 4) | (1 << 6), k)
                obs = signature(dp_tables(full), 1 << (len(full) - 1))
                need((obs[0], obs[2]) == (3 ** (m-1) * (756 * 2 ** k - 83),
                                         3 ** (m-1) * (756 * 2 ** k + 51)), 'actual root response')
        zero.append({'m': m, 'rows': rows, 'auxiliary': auxiliary, 'sigma': sig})
        rows, auxiliary = c2(rows, auxiliary)
    # Exact arithmetic necessary conditions, explicitly not realizability.
    resources = []
    for m in range(1, 9):
        N = 3 ** m
        for digit in (2, 4):
            candidates = []
            for b in range((2 * N - digit) // 180 + 1):
                h, f, s = 3 * N - 8 * b, 3 * N + b, N - 90 * b - digit // 2
                need(h > 0 and f < 2 * h and s >= 0, 'resource and port cardinality')
                parity = (b % 2 == (1 if digit == 2 else 0))
                if parity:
                    candidates.append({'b': b, 'h': h, 'f': f, 's': s, 'F_debt': -136 * b})
            if digit == 2 and m <= 4:
                need(not candidates, 'first four zero depths exclude digit2')
            if digit == 2 and m == 5:
                need(candidates == [{'b': 1, 'h': 721, 'f': 730, 's': 152, 'F_debt': -136}], 'first surviving target')
            resources.append({'m': m, 'digit': digit, 'necessary_candidates': candidates})
    receipt = {'status': 'FINITE-EXACT controls for an independently stated all-order proof',
               'universe': small, 'labelled_tournaments': cases, 'port_sets': port_cases,
               'root_cases': roots, 'checks': CHECKS, 'actual_zero_terminals': zero,
               'arithmetic_not_realizability': resources,
               'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    content = json.dumps(receipt, indent=2, sort_keys=True) + '\n'
    args.output.write_text(content)
    print(json.dumps({'checks': CHECKS, 'labelled_tournaments': cases, 'port_sets': port_cases,
                      'root_cases': roots, 'sha256': hashlib.sha256(content.encode()).hexdigest()}, sort_keys=True))


if __name__ == '__main__':
    main()
