#!/usr/bin/env python3
"""Small complete target-bounded probes of old-root access, with actual names."""
from math import isqrt
import argparse
import json


def need(test, message):
    if not test:
        raise RuntimeError(message)


def birth(rows, root):
    n = len(rows)
    return tuple(row | ((1 << n) if i == root else 0)
                 for i, row in enumerate(rows)) + (((1 << n)-1) ^ (1 << root),)


def count(rows):
    n = len(rows)
    full = (1 << n)-1
    dp = [[0]*n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1, full+1):
        for v, value in enumerate(dp[mask]):
            if not value:
                continue
            choices = rows[v] & (full ^ mask)
            while choices:
                bit = choices & -choices
                choices -= bit
                dp[mask | bit][bit.bit_length()-1] += value
    return sum(dp[-1])


def prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n)+1))


def probe(limit, degree_cap):
    triangle = (2, 4, 1)
    first = birth(triangle, 0)
    need(count(first) == 5, 'initial root gauge')
    states = [(first, 5, (0,))] if limit >= 5 else []
    hits = {3: (), 5: (0,)} if limit >= 5 else {3: ()}
    layers = []
    tested = 0
    while states:
        n = len(states[0][0])
        next_states = []
        for rows, h, word in states:
            for root in range(n):
                indegree = n-1-rows[root].bit_count()
                if degree_cap is not None and indegree > degree_cap:
                    continue
                child = birth(rows, root)
                value = count(child)
                tested += 1
                need(value >= h+2 and value % 2 == 1, 'strict odd growth')
                if value <= limit:
                    next_word = word+(root,)
                    next_states.append((child, value, next_word))
                    hits.setdefault(value, next_word)
        layers.append({'parent_order': n, 'parents': len(states),
                       'surviving_children': len(next_states)})
        states = next_states
    missing = [h for h in range(3, limit+1, 2) if h not in hits and h not in (7,21)]
    return {'cutoff': limit, 'indegree_cap': degree_cap,
            'normalization': 'first root0; all named later histories retained',
            'tested_children': tested, 'layers': layers,
            'missing_allowed_counts': missing,
            'missing_prime_targets': [h for h in missing if prime(h)],
            'missing_7p_targets': [h for h in missing if h % 7 == 0 and prime(h//7)],
            'witness_histories': {str(h): hits[h] for h in sorted(hits)}}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--limit', type=int, default=99)
    p.add_argument('--caps', default='1,2,all')
    args = p.parse_args()
    need(args.limit >= 5, 'cutoff at least5')
    caps = [None if s == 'all' else int(s) for s in args.caps.split(',')]
    need(all(d is None or d >= 1 for d in caps), 'degree caps at least1')
    output = [probe(args.limit, d) for d in caps]
    print(json.dumps({'scope': 'exact finite cutoff in specified root languages; no tail claim',
                      'probes': output}, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
