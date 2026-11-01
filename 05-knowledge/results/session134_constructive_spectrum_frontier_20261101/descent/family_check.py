#!/usr/bin/env python3
"""Actual source/cut and independent adjacency checks for adjacent infinite rays."""
from itertools import permutations
from hashlib import sha256
from pathlib import Path
import json


def require(ok, label):
    if not ok:
        raise RuntimeError(label)


def valid_tournament(rows):
    n = len(rows)
    for u in range(n):
        require(0 <= rows[u] < 1 << n and not(rows[u] >> u & 1), 'row range/loop')
        for v in range(u + 1, n):
            require(((rows[u] >> v) & 1) + ((rows[v] >> u) & 1) == 1, 'one arc per pair')


def path(rows, word):
    return all(rows[u] >> v & 1 for u, v in zip(word, word[1:]))


def literal_h(rows):
    return sum(path(rows, p) for p in permutations(range(len(rows))))


def backward_h(rows):
    """Recursive final-vertex DP, independent of the producer's forward table."""
    from functools import lru_cache
    n = len(rows)
    @lru_cache(None)
    def ending(mask, v):
        old = mask ^ (1 << v)
        if old == 0:
            return 1
        return sum(ending(old, u) for u in range(n)
                   if old >> u & 1 and rows[u] >> v & 1)
    return sum(ending((1 << n) - 1, v) for v in range(n))


def mark(parent, cut):
    n = len(parent)
    return [row | ((1 << n) if not(cut >> u & 1) else 0)
            for u, row in enumerate(parent)] + [cut]


def zero(rows, apex):
    n = len(rows)
    return [row | ((1 << n) if u != apex else 0)
            for u, row in enumerate(rows)] + [1 << apex]


def literal_state(parent, cut):
    h = g = q = 0
    n = len(parent)
    for word in permutations(range(n)):
        if path(parent, word):
            h += 1
            g += (cut >> word[-1]) & 1
        for split in range(n):
            P, Q = word[:split], word[split:]
            if path(parent, P) and path(parent, Q) and cut >> Q[0] & 1:
                q += 1
    return literal_h(mark(parent, cut)), q, g, h


def main():
    parent = [20, 1, 2, 7, 14]
    cut = 17
    source = mark(parent, cut)
    valid_tournament(parent)
    valid_tournament(source)
    require(source == [20, 33, 34, 39, 14, 17], 'source edge convention')
    state = literal_state(parent, cut)
    require(state == (33, 32, 2, 9), 'actual source coordinates')
    require(literal_state([0], 1) == (1, 1, 1, 1), 'plus-one source')
    controls = []
    for k in range(6, 13):
        low = mark([0], 1)
        high = source
        for unused in range(k):
            low = zero(low, 1)
        for unused in range(k - 5):
            high = zero(high, 5)
        valid_tournament(low)
        valid_tournament(high)
        low_H, high_H = backward_h(low), backward_h(high)
        require(low_H == (1 << k) + 1, 'plus-one actual adjacency')
        require(high_H == (1 << k) + 3, 'plus-three actual adjacency')
        require(high_H - low_H == 2 and low_H > 3, 'nonvacuous forbidden nesting pair')
        if k == 6:
            require(literal_h(low) == low_H and literal_h(high) == high_H,
                    'literal first-family controls')
        controls.append({'k': k, 'low_order': len(low), 'high_order': len(high),
                         'low_H': low_H, 'high_H': high_H})
    # r=0 is specifically not a member of the r>=1 count formula.
    require(state[0] != (1 << 5) + 3, 'zero-tail boundary hostile')
    try:
        require(literal_state(parent, cut ^ 1) == state, 'mutated cut rejected')
    except RuntimeError:
        pass
    else:
        raise RuntimeError('cut mutation was not detected')
    print(json.dumps({'status': 'PASS', 'source_parent_rows': parent, 'source_cut': cut,
                      'source_full_rows': source, 'state_order': ['F', 'q', 'g', 'h'],
                      'literal_source_state': state, 'controls': controls,
                      'infinite_claim': 'For every k>=6, both2^k+1 and2^k+3 are actual counts; non-nesting uses THEORY.md.',
                      'tail_scope': 'Literal source plus proved zero-tail recurrence gives all k, not extrapolation from controls.',
                      'hostiles': ['r=0 formula fails as required', 'changed cut rejected'],
                      'source_sha256': sha256(Path(__file__).read_bytes()).hexdigest()},
                     indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
