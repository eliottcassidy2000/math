#!/usr/bin/env python3
"""Literal controls for the compact strong oriented binary constructor."""
from functools import lru_cache
from hashlib import sha256
from itertools import permutations
from pathlib import Path
import json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def binary_graph(k, selected):
    """TT_k -> a -> b, with b -> selected chain vertices, no other arcs."""
    need(k >= 1, 'positive chain length')
    need(all(0 <= v < k for v in selected), 'selected vertex in chain')
    a, b = k, k + 1
    rows = [sum(1 << j for j in range(i + 1, k)) | (1 << a)
            for i in range(k)] + [1 << b, sum(1 << v for v in selected)]
    return tuple(rows)


def construct(h):
    need(isinstance(h, int) and h >= 0, 'nonnegative integer target')
    if h == 0:
        return (0, 0)
    if h == 1:
        return (0,)
    if h == 2:
        return (2, 4, 9, 2)  # inherited Session235 audit/AUDIT.md small witness
    k = (h - 2).bit_length()
    selected = [i for i in range(k) if (h - 2) >> (k - 1 - i) & 1]
    return binary_graph(k, selected)


def hamilton(rows):
    """Independent endpoint subset recurrence, with no construction formula."""
    n = len(rows)
    full = (1 << n) - 1
    @lru_cache(None)
    def walk(mask, v):
        if mask == full:
            return 1
        rest = rows[v] & (full ^ mask)
        total = 0
        while rest:
            bit = rest & -rest
            rest -= bit
            total += walk(mask | bit, bit.bit_length() - 1)
        return total
    return sum(walk(1 << v, v) for v in range(n)) if n else 1


def paths(rows):
    return [p for p in permutations(range(len(rows)))
            if all(rows[u] >> v & 1 for u, v in zip(p, p[1:]))]


def strong(rows):
    n = len(rows)
    for s in range(n):
        seen = todo = 1 << s
        while todo:
            bit = todo & -todo
            todo -= bit
            nxt = rows[bit.bit_length() - 1] & ~seen
            seen |= nxt
            todo |= nxt
        if seen != (1 << n) - 1:
            return False
    return bool(n)


def oriented(rows):
    n = len(rows)
    return (all(0 <= r < (1 << n) and not (r >> i & 1)
                for i, r in enumerate(rows))
            and all(not (rows[i] >> j & rows[j] >> i & 1)
                    for i in range(n) for j in range(i + 1, n)))


def missing(rows):
    return [(i, j) for i in range(len(rows)) for j in range(i + 1, len(rows))
            if not ((rows[i] >> j | rows[j] >> i) & 1)]


def no_zero_encoding(rows):
    """Canonical raw11-empty code: effective absent pairs are written11."""
    raw = list(rows)
    for i, j in missing(rows):
        raw[i] |= 1 << j
        raw[j] |= 1 << i
    return tuple(raw)


def cancel(rows):
    return tuple(sum(1 << j for j in range(len(rows)) if i != j
                     and (rows[i] >> j & 1) and not (rows[j] >> i & 1))
                 for i in range(len(rows)))


def matrix(rows):
    return [''.join(str(r >> j & 1) for j in range(len(rows))) for r in rows]


def main():
    table = []
    permutation_checks = 0
    # Every binary control word through k=6, including nonstrong hostiles.
    for k in range(1, 7):
        hits = []
        for mask in range(1 << k):
            selected = [i for i in range(k) if mask >> i & 1]
            rows = binary_graph(k, selected)
            actual = hamilton(rows)
            predicted = 1 + sum(1 << (k - 1 - i) for i in selected) + (0 in selected)
            need(actual == predicted, 'binary formula')
            need(oriented(rows), 'orientation')
            need(strong(rows) == (0 in selected), 'strongness guard')
            need(len(missing(rows)) == k - len(selected), 'hole count')
            if k <= 5:
                need(len(paths(rows)) == actual, 'permutation cross-check')
                permutation_checks += 1
            if 0 in selected:
                hits.append(actual)
        expected = list(range((1 << (k - 1)) + 2, (1 << k) + 2))
        need(sorted(hits) == expected, 'entire strong interval')
        table.append(dict(k=k, order=k+2, words=1 << k, strong_targets=expected))
    targets = sorted(set(range(131)) | {255, 256, 257, 258, 1023, 1024, 1025, 1026})
    for h in targets:
        rows = construct(h)
        need(hamilton(rows) == h and oriented(rows), 'constructed target')
        need(strong(rows) == (h > 0), 'constructed strongness')
        need(cancel(no_zero_encoding(rows)) == rows, 'semantic coding round trip')
        if h >= 3:
            need(len(rows) == (h - 2).bit_length() + 2, 'logarithmic size')
    witnesses = []
    for h in (0, 1, 2, 7, 21, 1025):
        rows = construct(h)
        raw = no_zero_encoding(rows)
        witnesses.append(dict(h=h, order=len(rows), rows=list(rows), matrix=matrix(rows),
            missing_pairs=missing(rows), strong=strong(rows), no00_raw_matrix=matrix(raw),
            raw_ordinary_H=hamilton(raw), cancelled_H=hamilton(cancel(raw)),
            paths=[list(p) for p in paths(rows)] if h in (7, 21) else None))
    for bad in (-1, 0.5):
        try:
            construct(bad)
        except RuntimeError:
            pass
        else:
            raise RuntimeError('invalid target accepted')
    here = Path(__file__).resolve()
    root = here.parents[4]
    sources = [
        '01-canon/theorems/THM-5071-weighted-cycle-formula-for-digraph-path-totals-and-one-reciprocal-pair-exclusions.md',
        '05-knowledge/results/session235_arc_pair_models_20261109/theory/PROOFS.md',
        '05-knowledge/results/session235_arc_pair_models_20261109/audit/AUDIT.md',
    ]
    print(json.dumps(dict(status='PASS', scope='Finite implementation controls; all-order proof in PROOF.md',
        all_binary_words=sum(x['words'] for x in table), permutation_checks=permutation_checks,
        independent_target_counts=len(targets), targets=targets, interval_table=table,
        witnesses=witnesses, source_sha256={p:sha256((root/p).read_bytes()).hexdigest() for p in sources},
        script_sha256=sha256(here.read_bytes()).hexdigest()), indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
