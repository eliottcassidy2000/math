"""Finite controls for the independent all-order proof; standard library only."""
from algorithm import *
from collections import Counter
from itertools import permutations
import json


def need(value, message):
    if not value:
        raise RuntimeError(message)


def literal(rows):
    return sum(all(rows[u] >> v & 1 for u, v in zip(p, p[1:]))
               for p in permutations(range(len(rows))))


def components(rows):
    n = len(rows)
    reach = [r | (1 << i) for i, r in enumerate(rows)]
    for k in range(n):
        for i in range(n):
            if reach[i] >> k & 1:
                reach[i] |= reach[k]
    unused = set(range(n))
    blocks = []
    while unused:
        v = min(unused)
        block = sorted(u for u in unused if reach[v] >> u & 1 and reach[u] >> v & 1)
        unused.difference_update(block)
        blocks.append(block)
    return blocks


def induce(rows, vertices):
    return tuple(sum(1 << j for j, w in enumerate(vertices) if rows[v] >> w & 1)
                 for v in vertices)


def main():
    total = literals = reductions = caps = 0
    strong_min = {}; reduced_min = {}; bank = {}; strong_hist = {}
    for n in range(1, 7):
        counts = Counter()
        for rows in tournaments(n):
            total += 1
            h = hp(rows)
            if n <= 5:
                need(h == literal(rows), 'literal path counter')
                literals += 1
                for cap in (2, 4, 6, 10, 24):
                    need(hp(rows, cap) == min(cap, h), 'saturation')
                    caps += 1
            cc = components(rows)
            product = 1
            kept = []
            for block in cc:
                product *= hp(induce(rows, block))
                if len(block) > 1:
                    kept.extend(block)
            need(product == h, 'SCC product')
            if h > 1:
                core = induce(rows, kept)
                need(hp(core) == h and 3 ** len(core) <= h ** 3, 'reduced bound')
                reductions += 1
            else:
                need(not kept, 'unit count reduction boundary')
            if len(kept) == n:
                reduced_min[n] = min(reduced_min.get(n, h), h)
            if len(cc) == 1 and n >= 3:
                need(is_strong(rows), 'independent SCC agreement')
                counts[h] += 1
                strong_min[n] = min(strong_min.get(n, h), h)
                if h <= 23:
                    bank.setdefault(h, rows)
        strong_hist[n] = dict(sorted(counts.items()))
    for n, value in strong_min.items():
        need(value == strong_floor(n), 'small strong minimum')
    for n, value in reduced_min.items():
        need(value == reduced_floor(n), 'small reduced minimum')
    # Independent multiplicative partition DP, allowing every published strong size.
    best = [None] * 201
    best[0] = 1
    for n in range(3, 201):
        values = [best[n-s] * strong_floor(s) for s in range(3, n+1)
                  if best[n-s] is not None]
        best[n] = min(values)
        need(best[n] == reduced_floor(n), 'exact reduced partition formula')
        need(strong_floor(n) ** 3 >= 5 ** (n-1), 'published floor simplification')
        need(strong_floor(n) ** 3 >= 3 ** n, 'per-component order inequality')
    rounding = 0
    for h in range(1, 10001):
        b = reduced_bound(h)
        need(3 ** b <= h ** 3 < 3 ** (b+1), 'integer rounding')
        rounding += 1
    decisions = {}
    for h in range(-1, 24):
        rows = from_bank(h, bank)
        expected = h > 0 and h % 2 == 1 and h not in (7, 21)
        need((rows is not None) == expected, 'small complete target answer')
        if rows is not None:
            need(hp(rows) == h, 'returned physical witness')
            need(h == 1 or len(rows) <= reduced_bound(h), 'witness size')
        decisions[h] = None if rows is None else list(rows)
    need(strong_bound(23) == 6, 'complete bank order bound')
    need(reduced_bound(5) == 4, 'floor location counterexample')
    need(reduced_bound(1) == 0 and decide(1) == (0,), 'unit boundary')
    need(decide(7) is None, 'standalone target invocation')
    return dict(status='PASS', tournaments=total, literal_counter_checks=literals,
                saturation_checks=caps, nonunit_reduction_checks=reductions,
                strong_minima=strong_min, reduced_minima=reduced_min,
                partition_orders_checked=200, integer_rounding_checks=rounding,
                complete_targets_through23=decisions, strong_histograms=strong_hist)


if __name__ == '__main__':
    print(json.dumps(main(), indent=2, sort_keys=True))
