#!/usr/bin/env python3
"""Recompute inherited regular-nine orbit distribution with eight-vertex tests.

Input contract: the pinned Session50 receipt's 85 almost-regular parents are
complete. This checker independently rebuilds all reroot rows and class blocks.
It never canonicalizes an order-nine or order-ten tournament.
"""
from collections import Counter
from itertools import permutations, product
from pathlib import Path
import hashlib
import json
import check

HERE = Path(__file__).resolve().parent
COUNTS = Counter()

def need(test, label):
    COUNTS[label] += 1
    if not test:
        raise RuntimeError(label)


def canonical_parent(rows):
    need(len(rows) == 8, 'canonical input has exactly eight vertices')
    cells = check.refine(rows)
    best = None
    for pieces in product(*(permutations(c) for c in cells)):
        order = sum(pieces, ())
        candidate = tuple(sum(check.arc(rows, v, w) << j for j, w in enumerate(order)) for v in order)
        COUNTS['eight vertex labellings'] += 1
        if best is None or candidate < best:
            best = candidate
    return best


def run():
    raw = (HERE / 'sources/session50-regular-nine-receipt.json').read_bytes()
    pins = json.loads((HERE / 'SOURCE_PINS.json').read_text())
    pin = next(x for x in pins['sources'] if x['archived_path'] == 'sources/session50-regular-nine-receipt.json')
    need(hashlib.sha256(raw).hexdigest() == pin['sha256'], 'inherited regular-sector receipt pin')
    inherited = json.loads(raw)
    parents = [tuple(s['parent']) for s in inherited['states']]
    keys = [canonical_parent(p) for p in parents]
    need(len(keys) == len(set(keys)), 'almost regular parents pairwise nonisomorphic')
    index = {key: i for i, key in enumerate(keys)}
    reroot = []
    for rows in parents:
        check.validate(rows)
        degrees = [r.bit_count() for r in rows]
        need(Counter(degrees) == {3: 4, 4: 4}, 'almost regular degree profile')
        child = list(rows) + [0]
        for v, d in enumerate(degrees):
            if d == 3:
                child[v] |= 1 << 8
            else:
                child[8] |= 1 << v
        need(all(r.bit_count() == 4 for r in child), 'unique extension regular')
        destinations = Counter()
        for deleted in range(9):
            order = tuple(v for v in range(9) if v != deleted)
            smaller = tuple(sum(check.arc(child, v, w) << j for j, w in enumerate(order)) for v in order)
            key = canonical_parent(smaller)
            need(key in index, 'all reroot destinations in inherited parent domain')
            destinations[index[key]] += 1
        reroot.append(destinations)
    remaining = set(range(len(parents)))
    histogram = Counter()
    orbit_histogram = Counter()
    while remaining:
        i = min(remaining)
        component = set(reroot[i])
        need(i in component and component <= remaining, 'regular component partition')
        need(all(reroot[j] == reroot[i] for j in component), 'regular reroot rows identical on component')
        remaining.difference_update(component)
        sizes = tuple(sorted(reroot[i].values()))
        need(sum(sizes) == 9, 'regular orbit sizes sum nine')
        e = sum(9 - s for s in sizes)
        need(e == 9 * (len(component) - 1), 'regular error from orbit count')
        histogram[e] += 1
        orbit_histogram[sizes] += 1
    need(histogram == Counter({0: 3, 18: 3, 36: 2, 72: 7}), 'regular nine independently rebuilt histogram')
    producer = json.loads((HERE / 'normal.json').read_text())
    strong10 = Counter()
    for sector in producer['five_cells']['sectors']:
        for e, count in sector['error_histogram']:
            if e:
                strong10[e] += count
    for record in producer['seven_cell_substitutions']['expansions']:
        if record['error']:
            strong10[record['error']] += 1
    for record in producer['three_triangle_cells']['class_histogram']:
        if record['error']:
            strong10[record['error']] += record['classes']
    need(strong10 == Counter({14: 4, 18: 32, 40: 84}), 'strong ten actual class distribution')
    strong9 = +histogram
    strong9.pop(0, None)
    strong9[14] += 1
    strong9[18] += 4
    strong7 = Counter({r['error']: 1 for r in producer['regular_seven']['classes'] if r['error']})
    total10 = strong10.copy()
    for e, count in strong9.items():
        total10[e] += 2 * count
    for e, count in strong7.items():
        total10[e] += 6 * count
    need(total10 == Counter({14: 12, 18: 46, 36: 4, 40: 84, 72: 14}), 'total ten exact defect distribution')
    need(sum(total10.values()) == 160, 'total ten errorful class count')
    need(sum(e * count for e, count in total10.items()) == producer['total_error_order_ten'], 'distribution first moment')
    return {'status': 'PASS', 'input_contract': 'complete almost-regular-eight parent slice certified in Session50',
            'regular_nine_error_histogram': sorted(histogram.items()),
            'regular_nine_vertex_orbit_multisets': sorted(orbit_histogram.items()),
            'strong_nine_positive_error_histogram': sorted(strong9.items()),
            'strong_ten_positive_error_histogram': sorted(strong10.items()),
            'total_ten_positive_error_histogram': sorted(total10.items()),
            'total_ten_errorful_classes': sum(total10.values()),
            'total_ten_error': sum(e * count for e, count in total10.items()),
            'total_ten_raw_second_moment': sum(e * e * count for e, count in total10.items()),
            'order_nine_or_ten_canonicalizations': 0,
            'checks': dict(sorted(COUNTS.items()))}

if __name__ == '__main__':
    print(json.dumps(run(), indent=2, sort_keys=True))
