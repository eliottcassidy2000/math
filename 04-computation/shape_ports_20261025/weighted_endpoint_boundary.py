#!/usr/bin/env python3
"""Exact first weighted-endpoint branch and pointed-C3 collapse.

No repository imports. Full weighted class quotient at k=2..6; all matching
span-two deletion/insertion checks through k=7; subset-DP feedback controls.
"""
from collections import defaultdict
from itertools import permutations, product
from math import prod
from pathlib import Path
import argparse
import hashlib
import json


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def matchings(vertices):
    if not vertices:
        yield ()
        return
    a = vertices[0]
    for i, b in enumerate(vertices[1:], 1):
        for tail in matchings(vertices[1:i] + vertices[i + 1:]):
            yield ((a, b),) + tail


def irreducible(edges, size):
    return all(any(a < cut <= b for a, b in edges) for cut in range(1, size))


def matching_count(k):
    return prod(range(1, 2 * k, 2))


def strong_counts(max_k):
    values = [0] * (max_k + 1)
    for k in range(1, max_k + 1):
        values[k] = matching_count(k) - sum(values[j] * matching_count(k - j) for j in range(1, k))
    return values


def build(edges, endpoint_weights, gap_weights, strong):
    count = len(endpoint_weights)
    order = []
    for i in range(count + 1):
        if gap_weights[i] and not (strong and i in (0, count)):
            order.append(('g', i, gap_weights[i]))
        if i < count:
            order.append(('e', i, endpoint_weights[i]))
    n = len(order)
    positions = {index: pos for pos, (kind, index, _) in enumerate(order) if kind == 'e'}
    reverse = {tuple(sorted((positions[a], positions[b]))) for a, b in edges}
    rows = [0] * n
    for a in range(n):
        for b in range(a + 1, n):
            if (a, b) in reverse:
                rows[b] |= 1 << a
            else:
                rows[a] |= 1 << b
    return tuple(rows), tuple(w for _, _, w in order), order


def has_clone(rows):
    full = (1 << len(rows)) - 1
    return any(not ((rows[a] ^ rows[b]) & (full ^ (1 << a) ^ (1 << b)))
               for a in range(len(rows)) for b in range(a + 1, len(rows)))


def strong_graph(rows):
    n = len(rows)
    if n < 2:
        return False
    reach = [row | (1 << i) for i, row in enumerate(rows)]
    for k in range(n):
        for i in range(n):
            if reach[i] >> k & 1:
                reach[i] |= reach[k]
    return all(row == (1 << n) - 1 for row in reach)


def canonical(rows, weights):
    """Weight/directed-colour refinement, then every unresolved permutation.

    This does not use ordered gaps, matchings, distinguished endpoints or
    the theorem's anchor decoder. Thus it can detect the k=2 owner collapse.
    """
    groups = defaultdict(list)
    for v, weight in enumerate(weights):
        groups[weight].append(v)
    cells = [groups[key] for key in sorted(groups)]
    while True:
        masks = [sum(1 << v for v in cell) for cell in cells]
        new = []
        for cell in cells:
            parts = defaultdict(list)
            for v in cell:
                signature = tuple((rows[v] & mask).bit_count() for mask in masks)
                parts[signature].append(v)
            new.extend(parts[key] for key in sorted(parts))
        if len(new) == len(cells):
            cells = new
            break
        cells = new
    best = None
    for blocks in product(*(permutations(cell) for cell in cells)):
        order = tuple(v for block in blocks for v in block)
        candidate = (tuple(weights[v] for v in order),
                     tuple(sum(1 << j for j, w in enumerate(order) if rows[v] >> w & 1) for v in order))
        if best is None or candidate < best:
            best = candidate
    return best


def feedback(rows, weights):
    n = len(rows)
    full = (1 << n) - 1
    incoming = [sum(1 << v for v in range(n) if rows[v] >> u & 1) for u in range(n)]
    weight_sum = [0] * (full + 1)
    for mask in range(1, full + 1):
        bit = mask & -mask
        weight_sum[mask] = weight_sum[mask ^ bit] + weights[bit.bit_length() - 1]
    dp = [0] * (full + 1)
    for mask in range(1, full + 1):
        best = 10 ** 20
        remaining = mask
        while remaining:
            bit = remaining & -remaining
            remaining ^= bit
            v = bit.bit_length() - 1
            # Place v first; the remaining incoming arcs are backward.
            best = min(best, dp[mask ^ bit] + weights[v] * weight_sum[incoming[v] & mask])
        dp[mask] = best
    return dp[-1]


def insert_span_two(edges, mark, size):
    def lift(v):
        return v if v < mark else v + 1 if v == mark else v + 2
    return tuple(sorted(tuple(sorted((lift(a), lift(b)))) for a, b in edges)) + ((mark, mark + 2),)


def delete_span_two(edges, marked_edge, size):
    a, b = marked_edge
    check(b == a + 2, 'not span two')
    remaining = [v for v in range(size) if v not in (a, b)]
    index = {v: i for i, v in enumerate(remaining)}
    reduced = tuple(sorted((index[x], index[y]) for x, y in edges if (x, y) != marked_edge))
    return reduced, index[a + 1]


def matching_bijection(max_k):
    indecomposable = strong_counts(max_k)
    output = []
    for k in range(2, max_k + 1):
        total = strong = 0
        for edges in matchings(tuple(range(2 * k))):
            original_strong = irreducible(edges, 2 * k)
            for pair in edges:
                if pair[1] - pair[0] != 2:
                    continue
                reduced, mark = delete_span_two(edges, pair, 2 * k)
                check(tuple(sorted(insert_span_two(reduced, mark, 2 * k - 2))) == tuple(sorted(edges)),
                      'span-two map not inverse')
                check(original_strong == irreducible(reduced, 2 * k - 2), 'cut irreducibility not preserved')
                total += 1
                strong += original_strong
        check(total == 2 * (k - 1) * matching_count(k - 1), 'H formula')
        check(strong == 2 * (k - 1) * indecomposable[k - 1], 'H strong formula')
        output.append({'k': k, 'all_marked_span_two': total, 'strong_marked_span_two': strong})
    return output


def boundary_classes(max_k):
    indecomposable = strong_counts(max_k)
    output = []
    dp_controls = 0
    relabel_controls = 0
    for k in range(2, max_k + 1):
        endpoint_count = 2 * k - 2
        for strong in (False, True):
            keys, owner_keys = set(), set()
            records = 0
            for edges in matchings(tuple(range(endpoint_count))):
                if strong and not irreducible(edges, endpoint_count):
                    continue
                for mark in range(endpoint_count):
                    endpoint_weights = [1] * endpoint_count
                    endpoint_weights[mark] = 2
                    rows, weights, order = build(edges, endpoint_weights, [k] * (endpoint_count + 1), strong)
                    check(not has_clone(rows), 'target clone')
                    check(not strong or strong_graph(rows), 'target not strong')
                    check(sum(weights) == (2 * k * k + (k if strong else 3 * k)) - (2 * k + 1), 'target mass')
                    key = canonical(rows, weights)
                    keys.add(key)
                    if records < 3:
                        reversed_rows = tuple(sum(1 << (len(rows) - 1 - j) for j in range(len(rows)) if rows[i] >> j & 1)
                                              for i in range(len(rows) - 1, -1, -1))
                        check(canonical(reversed_rows, tuple(reversed(weights))) == key, 'relabel covariance')
                        relabel_controls += 1
                    if k <= 4 or (k == 5 and records < 2):
                        check(feedback(rows, weights) == k, 'wrong exact feedback')
                        dp_controls += 1

                    # Reconstruct the previous pointed-C3 owner, retaining its port.
                    lifted = insert_span_two(edges, mark, endpoint_count)
                    source_gaps = [k] * (endpoint_count + 3)
                    source_gaps[mark + 1] = source_gaps[mark + 2] = 0
                    old_rows, old_weights, old_order = build(lifted, [1] * (endpoint_count + 2), source_gaps, strong)
                    check(not has_clone(old_rows), 'owner source clone')
                    check(not strong or strong_graph(old_rows), 'owner source not strong')
                    check(sum(old_weights) == sum(weights) + 1, 'collapse did not lose mass one')
                    check(len(old_rows) == len(rows) + 2, 'collapse did not lose two vertices')
                    owner_keys.add(canonical(old_rows, old_weights))
                    deleted = {i for i, item in enumerate(old_order) if item[0] == 'e' and item[1] in (mark, mark + 2)}
                    port = next(i for i, item in enumerate(old_order) if item[:2] == ('e', mark + 1))
                    survivors = [i for i in range(len(old_rows)) if i not in deleted]
                    collapsed_rows = tuple(sum(1 << j for j, v in enumerate(survivors) if old_rows[u] >> v & 1)
                                           for u in survivors)
                    collapsed_weights = tuple(2 if u == port else old_weights[u] for u in survivors)
                    check(collapsed_rows == rows and collapsed_weights == weights, 'literal port collapse failed')
                    records += 1
            expected_records = 2 * (k - 1) * (indecomposable[k - 1] if strong else matching_count(k - 1))
            check(records == expected_records, 'record count')
            check(len(keys) == (1 if k == 2 else records), 'weighted target owner count')
            check(len(owner_keys) == records, 'old pointed-owner records unexpectedly collide')
            output.append({'k': k, 'strong': strong, 'records': records,
                           'target_classes': len(keys), 'source_owner_classes': len(owner_keys),
                           'mass': (2 * k * k + (k if strong else 3 * k)) - (2 * k + 1)})
    return output, dp_controls, relabel_controls


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    bijection = matching_bijection(7)
    classes, dp_controls, relabel_controls = boundary_classes(6)
    result = {'status': 'PROVED scoped family; FINITE-EXACT controls',
              'matching_bijection': bijection, 'boundary_classes': classes,
              'exact_weighted_feedback_controls': dp_controls, 'relabel_controls': relabel_controls,
              'scope': 'Only r<k at mass deficit 2k+1; not the full mass layer or A000568'}
    payload = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(payload)
    print(json.dumps({'classes': classes, 'feedback_controls': dp_controls,
                      'sha256': hashlib.sha256(payload.encode()).hexdigest()}, sort_keys=True))


if __name__ == '__main__':
    main()
