#!/usr/bin/env python3
"""Exact parent-only canonical augmentation. Standard library; explicit checks."""
from collections import Counter
from itertools import permutations, product
import argparse
import json


def need(x, message):
    if not x:
        raise RuntimeError(message)


def matrix(n, code):
    rows = [0] * n
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (code >> k) & 1:
                rows[i] |= 1 << j
            else:
                rows[j] |= 1 << i
            k += 1
    return tuple(rows)


def encode(rows, order):
    value = 0
    k = 0
    for i, u in enumerate(order):
        for v in order[i + 1:]:
            value |= ((rows[u] >> v) & 1) << k
            k += 1
    return value


def move_mask(mask, order):
    return sum(((mask >> v) & 1) << i for i, v in enumerate(order))


def degree_orders(rows):
    degrees = [r.bit_count() for r in rows]
    groups = [[v for v, d in enumerate(degrees) if d == k] for k in sorted(set(degrees))]
    for choices in product(*(permutations(g) for g in groups)):
        yield tuple(v for block in choices for v in block)


class ParentCanonicalizer:
    """Caches full sets of attaining transports, never unaligned marginals."""
    def __init__(self, max_order=None):
        self.max_order = max_order
        self.cache = {}
        self.requests = 0
        self.searches = 0
        self.orders_tested = 0

    def maps(self, rows):
        need(self.max_order is None or len(rows) <= self.max_order,
             'production attempted child-size canonicalization')
        self.requests += 1
        if rows in self.cache:
            return self.cache[rows]
        self.searches += 1
        best = None
        transports = []
        for order in degree_orders(rows):
            self.orders_tested += 1
            code = encode(rows, order)
            if best is None or code < best:
                best = code
                transports = [order]
            elif code == best:
                transports.append(order)
        result = (best, tuple(transports))
        self.cache[rows] = result
        return result

    def rooted(self, rows, mask):
        code, transports = self.maps(rows)
        return code, min(move_mask(mask, p) for p in transports)


def extend(rows, mask):
    n = len(rows)
    return tuple(row | (((mask >> i & 1) ^ 1) << n) for i, row in enumerate(rows)) + (mask,)


def delete_root(rows, root):
    old = tuple(i for i in range(len(rows)) if i != root)
    return tuple(move_mask(rows[v], old) for v in old), move_mask(rows[root], old)


def orbit_masks(n, group):
    seen = set()
    result = []
    for mask in range(1 << n):
        if mask in seen:
            continue
        orbit = {move_mask(mask, p) for p in group}
        need(min(orbit) == mask, 'orbit traversal skipped an earlier representative')
        seen.update(orbit)
        result.append(mask)
    need(len(seen) == 1 << n, 'incomplete mask orbit coverage')
    return result


def produce(parents, n, audit=True, selector="degree"):
    """Produce size n+1 raw matrices using only <=n-vertex canon searches."""
    engine = ParentCanonicalizer(max_order=n)
    normalized = []
    for raw in parents:
        code, _ = engine.maps(raw)
        rows = matrix(n, code)
        canonical_code, group = engine.maps(rows)
        need(code == canonical_code, 'normalization is not idempotent')
        need(all(encode(rows, p) == code for p in group), 'not an actual automorphism')
        normalized.append((code, rows, group))
    normalized.sort()
    preparation_requests = engine.requests
    preparation_searches = engine.searches
    preparation_orders = engine.orders_tested
    count = Counter()
    tie_hostile = None
    children = []
    pointed = []
    accepted_keys = []
    child_audit = ParentCanonicalizer()  # separate, never supplies production data
    audit_all_codes = set()
    audit_accepted_codes = []
    for code, rows, group in normalized:
        masks = orbit_masks(n, group)
        count['raw_masks'] += 1 << n
        count['pointed_mask_orbits'] += len(masks)
        for mask in masks:
            child = extend(rows, mask)
            key = (code, mask)
            pointed.append((key, child))
            if audit:
                audit_all_codes.add(child_audit.maps(child)[0])
            degree = [r.bit_count() for r in child]
            if selector == 'degree':
                signatures = degree
            elif selector == 'neighbor_degrees':
                signatures = [(degree[v], tuple(sorted(degree[w] for w in range(n + 1)
                                      if child[v] >> w & 1)),
                               tuple(sorted(degree[w] for w in range(n + 1)
                                      if w != v and not (child[v] >> w & 1))))
                              for v in range(n + 1)]
            else:
                raise RuntimeError('unknown root selector')
            if signatures[n] != max(signatures):
                count['selector_rejected'] += 1
                continue
            count['selector_survivors'] += 1
            accepted = True
            # Every compared root belongs to the invariant maximum-degree locus.
            for v in range(n):
                if signatures[v] != signatures[n]:
                    continue
                parent, attachment = delete_root(child, v)
                count['competitor_key_queries'] += 1
                competitor = engine.rooted(parent, attachment)
                if competitor < key:
                    count['key_rejected'] += 1
                    if tie_hostile is None:
                        tie_hostile = {'raw_child_code': encode(child, tuple(range(n + 1))),
                                       'root': n, 'competitor': v,
                                       'equal_selected_signature': signatures[n],
                                       'rooted_key': key, 'smaller_competitor_key': competitor}
                    accepted = False
                    break
            if accepted:
                children.append(child)
                accepted_keys.append(key)
                if audit:
                    audit_accepted_codes.append(child_audit.maps(child)[0])
    need(len(children) == count['selector_survivors'] - count['key_rejected'], 'acceptance accounting')
    if audit:
        need(len(set(audit_accepted_codes)) == len(children), 'duplicate output class')
        need(set(audit_accepted_codes) == audit_all_codes, 'lost output class')
        # Full Q is an independently checked reference, not needed for generation.
        state_keys = {key for key, _ in pointed}
        reference_engine = ParentCanonicalizer()
        support = {}
        for key, child in pointed:
            row = Counter(reference_engine.rooted(*delete_root(child, v)) for v in range(n + 1))
            need(set(row) <= state_keys, 'rerooted state missing')
            signature = tuple(sorted(row.items()))
            support[key] = signature
        need(len(set(support.values())) == len(children), 'root-switch row count mismatch')
        for key, signature in support.items():
            need(sum(m for _, m in signature) == n + 1, 'Q row sum')
            need(all(support[j] == signature for j, _ in signature), 'Q block inconsistency')
        need(len({support[k] for k in accepted_keys}) == len(children), 'not one accepted root per Q block')
    report = {'order': n + 1, 'parent_classes': len(parents), 'root_selector': selector, **dict(sorted(count.items())),
              'accepted_classes': len(children),
              'reference_full_Q_key_queries': (n + 1) * len(pointed),
              'reference_raw_state_bank_key_queries': (1 << n) * len(parents),
              'parent_preparation_requests': preparation_requests,
              'parent_preparation_searches': preparation_searches,
              'parent_preparation_orders_tested': preparation_orders,
              'competitor_canon_searches_after_cache': engine.searches - preparation_searches,
              'competitor_orders_tested_after_cache': engine.orders_tested - preparation_orders,
              'production_max_canonicalized_order': n,
              'audit_full_Q_and_child_canon': audit,
              'equal_signature_distinct_root_orbit_hostile': tie_hostile}
    return children, report


def labelled_orbit_census(n):
    """Independent brute-force labelled universe + full S_n orbit marking."""
    universe = 1 << (n * (n - 1) // 2)
    unseen = set(range(universe))
    group = tuple(permutations(range(n)))
    histogram = Counter()
    reps = []
    while unseen:
        code = min(unseen)
        rows = matrix(n, code)
        orbit = {encode(rows, p) for p in group}
        need(orbit <= unseen, 'S_n orbits overlap')
        unseen.difference_update(orbit)
        histogram[len(orbit)] += 1
        reps.append(rows)
    need(sum(s * m for s, m in histogram.items()) == universe, 'labelled mass mismatch')
    return reps, {'order': n, 'labelled_universe': universe, 'class_count': len(reps),
                  'labelled_orbit_size_histogram': dict(sorted(histogram.items()))}


def permutation_cycles(p):
    unseen = set(range(len(p)))
    sizes = []
    while unseen:
        start = min(unseen)
        v = start
        length = 0
        while v in unseen:
            unseen.remove(v)
            length += 1
            v = p[v]
        need(v == start, 'not a permutation cycle')
        sizes.append(length)
    return tuple(sorted(sizes))


def blowup_hostile():
    # T=C3[A0,A0,A0], U=C3[A0,A0,A1], A1 literal converse of A0.
    def build(reversed_last):
        rows = [0] * 12
        for u in range(12):
            for v in range(u + 1, 12):
                b, q = divmod(u, 4)
                c, r = divmod(v, 4)
                if b != c:
                    arc = (c - b) % 3 == 1
                else:
                    arc = q == 0 or (q > 0 and (r - q) % 3 == 1)
                    if reversed_last and b == 2:
                        arc = not arc
                if arc:
                    rows[u] |= 1 << v
                else:
                    rows[v] |= 1 << u
        return tuple(rows)
    def explicit_group(outer):
        group = []
        for k in range(3 if outer else 1):
            for shifts in product(range(3), repeat=3):
                p = []
                for b in range(3):
                    for q in range(4):
                        target_q = 0 if q == 0 else 1 + (q - 1 + shifts[b]) % 3
                        p.append(4 * ((b + k) % 3) + target_q)
                group.append(tuple(p))
        return tuple(group)
    T, U = build(False), build(True)
    G, H = explicit_group(True), explicit_group(False)
    for rows, group in ((T, G), (U, H)):
        code = encode(rows, tuple(range(12)))
        need(all(encode(rows, p) == code for p in group), 'hostile generator not arc preserving')
    cyclic = lambda rows: {triple for triple in __import__('itertools').combinations(range(12), 3)
                          if all(sum((rows[v] >> w) & 1 for w in triple) == 1 for v in triple)}
    need(cyclic(T) == cyclic(U) and len(cyclic(T)) == 67, 'cyclic-triple hostile failed')
    need(all({p for p in G if p[v] == v} == {p for p in H if p[v] == v}
             for v in range(12)), 'literal point stabilizer mismatch')
    outputs = []
    for group in (G, H):
        masks = orbit_masks(12, group)
        burnside = sum(1 << len(permutation_cycles(p)) for p in group)
        need(burnside % len(group) == 0 and burnside // len(group) == len(masks), 'mask Burnside disagreement')
        outputs.append({'group_order': len(group), 'binary_attachment_orbits': len(masks),
                        'single_marked_port_orbits': sum(m.bit_count() == 1 for m in masks),
                        'transporter_0_to_4_size': sum(p[0] == 4 for p in group)})
    # Independent DP checks all named induced-subset H readouts agree.
    def subset_H(rows):
        N = 1 << len(rows)
        dp = [[0] * len(rows) for _ in range(N)]
        values = [0] * N
        values[0] = 1
        for v in range(len(rows)):
            dp[1 << v][v] = 1
        for mask in range(1, N):
            for end, ways in enumerate(dp[mask]):
                if not ways:
                    continue
                rest = rows[end] & ~mask
                while rest:
                    bit = rest & -rest
                    rest -= bit
                    dp[mask | bit][bit.bit_length() - 1] += ways
            values[mask] = sum(dp[mask])
        return values
    ht, hu = subset_H(T), subset_H(U)
    need(ht == hu, 'induced H hostile mismatch')
    return {'same_labelled_cyclic_triples': 67, 'same_all_subset_H_entries': len(ht),
            'whole_H': ht[-1], 'literal_point_stabilizers_equal': True,
            'attachments': outputs}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--max-order', type=int, default=7)
    parser.add_argument('--brute-order', type=int, default=6)
    parser.add_argument('--selector', choices=['degree', 'neighbor_degrees'], default='degree')
    parser.add_argument('--production-only', action='store_true',
                        help='Disable every child-canon/full-Q/labelled-universe audit.')
    args = parser.parse_args()
    need(1 <= args.max_order <= 8, 'explicit finite scope is orders 1..8')
    need(1 <= args.brute_order <= 6, 'independent labelled universe is limited to 6')
    parents = [()]
    table = []
    banks = {0: parents}
    for n in range(args.max_order):
        parents, row = produce(parents, n, audit=not args.production_only, selector=args.selector)
        table.append(row)
        banks[n + 1] = parents
    known = [1, 1, 2, 4, 12, 56, 456, 6880]
    need([r['accepted_classes'] for r in table] == known[:args.max_order], 'known counts audit mismatch')
    brute_rows = []
    for n in range(1, 1 if args.production_only else min(args.brute_order, args.max_order) + 1):
        brute, row = labelled_orbit_census(n)
        audit = ParentCanonicalizer()
        need({audit.maps(t)[0] for t in brute} == {audit.maps(t)[0] for t in banks[n]},
             'independent labelled-orbit oracle mismatch')
        brute_rows.append(row)
    return {'status': 'FINITE-EXACT', 'source_pin': 'b8672af00bcc6d440a9a00a6915325faf39fcdff',
            'generation': table, 'independent_labelled_orbit_census': brute_rows,
            'transport_hostile': blowup_hostile(),
            'production_uses_full_child_canon': False,
            'check_policy': 'Explicit exceptions, not removable assertions; all known counts are audit-only.'}

if __name__ == '__main__':
    print(json.dumps(main(), indent=2, sort_keys=True))
