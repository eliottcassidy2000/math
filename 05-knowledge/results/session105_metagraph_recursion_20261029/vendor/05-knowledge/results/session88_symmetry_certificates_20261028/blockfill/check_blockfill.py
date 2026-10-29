#!/usr/bin/env python3
"""Matched complete-Q algorithms, retaining every aligned rooted transporter.

All production canonicalization is on deletion parents. Known counts and child
canonicalization are audit-only. Run python3 -B check_blockfill.py --max-order 8.
"""
from collections import Counter
from hashlib import sha256
from itertools import permutations
from pathlib import Path
from statistics import median
from time import perf_counter
import argparse
import json
from baseline_sources import generation as old

PIN = '7a0bab9bfa7bb87c5a10c7ba239592d90cfafd51'
need = old.need


def inverse(p):
    q = [0] * len(p)
    for i, v in enumerate(p):
        q[v] = i
    return tuple(q)


def compose(p, q):
    return tuple(p[v] for v in q)


def aligned_root(engine, child, root):
    """Return complete key and ALL maps canonical pointed child -> child.

    A permutation is stored as destination label for each source label.
    Every map sends the distinguished last vertex to the specified root.
    """
    parent, mask = old.delete_root(child, root)
    code, orders = engine.maps(parent)
    best = min(old.move_mask(mask, p) for p in orders)
    remaining = tuple(v for v in range(len(child)) if v != root)
    maps = tuple(sorted(tuple(remaining[v] for v in p) + (root,)
                        for p in orders if old.move_mask(mask, p) == best))
    return (code, best), maps


def complete_row(engine, child):
    records = tuple(aligned_root(engine, child, v) for v in range(len(child)))
    signature = tuple(sorted(Counter(key for key, maps in records).items()))
    return signature, records


def prepare(parents, n):
    """Same literal cache, raw-mask bank, parent transports in both methods."""
    engine = old.ParentCanonicalizer(max_order=n)
    states = {}
    parent_data = []
    t0 = perf_counter()
    for raw in parents:
        code, incoming = engine.maps(raw)
        rows = old.matrix(n, code)
        code2, group = engine.maps(rows)
        need(code == code2, 'non-idempotent canonical parent')
        parent_data.append((code, raw, incoming, group))
        for mask in range(1 << n):
            key = engine.rooted(rows, mask)
            states.setdefault(key, old.extend(old.matrix(n, key[0]), key[1]))
    return engine, states, tuple(sorted(parent_data)), {
        'parent_preparation_requests': 2 * len(parents),
        'raw_state_classification_requests': len(parents) * (1 << n),
        'preparation_total_requests': engine.requests,
        'preparation_searches': engine.searches,
        'preparation_orders_tested': engine.orders_tested,
        'preparation_seconds': perf_counter() - t0}


def check_row(row, source, states, N):
    need(row and source in dict(row), 'row missing source state')
    need(all(j in states and type(w) is int and w > 0 for j, w in row),
         'row has absent state or nonpositive/noninteger weight')
    need(sum(w for j, w in row) == N, 'partial or incorrect-total row')


def fill_from_anchor(anchor_records, target, N):
    """Transport complete anchored data; no new classification calls."""
    candidates = [p for key, maps in anchor_records if key == target for p in maps]
    need(candidates, 'absent target port')
    to_anchor = min(candidates)
    from_anchor = inverse(to_anchor)
    records = [None] * N
    for v, (key, maps) in enumerate(anchor_records):
        records[from_anchor[v]] = (key, tuple(sorted(compose(from_anchor, p) for p in maps)))
    need(all(record is not None for record in records), 'transport omitted vertex')
    return tuple(records)


def build(parents, n, method):
    start = perf_counter()
    engine, states, parent_data, counts = prepare(parents, n)
    N = n + 1
    rows, roots, groups, representatives = {}, {}, {}, {}
    rows_built = 0
    fill_counts = Counter()
    for source, child in sorted(states.items()):
        if method == 'blockfill' and source in rows:
            continue
        row, records = complete_row(engine, child)
        rows_built += 1
        check_row(row, source, states, N)
        representatives.setdefault(row, child)
        destinations = [source] if method == 'full_reference' else [j for j, w in row]
        for target in destinations:
            need(target not in rows, 'pointed state filled more than once')
            target_records = records if target == source else fill_from_anchor(records, target, N)
            need(target_records[n][0] == target, 'distinguished port misaligned')
            rows[target] = row
            roots[target] = target_records
            groups[target] = tuple(sorted(p for key, maps in target_records if key == target for p in maps))
            fill_counts[target] += 1
    need(set(rows) == set(states), 'pointed states remain uncovered')
    need(set(fill_counts.values()) == {1}, 'incorrect exact coverage')
    # Every operation below is an intrinsic check on the materialized outputs.
    for source, row in rows.items():
        need(all(rows[j] == row for j, w in row), 'block inconsistency')
        need(tuple(sorted(Counter(key for key, maps in roots[source]).items())) == row,
             'root data disagree with Q weights')
        need(len(groups[source]) % 2 == 1, 'even tournament automorphism order')
        need(all(w % 2 == 1 for j, w in row), 'even root-orbit weight')
    counts.update({'method': method, 'order': N, 'pointed_states': len(states),
                   'ordinary_classes': len(representatives), 'rows_classified': rows_built,
                   'reroot_classification_requests': N * rows_built,
                   'all_canonicalizer_requests': engine.requests,
                   'reroot_searches_after_cache': engine.searches - counts['preparation_searches'],
                   'reroot_orders_tested_after_cache': engine.orders_tested - counts['preparation_orders_tested'],
                   'total_searches': engine.searches, 'total_orders_tested': engine.orders_tested,
                   'all_pointed_states_filled_once': True,
                   'production_max_canonicalized_order': n,
                   'seconds': perf_counter() - start})
    return {'states': states, 'Q': rows, 'root_transports': roots,
            'automorphisms': groups, 'parent_data': parent_data,
            'bank': tuple(representatives.values()), 'counts': counts}


def digest(value):
    return sha256(repr(value).encode()).hexdigest()


def compare(a, b):
    hashes = {}
    for field in ('states', 'Q', 'root_transports', 'automorphisms', 'parent_data', 'bank'):
        need(a[field] == b[field], 'full-output mismatch: ' + field)
        value = sorted(a[field].items()) if isinstance(a[field], dict) else a[field]
        hashes[field] = digest(value)
    return hashes


def independent_audit(result, n):
    """Full-child canon is audit only, plus unrestricted labels through six."""
    engine = old.ParentCanonicalizer()
    canon = set()
    orbit_hist = Counter()
    aut_hist = Counter()
    for child in result['bank']:
        code, transports = engine.maps(child)
        need(code not in canon, 'duplicate full-child canonical class')
        canon.add(code)
        # All automorphisms directly from all full-child canonical transports.
        base_inv = inverse(transports[0])
        group = tuple(sorted(compose(p, base_inv) for p in transports))
        key = old.ParentCanonicalizer().rooted(*old.delete_root(child, n - 1))
        need(group == result['automorphisms'][key], 'full-child automorphism mismatch')
        seen = set()
        orbit_weights = []
        for v in range(n):
            if v not in seen:
                orbit = {p[v] for p in group}
                seen.update(orbit)
                orbit_weights.append(len(orbit))
        expected = sorted(w for j, w in result['Q'][key])
        need(sorted(orbit_weights) == expected, 'independent orbit-weight mismatch')
        orbit_hist[tuple(expected)] += 1
        aut_hist[len(group)] += 1
    brute = None
    if n <= 6:
        labelled, brute = old.labelled_orbit_census(n)
        need({engine.maps(t)[0] for t in labelled} == canon, 'labelled universe mismatch')
    return {'independent_full_child_classes': len(canon),
            'automorphism_order_histogram': dict(sorted(aut_hist.items())),
            'orbit_weight_profiles': [{'weights': p, 'classes': c} for p, c in sorted(orbit_hist.items())],
            'labelled_universe': brute}


def hostile_tests():
    def rejects(row, source, states, N):
        try:
            check_row(row, source, states, N)
        except RuntimeError:
            return True
        return False
    states = {'a': None, 'b': None, 'c': None}
    need(rejects((('a', 1), ('b', 1)), 'a', states, 3), 'partial-row hostile undetected')
    need(rejects((('a', 1), ('missing', 2)), 'a', states, 3), 'absent-port hostile undetected')
    # Row-sum, positivity, Q^2=NQ and odd weights do NOT certify actual weights.
    true = ((1, 3), (1, 3))
    fake = ((3, 1), (3, 1))
    square = lambda A: tuple(tuple(sum(A[i][k] * A[k][j] for k in range(2))
                                  for j in range(2)) for i in range(2))
    need(square(true) == tuple(tuple(4*x for x in row) for row in true), 'true square law')
    need(square(fake) == tuple(tuple(4*x for x in row) for row in fake), 'fake square law')
    need(true != fake, 'weight hostile identical')
    # One real tournament realizes the profile (1,3): source -> C3.
    child = old.extend(old.matrix(3, 2), 7)
    E = old.ParentCanonicalizer()
    row, records = complete_row(E, child)
    need(sorted(w for j,w in row) == [1,3], 'real weight hostile not realized')
    # A proper subset of the 3 transitive roots falsely splits one class.
    # Renormalization can hide missing support while preserving all polynomial checks.
    split = ((3,0,0),(0,3,0),(0,0,3))
    need(all(sum(r)==3 for r in split), 'split row sum')
    return {'partial_row_rejected': True, 'absent_port_rejected': True,
            'wrong_odd_weights_pass_row_sum_and_Q_square': {'true': true, 'false': fake},
            'actual_tournament_profile': [1,3],
            'renormalized_missing_support_can_pass_Q_square': split,
            'repair': 'Enumerate every lawful root exactly once and preserve joint parent-mask transports; algebraic identities alone are insufficient.'}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--max-order', type=int, default=8)
    ap.add_argument('--benchmark', type=int, default=0, help='Alternating matched repeats at orders 7 and 8')
    args = ap.parse_args()
    need(1 <= args.max_order <= 8, 'finite scope 1..8')
    parents = [()]
    results, bench = [], []
    for n in range(args.max_order):
        a = build(parents, n, 'full_reference')
        b = build(parents, n, 'blockfill')
        hashes = compare(a,b)
        audit = independent_audit(b, n+1)
        results.append({'order': n+1, 'reference': a['counts'], 'blockfill': b['counts'],
                        'equal_complete_output_sha256': hashes, 'audit': audit})
        if args.benchmark and n >= 6:
            timings = {'full_reference': [], 'blockfill': []}
            preparation = {'full_reference': [], 'blockfill': []}
            for repeat in range(args.benchmark):
                methods = ('full_reference','blockfill') if repeat % 2 == 0 else ('blockfill','full_reference')
                for method in methods:
                    trial = build(parents, n, method)
                    compare(b,trial)
                    timings[method].append(trial['counts']['seconds'])
                    preparation[method].append(trial['counts']['preparation_seconds'])
                    del trial
            medians = {k: median(v) for k,v in timings.items()}
            bench.append({'order': n+1, 'repeat_count': args.benchmark,
                          'seconds': timings, 'preparation_seconds': preparation,
                          'median_seconds': medians,
                          'ratio_of_medians': medians['full_reference']/medians['blockfill']})
        parents = list(b['bank'])
        del a,b
    need([r['blockfill']['ordinary_classes'] for r in results] == [1,1,2,4,12,56,456,6880][:args.max_order],
         'post-construction known-count audit mismatch')
    return {'status': 'FINITE-EXACT', 'source_pin': PIN,
            'reference_file_sha256': sha256(Path(old.__file__).read_bytes()).hexdigest(),
            'scope': 'Ordinary tournament classes, all masks, child orders 1..8, identical complete sparse Q, raw representatives, all pointed transports, all automorphisms, incoming parent transports.',
            'benchmark_scope': 'Both cold per-stage literal caches, identical raw-mask state bank and full outputs. Includes parent preparation, all fills, transport reconstruction and intrinsic validation. Excludes lower levels, equality/hash comparison, independent child/labelled audits. Local timings, no asymptotic claim.',
            'results': results, 'hostiles': hostile_tests(), 'benchmarks': bench}

if __name__ == '__main__':
    print(json.dumps(main(), indent=2, sort_keys=True))
