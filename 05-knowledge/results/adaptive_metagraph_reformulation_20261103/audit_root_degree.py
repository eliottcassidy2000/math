#!/usr/bin/env python3
"""Independent small root-cap audit: adjacency matrices and recursive counts.

No producer code is imported. No isomorphism quotient or orbit census is used.
"""
import argparse
from collections import Counter
from functools import lru_cache
from itertools import permutations
from pathlib import Path
import hashlib
import json


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def seed():
    return ((False, True, False), (False, False, True), (True, False, False))


def add_vertex(a, root):
    n = len(a)
    require(0 <= root < n, 'root range')
    b = [list(row) + [False] for row in a] + [[False]*(n+1)]
    for v in range(n):
        b[v][n] = (v == root)
        b[n][v] = (v != root)
    return tuple(map(tuple, b))


def indegrees(a):
    return tuple(sum(row[v] for row in a) for v in range(len(a)))


def hcount(a):
    """Backtracking count memoized by current endpoint and unused vertex set."""
    n = len(a)
    @lru_cache(None)
    def walk(last, remaining):
        if not remaining:
            return 1
        total = 0
        for v in range(n):
            if (remaining >> v) & 1 and a[last][v]:
                total += walk(v, remaining ^ (1 << v))
        return total
    full = (1 << n)-1
    return sum(walk(v, full ^ (1 << v)) for v in range(n))


def literal_count(a):
    return sum(all(a[u][v] for u, v in zip(p, p[1:]))
               for p in permutations(range(len(a))))


def is_strong(a):
    n = len(a)
    for u in range(n):
        seen = {u}
        todo = [u]
        while todo:
            v = todo.pop()
            for w in range(n):
                if a[v][w] and w not in seen:
                    seen.add(w)
                    todo.append(w)
        if len(seen) != n:
            return False
    return True


def finite_search(cutoff, cap, first_root=0):
    start = add_vertex(seed(), first_root)
    require(hcount(start) == 5 and is_strong(start), 'first birth')
    states = [(start, 5)]
    hits = {3, 5}
    layers = []
    tested = 0
    checks = Counter()
    while states:
        n = len(states[0][0])
        children = []
        for a, old_h in states:
            ds = indegrees(a)
            for d in range(n+1):
                require(sum(x <= d for x in ds) <= 2*d, 'active bound')
                checks['active_bound'] += 1
            for r in range(n):
                if ds[r] > cap:
                    continue
                b = add_vertex(a, r)
                bs = indegrees(b)
                require(all(bs[v] == ds[v] + (v != r) for v in range(n)),
                        'old indegree update')
                require(bs[-1] == 1, 'new degree')
                require(all(ds[v] <= cap or bs[v] > cap for v in range(n)),
                        'inactive root return')
                checks['degree_update'] += 1
                h = hcount(b)
                require(h >= old_h + 2 and h % 2 == 1, 'strict odd growth')
                tested += 1
                if h <= cutoff:
                    children.append((b, h))
                    hits.add(h)
        layers.append({'parent_order': n, 'parents': len(states),
                       'surviving_children': len(children)})
        states = children
    return {'cutoff': cutoff, 'indegree_cap': cap,
            'first_root': first_root, 'tested_children': tested,
            'layers': layers, 'attained': sorted(hits),
            'missing_allowed_counts': [h for h in range(3, cutoff+1, 2)
                                       if h not in hits and h not in (7, 21)],
            'controls': dict(checks)}


def replay(word, cap):
    a = seed()
    ds = []
    for r in word:
        d = indegrees(a)[r]
        require(d <= cap, 'illegal witness cap')
        ds.append(d)
        a = add_vertex(a, r)
    return a, tuple(ds)


def verify_receipt(path):
    data = json.loads(path.read_text())
    rows = []
    for rec in data['probes']:
        cap = rec['indegree_cap']
        found = set()
        for target, word in rec['witness_histories'].items():
            a, ds = replay(word, cap if cap is not None else 10**9)
            value = hcount(a)
            require(value == int(target), 'receipt witness count')
            require(value <= rec['cutoff'] and is_strong(a), 'receipt source')
            found.add(value)
        misses = [v for v in range(3, rec['cutoff']+1, 2)
                  if v not in found and v not in (7,21)]
        require(misses == rec['missing_allowed_counts'], 'receipt missing list')
        rows.append({'cutoff': rec['cutoff'], 'cap': cap,
                     'witnesses_checked': len(found), 'missing': misses})
    return {'file': path.name, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
            'probes': rows}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--receipt-dir', type=Path, default=Path(__file__).parent)
    args = parser.parse_args()
    # Explicit root-rotation isomorphisms prove first-root normalization.
    root_iso = []
    base = add_vertex(seed(), 0)
    for r in range(3):
        other = add_vertex(seed(), r)
        perm = tuple((v+r) % 3 for v in range(3)) + (3,)
        require(all(base[u][v] == other[perm[u]][perm[v]]
                    for u in range(4) for v in range(4)), 'first root rotation')
        root_iso.append({'root': r, 'permutation': perm})
    word = (0, 0, 3, 3, 0)
    a, ds = replay(word, 3)
    rows = tuple(sum((1 << v) for v, bit in enumerate(row) if bit) for row in a)
    require(ds == (1, 1, 2, 2, 3), '127 degree trace')
    require(rows == (154, 4, 1, 102, 14, 23, 55, 126), '127 adjacency')
    literal = literal_count(a)
    require(literal == hcount(a) == 127 and is_strong(a), '127 literal count')
    # Repeat the complete cutoff199 cap2 search independently.
    search199 = finite_search(199, 2)
    require(search199['missing_allowed_counts'] == [127], 'cap2 boundary')
    require(7 not in search199['attained'] and 21 not in search199['attained'], 'known gap controls')
    require(all(127 % d for d in range(2, 12)), 'prime127')
    # Also use every first root at cutoff127: no reliance on the gauge for this hostile.
    gauge = [finite_search(127, 2, r) for r in range(3)]
    require(all(x['missing_allowed_counts'] == [127] for x in gauge), 'gauge hostile')
    require(all(x['layers'] == gauge[0]['layers'] for x in gauge), 'gauge layers')
    prod199 = json.loads((args.receipt_dir/'root_degree_199.json').read_text())['probes'][0]
    for key in ['layers','tested_children','missing_allowed_counts']:
        require(search199[key] == prod199[key], 'producer discrepancy: '+key)
    receipts = [verify_receipt(args.receipt_dir/f) for f in
                ['root_degree_99.json','root_degree_199.json','root_degree_3_199.json']]
    # Positive coverage is fully certified by actual witnesses, not by search completeness.
    cap3 = next(x for x in receipts[-1]['probes'] if x['cap'] == 3)
    require(cap3['missing'] == [] and cap3['witnesses_checked'] == 97,
            'cap3 finite coverage')
    print(json.dumps({'status': 'PASS',
       'scope': 'Independent cap2 cutoff199 search; all first-root gauges cutoff127; '
                'all saved witnesses verified recursively; literal8! count127. '
                'No all-target cap3 or unrestricted completeness claim.',
       'first_root_isomorphisms': root_iso,
       'witness127': {'history': word, 'root_indegrees_before_birth': ds,
                      'rows': rows, 'literal_permutations_examined': 40320,
                      'H': literal, 'strong': True},
       'cap2_cutoff199': search199,
       'cap2_cutoff127_all_first_roots': gauge,
       'receipt_witness_audits': receipts}, indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
