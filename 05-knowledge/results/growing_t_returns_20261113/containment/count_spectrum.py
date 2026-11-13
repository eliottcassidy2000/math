"""All induced count spectra of the complete first five ready-target heads.

Two independent subset recurrences and literal small controls.  The head's
all-order completeness is an inherited mathematical input, not inferred here.
"""
from pathlib import Path
from collections import Counter
from itertools import permutations
import argparse
import gzip
import hashlib
import json
import sys
import numpy as np
from numba import njit

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent.parent
sys.path.insert(0, str(RESULTS / 'berggren_target_return_20261111' / 'reset'))
import funnel_six_state as compiler


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


@njit
def backward(rows):
    n = len(rows)
    dp = np.zeros((1 << n, n), dtype=np.int64)
    for v in range(n):
        dp[1 << v, v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask >> v & 1):
                continue
            rest = mask ^ (1 << v)
            if rest:
                for u in range(n):
                    if rows[u] >> v & 1:
                        dp[mask, v] += dp[rest, u]
    return dp


@njit
def forward(rows):
    n = len(rows)
    dp = np.zeros((1 << n, n), dtype=np.int64)
    for v in range(n):
        dp[1 << v, v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if dp[mask, v]:
                for w in range(n):
                    if not (mask >> w & 1) and (rows[v] >> w & 1):
                        dp[mask | (1 << w), w] += dp[mask, v]
    return dp


def birth(rows, root):
    n = len(rows)
    return tuple(row | ((1 << n) if v == root else 0)
                 for v, row in enumerate(rows)) + (((1 << n) - 1) ^ (1 << root),)


def literal(rows, vertices):
    return sum(all(rows[a] >> b & 1 for a, b in zip(p, p[1:]))
               for p in permutations(vertices))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(output):
    source = RESULTS / 'multicore_metagraph_20261113' / 'next_target.json.gz'
    head = json.loads(gzip.decompress(source.read_bytes()))
    targets = [1691, 3311, 6551, 13031, 25991]
    records = []
    for hit in head['hits']:
        target = hit['H']
        if target not in targets:
            continue
        bound = next(b for b in head['bounds'] if b['H'] == target)
        need(bound['complete'], 'incomplete inherited target head')
        parent = compiler.a.unrank(hit['order'], hit['rank'])
        rows = birth(parent, 0)
        n = len(rows)
        full = (1 << n) - 1
        array = np.array(rows, dtype=np.int64)
        ends = backward(array)
        need(np.array_equal(ends, forward(array)), 'independent recurrence mismatch')
        h = ends.sum(axis=1)
        h[0] = 1
        need(h[full] == target, 'inherited compiler versus full subset count')
        need(all(value & 1 for value in h), 'odd Hamiltonian count control')
        root = n - 1
        start = backward(np.array([full ^ (1 << i) ^ row
                                   for i, row in enumerate(rows)], dtype=np.int64))
        F = sum(int(ends[mask, root]) * int(h[full ^ mask])
                for mask in range(1, full + 1))
        s = int(start[full, root])
        q = int(h[full ^ (1 << root)])
        need(F - target == q - s == 1, 'literal unit readiness')
        previous = (target + 71) // 2
        witnesses = [mask for mask, value in enumerate(h) if value == previous]
        need(not witnesses, 'unexpected induced predecessor count')
        # A positive hereditary control: the exact deletion parent is present.
        parent_h = int(h[(1 << len(parent)) - 1])
        need(parent_h == hit['state'][0], 'positive deletion-parent control')
        controls = []
        for vertices in [tuple(range(min(8, n))), tuple(range(n - min(7, n), n))]:
            mask = sum(1 << v for v in vertices)
            count = literal(rows, vertices)
            need(count == h[mask], 'literal permutation control')
            controls.append(dict(vertices=vertices, H=count))
        spectrum = Counter(map(int, h))
        below = max(value for value in spectrum if value < previous)
        above = min(value for value in spectrum if value > previous)
        by_size = []
        for size in range(n + 1):
            values = [int(h[mask]) for mask in range(full + 1) if mask.bit_count() == size]
            by_size.append(dict(size=size, subsets=len(values), distinct=len(set(values)),
                                minimum=min(values), maximum=max(values)))
        records.append(dict(parent_order=hit['order'], parent_rank=hit['rank'],
                            target=target, previous=previous, rows=rows, root=root,
                            state=[target, s, q, F], subsets=full + 1,
                            missing_previous=True, nearest_counts=[below, above],
                            spectrum=sorted(spectrum.items()), by_size=by_size,
                            literal_controls=controls,
                            subset_table_sha256=hashlib.sha256(h.astype('<i8').tobytes()).hexdigest()))
        print(target, hit['order'], hit['rank'], 'subsets', full + 1,
              'gap', below, previous, above, flush=True)
    counts = Counter(record['target'] for record in records)
    need([counts[x] for x in targets] == [3, 3, 4, 2, 2], 'inherited complete class counts')
    # A deliberately positive birth example prevents a blanket no-growth test.
    c3 = (2, 4, 1)
    child = birth(c3, 0)
    values = backward(np.array(child, dtype=np.int64)).sum(axis=1)
    need(values[7] == 3 and values[-1] == 5, 'positive induced birth control')
    pins = {str(source.relative_to(RESULTS.parent.parent)).replace('\\', '/'): digest(source)}
    for path in [Path(compiler.__file__), Path(compiler.a.__file__)] :
        pins[str(path.relative_to(RESULTS.parent.parent)).replace('\\', '/')] = digest(path)
    receipt = dict(status='FINITE-EXACT; all-order conclusion additionally uses the cited complete head and SCC reduction',
                   universe='All subsets, including empty, of all fourteen strong unit-ready targets at the first five successors on 405*2^k+71',
                   inherited_filters='Complete unique-start parent grammar and cited strong-order bound; no further subset filters',
                   source_pins=pins, total_subsets=sum(r['subsets'] for r in records),
                   independent_recurrences=True, positive_control=dict(rows=child, source_H=3, target_H=5),
                   records=records)
    (HERE / output).write_bytes(gzip.compress(json.dumps(receipt, sort_keys=True).encode('utf-8'), mtime=0))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default='receipt.json.gz')
    run(parser.parse_args().output)
