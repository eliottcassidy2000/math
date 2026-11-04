"""Exact endpoint-sidecar isomorphism/automorphism census on inherited banks."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
SOURCE = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def rows_from_code(n, code):
    rows = [0] * n
    for i, (u, v) in enumerate(itertools.combinations(range(n), 2)):
        a, b = (u, v) if code & (1 << i) else (v, u)
        rows[a] |= 1 << b
    return tuple(rows)


def ends(rows):
    n = len(rows)
    dp = [[0] * n for _ in range(1 << n)]
    for mask in range(1, 1 << n):
        for v in range(n):
            if not mask & (1 << v):
                continue
            rem = mask ^ (1 << v)
            dp[mask][v] = sum(dp[rem][u] for u in range(n) if rows[u] & (1 << v)) if rem else 1
    return dp


def sidecar(rows):
    n = len(rows)
    full = (1 << n) - 1
    e = ends(rows)
    s = ends(tuple(full ^ (1 << u) ^ row for u, row in enumerate(rows)))
    c = [[0] * n for _ in range(n)]
    for mask in range(1, full):
        left = [(u, a) for u, a in enumerate(e[mask]) if a]
        right = [(v, b) for v, b in enumerate(s[full ^ mask]) if b]
        for u, a in left:
            for v, b in right:
                c[u][v] += a * b
    return tuple(s[full]), tuple(e[full]), tuple(tuple(row) for row in c)


def graph_code(rows, order):
    return sum(((rows[u] >> v) & 1) << bit
               for bit, (u, v) in enumerate(itertools.combinations(order, 2)))


def canonical(rows, data):
    s, e, c = data
    n = len(rows)
    def compress(keys):
        index = {key: i for i, key in enumerate(sorted(set(keys)))}
        return tuple(index[key] for key in keys)
    colors = compress(tuple(zip(s, e)))
    while True:
        keys = [(colors[u], tuple(sorted((colors[v], c[u][v], c[v][u])
                                        for v in range(n) if v != u))) for u in range(n)]
        nxt = compress(keys)
        if nxt == colors:
            break
        colors = nxt
    cells = [tuple(u for u in range(n) if colors[u] == col) for col in sorted(set(colors))]
    best = None
    codes = set()
    best_count = leaves = 0
    for pieces in itertools.product(*(itertools.permutations(cell) for cell in cells)):
        order = sum(pieces, ())
        key = tuple(x for u in order for x in (s[u], e[u])) + tuple(c[u][v] for u in order for v in order)
        leaves += 1
        if best is None or key < best:
            best, codes, best_count = key, {graph_code(rows, order)}, 1
        elif key == best:
            codes.add(graph_code(rows, order))
            best_count += 1
    return best, codes, best_count, leaves


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--max-n', type=int, default=8)
    args = parser.parse_args()
    if not 1 <= args.max_n <= 9:
        raise RuntimeError('bank supports orders1..9')
    rel = '05-knowledge/results/hn_metagraph_information_20261026/enumeration/'
    bank = ROOT / (rel + ('representatives9.json.gz' if args.max_n == 9 else 'representatives.json'))
    raw = bank.read_bytes()
    data = json.loads(gzip.decompress(raw) if bank.suffix == '.gz' else raw)
    expected = (1, 1, 1, 2, 4, 12, 56, 456, 6880, 191536)
    receipt = {'source_sha256': SOURCE, 'bank_path': bank.relative_to(ROOT).as_posix(),
               'bank_sha256': hashlib.sha256(raw).hexdigest(), 'orders': [], 'collisions': []}
    for n in range(1, args.max_n + 1):
        members = data['members'][str(n)]
        if len(members) != expected[n]:
            raise RuntimeError(('bank count', n, len(members)))
        seen = {}
        total_leaves = max_leaves = automorphism_failures = 0
        for index, member in enumerate(members):
            rows = rows_from_code(n, member['code'])
            side = sidecar(rows)
            key, codes, aut, leaves = canonical(rows, side)
            total_leaves += leaves
            max_leaves = max(max_leaves, leaves)
            if len(codes) > 1:
                receipt['collisions'].append({'n': n, 'kind': 'labelled automorphism', 'codes': sorted(codes)})
            if aut != member['aut']:
                automorphism_failures += 1
            if key in seen:
                receipt['collisions'].append({'n': n, 'kind': 'class collision',
                                              'source_codes': [seen[key], member['code']]})
            else:
                seen[key] = member['code']
        row = {'n': n, 'classes': len(members), 'sidecar_classes': len(seen),
               'canonical_leaves': total_leaves, 'max_leaves': max_leaves,
               'automorphism_group_disagreements': automorphism_failures}
        receipt['orders'].append(row)
        print(json.dumps(row, sort_keys=True), file=sys.stderr, flush=True)
    if hashlib.sha256(Path(__file__).read_bytes()).hexdigest() != SOURCE:
        raise RuntimeError('source changed during run')
    receipt['status'] = 'COLLISION' if receipt['collisions'] else 'PASS_INJECTIVE'
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
