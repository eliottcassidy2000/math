"""Finite feature-menu minimality, sparse wedge repair, and a DRT obstruction.

Default compares a frozen receipt. --write regenerates intentionally.
All acceptance checks remain active with python -O. Standard library only.
"""
from pathlib import Path
from itertools import combinations
from collections import Counter, defaultdict
import argparse, csv, hashlib, json, math

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
CHECKS = 0


def require(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)


def digest(value):
    return hashlib.sha256(json.dumps(value, separators=(',', ':'), sort_keys=True).encode()).hexdigest()


def matrix_rows(a):
    return [sum(v << i for i, v in enumerate(row)) for row in a]


def matrices(n, inherited, hashes):
    if n == 8:
        relative = '04-computation/experiments/tiling_line_metagraph_census_20261020_classes_n8.tsv'
    else:
        relative = f'04-computation/experiments/electric_metagraph_20261021/electric_n{n}.json'
    path = ROOT / relative; raw = path.read_bytes(); h = hashlib.sha256(raw).hexdigest()
    inherited_sources = {k.replace('\\', '/'): v for k, v in inherited['sources'].items()}
    require(h == inherited_sources[relative], ('changed inherited universe', n))
    hashes[relative] = h
    if n == 8:
        for row in csv.DictReader(raw.decode().splitlines(), delimiter='\t'):
            m = int(row['canon']); a = [[0] * n for _ in range(n)]
            for b, (i, j) in enumerate(combinations(range(n), 2)):
                a[i][j] = (m >> (27 - b)) & 1; a[j][i] = 1 - a[i][j]
            yield m, a
    else:
        source = json.loads(raw)
        free = [(a - 1, b - 1) for b in range(1, n - 1) for a in range(b + 2, n + 1)]
        for node in source['nodes']:
            for m in [node['rep_mask']] + ([] if node['SC'] else [node['conv_mask']]):
                a = [[0] * n for _ in range(n)]
                for i in range(1, n):
                    a[i][i - 1] = 1
                for bit, (i, j) in enumerate(free):
                    u, v = (j, i) if m >> bit & 1 else (i, j)
                    a[u][v] = 1
                yield m, a


def incidence(a, independent=False):
    n = len(a); degrees = list(map(sum, a)); cv = [0] * n
    edge = Counter(); triangles = []
    for t in combinations(range(n), 3):
        i, j, k = t
        if a[i][j] == a[j][k] == a[k][i]:
            triangles.append(t)
            for v in t:
                cv[v] += 1
            for e in combinations(t, 2):
                edge[e] += 1
    # O(n^3) rooted strong-four formula: central-edge and outer-vertex cases.
    v4 = [0] * n
    for (i, j), count in edge.items():
        v4[i] += math.comb(count, 2); v4[j] += math.comb(count, 2)
    for t in triangles:
        for v in t:
            opposite = tuple(w for w in t if w != v)
            v4[v] += edge[opposite] - 1
    if independent:
        direct = [0] * n
        for t in combinations(range(n), 4):
            if sorted(sum(a[u][v] for v in t) for u in t) == [1, 1, 2, 2]:
                for v in t:
                    direct[v] += 1
        require(v4 == direct, ('rooted overlap formula', n, matrix_rows(a)))
    moments = tuple(sum(math.comb(c, k) for c in edge.values()) for k in (2, 3, 4))
    require(sum(v4) == 4 * moments[0], 'rooted/global strong-four consistency')
    return list(zip(degrees, cv, v4)), moments, edge


def records(a, tags, feature_bits, rounds):
    rec = [tuple(t[k] for k in range(3) if feature_bits >> k & 1) for t in tags]
    for _ in range(rounds):
        rec = [(rec[v], tuple(sorted(rec[w] for w in range(len(a)) if a[v][w]))) for v in range(len(a))]
    return rec


def signature(a, tags, feature_bits, rounds):
    return tuple(sorted(records(a, tags, feature_bits, rounds)))


def simple_records(a, tags):
    t = [x[2] for x in tags]
    return [(t[v], tuple(sorted(t[w] for w in range(len(a)) if a[v][w]))) for v in range(len(a))]


def wedge(a, rec, r, b, c):
    return sum(sum(a[v][w] for w in range(len(a)) if rec[w] == b)
               * sum(a[v][w] for w in range(len(a)) if rec[w] == c)
               for v in range(len(a)) if rec[v] == r)


WEDGE_TYPES = [
    ((20, (16, 18, 20)), (16, (16, 16, 20)), (20, (14, 16, 20))),
    ((20, (14, 16, 16, 20)), (16, (14, 16, 18, 20)), (20, (14, 16, 16, 20))),
]

# The two ambiguous first-round record multisets, in ordinary unlabelled gauge.
GUARD_KEYS = {
    tuple(sorted([(14, (16, 16, 16, 20, 20)), (16, (16, 16, 20)),
                  (16, (16, 20, 20)), (16, (20, 20, 20)), (18, (14, 16, 16, 16, 20)),
                  (20, (14, 16, 20)), (20, (16, 18, 20)), (20, (16, 18, 20))])),
    tuple(sorted([(14, (18, 20)), (16, (14, 16, 16, 18)), (16, (14, 16, 18, 20)),
                  (16, (14, 18, 20, 20)), (18, (20, 20)),
                  (20, (14, 16, 16, 20)), (20, (14, 16, 16, 20)), (20, (16, 16, 18, 20))])),
}


def adaptive_wedge_signature(a):
    """Invariant for actual tournaments; the one-bit repair is scoped to n=8.

    Returns (first-round S4-only record multiset, bit). The bit is measured
    only on the two authenticated ambiguous multisets; it is zero elsewhere.
    """
    tags, _, _ = incidence(a)
    rec = simple_records(a, tags); base = tuple(sorted(rec))
    bit = sum(wedge(a, rec, *selector) for selector in WEDGE_TYPES) if len(a) == 8 and base in GUARD_KEYS else 0
    require(bit in (0, 1), 'guarded adaptive measurement is not a bit')
    return base, bit


def paley(q):
    squares = {i * i % q for i in range(1, q)}
    return [[int((j - i) % q in squares) for j in range(q)] for i in range(q)]


def double(a):
    n = len(a); b = [[0] * (2 * n + 1) for _ in range(2 * n + 1)]
    for i in range(n):
        b[n][i] = 1; b[n + 1 + i][n] = 1
        for j in range(n):
            b[i][j] = a[i][j]; b[n + 1 + i][n + 1 + j] = a[j][i]
            b[i][n + 1 + j] = int(i == j or a[i][j]); b[n + 1 + j][i] = 1 - b[i][n + 1 + j]
    return b


def is_drt(a):
    n = len(a)
    return (n % 4 == 3 and all(sum(row) == (n - 1) // 2 for row in a)
            and all(sum(a[i][k] * a[j][k] for k in range(n)) == (n - 3) // 4
                    for i, j in combinations(range(n), 2)))


def drt_out_count(a):
    n = len(a)
    return sum(is_drt([[a[i][j] for j in range(n) if a[v][j]]
                       for i in range(n) if a[v][i]]) for v in range(n))


def rank_three_profiles(a):
    """One paired-neighbor refinement of diagonal/arc/reverse-arc colors."""
    n = len(a)
    color = [[0 if i == j else 1 if a[i][j] else 2 for j in range(n)] for i in range(n)]
    by = defaultdict(set)
    for i in range(n):
        for j in range(n):
            profile = tuple(sorted(Counter((color[i][k], color[k][j]) for k in range(n)).items()))
            by[color[i][j]].add(profile)
    require(all(len(values) == 1 for values in by.values()), 'rank-three relation splits')
    return [(c, next(iter(by[c]))) for c in sorted(by)]


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--write', action='store_true'); args = parser.parse_args()
    inherited_path = ROOT / '05-knowledge/results/metagraph_bucket_obstruction_20261023/receipt.json'
    inherited_raw = inherited_path.read_bytes(); inherited = json.loads(inherited_raw)
    hashes = {str(inherited_path.relative_to(ROOT)): hashlib.sha256(inherited_raw).hexdigest()}
    out = {'scope': 'FINITE-EXACT complete inherited ordinary classes n=4..8; 192 explicit feature menus per order',
           'feature_bits': ['degree', 'cyclic_triangles_at_vertex', 'strong_fours_at_vertex'],
           'moment_bits': ['M2', 'M3', 'M4'], 'orders': {}, 'input_hashes': hashes}
    n8data = None
    for n in range(4, 9):
        mats = list(matrices(n, inherited, hashes)); expected = {4: 4, 5: 12, 6: 56, 7: 456, 8: 6880}[n]
        require(len(mats) == expected and len({m for m, a in mats}) == expected, ('inherited class list', n))
        data = [incidence(a, independent=True) for m, a in mats]; rows = []
        for f in range(8):
            for r in range(3):
                base = [signature(a, d[0], f, r) for (m, a), d in zip(mats, data)]
                for moment_bits in range(8):
                    by = defaultdict(list)
                    for (m, a), b, d in zip(mats, base, data):
                        by[b, tuple(v for k, v in enumerate(d[1]) if moment_bits >> k & 1)].append(m)
                    partition = sorted(sorted(v) for v in by.values()); bad = [v for v in partition if len(v) > 1]
                    rows.append({'features': f, 'rounds': r, 'moments': moment_bits, 'count': len(by),
                                 'collision_groups': len(bad), 'first_collision': bad[0] if bad else [],
                                 'partition_sha256': digest(partition)})
        wins = [x for x in rows if x['count'] == expected]
        minimal = [x for x in wins if not any(y is not x and y['features'] & x['features'] == y['features']
                   and y['rounds'] <= x['rounds'] and y['moments'] & x['moments'] == y['moments'] for y in wins)]
        out['orders'][str(n)] = {'classes': expected, 'menus': rows,
                                'minimal': [[x['features'], x['rounds'], x['moments']] for x in minimal]}
        require(len(rows) == 192, 'menu coverage')
        print(n, out['orders'][str(n)]['minimal'], flush=True)
        if n == 8:
            n8data = mats, data
    require(out['orders']['8']['minimal'] == [[2, 2, 1], [4, 2, 0]], 'minimal menu frontier changed')
    mats, data = n8data; groups = defaultdict(list)
    for (m, a), d in zip(mats, data):
        groups[tuple(sorted(simple_records(a, d[0])))].append((m, a, d))
    collision = [v for v in groups.values() if len(v) > 1]
    require(len(collision) == 2 and all(len(v) == 2 for v in collision), 'adaptive universe changed')
    require({k for k, v in groups.items() if len(v) > 1} == GUARD_KEYS, 'adaptive guard keys changed')
    adaptive = set(); witness = []
    for record_key, members in groups.items():
        if len(members) == 1:
            adaptive.add((record_key, 0)); continue
        pair = []
        for m, a, d in members:
            rec = simple_records(a, d[0]); types = sorted(set(rec)); ids = [types.index(x) for x in rec]
            edge_matrix = [[sum(a[i][j] for i in range(8) for j in range(8) if ids[i] == u and ids[j] == v)
                            for v in range(len(types))] for u in range(len(types))]
            bit = sum(wedge(a, rec, *selector) for selector in WEDGE_TYPES)
            require(bit in (0, 1), 'repair is not one bit on guarded collision family')
            adaptive.add((record_key, bit))
            pair.append({'canon': m, 'rows': matrix_rows(a), 'bit': bit, 'type_edge_matrix': edge_matrix})
        require(pair[0]['type_edge_matrix'] == pair[1]['type_edge_matrix'], 'aggregate edge matrix separates already')
        require({p['bit'] for p in pair} == {0, 1}, 'wedge bit failed')
        witness.append(pair)
    require(len(adaptive) == 6880, 'adaptive signature not complete at eight')
    require({adaptive_wedge_signature(a) for m, a in mats} == adaptive, 'public adaptive API differs')
    out['adaptive_n8'] = {'one_round_signatures': len(groups), 'repaired_signatures': len(adaptive),
                          'measure_only_on_classes': 4, 'selectors': WEDGE_TYPES, 'witness_pairs': witness}
    hostile = []; pair_profiles = []
    for name, a in [('Paley31', paley(31)), ('double_double_Paley7', double(double(paley(7))))]:
        require(all(a[i][j] + a[j][i] == 1 for i, j in combinations(range(31), 2)), 'hostile not tournament')
        require(is_drt(a), 'hostile not doubly regular')
        tags, moments, edge = incidence(a, independent=True)
        require(set(tags) == {(15, 120, 1680)}, 'DRT rooted tags differ')
        require(len(edge) == math.comb(31, 2) and set(edge.values()) == {8}, 'DRT edge incidence differs')
        pair_profiles.append(rank_three_profiles(a))
        if name == 'Paley31':
            require(all(a[(-i) % 31][(-j) % 31] == a[j][i] for i in range(31) for j in range(31)),
                    'Paley negation is not converse isomorphism')
        hostile.append({'name': name, 'rows': matrix_rows(a), 'tag': list(tags[0]), 'edge_incidence': 8,
                        'M2_M3_M4': moments, 'doubly_regular_outneighborhoods': drt_out_count(a)})
    require([h['doubly_regular_outneighborhoods'] for h in hostile] == [0, 9], 'nonisomorphism witness changed')
    require(pair_profiles[0] == pair_profiles[1], 'rank-three pair profiles differ')
    out['all_round_all_edge_moment_hostile'] = hostile
    out['rank_three_pair_profiles'] = pair_profiles[0]
    out['checks'] = CHECKS
    out['source_sha256_lf'] = hashlib.sha256(Path(__file__).read_bytes().replace(b'\r\n', b'\n')).hexdigest()
    if args.write:
        (HERE / 'receipt.json').write_text(json.dumps(out, indent=2, sort_keys=True) + '\n')
    else:
        require(digest(json.loads((HERE / 'receipt.json').read_text())) == digest(out), 'frozen receipt mismatch')
    print(json.dumps({'checks': CHECKS, 'receipt_sha256_semantic': digest(out), 'adaptive_classes': len(adaptive),
                      'DRT_out_counts': [h['doubly_regular_outneighborhoods'] for h in hostile]}))


if __name__ == '__main__':
    main()
