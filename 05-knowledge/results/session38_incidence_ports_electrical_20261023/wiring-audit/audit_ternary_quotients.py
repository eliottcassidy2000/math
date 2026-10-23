"""Independent ternary-sum forest and exhaustive 255 binary E33 quotients.

Uses only the inherited, already independently certified canonical class list.
No producer imports, pruning, search transcript, or converse IDs are used.
"""
import argparse
import collections
import csv
import hashlib
import itertools
import json
import pathlib

HERE = pathlib.Path(__file__).resolve().parent


def need(test, message):
    if not test:
        raise ValueError(message)


def forest(left, right, codes, detailed=False):
    lv = {s: i for i, s in enumerate(sorted(set(left)))}
    rv = {s: len(lv) + i for i, s in enumerate(sorted(set(right)))}
    parent = list(range(len(lv) + len(rv)))
    degree = [0] * len(parent)
    edgecount = [0] * len(parent)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    cycles = []
    endpoints = []
    for code, a, b in zip(codes, left, right):
        a, b = lv[a], rv[b]
        degree[a] += 1
        degree[b] += 1
        endpoints.append((a, b))
        ra, rb = find(a), find(b)
        if ra == rb:
            cycles.append(code)
        else:
            parent[rb] = ra
    for a, b in endpoints:
        edgecount[find(a)] += 1
    components = collections.Counter(edgecount[i] for i in range(len(parent)) if find(i) == i)
    out = {
        'ordinary_signatures': len(lv),
        'converse_signatures': len(rv),
        'joint_signatures': len(set(zip(left, right))),
        'vertices': len(parent),
        'edges': len(codes),
        'components': sum(components.values()),
        'cycle_rank': len(cycles),
        'first_cycle_closing_canon': cycles[0] if cycles else None,
    }
    need(out['edges'] - out['vertices'] + out['components'] == len(cycles), 'Euler cycle rank')
    if detailed:
        out['component_edge_count_histogram'] = sorted(components.items())
        out['maximum_degree'] = max(degree)
        out['both_nonsingleton_canons'] = [code for code, (a, b) in zip(codes, endpoints) if degree[a] > 1 and degree[b] > 1]
        out['components'] = sum(components.values())
        out['leaf_peeling_rounds'] = []
        active = set(range(len(codes)))
        while active:
            take = {i for i in active if degree[endpoints[i][0]] == 1 or degree[endpoints[i][1]] == 1}
            if not take:
                break
            out['leaf_peeling_rounds'].append(len(take))
            for i in take:
                a, b = endpoints[i]
                degree[a] -= 1
                degree[b] -= 1
            active.difference_update(take)
        out['unpeeled_classes'] = len(active)
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    source = HERE / 'classes_input_n8.tsv'
    codes = sorted(int(r['canon']) for r in csv.DictReader(source.open(), delimiter='\t'))
    need(len(codes) == len(set(codes)) == 6880, 'complete certified class list')
    all_rows = []
    root_types = {}
    local_checks = 0
    for code in codes:
        A = [[0] * 8 for _ in range(8)]
        for k, (u, v) in enumerate(itertools.combinations(range(8), 2)):
            A[u][v] = (code >> (27 - k)) & 1
            A[v][u] = 1 - A[u][v]
        c = [0] * 8
        e = [[0] * 8 for _ in range(8)]
        triangles = 0
        # Direct cyclic triples, unlike the producer's two-step intersections.
        for u, v, w in itertools.combinations(range(8), 3):
            if A[u][v] == A[v][w] == A[w][u]:
                triangles += 1
                for x in (u, v, w):
                    c[x] += 1
                for x, y in ((u, v), (u, w), (v, w)):
                    e[x][y] += 1
                    e[y][x] += 1
        need(sum(c) == 3 * triangles, 'triple incidence mass')
        halves = []
        sums = []
        for reverse in (False, True):
            rowids = []
            weighted = []
            for u in range(8):
                hist = [0] * 9
                for v in range(8):
                    if A[v][u] if reverse else A[u][v]:
                        hist[3 * (c[v] % 3) + (e[u][v] % 3)] += 1
                hist = tuple(hist)
                rowids.append(root_types.setdefault(hist, len(root_types)))
                weighted.append(sum(((i // 3 + i % 3) % 3) * h for i, h in enumerate(hist)) % 3)
            halves.append(tuple(rowids))
            sums.append(weighted)
        all_rows.append(tuple(halves))
        for u in range(8):
            need((sums[0][u] + sums[1][u]) % 3 == c[u] % 3, 'aligned converse decoding')
            local_checks += 1
            for v in range(8):
                if A[u][v]:
                    color = (c[v] + e[u][v]) % 3
                    opposite = (c[u] + e[u][v]) % 3
                    need((color - c[v]) % 3 == e[u][v] % 3, 'edge decoding given aligned head')
                    need((color - opposite) % 3 == (c[v] - c[u]) % 3, 'original minus opposite gradient sign')
                    local_checks += 2
    types = list(root_types)

    def signatures(mapping):
        rows = []
        for hist in types:
            counts = [0] * (max(mapping) + 1)
            for symbol, count in enumerate(hist):
                counts[mapping[symbol]] += count
            rows.append(bytes(counts))
        left, right = [], []
        for forward, backward in all_rows:
            left.append(b''.join(sorted(rows[t] for t in forward)))
            right.append(b''.join(sorted(rows[t] for t in backward)))
        return left, right

    selected = [(i // 3 + i % 3) % 3 for i in range(9)]
    need(''.join(map(str, selected)) == '012120201', 'sum partition encoding')
    ternary = forest(*signatures(selected), codes, detailed=True)
    need(ternary['ordinary_signatures'] == 6864, 'ternary ordinary count')
    need(ternary['cycle_rank'] == 0 and ternary['unpeeled_classes'] == 0, 'ternary forest')
    need(ternary['component_edge_count_histogram'] == [(1, 6820), (2, 24), (3, 4)], 'component sizes')
    need(ternary['maximum_degree'] == 2, 'size-three components are paths')
    need(len(ternary['both_nonsingleton_canons']) == 4, 'four central path edges')
    print('PASS selected ternary sum forest; now directly checking all 255 binary quotients.', flush=True)
    binary = []
    # Canonical two-block partitions have residue 0 in block 0; at least one
    # of the remaining eight residues is in block 1. This is exactly 255.
    for bits in range(1, 256):
        mapping = [0] + [(bits >> i) & 1 for i in range(8)]
        info = forest(*signatures(mapping), codes)
        need(info['cycle_rank'] > 0, 'binary forest success contradicts lower bound')
        binary.append({'partition': ''.join(map(str, mapping)), **info})
    unary = forest(*signatures([0] * 9), codes)
    need(unary['cycle_rank'] > 0, 'unary failure')
    result = {
        'status': 'PASS',
        'scope': 'All 6880 certified order-eight classes; exact minimum of three symbols for the forest criterion among every set partition of the fixed nine E33 symbols. No claim that two symbols fail actual port-table uniqueness.',
        'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
        'method': 'Independent direct cyclic triples; direct converse orientation, no converse ID oracle; all 255 binary partitions evaluated without pruning or producer imports.',
        'selected_partition': '012120201',
        'selected': ternary,
        'aligned_vertex_local_identity_checks': local_checks,
        'distinct_nine_symbol_root_histograms': len(types),
        'binary_partitions_checked': len(binary),
        'binary_forest_successes': 0,
        'minimum_binary_cycle_rank': min(r['cycle_rank'] for r in binary),
        'maximum_binary_cycle_rank': max(r['cycle_rank'] for r in binary),
        'binary_results': binary,
        'unary': unary,
    }
    result = json.loads(json.dumps(result))
    dest = HERE / 'independent_ternary_quotient_receipt.json'
    if args.write:
        dest.write_text(json.dumps(result, indent=2) + '\n')
    else:
        need(json.loads(dest.read_text()) == result, 'frozen receipt mismatch')
    print(json.dumps({k: v for k, v in result.items() if k != 'binary_results'}, indent=2))


if __name__ == '__main__':
    main()
