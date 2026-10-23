"""Canonicalize supplied named equitable cells of sizes1,3,5,7.

No target inventory is supplied. Exact matrices and all actual cell owners
are inputs. A canonicalizing chart is a gauge, not an intrinsic owner choice.
"""
from collections import Counter
from functools import lru_cache
from itertools import combinations, permutations
from pathlib import Path
import argparse
import gzip
import hashlib
import importlib.util
import json
import random

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('seven_kernel_checks', HERE / 'check.py')
kernel_checks = importlib.util.module_from_spec(spec)
spec.loader.exec_module(kernel_checks)
CHECKS = 0


def need(b, message):
    global CHECKS
    CHECKS += 1
    if not b:
        raise ValueError(message)


@lru_cache(None)
def templates(A):
    n = len(A)
    need(n in (1, 3, 5, 7), 'unsupported cell size')
    need(all(sum(row) == (n - 1) // 2 for row in A), 'cell is not regular')
    best = None
    charts = []
    for p in permutations(range(n)):
        key = tuple(A[p[i]][p[j]] for i, j in combinations(range(n), 2))
        if best is None or key < best:
            best, charts = key, [p]
        elif key == best:
            charts.append(p)
    return best, tuple(charts)


def canonicalize(A, cells):
    cells = tuple(tuple(D) for D in cells)
    need(sorted(v for D in cells for v in D) == list(range(len(A))), 'cells do not partition owners')
    options = []
    for D in cells:
        _, local = templates(tuple(tuple(A[u][v] for v in D) for u in D))
        options.append(tuple(tuple(D[i] for i in p) for p in local))
    support = {i: [] for i in range(len(cells))}
    for i, j in combinations(range(len(cells)), 2):
        B = [[A[u][v] for v in cells[j]] for u in cells[i]]
        need(len({sum(row) for row in B}) == 1, 'not row equitable')
        need(len({sum(row[k] for row in B) for k in range(len(cells[j]))}) == 1, 'not column equitable')
        if 0 < sum(B[0]) < len(cells[j]):
            support[i].append(j)
            support[j].append(i)
    visited = set()
    chosen = {}
    metadata = []
    for root in range(len(cells)):
        if root in visited:
            continue
        tree_order = [root]
        parent = {root: None}
        visited.add(root)
        for i in tree_order:
            for j in sorted(support[i]):
                if j not in visited:
                    visited.add(j)
                    parent[j] = i
                    tree_order.append(j)
        component = sorted(tree_order)
        candidates = []
        for root_chart in options[root]:
            maps = {root: root_chart}
            for j in tree_order[1:]:
                i = parent[j]
                costs = [(tuple(A[u][v] for u in maps[i] for v in q), q) for q in options[j]]
                best_cost = min(cost for cost, _ in costs)
                winners = [q for cost, q in costs if cost == best_cost]
                need(len(winners) == 1, 'one-sided kernel or hypothesis failure')
                maps[j] = winners[0]
            owners = tuple(v for i in component for v in maps[i])
            key = tuple(A[owners[i]][owners[j]] for i, j in combinations(range(len(owners)), 2))
            candidates.append((key, owners, maps))
        minimum = min(key for key, owners, maps in candidates)
        winners = [(owners, maps) for key, owners, maps in candidates if key == minimum]
        # A literal representative may break residual symmetry by input labels;
        # its MATRIX is canonical. Keep the full canonicalizer fibre separately.
        owners, maps = min(winners, key=lambda item: item[0])
        chosen.update(maps)
        metadata.append({'cells': component, 'root': root,
                         'root_charts_tried': len(options[root]),
                         'canonicalizer_count': len(winners),
                         'canonicalizers': [list(o) for o, _ in winners],
                         'local_group_orders': [len(options[i]) for i in component]})
    order = tuple(v for i in range(len(cells)) for v in chosen[i])
    key = (tuple(map(len, cells)), tuple(A[order[i]][order[j]] for i, j in combinations(range(len(A)), 2)))
    return key, order, metadata


def canonical_port(A, cells, path):
    need(sorted(path) == list(range(len(A))) and all(A[u][v] for u, v in zip(path, path[1:])), 'not actual Hamiltonian path')
    key, order, data = canonicalize(A, cells)
    offsets = []
    at = 0
    for D in cells:
        offsets.append(at)
        at += len(D)
    chosen = {}
    for component in data:
        coords = tuple(offsets[i] + j for i in component['cells'] for j in range(len(cells[i])))
        choices = []
        for owners in component['canonicalizers']:
            inverse = dict(zip(owners, coords))
            word = tuple(inverse[v] for v in path if v in inverse)
            choices.append((word, tuple(owners)))
        word, owners = min(choices)
        need(sum(candidate == word for candidate, _ in choices) == 1, 'full path did not kill canonicalizer symmetry')
        chosen.update(zip(coords, owners))
    order = tuple(chosen[i] for i in range(len(A)))
    inverse = {v: i for i, v in enumerate(order)}
    encoded = tuple(inverse[v] for v in path)
    need(tuple(A[order[i]][order[j]] for i, j in combinations(range(len(A)), 2)) == key[1], 'port gauge changed canonical source')
    return (key, encoded), order


def double(A):
    n = len(A)
    B = [[0] * (2 * n + 1) for _ in range(2 * n + 1)]
    for i in range(n):
        B[n][i] = B[n + 1 + i][n] = 1
        for j in range(n):
            B[i][j] = A[i][j]
            B[n + 1 + i][n + 1 + j] = A[j][i]
            B[i][n + 1 + j] = int(i == j or A[i][j])
            B[n + 1 + j][i] = 1 - B[i][n + 1 + j]
    return B


def owned_path(A):
    P = []
    for v in range(len(A)):
        i = next((i for i, u in enumerate(P) if A[v][u]), len(P))
        P.insert(i, v)
    return tuple(P)


def relabel(A, cells, p):
    B = [[0] * len(A) for _ in A]
    for i, j in combinations(range(len(A)), 2):
        B[p[i]][p[j]], B[p[j]][p[i]] = A[i][j], A[j][i]
    return B, [tuple(p[v] for v in D) for D in cells]


def path_complement(A, P):
    edges = set(zip(P, P[1:]))
    return [[0 if i == j else A[i][j] if (i, j) in edges or (j, i) in edges else 1 - A[i][j]
             for j in range(len(A))] for i in range(len(A))]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    rng = random.Random(731001)
    P7 = [[int((j - i) % 7 in (1, 2, 4)) for j in range(7)] for i in range(7)]
    A = double(double(P7))
    P = owned_path(A)
    histogram = Counter()
    port_cases = 0
    chart_trials = 0
    for u, v in ((u, v) for u in range(31) for v in range(31) if u != v):
        cells, _ = kernel_checks.refine(A, kernel_checks.marked_cells(31, (u, v)))
        key, order, data = canonicalize(A, cells)
        group = 1
        for record in data:
            group *= record['canonicalizer_count']
            chart_trials += record['root_charts_tried']
        histogram[max(map(len, cells)), group] += 1
        p = list(range(31))
        rng.shuffle(p)
        B, new_cells = relabel(A, cells, p)
        new_key, new_order, new_data = canonicalize(B, new_cells)
        need(new_key == key, 'canonical source changed under literal relabelling')
        need([x['canonicalizer_count'] for x in data] == [x['canonicalizer_count'] for x in new_data], 'canonicalizer group count changed')
        f = dict(zip(order, new_order))
        need(all(A[i][j] == B[f[i]][f[j]] for i in range(31) for j in range(31)), 'canonical code match did not construct isomorphism')
        if max(map(len, cells)) == 7:
            port, source_order = canonical_port(A, cells, P)
            new_port, target_order = canonical_port(B, new_cells, tuple(p[x] for x in P))
            need(port == new_port, 'actual port code changed under relabelling')
            f = dict(zip(source_order, target_order))
            need(tuple(f[x] for x in P) == tuple(p[x] for x in P), 'canonical port map lost actual path')
            K, L = path_complement(A, P), path_complement(B, tuple(p[x] for x in P))
            need(all(K[i][j] == L[f[i]][f[j]] for i in range(31) for j in range(31)), 'actual complement did not transport')
            port_cases += 1
    expected = {(1, 1): 504, (3, 1): 126, (3, 3): 168, (7, 1): 42, (7, 3): 84, (7, 21): 6}
    need(dict(histogram) == expected, 'order31 canonicalizer count differs from inherited automorphism census')
    inherited = json.loads(gzip.decompress((HERE.parents[1] / 'endpoint_overlap_ports_20261023' / 'receipt.json.gz').read_bytes()))
    need(inherited['matrix_outmasks'] == [sum(x << j for j, x in enumerate(row)) for row in A], 'actual port source matrix differs')
    native_cases = 0
    native_classes = Counter()
    for record in inherited['actual_ports']['records']:
        path = tuple(record['path'])
        cells, _ = kernel_checks.refine(A, kernel_checks.marked_cells(31, (path[0], path[-1])))
        port, source_order = canonical_port(A, cells, path)
        p = list(range(31))
        rng.shuffle(p)
        B, new_cells = relabel(A, cells, p)
        new_port, target_order = canonical_port(B, new_cells, tuple(p[x] for x in path))
        need(port == new_port, 'actual endpoint port changed under relabelling')
        f = dict(zip(source_order, target_order))
        need(tuple(f[x] for x in path) == tuple(p[x] for x in path), 'actual endpoint port map lost path')
        K, L = path_complement(A, path), path_complement(B, tuple(p[x] for x in path))
        need([sum(x << j for j, x in enumerate(row)) for row in K] == record['target_outmasks'], 'inherited actual complement differs')
        need(all(K[i][j] == L[f[i]][f[j]] for i in range(31) for j in range(31)), 'actual endpoint complement did not transport')
        native_cases += 1
        native_classes[max(map(len, cells))] += 1
    # Independent actual symmetry control: a mixed Aut3/Aut7 component is rigid.
    regular = kernel_checks.regular_sevens()
    by_size = {}
    for C in regular:
        G = kernel_checks.autos(C)
        by_size.setdefault(len(G), C)
        if len(by_size) == 3:
            break
    control_A = by_size[3]
    control_B = by_size[7]
    M = tuple(tuple(int(i == j) for j in range(7)) for i in range(7))
    T = kernel_checks.join(control_A, control_B, M)
    key, order, data = canonicalize(T, [tuple(range(7)), tuple(range(7, 14))])
    need(data[0]['canonicalizer_count'] == 1, 'coprime local group orders did not force rigidity')
    out = {'status': 'FINITE-EXACT source-only canonicalization and actual-port controls',
           'ordered_two_mark_sources': 930,
           'largest_cell_and_canonicalizer_counts': [[a, b, n] for (a, b), n in sorted(histogram.items())],
           'independent_random_relabellings': 930,
           'arbitrary_two_mark_path_complement_controls': port_cases,
           'actual_endpoint_port_controls': native_cases,
           'actual_endpoint_largest_cell_counts': {str(k): v for k, v in sorted(native_classes.items())},
           'inherited_distinct_literal_paths': inherited['actual_ports']['distinct_literal_paths'],
           'retained_mod8_receipt_sha256_lf': inherited['retained_mod8_receipt_sha256_lf'],
           'root_charts_tried_in_first_census': chart_trials,
           'coprime_internal_group_control': {'local_orders': [3, 7], 'canonicalizers': 1},
           'checks': CHECKS,
           'source_sha256_lf': hashlib.sha256(Path(__file__).read_bytes().replace(b'\r\n', b'\n')).hexdigest()}
    text = json.dumps(out, indent=2, sort_keys=True) + '\n'
    target = HERE / 'canonical.json'
    if args.write:
        target.write_text(text, encoding='utf8')
    else:
        need(json.loads(target.read_text(encoding='utf8')) == out, 'canonical receipt differs')
    print(json.dumps({'checks': CHECKS, 'census': out['largest_cell_and_canonicalizer_counts'],
                      'actual_endpoint_port_controls': native_cases, 'root_charts': chart_trials}, sort_keys=True))


if __name__ == '__main__':
    main()
