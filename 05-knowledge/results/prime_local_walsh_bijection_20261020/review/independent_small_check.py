#!/usr/bin/env python3
"""Independent labelled controls through n=5 and inverse-minor checks.
No imports from the target checker; no n=6 enumeration.
"""
from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations
from math import factorial, gcd
from pathlib import Path
import json
import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
TARGET = ROOT


def require(ok, msg):
    if not ok:
        raise RuntimeError(msg)


def atlas(n):
    pairs = list(combinations(range(n), 2))
    ix = {e: k for k, e in enumerate(pairs)}
    perms = list(permutations(range(n)))
    universe = set(range(1 << len(pairs)))
    tours = []
    remaining = set(universe)
    while remaining:
        t = min(remaining)
        directed = [(v, u) if t >> k & 1 else (u, v)
                    for k, (u, v) in enumerate(pairs)]
        images = []
        for g in perms:
            backward_edges = [tuple(sorted((g[u], g[v])))
                              for u, v in directed if g[u] > g[v]]
            images.append(sum(1 << ix[e] for e in backward_edges))
        orbit = set(images)
        tours.append({'mask': t, 'aut': images.count(t), 'orbit': len(orbit)})
        require(len(orbit) * images.count(t) == factorial(n), 'tournament orbit-stabilizer')
        remaining.difference_update(orbit)
    graphs = []
    graph_count = 0
    remaining = set(universe)
    while remaining:
        x = min(remaining)
        edges = [e for k, e in enumerate(pairs) if x >> k & 1]
        signed_images = {}
        stabilizer = 0
        for g in perms:
            y = sum(1 << ix[tuple(sorted((g[u], g[v])))] for u, v in edges)
            sign = -1 if sum(g[u] > g[v] for u, v in edges) % 2 else 1
            signed_images.setdefault(y, set()).add(sign)
            stabilizer += y == x
        remaining.difference_update(signed_images)
        graph_count += 1
        if all(len(s) == 1 for s in signed_images.values()):
            signs = {y: next(iter(s)) for y, s in signed_images.items()}
            graphs.append({'mask': x, 'aut': stabilizer,
                           'orbit': len(signs), 'signs': signs})
    matrix = []
    for t in tours:
        row = []
        for g in graphs:
            value = 0
            for y, sign in g['signs'].items():
                character = 1
                for k in range(len(pairs)):
                    if (y >> k & 1) and (t['mask'] >> k & 1):
                        character = -character
                value += sign * character
            row.append(value)
        matrix.append(row)
    return tours, graphs, matrix, graph_count


def inverse(A, p):
    if not A.rows:
        return sp.zeros(0)
    return A.inv() if p is None else A.inv_mod(p)


def cofactor_replay(matrix, p):
    A = sp.Matrix(matrix)
    rows, cols = list(range(A.rows)), list(range(A.cols))
    Q = inverse(A, p)
    pairs, trace = [], []
    while rows:
        nz = (lambda x: bool(x)) if p is None else (lambda x: bool(x % p))
        j = next(j for j in range(A.cols) if nz(A[0, j]) and nz(Q[j, 0]))
        pairs.append((rows[0], cols[j]))
        trace.append([rows[0], cols[j], str(Q[j, 0])])
        keep = [k for k in range(A.cols) if k != j]
        retained = A.extract(range(1, A.rows), keep)
        denominator_inverse = 1 / Q[j, 0] if p is None else pow(int(Q[j, 0]), -1, p)
        updated = sp.Matrix(len(keep), len(keep), lambda k, l:
                            Q[keep[k], l + 1] - Q[keep[k], 0] * Q[j, l + 1] * denominator_inverse)
        if p is not None:
            updated = updated.applyfunc(lambda x: x % p)
        actual = inverse(retained, p)
        require(updated == actual, f'inverse-minor identity p={p}, step={len(pairs)}')
        A, Q = retained, actual
        rows, cols = rows[1:], [cols[k] for k in keep]
    return pairs, trace


def valuation(a, p):
    v = 0
    while a % p == 0:
        a //= p
        v += 1
    return v


def main():
    source_hashes = {}
    for line in (TARGET / 'CHECKPOINT_SHA256.txt').read_text().splitlines():
        expected, name = line.split()
        actual = sha256((TARGET / name).read_bytes()).hexdigest()
        require(expected == actual, f'frozen checksum {name}')
        source_hashes[name] = actual
    supplied = json.loads((TARGET / 'walsh_support.json').read_text())
    controls = []
    for n in range(6):
        T, G, matrix, graph_count = atlas(n)
        stored = supplied[n]
        require(matrix == stored['matrix'], f'n={n} independent Walsh matrix')
        for own, recorded in ((T, stored['tournaments']), (G, stored['even_graphs'])):
            require([{k: x[k] for k in ('mask', 'aut', 'orbit')} for x in own] ==
                    [{k: x[k] for k in ('mask', 'aut', 'orbit')} for x in recorded],
                    f'n={n} independent class catalogue')
        A = sp.Matrix(matrix)
        Da = sp.diag(*[t['aut'] for t in T])
        Db = sp.diag(*[g['aut'] for g in G])
        require(A * Db * A.T == 2**(n*(n-1)//2) * Da, 'weighted form identity')
        require(int(A.det()) == stored['det'], 'exact small determinant')
        checked_steps = 0
        for p in (None, 3, 5, 7):
            matching, trace = cofactor_replay(matrix, p)
            encoded = [[T[i]['mask'], G[j]['mask'], matrix[i][j]] for i, j in matching]
            expected = stored['rational_matching'] if p is None else stored['prime_data'][str(p)]['matching']
            expected_trace = stored['rational_cofactor_trace'] if p is None else stored['prime_data'][str(p)]['cofactor_trace']
            require(encoded == expected and trace == expected_trace, f'n={n} cofactor construction p={p}')
            if p is not None:
                require(all(valuation(T[i]['aut'], p) == valuation(G[j]['aut'], p)
                            for i, j in matching), 'matching valuations')
            checked_steps += len(matching)
        odd = lambda x: x // (x & -x)
        joint_equal = Counter(odd(t['aut']) for t in T) == Counter(odd(g['aut']) for g in G)
        require(joint_equal == (n < 5), 'first joint-prime failure')
        controls.append({'n': n, 'tournament_and_even_graph_classes': len(T),
                         'all_graph_classes': graph_count, 'determinant': int(A.det()),
                         'cofactor_steps_checked_by_fresh_inverse': checked_steps,
                         'joint_odd_parts_equal': joint_equal})
    n5 = supplied[5]
    row_indices = [next(i for i, t in enumerate(n5['tournaments']) if t['mask'] == mask)
                   for mask in (5, 76)]
    neighborhoods = [[j for j, x in enumerate(n5['matrix'][i]) if gcd(x, 15) == 1]
                     for i in row_indices]
    require(neighborhoods == [[0], [0]], 'two-row Hall obstruction')
    require(n5['even_graphs'][0]['mask'] == 0, 'Hall neighbor is empty graph')
    receipt = {'status': 'PASS', 'source_hashes': source_hashes,
               'independent_atlas_controls': controls,
               'hall_obstruction': {'tournament_masks': [5, 76],
                                    'automorphism_orders': [3, 5],
                                    'sole_coprime15_neighbor_mask': 0,
                                    'neighbor_automorphism_order': 120,
                                    'rows': [n5['matrix'][i] for i in row_indices]},
               'scope': 'Complete independent labelled atlases only through n=5. No n=6 enumeration. No imported target code.'}
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
