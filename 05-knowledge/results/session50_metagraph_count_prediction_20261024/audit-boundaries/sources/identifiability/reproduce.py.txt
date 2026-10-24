#!/usr/bin/env python3
"""Construct exact abstract C,R countermodels; no third-party dependencies."""
from collections import Counter
from fractions import Fraction
from math import lcm
from pathlib import Path
import argparse, json

ROOT = Path(__file__).resolve().parent

def require(ok, message):
    if not ok:
        raise RuntimeError(message)

def collision(degree):
    odds = list(range(1, 2 * degree + 2, 2))
    weights = []
    for a in odds:
        z = Fraction(1)
        for b in odds:
            if b != a:
                z *= Fraction(2 - b, a - b)
        weights.append(z)
    denominator = lcm(*(z.denominator for z in weights))
    relation = {2: denominator}
    relation.update({a: -int(z * denominator) for a, z in zip(odds, weights)})
    positive = sorted(a for a, c in relation.items() for _ in range(max(c, 0)))
    negative = sorted(a for a, c in relation.items() for _ in range(max(-c, 0)))
    # Make the common leaf mass odd, so the shared center is an SC node.
    if sum(positive) % 2 == 0:
        positive.append(1)
        negative.append(1)
        positive.sort()
        negative.sort()
    require(len(positive) == len(negative), 'zeroth moment')
    for q in range(degree + 1):
        require(sum(d**q for d in positive) == sum(d**q for d in negative), f'moment {q}')
    return denominator, relation, positive, negative

def cube_parameters(n):
    m = (n - 1) * (n - 2) // 2
    f = (n - 1) // 2
    return m, f, 1 << m, 1 << ((m + f) // 2)

def choose_order(leaves):
    base = sum(leaves[0])
    maximum_s = max(sum(d % 2 for d in row) for row in leaves)
    for n in range(3, 100):
        _, _, total, blue = cube_parameters(n)
        if blue > 2 * maximum_s and total > blue + 2 * base:
            return n
    raise RuntimeError('no order')

def star_graph(leaves, n):
    _, _, total, blue_mass = cube_parameters(n)
    center = len(leaves)
    center_degree = total - sum(leaves)
    s_leaves = sum(d % 2 for d in leaves)
    b_center = blue_mass - s_leaves
    blue = Counter()
    black_half = Counter()
    for i, d in enumerate(leaves):
        b = d % 2
        if b:
            blue[i, center] = 1
        h = (d - b) // 2
        if h:
            black_half[i, center] = h
    blue[center, center] = (b_center - s_leaves) // 2
    h_center = (center_degree - b_center) // 2
    black_half[center, center] = (h_center - sum(black_half.values())) // 2
    require(blue[center, center] > 0 and black_half[center, center] > 0, 'strict padding')
    d = list(leaves) + [center_degree]
    s = [bool(x % 2) for x in d]
    require(sum(blue.values()) * 2 == blue_mass, 'blue mass')
    require(sum(black_half.values()) * 4 == total - blue_mass, 'black mass')
    return {
        'n': n, 'degrees': d, 'self_converse': s,
        'blue_edges': [[u, v, w] for (u, v), w in sorted(blue.items())],
        'black_half_edges': [[u, v, w] for (u, v), w in sorted(black_half.items())],
    }

def lift_actual_cube(graph):
    n = graph['n']
    coords = [(i, j) for i in range(n) for j in range(i + 2, n)]
    loc = {p: i for i, p in enumerate(coords)}
    perm = [loc[n - 1 - j, n - 1 - i] for i, j in coords]
    total = 1 << len(coords)
    full = total - 1
    def reflect(x):
        return sum(((x >> j) & 1) << perm[j] for j in range(len(coords)))
    seen = set()
    blue_orbits, black_orbits = [], []
    for x in range(total):
        if x in seen:
            continue
        rx = reflect(x)
        orbit = {x, x ^ full, rx, rx ^ full}
        seen.update(orbit)
        if rx == x:
            require(len(orbit) == 2, 'blue orbit size')
            blue_orbits.append(x)
        else:
            require(len(orbit) == 4, 'black orbit size')
            black_orbits.append(x)
    merged = [-1] * total
    ordinary = [None] * total
    def assign(x, node, sheet):
        require(merged[x] == -1, 'assignment overlap')
        merged[x] = node
        ordinary[x] = [node, 0 if graph['self_converse'][node] else sheet]
    bpos = 0
    for u, v, w in graph['blue_edges']:
        for _ in range(w):
            x = blue_orbits[bpos]
            bpos += 1
            assign(x, u, 0)
            assign(x ^ full, v, 0)
    kpos = 0
    for u, v, w in graph['black_half_edges']:
        for _ in range(w):
            x = black_orbits[kpos]
            kpos += 1
            rx = reflect(x)
            assign(x, u, 0)
            assign(rx, u, 1)
            assign(x ^ full, v, 0)
            assign(rx ^ full, v, 1)
    require(bpos == len(blue_orbits) and kpos == len(black_orbits), 'orbit exhaustion')
    require(-1 not in merged, 'word exhaustion')
    return dict(graph, merged_by_word=merged, ordinary_by_word=ordinary)

def connect_degrees(degrees, loops):
    """Loopless greedy pairing, then degree-preserving component joins."""
    edges = Counter({(i, i): q for i, q in loops.items() if q})
    residual = [d - 2 * loops.get(i, 0) for i, d in enumerate(degrees)]
    require(all(d >= 0 for d in residual), 'loop budget')
    require(sum(residual) % 2 == 0, 'degree parity')
    require(not residual or max(residual) * 2 <= sum(residual), 'loopless degree bound')
    while any(residual):
        order = sorted(range(len(residual)), key=lambda i: (-residual[i], i))
        u, v = order[:2]
        require(residual[v] > 0, 'unpaired stub')
        edges[min(u, v), max(u, v)] += 1
        residual[u] -= 1
        residual[v] -= 1
    active = {i for i, d in enumerate(degrees) if d}
    def components(skip=None):
        unseen = set(active)
        out = []
        while unseen:
            stack = [unseen.pop()]
            comp = set(stack)
            while stack:
                u = stack.pop()
                for (a, b), w in edges.items():
                    if w - int(skip == (a, b)) <= 0 or a == b:
                        continue
                    v = b if a == u else a if b == u else None
                    if v in unseen:
                        unseen.remove(v);comp.add(v);stack.append(v)
            out.append(comp)
        return out
    while len(components()) > 1:
        comps = components()
        cyclic = None
        for e, w in sorted(edges.items()):
            if e[0] != e[1] and w and len(components(skip=e)) == len(comps):
                cyclic = e
                break
        require(cyclic is not None, 'insufficient cycle reserve for connected realization')
        u, v = cyclic
        own = next(c for c in comps if u in c)
        other = next(c for c in comps if c != own)
        target = next((e for e, w in sorted(edges.items()) if w and e[0] != e[1] and e[0] in other), None)
        require(target is not None, 'isolated loop-only component')
        a, b = target
        edges[cyclic] -= 1;edges[target] -= 1
        edges[min(u, a), max(u, a)] += 1
        edges[min(v, b), max(v, b)] += 1
    return [[u, v, w] for (u, v), w in sorted(edges.items()) if w]

def bounded_graph(pairs):
    degrees = [d for d, b in pairs]
    blue = [b for d, b in pairs]
    half = [(d-b)//2 for d, b in pairs]
    blue_loops = {}
    for i, b in enumerate(blue):
        if b >= 3 and len(blue_loops) < 2:
            blue_loops[i] = 1
    black_loops = {}
    for i, h in enumerate(half):
        if blue[i] == 0 and h >= 7 and len(black_loops) < 4:
            black_loops[i] = 3
    require(len(blue_loops) == 2 and len(black_loops) == 4, 'bounded loop allocation')
    return {'n': 6, 'degrees': degrees, 'self_converse': [bool(b) for b in blue],
            'blue_edges': connect_degrees(blue, blue_loops),
            'black_half_edges': connect_degrees(half, black_loops)}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    hierarchy = []
    for degree in range(1, 9):
        denominator, relation, a, b = collision(degree)
        n = choose_order([a, b])
        rows = [star_graph(x, n) for x in [a, b]]
        hierarchy.append({
            'moment_degree': degree, 'denominator': denominator,
            'signed_relation': {str(k): v for k, v in sorted(relation.items())},
            'n': n, 'V': len(a) + 1,
            'S': [sum(g['self_converse']) for g in rows],
            'A_abstract': [2 * len(g['degrees']) - sum(g['self_converse']) for g in rows],
            'degree_moments': [[sum(d**q for d in g['degrees']) for q in range(degree + 1)] for g in rows],
        })
    _, _, a, b = collision(2)
    n = choose_order([a, b])
    require(n == 6, 'explicit order')
    native = json.loads((ROOT / 'native_reference.json').read_text())[-1]
    native_pairs = sorted((z['degree'], z['blue_degree']) for z in native['nodes'])
    alternate_hist = [(1,1,2),(2,0,1),(10,0,1),(15,3,16),(29,7,1),(33,7,1),(38,0,2),(50,0,1),(58,0,1),(62,0,1),(66,0,7)]
    alternate_pairs = [(d,b) for d,b,count in alternate_hist for _ in range(count)]
    require(set(alternate_pairs) <= set(native_pairs), 'native joint support')
    bounded_models = [lift_actual_cube(bounded_graph(pairs)) for pairs in [native_pairs, alternate_pairs]]
    result = {
        'bounded_models': bounded_models,
        'schema': 'abstract-cube-identifiability-v1',
        'scope': 'Actual fixed C,R action; arbitrary odd-fibre partition, not tournament isomorphism fibres.',
        'hierarchy': hierarchy,
        'models': [lift_actual_cube(star_graph(x, n)) for x in [a, b]],
    }
    encoded = json.dumps(result, indent=2, sort_keys=True) + '\n'
    path = ROOT / 'certificate.json'
    if args.check:
        require(path.read_text() == encoded, 'producer replay differs')
        print(json.dumps({'producer_replay': 'PASS', 'models': 2, 'finite_moment_orders': list(range(1, 9))}, sort_keys=True))
    else:
        path.write_text(encoded)
        print(json.dumps({'written': path.name, 'n': n, 'models': 2}, sort_keys=True))

if __name__ == '__main__':
    main()
