"""Exact ellipse geometry and bounded necessary-condition audit; no HN claim."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import networkx as nx

HERE = Path(__file__).resolve().parent
CHECKS = 0


def need(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(label)


# Basis 1,sqrt(2),sqrt(3),sqrt(6), with bit masks for the square roots.
def field(a=0, b=0, c=0, d=0):
    return tuple(map(F, (a, b, c, d)))


def add(a, b):
    return tuple(x+y for x, y in zip(a, b))


def neg(a):
    return tuple(-x for x in a)


def mul(a, b):
    out = [F(0)] * 4
    for i in range(4):
        for j in range(4):
            common = i & j
            out[i ^ j] += a[i]*b[j]*(2 if common & 1 else 1)*(3 if common & 2 else 1)
    return tuple(out)


def normdiff(p, q):
    v = [add(p[i], neg(q[i])) for i in range(2)]
    return add(mul(v[0], v[0]), mul(v[1], v[1]))


def gallai(g):
    for b in nx.biconnected_components(g):
        h = g.subgraph(b)
        n, e = len(h), h.number_of_edges()
        if e == n*(n-1)//2:
            continue
        if n % 2 and all(d == 2 for _, d in h.degree()):
            continue
        return False
    return True


def structural_candidate(g, special):
    light = set(g) - set(special)
    if not light or any(g.degree(v) != 4 for v in light):
        return False
    if any(g.degree(v) < 3 for v in special):
        return False
    return all(gallai(g.subgraph(c)) for c in nx.connected_components(g.subgraph(light)))


def run():
    one = field(1)
    points = [(field(), field(F(1, 4))),
              (field(d=F(2, 5)), field(F(1, 20))),
              (field(d=F(-2, 5)), field(F(1, 20))),
              (field(b=F(2, 3)), field(F(-1, 12))),
              (field(b=F(-2, 3)), field(F(-1, 12)))]
    need(len(set(points)) == 5, 'five distinct exact points')
    for x, y in points:
        need(add(mul(x, x), mul(field(16), mul(y, y))) == one, 'ellipse membership')
    edges = [(i, j) for i, j in combinations(range(5), 2)
             if normdiff(points[i], points[j]) == one]
    need(edges == [(0, i) for i in range(1, 5)], 'sharp four-neighbour star')
    # Independent scalar equations for the two y roots of the intersection.
    for y in (F(1, 20), F(-1, 12)):
        need(240*y*y+8*y-1 == 0, 'eliminated quadratic')
        need(1-16*y*y+(y-F(1, 4))**2 == 1, 'unit circle scalar check')

    counts = {'atlas_graphs': 0, 'geometric_necessary_graphs': 0,
              'special_choices': 0, 'failed_minimal_candidates': 0,
              'codegree_removed_candidates': 0}
    near = None
    for g in nx.graph_atlas_g():
        counts['atlas_graphs'] += 1
        if not g or not nx.is_connected(g) or max(dict(g.degree()).values()) > 4:
            continue
        if any(len(c) >= 4 for c in nx.find_cliques(g)):
            continue
        codegree = max((len(set(g[u]) & set(g[v])) for u, v in combinations(g, 2)), default=0)
        if codegree <= 2:
            counts['geometric_necessary_graphs'] += 1
        for k in range(min(2, len(g))+1):
            for special in combinations(g, k):
                counts['special_choices'] += 1
                if not structural_candidate(g, special):
                    continue
                if codegree <= 2:
                    counts['failed_minimal_candidates'] += 1
                else:
                    counts['codegree_removed_candidates'] += 1
                    if near is None:
                        near = {'vertices': len(g), 'edges': sorted(map(list, g.edges())),
                                'special': list(special), 'maximum_codegree': codegree}
    need(counts['atlas_graphs'] == 1253, 'full standard atlas through seven vertices')
    need(counts['failed_minimal_candidates'] == 0, 'no necessary critical candidate')
    need(near is not None and near['vertices'] == 7, 'retain codegree hostile mechanism')
    return {'status': 'FINITE-EXACT bounded audit; general proof in README',
            'checks': CHECKS, 'unit_edges': edges, 'counts': counts,
            'codegree_guard_near_miss': near,
            'checker_sha256_lf': sha256(Path(__file__).read_text(encoding='utf-8').encode()).hexdigest()}


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    result = run()
    path = HERE / 'receipt.json'
    if a.verify:
        need(json.loads(path.read_text(encoding='utf-8')) == json.loads(json.dumps(result)), 'saved receipt replay')
    else:
        path.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, indent=2))
