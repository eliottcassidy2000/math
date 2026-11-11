#!/usr/bin/env python3
"""Exact radical geometry and finite conditional-colouring response controls."""
import sympy as s
import itertools, json

def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)
from pathlib import Path
R = s.Rational
core = [s.Matrix(v) for v in [(0, 0), (1, 0), (1, 1), (0, 1)]]
ports = [s.Matrix([R(1, 3), R(2, 5)]), s.Matrix([R(7, 5), R(4, 3)])]
pts = core + ports
expected = {(0, 1), (1, 2), (2, 3), (0, 3)}
for j, p in enumerate(ports):
    for i, w in enumerate(core):
        v = w - p
        d = v.dot(v)
        require(0 < d < 4, 'Check failed: 0 < d < 4')
        x = (p + w) / 2 + s.sqrt(1 / d - R(1, 4)) * s.Matrix([-v[1], v[0]])
        k = len(pts)
        pts.append(x)
        expected |= {tuple(sorted((k, i))), tuple(sorted((k, 4 + j)))}
actual = set()
distinct = True
nonedge_bounds = []
for i, j in itertools.combinations(range(len(pts)), 2):
    d = s.expand((pts[i] - pts[j]).dot(pts[i] - pts[j]))
    d = s.radsimp(s.simplify(d))
    require(d != 0 and d.is_zero is False, (i, j, d))
    if s.simplify(d - 1) == 0:
        actual.add((i, j))
    else:
        require((d - 1).is_positive is True or (d - 1).is_negative is True, (i, j, d))
require(actual == expected, (actual - expected, expected - actual))

def response(n, edges, q, a, b):
    total = 0
    proper = 0
    for cs in itertools.product(range(q), repeat=n):
        if any((cs[i] == cs[j] for i, j in edges)):
            continue
        proper += 1
        v = 1
        for c in cs:
            v *= (q - 2 + (a == c)) * (q - 2 + (b == c))
        total += v
    return (total, proper)
edges = [(0, 1), (1, 2), (2, 3), (3, 0)]
rows = {}
for q in range(2, 7):
    eq, pc = response(4, edges, q, 0, 0)
    ne, _ = response(4, edges, q, 0, 1)
    rows[str(q)] = {'equal': eq, 'unequal': ne, 'core_colourings': pc}
    require(q == 2 or (eq > 0 and ne > 0), 'Check failed: q == 2 or (eq > 0 and ne > 0)')
require(response(3, [(0, 1), (1, 2), (0, 2)], 3, 0, 0)[0] == 24, 'Check failed: response(3, [(0, 1), (1, 2), (0, 2)], 3, 0, 0)[0] == 24')
require(response(3, [(0, 1), (1, 2), (0, 2)], 3, 0, 1)[0] == 24, 'Check failed: response(3, [(0, 1), (1, 2), (0, 2)], 3, 0, 1)[0] == 24')
require(response(4, list(itertools.combinations(range(4), 2)), 3, 0, 1)[0] == 0, 'Check failed: response(4, list(itertools.combinations(range(4), 2)), 3, 0, 1)[0] == 0')
require(response(1, [], 6, 0, 0)[0] == 105, 'Check failed: response(1, [], 6, 0, 0)[0] == 105')
require(response(1, [], 6, 0, 1)[0] == 104, 'Check failed: response(1, [], 6, 0, 1)[0] == 104')
triangle_equal, triangle_proper = response(3, [(0, 1), (1, 2), (0, 2)], 3, 0, 0)
triangle_unequal, _ = response(3, [(0, 1), (1, 2), (0, 2)], 3, 0, 1)
require(rows['2'] == {'equal': 0, 'unequal': 0, 'core_colourings': 2}, 'q=2 support boundary')
require(all(v['equal'] > v['unequal'] for k,v in rows.items() if int(k) >= 3), 'Square-core strict Gram bias')
out = {'status': 'FINITE-EXACT', 'vertices': len(pts), 'unit_edges': len(actual), 'injective': True, 'all_pairs_exact': 91, 'core': 'unit square C4', 'port_coordinates': [[str(t) for t in p] for p in ports], 'response': rows, 'single_core_vertex_q6': [105, 104], 'balanced_triangle_gram_control': {'q': 3, 'equal': triangle_equal, 'unequal': triangle_unequal, 'core_colourings': triangle_proper, 'all_class_sizes': [1, 1, 1]}}
Path(__file__).with_name('receipt.json').write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps(out, indent=2))
