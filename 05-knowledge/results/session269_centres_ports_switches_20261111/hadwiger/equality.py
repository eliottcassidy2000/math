#!/usr/bin/env python3
import sympy as s, itertools, json
from pathlib import Path
r = s.sqrt

def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)
R = s.Rational
base = [s.Matrix([r(3) / 2, R(1, 2)]), s.Matrix([r(3) / 2, -R(1, 2)]), s.Matrix([r(3), 0])]
rot = s.Matrix([[R(107, 108), -r(215) / 108], [r(215) / 108, R(107, 108)]])
p = [s.zeros(2, 1)] + base + [rot * v for v in base]
expected = {(0, 1), (0, 2), (1, 2), (1, 3), (2, 3), (0, 4), (0, 5), (4, 5), (4, 6), (5, 6)}
edges = set()
for i, j in itertools.combinations(range(7), 2):
    d = s.simplify((p[i] - p[j]).dot(p[i] - p[j]))
    require(d.is_zero is False, 'Check failed: d.is_zero is False')
    if s.simplify(d - 1) == 0:
        edges.add((i, j))
    else:
        require((d - 1).is_positive is True or (d - 1).is_negative is True, 'Check failed: (d - 1).is_positive is True or (d - 1).is_negative is True')
require(edges == expected, 'Check failed: edges == expected')
require(s.simplify((p[3] - p[6]).dot(p[3] - p[6])) == R(1, 18), 'Check failed: s.simplify((p[3] - p[6]).dot(p[3] - p[6])) == R(1, 18)')
counts = {}
for q in [3, 4]:
    states = [0, 0]
    for c in itertools.product(range(q), repeat=7):
        if c[3] != 0 or c[6] not in (0, 1):
            continue
        if any((c[i] == c[j] for i, j in edges)):
            continue
        states[c[6]] += 1
    counts[str(q)] = states
require(counts['3'] == [4, 0], "Check failed: counts['3'] == [4, 0]")
require(counts['4'] == [48, 32], 'Higher-colour hostile control')
out = {'status': 'FINITE-EXACT', 'vertices': 7, 'edges': 10, 'all_pairs_exact': 21, 'squared_port_distance': '1/18', 'port_response_fixed_equal_unequal': counts, 'rotation_cos': '107/108', 'rotation_sin': 'sqrt(215)/108'}
Path(__file__).with_name('equality_receipt.json').write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps(out, indent=2))
