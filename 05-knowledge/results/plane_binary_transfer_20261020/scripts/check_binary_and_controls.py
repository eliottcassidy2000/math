"""Exact finite sanity checks for Q, Fourier degree and SAT positive controls."""
from fractions import Fraction as F
from itertools import product, combinations
from pathlib import Path
import json
Q = lambda a, b, c: 1 + a * b + a * c + b * c
mono = lambda xs: len(set(xs)) == 1
for xs in product([-1, 1], repeat=3):
    if not Q(*xs) == 4 * mono(xs):
        raise RuntimeError('Exact verification gate failed: scripts/check_binary_and_controls.py:8')
grid = [F(k, 4) for k in range(-4, 5)]
for a, b, c in product(grid, repeat=3):
    q = Q(a, b, c)
    if not q >= 0:
        raise RuntimeError('Exact verification gate failed: scripts/check_binary_and_controls.py:12')
    if not int(mono([x >= 0 for x in [a, b, c]])) <= q:
        raise RuntimeError('Exact verification gate failed: scripts/check_binary_and_controls.py:13')
    if not q / 4 == ((1 + a) * (1 + b) * (1 + c) + (1 - a) * (1 - b) * (1 - c)) / 8:
        raise RuntimeError('Exact verification gate failed: scripts/check_binary_and_controls.py:14')
rows = []
for k in range(2, 9):
    coeff = {}
    for mask in range(1 << k):
        s = [i for i in range(k) if mask >> i & 1]
        coeff[mask] = sum((F(mono(x)) * __import__('math').prod((x[i] for i in s)) for x in product([-1, 1], repeat=k))) / 2 ** k
        if not coeff[mask] == (F(1, 2 ** (k - 1)) if len(s) % 2 == 0 else 0):
            raise RuntimeError('Exact verification gate failed: scripts/check_binary_and_controls.py:21')
    rows.append({'k': k, 'degree': max((mask.bit_count() for mask, c in coeff.items() if c))})
even = [x for x in product([-1, 1], repeat=4) if __import__('math').prod(x) == 1]
odd = [x for x in product([-1, 1], repeat=4) if __import__('math').prod(x) == -1]
for k in range(4):
    for s in combinations(range(4), k):
        mean = lambda pop: sum((F(__import__('math').prod((x[i] for i in s)), len(pop)) for x in pop))
        if not mean(even) == mean(odd):
            raise RuntimeError('Exact verification gate failed: scripts/check_binary_and_controls.py:28')
if not (sum(map(mono, even)) == 2 and sum(map(mono, odd)) == 0):
    raise RuntimeError('Exact verification gate failed: scripts/check_binary_and_controls.py:29')
controls = []
for line in (Path(__file__).resolve().parents[1] / 'results/finite_triangle_probe.jsonl').read_text().splitlines():
    row = json.loads(line)
    if row['status'] != 'SAT':
        continue
    p, d = (row['p'], row['d'])
    colors = row['coloring']
    pts = list(product(range(p), repeat=2))
    norm = lambda i, j: ((pts[i][0] - pts[j][0]) ** 2 - d * (pts[i][1] - pts[j][1]) ** 2) % p
    adj = [{j for j in range(p * p) if i != j and norm(i, j) == 1} for i in range(p * p)]
    es = [(i, j, k) for i in range(p * p) for j in adj[i] if i < j for k in adj[i] & adj[j] if j < k]
    if not len(es) == row['triangle_count']:
        raise RuntimeError('Exact verification gate failed: scripts/check_binary_and_controls.py:39')
    if not all((not mono([colors[v] for v in e]) for e in es)):
        raise RuntimeError('Exact verification gate failed: scripts/check_binary_and_controls.py:40')
    controls.append({'p': p, 'triangles': len(es), 'result': 'PASS'})
print(json.dumps({'result': 'PASS', 'rational_Q_grid_points': len(grid) ** 3, 'Boolean_degrees': rows, 'four_bit_parity_mono_counts': [2, 0], 'positive_controls': controls}, indent=2))
