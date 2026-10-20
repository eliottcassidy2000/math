"""Certificate excluding <=2 monochromatic triangles with one fixed mono root.
Root symmetry is justified using the separately certified impossibility of zero.
"""
from pathlib import Path
from collections import defaultdict
import gzip, hashlib, json
from probe_finite_triangles import field
OUT = Path(__file__).resolve().parents[1] / 'results'
p = 17
d, zeta, pts, U, edges = field(p)
n = p * p
incid = defaultdict(list)
for i, e in enumerate(edges):
    for v in e:
        incid[v].append(i)
nodes = 0

def scan(a):
    bad = []
    units = []
    for i, e in enumerate(edges):
        known = [v for v in e if a[v] >= 0]
        if len(known) == 3 and a[e[0]] == a[e[1]] == a[e[2]]:
            bad.append(i)
        elif len(known) == 2 and a[known[0]] == a[known[1]]:
            units.append((next((v for v in e if a[v] < 0)), 1 - a[known[0]], i))
    return (bad, units)

def dfs(a):
    global nodes
    nodes += 1
    forces = []
    while True:
        bad, units = scan(a)
        if len(bad) > 2:
            return {'forces': forces, 'conflicts': bad[:3]}
        if len(bad) < 2 or not units:
            break
        for v, c, eid in units:
            if a[v] >= 0:
                if a[v] != c:
                    return {'forces': forces, 'conflicts': [bad[0], bad[1], eid]}
                continue
            a[v] = c
            forces.append([v, c, eid, bad[0], bad[1]])
    if all((v >= 0 for v in a)):
        raise AssertionError('Found a <=2-defect coloring')
    if units:
        v = units[0][0]
    else:
        rem = [v for v in range(n) if a[v] < 0]
        v = max(rem, key=lambda v: sum((10 if sum((a[w] >= 0 for w in edges[i])) == 1 else 1 if all((a[w] < 0 for w in edges[i])) else 0 for i in incid[v])))
    children = []
    for color in [0, 1]:
        b = a[:]
        b[v] = color
        children.append(dfs(b))
    return {'forces': forces, 'branch': v, 'children': children}
a = [-1] * n
for v in edges[0]:
    a[v] = 0
cert = {'p': p, 'budget': 2, 'root_monochromatic_edge': edges[0], 'proof': dfs(a)}
raw = json.dumps(cert, separators=(',', ':')).encode()
compressed = gzip.compress(raw, compresslevel=9, mtime=0)
path = OUT / 'q17_no_two_defects_certificate.json.gz'
path.write_bytes(compressed)
print(json.dumps({'nodes': nodes, 'raw_bytes': len(raw), 'compressed_bytes': len(compressed), 'sha256': hashlib.sha256(compressed).hexdigest()}, sort_keys=True))
