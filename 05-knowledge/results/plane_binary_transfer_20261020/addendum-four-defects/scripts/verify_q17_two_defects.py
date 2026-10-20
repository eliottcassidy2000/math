"""Independent exact verifier of <=2-mono refutation and root symmetry.
Uses direct norm adjacency, never the certificate generator's solver.
"""
from pathlib import Path
import gzip, json, hashlib, sys
P = Path(__file__).resolve().parents[1] / 'results'
certificate_path = Path(sys.argv[1]) if len(sys.argv) > 1 else P / 'q17_no_two_defects_certificate.json.gz'
raw = certificate_path.read_bytes()
c = json.loads(gzip.decompress(raw))
p = 17
d = 3
if not (type(c['p']) is int and c['p'] == p and type(c['budget']) is int and c['budget'] == 2):
    raise AssertionError('Certificate check failed')
pts = [(i, j) for i in range(p) for j in range(p)]
ix = {z: i for i, z in enumerate(pts)}
norm = lambda i, j: ((pts[i][0] - pts[j][0]) ** 2 - d * (pts[i][1] - pts[j][1]) ** 2) % p
adj = [{j for j in range(p * p) if j != i and norm(i, j) == 1} for i in range(p * p)]
edges = sorted(((i, j, k) for i in range(p * p) for j in adj[i] if i < j for k in adj[i] & adj[j] if j < k))
if not (type(c['root_monochromatic_edge']) is list and all(type(v) is int for v in c['root_monochromatic_edge']) and tuple(c['root_monochromatic_edge']) == edges[0]):
    raise AssertionError('Certificate check failed')
add = lambda z, w: ((z[0] + w[0]) % p, (z[1] + w[1]) % p)
mul = lambda z, w: ((z[0] * w[0] + d * z[1] * w[1]) % p, (z[0] * w[1] + z[1] * w[0]) % p)
units = [z for z in pts if (z[0] * z[0] - d * z[1] * z[1]) % p == 1]
orbit = {tuple(sorted((ix[add(t, mul(u, pts[v]))] for v in edges[0]))) for t in pts for u in units}
if not orbit == set(edges):
    raise AssertionError('Certificate check failed')
nodes = leaves = forces = 0

def index(value, size):
    if type(value) is not int or not 0 <= value < size:
        raise ValueError('Noncanonical certificate index')

def mono(a, eid):
    index(eid, len(edges))
    v = [a[x] for x in edges[eid]]
    return v[0] is not None and v[0] == v[1] == v[2]

def verify(t, a):
    global nodes, leaves, forces
    nodes += 1
    if type(t) is not dict or type(t.get('forces')) is not list:
        raise ValueError('Malformed proof node')
    for record in t['forces']:
        if type(record) is not list or len(record) != 5:
            raise ValueError('Malformed force record')
        v, value, eid, b1, b2 = record
        index(v, len(a))
        index(value, 2)
        for edge_id in [eid, b1, b2]:
            index(edge_id, len(edges))
        if not (value in [0, 1] and a[v] is None and (len({eid, b1, b2}) == 3)):
            raise AssertionError('Certificate check failed')
        if not (mono(a, b1) and mono(a, b2)):
            raise AssertionError('Certificate check failed')
        e = edges[eid]
        if not v in e:
            raise AssertionError('Certificate check failed')
        if not all((a[w] == 1 - value for w in e if w != v)):
            raise AssertionError('Certificate check failed')
        a[v] = value
        forces += 1
    if 'conflicts' in t:
        if type(t['conflicts']) is not list or len(t['conflicts']) != 3:
            raise ValueError('Malformed conflict list')
        for edge_id in t['conflicts']:
            index(edge_id, len(edges))
        if not (set(t) == {'forces', 'conflicts'} and len(set(t['conflicts'])) == 3):
            raise AssertionError('Certificate check failed')
        if not all((mono(a, eid) for eid in t['conflicts'])):
            raise AssertionError('Certificate check failed')
        leaves += 1
        return
    if not (set(t) == {'forces', 'branch', 'children'} and type(t['children']) is list and len(t['children']) == 2):
        raise AssertionError('Certificate check failed')
    v = t['branch']
    index(v, len(a))
    if not a[v] is None:
        raise AssertionError('Certificate check failed')
    for color, child in enumerate(t['children']):
        b = a[:]
        b[v] = color
        verify(child, b)
a = [None] * (p * p)
for v in edges[0]:
    a[v] = 0
verify(c['proof'], a)
print(json.dumps({'result': 'PASS', 'excluded_budget': 2, 'root_orbit_triangles': len(orbit), 'nodes': nodes, 'leaves': leaves, 'forced_assignments': forces, 'certificate_sha256': hashlib.sha256(raw).hexdigest(), 'consequence': 'Together with the separately verified zero-defect exclusion and even cut parity, every coloring has >=4 mono triangles.'}, sort_keys=True))
