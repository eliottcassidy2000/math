"""Independent exact verifier: reconstruct all norm triangles, then check tree."""
from pathlib import Path
import json, hashlib, argparse
parser=argparse.ArgumentParser()
parser.add_argument("--certificate",type=Path)
args=parser.parse_args()
P = Path(__file__).resolve().parents[1] / 'results/q17_nae_certificate.json'
if args.certificate is not None:
    P=args.certificate
raw = P.read_bytes()
c = json.loads(raw)
p = c['p']
d = c['d']
if not ((p, d) == (17, 3) and pow(d, (p - 1) // 2, p) == p - 1):
    raise RuntimeError('Exact verification gate failed: scripts/verify_q17.py:6')
vertices = [tuple(z) for z in c['vertices']]
if not vertices == [(i, j) for i in range(p) for j in range(p)]:
    raise RuntimeError('Exact verification gate failed: scripts/verify_q17.py:8')

def normdiff(i, j):
    a, b = vertices[i]
    x, y = vertices[j]
    return ((a - x) ** 2 - d * (b - y) ** 2) % p
adj = [set((j for j in range(p * p) if j != i and normdiff(i, j) == 1)) for i in range(p * p)]
edges = sorted(((i, j, k) for i in range(p * p) for j in adj[i] if i < j for k in adj[i] & adj[j] if j < k))
if not edges == [tuple(e) for e in c['edges']]:
    raise RuntimeError('Exact verification gate failed: scripts/verify_q17.py:15')
if not (len(edges) == 1734 and all((len(a) == 18 for a in adj))):
    raise RuntimeError('Exact verification gate failed: scripts/verify_q17.py:16')
def checked_index(value, bound, label):
    if type(value) is not int or not 0 <= value < bound:
        raise RuntimeError('Invalid certificate '+label)

def validate_tree_shape(t):
    if type(t) is not dict or type(t.get('forces')) is not list:
        raise RuntimeError('Invalid certificate node/forces')
    for item in t['forces']:
        if type(item) is not list or len(item)!=3:
            raise RuntimeError('Invalid certificate force record')
        checked_index(item[0],p*p,'force vertex')
        checked_index(item[1],2,'force color')
        checked_index(item[2],len(edges),'force edge')
    if 'conflict' in t:
        if set(t)!={'forces','conflict'}:
            raise RuntimeError('Invalid conflict keys')
        checked_index(t['conflict'],len(edges),'conflict edge')
    else:
        if set(t)!={'forces','branch','children'} or type(t['children']) is not list or len(t['children'])!=2:
            raise RuntimeError('Invalid branch structure')
        checked_index(t['branch'],p*p,'branch vertex')
        for child in t['children']:
            validate_tree_shape(child)

validate_tree_shape(c['proof'])

nodes = leaves = forces = 0

def verify(t, a):
    global nodes, leaves, forces
    nodes += 1
    for v, value, eid in t['forces']:
        if not (value in [0, 1] and a[v] is None):
            raise RuntimeError('Exact verification gate failed: scripts/verify_q17.py:23')
        e = edges[eid]
        if not v in e:
            raise RuntimeError('Exact verification gate failed: scripts/verify_q17.py:24')
        other = [w for w in e if w != v]
        if not all((a[w] == 1 - value for w in other)):
            raise RuntimeError('Exact verification gate failed: scripts/verify_q17.py:26')
        a[v] = value
        forces += 1
    if 'conflict' in t:
        if not set(t) == {'forces', 'conflict'}:
            raise RuntimeError('Exact verification gate failed: scripts/verify_q17.py:29')
        vals = [a[w] for w in edges[t['conflict']]]
        if not (vals[0] is not None and vals[0] == vals[1] == vals[2]):
            raise RuntimeError('Exact verification gate failed: scripts/verify_q17.py:31')
        leaves += 1
        return
    if not (set(t) == {'forces', 'branch', 'children'} and len(t['children']) == 2):
        raise RuntimeError('Exact verification gate failed: scripts/verify_q17.py:33')
    v = t['branch']
    if not a[v] is None:
        raise RuntimeError('Exact verification gate failed: scripts/verify_q17.py:34')
    for val, child in enumerate(t['children']):
        b = a[:]
        b[v] = val
        verify(child, b)
if not c['root_fixed'] == [0, 0]:
    raise RuntimeError('Exact verification gate failed: scripts/verify_q17.py:38')
a = [None] * (p * p)
a[0] = 0
verify(c['proof'], a)
print(json.dumps({'result': 'PASS', 'field': 'F_17[sqrt(3)]', 'vertices': len(vertices), 'undirected_unit_edges': sum(map(len, adj)) // 2, 'triangle_edges': len(edges), 'nodes': nodes, 'leaves': leaves, 'forced_assignments': forces, 'certificate_sha256': hashlib.sha256(raw).hexdigest()}, sort_keys=True))
