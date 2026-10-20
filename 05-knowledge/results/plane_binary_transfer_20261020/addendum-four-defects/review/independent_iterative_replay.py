#!/usr/bin/env python3
"""Retained independent iterative replay used in the addendum review.

No search and no imports from the generator or the primary verifier. Reconstructs
norm triangles directly and checks every recorded proof inference. The default
certificate is ../results/q17_no_two_defects_certificate.json.gz; an optional
first argument selects another certificate. The primary verifier separately
checks the root orbit; the frozen base supplies zero-defect exclusion.

This is the reviewed one-off replay saved as a reusable script, with a file-path
wrapper, certificate identity reporting, and explicit metadata/root checks.
"""
from pathlib import Path
import gzip
import hashlib
import json
import sys


def require(condition):
    if not condition:
        raise ValueError('Independent replay rejected')


p = 17
pts = [(x, y) for x in range(p) for y in range(p)]
adj = []
for x, y in pts:
    adj.append({j for j, (a, b) in enumerate(pts)
                if ((x-a)**2 - 3*(y-b)**2) % 17 == 1})
triangles = sorted((i, j, k) for i in range(289) for j in adj[i] if i < j
                   for k in (adj[i] & adj[j]) if j < k)
require(len(triangles) == 1734)

certificate = (Path(sys.argv[1]) if len(sys.argv) > 1 else
               Path(__file__).resolve().parents[1] / 'results' /
               'q17_no_two_defects_certificate.json.gz')
raw = certificate.read_bytes()
c = json.loads(gzip.decompress(raw))
require(type(c) is dict and set(c) ==
        {'p', 'budget', 'root_monochromatic_edge', 'proof'})
require(type(c['p']) is int and c['p'] == 17)
require(type(c['budget']) is int and c['budget'] == 2)
require(type(c['root_monochromatic_edge']) is list and
        all(type(v) is int for v in c['root_monochromatic_edge']) and
        tuple(c['root_monochromatic_edge']) == triangles[0])
root = [None] * 289
for v in triangles[0]:
    root[v] = 0
stack = [(c['proof'], root)]
nodes = leaves = forces = 0


def is_mono(a, eid):
    values = [a[v] for v in triangles[eid]]
    return values[0] in (0, 1) and values[0] == values[1] == values[2]


while stack:
    t, a = stack.pop()
    nodes += 1
    require(type(t) is dict and type(t.get('forces')) is list)
    for record in t['forces']:
        require(type(record) is list and len(record) == 5)
        v, value, eid, b1, b2 = record
        require(type(v) is int and 0 <= v < 289 and a[v] is None)
        require(type(value) is int and value in (0, 1))
        require(all(type(e) is int and 0 <= e < 1734
                    for e in (eid, b1, b2)))
        require(len({eid, b1, b2}) == 3 and
                is_mono(a, b1) and is_mono(a, b2))
        require(v in triangles[eid])
        others = [w for w in triangles[eid] if w != v]
        require(len(others) == 2 and all(a[w] == 1-value for w in others))
        a[v] = value
        forces += 1
    if 'conflicts' in t:
        require(set(t) == {'forces', 'conflicts'} and
                type(t['conflicts']) is list and len(t['conflicts']) == 3)
        require(all(type(e) is int and 0 <= e < 1734
                    for e in t['conflicts']))
        require(len(set(t['conflicts'])) == 3 and
                all(is_mono(a, e) for e in t['conflicts']))
        leaves += 1
    else:
        require(set(t) == {'forces', 'branch', 'children'} and
                type(t['children']) is list and len(t['children']) == 2)
        v = t['branch']
        require(type(v) is int and 0 <= v < 289 and a[v] is None)
        for value, child in enumerate(t['children']):
            b = a.copy()
            b[v] = value
            stack.append((child, b))

print(json.dumps({
    'independently_implemented_iterative_replay': 'PASS',
    'triangles': len(triangles),
    'nodes': nodes,
    'leaves': leaves,
    'forced_assignments': forces,
    'certificate_sha256': hashlib.sha256(raw).hexdigest(),
}, sort_keys=True))
