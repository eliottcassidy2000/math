from itertools import permutations, product
from pathlib import Path
import json

def need(ok, message):
    if not ok:
        raise RuntimeError(message)

def words(rows, vertices=None):
    vertices = list(range(len(rows))) if vertices is None else list(vertices)
    return [list(p) for p in permutations(vertices)
            if all(rows[a] >> b & 1 for a,b in zip(p,p[1:]))]

q = [8,1,3,6]
a, b = words(q,[0,1,3]), words(q,[0,2,3])
need(len(a)==len(b)==3 and len(words(q))==5, 'two-port hostile')
need(q[0] >> 3 & 1, 'shared arc orientation')
grafts=[]
for orientations in product((0,1),repeat=4):
    rows=[0]*5
    for u,v in ((0,1),(1,2),(2,0),(0,3),(3,4),(4,0)):
        rows[u] |= 1<<v
    for direction,(u,v) in zip(orientations,((1,3),(1,4),(2,3),(2,4))):
        if direction:
            u,v=v,u
        rows[u] |= 1<<v
    h=len(words(rows))
    need(h>=9,'one-port positive control')
    grafts.append({'cross_bits':orientations,'H':h})
report={'status':'FINITE_EXACT_CONTROLS_FOR_SCOPED_MERGE_PROOF',
        'Q':q,'Q_paths':words(q),'triangle_paths':[a,b],
        'one_port_grafts':grafts,
        'scope':'All16 orientations between the two disjoint vertex pairs of two fixed triangles meeting at one vertex; one literal two-port hostile.'}
Path(__file__).with_name('overlap.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
print('PASS:16 one-port grafts and exact two-port hostile')
