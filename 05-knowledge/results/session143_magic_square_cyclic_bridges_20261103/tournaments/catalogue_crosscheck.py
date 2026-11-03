#!/usr/bin/env python3
"""Independent matrix-relabeling comparison with McKay's primary catalogues."""
from itertools import permutations,combinations
import json,pathlib
root=pathlib.Path(__file__).resolve().parent
catalogue=json.loads((root/'mckay_catalogue.json').read_text())
ours=json.loads((root/'classification.json').read_text())['classes']
edges=list(combinations(range(5),2))
def canon(bits):
    if len(bits)!=10 or not set(bits)<=set('01'):
        raise RuntimeError('Invalid upper-triangle catalogue record')
    a=[[0]*5 for _ in range(5)]
    for bit,(i,j) in zip(bits,edges):
        a[i][j]=int(bit)
        a[j][i]=1-int(bit)
    return min(sum(a[p[i]][p[j]]<<k for k,(i,j) in enumerate(edges))
               for p in permutations(range(5)))
expected={'all':{r['mask'] for r in ours},
          'self_converse':{r['mask'] for r in ours if r['self_dual']},
          'regular':{r['mask'] for r in ours if r['score']==[2]*5}}
for name in expected:
    records=catalogue[name]['records']
    got={canon(bits) for bits in records}
    if len(got)!=len(records) or got!=expected[name]:
        raise RuntimeError('Catalogue mismatch: '+name)
print('PASS: exact canonical match with McKay catalogues; all12, self-converse8, regular1')
