#!/usr/bin/env python3
"""Read-only standalone verifier of the explicit 210-vertex certificates."""
from collections import Counter
from itertools import combinations
from pathlib import Path
import hashlib
import json

path=Path(__file__).resolve().parent/'constructor-210.json'
raw=path.read_bytes()
data=json.loads(raw)
def require(ok, message):
    if not ok:
        raise RuntimeError(message)
ports=[tuple(x) for x in data['vertex_order']]
p=data['permutation']
n=len(ports)
require(n==210 and len(set(ports))==210,'carrier count')
require(set(ports)=={(r,b) for r in range(1,210,2) for b in (-1,1)},'carrier values')
require(sorted(p)==list(range(n)),'vertex permutation')
for i,(r,b) in enumerate(ports):
    require(ports[p[i]]==((2*r+b)%210,-b),'actual L action')

def adjacency(bits, size):
    require(len(bits)==size*(size-1)//2 and set(bits)<={'0','1'},'bit encoding')
    a=[[False]*size for _ in range(size)]
    for bit,(u,v) in zip(bits,combinations(range(size),2)):
        a[u][v]=bit=='1'
        a[v][u]=bit=='0'
    return a
def failures(a, permutation):
    return [(u,v) for u,v in combinations(range(len(a)),2)
            if a[u][v] != a[permutation[u]][permutation[v]]]

free=adjacency(data['unrestricted_orientation_bits'],n)
module=adjacency(data['module_orientation_bits'],n)
blocks=[tuple(x) for x in data['blocks']]
require(len(blocks)==105 and sorted(v for block in blocks for v in block)==list(range(n)), 'pair partition')
where={v:i for i,block in enumerate(blocks) for v in block}
for u,v in blocks:
    require(u<v,'block orientation convention')
    x=u
    orbit=[]
    while x not in orbit:
        orbit.append(x)
        x=p[x]
    require(len(orbit)%2==0 and orbit[len(orbit)//2]==v,'antipodal block')
    for w in range(n):
        if w not in (u,v):
            require(module[u][w]==module[v][w],'homogeneous pair')
pi=data['quotient_permutation']
q=adjacency(data['quotient_orientation_bits'],len(blocks))
for i,(u,v) in enumerate(blocks):
    require(where[p[u]]==where[p[v]]==pi[i], 'quotient action')
    for j,(x,y) in enumerate(blocks):
        if i!=j:
            require(q[i][j]==module[u][x]==module[u][y]==module[v][x]==module[v][y], 'quotient lift')
free_fail=failures(free,p)
module_fail=failures(module,p)
q_fail=failures(q,pi)
inner=sum(where[u]==where[v] for u,v in module_fail)
external=len(module_fail)-inner
require((len(free_fail),len(module_fail),len(q_fail),inner,external)==(21,77,14,21,56),'actual defect totals')
require(external==4*len(q_fail),'external weight')
print(json.dumps({'status':'PASS','vertices':n,'blocks':len(blocks),
                  'direct_unrestricted_defect':len(free_fail),'direct_module_defect':len(module_fail),
                  'direct_quotient_defect':len(q_fail),'internal_failures':inner,
                  'external_failures':external,'constructor_sha256':hashlib.sha256(raw).hexdigest()},
                 indent=2,sort_keys=True))
