"""Independent literal-state audit of the frozen selective-star kernel."""
from pathlib import Path
import hashlib
import json
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
def need(ok,msg):
    if not ok: raise RuntimeError(msg)

def paths(rows):
    n=len(rows); full=(1<<n)-1
    e=[[0]*n for _ in range(full+1)]; h=[1]+[0]*full
    for mask in range(1,full+1):
        for v in range(n):
            if mask>>v&1:
                rest=mask^(1<<v)
                e[mask][v]=sum(e[rest][u] for u in range(n) if rows[u]>>v&1) if rest else 1
        h[mask]=sum(e[mask])
    return e,h

def state(rows,r):
    full=(1<<len(rows))-1
    e,h=paths(rows)
    s,_=paths([full^(1<<v)^rows[v] for v in range(len(rows))])
    f=sum(e[m][r]*h[full^m] for m in range(1,full+1))
    return [h[full],s[full][r],h[full^(1<<r)],f]

a=[518,388,456,1011,487,455,387,257,1,502]
p=[130,388,457,1011,487,455,387,256,1,503]
t=[1154,388,457,1011,487,455,387,256,1,503,1022]
patch=a[:]
for v in [2,7,9]:
    patch[0]^=1<<v
    patch[v]^=1
need(patch==p,'exact three-star patch')
need(t==[row|((1<<10) if v==0 else 0) for v,row in enumerate(p)]+[1022], 'physical newborn')
states=[state(a,0),state(p,0),state(t,10)]
need(states[0][0]==881,'source count')
need(states[1]==[153,1,27,1539],'literal parent')
need(states[2]==[1691,152,153,1692],'literal child')
need(states[1][0]+states[1][3]-states[1][1]==states[2][0], 'full target')
need(states[2][3]-states[2][0]==states[1][1]==1,'readiness')
base=HERE.parent/'star_kernel'
normal=json.loads((base/'receipt.json').read_text(encoding='utf-8'))
optimized=json.loads((base/'receipt.optimized.json').read_text(encoding='utf-8'))
need(normal==optimized,'inherited normal optimized pair')
pins=[{'path':str(p.relative_to(ROOT)).replace('\\','/'),
       'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
      for p in [base/'PROOF.md',base/'check.py',base/'receipt.json',base/'receipt.optimized.json']]
out={'status':'VERIFIED proof review plus FINITE-EXACT independent literal controls',
     'scope':'three literal graphs and exact repair; exhaustive star menus inherited from reviewed root checker',
     'literal_states':states,'source_pins':pins}
suffix='.optimized' if sys.flags.optimize else ''
(HERE/('star_kernel_audit'+suffix+'.json')).write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(states))
