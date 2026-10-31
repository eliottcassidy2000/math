#!/usr/bin/env python3
"""Independent exact factor matrix and fraction-free determinant certificate."""
from pathlib import Path
import json,sys
from audit_fractions import factor,require
ROOT=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent.parent
reported=json.loads((ROOT/'pure_blue_fractran/geometry.json').read_text())
rows=[]
for line in (ROOT/'pure_blue_fractran/pure_blue.fractran').read_text().splitlines():
    line=line.split('#')[0].strip()
    if line:
        a,b=map(int,line.split('/'));rows.append((factor(a),factor(b)))
ps=reported['matrix_row_primes']
require(set(ps)==set().union(*(set(a)|set(b) for a,b in rows)),'Unknown/missing prime rows')
M=[[a[p]-b[p] for a,b in rows] for p in ps]
require(M==reported['matrix'],'Valuation matrix mismatch')
controls=[int(p>=23) for p in ps]
require(all(sum(w*M[i][j] for i,w in enumerate(controls))==0 for j in range(22)),'Control invariant fails')
B=[list(row) for p,row in zip(ps,M) if p!=23]
prev=1;sign=1
for k in range(len(B)-1):
    pivot=next((j for j in range(k,len(B)) if B[j][k]),None)
    require(pivot is not None,'Singular minor')
    if pivot!=k:B[pivot],B[k]=B[k],B[pivot];sign=-sign
    x=B[k][k]
    for i in range(k+1,len(B)):
        for j in range(k+1,len(B)):
            num=x*B[i][j]-B[i][k]*B[k][j]
            require(num%prev==0,'Nonexact Bareiss division')
            B[i][j]=num//prev
        B[i][k]=0
    prev=x
det=sign*B[-1][-1]
require(det==reported['minor_determinant']==1,'Wrong unimodular minor')

def strong(vertices,edges):
    for start in vertices:
        seen={start};todo=[start]
        while todo:
            v=todo.pop()
            for a,b in edges:
                if a==v and b not in seen:seen.add(b);todo.append(b)
        require(seen==set(vertices),'Not strongly connected')
controls_set={p for p in ps if p>=23}
control_edges=[];data_edges=[]
for pos,neg in rows:
    qa=next(p for p in pos if p in controls_set);qb=next(p for p in neg if p in controls_set)
    control_edges.append((qb,qa))
    for a in neg:
        for b in pos:
            if a<19 and b<19:data_edges.append((a,b))
strong(controls_set,control_edges);strong([2,3,5,7,11,13,17],data_edges)
weights={2:12,3:10,5:8,7:6,11:18,13:16,17:14,19:8,
         23:0,29:-7,31:-8,37:-6,41:-7,43:-5,47:-6,53:-4,59:-5,
         61:-3,67:-4,71:-2,73:-3,79:-1,83:-2}
require(all(p**12>=2**w for p,w in weights.items() if w>=0), 'Magnitude envelope failed')
require(all(sum(weights[p]*M[i][j] for i,p in enumerate(ps))==1 for j in range(22)), 'Microstep clock mismatch')
checkpoints=json.loads((ROOT/'pure_blue_fractran/verification.json').read_text())['checkpoints']
for c in checkpoints:
    w=c['window'];n=c['n']
    require(sum(x*y for x,y in zip(w,[12,10,8,6,18,16]))+8*n-108==c['rational_steps_from_seed'],'Cumulative clock mismatch')
print(json.dumps(dict(status='PASS',integer_magnitude_envelope_power=12,microstep_clock_checked_fractions=22,clock_checkpoint_checks=len(checkpoints),minor_deleted_prime=23,minor_determinant=det,
                     rank=22,left_kernel_dimension=1,right_kernel_dimension=0,
                     no_nonempty_periodic_orbit_from_any_positive_integer=True),indent=2))
