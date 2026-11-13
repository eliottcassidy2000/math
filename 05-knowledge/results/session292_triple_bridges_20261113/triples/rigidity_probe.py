#!/usr/bin/env python3
"""Exhaustive physical single-entry seam search, fixed ports0,1 and fixed exterior incidences."""
import itertools,json
from pathlib import Path

def add(a,b):
    c=a.copy()
    for m,v in b.items():c[m]=c.get(m,0)+v
    return c
def mul(a,b):
    c={}
    for m,u in a.items():
        for n,v in b.items():
            if not m&n:c[m|n]=c.get(m|n,0)+u*v
    return c
rows=[]
for n in range(3,7):
    pairs=list(itertools.combinations(range(n),2));groups={};fiber_forward={};fiber_reverse={};h=(1<<n)-1
    target=h^1 # x_I x_b, a=0,b=1
    samples=[];collisions=0;distinct_coeff=0
    for code in range(1<<len(pairs)):
        A=[[0]*n for _ in range(n)]; q=[[{} for _ in range(n+2)] for _ in range(n+2)]
        exterior=0
        for k,(i,j) in enumerate(pairs):
            bit=code>>k&1
            if i<2:exterior|=bit<<k
            if bit:A[i][j]=1;q[i][j]={1<<j:1}
            else:A[j][i]=1;q[j][i]={1<<i:1}
        for i in range(n):q[n][i]={1<<i:1};q[i][n+1]={0:1}
        live=list(range(n+2))
        for k in range(2,n):
            live.remove(k)
            for i in live:
                for j in live:q[i][j]=add(q[i][j],mul(q[i][k],q[k][j]))
        coeff=q[0][1].get(target,0)
        signature=[]
        for i in live:
            for j in live:
                signature.append(tuple(sorted((m,c) for m,c in q[i][j].items() if True)))
        key=(exterior,tuple(sorted((m,c) for m,c in q[n][n+1].items() if m.bit_count()==3)),tuple((p,tuple(sorted((m,c) for m,c in q[p][n+1].items() if m.bit_count()==2))) for p in (0,1)))
        cyclic=tuple((i,j,k) for i,j,k in itertools.combinations(range(n),3) if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[k][j] and A[i][k]))
        fiber=(exterior,cyclic)
        if key in fiber_forward and fiber_forward[key]!=fiber: raise RuntimeError('observer does not fix relative C3 fibre')
        if fiber in fiber_reverse and fiber_reverse[fiber]!=tuple(signature): raise RuntimeError('relative C3 does not fix boundary')
        fiber_forward[key]=fiber;fiber_reverse[fiber]=tuple(signature)
        if key in groups:
            oldcode,oldc,oldsig=groups[key];collisions+=1
            if oldsig!=tuple(signature): raise RuntimeError("boundary reconstruction failed")
            if coeff!=oldc:
                distinct_coeff+=1
                if len(samples)<10:samples.append({'codes':[oldcode,code],'coefficients':[oldc,coeff]})
        else:groups[key]=(code,coeff,tuple(signature))
    rows.append({'n':n,'tournaments':1<<len(pairs),'distinct_signatures':len(groups),'relative_C3_fibres':len(fiber_reverse),'signature_collisions':collisions,'collisions_with_changed_target_coefficient':distinct_coeff,'examples':samples})
    print(rows[-1],flush=True)
Path(__file__).with_name('rigidity_probe.json').write_text(json.dumps({'scope':'All labelled tournaments n=3..6; ports0,1; fixed exterior incidences/port edge; group by degree3 Q_st and degree2 Q_0t/Q_1t, check complete boundary and augmented relative C3 iff.','rows':rows},indent=2)+'\n')
