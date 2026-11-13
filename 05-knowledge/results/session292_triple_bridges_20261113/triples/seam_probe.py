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
    pairs=list(itertools.combinations(range(n),2));groups={};h=(1<<n)-1
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
                signature.append(tuple(sorted((m,c) for m,c in q[i][j].items() if not (i==0 and j==1 and m==target))))
        key=(exterior,tuple(signature))
        if key in groups:
            oldcode,oldc=groups[key];collisions+=1
            if coeff!=oldc:
                distinct_coeff+=1
                if len(samples)<10:samples.append({'codes':[oldcode,code],'coefficients':[oldc,coeff]})
        else:groups[key]=(code,coeff)
    rows.append({'n':n,'tournaments':1<<len(pairs),'distinct_signatures':len(groups),'signature_collisions':collisions,'collisions_with_changed_target_coefficient':distinct_coeff,'examples':samples})
    print(rows[-1],flush=True)
Path(__file__).with_name('seam_probe.json').write_text(json.dumps({'scope':'All labelled tournaments n=3..6; ports0,1; fixed exterior incidences and port edge; source/sink included in boundary equality.','rows':rows},indent=2)+'\n')
