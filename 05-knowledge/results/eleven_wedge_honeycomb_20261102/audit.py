#!/usr/bin/env python3
"""Independent pattern counts plus exact rational matrix consequences."""
import json,subprocess,sys,random
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from math import prod
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'square_sum_archive_stability_20261102'))
from verify_reduction import det,hp,cycles,matrix_from_pairs
def need(b,s):
    if not b:raise RuntimeError(s)
def decode(s):
    n=ord(s[0])-63;bits=[]
    for c in s[1:]:bits.extend((ord(c)-63>>k)&1 for k in range(5,-1,-1))
    edges=[];k=0
    for j in range(1,n):
        for i in range(j):
            if bits[k]:edges.append((i,j))
            k+=1
    return n,edges
def counts(n,edges):
    # Enumerate actual edge subsets, not the producer's vertex recurrence.
    masks=[(1<<i)|(1<<j) for i,j in edges];m=[1,len(edges),0,0,0,0]
    def rec(start,used,k):
        if k>=2:m[k]+=1
        if k==5:return
        for e in range(start,len(masks)):
            if not masks[e]&used:rec(e+1,used|masks[e],k+1)
    rec(0,0,0)
    wedges=[]
    for e,f in combinations(masks,2):
        if e&f:wedges.append(e|f)
    p21=sum(not (e&w) for e in masks for w in wedges)
    return m,len(wedges),p21
def triangle8():
    stream=subprocess.check_output(['geng','-q','8'],text=True).splitlines()
    minimum=None;witness=[]
    for s in stream:
        n,e=decode(s);m,w,p=counts(n,e);v=512+33*w-32*m[2]+7*m[3]
        if minimum is None or v<minimum:minimum=v;witness=[s]
        elif v==minimum:witness.append(s)
    need(len(stream)==12346 and minimum==212,'triangle8 independent census')
    return {'universe':len(stream),'minimum':minimum,'equality_graph6':witness,'kernel':'explicit subsets of edges, wedge intersection masks'}
def exceptions():
    data=json.loads((HERE/'extension11.json').read_text());rows=[]
    for s in data['exceptional_bases']:
        n,e=decode(s);lo=None;zero=[]
        for mask in range(1<<n):
            f=e+[(i,n) for i in range(n) if mask>>i&1]
            m,w,p=counts(n+1,f)
            val=499*m[1]-320*m[2]+80*m[3]-20*m[4]+5*m[5]+80*p+56*w
            need(val>=0,'exception hostile '+s)
            lo=val if lo is None else min(lo,val)
            if val==0:zero.append(mask)
        deg=[sum(i in edge for edge in e) for i in range(n)]
        rows.append({'graph6':s,'edges':len(e),'degree_sequence':sorted(deg),'all_1024_minimum':lo,'zero_neighborhoods':zero})
    need(len(rows)==19,'exception count')
    return rows
def matrices():
    n=11;es=list(combinations(range(n),2));rng=random.Random(11020740)
    fair=[F(1,2) if i//2==j//2 and j<10 else F(1) for i,j in es]
    cases=[('transitive',[F(1)]*55),('all_fair',[F(1,2)]*55),('five_fair_pairs',fair)]
    for t in range(8):cases.append(('dyadic_'+str(t),[F(rng.randrange(5),4) for _ in es]))
    for t in [0,1,5,9,33,54]:
        p=fair.copy();p[t]=F(7,16) if p[t]==F(1,2) else F(15,16);cases.append(('near_fair_'+str(t),p))
    hostile={(0,1),(1,2),(3,4),(5,6),(7,8),(9,10)}
    cases.append(('P3_plus_4K2',[F(1,2) if e in hostile else F(1) for e in es]))
    rows=[]
    for name,p in cases:
        a=matrix_from_pairs(p,n);H=hp(a);u=det([[a[i][j]+int(i==j) for j in range(n)] for i in range(n)]);v=-det(a)
        V=sum((x*(1-x) for x in p),F(0));C=cycles(a)
        rhs=F(781,1280)*V+F(153,640)*C[3]+C[4]/64+C[5]+C[6]/16+C[7]+C[8]/4+C[9]+C[10]
        slack=H-u-abs(v)-rhs;need(slack>=0,'matrix consequence '+name)
        if name=='five_fair_pairs':need(slack==0 and H==1 and v==0 and V==F(5,4),'fair equality')
        rows.append({'name':name,'H':str(H),'u':str(u),'v':str(v),'V':str(V),'tau':str(C[3]),'slack':str(slack)})
    return rows
if __name__=='__main__':
    mode=sys.argv[1] if len(sys.argv)>1 else 'all'
    out={'status':'passed','scope':'Independent same-session kernels; not an independent-agent or formal audit.'}
    if mode in ['all','graphs']:out.update(triangle8=triangle8(),exceptional_bases=exceptions())
    if mode in ['all','matrices']:out['rational_matrix_controls']=matrices()
    print(json.dumps(out,indent=2,sort_keys=True))
