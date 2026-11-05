#!/usr/bin/env python3
import sys,json,random
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from verify import *
from verify_extension import clone,skeleton_operator
from verify_rooted_composition import scalar_shifts

def converse(a):return tuple(sum(1<<j for j in range(len(a)) if a[j]>>i&1) for i in range(len(a)))
def source(S,a):
    t=len(S);T=list(join(converse(S),join((0,),S)));u=a;v=t+1+a
    require(T[u]>>v&1,'initial forward crossing');T[u]^=1<<v;T[v]|=1<<u
    validate(T);return tuple(T),t

def response(S,a):
    n=len(S);d=[[0]*n for _ in range(1<<n)]
    for v in range(n):d[1<<v][v]=1
    for mask in range(1,1<<n):
        for v in range(n):
            if not mask>>v&1:continue
            prior=mask^(1<<v)
            if prior:d[mask][v]=sum(d[prior][u] for u in range(n) if S[u]>>v&1)
    H=[sum(row) for row in d];H[0]=1
    return H[-1],sum(d[mask][a]*H[((1<<n)-1)^mask] for mask in range(1<<n))

def run():
    rng=random.Random(197871);cases=[]
    for n in range(1,4):
        for mask in range(1<<(n*(n-1)//2)):
            for a in range(n):cases.append((tournament(n,mask),a))
    for n,reps in ((4,12),(5,3)):
        for _ in range(reps):cases.append((tournament(n,rng.getrandbits(n*(n-1)//2)),rng.randrange(n)))
    physical=[];operator_checks=0
    for S,a in cases:
        h,q=response(S,a);T,r=source(S,a)
        if len(S)<=3:
            op=skeleton_operator(T,r)
            require({k:v for k,v in scalar_shifts(op,1).items() if k}=={1:h*h,2:q*q},('full scalar operator',S,a,h,q));operator_checks+=1
        for m in (1,2):
            A=clone(T,r,m);H=value(packing_polynomial(A),2)
            require(H==h*h+q*q*2**m,'independent odd-cycle physical count')
            physical.append((len(S),a,m,H,h,q))
    matrices=0
    for n in range(5):
        entries=list((i,j) for i in range(n) for j in range(i,n))
        for mask in range(1<<len(entries)):
            edges=set()
            for b,(i,j) in enumerate(entries):
                if mask>>b&1:edges.update(((i,j),(j,i)))
            rankone=bool(edges) and not any(i!=u and j!=v for i,j in edges for u,v in edges)
            require(rankone==(len(edges)==1 and next(iter(edges))[0]==next(iter(edges))[1]),'independent symmetric rank1');matrices+=1
    classified=0
    for n in range(1,6):
        S=transitive(n)
        for j in range(n):
            h,q=response(S,j);T,r=source(S,j)
            require(h==1 and q==2**j,'unit response')
            order=list(range(n-1,j,-1))+list(range(j,-1,-1))+[r]+list(range(n+1,n+j+2))+list(range(n+j+2,2*n+1))
            first=n-1-j;last=first+2*j+2
            require(sorted(order)==list(range(2*n+1)),'core order')
            for i in range(len(order)):
                for k in range(i+1,len(order)):
                    require(bool(T[order[i]]>>order[k]&1)!=(i==first and k==last),'single backward core')
            classified+=1
    residue_checks=0
    for h in range(7):
        for q in range(7):
            for m in range(1,7):
                if (h*h+2**m*q*q)%7==0:require(h==q==0,'binary norm obstruction')
                residue_checks+=1
    print(json.dumps({'rooted_side_sources':len(cases),'independent_odd_cycle_graphs':len(physical),'scalar_operator_checks':operator_checks,'symmetric_cross_matrices':matrices,'unit_core_controls':classified,'residue_controls':residue_checks,'physical_controls':physical,'status':'PASS'},sort_keys=True,indent=2))
if __name__=='__main__':run()
