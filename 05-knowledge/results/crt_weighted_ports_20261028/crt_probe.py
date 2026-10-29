#!/usr/bin/env python3
"""Exact squarefree-odd-factor CRT controls; explicit universe, no assertions."""
from array import array
from collections import Counter
import itertools
import json
import math
import random

CHECKS=0
def ck(ok,msg):
    global CHECKS
    CHECKS+=1
    if not ok:raise RuntimeError(msg)

def lucas(n):
    a,b=2,1
    for _ in range(n):a,b=b,a+b
    return a

def factor(n):
    out=[]
    for p in range(2,math.isqrt(n)+1):
        if n%p==0:
            e=0
            while n%p==0:n//=p;e+=1
            out.append((p,e))
    if n>1:out.append((n,1))
    return out

def make(n):
    N=lucas(n); fields=[]
    for p,e in factor(N):
        if p==2:
            ck(e==2,'F4 multiplicity');fields.append((4,2))
        else:
            ck(e==1,'this finite CRT implementation excludes repeated odd prime factors')
            roots=[a for a in range(p) if (a*a-a-1)%p==0 and pow(a,n,p)==1]
            ck(len(roots)==1,'unique period-compatible split root')
            fields.append((p,roots[0]))
    qs=[q for q,a in fields]
    strides=[];v=1
    for q in qs:strides.append(v);v*=q
    ck(v==N,'CRT cardinality')
    powers=[[([1,2,3][i%3] if q==4 else pow(a,i,q)) for i in range(n)] for q,a in fields]
    reverse=[[row[-i%n] for i in range(n)] for row in powers]
    R=array('i',[-1])*N; W=array('I',[0])*N
    def put(i,first,previous,x,y,w):
        if i==n:
            ix=sum(a*s for a,s in zip(x,strides));iy=sum(a*s for a,s in zip(y,strides))
            ck(R[ix]==-1,'no golden-reader collision')
            R[ix]=iy;W[ix]=w;return
        put(i+1,first,False,x,y,w)
        if not previous and not(i==n-1 and first):
            xx=[(a^powers[j][i]) if qs[j]==4 else (a+powers[j][i])%qs[j] for j,a in enumerate(x)]
            yy=[(a^reverse[j][i]) if qs[j]==4 else (a+reverse[j][i])%qs[j] for j,a in enumerate(y)]
            put(i+1,first or i==0,True,xx,yy,w|(1<<i))
    put(0,False,False,[0]*len(qs),[0]*len(qs),0)
    ck(all(v>=0 for v in R),'reader surjectivity')
    ck(all(R[R[x]]==x for x in range(N)),'literal word reflection involution')
    return qs,strides,R,W

def projection(x,chosen,qs,strides):
    v=0;s=1
    for j in chosen:v+=(x//strides[j]%qs[j])*s;s*=qs[j]
    return v

def census(n):
    qs,strides,R,W=make(n);N=len(R);stable=[];rows=[]
    for j,q in enumerate(qs):
        out=[set() for _ in range(q)]
        for x,y in enumerate(R):out[x//strides[j]%q].add(y//strides[j]%q)
        rows.append(dict(field=q,output_ambiguity_histogram=dict(sorted(Counter(map(len,out)).items()))))
    for size in range(1,len(qs)+1):
        for chosen in itertools.combinations(range(len(qs)),size):
            d=math.prod(qs[j] for j in chosen); table=array('i',[-1])*d;ok=True
            for x,y in enumerate(R):
                a=projection(x,chosen,qs,strides);b=projection(y,chosen,qs,strides)
                if table[a]==-1:table[a]=b
                elif table[a]!=b:ok=False;break
            if ok:stable.append([qs[j] for j in chosen])
    # CRT factor influence is tested on every fibre, not random perturbations.
    influence=[]
    for j,q in enumerate(qs):
        other=[i for i in range(len(qs)) if i!=j];d=N//q
        for i,p in enumerate(qs):
            table=array('i',[-1])*d;changes=False
            for x,y in enumerate(R):
                a=projection(x,other,qs,strides);b=y//strides[i]%p
                if table[a]==-1:table[a]=b
                elif table[a]!=b:changes=True;break
            if changes:influence.append([q,p])
    return dict(n=n,states=N,fields=qs,stable_factor_projections=stable,
                factor_output_ambiguities=rows,dependency_edges=influence)

def probes(n,q,shifts):
    qs,strides,R,W=make(n);j=qs.index(q)
    # Integer translation in every residue field (F4 prime subfield addition).
    def translated(x,t):
        out=0
        for a,p,s in zip([x//s%p for p,s in zip(qs,strides)],qs,strides):
            out+=((a^(t%2)) if p==4 else (a+t)%p)*s
        return out
    signatures=[]
    for x in range(len(R)):
        signatures.append(tuple([x//strides[j]%q]+[R[translated(x,t)]//strides[j]%q for t in shifts]))
    histogram=Counter(Counter(signatures).values())
    return dict(n=n,observed_field=q,integer_shifts=shifts,signature_count=len(set(signatures)),
                fibre_histogram=dict(sorted(histogram.items())))

def decoder(n,q):
    qs,strides,R,W=make(n);j=qs.index(q);N=len(R)
    rng=random.Random(243+n+q)
    candidates=[tuple((t%2 if p==4 else t%p) for p in qs) for t in range(8)]
    candidates+= [tuple(rng.randrange(p) for p in qs) for _ in range(40)]
    candidates=list(dict.fromkeys(candidates))
    columns=[]
    for shift in candidates:
        col=array('I')
        for x in range(N):
            y=sum((((x//s%p)^a) if p==4 else (x//s%p+a)%p)*s for a,p,s in zip(shift,qs,strides))
            col.append(R[y]//strides[j]%q)
        columns.append(col)
    ids=[x//strides[j]%q for x in range(N)];used=[];steps=[]
    while len(set(ids))<N:
        best=None
        for c,col in enumerate(columns):
            if c in used:continue
            count=Counter(zip(ids,col));score=sum(v*(v-1)//2 for v in count.values())
            candidate=(score,-len(count),c)
            if best is None or candidate<best:best=candidate
        ck(best is not None,'probe menu not exhausted')
        score,negative,c=best;used.append(c);table={};new=[]
        for pair in zip(ids,columns[c]):
            if pair not in table:table[pair]=len(table)
            new.append(table[pair])
        ck(len(table)>len(set(ids)),'probe menu makes strict progress')
        ids=new;steps.append(dict(shift=list(candidates[c]),cells=len(table),unresolved_pairs=score))
    return dict(n=n,field=q,candidate_count=len(candidates),steps=steps,
                status='FINITE-EXACT feasible decoder; no optimality claim')

def main():
    rows=[census(n) for n in [9,15,21,23,25,27]]
    probe_rows=[probes(n,q,ss) for n,q in [(9,19),(15,11),(15,31),(21,29)]
                for ss in [[0],[0,1],[0,1,2],[0,1,2,3]]]
    decoders=[decoder(n,q) for n,q in [(15,11),(15,31),(21,29)]]
    print(json.dumps(dict(scope='all legal cyclic words at the listed odd lengths; factors exactly as listed',
                         rows=rows,probes=probe_rows,decoders=decoders,active_checks=CHECKS),indent=2,sort_keys=True))

if __name__=='__main__':main()
