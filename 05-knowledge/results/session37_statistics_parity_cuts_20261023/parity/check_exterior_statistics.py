#!/usr/bin/env python3
"""Exact, standard-library-only controls for the Session37 exterior construction."""
from collections import Counter
from itertools import combinations, permutations
from math import factorial, comb, gcd
from fractions import Fraction
import json

def need(p, message):
    if not p:
        raise RuntimeError(message)

def sign(p):
    return -1 if sum(p[i] > p[j] for i in range(len(p)) for j in range(i+1,len(p))) % 2 else 1

def pair_action(n, p):
    pairs=list(combinations(range(n),2)); lookup={e:i for i,e in enumerate(pairs)}
    target=[]; flips=[]
    for u,v in pairs:
        a,b=p[u],p[v]
        target.append(lookup[tuple(sorted((a,b)))])
        flips.append(1 if a<b else -1)
    return target,flips

def cycle_data(target, flips):
    seen=set(); out=[]
    for i in range(len(target)):
        if i in seen: continue
        j=i; length=0; s=1
        while j not in seen:
            seen.add(j);length+=1;s*=flips[j];j=target[j]
        out.append((length,s))
    return out

def det_poly(target,flips):
    """det(I+tR) as the exact product over signed permutation cycles."""
    a=[1]+[0]*len(target)
    for k,s in cycle_data(target,flips):
        c=-s*((-1)**k)
        for j in range(len(target),k-1,-1): a[j]+=c*a[j-k]
    return a

def vertex_odd_order(p):
    return all(k%2 for k,_ in cycle_data(p,[1]*len(p)))

def transform(mask,target):
    out=0
    while mask:
        b=mask&-mask;i=b.bit_length()-1;mask-=b
        out|=1<<target[i]
    return out

def connected(n,mask):
    if n==0:return False
    pairs=list(combinations(range(n),2)); seen={0}
    while True:
        old=len(seen)
        for i,(u,v) in enumerate(pairs):
            if (mask>>i)&1 and ((u in seen) or (v in seen)):seen.update((u,v))
        if len(seen)==old:return len(seen)==n

def graph_atlas(n):
    m=comb(n,2);actions=[pair_action(n,p) for p in permutations(range(n))]
    unseen=bytearray([1])*(1<<m); rows=[]
    for mask in range(1<<m):
        if not unseen[mask]:continue
        edges=[i for i in range(m) if (mask>>i)&1]
        ext=True;royle=True;aut=0
        for target,flips in actions:
            other=transform(mask,target);unseen[other]=0
            if other!=mask:continue
            aut+=1
            rev=1
            for e in edges:rev*=flips[e]
            edge_order=[target[e] for e in edges]
            if rev<0:royle=False
            if rev*sign(edge_order)<0:ext=False
        rows.append(dict(mask=mask,edges=len(edges),connected=connected(n,mask),exterior=ext,royle=royle,aut=aut))
    return rows

def mul(a,b,N):
    out=[0]*(N+1)
    for i,x in enumerate(a):
        for j,y in enumerate(b[:N+1-i]):out[i+j]+=x*y
    return out

def main():
    N=7;character_rows=[];polys={};permutation_checks=0
    for n in range(N+1):
        accum=[0]*(comb(n,2)+1);tour=0
        for p in permutations(range(n)):
            target,flips=pair_action(n,p);poly=det_poly(target,flips)
            fixed=2**len(cycle_data(target,flips)) if vertex_odd_order(p) else 0
            need(sum(poly)==fixed,('character mismatch',n,p,poly,fixed))
            if n>=3:need(sum((-1)**k*x for k,x in enumerate(poly))==0,('missing +1 eigenvalue',n,p))
            accum=[a+b for a,b in zip(accum,poly)];tour+=fixed;permutation_checks+=1
        need(all(x%factorial(n)==0 for x in accum),'nonintegral invariant dimension')
        poly=[x//factorial(n) for x in accum];polys[n]=poly
        character_rows.append(dict(n=n,tournament_types=tour//factorial(n),invariant_exterior_polynomial=poly,superdimension=sum((-1)**k*x for k,x in enumerate(poly))))
    a=[r['tournament_types'] for r in character_rows]
    need(a==[1,1,1,2,4,12,56,456],'A000568 control')
    atlas_rows=[];B=[0]*7;F=[0]*7;all_graphs={}
    for n in range(7):
        atlas=graph_atlas(n);all_graphs[n]=atlas
        ext=Counter(x['edges'] for x in atlas if x['exterior'])
        need([ext[k] for k in range(comb(n,2)+1)]==polys[n],('graph-vs-character',n))
        need(sum(x['royle'] for x in atlas)==a[n],('Royle-vs-tournament',n))
        B[n]=sum(x['exterior'] and x['connected'] and x['edges']%2==0 for x in atlas)
        F[n]=sum(x['exterior'] and x['connected'] and x['edges']%2==1 for x in atlas)
        need(B[n]-F[n]==int(n==1)-int(n==3),('connected balance',n))
        atlas_rows.append(dict(n=n,graph_types=len(atlas),exterior_admissible=sum(x['exterior'] for x in atlas),royle_even=sum(x['royle'] for x in atlas),both=sum(x['exterior'] and x['royle'] for x in atlas),connected_bosonic=B[n],connected_fermionic=F[n]))
    need([r['graph_types'] for r in atlas_rows]==[1,1,2,4,11,34,156],'graph atlas completeness')
    # Direct product by actual connected graph types, including the edge grading.
    prod={(0,0):1};top=6
    for n,atlas in all_graphs.items():
        if n==0:continue
        for row in atlas:
            if not(row['connected'] and row['exterior']):continue
            m=row['edges'];nextprod=Counter()
            for (v,e),coefficient in prod.items():
                maxr=min((top-v)//n,1) if m%2 else (top-v)//n
                for r in range(maxr+1):nextprod[v+r*n,e+r*m]+=coefficient
            prod=dict(nextprod)
    for n in range(top+1):
        need([prod.get((n,m),0) for m in range(comb(n,2)+1)]==polys[n],('bigraded product',n))
    # Exact small witnesses: P3 is Royle-even but exterior-odd, C3 vice versa.
    p3=next(x for x in all_graphs[3] if x['edges']==2)
    c3=next(x for x in all_graphs[3] if x['edges']==3)
    need(p3['royle'] and not p3['exterior'],'P3 sign hostile')
    need(c3['exterior'] and not c3['royle'],'C3 sign hostile')
    pairmask=sum(1<<i for i,e in enumerate(combinations(range(6),2)) if set(e)<=set(range(3)) or set(e)<=set(range(3,6)))
    pairrow=next(x for x in all_graphs[6] if x['mask']==min(transform(pairmask,pair_action(6,p)[0]) for p in permutations(range(6))))
    need(not pairrow['exterior'],'two triangles must be excluded')
    # Recover e_n, the ordinary Euler exponents, independently by recursive division.
    residual=a[:];E=[0]*len(a)
    for n in range(1,len(a)):
        E[n]=residual[n]
        for _ in range(E[n]):
            for j in range(len(a)-1,n-1,-1):residual[j]-=residual[j-n]
    need(E[1:]==[1,0,1,2,8,43,398],'Euler exponent control')
    for n in range(1,7):need(E[n]==B[n]+F[n]-(F[n//2] if n%2==0 else 0),'Euler/fermion recurrence')
    report=dict(status='PASS',universe='all permutations n<=7; all graph isomorphism types n<=6 via complete labelled masks and permutation orbits',permutation_checks=permutation_checks,character_rows=character_rows,graph_rows=atlas_rows,connected_royle_even=E[1:],bosonic=B[1:],fermionic=F[1:],hostiles=dict(P3=p3,C3=c3,two_triangles=pairrow),bigraded_component_product='verified through vertex degree 6')
    print(json.dumps(report,sort_keys=True,indent=2))

if __name__=='__main__':main()
