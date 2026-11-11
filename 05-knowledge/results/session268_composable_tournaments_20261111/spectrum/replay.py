#!/usr/bin/env python3
"""Exact cyclic tomography, inverse compilation, and physical seam controls."""
from pathlib import Path
from math import factorial
from itertools import combinations, permutations
from functools import lru_cache
import json, random, heapq
OUT=Path(__file__).resolve().parent

def transitive(n): return tuple(sum(1<<j for j in range(i+1,n)) for i in range(n))
def mask_graph(n,m):
    rows=[0]*n
    for k,(i,j) in enumerate(combinations(range(n),2)):
        a,b=(i,j) if (m>>k)&1 else (j,i)
        rows[a]|=1<<b
    return tuple(rows)
def cyclic(A,B,C):
    parts=[A,B,C]; offsets=[0,len(A),len(A)+len(B)]; rows=[]
    for k,P in enumerate(parts):
        target=(k+1)%3; cross=((1<<len(parts[target]))-1)<<offsets[target]
        rows.extend((r<<offsets[k])|cross for r in P)
    return tuple(rows)
def path_data(rows):
    n=len(rows); end=[[0]*n for _ in range(1<<n)]; hs=[0]*(1<<n); hs[0]=1
    for v in range(n): end[1<<v][v]=1
    for mask in range(1,1<<n):
        hs[mask]=sum(end[mask])
        for v,x in enumerate(end[mask]):
            if not x: continue
            avail=rows[v]&~mask
            while avail:
                bit=avail&-avail; avail-=bit
                end[mask|bit][bit.bit_length()-1]+=x
    return hs,end

def count(rows): return path_data(rows)[0][-1]
def cover_vector(rows):
    n=len(rows); hs,_=path_data(rows); old=[0]*(1<<n);old[0]=1;ans=[]
    for k in range(1,n+1):
        new=[0]*(1<<n)
        for mask in range(1,1<<n):
            if mask.bit_count()<k:continue
            sub=mask
            while sub:
                new[mask]+=hs[sub]*old[mask^sub];sub=(sub-1)&mask
        ans.append(new[-1]);old=new
    return tuple(ans)
@lru_cache(None)
def F(d,r):
    if r<1 or r>d:return 0
    if d==1:return 1
    return r*(F(d-1,r)+F(d-1,r-1)) if r>1 else 1

def coeff(d,r):
    a,b,c=F(d,r-1),F(d,r),F(d,r+1)
    return 3*b*b+2*b*c+c*c+a*a+2*a*b

def probe(q,d):return sum(coeff(d,r)*q[r-1] for r in range(1,min(len(q),d+1)+1))
def decode(h,probes):
    q=[h]
    for d,p in enumerate(probes,1):
        v=p-sum(coeff(d,r)*q[r-1] for r in range(1,d+1));den=factorial(d)**2
        if v%den:raise ValueError('not a physical tomography vector')
        q.append(v//den)
    return tuple(q)
def atom(b,c):return 2+((1<<b)-1)*((1<<c)-1)
def inverse_atom(h):
    if h<3:return None
    M=h-2
    for b in range(1,(M+1).bit_length()):
        f=(1<<b)-1
        if M%f:continue
        q=M//f+1
        if q&(q-1)==0:
            c=q.bit_length()-1
            if b<=c:return b,c
    return None

def inverse_lift(N,q):
    h=q[0];Q=q[1] if len(q)>1 else 0; num=N-h+Q; den=h+Q
    if num<2*den or num%den:return None
    ratio=num//den
    return ratio.bit_length()-1 if ratio&(ratio-1)==0 else None

def seam(A,B,a=0,b=0):
    na=len(A);nb=len(B);rows=[r|(((1<<nb)-1)<<na) for r in A]
    rows += [r<<na for r in B]
    rows[a]&=~(1<<(na+b));rows[na+b]|=1<<a
    return tuple(rows)

def monoid_census(bound):
    atoms={}
    for b in range(1,bound.bit_length()+1):
        for c in range(b,bound.bit_length()+1):
            h=atom(b,c)
            if h<=bound:
                if h in atoms and atoms[h]!=(b,c):raise AssertionError('atom collision')
                atoms[h]=(b,c)
    gens=sorted(atoms);seen={1};heap=[1]
    while heap:
        x=heapq.heappop(heap)
        for g in gens:
            y=x*g
            if y>bound:break
            if y not in seen:seen.add(y);heapq.heappush(heap,y)
    return {'bound':bound,'atoms':len(atoms),'product_counts':len(seen),
            'odd_counts_in_interval':(bound+1)//2,'first_absent_odds':[x for x in range(1,100,2) if x not in seen],
            'max_atom_vertices':max(b+c+1 for b,c in atoms.values())}

def main():
    rng=random.Random(268);seeds=[]
    for n in range(1,5):
        seeds += [mask_graph(n,m) for m in range(1<<(n*(n-1)//2))]
    for n in range(5,8):
        seeds += [mask_graph(n,rng.getrandbits(n*(n-1)//2)) for _ in range(10)]
    probes=0; affine=0
    for T in seeds:
        q=cover_vector(T);n=len(T)
        if q[-1]!=factorial(n):raise AssertionError('cover normalization')
        ps=[]
        for d in range(1,n):
            value=probe(q,d);ps.append(value)
            if n+2*d<=10:
                if count(cyclic(transitive(d),transitive(d),T))!=value:raise AssertionError('physical probe')
                probes+=1
        if decode(q[0],ps)!=q:raise AssertionError('tomography')
        for b in range(1,4):
            h=q[0];Q=q[1] if n>1 else 0;target=(h+Q)*(1<<b)+h-Q
            if count(cyclic((0,),transitive(b),T))!=target:raise AssertionError('affine lift')
            if inverse_lift(target,q)!=b:raise AssertionError('lift inverse')
            affine+=1
    atom_controls=0;seam_controls=0
    for b in range(1,5):
        for c in range(b,6):
            A=cyclic((0,),transitive(b),transitive(c));h=atom(b,c)
            if count(A)!=h or inverse_atom(h)!=(b,c):raise AssertionError('atom')
            atom_controls+=1
            if len(A)<=7:
                for m,a0,r0 in [(8,24,15),(10,20,11)]:
                    value=a0*h+r0*((1<<b)-2)-1
                    if count(seam(A,mask_graph(5,m)))!=value:raise AssertionError('seam ray')
                    seam_controls+=1
    qa=(9,150,756,1704,1800,720);qb=(9,126,684,1656,1800,720)
    out={'scope':'physical exhaustive labelled seeds through n4 plus ten seeded random seeds at each n5,n6,n7',
         'seed_count':len(seeds),'physical_tomography_checks':probes,'physical_affine_lift_checks':affine,
         'atom_controls':atom_controls,'physical_seam_ray_checks':seam_controls,
         'coefficient_triangle':[[coeff(d,r) for r in range(1,d+2)] for d in range(1,7)],
         'hidden_seed_probe_pairs':[[probe(qa,d),probe(qb,d)] for d in range(1,6)],
         'atom_monoid_censuses':[monoid_census(x) for x in [1000,10000,100000,1000000]]}
    OUT.joinpath('results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
