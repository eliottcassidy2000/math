#!/usr/bin/env python3
"""Authored independent divisor/carrier checks, stdlib only, assertions not used."""
from itertools import combinations, product
from math import isqrt, prod, comb
import json

def ck(b, msg):
    if not b: raise RuntimeError(msg)
def divisors(n):
    lo=[d for d in range(1,isqrt(n)+1) if n%d==0]
    return sorted(set(lo+[n//d for d in lo]))
def factor(n):
    fs=[]; p=2
    while p*p<=n:
        a=0
        while n%p==0: n//=p; a+=1
        if a: fs.append((p,a))
        p+=1
    if n>1: fs.append((n,1))
    return fs

def minimal_orientation(ds,N):
    # Preserve each nonfixed unordered pair orbit, except unavoidable swapped endpoints.
    arcs=set(); seen=set(); defects=0
    for u,v in combinations(ds,2):
        e=frozenset((u,v))
        if e in seen: continue
        su,sv=N//u,N//v; se=frozenset((su,sv))
        arcs.add((u,v)); seen.add(e)
        if se!=e:
            arcs.add((su,sv)); seen.add(se)
    for u,v in arcs:
        if (N//u,N//v) not in arcs: defects+=1
    ck(len(arcs)==comb(len(ds),2),'complete orientation')
    return defects,arcs

cases=0; antichains=0
for N in range(2,2001):
    ds=divisors(N); mid=ds[1:-1]; fs=factor(N)
    tuples=list(product(*(range(a+1) for p,a in fs)))
    weighted=[prod(p**e for (p,a),e in zip(fs,t)) for t in tuples]
    ck(sorted(weighted)==ds,'monomial evaluation')
    ck(len(tuples)-1==len(ds)-1,'augmentation')
    ck(sum(weighted)-N==sum(ds[:-1]),'aliquot evaluation')
    f=int(isqrt(N)**2==N); b=(len(mid)-f)//2
    d,arcs=minimal_orientation(mid,N)
    ck(d==b,'min defect witness')
    ck(sum(N//x==x for x in mid)==f,'fixed point')
    ck(sum((N//u)>(N//v) for u,v in combinations(mid,2))==comb(len(mid),2),'numeric orientation is converse')
    pairs=[(x,N//x) for x in mid if x<N//x]
    ck(1+sum(x+y for x,y in pairs)+(isqrt(N) if f else 0)==sum(ds[:-1]),'pair sum')
    cases+=1
for n in range(1,251):
    sn=sum(divisors(n)[:-1])
    for k in range(2,41):
        sm=sum(divisors(n*k)[:-1])
        ck(sm>=sn+n,'strict divisibility growth')
        ck(sm>=k*sn+1,'scaled divisor injection')
        antichains+=1
# Exhaustive all orientations for small divisor interiors, independently derive true minimum.
exhaustive=0
for N in [2,4,6,12,24,30]:
    ds=divisors(N)[1:-1]; pairs=list(combinations(ds,2)); best=len(pairs)
    for mask in range(1<<len(pairs)):
        arcs={(v,u) if mask>>i&1 else (u,v) for i,(u,v) in enumerate(pairs)}
        d=sum((N//u,N//v) not in arcs for u,v in arcs)
        best=min(best,d); exhaustive+=1
    ck(best==(len(ds)-int(isqrt(N)**2==N))//2,'exhaustive minimum')
N=210; ds=divisors(N); mid=ds[1:-1]; pairs=[(d,N//d) for d in mid if d<N//d]
ck(len(ds)==16 and len(ds[:-1])==15 and len(mid)==14,'counts')
ck(min(pairs,key=lambda xy:xy[1]-xy[0])==(14,15),'closest pair')
ck(210 not in ds[:-1] and 1 in ds[:-1],'boundary hostile')
ck(sum(divisors(4)[:-1])!=sum(divisors(9)[:-1]),'prime relabel hostile')
# Complement pairs are unordered; forcing invariant trace to orient them has seven ties.
ck(sum(x+N//x==y+N//y for x,y in combinations(mid,2))==7,'trace ties')
print(json.dumps({'status':'PASS','divisor_carriers_N_2_to_2000':cases,'divisibility_controls':antichains,'small_orientations':exhaustive,'N210':{'divisors':ds,'complement_pairs':pairs,'pair_sums':[a+b for a,b in pairs],'aliquot':sum(ds[:-1]),'orientation_minimum_defect':7,'numeric_orientation_defect':91,'minimum_orientation_arcs':sorted(minimal_orientation(mid,N)[1])},'F1_hostile':{'s4':3,'s9':4,'same_unlabelled_divisor_chain':True}},sort_keys=True,indent=2))
