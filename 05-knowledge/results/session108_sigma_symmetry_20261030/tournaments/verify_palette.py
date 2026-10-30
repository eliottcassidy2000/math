#!/usr/bin/env python3
"""Exact independent orbit and operator checks; standard library + sympy."""
from itertools import product, permutations
from collections import Counter
from fractions import Fraction
from math import gcd, factorial
import json
import sympy as sp

def require(x,msg):
    if not x: raise RuntimeError(msg)
def lin(cols,x):
    y=0
    for i in range(3):
        if x>>i&1:y^=cols[i]
    return y
GL=sorted(set(tuple(lin(c,x) for x in range(8)) for c in product(range(1,8),repeat=3) if len(set(lin(c,x) for x in range(8)))==8))
AGL=sorted(set(tuple(v^b for v in g) for g in GL for b in range(8)))
# Named nonzero port is the functional z -> z_0. Preserve all its answers.
H=[g for g in GL if all((g[x]&1)==(x&1) for x in range(8))]
groups={'named':[tuple(range(8))],'GL':GL,'port_stabilizer':H,'AGL':AGL}
require([len(GL),len(H),len(AGL)]==[168,24,1344],'group sizes')
def parts(n,lo=1):
    if n==0:yield ()
    for j in range(lo,n+1,2):
        for t in parts(n-j,j):yield (j,)+t

def fixed(g,k):
    return sum(iterate(g,x,k)==x for x in range(8))
def iterate(g,x,k):
    for _ in range(k):x=g[x]
    return x
def burn(n,G):
    ans=Fraction(0)
    for lam in parts(n):
        e=sum((x-1)//2 for x in lam)+sum(gcd(lam[i],lam[j]) for i in range(len(lam)) for j in range(i))
        z=1
        for m,c in Counter(lam).items():z*=m**c*factorial(c)
        w=0
        for g in G:
            f=1
            for l in lam:f*=fixed(g,l)
            w+=f
        ans+=Fraction(2**e*w,z*len(G))
    require(ans.denominator==1,'Burnside integral')
    return int(ans)

def pairs(n):return [(i,j) for i in range(n) for j in range(i+1,n)]
def transform(n,mask,colors,p,g):
    edges=pairs(n);out=0;d={e:i for i,e in enumerate(edges)}
    for k,(i,j) in enumerate(edges):
        u,v=p[i],p[j];bit=(mask>>k)&1
        if u>v:u,v=v,u;bit^=1
        out|=bit<<d[u,v]
    cs=[0]*n
    for i in range(n):cs[p[i]]=g[colors[i]]
    return out,tuple(cs)

def banks(G,N=3):
    result=[];lookups=[]
    for n in range(N+1):
        seen={};reps=[]
        for mask in range(2**len(pairs(n))):
            for cs in product(range(8),repeat=n):
                key=mask,cs
                if key in seen:continue
                idx=len(reps);reps.append(key)
                for p in permutations(range(n)):
                    for g in G:seen[transform(n,mask,cs,p,g)]=idx
        require(len(reps)==burn(n,G),'direct orbit count')
        result.append(reps);lookups.append(seen)
    return result,lookups

def delete(n,state,v):
    mask,cs=state;remaining=[i for i in range(n) if i!=v];d={e:k for k,e in enumerate(pairs(n))}
    out=sum(((mask>>d[remaining[i],remaining[j]])&1)<<k for k,(i,j) in enumerate(pairs(n-1)))
    return out,tuple(cs[i] for i in remaining)
def attach(n,state,c,bits):
    mask,cs=state;d={e:k for k,e in enumerate(pairs(n))}
    out=0
    for k,(i,j) in enumerate(pairs(n+1)):
        bit=((bits>>i)&1) if j==n else ((mask>>d[i,j])&1)
        out|=bit<<k
    return out,cs+(c,)
def source(n,state):
    mask,cs=state;d={e:k for k,e in enumerate(pairs(n))}
    for v in range(n):
        if all((((mask>>d[min(v,w),max(v,w)])&1)==(v<w)) for w in range(n) if v!=w):return v
    return None

def cycle_type(g):
    unseen=set(range(8));ls=[]
    while unseen:
        x=min(unseen);y=g[x];l=1;unseen.remove(x)
        while y!=x:unseen.remove(y);l+=1;y=g[y]
        ls.append(l)
    return tuple(sorted(ls))
output={}
for name,G in groups.items():
    counts=[burn(n,G) for n in range(9)]
    item={'group_order':len(G),'counts_n0_to_n8':counts,'harmonic_dimensions':[counts[0]]+[counts[n]-counts[n-1] for n in range(1,9)],'palette_cycle_types':{','.join(map(str,k)):v for k,v in sorted(Counter(cycle_type(g) for g in G).items())}}
    # Full operators for the symmetry quotients; named counts inherit prior q-color checks.
    if name!='named':
        B,L=banks(G);D=[None];U=[]
        for n in range(1,4):
            mat=sp.zeros(len(B[n-1]),len(B[n]))
            for i,s in enumerate(B[n]):
                for v in range(n):mat[L[n-1][delete(n,s,v)],i]+=1
            D.append(mat)
        for n in range(3):
            mat=sp.zeros(len(B[n+1]),len(B[n]))
            for i,s in enumerate(B[n]):
                for c in range(8):
                    for bits in range(2**n):mat[L[n+1][attach(n,s,c,bits)],i]+=1
            U.append(mat)
            rhs=8*2**n*sp.eye(len(B[n]))
            if n:rhs+=2*U[n-1]*D[n]
            require(D[n+1]*U[n]==rhs,'commutator')
        for n in range(1,4):
            K=U[n-1]*D[n]/(n*8*2**(n-1))
            expected={sp.Rational(n-j,n):counts[j]-(counts[j-1] if j else 0) for j in range(n+1)}
            expected={k:v for k,v in expected.items() if v}
            require(K.eigenvals()==expected,'spectrum')
            weights=Counter(L[n].values());pi=sp.Matrix([sp.Rational(weights[i],len(L[n])) for i in range(len(B[n]))])
            require(K*pi==pi,'stationary')
            require(K*sp.diag(*pi)==sp.diag(*pi)*K.T,'detailed balance')
        reps=[];seen=set()
        for c in range(8):
            if c not in seen:reps.append(c);seen.update(g[c] for g in G)
        source_counts=[sum(source(n,s) is not None for s in B[n]) for n in range(1,4)]
        predicted=[sum(burn(n-1,[g for g in G if g[c]==c]) for c in reps) for n in range(1,4)]
        require(source_counts==predicted,'stabilizer source formula')
        item.update(direct_orbit_counts=[len(b) for b in B],source_bearing_n1_to_n3=source_counts,checks='orbit counts n<=3; commutators input n<=2; spectra, stationarity, detailed balance n<=3; stabilizer source formula n<=3')
    output[name]=item
print(json.dumps(output,indent=2,sort_keys=True))
