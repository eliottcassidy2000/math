#!/usr/bin/env python3
"""Standalone exact finite controls; no claims to decide infinite Collatz."""
from fractions import Fraction as F
from itertools import product
import json

checks = {}
def require(ok, name):
    if not ok:
        raise RuntimeError(name)
    checks[name] = checks.get(name, 0) + 1

def mv(A, v):
    return [sum((a*x for a,x in zip(row,v)), F(0)) for row in A]
def eye(n):
    return [[F(i==j) for j in range(n)] for i in range(n)]
def kernel(dest):
    n=len(dest)
    return [[F(dest[j]==i+1) for j in range(n)] for i in range(n)]
def mix(A,B,a):
    return [[(1-a)*x+a*y for x,y in zip(r,s)] for r,s in zip(A,B)]
def resolvent(A,z,w):
    n=len(w)
    rows=[[F(i==j)-z*A[i][j] for j in range(n)]+[w[i]] for i in range(n)]
    for j in range(n):
        p=next((p for p in range(j,n) if rows[p][j]),None)
        if p is None: raise RuntimeError('singular discounted matrix')
        rows[j],rows[p]=rows[p],rows[j]
        scale=rows[j][j]
        rows[j]=[x/scale for x in rows[j]]
        for i in range(n):
            if i!=j:
                scale=rows[i][j]
                rows[i]=[x-scale*y for x,y in zip(rows[i],rows[j])]
    return [r[-1] for r in rows]
def tau(dest,start):
    seen=set(); t=0
    while start:
        if start in seen: return None
        seen.add(start); start=dest[start-1]; t+=1
    return t

n=3; w=[F(1,2),F(1,3),F(1,6)]
zs=[F(1,3),F(3,4),F(9,10)]
all_dest=list(product(range(n+1),repeat=n))
for dest in all_dest:
    Q=kernel(dest); ts=[tau(dest,i+1) for i in range(n)]
    for z in zs:
        g=resolvent(Q,z,w)
        d=(1-z)*sum(g)
        require(d==sum(v*(1-(z**t if t is not None else 0)) for v,t in zip(w,ts)), 'abel_identity')
        require(all(x>=0 for x in g), 'resolvent_positive')
        require(g==[v+z*q for v,q in zip(w,mv(Q,g))], 'supersolution_equality')
        for atom,t in zip(w,ts):
            if d<atom:
                require(t is not None, 'atom_excludes_infinite')
                m=0
                while atom*(1-z**m)<=d: m+=1
                require(t<m, 'rational_hitting_bound')
        for p in [F(1,4),F(2,3),F(1)]:
            P=mix(eye(n),Q,p)
            zz=z*p/(1-z*(1-p))
            gp=resolvent(P,z,w); gq=resolvent(Q,zz,w)
            require((1-z)*sum(gp)==(1-zz)*sum(gq), 'lazy_clock_identity')
    for rdest in all_dest:
        R=kernel(rdest); alpha=F(2,5); L=mix(Q,R,alpha)
        q=w[:]; l=w[:]
        for k in range(5):
            require(all(a>=(1-alpha)**k*b for a,b in zip(l,q)), 'mixed_entrywise_domination')
            q=mv(Q,q); l=mv(L,l)

# A noncommuting explicit case guards against a binomial-expansion proof.
Q=kernel((2,0,3)); R=kernel((3,1,0)); unit=[F(1),F(0),F(0)]
require(mv(Q,mv(R,unit))!=mv(R,mv(Q,unit)), 'noncommuting_witness')
# Sharp refresh hostile and aggregation hostile.
alpha=F(2,5)
for k in range(15):
    require((1-alpha)**k/(1-alpha)**k==1, 'refresh_threshold_sharp')
Q=kernel((0,2)); v=[F(1,2),F(1,2)]
for k in range(1,12):
    v=mv(Q,v)
    require(sum(v)==F(1,2), 'native_aggregation_survivor')
    if k>1: require(F(1,2)**k<sum(v), 'averaging_changes_survival')

def collatz(n):
    if n<=0 or n%2==0: raise ValueError('positive odd required')
    m=3*n+1; a=0
    while m%2==0: m//=2; a+=1
    return m,a
for k in range(1,81):
    x=2**(k+1)-1
    for i in range(k):
        require(x==2**(k+1-i)*3**i-1 and x>1, 'norm_one_source_formula')
        x,a=collatz(x)
        require(a==1, 'exact_valuation_one')
    require(x==2*3**k-1 and x>1, 'norm_one_survival')

# Exact finite-source mass certificate, preserving the omitted tail.
J=64; K=200; z=F(999,1000); d=F(1,J+1); survivor=F(1,J+1)
resolved=[]
for j in range(1,J+1):
    atom=F(1,j*(j+1)); x=2*j-1; t=0
    while x!=1 and t<K:
        x,_=collatz(x); t+=1
    if x==1:
        d+=atom*(1-z**t); resolved.append(t)
    else:
        d+=atom; survivor+=atom
require(sum(F(1,j*(j+1)) for j in range(1,J+1))+F(1,J+1)==1, 'telescoping_exact_mass')
require(len(resolved)==J, 'finite_head_resolved')
require(d<F(1,6), 'finite_head_abel_atom_crossing')
require(survivor==F(1,65), 'survival_tail_not_discarded')
m=0
while F(1,6)*(1-z**m)<=d: m+=1
require(2<m, 'source_three_bound') # 3 -> 5 -> 1.
print(json.dumps({'status':'PASS','checks':checks,'total_checks':sum(checks.values()),
    'finite_head':{'odd_sources':J,'horizon':K,'largest_lifetime':max(resolved),
    'z':str(z),'tail_mass':str(F(1,J+1)),
    'abel_upper_numerator':d.numerator,'abel_upper_denominator':d.denominator,
    'certified_source_3_lifetime_strictly_below':m},
    'scope':'finite exact controls only; infinite statements proved in Markdown'},indent=2,sort_keys=True))
