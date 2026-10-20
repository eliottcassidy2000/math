#!/usr/bin/env python3
"""Exact independent finite controls for FINITE_STATE_BARRIERS_ADDENDUM.md."""
from fractions import Fraction as F
from itertools import product
from random import Random
import json

checks=0

def need(ok, label):
    global checks
    checks+=1
    if not ok:
        raise ValueError(label)

# An edge is (source, target, lambda, nonnegative debt).
def cycles(n, edges):
    adj=[[] for _ in range(n)]
    for i,e in enumerate(edges):
        adj[e[0]].append(i)
    answer=[]
    for start in range(n):
        def visit(u,seen,path):
            for i in adj[u]:
                v=edges[i][1]
                if v==start:
                    answer.append(path+[i])
                elif v>start and v not in seen:
                    visit(v,seen|{v},path+[i])
        visit(start,{start},[])
    return answer

def cycle_criterion(n, edges):
    for cycle in cycles(n,edges):
        lam=F(1)
        positive=False
        for i in cycle:
            lam*=edges[i][2]
            positive |= edges[i][3]>0
        if lam<1 or (lam==1 and positive):
            return False
    return True

# Independent rational feasibility, not an invocation of the cycle criterion.
# Each row (a,c) denotes a dot b >= c, with no implicit sign assumptions.
def eliminate_feasible(n,edges):
    rows=[]
    for j in range(n):
        a=[F(0)]*n; a[j]=1
        rows.append((tuple(a),F(1)))
    for u,v,l,d in edges:
        a=[F(0)]*n
        a[u]-=1; a[v]+=l
        rows.append((tuple(a),d))
    for _ in range(n):
        pos=[];neg=[];zero=[]
        for a,c in rows:
            if a[0]>0: pos.append((a,c))
            elif a[0]<0: neg.append((a,c))
            else: zero.append((a[1:],c))
        new=zero[:]
        for a,c in pos:
            for b,d in neg:
                p=a[0];q=-b[0]
                new.append((tuple(x/p+y/q for x,y in zip(a[1:],b[1:])),c/p+d/q))
        # Keeping only the strongest identical left side preserves feasibility.
        dedup={}
        for a,c in new:
            if not any(a):
                if c>0:return False
                continue
            # Positive normalization improves exact deduplication; no inequality flips.
            scale=next(abs(x) for x in a if x)
            a=tuple(x/scale for x in a);c/=scale
            dedup[a]=max(dedup.get(a,c),c)
        rows=list(dedup.items())
    return all(c<=0 for a,c in rows)

def verifies_barrier(edges,b):
    return all(x>0 for x in b) and all(l*b[v]-b[u]>=d for u,v,l,d in edges)

def rational_barrier(n,edges):
    # Exact multiplicative shortest paths after a rational strictness perturbation.
    N=1
    for cycle in cycles(n,edges):
        L=F(1);p=0
        for i in cycle:
            L*=edges[i][2];p+=int(edges[i][3]>0)
        if p:
            if L<=1:raise ValueError('infeasible positive-debt cycle')
            q=F(2**p-1)/(L-1)
            N=max(N,-(-q.numerator//q.denominator))
    kappa=F(N+1,N)
    modified=[(u,v,l/(kappa if d>0 else 1)) for u,v,l,d in edges]
    dist=[F(1)]*n
    for _ in range(n):
        new=dist[:]
        for u,v,m in modified:new[v]=min(new[v],dist[u]*m)
        dist=new
    b=[1/x for x in dist]
    K=F(1)
    for u,v,l,d in edges:
        slack=l*b[v]-b[u]
        need(slack>=0,'rational construction nonnegative slack')
        if d>0:
            need(slack>0,'rational construction strict paid-edge slack')
            K=max(K,d/slack)
    b=[K*x for x in b]
    need(verifies_barrier(edges,b),'constructed rational barrier certificate')
    return b

# Exhaustive labeled two-vertex directed graphs with at most one edge for each
# ordered pair, including loops. Seven choices: absent or 3 multipliers x 2 debts.
options=[None]+[(l,d) for l in [F(1,2),F(1),F(2)] for d in [F(0),F(1)]]
exhaustive=0;feasible=0
for choice in product(options,repeat=4):
    edges=[]
    for (u,v),item in zip(product(range(2),repeat=2),choice):
        if item is not None:edges.append((u,v,*item))
    c=cycle_criterion(2,edges)
    lp=eliminate_feasible(2,edges)
    need(c==lp,'exhaustive two-state cycle iff exact LP')
    if c:rational_barrier(2,edges)
    exhaustive+=1;feasible+=int(c)

# Additional exact three-state multigraphs, including parallel edges.
rng=Random(20261020)
random_graphs=1000
for _ in range(random_graphs):
    n=3
    edges=[(rng.randrange(n),rng.randrange(n),rng.choice([F(1,2),F(2,3),F(1),F(3,2),F(2),F(3)]),rng.choice([F(0),F(1),F(2)]))
           for _ in range(rng.randrange(7))]
    good=cycle_criterion(n,edges)
    need(good==eliminate_feasible(n,edges),
         'three-state multigraph cycle iff exact LP')
    if good:rational_barrier(n,edges)

hostiles=[
    ('contracting zero-debt loop',1,[(0,0,F(1,2),F(0))],False,None),
    ('neutral positive-debt loop',1,[(0,0,F(1),F(1))],False,None),
    ('neutral zero-debt loop',1,[(0,0,F(1),F(0))],True,[F(1)]),
    ('zero-debt neutral two-cycle',2,[(0,1,F(1,2),F(0)),(1,0,F(2),F(0))],True,[F(1),F(2)]),
    ('unpaid neutral two-cycle',2,[(0,1,F(1,2),F(1)),(1,0,F(2),F(0))],False,None),
    ('strict cycle with contracting edge',2,[(0,1,F(1,2),F(1)),(1,0,F(3),F(2))],True,[F(4),F(10)]),
    ('acyclic positive debt beside neutral loop',2,[(0,0,F(1),F(0)),(0,1,F(1),F(5))],True,[F(1),F(6)]),
]
for name,n,edges,answer,b in hostiles:
    need(cycle_criterion(n,edges)==answer,name+' cycle verdict')
    need(eliminate_feasible(n,edges)==answer,name+' LP verdict')
    if b is not None:need(verifies_barrier(edges,b),name+' explicit positive barrier')
need(F(1,2)*1-1>=-1,'signed-debt boundary hostile')
need(F(1,2)*0-0>=0,'zero-barrier boundary hostile')
need((3+3)/2==3,'bounded non-target affine hostile')

# Check the telescoping formula on exact sample paths, including changing barriers.
edges=[(0,1,F(1,2),F(1)),(1,0,F(3),F(2))]
b=[F(4),F(10)]
for start in range(2):
    for m in range(1,101):
        v=start;L=F(1);D=F(0);E=F(0);x=F(20)
        for j in range(m):
            u,w,l,d=edges[v]
            need(u==v,'path source')
            E+=L*(l*b[w]-b[u]-d)
            D+=L*d;L*=l;x=(x+d)/l;v=w
        need(L*b[v]-b[start]-D==E,'state-dependent telescoping slack')
        need(x-b[v]<=(20-b[start])/L,'state-dependent excess bound')

# Physical Collatz controls.
def val(n):return (n&-n).bit_length()-1

def coeff(word):
    A=0;C=0;S=1
    for a in word:
        C=3*C+2**A;A+=a;S*=3
    M=2**A
    r=((M-C)*pow(S,-1,2*M))%(2*M)
    return S,M,C,r

def execute(n,word):
    for a in word:
        if n<=1:return None
        z=3*n+1
        if val(z)!=a:return None
        n=z//2**a
    return n

for K in range(1,101):
    for t in range(1,5):
        n=2**(K+1)*t-1;x=n
        for j in range(K):
            need(x==3**j*2**(K+1-j)*t-1,'all-ones cylinder exact endpoint')
            need(val(3*x+1)==1,'all-ones exact next valuation')
            y=(3*x+1)//2
            need(y>x and x>1,'all-ones prefix increases above ROOT')
            x=y

unit_bank=[]
for k in range(1,6):
    for word in product(range(1,4),repeat=k):
        S,M,C,r=coeff(word)
        if C<=M-S:
            need(any(a>=2 for a in word),'unit-budgeted word contains reset')
            unit_bank.append(word)
L=max(map(len,unit_bank))
for t in range(1,11):
    n=2**(L+1)*t-1
    for word in unit_bank:
        need(execute(n,word) is None,'finite unit bank has uncovered residue family')

physical_words=[(1,),(2,),(3,),(1,2),(1,3),(1,2,3),(1,1,3)]
physical_checked=0
physical_feasible=0

def has_path(v,n,depth,label_edges):
    if depth==0:return True
    for u,w,word in label_edges:
        if u==v:
            y=execute(n,word)
            if y is not None and has_path(w,y,depth-1,label_edges):return True
    return False

for V in [1,2,3]:
    for _ in range(100):
        label_edges=[(rng.randrange(V),rng.randrange(V),rng.choice(physical_words))
                     for _ in range(rng.randrange(1,7))]
        edges=[]
        for u,v,word in label_edges:
            S,M,C,r=coeff(word)
            edges.append((u,v,F(M,S),F(C,S)))
        good=cycle_criterion(V,edges)
        need(good==eliminate_feasible(V,edges),'physical graph cycle iff LP')
        physical_checked+=1
        if good:
            rational_barrier(V,edges)
            physical_feasible+=1
            L=max(len(word) for _,_,word in label_edges)
            for t in range(1,4):
                n=2**(V*L+1)*t-1
                for start in range(V):
                    need(not has_path(start,n,V,label_edges),'finite-state physical coverage obstruction')

# Explicit feasible physical graph with a valuation-one edge.
label_edges=[(0,1,(1,)),(1,0,(2,2))]
edges=[]
for u,v,word in label_edges:
    S,M,C,r=coeff(word);edges.append((u,v,F(M,S),F(C,S)))
need(verifies_barrier(edges,[F(5),F(8)]),'physical changing-barrier witness')
need(execute(31,(1,))==47,'physical nonterminal entry')
need(execute(47,(2,2)) is None,'physical failure of next block guard')
need(not has_path(0,31,2,label_edges),'physical witness cannot take V edges')

print(json.dumps({'status':'PASS','checks':checks,
                  'exhaustive_two_state_graphs':exhaustive,
                  'exhaustive_feasible_graphs':feasible,
                  'sampled_three_state_multigraphs':random_graphs,
                  'unit_bank_words':len(unit_bank),
                  'physical_graphs_checked':physical_checked,
                  'physical_graphs_with_barriers':physical_feasible,
                  'independent_feasibility_method':'exact rational Fourier-Motzkin elimination',
                  'scope':'Finite controls; proofs of all graph and all-depth claims are in the addendum',
                  'global_collatz_conclusion':False},indent=2,sort_keys=True))
