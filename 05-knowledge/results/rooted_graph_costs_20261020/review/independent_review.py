#!/usr/bin/env python3
"""Independent exact tests of star tomography and integrated finite patching."""
from fractions import Fraction as F
from itertools import product, combinations
from math import comb
from random import Random
import json

def require(x,label):
    if not x:raise RuntimeError(label)

def solve3(rows):
    A=[[F(z) for z in row] for row in rows]
    for i in range(3):
        pivot=next((j for j in range(i,3) if A[j][i]),None)
        if pivot is None:return None
        A[i],A[pivot]=A[pivot],A[i]
        d=A[i][i];A[i]=[x/d for x in A[i]]
        for j in range(3):
            if j!=i:
                d=A[j][i];A[j]=[x-d*y for x,y in zip(A[j],A[i])]
    return tuple(row[-1] for row in A)

def exact_conditional_optimum(costs,sizes,C):
    # q0,q1 are the independent conditional probabilities of action 0 at the
    # two root words, whose probabilities are 1/3 and 2/3. Minimize t.
    planes=[(1,0,0,0),(1,0,0,1),(0,1,0,0),(0,1,0,1)]
    lines=[]
    for f,s in zip(costs,sizes):
        a=F(f[0][0]-f[0][1],3)
        b=F(2*(f[1][0]-f[1][1]),3)
        c=F(f[0][1]+2*f[1][1],3)-C*s
        lines.append((a,b,c));planes.append((-a,-b,1,c))
    feasible=[]
    for rows in combinations(planes,3):
        sol=solve3(rows)
        if sol is None:continue
        q0,q1,t=sol
        if 0<=q0<=1 and 0<=q1<=1 and all(a*q0+b*q1+c<=t for a,b,c in lines):feasible.append(sol)
    require(feasible,'LP has a feasible minimizing vertex')
    return min(feasible,key=lambda z:z[2])

def hypothesis(costs,sizes,C,k):
    r=len(sizes)
    for m in product(range(k+1),repeat=r):
        if sum(m)>k+r-1:continue
        a=min(sum(m[i]*costs[i][0][b] for i in range(r)) for b in range(2))
        b=min(sum(m[i]*costs[i][1][b] for i in range(r)) for b in range(2))
        if F(a+2*b,3)>C*(1+sum(x*y for x,y in zip(m,sizes))):return False
    return True

rng=Random(20261020)
passed=failed_assumption=0;strict=0;examples=[]
for trial in range(600):
    costs=[[[rng.randrange(4) for b in range(2)] for root in range(2)] for i in range(3)]
    sizes=[rng.randrange(3) for i in range(3)]
    C=rng.choice([F(0),F(1,3),F(1,2),F(1),F(2)])
    k=1+trial%5
    if not hypothesis(costs,sizes,C,k):failed_assumption+=1;continue
    q0,q1,opt=exact_conditional_optimum(costs,sizes,C)
    eps=C*F(1+sum(sizes)-min(sizes),k)
    require(opt<=eps,f'integrated finite patching fails {costs,sizes,C,k,opt,eps}')
    passed+=1
    if opt>0:
        strict+=1
        if len(examples)<3:examples.append(dict(costs=costs,sizes=sizes,C=str(C),k=k,optimum_excess=str(opt),tariff=str(eps),conditional_probabilities=[str(q0),str(q1)]))

# Star top coefficient by a third computation: discrete differences of |sum|.
star=[]
for n in range(2,21,2):
    vals=[F(n-abs(n-2*k),2) for k in range(n+1)]
    while len(vals)>1:vals=[b-a for a,b in zip(vals,vals[1:])]
    top=vals[0]/(1 << n)
    desired=F((-1)**(n//2)*comb(n-2,n//2-1),1 << (n-1))
    require(top==desired,'star coefficient via finite differences')
    star.append([n,str(top)])

# The five-path: enumerate all boundary subsets of endpoints and all possible
# interior subsets and edges. An omitted endpoint is absent, not a free vertex.
V=range(5);edges={(0,1),(1,2),(2,3),(3,4)};count=0
for rootmask in range(4):
    roots={v for j,v in enumerate([0,4]) if (rootmask>>j)&1}
    for imask in range(8):
        interior={v for j,v in enumerate([1,2,3]) if (imask>>j)&1}
        allowed=sorted(e for e in edges if set(e)<=roots|interior and set(e)&interior)
        for emask in range(1 << len(allowed)):
            used=[e for j,e in enumerate(allowed) if (emask>>j)&1]
            best=100
            for bits in product([0,1],repeat=len(interior)):
                x={0:0,4:1};x.update(zip(sorted(interior),bits))
                best=min(best,sum(x[a]==x[b] for a,b in used))
            require(F(best)<=F(len(interior),3),'five path all root subsets')
            count+=1

result=dict(status='PASS',integrated_patching_random_seed=20261020,integrated_patching_cases=600,assumption_passes=passed,assumption_failures=failed_assumption,positive_optimum_excess_cases=strict,positive_excess_examples=examples,star_top_coefficients=star,five_path_all_root_subsets_and_interior_patterns=count)
print(json.dumps(result,indent=2))
