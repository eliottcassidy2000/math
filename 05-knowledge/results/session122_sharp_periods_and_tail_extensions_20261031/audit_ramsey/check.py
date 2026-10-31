#!/usr/bin/env python3
"""Independent finite controls for the all-budget mathematical period proof.
No producer import. Exact Fraction arithmetic and explicit failures survive -O.
"""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
import json

def require(ok, message):
    if not ok:
        raise RuntimeError(message)

@lru_cache(None)
def weight(a):
    return sum((Q(3,2)**j.bit_count() for j in range(a)), Q(0))

def solve(matrix, rhs):
    m=[[Q(x) for x in row]+[Q(y)] for row,y in zip(matrix,rhs)]
    for col in range(5):
        pivot=next((i for i in range(col,5) if m[i][col]),None)
        if pivot is None:return None
        m[col],m[pivot]=m[pivot],m[col]
        z=m[col][col];m[col]=[v/z for v in m[col]]
        for i in range(5):
            if i!=col:
                z=m[i][col];m[i]=[a-z*b for a,b in zip(m[i],m[col])]
    return tuple(row[-1] for row in m)

def vertices():
    inequalities=[]
    for i in range(5):
        row=[0]*5;row[i]=-1;inequalities.append((row,0))
    for i in range(5):
        row=[0]*5;row[i]=row[(i+2)%5]=1;inequalities.append((row,1))
    out=set()
    for ids in combinations(range(10),5):
        p=solve([inequalities[i][0] for i in ids],[inequalities[i][1] for i in ids])
        if p is not None and all(sum(x*y for x,y in zip(row,p))<=b for row,b in inequalities):out.add(p)
    return out

@lru_cache(None)
def patterns(budget,step,maximal=False):
    out=[]
    # Straight full Cartesian universe; maximality is optional and audited.
    for p in product(range(1,budget), repeat=5):
        if any(p[i]+p[(i+step)%5]>budget for i in range(5)):continue
        if maximal and any(p[i]+p[(i+step)%5]<budget and p[i]+p[(i-step)%5]<budget for i in range(5)):continue
        out.append(p)
    return tuple(out)

def capacities(amax,bmax,maximal):
    values={};candidates=0
    for a in range(1,amax+1):
        for b in range(1,bmax+1):
            best=1
            best=max([best]+[values[x,b]+values[a-x,b] for x in range(1,a)]+[values[a,y]+values[a,b-y] for y in range(1,b)])
            for aa in patterns(a,1,maximal):
                for bb in patterns(b,2,maximal):
                    best=max(best,sum(values[x,y] for x,y in zip(aa,bb)));candidates+=1
            values[a,b]=best
    return values,candidates

@lru_cache(None)
def sharp_statistics(a,b):
    require(b>0 and b%(1<<(a-1).bit_length())==0,'construction hypothesis')
    if a==1:return (b,1,b)
    h=a//2
    aa=(h+1,h,h+1,h,h) if a%2 else (h,)*5
    children=[sharp_statistics(x,b//2) for x in aa]
    n=sum(z[0] for z in children)
    om=max(children[i][1]+children[(i+1)%5][1] for i in range(5))
    al=max(children[i][2]+children[(i+2)%5][2] for i in range(5))
    return n,om,al

def main():
    recurrence=0
    for a in range(2,513):
        h=a//2
        expected=Q(5,2)*weight(h) if a%2==0 else Q(3,2)*weight(h)+weight(h+1)
        require(weight(a)==expected,'binary recurrence vs digit sum');recurrence+=1
    strictness=0
    for total in range(2,257):
        for x in range(1,total):
            require(weight(x)+weight(total-x)<weight(total),'strict superadditivity');strictness+=1
    vv=vertices()
    indicators={tuple(Q(x) for x in p) for p in product((0,1),repeat=5) if all(p[i]+p[(i+2)%5]<=1 for i in range(5))}
    require(vv==indicators|{(Q(1,2),)*5},'exact polytope vertices')
    patterns_checked=0;equal_patterns=0;indicator_tests=0;levels=[]
    for a in range(2,13):
        counts=0
        for aa in patterns(a,1):
            ww=[weight(x) for x in aa]
            require(sum(ww)<=2*weight(a),'C5 weight upper bound')
            for v in indicators:
                require(sum(c*t for c,t in zip(ww,v))<weight(a),'indicator strictness')
                indicator_tests+=1
            if sum(ww)==2*weight(a):
                counts+=1;equal_patterns+=1
                if a>=3:
                    k=(a-1).bit_length()
                    require(max(aa)>1<<(k-2),'deep-child dyadic witness')
            patterns_checked+=1
        levels.append([a,counts])
    f,nmax=capacities(6,16,True)
    full,nfull=capacities(4,6,False)
    require(all(f[key]==val for key,val in full.items()),'maximal/full finite capacities')
    for (a,b),n in f.items():
        L=1<<(a-1).bit_length()
        require((Q(n)==weight(a)*b)==(b%L==0),'finite sharpness iff lattice')
    require(f[3,2]==7 and weight(3).denominator==1,'integral slope hostile')
    require(f[6,4]<4*weight(6) and weight(6).denominator==1,'higher integral slope hostile')
    construction_count=0
    for a in range(1,65):
        L=1<<(a-1).bit_length()
        for copies in (1,2,3):
            b=copies*L;n,om,al=sharp_statistics(a,b)
            require((n,om,al)==(weight(a)*b,a,b),'actual construction statistics');construction_count+=1
    out={'status':'PASS','scope':'independent finite controls; all-budget result rests on audited proof',
         'weight_recurrence_controls':recurrence,'strict_superadditivity_pairs':strictness,
         'polytope_vertices':len(vv),'indicator_vertices':len(indicators),
         'positive_clique_patterns':patterns_checked,'indicator_strictness_tests':indicator_tests,
         'sharp_clique_patterns':equal_patterns,'sharp_patterns_by_A':levels,
         'capacity_rectangle':[6,16],'maximal_candidates':nmax,
         'full_capacity_rectangle':[4,6],'full_candidates':nfull,
         'capacities':[[f[a,b] for b in range(1,17)] for a in range(1,7)],
         'sharp_construction_statistics':construction_count,
         'hostiles':{'F(3,2)':f[3,2],'F(6,4)':f[6,4]},
         'literal_graph_census_claimed':False}
    print(json.dumps(out,sort_keys=True,indent=2))

if __name__=='__main__':main()
