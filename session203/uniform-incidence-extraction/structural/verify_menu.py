#!/usr/bin/env python3
"""Exact pinned-path and finite-colored-rectangle checks for the menu addendum."""
from itertools import product, combinations
import json
import random
from verify import hp, require

def graph(bits,pa,pc,arc):
    m=len(bits);rows=[0]*(m+2)
    if arc:rows[0]|=2
    else:rows[1]|=1
    for j,b in enumerate(bits,2):
        for a,p in [(0,pa),(1,pc)]:
            if b^p:rows[a]|=1<<j
            else:rows[j]|=1<<a
        for i in range(2,j):rows[i]|=1<<j
    return rows

def pinned_paths(free,pa,pc,arc):
    bits=[0,0,1,1]+list(free)+[0,0,1,1];m=len(bits)
    prefix=list(range(4));suffix=list(range(m-4,m))
    pM=next(i for i in prefix if bits[i]^pa)
    pR=next(i for i in prefix if i!=pM and bits[i]^pc)
    sL=next(i for i in suffix if not bits[i]^pa)
    sM=next(i for i in suffix if i!=sL and not bits[i]^pc)
    rows=graph(bits,pa,pc,arc);seen=set()
    for assign in product(range(3),repeat=len(free)):
        bins=[[i for i in prefix if i not in [pM,pR]],[pM],[pR]]
        for i,c in enumerate(assign,4):bins[c].append(i)
        bins[0].append(sL);bins[1].append(sM)
        bins[2]+=[i for i in suffix if i not in [sL,sM]]
        bins=[sorted(b) for b in bins]
        require(bins[0][-1]==sL and bins[1][0]==pM and bins[1][-1]==sM and bins[2][0]==pR,'prescribed extrema')
        path=tuple([i+2 for i in bins[0]]+[0]+[i+2 for i in bins[1]]+[1]+[i+2 for i in bins[2]])
        require(len(path)==len(set(path))==m+2,'permutation')
        require(all(rows[i]>>j&1 for i,j in zip(path,path[1:])),'pinned path arcs')
        seen.add(path)
    require(len(seen)==3**len(free),'three-slot count')
    if len(free)<=3:require(hp(rows)>=len(seen),'independent subset HP count')
    return len(seen)

def compositions(n,k):
    if k==1:
        yield(n,);return
    for a in range(n+1):
        for tail in compositions(n-a,k-1):yield(a,)+tail

def rectangle(colors,q):
    old={}
    for c in range(len(colors[0])):
        occurrences=0
        for a,b in combinations(range(2*q),2):
            if colors[a][c]!=colors[b][c]:continue
            occurrences+=1;key=(a,b,colors[a][c])
            if key in old:
                d=old[key]
                require(d<c,'different columns')
                require(colors[a][d]==colors[b][d]==colors[a][c]==colors[b][c],'monochromatic rectangle')
                return a,b,d,c
            old[key]=c
        require(occurrences>=q,'column pair lower bound')
    raise RuntimeError('rectangle absent')

def main():
    words=paths=dpchecks=0
    for t in range(6):
        for free in product(range(2),repeat=t):
            for pa,pc,arc in product(range(2),repeat=3):
                paths+=pinned_paths(free,pa,pc,arc);words+=1
                if t<=3:dpchecks+=1
    compchecks=0
    for q in range(1,6):
        for sizes in compositions(2*q,q):
            require(sum(n*(n-1)//2 for n in sizes)>=q,'pair convexity')
            compchecks+=1
    rng=random.Random(203);rectchecks=0
    for q in range(1,9):
        K=q*(2*q-1)
        for _ in range(100):
            colors=[[rng.randrange(q) for c in range(K+1)] for r in range(2*q)]
            a,b,c,d=rectangle(colors,q)
            # Independently test the negative quadratic rectangle identity.
            # Bit increment groups have sizes b-a and d-c respectively.
            coeff=[[-rng.randrange(1,20) for j in range(K)] for i in range(2*q-1)]
            def f(r,s):return sum(coeff[i][j] for i in range(r) for j in range(s))+7*r-11*s+29
            mixed=f(b,d)-f(b,c)-f(a,d)+f(a,c)
            require(mixed==sum(coeff[i][j] for i in range(a,b) for j in range(c,d))<0,'rectangle quadratic identity')
            rectchecks+=1
    growth=3**20;target=6*2**29+1
    require(growth==3486784401 and target==3221225473 and growth>target,'m28 bound')
    require(8+(3+1)*(2*3-1)==28,'three-template bound')
    print(json.dumps({'status':'PASS','pinned_word_polarity_arc_cases':words,'explicit_distinct_pinned_paths':paths,'independent_subset_DP_cases':dpchecks,'all_compositions_checked_q_at_most_5':compchecks,'deterministic_random_rectangle_controls':rectchecks,'m28_lower_bound':growth,'m28_larger_target_upper_bound':target,'menu_bound_at_m28':3,'scope':'pinned lower bound and quadratic finite-menu rectangle mechanism'},sort_keys=True,indent=2))
if __name__=='__main__':main()
