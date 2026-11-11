#!/usr/bin/env python3
"""Exact controls: directed odd-path capacities versus independent subset DP."""
from functools import lru_cache
from itertools import combinations, permutations
from collections import Counter
import json
from pathlib import Path

def require(ok,msg):
    if not ok: raise RuntimeError(msg)

def rows(n,code):
    a=[0]*n
    for k,(i,j) in enumerate(combinations(range(n),2)):
        u,v=(i,j) if code>>k&1 else (j,i)
        a[u]|=1<<v
    return tuple(a)

def counter(a):
    @lru_cache(None)
    def end(mask,v):
        rest=mask^(1<<v)
        if not rest:return 1
        return sum(end(rest,u) for u in range(len(a)) if rest>>u&1 and a[u]>>v&1)
    @lru_cache(None)
    def h(mask):
        return sum(end(mask,v) for v in range(len(a)) if mask>>v&1) if mask else 1
    return h

def capacities(a):
    n=len(a); full=(1<<n)-1; h=counter(a); w=[[0]*n for _ in a]
    def walk(start,v,used,length):
        if length%2==0:w[start][v]+=h(full^used)
        available=a[v]&~used
        while available:
            b=available&-available; available-=b
            walk(start,b.bit_length()-1,used|b,length+1)
    for s in range(n):walk(s,s,1<<s,1)
    return w

def extend(a,x):
    n=len(a); b=list(a)+[x]
    for i in range(n):
        if not x>>i&1:b[i]|=1<<n
    return tuple(b)

def cut(w,x):
    return sum(w[i][j] for i in range(len(w)) for j in range(len(w)) if x>>i&1 and not x>>j&1)

def unit_candidates(w):
    n=len(w); found=[]
    for a in range(n):
        for b in range(n):
            if w[a][b]!=1:continue
            seen=1<<a
            while True:
                old=seen
                for u in range(n):
                    if seen>>u&1:
                        for v in range(n):
                            if w[u][v] and (u,v)!=(a,b):seen|=1<<v
                if old==seen:break
            if not seen>>b&1:found.append((a,b,seen))
    return found

def seed_coordinates(a):
    n=len(a); full=(1<<n)-1
    starts=[[0]*n for _ in range(full+1)]; ends=[[0]*n for _ in range(full+1)]
    for mask in range(1,full+1):
        for p in permutations([v for v in range(n) if mask>>v&1]):
            if all(a[u]>>v&1 for u,v in zip(p,p[1:])):
                starts[mask][p[0]]+=1;ends[mask][p[-1]]+=1
    K=[[0]*n for _ in a]
    for mask in range(1,full):
        for i in range(n):
            for j in range(n):K[i][j]+=ends[mask][i]*starts[full^mask][j]
    return starts[full],ends[full],K

def run():
    report={'status':'FINITE-EXACT','universe':'all labelled tournaments of orders 1..5, all newborn cuts','tournaments':0,'births':0,'unit_parent_counts':{},'controls':{},'seed_budget_checks':'all selected subsets; independently enumerated ordered two-path covers'}
    unit=Counter(); examples={}
    for n in range(1,6):
        for code in range(1<<(n*(n-1)//2)):
            a=rows(n,code); h=counter(a)((1<<n)-1); w=capacities(a)
            s,e,K=seed_coordinates(a)
            for O in range(1<<n):
                inside=sum(w[i][j] for i in range(n) for j in range(n) if O>>i&1 and O>>j&1)
                volume=sum(w[i][j] for i in range(n) for j in range(n) if O>>i&1)
                original_W=sum(K[i][j] for i in range(n) for j in range(n) if O>>i&1 and O>>j&1)
                original_U=sum(s[j]+sum(K[i][j] for i in range(n)) for j in range(n) if O>>j&1)
                require(original_W==2*inside,'seed_internal_budget')
                require(original_U==sum(e[j] for j in range(n) if O>>j&1)+2*volume,'seed_volume_budget')
                if original_W==h-1 and h==1:
                    require(original_U==0 or original_U&(original_U-1)==0,'transitive_slope')
            vals=[]
            for x in range(1<<n):
                d=cut(w,x); actual=counter(extend(a,x))((1<<(n+1))-1)
                require(actual==h+2*d,('cut_identity',a,x,actual,h,d))
                if d==0:
                    require(all(not(a[i]>>j&1) for i in range(n) for j in range(n) if x>>i&1 and not x>>j&1),'zero_cut')
                vals.append(d);report['births']+=1
            candidates=unit_candidates(w)
            require(bool(candidates)==(1 in vals),('unit_classification',a,w,vals))
            for u,v,x in candidates:require(cut(w,x)==1,'unit_constructor')
            if candidates:require(h in (1,3),'sharp_plus_two_obstruction');unit[h]+=1;examples.setdefault(h,{'parent':a,'cut':candidates[0][2],'child_count':h+2})
            # Exact heavy-arc implication pruning, tested at every observed target.
            for r in set(vals):
                feasible=[x for x in range(1<<n) if all(not(x>>i&1) or x>>j&1 for i in range(n) for j in range(n) if w[i][j]>r)]
                require({x for x in feasible if cut(w,x)==r}=={x for x,d in enumerate(vals) if d==r},'heavy_pruning')
            report['tournaments']+=1
    report['unit_parent_counts']=dict(sorted(unit.items()));report['unit_witnesses']=examples
    tri=(2,4,1); w=capacities(tri)
    require(sorted(set(cut(w,x) for x in range(8)))==[0,1],'triangle_positive')
    # Strong four-vertex tournament with H=5: no +2 birth, despite C1=1.
    hostile=next(rows(4,c) for c in range(64) if counter(rows(4,c))(15)==5)
    spectra=sorted(set(counter(extend(hostile,x))(31) for x in range(16)))
    from math import gcd
    from functools import reduce
    require(reduce(gcd,spectra)==1 and 7 not in spectra,'gcd_hostile')
    report['controls']={'triangle':{'parent':tri,'half_increment_spectrum':[0,1]},'gcd_one_not_support':{'parent':hostile,'H':5,'extension_spectrum':spectra,'minimum_positive_half_increment':min(cut(capacities(hostile),x) for x in range(16) if cut(capacities(hostile),x)>0)}}
    return report
if __name__=='__main__':
    result=run();print(json.dumps(result,indent=2))
