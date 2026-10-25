#!/usr/bin/env python3
"""Exact independent score/forest/partition checks; Python standard library only."""
from functools import lru_cache
from itertools import combinations, product
from math import comb, factorial
import json
from time import perf_counter
from pathlib import Path

def require(condition, detail="exact check failed"):
    if not condition: raise RuntimeError(detail)

def edges(n): return [(i,j) for i in range(n) for j in range(i+2,n)]
def direct(n):
    es=edges(n); scores=set(); counts=[0]*n
    for mask in range(1<<len(es)):
        d=[0]*n; p=list(range(n)); forest=True; k=0
        def root(u):
            while p[u]!=u: u=p[u]
            return u
        for j,(u,v) in enumerate(es):
            d[v if mask>>j&1 else u]+=1
            if mask>>j&1:
                k+=1; a,b=root(u),root(v)
                if a==b: forest=False
                else: p[a]=b
        scores.add(tuple(d))
        if forest: counts[k]+=1
    return len(scores),counts

def trim(a):
    a=list(a)
    while len(a)>1 and a[-1]==0: a.pop()
    return tuple(a)
@lru_cache(None)
def weighted(parts):
    """Weighted forest polynomial on distinct blocks with sizes parts, sorted."""
    if not parts: return (1,)
    r=parts[0]; rest=parts[1:]
    kinds=sorted(set(rest)); mult=[rest.count(i) for i in kinds]
    out=[0]*len(parts)
    for b in product(*(range(m+1) for m in mult)):
        selected=[i for i,c in zip(kinds,b) for _ in range(c)]
        remaining=tuple(i for i,m,c in zip(kinds,mult,b) for _ in range(m-c))
        h=len(selected); ways=1
        for m,c in zip(mult,b): ways*=comb(m,c)
        w=1
        if h:
            w=(r+sum(selected))**(h-1)*r
            for i in selected:w*=i
        for k,v in enumerate(weighted(remaining)): out[k+h]+=ways*w*v
    return trim(out)

def partitions(n,minimum=1):
    if n==0: yield (); return
    for x in range(minimum,n+1):
        for rest in partitions(n-x,x): yield (x,)+rest

def ie(n):
    out=[0]*max(n,1)
    for parts in partitions(n):
        k=len(parts); mult=factorial(k)
        for x in set(parts): mult//=factorial(parts.count(x))
        shift=n-k; sign=(-1)**shift
        for j,v in enumerate(weighted(parts)): out[shift+j]+=sign*mult*v
    return trim(out)

def det(a):
    a=[r[:] for r in a]; n=len(a)
    if not n:return 1
    prev=1; sign=1
    for k in range(n-1):
        if not a[k][k]:
            j=next((j for j in range(k+1,n) if a[j][k]),None)
            if j is None:return 0
            a[k],a[j]=a[j],a[k];sign=-sign
        pivot=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                v=a[i][j]*pivot-a[i][k]*a[k][j]
                require(v%prev==0)
                a[i][j]=v//prev
        for i in range(k+1,n):a[i][k]=0
        prev=pivot
    return sign*a[-1][-1]
def subset(n, fast=False):
    @lru_cache(None)
    def tree(U):
        vs=[i for i in range(n) if U>>i&1]; z=len(vs)
        if z==1:return 1
        if fast:
            runs=[]; r=1
            for a,b in zip(vs,vs[1:]):
                if b==a+1:r+=1
                else:runs.append(r);r=1
            runs.append(r); numerator=1
            for r in runs:
                if r==1: block=z
                else:
                    prev,current=1,z-1
                    for j in range(1,r):
                        d=z-1 if j==r-1 else z-2
                        prev,current=current,d*current-prev
                    block=current
                numerator*=block
            require(numerator%(z*z)==0)
            return numerator//(z*z)
        L=[[0]*z for _ in vs]
        for i,j in combinations(range(z),2):
            if abs(vs[i]-vs[j])!=1:
                L[i][i]+=1;L[j][j]+=1;L[i][j]-=1;L[j][i]-=1
        return det([r[:-1] for r in L[:-1]])
    @lru_cache(None)
    def calc(U):
        if not U:return (1,)
        root=U&-U; rest=U^root; V=rest; out=[0]*U.bit_count()
        while True:
            B=V|root; w=tree(B); shift=B.bit_count()-1
            for j,v in enumerate(calc(U^B)):out[j+shift]+=w*v
            if not V:break
            V=(V-1)&rest
        return trim(out)
    return calc((1<<n)-1)

def main():
    rows=[]; start=perf_counter()
    for n in range(1,15):
        p=ie(n); q=subset(n,fast=True); require(p==q,(n,p,q))
        if n<=12: require(p==subset(n))
        row={'n':n,'forest_coefficients':p,'score_fibres':sum(p),'methods':['weighted_composition_inclusion_exclusion','subset_tree_partition_recurrence_with_path_run_determinants']}
        if n<=7:
            score,f=direct(n); require(score==sum(p) and trim(f)==p)
            row['direct_orientation_scores']=score;row['direct_forest_coefficients']=f
        rows.append(row)
    require([rows[n-1]['score_fibres'] for n in (5,6,7)]==[52,478,5706])
    triangle=weighted((1,1,1)); require(triangle==(1,3,3))
    require(det([[3,-1,-1],[-1,3,-1],[-1,-1,3]])==16)
    result={'rows':rows,'triangle_control':{'unrooted_forest_polynomial':triangle,'unrooted_total':sum(triangle),'det_I_plus_L':16},'weighted_cache_states':weighted.cache_info().currsize,'elapsed_seconds':perf_counter()-start}
    Path(__file__).with_name('certificate.json').write_text(json.dumps(result,indent=2)+'\n')
    for row in rows: print(row['n'],row['score_fibres'],list(row['forest_coefficients']))
    print('All exact checks passed; triangle 7 versus 16.')
if __name__=='__main__':main()
