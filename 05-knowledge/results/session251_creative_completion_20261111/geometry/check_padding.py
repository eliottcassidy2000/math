#!/usr/bin/env python3
"""Exact new padding normal form; exhaustive labelled bases through order 5."""
from itertools import combinations
import json

def tournament(n,code):
    a=[0]*n
    for k,(i,j) in enumerate(combinations(range(n),2)):
        if code>>k&1: a[i]|=1<<j
        else: a[j]|=1<<i
    return a

def counts(a):
    n=len(a); dp=[[0]*n for _ in range(1<<n)]; h=[0]*(1<<n);h[0]=1
    for i in range(n):dp[1<<i][i]=1
    for S in range(1,1<<n):
        for i in range(n):
            if not (S>>i&1):continue
            if S!=1<<i:
                dp[S][i]=sum(dp[S^(1<<i)][j] for j in range(n) if a[j]>>i&1)
        h[S]=sum(dp[S])
    return h

def scc_order(a):
    n=len(a); reach=[a[i]|(1<<i) for i in range(n)]
    for k in range(n):
        for i in range(n):
            if reach[i]>>k&1: reach[i]|=reach[k]
    unused=set(range(n)); cs=[]
    while unused:
        i=min(unused); C={j for j in unused if reach[i]>>j&1 and reach[j]>>i&1}
        cs.append(C);unused-=C
    old=cs[:]
    cs.sort(key=lambda C: -sum(bool(a[next(iter(C))]>>next(iter(D))&1) for D in old if C is not D))
    return cs

def extension(a,P):
    n=len(a);b=a[:]+[P]
    for i in range(n):
        if not(P>>i&1):b[i]|=1<<n
    return b

def capacities(a,h):
    n=len(a);full=(1<<n)-1;c=[[0]*n for _ in range(n)]
    def visit(i,j,S,k):
        if k%2==0:c[i][j]+=h[full^S]
        for v in range(n):
            if not(S>>v&1) and a[j]>>v&1:visit(i,v,S|1<<v,k+1)
    for i in range(n):visit(i,i,1<<i,1)
    return c

def run():
    total=low=inv_cases=branches=capacity_checks=flank_checks=0;first=None
    for n in range(1,6):
        for code in range(1<<(n*(n-1)//2)):
            a=tournament(n,code);h=counts(a);H=h[-1];C=scc_order(a)
            singles=[next(iter(c)) for c in C if len(c)==1]
            positions={next(iter(c)):k for k,c in enumerate(C) if len(c)==1}
            cap=capacities(a,h)
            for u,v in combinations(singles,2):
                middle_set=set().union(*C[positions[u]+1:positions[v]])
                middle=len(middle_set); M=sum(1<<x for x in middle_set)
                even_sum=sum(h[S]*h[M^S] for S in range(1<<n) if S&M==S and S.bit_count()%2==0)
                assert cap[u][v]*h[M]==H*even_sum
                assert cap[u][v]>=(middle//2+1)*H
                assert cap[v][u]==0
                capacity_checks+=1
            for P in range(1<<n):
                N=counts(extension(a,P))[-1];total+=1
                cut=sum(cap[i][j] for i in range(n) for j in range(n) if P>>i&1 and not(P>>j&1))
                assert N==H+2*cut
                inv=[(u,v) for u,v in combinations(singles,2) if P>>u&1 and not(P>>v&1)]
                assert N>=H*(1+2*len(inv))
                if N<5*H:
                    low+=1;assert len(inv)<=1
                    nontrivial=[c for c in C if len(c)>1]
                    if len(nontrivial)==1:
                        core=nontrivial[0]; ci=C.index(core)
                        root_component=next(c for c in scc_order(extension(a,P)) if n in c)
                        if core<=root_component:
                            flank_checks+=1
                            left=set().union(*C[:ci]);right=set().union(*C[ci+1:])
                            assert not(left&root_component and right&root_component)
                    if inv:
                        inv_cases+=1;u,v=inv[0]
                        assert positions[v]==positions[u]+1
                        assert N==3*H
                        left=set().union(*C[:positions[u]]);right=set().union(*C[positions[v]+1:])
                        assert all(not(P>>i&1) for i in left)
                        assert all(P>>i&1 for i in right)
                        if H>=3 and first is None:first={'rows':a,'outset':P,'h':H,'child':N,'pair':[u,v]}
                if H>=3 and N in (3*H+2,3*H+4):
                    branches+=1;assert not inv
    # Scaled bound is sharp, and strict threshold cannot be weakened:
    # TT4, beaten set {0} has two inversions when rows oriented ascending? choose directly.
    a=[14,12,8,0]; P=2 # 0100: two inversions, child 5.
    assert counts(extension(a,P))[-1]==5
    return dict(labelled_attachment_checks=total,low_ratio_checks=low,
                low_ratio_unique_inversion_cases=inv_cases,branch_cases=branches,
                singleton_interval_checks=capacity_checks,one_sided_core_checks=flank_checks,nontrivial_equality_example=first,
                threshold_hostile={'rows':a,'outset':P,'h':1,'child':5,'inversions':2})
if __name__=='__main__':print(json.dumps(run(),indent=2,sort_keys=True))
