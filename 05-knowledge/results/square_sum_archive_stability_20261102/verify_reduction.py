#!/usr/bin/env python3
"""Exact cycle ledger and independent weighted-matrix consequence controls."""
from fractions import Fraction as F
from collections import defaultdict
from itertools import combinations
from math import factorial
import random,json

def need(b,m):
    if not b:raise RuntimeError(m)
def parts(n,m=2):
    if n==0:yield ()
    for a in range(m,n+1):
        for t in parts(n-a,a):yield (a,)+t

def ledger(n,branch):
    out={}
    for s in range(2,n+1):
        for t in parts(s):
            h=2**len(t) if all(x%2 for x in t) else 0
            out[t]=h-(1 if s<n else branch)*(-1)**sum(x%2==0 for x in t)
    return out

def complement_det(n,sign):
    # u + sign*v, with v=(-1)^n det A.
    out={():1}
    for s in range(2,n+1):
        for t in parts(s):out[t]=(-1)**sum(x%2==0 for x in t)+(sign*(-1)**len(t) if s==n else 0)
    return out

def lift(k,p):
    out=defaultdict(F)
    for t,c in p.items():out[tuple(sorted(t+(k,)))]+=c*(t.count(k)+1)
    return dict(out)
def add(*ps):
    out=defaultdict(F)
    for p in ps:
        for t,c in p.items():out[t]+=c
    return {t:c for t,c in out.items() if c}
def scale(c,p):return {t:c*v for t,v in p.items()}
def fmt(p):return {','.join(map(str,t)):str(c) for t,c in sorted(p.items()) if c}

def symbolic():
    P={(2,):1,(2,2):-1,(2,2,2):1,(2,2,2,2):-1}
    tri={(3,):1,(3,3):3,(2,2,3):-1,(2,2,3,3):-2}
    six={(6,):1,(2,6):-1,(3,6):1}
    five={(5,):1,(2,2,5):-1}
    eight={(8,):1,(2,8):-2}
    leftover={(2,3):1,(7,):1,(9,):1}
    rows=[]
    for branch in [0,2]:
        four=add(lift(4,complement_det(6,-1 if branch==0 else 1)),{(4,4):1},({(2,4,4):-2} if branch==2 else {}))
        remainder=add(ledger(10,branch),scale(-1,add(P,tri,six,five,eight,four,leftover)))
        need(all(c>=0 for c in remainder.values()),'negative unallocated cycle coefficient')
        rows.append({'full_cover_multiplier':branch,'four_cycle_group':fmt(four),'nonnegative_remainder':fmt(remainder),'complete_ledger':fmt(ledger(10,branch))})
    a=F(21,16);quad=F(1,2)-(1+a)*F(2,7);lin=(1+a)/4
    max_loss=-lin*lin/(4*quad)
    need(max_loss==F(9583,18432),'triangle square completion')
    need(1-max_loss==F(8849,18432),'positive triangle reserve')
    return rows

def det(a):
    a=[list(map(F,r)) for r in a];n=len(a);out=F(1)
    for k in range(n):
        i=next((i for i in range(k,n) if a[i][k]),None)
        if i is None:return F(0)
        if i!=k:a[k],a[i]=a[i],a[k];out=-out
        pivot=a[k][k];out*=pivot
        for i in range(k+1,n):
            q=a[i][k]/pivot
            for j in range(k+1,n):a[i][j]-=q*a[k][j]
    return out

def hp(a):
    n=len(a);dp=[{} for _ in range(1<<n)]
    for v in range(n):dp[1<<v][v]=F(1)
    for mask in range(1,1<<n):
        for v,c in dp[mask].items():
            for w in range(n):
                if not mask>>w&1 and a[v][w]:
                    r=dp[mask|1<<w];r[w]=r.get(w,F(0))+c*a[v][w]
    return sum(dp[-1].values(),F(0))

def cycles(a):
    n=len(a);out=[F(0) for _ in range(n+1)]
    # Each cycle is counted at its least labelled vertex, with direction retained.
    for root in range(n):
        others=list(range(root+1,n));m=len(others);dp=[{} for _ in range(1<<m)]
        for j,v in enumerate(others):
            if a[root][v]:dp[1<<j][j]=a[root][v]
        for mask,row in enumerate(dp):
            size=mask.bit_count()+1
            for j,c in row.items():
                v=others[j]
                if size>=3:out[size]+=c*a[v][root]
                for k,w in enumerate(others):
                    if not mask>>k&1 and a[v][w]:
                        rr=dp[mask|1<<k];rr[k]=rr.get(k,F(0))+c*a[v][w]
    return out

def matrix_from_pairs(ps,n=10):
    a=[[F(0) for _ in range(n)] for _ in range(n)]
    for (i,j),p in zip(combinations(range(n),2),ps):a[i][j]=F(p);a[j][i]=1-F(p)
    return a

def weighted_controls():
    n=10;edges=list(combinations(range(n),2));rng=random.Random(20261102)
    inputs=[('all_fair',[F(1,2)]*45),('transitive',[F(1)]*45)]
    fair=[F(1,2) if i//2==j//2 else F(1) for i,j in edges];inputs.append(('five_fair_pairs',fair))
    for j in range(16):inputs.append(('dyadic_'+str(j),[F(rng.randrange(5),4) for _ in edges]))
    for idx in [0,1,8,17,30,44]:
        x=fair.copy();x[idx]=F(15,16) if x[idx]==1 else F(7,16);inputs.append(('near_fair_'+str(idx),x))
    rows=[]
    for name,ps in inputs:
        a=matrix_from_pairs(ps);H=hp(a);u=det([[a[i][j]+int(i==j) for j in range(n)] for i in range(n)]);v=det(a)
        V=sum((p*(1-p) for p in ps),F(0));C=cycles(a)
        rhs=F(39,64)*V+F(8849,18432)*C[3]+F(1,32)*C[4]+F(1,16)*C[5]+F(1,2)*C[6]+C[7]+F(1,2)*C[8]+C[9]
        slack=H-u-abs(v)-rhs;need(slack>=0,'weighted consequence '+name)
        if name=='five_fair_pairs':need(slack==0 and V==F(5,4) and H-u-abs(v)==F(195,256),'sharp witness')
        rows.append({'name':name,'H':str(H),'u':str(u),'v':str(v),'V':str(V),'tau':str(C[3]),'claimed_inequality_slack':str(slack)})
    return rows

def main():
    print(json.dumps({'status':'passed','symbolic_ledgers':symbolic(),'weighted_controls':weighted_controls(),'scope':'Cycle grouping is exact for all matrices. Matrix examples are controls; the full real-box proof additionally uses the exhaustive graph endpoint certificate.'},indent=2,sort_keys=True))
if __name__=='__main__':main()
