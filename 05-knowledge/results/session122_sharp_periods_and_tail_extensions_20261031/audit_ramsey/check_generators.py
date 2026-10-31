#!/usr/bin/env python3
"""Independent finite-generator controls; exact already-proved rows2..8.
Uses auditor routines only. Bellman-Ford differs from producer heap search.
"""
from fractions import Fraction as Q
from pathlib import Path
import json
from check import weight, patterns, vertices, solve, require

ROWS={1:(1,(0,)),2:(2,(0,2)),3:(4,(0,3,7,11)),4:(4,(0,4,10,16)),
      5:(8,(0,5,12,19,29,35,43,51)),6:(8,(0,6,15,23,35,43,55,65)),
      7:(8,(0,7,17,27,41,51,65,79)),8:(8,(0,8,20,32,50,62,80,98))}

def row(a,b):
    L,h=ROWS[a]
    z=weight(a)*L*(b//L)+h[b%L]
    require(z.denominator==1,'row integral')
    return int(z)

def bound(a):
    s=weight(a);L=1<<(a-1).bit_length();D=(s-a)*(L-1);M=max(1,L)
    details=[];vv=vertices()
    for x in range(1,a):
        gap=s-weight(x)-weight(a-x);require(gap>0,'join gap')
        m=int(D//gap);M=max(M,m);details.append(('join',m))
    for aa in patterns(a,1,True):
        co=[weight(x) for x in aa]
        half=sum(co)/2
        if half==s:
            mat=[]
            for i in range(5):
                r=[0]*5;r[i]=r[(i-2)%5]=1;mat.append(r)
            dual=solve(mat,co);require(dual and min(dual)>0,'positive sharp dual')
            S=int(D//min(dual));N=1
            m=2*N+S+2*L-1;details.append(('sharp',m))
        else:
            v=max(sum(c*x for c,x in zip(co,v)) for v in vv)
            require(v<s,'nonsharp vertex gap')
            m=int(D//(s-v));details.append(('other',m))
        M=max(M,m)
    return M,details

def shortest(L,p,q,items):
    # Finite simple-path Bellman-Ford; all lexicographic edge weights positive.
    d=[None]*L;paths=[None]*L;d[0]=(0,0);paths[0]=()
    for _ in range(L-1):
        old=d[:];oldpaths=paths[:]
        for r in range(L):
            if old[r] is None:continue
            for b,n in items:
                cost=p*b-q*n;require(cost>=0,'nonnegative deficit')
                dest=(r+b)%L;value=(old[r][0]+cost,old[r][1]+b)
                if d[dest] is None or value<d[dest]:
                    d[dest]=value;paths[dest]=oldpaths[r]+(b,)
    require(all(z is not None for z in d),'residue reachability')
    # Check no edge can improve a completed distance, including zero-cost edges.
    for r in range(L):
        for b,n in items:
            val=(d[r][0]+p*b-q*n,d[r][1]+b)
            require(d[(r+b)%L]<=val,'Bellman optimality')
    return d,paths

def knapsack(items,limit):
    f=[0]+[-1]*limit
    for n in range(1,limit+1):
        f[n]=max((f[n-b]+profit for b,profit in items if b<=n and f[n-b]>=0),default=-1)
    return f

def main():
    rows=[];total=0
    for a in range(2,9):
        s=weight(a);p,q=s.numerator,s.denominator;L=ROWS[a][0]
        M,details=bound(a);items=[(b,row(a,b)) for b in range(1,M+1)]
        d,paths=shortest(L,p,q,items)
        for r,(cost,length) in enumerate(d):
            require(length<=max(0,(L-1)*M),'simple path quantitative bound')
            require(sum(paths[r])==length and length%L==r,'actual shortest path length')
            require(sum(row(a,b) for b in paths[r])==Q(p*length-cost,q),'actual path profit')
            require(cost==p*r-q*row(a,r),'correct eventual residue deficit')
        limit=max(64,(L-1)*M+L)
        f=knapsack(items,limit)
        require(all(f[b]==row(a,b) for b in range(1,limit+1)),'entire finite coin closure')
        require(all(f[b+L]-f[b]==s*L for b in range(1,limit-L+1)),'known threshold1')
        rows.append({'a':a,'M':M,'bound_types':{kind:sum(k==kind for k,_ in details) for kind in ('join','sharp','other')},
                     'path_pairs':d,'path_lengths':list(map(list,paths)),
                     'verified_coin_prefix':limit,'safe_tail_threshold':max(1,max(z[1] for z in d)),'least_threshold':1})
        total+=limit
    hostile=[]
    for late in (3,11,101):
        items=[(1,1),(2,4),(late,2*late)]
        d,paths=shortest(2,2,1,items)
        require(d==[(0,0),(0,late)],'late odd generator distance')
        f=knapsack(items,late+8)
        require(all(f[b]==2*b-int(b%2==1 and b<late) for b in range(late+9)),'delayed plateau actual knapsack')
        lastbad=max((b for b in range(1,late+1) if f[b]!=2*b),default=0)
        require(lastbad+1==late-1,'delayed least threshold')
        hostile.append({'late_generator':late,'odd_first_zero_cost_weight':d[1][1],'least_threshold':lastbad+1})
    print(json.dumps({'status':'PASS','scope':'independent finite-generator algorithm controls for certified rows2..8; all-A proof separately audited',
                     'rows':rows,'total_prefix_values':total,'late_generator_hostiles':hostile},sort_keys=True,indent=2))

if __name__=='__main__':main()
