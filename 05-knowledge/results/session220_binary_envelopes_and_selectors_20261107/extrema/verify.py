#!/usr/bin/env python3
"""Exact bounded checks of the all-length extremal theorems; no coverage claim."""
import json
from pathlib import Path
from collections import deque
states=[tuple(map(int,s.split())) for s in Path(__file__).with_name('CORE_STATES.txt').read_text().splitlines() if s.strip()]
def f(n):
    out=0
    while n:
        p=n//2
        if n%2: out+=p+p%2
        n=p
    return out
def fp(n):
    v=[i for i in range(n.bit_length()) if n>>i&1]
    return sum((1<<(j-i-1))+(j==i+1) for i in v for j in v if i<j)
def rev(n,b): return sum(((n>>i)&1)<<(b-1-i) for i in range(b))
def value(c,p,r,b):
    h,q,A,G,T,S,E=c
    return (q*p+A)*((1<<(b+1))-rev(r,b))+h*(2*r+1-f(r))-p*G-h*f(p)-T-(p%2)*S-E*(r>>(b-1) if b else 0)
def compiler(c,p,a,r,b):
    h,q,A,G,T,S,E=c
    # Prepend prefix bits from right to left.
    for i in range(a):
        z=p>>i&1
        q,A,G,T,S,E=2*q,A+z*q,2*G+z*h,T+z*(G+S),z*h,E
    for i in range(b-1,-1,-1):
        z=r>>i&1
        q,A,G,T,S,E=2*q,2*A+z*h,G+z*q,T+z*(A+E),S,z*h
    return 2*A+h-T

def mincut(c,prefix_bits,suffix_bits):
    # Values None are free, 0/1 fixed. Index order is low-to-high per word.
    a,b=len(prefix_bits),len(suffix_bits)
    fixed=prefix_bits+suffix_bits
    free=[i for i,v in enumerate(fixed) if v is None]
    base=sum((v or 0)<<i for i,v in enumerate(prefix_bits))
    baser=sum((v or 0)<<i for i,v in enumerate(suffix_bits))
    def fun(ones):
        p,r=base,baser
        for i in ones:
            if i<a:p|=1<<i
            else:r|=1<<(i-a)
        return value(c,p,r,b)
    C=fun([]);L=[fun([i])-C for i in free]
    n=len(free);cs=[];d=[2*l for l in L]
    for i in range(n):
        for j in range(i):
            w=-(fun([free[i],free[j]])-C-L[i]-L[j])
            assert w>=0
            if w:cs.append((i,j,w));d[i]-=w;d[j]-=w
    s,t=n,n+1; cap=[[0]*(n+2) for _ in range(n+2)]
    for i,j,w in cs:cap[i][j]+=w;cap[j][i]+=w
    for i,v in enumerate(d):
        cap[s][i]+=max(-v,0);cap[i][t]+=max(v,0)
    flow=0
    while True:
        pred=[-1]*(n+2);pred[s]=s;todo=deque([s])
        while todo and pred[t]<0:
            u=todo.popleft()
            for v in range(n+2):
                if cap[u][v] and pred[v]<0:pred[v]=u;todo.append(v)
        if pred[t]<0:break
        aug=sum(cap[s]);v=t
        while v!=s:u=pred[v];aug=min(aug,cap[u][v]);v=u
        v=t
        while v!=s:u=pred[v];cap[u][v]-=aug;cap[v][u]+=aug;v=u
        flow+=aug
    ans=2*C+sum(min(v,0) for v in d)+flow
    assert ans%2==0
    return ans//2

def main():
    checks=0; extrema=0;cuts=0;cylinders=0;scalar=0
    def g(n):
        if not n:return 1
        u=n//2
        return g(u)+(u+u%2 if n%2==0 else 0)
    for n in range(1<<13): assert g(n)+f(n)==n+1-n%2
    for m in range(1,10):
        for K in (1,3,9,32):
            for h in (1,5,9):
                for S in (0,2,20):
                    fun=lambda n:K*n-h*f(n)-S*(n%2)
                    cand=[0]+[(1<<m)-(1<<j) for j in range(m)]
                    assert max(map(fun,range(1<<m)))==max(map(fun,cand))
                    scalar+=1
    for n in range(1<<13):
        assert f(n)==fp(n)
        # Every applicable internal-zero rounding, not just a chosen rounding.
        for k in range(1,n.bit_length()):
            U=n>>(k+1);V=n% (1<<k)
            if U and V and not(n>>k&1):
                new=(2*U+1)<<k
                assert new>n and f(new)<=f(n)
                assert rev(new,n.bit_length())<rev(n,n.bit_length())
                checks+=1
    for ci,c in enumerate(states):
        for a in range(9):
            for b in range(9-a):
                if not a and ci in (0,16,48):continue
                pp=range(1<<(a-1),1<<a) if a else [0]
                rr=range(1<<b)
                vals=[value(c,p,r,b) for p in pp for r in rr]
                pcan=[(1<<a)-(1<<j) for j in range(a)] if a else [0]
                rcan=[0]+[((1<<l)-1)<<j for l in range(1,b+1) for j in range(b-l+1)]
                assert max(vals)==max(value(c,p,r,b) for p in pcan for r in rcan)
                extrema+=1
                for p in pp:
                    for r in rr:
                        assert value(c,p,r,b)==compiler(c,p,a,r,b)
                if a+b<=6:
                    pb=[None]*a
                    if a:pb[-1]=1
                    assert min(vals)==mincut(c,pb,[None]*b)
                    cuts+=1
                    rb=[None]*b
                    # Interior bits fixed on both sides; leading bit stays legal.
                    for i in range(a-1):
                        if (i+ci)%3==0:pb[i]=(i+ci)%2
                    for i in range(b):
                        if (i+ci)%3==1:rb[i]=(i+ci)%2
                    sub=[value(c,p,r,b) for p in pp for r in rr
                         if all(v is None or (p>>i&1)==v for i,v in enumerate(pb))
                         and all(v is None or (r>>i&1)==v for i,v in enumerate(rb))]
                    assert min(sub)==mincut(c,pb,rb)
                    cylinders+=1
    # Capacities exceed 10^100: augmentation must use the true path bottleneck.
    pb=[0]*400;pb[-1]=1;pb[398]=None;pb[0]=None
    large_values=[value(states[0],(1<<399)+(x<<398)+y,0,0)
                  for x in (0,1) for y in (0,1)]
    assert mincut(states[0],pb,[])==min(large_values)
    # Negative control: min/max interval does not imply occupancy.
    c=states[0];outs={value(c,p,0,0) for p in range(8,16)}
    holes=[n for n in range(min(outs),max(outs)+1,2) if n not in outs]
    assert holes
    print(json.dumps({'status':'PASS','core_states':len(states),'rounding_checks':checks,'fixed_length_maximum_checks':extrema,'minimum_cut_checks':cuts,'fixed_interior_cylinder_checks':cylinders,'general_scalar_maximum_checks':scalar,'join_bridge_checks':8192,'large_capacity_400_bit_check':'PASS','hostile':{'core_index':0,'a':4,'b':0,'support':sorted(outs),'first_missing_odd_in_hull':holes[0]},'scope':'Extremal identities only; all-target coverage remains open.'},indent=2))
if __name__=='__main__':main()
