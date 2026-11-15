#!/usr/bin/env python3
"""Literal directed subset DP for the exact swaps and interface certificates."""
from itertools import combinations, permutations
from functools import lru_cache
from math import prod
from random import Random

def strong(A):
    full=(1<<len(A))-1
    for i in range(len(A)):
        s=1<<i
        for _ in A:
            for j in range(len(A)):
                if s>>j&1:s|=A[j]
        if s!=full:return False
    return True

def paths(A):
    n=len(A); N=1<<n; end=[[0]*n for _ in range(N)]; H=[1]+[0]*(N-1)
    for i in range(n):end[1<<i][i]=1
    for s in range(1,N):
        row=end[s]
        H[s]=sum(row)
        for i,w in enumerate(row):
            if not w:continue
            t=A[i]&~s
            while t:
                b=t&-t;t^=b
                end[s|b][b.bit_length()-1]+=w
    return H

def slot(A,k):
    h=paths(A); prev=[1]+[0]*(len(h)-1)
    out=[prev[-1]]
    for _ in range(k):
        new=[0]*len(h)
        for s in range(len(h)):
            t=s
            while True:
                new[s]+=prev[t]*h[s^t]
                if not t:break
                t=(t-1)&s
        prev=new
        out.append(prev[-1])
    return out

def compose(Q,modules):
    v=[(i,j) for i,a in enumerate(modules) for j in range(len(a))]
    return tuple(sum(1<<t for t,(c,d) in enumerate(v)
                     if ((modules[a][b]>>d)&1 if c==a else (Q[a]>>c)&1))
                 for a,b in v)

def cycle_lengths(A, supports=False):
    n=len(A); C=[0]*(n+1); D={}
    for root in range(n):
        def step(u,s,d):
            if d>=3 and A[u]>>root&1:
                C[d]+=1
                if d%2:D[s]=D.get(s,0)+1
            t=A[u]&~s&~((1<<(root+1))-1)
            while t:
                b=t&-t;t^=b
                step(b.bit_length()-1,s|b,d+1)
        step(root,1<<root,1)
    return (C,D) if supports else C

@lru_cache(None)
def schedules(Q,m,last=-1):
    if not any(m):return 1
    z=0
    for i,r in enumerate(m):
        if r and (last<0 or Q[last]>>i&1):
            q=list(m);q[i]-=1;z+=schedules(Q,tuple(q),i)
    return z

def main():
    U=(120,85,89,114,192,86,128,47);V=(158,92,64,116,196,23,161,46)
    assert strong(U) and strong(V)
    assert slot(U,3)==[0,345,5640,44505] and slot(V,3)==[0,345,5640,45225]
    # z=0, O_1=1,O_2=2, I_1=3,I_2=4.
    Q=(6,12,16,21,3)
    a,b=[compose(Q,(W,(0,),(0,),(0,),(0,))) for W in (U,V)]
    assert strong(a) and strong(b)
    low=paths(a)[-1];high=paths(b)[-1];assert high-low==1440
    U=(24,17,3,6,12);V=(16,9,3,5,14)
    assert strong(U) and strong(V)
    assert paths(U)[-1]==paths(V)[-1]==15
    assert (cycle_lengths(U)[5],cycle_lengths(V)[5])==(2,3)
    Q=(6,12,8,1)
    a,b=[compose(Q,(*x,(2,0),(2,0))) for x in ((U,V),(V,U))]
    assert strong(a) and strong(b)
    assert (paths(a)[-1],paths(b)[-1])==(2270565,2286945)
    rng=Random(91)
    for n in range(3,9):
        for _ in range(10):
            A=[0]*n
            for i in range(n):
                for j in range(i+1,n):
                    if rng.randrange(2):A[i]|=1<<j
                    else:A[j]|=1<<i
            _,H,W,P3,P4,P5=slot(A,5); c,D=cycle_lengths(A,True)
            drift=sum((u.bit_count()-1)*v for u,v in D.items())
            drift+=sum(2*(u.bit_count()+v.bit_count()-2)*a*b
                       for (u,a),(v,b) in combinations(D.items(),2) if not u&v)
            assert (P3-3**n*(H-8*drift))%64==0
            assert (P5-5**n*H)%16==0
            if n>=5:assert P4%64==((8*c[n]) if n%2 else (32*c[n-1]))%64
            J=sum(c[k] for k in range(3,n+1,4))%2
            assert (P3-3**n*(H+16*J))%32==0
            if n>=3:assert (6*(P4-4*P3+6*W-4*H)-16*(1-n%2))%32==0
            if n>=5 and n%2:assert W%16==(4*c[n])%16
            if n>=4 and n%2==0:assert W%16==(8*c[n-1])%16
    # Full simultaneous formula at a few mixed module sizes.
    for _ in range(16):
        d=3+(_%2); Q=[0]*d
        for i in range(d):
            for j in range(i+1,d):
                if rng.randrange(2):Q[i]|=1<<j
                else:Q[j]|=1<<i
        Q=tuple(Q)
        modules=[rng.choice(((2,0),(6,4,0),(2,4,1))) for i in range(d)]
        v=paths(compose(Q,modules))[-1]; h=[paths(u)[-1] for u in modules]
        A=[schedules(Q,tuple(2 if i==j else 1 for i in range(d))) for j in range(d)]
        B=[schedules(Q,tuple(3 if i==j else 1 for i in range(d))) for j in range(d)]
        C=sum(schedules(Q,tuple(2 if i in (j,k) else 1 for i in range(d)))
              for j,k in combinations(range(d),2))
        eps=[len(u)%2 for u in modules];chi=[slot(u,2)[2]//4%2 for u in modules]
        assert all(b%2==0 for b in B)
        assert v%8==((paths(Q)[-1]-2*sum(A))*prod(h)+4*(sum(a*c+(b//2)*e for a,b,c,e in zip(A,B,chi,eps))+C))%8
        shadow=tuple((0,) if e else (2,0) for e in eps)
        assert v%3==(prod(h)*paths(compose(Q,shadow))[-1])%3
    print('exact probe outputs',low,high,'; swaps (2270565,2286945); 60 dyadic slot checks')
if __name__=='__main__':main()
