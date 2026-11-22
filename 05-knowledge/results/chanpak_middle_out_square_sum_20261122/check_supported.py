"""Exact outward paired-cut search, specializing quickly on supported targets."""
from collections import deque
from itertools import combinations
from random import Random
from pathlib import Path
import sys

# The two earlier proof packets use the same row numbering and literal checker.
sys.path.insert(0,str(Path(__file__).resolve().parent.parent /
                     'chanpak_mirimanoff_unit_lifts_20261122'))
try:
    from check_seams import flipped, subpaths
except ImportError:  # standalone companion for replay outside the repository
    from fixture_seams import flipped, subpaths

def capacities(rows):
    n=len(rows);N=1<<n;h=subpaths(rows);h[0]=1
    W=[[0]*n for _ in rows]
    for u in range(n):
        D=[[0]*n for _ in range(N)]
        D[1<<u][u]=1
        for s in range(1<<u,N):
            for v,x in enumerate(D[s]):
                if not x:continue
                if s.bit_count()%2==0: W[u][v]+=x*h[(N-1)^s]
                use=rows[v]&~s
                while use:
                    k=use&-use;use^=k
                    D[s|k][k.bit_length()-1]+=x
    return W

def value(W,S):
    return sum(w for i,row in enumerate(W) if S>>i&1
               for j,w in enumerate(row) if not S>>j&1)

def mincut(W,inside,outside):
    n=len(W);src,sink=n,n+1;B=sum(map(sum,W))+1
    R=[r.copy()+[0,0] for r in W]+[[0]*(n+2) for _ in range(2)]
    for i in range(n):
        if inside>>i&1:R[src][i]=B
        if outside>>i&1:R[i][sink]=B
    m=0
    while True:
        prev=[-1]*(n+2);prev[src]=src;dq=deque([src])
        while dq and prev[sink]<0:
            a=dq.popleft()
            for b,c in enumerate(R[a]):
                if c and prev[b]<0:
                    prev[b]=a;dq.append(b)
        if prev[sink]<0:
            return m,sum(1<<i for i in range(n) if prev[i]>=0)
        b=sink;cap=B
        while b!=src:
            a=prev[b];cap=min(cap,R[a][b]);b=a
        b=sink
        while b!=src:
            a=prev[b];R[a][b]-=cap;R[b][a]+=cap;b=a
        m+=cap

def minimum_chain(W,inside,outside,minimum,first):
    n=len(W);S=first
    while True:
        yield S
        new=[]
        for i in range(n):
            if (S|outside)>>i&1: continue
            val,U=mincut(W,S|1<<i,outside)
            if val==minimum:new.append(U)
        if not new:return
        S=min(new,key=int.bit_count)

def signed_bound(W,inside,outside):
    n=len(W);lo=hi=0
    def bits(i):
        return [1] if inside>>i&1 else [0] if outside>>i&1 else [0,1]
    for i,j in combinations(range(n),2):
        poss=[W[i][j]*a*(1-b)+W[j][i]*b*(1-a)
              for a in bits(i) for b in bits(j)]
        lo+=min(poss);hi+=max(poss)
    return lo,hi

def solve(W0,W1,target,order,inside=0,outside=0,stats=None):
    n=len(W0)
    assert sorted(order)==list(range(n))
    if stats is None:stats={}
    def hit(k):stats[k]=stats.get(k,0)+1
    def weighted(a,b):
        return [[a*W0[i][j]+b*W1[i][j] for j in range(n)] for i in range(n)]
    projections=[(weighted(a,b),a*target[0]+b*target[1]) for a,b in
                 ((1,0),(0,1),(8,9),(9,8),(1,1))]
    D=weighted(-1,1)
    delta=target[1]-target[0]
    def branch(I,O):
        hit('branch')
        if I&O:return
        lo,hi=signed_bound(D,I,O)
        if not lo<=delta<=hi:hit('signed');return
        for W,g in projections:
            m,S=mincut(W,I,O)
            if m>g:hit('mincut');return
            if m==g:
                hit('chain')
                for A in minimum_chain(W,I,O,m,S):
                    if (value(W0,A),value(W1,A))==target:return A
                return
        for i in order:
            if not (I|O)>>i&1:
                a=branch(I|1<<i,O)
                return a if a is not None else branch(I,O|1<<i)
    return branch(inside,outside)

def graph_with_contact(rows,S):
    n=len(rows);return tuple(r|((not(S>>i&1))<<n) for i,r in enumerate(rows))+(S,)

def verify(a,ij,tests):
    n=len(a);b=flipped(a,*ij)
    old=(subpaths(a)[-1],subpaths(b)[-1]);assert old[1]==old[0]+2
    W0,W1=capacities(a),capacities(b)
    data=[(value(W0,S),value(W1,S)) for S in range(1<<n)]
    for S,(u,v) in enumerate(data):
        assert (subpaths(graph_with_contact(a,S))[-1],
                subpaths(graph_with_contact(b,S))[-1])==(old[0]+2*u,old[1]+2*v)
    for k,I,O,expected in tests:
        stats={}
        result=solve(W0,W1,(k,k),list(ij)+[i for i in range(n) if i not in ij],I,O,stats)
        truth=[s for s,t in enumerate(data) if t==(k,k) and s&I==I and not s&O]
        assert (result is None)==(not truth)
        if expected is not None:assert result==expected
        print(old,k,I,O,result,stats)
    rng=Random(85)
    for _ in range(100):
        s=rng.randrange(1<<n);p=rng.randrange(1<<n);I=s&p;O=p&~s
        t=data[s] if rng.randrange(2) else (rng.randrange(15),rng.randrange(15))
        order=list(range(n));rng.shuffle(order)
        r=solve(W0,W1,t,order,I,O)
        assert (r is not None)==any(x==t and j&I==I and not j&O for j,x in enumerate(data))
        if r is not None:assert data[r]==t and r&I==I and not r&O

if __name__=='__main__':
    verify((16,1,3,7,14),(0,3),[(4,1,16,15)])
    verify((4,1,2,7),(0,3),[(2,1,2,None)])
    verify((4,1,2,7,15),(0,3),[(6,0,0,None)])
    print('literal weighted cuts agree with direct DP over all contact sets')
