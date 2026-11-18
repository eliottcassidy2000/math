#!/usr/bin/env python3
"""Finite target decision, least-order witnesses and exact arc loads (stdlib)."""
import argparse
import itertools
import json
import math
from functools import lru_cache

def dp(rows):
    n=len(rows)
    d=[[0]*n for _ in range(1<<n)]
    for v in range(n): d[1<<v][v]=1
    for s in range(1,(1<<n)-1):
        for v,p in enumerate(d[s]):
            if not p: continue
            avail=rows[v]&~s
            while avail:
                bit=avail&-avail
                avail^=bit
                d[s|bit][bit.bit_length()-1]+=p
    return d

def check_rows(rows):
    n=len(rows)
    if not n or any(x<0 or x>=1<<n or x>>i&1 for i,x in enumerate(rows)):
        raise ValueError('outgoing masks must encode a nonempty loopless graph')
    for u,v in itertools.combinations(range(n),2):
        if ((rows[u]>>v&1)+(rows[v]>>u&1))!=1:
            raise ValueError('expected exactly one orientation per unordered pair')

def profile(rows):
    check_rows(rows)
    n=len(rows); f=(1<<n)-1
    out=dp(rows)
    rev=[sum(1<<u for u in range(n) if rows[u]>>v&1) for v in range(n)]
    inc=dp(rev)
    H=sum(out[f]); weights=[[0]*n for _ in range(n)]
    two_cover=0
    for s in range(1,f):
        t=f^s
        two_cover+=sum(out[s])*sum(inc[t])
        for u,p in enumerate(out[s]):
            if not p: continue
            for v,q in enumerate(inc[t]):
                if q: weights[u][v]+=p*q
    arcs=[(u,v,weights[u][v],weights[v][u])
          for u in range(n) for v in range(n) if rows[u]>>v&1]
    assert sum(w for u,v,w,z in arcs)==(n-1)*H
    assert sum(w+z for u,v,w,z in arcs)==two_cover
    assert all((z-w)%2==0 for u,v,w,z in arcs)
    for i in range(n):
        assert sum(w for u,v,w,z in arcs if u==i)==H-out[f][i]
        assert sum(w for u,v,w,z in arcs if v==i)==H-inc[f][i]
    return dict(order=n,H=H,rows=rows,start=inc[f],end=out[f],
                ordered_nonempty_two_covers=two_cover,
                unused=sum(w==0 for u,v,w,z in arcs),
                inert=sum(w==z==0 for u,v,w,z in arcs),
                zero_arcs=[[u,v,z] for u,v,w,z in arcs if w==0],
                arcs_uv_load_reverse_load=arcs)

def connected(rows):
    n=len(rows);f=(1<<n)-1
    def reach(arcs):
        s=1
        while s!=f:
            t=s
            for i,a in enumerate(arcs):
                if s>>i&1:t|=a
            if s==t:return False
            s=t
        return True
    return reach(rows) and reach([sum(1<<u for u,x in enumerate(rows) if x>>v&1)
                                 for v in range(n)])

def busch(n):
    r=n%3
    return ((3,5,9)[r])*5**((n-(3,4,5)[r])//3)

def frame(n):
    """Every n-tournament has a labeling whose path runs 0..n-1."""
    chords=[(u,v) for u in range(n) for v in range(u+2,n)]
    for b in range(1<<len(chords)):
        rows=[1<<(i+1) for i in range(n-1)]+[0]
        for j,(u,v) in enumerate(chords):
            if b>>j&1:rows[v]|=1<<u
            else:rows[u]|=1<<v
        yield tuple(rows)

def canonical(rows):
    n=len(rows)
    blocks=sorted({x.bit_count() for x in rows})
    perms=[list(itertools.permutations([u for u in range(n)
                                      if rows[u].bit_count()==b])) for b in blocks]
    def relabel(p):
        return tuple(sum(1<<v for v,b in enumerate(p) if rows[a]>>b&1)
                     for a in p)
    return min(relabel(sum(x,())) for x in itertools.product(*perms))

def join(A,B):
    n,m=len(A),len(B)
    fill=((1<<m)-1)<<n
    return tuple(x|fill for x in A)+tuple(x<<n for x in B)

def solve(N, all_classes=False, cap=None):
    if N<1 or N%2==0: return []
    if N==1: return [(0,)]
    if N in (7,21): return []
    @lru_cache(None)
    def rec(h,lim):
        best=lim+1
        sources=set()
        def admit(rows):
            nonlocal best
            n=len(rows)
            if n>best or n>lim:return
            if n<best:best=n;sources.clear()
            sources.add(canonical(rows) if all_classes else rows)

        # Split into a first nontrivial SCC and what follows, then minimize.
        for d in range(3,math.isqrt(h)+1,2):
            if h%d:continue
            e=h//d
            for x in rec(d,lim-3):
                for y in rec(e,lim-len(x)):
                    admit(join(x,y));admit(join(y,x))

        # Every strongly connected birth has order bounded by Busch's minimum.
        for n in itertools.count(3):
            if n>min(best,lim) or busch(n)>h:break
            if math.factorial(n)<h:continue
            for rows in frame(n):
                if connected(rows) and sum(dp(rows)[-1])==h:
                    admit(rows)
                    if not all_classes:break
        return tuple(sorted(sources)) if all_classes else tuple(sorted(sources)[:1])
    # Each nontrivial strong component consumes at least 3^(n/3) paths.
    # Integer-only normalization cap: floor(log_3(N^3)).
    t=N**3; bound=0
    while t>=3:t//=3;bound+=1
    return list(rec(N,bound if cap is None else min(bound,cap)))

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('target',nargs='?',type=int)
    p.add_argument('--rows',help='comma-separated outgoing masks; inspect a supplied witness')
    p.add_argument('--all-classes',action='store_true')
    p.add_argument('--max-order',type=int)
    a=p.parse_args()
    if a.rows:
        print(json.dumps(profile(list(map(int,a.rows.split(',')))),indent=2))
    elif a.target is not None:
        witnesses=solve(a.target,a.all_classes,a.max_order)
        print(json.dumps(dict(target=a.target,status='found' if witnesses else
             ('globally absent' if a.max_order is None else 'absent through order cap'),
             witnesses=[profile(w) for w in witnesses]),indent=2))
    else:p.error('supply a target or --rows')
