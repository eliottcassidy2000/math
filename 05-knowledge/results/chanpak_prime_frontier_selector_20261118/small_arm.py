#!/usr/bin/env python3
"""Literal small-arm F1 witness, profiled without exponential path DP."""
from collections import Counter
from functools import lru_cache
import argparse
import json
from selector import profile

def g(v):
    if v<2:return 1
    return g(v>>1)+(v>>1)+(v>>1&1) if not v&1 else g(v>>1)

def witness(e):
    assert e>=1
    v=e-1; words=bin(v)[2:][::-1] if v else ''; n=len(words)
    # y is right apex n, a is n+1 and x is n+2.
    y,a,x=n,n+1,n+2
    R=[0]*(n+3)
    for i,c in enumerate(words):
        for j in range(i): R[i]|=1<<j
        if c=='0': R[i]|=1<<y
        else:R[y]|=1<<i
    for i in range(n+1):
        R[i]|=1<<a
        if i!=y:R[i]|=1<<x
    R[a]|=1<<x;R[x]|=1<<y
    return R,words

def analyze(e):
    R,word=witness(e);n=len(word);y,a,x=n,n+1,n+2
    ones={i for i,b in enumerate(word) if b=='1'}
    def scan(kind):
        """DAG of two ordered subsequences in ORIGINAL spine orientation.
        Count completions first, then charge each actual reversed arc."""
        def hops(i,last,b):
            u,v=last; before=u if b==0 else v
            z=(i,v) if b==0 else (u,i)
            es=[(i,before)] if before>=0 else []
            if before<0:
                if kind==0:
                    if b==1:
                        if i in ones:return None
                        es.append((i,y))
                    if b==0:es.append((i,a))
                else:
                    if b==0:es.extend(((i,a),(i,x)))
                    else:es.append((i,a))
            return z,es
        def done(last):
            u,v=last
            if kind==0:
                if u>=0 and u not in ones:return None
                return [(y,u) if u>=0 else (y,a),(a,x)]
            if v>=0 and v not in ones:return None
            return [(x,y),(x,y),(a,x)]+([(y,v)]*2 if v>=0 else [(y,a)])
        @lru_cache(None)
        def back(i,last):
            if i==n:return int(done(last) is not None)
            return sum(back(i+1,z[0]) for b in (0,1) if (z:=hops(i,last,b)) is not None)
        loads=Counter(); front={(-1,-1):1}
        for i in range(n):
            nxt=Counter()
            for last,p in front.items():
                for b in (0,1):
                    z=hops(i,last,b)
                    if z is None:continue
                    state,es=z;q=back(i+1,state)
                    if not q:continue
                    nxt[state]+=p
                    for arc in es:loads[arc]+=p*q*(2 if kind and last[b]>=0 else 1)
            front=nxt
        for last,p in front.items():
            es=done(last)
            if es is not None:
                for arc in es:loads[arc]+=p
        return back(0,(-1,-1)), loads
    g0,L=scan(0);e0,M=scan(1);L.update(M)
    H=g0+2*e0
    assert (g0,e0)==(g(e-1),e)
    arcs=[(u,v,L[u,v]) for u in range(n+3) for v in range(n+3) if R[u]>>v&1]
    assert all(R[u]>>v&1 for u,v in L) and sum(w for u,v,w in arcs)==(n+2)*H
    return dict(order=n+3,H=H,e=e,word_lsb_first=word,rows=R,unused=sum(w==0 for u,v,w in arcs),
                zero_arcs=[[u,v] for u,v,w in arcs if not w],
                arcs_uv_load=arcs)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('e',nargs='*',type=int);p.add_argument('--check',action='store_true')
    a=p.parse_args()
    if a.check:
        for e in range(1,129):
            P=analyze(e);Q=profile(P['rows'])
            assert P['H']==Q['H']
            assert P['arcs_uv_load']==[arc[:3] for arc in Q['arcs_uv_load_reverse_load']], (e,P['zero_arcs'])
        print('128 small-arm full load tables agree with ambient subset DP')
    for e in a.e:print(json.dumps(analyze(e)))
