"""Exploratory bounded block-word search, never a theorem checker.

Modes 7, 8 and end change the fixed block model. Surviving sampled templates
need a separate all-order proof. No survivor in this finite template bank
does not mean no detour. The production proof/replay is endpoint_check.py.
"""
from itertools import combinations, permutations
from random import Random
from collections import Counter
import json
import sys
ARG=sys.argv[1] if len(sys.argv)>1 else '7'
MODE=8 if ARG=='end' else int(ARG)
END=ARG=='end'

def need(ok,why):
    if not ok: raise RuntimeError(why)
def pe(q): return frozenset(frozenset((u,v)) for u,v in zip(q,q[1:]))
def path(T,q): return all((u,v) in T for u,v in zip(q,q[1:]))
def blue(T,q):
    r=dict(zip(q,q[::-1]))
    return all((r[v],r[u]) in T for u,v in T)
def move(T,q):
    e=pe(q)
    return frozenset((u,v) if frozenset((u,v)) in e else (v,u) for u,v in T)
def source(n,h,rng,cy=True):
    I={h,n-3,n-2}; J={n-1-v for v in I}; T=set(); seen=set()
    forced={(n-2,1),(n-2,2)}
    if MODE>=8 and not END:forced.add((n-2,3))
    if cy is not None: forced.add((n-3,2) if cy else (2,n-3))
    for u,v in combinations(range(1,n-1),2):
        if (u,v) in seen: continue
        orb={(u,v),(n-1-v,n-1-u)};seen|=orb
        if v==u+1: bit=True
        elif any((a,b) in forced for a,b in orb): bit=True
        elif any((b,a) in forced for a,b in orb): bit=False
        else: bit=bool(rng.getrandbits(1))
        T.update((a,b) if bit else (b,a) for a,b in orb)
    for v in range(1,n-1):
        T.add((v,0) if v in I else (0,v))
        T.add((n-1,v) if v in J else (v,n-1))
    T.add((n-1,0));T=frozenset(T)
    need(blue(T,tuple(range(n))) and path(T,tuple(range(n))), 'source')
    need(forced<=T,'forced guards')
    return T
def word(n,code):
    blocks=[(0,),(1,),(2,),tuple(range(3,n-3)),(n-3,),(n-2,),(n-1,)] if MODE==7 else ([(0,),(1,),(2,),tuple(range(3,n-4)),(n-4,),(n-3,),(n-2,),(n-1,)] if END else [(0,),(1,),(2,),(3,),tuple(range(4,n-3)),(n-3,),(n-2,),(n-1,)])
    out=()
    for c in code:
        b=blocks[abs(c)-1]
        out+=b if c>0 else b[::-1]
    return out
def recover(e,n):
    adj=[[] for _ in range(n)]
    for pair in e:
        u,v=tuple(pair);adj[u].append(v);adj[v].append(u)
    if len(e)!=n-1 or any(len(a) not in (1,2) for a in adj):return None
    ends=[i for i in range(n) if len(adj[i])==1]
    if len(ends)!=2:return None
    out=[ends[0]];prev=-1
    while True:
        nxt=[v for v in adj[out[-1]] if v!=prev]
        if not nxt:break
        prev=out[-1];out.append(nxt[0])
        if len(out)>n:return None
    return tuple(out) if len(out)==n else None
def evaluate(T,qs):
    n=max(max(q) for q in qs)+1;now=T
    for q in qs:
        if not path(now,q) or blue(now,q):return False
        now=move(now,q)
    return now==move(T,tuple(range(n)))
def main():
    n=15;rng=Random(20261025);T=source(n,n-4 if END else 2,rng,True);B=move(T,tuple(range(n)))
    codes=[tuple(-c if c==(4 if END else MODE-3) and rev else c for c in p) for p in permutations(range(1,MODE+1)) for rev in (0,1)]
    lookup={word(n,c):c for c in codes}
    first=[(q,pe(q),c) for q,c in lookup.items() if path(T,q) and not blue(T,q)]
    last=[(q,pe(q),c) for q,c in lookup.items() if path(B,q) and not blue(B,q)]
    print('seven-block search',len(first),len(last),flush=True)
    e=pe(tuple(range(n)));cands=[]
    for q,eq,cq in first:
        S=move(T,q)
        for s,es,cs in last:
            r=recover(e^eq^es,n)
            if r is None:continue
            if not path(S,r):r=r[::-1]
            if path(S,r) and not blue(S,r):
                if r not in lookup:continue
                cands.append((cq,lookup[r],cs))
    print('candidate triples',len(cands),json.dumps(cands[:20]),flush=True)
    for sample in range(40):
        for cy in ((True,) if END else (True,False)):
            n=15+sample%7;T=source(n,n-4 if END else 2,rng,cy)
            cands=[cs for cs in cands if evaluate(T,[word(n,c) for c in cs])]
        if not cands:break
    print('uniform survivors',len(cands),flush=True)
    print(json.dumps(cands[:20]))
if __name__=='__main__':main()
