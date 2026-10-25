"""Explore sufficient affine contradiction chains for fixed ledgers.

No completeness or shortest-chain claim: a missing chain does not establish
satisfiability. The final proof uses direct witnesses, and endpoint_check.py
independently checks those witnesses rather than trusting this exploration.
"""
from itertools import combinations
from collections import deque

def canonical(u,v):return ((u,v),0) if u<v else ((v,u),1)
def edges(q):return {tuple(sorted((u,v))) for u,v in zip(q,q[1:])}
def words(n):
    a,b,x,d,y,c,z=0,1,2,3,n-3,n-2,n-1;U=tuple(range(4,n-3))
    return ((b,x,a,d)+U+(y,c,z),(b,x,a,d,c,z)+U+(y,),(a,b,x,d,c,z)+U+(y,))
def contradiction(n,stage):
    P=tuple(range(n));qs=words(n);q=qs[stage];prior=set()
    for prev in qs[:stage]:prior^=edges(prev)
    rp=dict(zip(P,P[::-1]));rq=dict(zip(q,q[::-1]));pairs=list(combinations(range(n),2));adj={p:[] for p in pairs}
    for pair in pairs:
        u,v=pair
        for r,label in ((rp,'P'),(rq,'Q')):
            other,rev=canonical(r[v],r[u]); parity=rev
            if label=='Q':parity^=(pair in prior)^(other in prior)
            adj[pair].append((other,parity,label));adj[other].append((pair,parity,label))
    known={}
    def setarc(u,v,why):
        pair,rev=canonical(u,v);val=1^rev
        if pair in known and known[pair][0]!=val:raise RuntimeError('input conflict')
        known[pair]=(val,why)
    for u,v in zip(P,P[1:]):setarc(u,v,'path')
    I={2,n-3,n-2};J={n-1-i for i in I}
    for v in range(1,n-1):
        setarc(v,0,'port') if v in I else setarc(0,v,'port')
        setarc(n-1,v,'port') if v in J else setarc(v,n-1,'port')
    setarc(n-1,0,'closed');setarc(n-2,1,'closed middle');setarc(n-2,2,'chord');setarc(n-2,3,'next chord')
    best=None
    for start,(sv,sw) in known.items():
        seen={start:(sv,[])};que=deque([start])
        while que:
            p=que.popleft();val,chain=seen[p]
            if p in known and known[p][0]!=val:
                if best is None or len(chain)<len(best[-1]):best=(start,sv,sw,p,known[p],chain)
                break
            if best is not None and len(chain)>=len(best[-1]):continue
            for v,d,l in adj[p]:
                if v not in seen:seen[v]=(val^d,chain+[(l,p,v,d)]);que.append(v)
    return best
if __name__=='__main__':
    for n in range(15,31):
        print(n,contradiction(n,1))
