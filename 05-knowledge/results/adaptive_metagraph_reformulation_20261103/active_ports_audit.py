#!/usr/bin/env python3
"""Independent endpoint-state lumpability by explicit successor-edge covers.
No producer imports, no permutation/cut cover enumerator, no subset DP.
"""
from itertools import product,permutations
from collections import Counter,defaultdict,deque
from functools import lru_cache
import json


def need(p,m):
    if not p:raise RuntimeError(m)


def birth(T,r):
    n=len(T)
    return tuple(row|((1<<n) if i==r else 0) for i,row in enumerate(T))+(((1<<n)-1)^(1<<r),)


def ports(T,d):return tuple(v for v,r in enumerate(T) if len(T)-1-r.bit_count()<=d)


@lru_cache(None)
def covers(T):
    """Choose directed successor edges, each at most once in/out, forbid cycles."""
    n=len(T);out={}
    options=[(-1,)+tuple(j for j in range(n) if row>>j&1) for row in T]
    for succ in product(*options):
        incoming=[v for v in succ if v!=-1]
        if len(set(incoming))!=len(incoming):continue
        starts=[v for v in range(n) if v not in incoming]
        seen=set();pairs=[]
        for s in starts:
            v=s
            while True:
                if v in seen:break
                seen.add(v)
                if succ[v]==-1:pairs.append((s,v));break
                v=succ[v]
        if len(seen)!=n:continue
        out[succ]=tuple(sorted(pairs))
    return out


def signature(pairs,A):
    pos={v:i for i,v in enumerate(A)}
    s=tuple(sorted((pos.get(a,-1),pos.get(b,-1)) for a,b in pairs))
    return None if sum(b==-1 for _,b in s)>1 else s


def direct_H(T):
    return sum(all(T[p[i]]>>p[i+1]&1 for i in range(len(T)-1)) for p in permutations(range(len(T))))


def check_step(T,r,d):
    U=birth(T,r);old=covers(T);new=covers(U);A=ports(T,d);B=ports(U,d)
    n=len(T);lifts=defaultdict(Counter)
    for succ,pairs in new.items():
        sig=signature(pairs,B)
        if sig is None:continue
        parent=tuple(-1 if e==n else e for e in succ[:n])
        need(parent in old,'deleting newborn gives a literal parent path cover')
        need(signature(old[parent],A) is not None,'two forgotten ends cannot be rescued')
        lifts[parent][sig]+=1
    bysig={};max_mult=0;positive=0;forgotten=0
    for edges,pairs in old.items():
        sig=signature(pairs,A)
        if sig is None:
            need(not lifts[edges],'discarded parent permanently dead');continue
        terms=lifts[edges]
        if sig in bysig:need(bysig[sig]==terms,'same endpoint signature has identical full lift kernel')
        else:bysig[sig]=terms
        max_mult=max(max_mult,max(terms.values(),default=0))
        positive+=bool(terms);forgotten+=any(a==-1 or b==-1 for a,b in sig)
    # Covers are edge sets, so one component means exactly one Hamilton word.
    need(sum(len(p)==1 for p in new.values())==direct_H(U),'independent edge-cover / word count')
    need(len(B)<=2*d,'actual temporal boundary')
    return {'source_covers':len(old),'target_covers':len(new),'signature_classes':len(bysig),
            'positive_parent_covers':positive,'forgotten_endpoint_parent_covers':forgotten,
            'largest_transition_multiplicity':max_mult}


def control_closure(d):
    T=birth((2,4,1),0);seed=tuple(3-T[v].bit_count() for v in ports(T,d))
    seen={seed};q=deque([seed])
    while q:
        a=q.popleft()
        for j in range(len(a)):
            b=tuple(x if i==j else x+1 for i,x in enumerate(a))
            b=tuple(x for x in b if x<=d)+(1,)
            if b not in seen:seen.add(b);q.append(b)
    return len(seen)


def main():
    result={};seed=birth((2,4,1),0)
    for d in (1,2,3):
        todo=[seed];steps=[]
        while todo:
            T=todo.pop()
            if len(T)>=6:continue
            for r in ports(T,d):
                steps.append(check_step(T,r,d));todo.append(birth(T,r))
        count=control_closure(d);need(count=={1:1,2:6,3:24}[d],'complete finite degree control')
        result[str(d)]={'named_birth_edges':len(steps),'degree_controls':count,
            'source_cover_observations':sum(s['source_covers'] for s in steps),
            'forgotten_endpoint_parent_observations':sum(s['forgotten_endpoint_parent_covers'] for s in steps),
            'largest_transition_multiplicity':max(s['largest_transition_multiplicity'] for s in steps)}
    T=(2,4,1)
    for r in (0,0,0):T=birth(T,r)
    extra=check_step(T,0,3)
    need(extra['forgotten_endpoint_parent_covers']>0,'cap-three forgotten-port hostile')
    return {'status':'PASS independent cover-by-cover endpoint lumpability', 'controls':result,
            'cap3_extra_forgotten_port_step':extra,
            'scope':'All actual cap1/2/3 histories from normalized T4 through order6; all path-cover lifts per birth edge. Exact support dimensions15/460/9672 remain producer closure receipts, not independently recomputed here.'}


if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
