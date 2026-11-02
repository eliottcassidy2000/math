#!/usr/bin/env python3
"""Independent SCC reduction consumer and exact small-height dyadic face banks."""
from itertools import permutations
import json
from check import tournament,state,pos,HP,E,chain,need,induced

def strong(A):
    n=len(A)
    if n<=1:return True
    for root in range(n):
        seen=1<<root
        while True:
            new=seen
            for i in range(n):
                if seen>>i&1:new|=A[i]
            if new==seen:break
            seen=new
        if seen!=(1<<n)-1:return False
    return True

def join(A,B):
    n=len(A);m=len(B);bits=((1<<m)-1)<<n
    return tuple(x|bits for x in A)+tuple(x<<n for x in B)

def canonical_bases(h):
    sources={}
    if h==1:
        A=chain(2);sources[1]=(A,1)
    else:
        # h in {3,5} is prime, so its nonsingleton core is one strong block.
        # Every such block has order<=h by its spanning Hamilton cycle.
        for n in range(3,h+1):
            for code in range(1<<(n*(n-1)//2)):
                B=tournament(n,code)
                if not strong(B) or HP(B)!=h:continue
                candidates=[]
                for v in range(n):candidates.append((join(B,(0,)),1<<v))
                candidates.append((join((0,),B),1))
                candidates.append((join(B,chain(2)),1<<n))
                for A,cut in candidates:
                    C,L,g,t=pos(state(A,cut));need(L==g==0 and t==h,'invalid base recipe')
                    sources.setdefault(C,(A,cut))
    # Keep the smallest representative for each odd part.
    best={}
    for b in sorted(sources):
        odd=b
        while odd%2==0:odd//=2
        best.setdefault(odd,b)
    return {b:sources[b] for b in sorted(best.values())}

def sccs(A):
    n=len(A);reach=[]
    for i in range(n):
        R=1<<i
        while True:
            Q=R
            for v in range(n):
                if R>>v&1:Q|=A[v]
            if Q==R:break
            R=Q
        reach.append(R)
    todo=set(range(n));parts=[]
    while todo:
        v=min(todo);part={u for u in todo if reach[v]>>u&1 and reach[u]>>v&1};parts.append(part);todo-=part
    # Strict dominance is a total order on SCCs.
    snapshot=tuple(parts)
    parts.sort(key=lambda P:-sum(A[next(iter(P))]>>next(iter(Q))&1 for Q in snapshot if P is not Q))
    return parts

def reduce_face(A,v):
    parts=sccs(A);i=next(i for i,P in enumerate(parts) if v in P)
    need(i<len(parts)-1,'selected vertex terminal')
    keep=set().union(*(P for P in parts if len(P)>1));keep.add(v)
    if not any(len(P)>1 for P in parts[i+1:]):keep.add(next(iter(parts[-1])))
    suffix_deleted=sum(len(P)==1 and not(P<=keep) for P in parts[i+1:])
    labels=sorted(keep);B=induced(A,labels);cut=1<<labels.index(v)
    return B,cut,suffix_deleted

def main():
    result={};total_core=0
    for h,expected in [(1,[1]),(3,[3]),(5,[5,6,7,9,16])]:
        B=canonical_bases(h);need(list(B)==expected,'small-height basis changed')
        entries=[]
        for b,(A,cut) in B.items():
            entries.append({'C_base':b,'parent':A,'cut':cut,'native_state':state(A,cut)})
        result[str(h)]=entries
    # Complete independent labelled parent universe through order five; all singleton
    # cuts with g=0, and exact SCC deletion gives the dyadic factor in every case.
    tests=0;small_membership=0
    banks={h:canonical_bases(h) for h in (1,3,5)}
    for n in range(2,6):
        for code in range(1<<(n*(n-1)//2)):
            A=tournament(n,code)
            for v in range(n):
                C,L,g,h=pos(state(A,1<<v))
                need(L==0,'singleton interior weight')
                if g:continue
                B,cut,k=reduce_face(A,v);c,l,G,t=pos(state(B,cut))
                need((l,G,t)==(0,0,h) and C==c*2**k,'exact SCC face reduction');tests+=1
                if h in banks:
                    need(any(C%b==0 and (C//b)&((C//b)-1)==0 for b in banks[h]),'small height face membership');small_membership+=1
    print(json.dumps({'status':'PASS','exact_all_order_small_height_bases':result,
       'complete_labelled_parent_orders':[2,3,4,5],'literal_SCC_reductions':tests,
       'small_height_membership_controls':small_membership,
       'scope':'All-order conclusion uses proved finite-core reduction; controls enumerate the stated exact finite universes.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
