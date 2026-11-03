#!/usr/bin/env python3
"""Independent all-labelled strong n<=5 peel audit and two primitive-core orbits.
No producer imports; Burnside controls enumerate root words, not tournaments.
"""
from itertools import combinations, permutations, product
from collections import Counter
from math import prod
import json


def need(p,m):
    if not p:raise RuntimeError(m)


def decode(n,code):
    R=[0]*n
    for k,(a,b) in enumerate(combinations(range(n),2)):
        u,v=(a,b) if code>>k&1 else (b,a);R[u]|=1<<v
    return tuple(R)


def strong(T,V):
    if len(V)<3:return False
    allowed=sum(1<<v for v in V)
    for start in V:
        seen=1<<start
        while True:
            reach=seen
            for v in V:
                if seen>>v&1:reach|=T[v]&allowed
            if reach==seen:break
            seen=reach
        if seen!=allowed:return False
    return True


def automorphisms(T,V=None):
    if V is None:V=tuple(range(len(T)))
    return tuple(dict(zip(V,p)) for p in permutations(V)
        if all(((T[a]>>b)&1)==((T[x]>>y)&1) for (a,b),(x,y) in zip(combinations(V,2), combinations(p,2))))


def peel(T):
    V=tuple(range(len(T)));reverse=[]
    while len(V)>3:
        candidates=[]
        for x in V:
            parent=tuple(v for v in V if v!=x)
            incoming=[v for v in parent if T[v]>>x&1]
            if len(incoming)==1 and strong(T,parent):candidates.append((x,incoming[0]))
        need(len(candidates)<=1,'unique legal inverse')
        if not candidates:break
        x,r=candidates[0];reverse.append((x,r));V=tuple(v for v in V if v!=x)
    return V,tuple(reversed(reverse))


def add(T,r):
    n=len(T)
    return tuple(row|((1<<n) if i==r else 0) for i,row in enumerate(T))+(((1<<n)-1)^(1<<r),)


def code(T,p):
    return sum(((T[p[i]]>>p[j])&1)<<k for k,(i,j) in enumerate(combinations(range(len(T)),2)))


def canonical(T):
    return min(code(T,p) for p in permutations(range(len(T))))


def main():
    census=[];audited=0
    for n in range(3,6):
        total=primitive=0
        for mask in range(1<<(n*(n-1)//2)):
            T=decode(n,mask)
            if not strong(T,tuple(range(n))):continue
            total+=1;K,history=peel(T);primitive+=not history
            core_roots={r for _,r in history if r in K}
            actual={tuple(g[i] for i in range(n)) for g in automorphisms(T)}
            predicted=set()
            for g in automorphisms(T,K):
                if any(g[r]!=r for r in core_roots):continue
                predicted.add(tuple(g[i] if i in g else i for i in range(n)))
            need(actual==predicted,'exact residual automorphism group')
            V=set(K)
            for x,r in history:
                need(r in V and all(bool(T[v]>>x&1)==(v==r) for v in V),'named history replay')
                V.add(x)
            need(V==set(range(n)),'complete names')
            audited+=1
        census.append({'n':n,'labelled_strong':total,'irreducible_labelled':primitive})
    need([(x['labelled_strong'],x['irreducible_labelled']) for x in census]==[(2,2),(24,0),(544,64)],'small labelled controls')
    R5=tuple(sum(1<<j for j in range(5) if (j-i)%5 in (1,2)) for i in range(5))
    B=(2|8,4|8,1|8,16,1|2|4) # C3 -> b -> c -> C3
    records=[]
    for name,K,expected in (('regular5',R5,[1,6]),('C3_C3_1_1',B,[3,14])):
        need(strong(K,tuple(range(5))) and not peel(K)[1],'actual irreducible core')
        G=automorphisms(K);fix=[sum(g[i]==i for i in range(5)) for g in G]
        depths=[]
        for m in (1,2):
            seen={};word_counts=Counter()
            for word in product(*(range(5+i) for i in range(m))):
                T=K
                for r in word:T=add(T,r)
                c=canonical(T);seen.setdefault(c,T);word_counts[c]+=1
            numerator=sum(prod(f+i for i in range(m)) for f in fix)
            need(numerator%len(G)==0 and numerator//len(G)==len(seen)==expected[m-1],'Burnside orbit count')
            for c,T in seen.items():
                a=len(automorphisms(T))
                need(word_counts[c]*a==len(G),'source word orbit/stabilizer weight')
            depths.append({'depth':m,'root_words':sum(word_counts.values()),'ordinary_classes':len(seen),
                           'automorphism_order_histogram':dict(Counter(len(automorphisms(T)) for T in seen.values()))})
        records.append({'core':name,'rows':K,'automorphism_order':len(G),'fixed_points_by_automorphism':fix,'depths':depths})
    converse_controls=[]
    for depth in (1,2,3):
        exits=[]
        for word in product(*(range(3+i) for i in range(depth))):
            T=(2,4,1)
            for r in word:T=add(T,r)
            U=tuple(((1<<len(T))-1)^row^(1<<i) for i,row in enumerate(T))
            if len(peel(U)[0])>3:exits.append((word,T,U))
        need(bool(exits)==(depth==3),'first adaptive converse exit')
        converse_controls.append({'order':depth+3,'named_histories':prod(range(3,3+depth)),
                                  'converse_exit_histories':len(exits)})
    word=(0,1,2);T=(2,4,1)
    for r in word:T=add(T,r)
    U=tuple(((1<<len(T))-1)^row^(1<<i) for i,row in enumerate(T))
    need(T==(10,20,33,6,13,27) and U==(52,41,26,49,34,4),'literal converse hostile')
    need(len(peel(T)[0])==3 and len(peel(U)[0])==6,'ancestry does not descend through converse merge')
    return {'status':'PASS','all_labelled_strong_controls':census,'sources_audited':audited,'primitive_core_controls':records,
            'converse_controls':converse_controls,'converse_hostile':{'history':word,'source_rows':T,'converse_rows':U,'core_orders':[3,6]},
            'scope':'Canonical peeling, exact residual automorphism groups, literal named root-word realizations; not a primitive-core inventory at all orders.'}


if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
