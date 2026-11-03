#!/usr/bin/env python3
"""Independent literal Hamilton-word / broken-cover audit; no producer imports."""
from itertools import combinations, permutations
from collections import Counter
import json


def need(p,m):
    if not p: raise RuntimeError(m)


def decode(n,code):
    R=[0]*n
    for k,(a,b) in enumerate(combinations(range(n),2)):
        u,v=(a,b) if code>>k&1 else (b,a);R[u]|=1<<v
    return tuple(R)


def profile(T):
    n=len(T);h=0;s=[0]*n;F=[0]*n
    for p in permutations(range(n)):
        bad=[i for i in range(n-1) if not (T[p[i]]>>p[i+1]&1)]
        if not bad:
            h+=1;s[p[0]]+=1
            for v in p:F[v]+=1
        elif len(bad)==1:
            F[p[bad[0]]]+=1
    q=[sum(all(T[p[i]]>>p[i+1]&1 for i in range(n-2))
           for p in permutations([w for w in range(n) if w!=v])) for v in range(n)]
    return tuple((h,s[v],q[v],F[v]) for v in range(n))


def birth(T,v):
    n=len(T)
    return tuple(r|((1<<n) if i==v else 0) for i,r in enumerate(T))+(((1<<n)-1)^(1<<v),)


def R(z):
    h,s,q,F=z;return h+F-s,q,q,2*F-q


def N(z):
    h,s,q,F=z;return h+F-s,h-s,h,h+F


def delta(z):
    h,s,q,F=z;return F-h-q+s


def main():
    checked=0
    for n in range(1,5):
        for c in range(1<<(n*(n-1)//2)):
            T=decode(n,c)
            for v,z in enumerate(profile(T)):
                Z=profile(birth(T,v))
                need(Z[v]==R(z),'literal keep-old-root')
                need(Z[n]==N(z),'literal choose-newborn')
                need(delta(z)>=0 and delta(R(z))==delta(z) and delta(N(z))==0,'broken-cover defect')
                d=n-1-T[v].bit_count()
                need(delta(z)>=d*(d-1) and (delta(z)==0)==(d<=1),'sharp indegree detector')
                h,s,q,F=z
                need(R(N(z))==(h+2*F,h,h,h+2*F),'balanced reset')
                checked+=1
    A=(40,17,3,6,13,30);B=(34,24,3,5,13,30);z=(37,5,13,71)
    need(profile(A)[2]==profile(B)[3]==z,'literal same-profile witness')
    receipts=[]
    for T,v,expected in ((A,2,[287,321,219,329,275,233,179]),(B,3,[289,307,339,219,271,227,179])):
        U=birth(T,v)
        vals=[R(p)[0] for p in profile(U)]
        need(vals==expected,'literal arbitrary-root continuation')
        # Recompute actual size-eight children for each root; no profile formula.
        actual=[profile(birth(U,w))[0][0] for w in range(7)]
        need(actual==expected,'actual size-eight Hamilton words')
        receipts.append({'source_rows':T,'root':v,'profile':z,'next_H':actual})
    return {'status':'PASS independent literal permutation audit','all_labelled_n1_to4_root_checks':checked,
        'profile_nonclosure':receipts,'size8_literal_children':14,
        'boundary':'Minimal order six rests on producer complete ordinary augmentation; this independent audit validates the exact witnesses, not that census.'}


if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
