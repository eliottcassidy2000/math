#!/usr/bin/env python3
"""Independent exhaustive p=5 direction census and p=2,3 frame census.
Uses coordinate tuples for elementary abelian additive groups, never Z/p^n.
"""
import itertools, collections, json, math

def directions(s,p):
    out=set()
    for a,b in itertools.combinations(s,2):
        x,y=((a[i]-b[i])%p for i in range(2))
        out.add(p if x==0 else y*pow(x,-1,p)%p)
    return out

def census(p):
    pts=list(itertools.product(range(p),repeat=2))
    hist=collections.Counter(len(directions(s,p)) for s in itertools.combinations(pts,p))
    assert sum(hist.values())==math.comb(p*p,p)
    assert hist[1]==p*(p+1)
    assert all(k==1 or k>=(p+3)//2 for k in hist)
    return dict(sorted(hist.items()))

def frame_census(p,n,c):
    pts=list(itertools.product(range(p),repeat=n))
    add=lambda a,b:tuple((x+y)%p for x,y in zip(a,b))
    sub=lambda a,b:tuple((x-y)%p for x,y in zip(a,b))
    A={add(c,tuple(x)+(0,)*(n-2)) for x in itertools.product(range(p),repeat=2)}
    # Completeness: one factor lies in a coset of U. Enumerate all p-subsets
    # within every such coset (a necessary consequence of fixed-plane sums).
    subs=[]
    for tail in itertools.product(range(p),repeat=n-2):
        plane=[tuple(x)+tail for x in itertools.product(range(p),repeat=2)]
        subs.extend(frozenset(s) for s in itertools.combinations(plane,p))
    cnt=disj=both=nonline=0
    for B,D in itertools.combinations(subs,2):
        if {add(b,d) for b in B for d in D}!=A: continue
        cnt+=1; disj+=not(B&D)
        db=directions([b[:2] for b in B],p)
        dd=directions([d[:2] for d in D],p)
        assert db.isdisjoint(dd)
        both+=(len(db)==len(dd)==1)
        nonline+=(len(db)>1 and len(dd)>1)
        assert len({sub(b,d) for b in B for d in D})==p*p
    q=p**n
    cp=(p+1)*p**(p-1)-p*(p+1)//2
    assert cnt==cp*q and both==p*(p+1)*q//2 and nonline==0
    expected=cp*(q-p*p) if p%2 else (0 if (0,)*n in A else cp*q)
    assert disj==expected
    return dict(p=p,n=n,C=c,unordered_frames=cnt,disjoint=disj,both_line=both)

out={"direction_census":{str(p):census(p) for p in (2,3,5)},"frame_census":[frame_census(p,n,c) for p,n,c in [(2,2,(0,0)),(2,3,(0,0,0)),(2,3,(0,0,1)),(3,2,(0,0)),(3,3,(0,0,0)),(3,3,(0,0,1))]]}
print(json.dumps(out,indent=2))
