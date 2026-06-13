#!/usr/bin/env python3
"""TRANSLATOR: all-short interval-cover assignments -> endpoint-owner congruences.
For S = S' u {v=nw}, n|v, the config is tight (G(S') subset D_v) iff every component
(a,b) of G(S') fits inside one v-arc center j/(nw) radius 1/(n^2 w). Endpoints carry
an OWNER (u,k,eps): left a=(k n+1)/(n u) [eps=+1], right b=(k n-1)/(n u) [eps=-1].
Cover-by-arc-j  <=>  |w(k n + eps) - j u| < u/n  (clear denominators by n w u).
For owner u<n this is the EXACT congruence w(k n+eps) = j u  (u | w(k n+eps)); matching
j across the two endpoints cancels w: u_b(k_a n+1) = u_a(k_b n-1).
We VERIFY: tight <=> every interval's two endpoint-owner congruence windows intersect.
opus-2026-06-03-S574."""
from fractions import Fraction as F
from math import gcd, ceil, floor
import random
def dist(x): x%=1; return min(x,1-x)

def G_components(Sp,n):
    """exact open components of {t: ||u t||>1/n for all u in Sp}, each with owner data."""
    THR=F(1,n); pts={}
    # boundary points with (owner u, k, eps); eps=+1 left edge, -1 right edge
    for u in Sp:
        # u-safe left edges (kn+1)/(nu), right edges (kn-1)/(nu)
        for k in range(u+1):
            for eps in (1,-1):
                t=F(k*n+eps,n*u)%1
                pts.setdefault(t,[]).append((u,k,eps))
    order=sorted(pts)
    comps=[]; L=len(order)
    for i in range(L):
        a=order[i]; b=order[(i+1)%L]; ln=(b-a) if b>a else (b-a+1); mid=(a+ln/2)%1
        if all(dist(u*mid)>THR for u in Sp):
            comps.append((a,b,ln,pts[a],pts[b]))
    return comps

def owner_pick(owners, eps_want):
    """choose an owner tuple matching eps_want (left endpoint wants a runner whose
    +1 edge is here; right endpoint wants -1 edge)."""
    cand=[o for o in owners if o[2]==eps_want]
    return cand[0] if cand else (owners[0] if owners else None)

def interval_coverable(a_owner,b_owner,w,n):
    """does a v=nw arc cover both endpoints?  via the congruence window on j."""
    if a_owner is None or b_owner is None: return None
    (ua,ka,ea)=a_owner; (ub,kb,eb)=b_owner
    Aa=w*(ka*n+ea); Ab=w*(kb*n+eb)
    # need integer j with |Aa - j ua|<ua/n  AND  |Ab - j ub|<ub/n
    # j-window from a: (Aa/ua - 1/n, Aa/ua + 1/n); from b similarly
    loa=F(Aa,ua)-F(1,n); hia=F(Aa,ua)+F(1,n)
    lob=F(Ab,ub)-F(1,n); hib=F(Ab,ub)+F(1,n)
    lo=max(loa,lob); hi=min(hia,hib)
    # integer j strictly in (lo,hi)?
    import math
    jmin=floor(lo)+1
    while jmin<=hi:
        if lo<jmin<hi: return jmin
        jmin+=1
        if jmin>hi+2: break
    # robust integer search
    for j in range(floor(lo),ceil(hi)+1):
        if lo<j<hi: return j
    return None

def prim(V):
    g=0
    for v in V: g=gcd(g,v)
    return tuple(sorted(v//g for v in V))

def tight_direct(V,n):
    """G(V) empty (open)? i.e. measure-0 / no open safe interval."""
    THR=F(1,n); eps=set([F(0)])
    for v in V:
        for k in range(v+1):
            for s in(-1,1): eps.add(F(k*n+s,n*v)%1)
    pts=sorted(eps); L=len(pts)
    for i in range(L):
        a=pts[i]; b=pts[(i+1)%L]; ln=(b-a) if b>a else (b-a+1); mid=(a+ln/2)%1
        if all(dist(v*mid)>THR for v in V): return False
    return True

def main():
    rng=random.Random(23)
    print("Translator verification: config TIGHT <=> every G(S') component coverable by congruence")
    for n in [6,8,10,12,14]:
        m=n-1; tot=0; agree=0; small_owner=0; owner_tot=0; tight_cnt=0
        for _ in range(2500):
            others=set(rng.sample([x for x in range(1,n+8) if x%n!=0],m-1))
            w=rng.randint(1,3); v=n*w
            if v in others: continue
            V=prim(tuple(sorted(others|{v})))
            if len(V)!=m or not any(x%n==0 for x in V): continue
            mults=[x for x in V if x%n==0]; vv=mults[0]; ww=vv//n
            Sp=tuple(x for x in V if x!=vv)
            comps=G_components(Sp,n)
            if not comps: continue
            tot+=1
            # translator prediction: tight iff every component coverable
            pred_tight=True
            for (a,b,ln,oa,ob) in comps:
                ao=owner_pick(oa,1); bo=owner_pick(ob,-1)
                for o in (ao,bo):
                    if o is not None:
                        owner_tot+=1
                        if o[0]<n: small_owner+=1
                if interval_coverable(ao,bo,ww,n) is None:
                    pred_tight=False; break
            actual=tight_direct(V,n)
            if pred_tight==actual: agree+=1
            if actual: tight_cnt+=1
        print(f"  n={n:2d}: {tot} configs; translator agrees with direct tight/loose = {agree}/{tot}; "
              f"endpoints owned by speed<n: {small_owner}/{owner_tot} ({100*small_owner/max(owner_tot,1):.0f}%); "
              f"actually tight={tight_cnt}")
if __name__=='__main__': main()

def cross_relation_audit():
    """For components with BOTH owners < n, the w-free necessary cover condition is
    u_b(k_a n+1) = u_a(k_b n-1). Measure how often it holds vs fails (failure => loose
    independent of w)."""
    import random
    rng=random.Random(91)
    print()
    print("w-FREE cross-relation u_b(k_a n+1)=u_a(k_b n-1) on small-owner (<n) components")
    for n in [6,8,10,12,14]:
        m=n-1; hold=0; fail=0
        for _ in range(3000):
            others=set(rng.sample([x for x in range(1,n+8) if x%n!=0],m-1))
            Sp=tuple(sorted(others))
            for (a,b,ln,oa,ob) in G_components(Sp,n):
                ao=owner_pick(oa,1); bo=owner_pick(ob,-1)
                if ao and bo and ao[0]<n and bo[0]<n:
                    (ua,ka,_)=ao; (ub,kb,_)=bo
                    if ub*(ka*n+1)==ua*(kb*n-1): hold+=1
                    else: fail+=1
        tot=hold+fail
        print(f"  n={n:2d}: small-owner components: cross-relation HOLDS {hold}/{tot} "
              f"({100*hold/max(tot,1):.2f}%), FAILS {fail} (=> those components uncoverable for ANY w)")
if __name__=='__main__': cross_relation_audit()
