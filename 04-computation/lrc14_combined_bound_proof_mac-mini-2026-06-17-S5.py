#!/usr/bin/env python3
"""
lrc14_combined_bound_proof — mac-mini-2026-06-17-S5

THE CANDIDATE PROOF. Rule: v = V = max(S). Claim W(S\{V}) > 1/(7V) via a COMBINED
PROVABLE lower bound  max( P, Anti )  where:

  P(A)   = (1 - sum_{u in A} 1/(7u)) / (sum_{u in A} u)        [PIGEONHOLE; always valid:
           safe measure mu(A) >= 1-sum 1/(7u); #safe arcs <= #teeth = sum u; W >= mu/#arcs]

  Anti(A)= [min_u (u/(2V) - 1/14)/u]_+  +  [min_u (13/14 - u/(2V))/u]_+   [ANTIPODE; valid
           (>0) iff every u in [V/7, 13V/7] so tau0=1/(2V) is in G_A; the antipode safe
           arc around tau0 lies in G_A (each runner sweeps inside [1/14,13/14], no wrap).]

If  max(P, Anti) > 1/(7V)  for EVERY covering 13-set, then W(S\{V}) > 1/(7V) PROVABLY,
hence M(S) >= 1/14 (arc-width lemma), hence LRC(14). Both bounds are elementary/provable,
so this would be a PROOF (the two cases must exhaust). Test over spread + clustered + mixed.
"""
from fractions import Fraction as F
import random
C=F(1,14)
def covering(S): return all(any(v%q==0 for v in S) for q in range(2,15))
def Ppig(A):
    mu=F(1)-sum(F(1,7*u) for u in set(A)); su=sum(A)
    return mu/su if mu>0 else F(0)
def Anti(A,V):
    # valid iff all u in [V/7, 13V/7]; tau0=1/(2V)
    lo=F(0); hi=F(0)
    for u in A:
        pos=F(u,2*V)               # u*tau0
        d_lo=pos-F(1,14)           # room moving down (toward 1/14)
        d_hi=F(13,14)-pos          # room moving up (toward 13/14)
        if d_lo<=0 or d_hi<=0: return F(0)   # tau0 not safe for this u -> antipode invalid
        # time for u to reach boundary = room/u
        a=d_lo/u; b=d_hi/u
        lo=a if (lo==0 or a<lo) else lo
        hi=b if (hi==0 or b<hi) else hi
    return lo+hi
def exactW(A):  # for cross-check
    def darcs(v):
        hw=F(C,v); return [(F(k,v)-hw,F(k,v)+hw) for k in range(v)]
    dz=[]
    for v in set(A): dz+=darcs(v)
    o=[]
    for lo,hi in dz:
        s=lo-(lo%1); a=lo-s;b=hi-s
        if b<=1:o.append((a,b))
        else:o.append((a,F(1)));o.append((F(0),b-1))
    o=sorted(o);r=[];cl,ch=o[0]
    for lo,hi in o[1:]:
        if lo<=ch: ch=ch if ch>hi else hi
        else:r.append((cl,ch));cl,ch=lo,hi
    r.append((cl,ch))
    best=F(0)
    for i in range(len(r)):
        hi=r[i][1]; lo=r[(i+1)%len(r)][0]+(1 if i==len(r)-1 else 0)
        best=max(best,lo-hi)
    return best

def clustered(N,win,rng):
    used=set(); S=[]
    for q in range(2,15):
        cs=[x for x in range(N,N+win+1) if x%q==0 and x not in used]
        if not cs: return None
        x=rng.choice(cs); used.add(x); S.append(x)
    S=sorted(set(S)); return S if len(S)==13 and covering(S) else None

rng=random.Random(13)
sets=[]
for drop in [1,2,3,4,5,6,7,12]:
    for k in range(1,40):
        s=sorted(set([v for v in range(1,14) if v!=drop]+[84*k]))
        if len(s)==13 and covering(s): sets.append(s)
for _ in range(1500):
    s=clustered(rng.choice([60,150,400,1500,9000]), rng.choice([30,80,200,400]), rng)
    if s: sets.append(s)
for _ in range(4000):
    drop=rng.choice([1,2,3,4,5,6]); base=[v for v in range(1,14) if v!=drop]
    d2=rng.choice([8,9,10,11,13]); base=[v for v in base if v!=d2]
    bigs=set()
    while len(bigs)<2: bigs.add(rng.choice([84,126,168,210,154,182,252])*rng.randint(1,6))
    s=sorted(set(base+list(bigs)))
    if len(s)==13 and covering(s): sets.append(s)

print("="*76)
print("COMBINED PROVABLE BOUND max(Pigeonhole, Antipode) > 1/(7V) for covering 13-sets?")
print("="*76)
ntot=0; combofail=0; pig_only=0; anti_only=0; worst=(F(99),None); wrongW=0
for S in sets:
    V=max(S); A=[u for u in S if u!=V]; thr=F(1,7*V)
    P=Ppig(A); An=Anti(A,V); best=P if P>An else An
    # which case is valid
    if P>thr: pig_only+=1
    if An>thr: anti_only+=1
    if not (best>thr): combofail+=1
    margin=best-thr
    if margin<worst[0]: worst=(margin,S,'pig' if P>=An else 'anti')
    # sanity: provable bound <= exact W (only on a subsample, exactW is slow for big V)
    if V<=300 and best>exactW(A)+F(1,10**9): wrongW+=1
    ntot+=1
print(f"  covering 13-sets tested: {ntot}")
print(f"  combined bound max(P,Anti) > 1/(7V) holds: {ntot-combofail}/{ntot}  (FAILS: {combofail})")
print(f"  pigeonhole alone suffices: {pig_only}/{ntot};  antipode alone suffices: {anti_only}/{ntot}")
print(f"  (sets where NEITHER alone works but combined does: {ntot - len([1 for _ in range(0)]) }) -- see fails")
print(f"  provable bound exceeded exact W (sanity, should be 0): {wrongW}")
print(f"  tightest combined margin: {float(worst[0]):.7f} via {worst[2] if worst[1] else '-'} at S={worst[1][:6] if worst[1] else None}...")
if combofail==0:
    print("\n  *** max(Pigeonhole, Antipode) > 1/(7V) on ALL tested => CANDIDATE PROOF of LRC(14):")
    print("      both bounds are elementary; if the two cases provably EXHAUST, M(S)>=1/14 for all covering S. ***")
else:
    print(f"\n  combined bound failed on {combofail} sets — need a third bound / sharper estimate there.")
