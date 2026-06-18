#!/usr/bin/env python3
"""
lrc14_three_bound_proof — mac-mini-2026-06-17-S5

CANDIDATE PROOF (3 elementary provable lower bounds on W(S\{V}), V=max(S)):
  P    = pigeonhole  mu/Sum(u)
  AntiV= antipode safe-arc at tau=1/(2V)
  AntiU= max over core u0 of antipode safe-arc at tau=1/(2 u0)   (generalizes AntiV)
Claim: max(P, AntiU) > 1/(7V) for EVERY covering 13-set  =>  W(S\{V})>1/(7V)
       => M(S)>=1/14 (arc-width lemma) => LRC(14). All bounds elementary => proof if universal.
Aggressive sampling (spread, clustered, mixed, two-equal-large, scaled).
"""
from fractions import Fraction as F
import random
C=F(1,14)
def covering(S): return all(any(v%q==0 for v in S) for q in range(2,15))
def Ppig(A):
    mu=F(1)-sum(F(1,7*u) for u in set(A)); su=sum(A); return mu/su if mu>0 else F(0)
def antipode_at(A,u0):
    """safe-arc width of A around tau=1/(2 u0) (0 if not safe there)."""
    tau=F(1,2*u0); lo=hi=None
    for u in A:
        fr=(u*tau)%1; d=min(fr,1-fr)
        if d<C: return F(0)
        dl=(fr-F(1,14))/u; dh=(F(13,14)-fr)/u
        if dl<0 or dh<0: return F(0)
        lo=dl if lo is None or dl<lo else lo
        hi=dh if hi is None or dh<hi else hi
    return (lo or F(0))+(hi or F(0))
def AntiU(A):
    return max([antipode_at(A,u0) for u0 in set(A)]+[F(0)])
def bound(S):
    V=max(S); A=[u for u in S if u!=V]
    return max(Ppig(A), AntiU(A)), F(1,7*V), A, V

def clustered(N,win,rng):
    used=set(); S=[]
    for q in range(2,15):
        cs=[x for x in range(N,N+win+1) if x%q==0 and x not in used]
        if not cs: return None
        x=rng.choice(cs); used.add(x); S.append(x)
    S=sorted(set(S)); return S if len(S)==13 and covering(S) else None

rng=random.Random(99); sets=[]
for drop in [1,2,3,4,5,6,7,12]:
    for k in range(1,60):
        s=sorted(set([v for v in range(1,14) if v!=drop]+[84*k]))
        if len(s)==13 and covering(s): sets.append(s)
for _ in range(2500):
    s=clustered(rng.choice([40,90,200,600,2000,9000,40000]), rng.choice([26,60,140,320,700]), rng)
    if s: sets.append(s)
for _ in range(9000):
    nd=rng.choice([2,3]); base=list(range(1,14))
    for _ in range(nd):
        if base: base.remove(rng.choice(base))
    bigs=set()
    while len(bigs)<nd: bigs.add(rng.choice([84,126,168,210,154,182,252,330,390,420])*rng.randint(1,8))
    s=sorted(set(base+list(bigs)))
    if len(s)==13 and covering(s): sets.append(s)
# two-equal-large
for _ in range(2000):
    drop=rng.choice([1,2,3,4,5]); base=[v for v in range(1,14) if v!=drop]
    d2=rng.choice([9,10,11,13]); base=[v for v in base if v!=d2]
    b=rng.choice([84,168,252,420]); s=sorted(set(base+[b,b+rng.choice([14,28,42,56,84])]))
    if len(s)==13 and covering(s): sets.append(s)

print("="*76)
print("THREE-BOUND max(Pigeonhole, AntipodeAny) > 1/(7V) for covering 13-sets?")
print("="*76)
ntot=0; fail=0; worst=(F(99),None); byfam={}
for S in sets:
    b,thr,A,V=bound(S); ntot+=1
    if not (b>thr):
        fail+=1
        if fail<=5: print("  FAIL:",S,"  bound-thr=",float(b-thr))
    m=b-thr
    if m<worst[0]: worst=(m,S)
print(f"  covering 13-sets tested: {ntot}")
print(f"  max(P,AntiU) > 1/(7V) holds: {ntot-fail}/{ntot}   (FAILS: {fail})")
print(f"  tightest margin: {float(worst[0]):.8f} at S={worst[1][:7] if worst[1] else None}... V={max(worst[1]) if worst[1] else '-'}")
if fail==0:
    print("\n  *** UNIVERSAL on this aggressive sample => CANDIDATE PROOF of LRC(14): ***")
    print("      v=max(S); W(S\\{V}) >= max(pigeonhole, antipode-any) > 1/(7V); arc-width lemma => M>=1/14.")
    print("      Remaining: prove this 2-bound max > 1/(7V) for ALL covering 13-sets (the two cases exhaust).")
else:
    print(f"\n  {fail} residual sets need a further bound.")
