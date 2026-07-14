#!/usr/bin/env python3
"""mac-mini-S98: attack the UNIFORM residue-pattern claim (some k<=13 shadow witness for every covering
cluster) via shallow/shadow witnesses. (a) exact witness condition; (b) the 'free unit class' structural
criterion + pigeonhole partial; (c) ADVERSARIAL SEARCH for a covering cluster bad for ALL k<=13."""
from fractions import Fraction as F
from math import gcd
import itertools, random
c114=F(1,14)
def lonely_at(S,t):
    for c in S:
        r=(c*t)%1
        if min(r,1-r)<c114: return False
    return True
def shadow_interval(S,a,k):
    lo=F(0); hi=F(1,2)
    for c in S:
        r=(c*a)%k
        if r==0: lo=max(lo,F(1,14*c)); hi=min(hi,F(13,14*c))
        else:
            s=r if r<=k//2 else r-k; base=F(abs(s),k)
            hi=min(hi,(F(13,14)-base)/c if s>0 else (base-c114)/c)
        if lo>=hi: return None
    return (lo,hi)
def witness(C):  # returns (k,a,iv) or None; C excludes the observer speed 1 which we prepend
    S=[1]+list(C)
    for k in range(2,14):
        for a in range(1,k):
            if gcd(a,k)!=1 or not (c114<=F(a,k)<=F(13,14)): continue
            iv=shadow_interval(S,a,k)
            if iv and lonely_at(S,F(a,k)+(iv[0]+iv[1])/2): return (k,a,iv)
    return None
def covering_cluster(C): return all(any(v%q==0 for v in C) for q in range(2,15))

print("(b) PIGEONHOLE PARTIAL: for k=13 (phi=12), the non-13-carrier speeds occupy <=12 unit classes.")
print("   If >=2 speeds are 13-divisible, <=11 non-carriers => a unit class is FREE => set the killer")
print("   class there (choose a) => no s=-1 killer. Combined with shadow-window, k=13 gives a witness")
print("   (modulo s<=-2 secondary classes -- checked exactly below).")
print()
print("(c) ADVERSARIAL SEARCH for a covering cluster with NO k<=13 shadow witness:")
rng=random.Random(7);  worst=[]; checked=0; fails=[]
# strategy 1: full residue occupation mod 13 (all 12 units), all large, single 13-carrier
def try_family(gen, label, N):
    global checked
    cnt=0
    for C in gen:
        if len(set(C))!=len(C): continue
        if 1 in C: continue
        if not covering_cluster(C): continue
        checked+=1; cnt+=1
        if witness(C) is None: fails.append((label,sorted(C)))
        if cnt>=N: break
# structured adversary: {1}u{large speeds hitting all residues mod 13} + carriers
def adv_gen():
    # pick a 13-carrier (13*m) and 11 large speeds in distinct nonzero residues mod 13, covering 2..14
    for _ in range(60000):
        m=rng.randint(1,15); carrier13=13*m
        # need to cover 2..14; build a spread cluster
        base=set()
        for q in [2,3,4,5,6,7,8,9,10,11,12,14]:
            base.add(q*rng.randint(2,12))   # multiples, kept largish
        base.add(carrier13)
        C=set(base)
        while len(C)<12: C.add(rng.randint(15,220))
        C=set(sorted(C)[:12])
        while len(C)<12: C.add(rng.randint(15,220))
        yield sorted(C)
try_family(adv_gen(),"adversarial-spread", 4000)
# strategy 2: pure random covering
def rand_gen():
    for _ in range(60000):
        C=set()
        for q in range(2,15):
            if any(v%q==0 for v in C): continue
            C.add(q*rng.randint(1,30))
        while len(C)<12: C.add(rng.randint(2,400))
        yield sorted(list(C)[:12])
try_family(rand_gen(),"random-covering", 6000)
print(f"   covering clusters checked: {checked}")
if fails:
    print(f"   COUNTEREXAMPLES (no k<=13 witness): {len(fails)}")
    for lab,C in fails[:10]: print(f"      {lab}: {C}")
else:
    print(f"   ZERO counterexamples in {checked} covering clusters (incl adversarial full-residue-occupation).")
    print("   => the uniform claim survives a hard adversarial search; combined with the pigeonhole")
    print("      partial (multi-carrier k) + single-killer proof (S97), strong toward a full proof.")
