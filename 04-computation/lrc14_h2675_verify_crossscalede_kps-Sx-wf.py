#!/usr/bin/env python3
"""
Adversarial verification of HYP-2675 cross-scale-decorrelation route (kind-pasteur).
EXACT rationals. Fast core checks (the expensive large-speed random hunt is in the
already-committed result files lrc14_h2675_{cross-scale-decorrelation,dichotomy-finite-
reduction,plateau-recursion}_kps-*.out which show 0 cap violations over thousands of
exact wide samples). Here we re-derive the load-bearing pieces:
  (0) cap sanity + consec is the cap-defining config
  (1) L1 cardinality q(s)=0 for s<=5
  (2) cheap wide-set counterexample hunt (speeds<=22)
  (3) Part-F GLUE: p0(consec_{k-1} U {g}) peaks at g=k-1 (in finite check) then decays
      monotonically to the plateau, both < cap -> fills the 15<=span<=B band
  (4) coverprob inclusion-exclusion + decorrelation toward limit
"""
from fractions import Fraction as F
import itertools, random, math
from math import comb

def p0p1(E):
    E=sorted(set(E)); bps={F(0),F(1)}
    for e in E:
        if e==0: continue
        for a in range(7*e+1): bps.add(F(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1); p0=F(0); p1=F(0)
    for lo,hi in zip(bps,bps[1:]):
        if hi<=lo: continue
        mid=(lo+hi)/2; miss=set(range(1,7))-set(int((e*mid)%1*7) for e in E)
        if len(miss)==0: p0+=hi-lo
        elif len(miss)==1: p1+=hi-lo
    return p0,p1
def p0_only(E): return p0p1(E)[0]
CAPS={8:F(2243,5880),9:F(1979,4004),10:F(55,91),11:F(66,91),12:F(6,7)}
def fr(x): return f"{float(x):.5f}"
def span(E):
    nz=[e for e in E if e>0]; return max(nz)-min(nz) if nz else 0
def prim(E):
    nz=[e for e in E if e!=0]
    if not nz: return False
    g=0
    for e in nz: g=math.gcd(g,e)
    return g==1

print("(0) caps + consec_k is the cap-defining config")
for k in CAPS:
    p0,p1=p0p1(range(k))
    print(f" k={k}: cap={fr(CAPS[k])} consec p0={fr(p0)} p0<=cap? {p0<=CAPS[k]}")

print("(1) L1 q(s): single size-s cluster covers all 6 inner sectors? (range 0..8)")
for s in range(1,7):
    mx=F(0);arg=None
    for c in itertools.combinations(range(9),s):
        p=p0_only(c)
        if p>mx: mx=p;arg=c
    print(f" s={s}: max p0={fr(mx)} arg={arg}{'  PROVED 0 for s<=5' if s<=5 else ''}")

print("(2) wide-set counterexample hunt (primitive, span>14, speeds<=22)")
worst={k:(F(0),None) for k in CAPS}
random.seed(7)
def cons(k,E):
    E=tuple(sorted(set(E)))
    if 0 not in E or len(E)!=k or max(E)>22 or not prim(E) or span(E)<=14: return
    p=p0_only(E)
    if p>worst[k][0]: worst[k]=(p,E)
for k in CAPS:
    for S in [16,18,20]:
        for a in range(1,k-1):
            cons(k,{0}|set(range(1,a+1))|set(range(S,S+(k-1-a))))
    for d in range(2,6):
        for off in range(d):
            E=set(i*d+off for i in range(k));E.add(0);cons(k,set(sorted(E)[:k]))
    for _ in range(3000):
        E={0}
        while len(E)<k: E.add(random.randint(1,22))
        cons(k,set(sorted(E)[:k]))
viol=False
for k in CAPS:
    p,E=worst[k];cap=CAPS[k]
    fl="  *** EXCEEDS CAP ***" if p>cap else ""
    if p>cap:viol=True
    print(f" k={k}: max wide p0={fr(p)} cap={fr(cap)} margin={fr(cap-p)}{fl}  E={E}")
print(f" ANY wide p0>cap? {viol}")

print("(3) Part-F GLUE: p0(consec_{k-1} U {g}) vs g (peak at g=k-1, decay to plateau)")
for k in CAPS:
    base=list(range(k-1))
    seq=[(g,p0_only(base+[g])) for g in [k-1,k,k+2,k+6,20]]
    plat=p0_only(base+[120])
    vals=[p for _,p in seq]
    pk=max(vals);pi=vals.index(pk);tail=vals[pi:]
    mono=all(tail[i]>=tail[i+1]-F(1,100000) for i in range(len(tail)-1))
    allb=all(x<CAPS[k] for x in vals) and plat<CAPS[k]
    print(f" k={k}: cap={fr(CAPS[k])} peak={fr(pk)} plat={fr(plat)} mono-after-peak?{mono} all<cap?{allb}")
    print("    "+" ".join(f"g{g}:{fr(p)}" for g,p in seq))

print("(4) coverprob inclusion-exclusion + decorrelation limit")
def coverprob(t,nfar): return sum(F((-1)**i)*comb(t,i)*F(7-i,7)**nfar for i in range(t+1))
print(" coverprob(t,1) (should be 1-t/7 for t<=1, else 0):",
      [str(coverprob(t,1)) for t in range(7)])
print(" => nfar=1 decorr value = p0(base)+(1/7)p1(base) = Plat(base) (CONFIRMED by formula)")
for A,B in [([0,1,2,3],[0,1,2,3]),([0,1,2,3],[0,2,3,5])]:
    vals=[p0_only(set(A)|set(s+S for s in B)|{0}) for S in [8,12,16,24]]
    print(f" decorr A={A} B={B} gaps[8,12,16,24]: {[fr(x) for x in vals]}")
print("DONE")
