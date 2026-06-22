#!/usr/bin/env python3
"""Test two strategic levers: (A) Node-2 slack-split (binding only k=8,9,10?),
(B) Node-3 quasi-independence ratio R' (exact, sector-cover proxy). kind-mendel-S2."""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
def lcm(a,b): return a*b//gcd(a,b)
def gl(xs): return reduce(lcm,[x for x in xs if x],1)
def gall(xs): return reduce(gcd,[x for x in xs if x],0)

def p0_fast(E):
    pos=[e for e in E if e!=0]
    if not pos: return F(0)
    D=7*gl(pos); bset=set([0,D])
    for e in pos:
        step=D//(7*e); x=0
        while x<=D: bset.add(x); x+=step
    bps=sorted(bset); num=0
    for a,b in zip(bps,bps[1:]):
        if b<=a: continue
        mid2=a+b; sec=set()
        for e in E:
            sec.add((7*((e*mid2)%(2*D)))//(2*D))
            if len(sec)==7: break
        if len(sec)==7: num+=b-a
    return F(num,D)

capvals={8:F(2243,5880),9:F(1979,4004),10:F(55,91),11:F(66,91),12:F(6,7),13:F(1)}

print("=== LEVER A (Node 2): max p0 over bounded-spread E, and wide asymptote, vs cap_k ===")
print("    (if margin huge at k>=11, those need only a CRUDE bound, not exact extremality)")
for k in [8,9,10,11,12,13]:
    cap=capvals[k]
    # bounded-spread search (spread<=k+3, modest)
    S=k+3; best=(F(-1),None)
    if k<=11:
        for sub in combinations(range(1,S+1),k-1):
            E=[0]+list(sub)
            if gall(E)!=1: continue
            v=p0_fast(E)
            if v>best[0]: best=(v,E)
    else:
        # k=12,13 too big to search; test reported maximizer candidates
        cands=[list(range(k)), [0]+list(range(2,k+1)), list(range(k-1))+[k]]
        for E in cands:
            if gall(E)!=1: continue
            v=p0_fast(E)
            if v>best[0]: best=(v,E)
    # wide asymptote: k-1 consecutive + one far runner pushed out
    Wfar=[0]+list(range(1,k-1))+[97]   # last element far & coprime-ish
    pw=p0_fast(Wfar)
    consec=p0_fast(list(range(k)))
    print(f"k={k:2d}: cap={float(cap):.4f} | consec p0={float(consec):.4f} | "
          f"bounded-max p0={float(best[0]):.4f} (margin {float(cap-best[0]):+.4f}) | "
          f"one-far p0={float(pw):.4f}")

print("\n=== LEVER B (Node 3): quasi-independence R' = meas(coverSet^c∩G_P)/(measGP*(1-p0)) ===")
print("    R'≈1 => decorrelation; witnessG2 >= measGP*(1-p0)*R' = floor. (exact rationals)")
def measGP_and_inter(E,P):
    "return measGP, p0=meas(coverSet), meas(coverSet^c ∩ GP), all exact"
    D=14*gl(list(E)+list(P)); bset=set([0,D])
    for e in E:
        if e==0: continue
        step=D//(7*e); x=0
        while x<=D: bset.add(x); x+=step
    for p in P:
        for n in range(p):
            c=n*D//p; bset.add((c-D//(14*p))%D); bset.add((c+D//(14*p))%D)
    bps=sorted(bset); ng=cov=covc=0; p0n=0
    for a,b in zip(bps,bps[1:]):
        if b<=a: continue
        mid2=a+b
        sec={(7*((e*mid2)%(2*D)))//(2*D) for e in E}
        allhit=(len(sec)==7)
        if allhit: p0n+=b-a
        inGP=all(min((p*mid2)%(2*D),2*D-(p*mid2)%(2*D))>=2*D//14 for p in P)
        if inGP:
            ng+=b-a
            if allhit: cov+=b-a
            else: covc+=b-a
    return F(ng,D),F(p0n,D),F(covc,D)
tests=[
 (8,[0,1,2,3,4,5,6,7],[1,5,7,8,9]),
 (8,[0,2,4,6,8,10,12,14],[1,3,5,7,9]),     # dilated even base
 (8,[0,1,2,3,4,5,6,7],[1,2,3,11,13]),
 (10,list(range(10)),[1,12,13]),
 (10,[0,1,2,3,4,5,6,7,8,9],[2,3,5]),
 (6,[0,1,2,3,4,5],[1,2,3,5,7,8,9]),
]
for k,E,P in tests:
    mg,p0v,covc=measGP_and_inter(E,P)
    indep=mg*(1-p0v)
    Rp = covc/indep if indep>0 else F(0)
    print(f"k={k} P={P}: measGP={float(mg):.4f} p0={float(p0v):.4f} | true wG2-proxy={float(covc):.4f} "
          f"| indep={float(indep):.4f} | R'={float(Rp):.4f} | Bonf(measGP-p0)={float(max(F(0),mg-p0v)):.4f}")
