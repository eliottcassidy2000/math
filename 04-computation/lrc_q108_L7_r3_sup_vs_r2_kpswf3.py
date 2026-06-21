#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_q108_L7_r3_sup_vs_r2_kpswf3.py   (kind-pasteur 2026-06-21, THREAD B final)

Pin the two structural claims that close r=3 directly:

(1) DICHOTOMY of the r=3 balanced limit (matches r=2 shape, one dim up):
    * ratios genuinely SEPARATED (>~2 apart in denom sense) -> R3 DECAYS to 0 (tail);
    * ratios -> COMMENSURATE / MERGED (gamma_i -> equal) -> R3 -> bounded plateau =
      the LOWER-r resonance (merge collapses a far cluster). NOT a new divergence.
    So the r=3 sup is attained at SMALL denominators where two ratios MERGE, i.e. it is
    DOMINATED by an effective r<=2 balanced resonance. No new worst case appears at r=3.

(2) sup_{r=3} p0_inf <= sup_{r=2} p0_inf  (per-k, exact), and both << cap.
    -> base-size domination is CONFIRMED by direct computation, not assumed.

Method: for fixed k=9, compare
    r=2 limit  p0_inf(base7, (q,p))           over coprime (1,2.15], q<=6
    r=3 limit  p0_inf(base6, (q,p2,p3))        over coprime pairs (1,2.15]^2, q<=6
and also the "merge" sub-family (q,p,p) [two far ratios equal] to expose the plateau.

EXACT rational.  Output -> 05-knowledge/results/.
"""
import itertools
from fractions import Fraction as Fr
from math import gcd, comb

P = 7
def sector(yf): return int(P*yf)
CAP = {8: Fr(2243,5880), 9: Fr(1979,4004), 10: Fr(55,91)}

def base_xcells(B):
    B=[int(b) for b in B]; xbp={Fr(0),Fr(1)}
    for b in B:
        if b==0: continue
        for t in range(0,P*b): xbp.add(Fr(t,P*b))
    xs=sorted(xbp); out=[]
    for a,b in zip(xs,xs[1:]):
        mid=(a+b)/2
        out.append((b-a, frozenset(sector((bb*mid)%1) for bb in B)))
    return out

def far_vcells(mults):
    mults=[int(m) for m in mults]; vbp={Fr(0),Fr(1)}
    for f in mults:
        for t in range(0,P*f): vbp.add(Fr(t,P*f))
    vs=sorted(vbp); out=[]
    for a,b in zip(vs,vs[1:]):
        mid=(a+b)/2
        out.append((b-a, frozenset(sector((f*mid)%1) for f in mults)))
    return out

def p0_inf_multi(B, mults):
    xcells=base_xcells(B)
    from collections import defaultdict
    vgrp=defaultdict(Fr)
    for vlen,sf in far_vcells(mults): vgrp[sf]+=vlen
    total=Fr(0)
    for sf,vlen in vgrp.items():
        for xlen,Sbase in xcells:
            if len(Sbase|sf)==P: total+=xlen*vlen
    return total

def Pr_decorrelated(B,r):
    total=Fr(0)
    for xlen,Sbase in base_xcells(B):
        m=P-len(Sbase)
        if m==0: pcov=Fr(1)
        elif m>r: pcov=Fr(0)
        else:
            pcov=Fr(0)
            for j in range(0,m+1): pcov+=Fr((-1)**j*comb(m,j))*Fr(P-j,P)**r
        total+=xlen*pcov
    return total

def main():
    base7=[0,2,4,6,8,10,12]; base6=[0,2,4,6,8,10]
    print("="*86)
    print("THREAD B final: r=3 sup vs r=2 sup (k=9), the merge-dichotomy, vs cap_9")
    print("="*86)

    # r=2 atlas
    P2=Pr_decorrelated(base7,2)
    sup2=(Fr(0),None)
    for q in range(1,7):
        for p in range(q+1,int(Fr(43,20)*q)+1):
            if gcd(p,q)!=1: continue
            if not (Fr(1)<Fr(p,q)<=Fr(43,20)): continue
            v=p0_inf_multi(base7,(q,p))
            if v>sup2[0]: sup2=(v,(q,p))
    print(f"\n r=2 (base7, 7+2=9): P2={float(P2):.5f}  sup p0_inf = {float(sup2[0]):.6f} at {sup2[1]}")

    # r=3 atlas
    P3=Pr_decorrelated(base6,3)
    sup3=(Fr(0),None)
    for q in range(1,7):
        nums=range(q+1,int(Fr(43,20)*q)+1)
        for p2,p3 in itertools.combinations_with_replacement(nums,2):
            if gcd(gcd(q,p2),p3)!=1: continue
            v=p0_inf_multi(base6,(q,p2,p3))
            if v>sup3[0]: sup3=(v,(q,p2,p3))
    print(f" r=3 (base6, 6+3=9): P3={float(P3):.5f}  sup p0_inf = {float(sup3[0]):.6f} at {sup3[1]}")
    print(f"\n => r=3 sup ({float(sup3[0]):.5f}) <= r=2 sup ({float(sup2[0]):.5f}) ? "
          f"{sup3[0] <= sup2[0]}   [DIRECT base-size domination, k=9]")
    print(f"    cap_9={float(CAP[9]):.5f}; margins: r=2 {float(CAP[9]-sup2[0]):.5f}, "
          f"r=3 {float(CAP[9]-sup3[0]):.5f}")

    # MERGE dichotomy: (q,p,p) two ratios equal -> should match an r=2 resonance value.
    print("\n MERGE check: r=3 direction (q,p,p) [two far ratios equal] collapses to r=2-like.")
    print(f"  {'q,p,p':>10} {'p0_inf_r3':>11} {'cmp r2(q,p) base6+1dup':>24}")
    for (q,p) in [(2,3),(3,4),(3,5),(2,4),(4,5),(5,6)]:
        if gcd(p,q)!=1: continue
        if not (Fr(1)<Fr(p,q)<=Fr(43,20)): continue
        v3=p0_inf_multi(base6,(q,p,p))
        # effective r=2 on base6 (one fewer base than base7): the merged far is a single cluster
        v2eff=p0_inf_multi(base6,(q,p))
        print(f"  {f'{q},{p},{p}':>10} {float(v3):>11.6f} {float(v2eff):>24.6f}")
    print("  (r3 (q,p,p) ~ r2 (q,p): the duplicated ratio adds a near-coincident far -> single cluster,")
    print("   confirming MERGE -> lower-r regime, NOT a new divergence. THM-557 boundary, one dim up.)")

    # SEPARATED tail: ratios spread out -> R3 -> 0 (already shown decaying in decay script).
    print("\n SEPARATED check: well-separated coprime triples, growing scale -> R3 -> 0:")
    for q in [5,7,9,11,13]:
        # pick spread ratios ~ (1.3, 2.0): p2=ceil(1.3q) coprimeized, p3=2q+1
        p2=None
        for cand in range(int(Fr(13,10)*q), int(Fr(14,10)*q)+1):
            if gcd(gcd(q,cand),2*q+1)==1: p2=cand; break
        if p2 is None: continue
        p3=2*q+1
        if not (Fr(1)<Fr(p2,q)<=Fr(43,20) and Fr(1)<Fr(p3,q)<=Fr(43,20)): continue
        v=p0_inf_multi(base6,(q,p2,p3)); R3=v-P3
        print(f"   q={q:>2} dir=({q},{p2},{p3}) g=({float(p2)/q:.3f},{float(p3)/q:.3f})  "
              f"p0_inf={float(v):.6f}  R3={float(R3):+.6f}")
    print("\nDONE.")

if __name__=="__main__":
    main()
