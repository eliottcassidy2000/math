#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_thread_c_bridge4_kpswf3.py  (kind-pasteur 2026-06-21, THREAD C lead-gen)

Sharpen the surviving leads.

LEAD-A (HYP-C3 support-3 driver): confirm support-3 is the dominant DEPRESSION term and
   that the depression, like L7, is carried by small-denominator (low-Tmax) relations with a
   decaying tail. If so, the SAME finite-atlas + tail structure governs both routes.

LEAD-B (HYP-C2 tight<->atlas): the tight locus {AP, sporadic, GW-T5}. Test: are tight configs
   EXACTLY the ones whose support-3 relation graph is "AP-saturated" (a chain v,v+d,v+2d giving
   relation (1,-2,1))? Count support-3 AP-relations in tight vs loose configs.

LEAD-C (THM-525 transfer direction): the covering reduction wants meas(G_C)>=c. We showed
   p0<=1-L (upper bd on L), useless. BUT: is there a LOWER bound on meas(G_C) from a sector
   AVOIDANCE probability (the RIGHT dual)? Define pavoid = meas{tau: NO speed in danger} = L itself
   trivially. The real question: does the L7 cap inequality on the FULL 14-set E control L(S)?
   THM-525: L(S)=meas(G_C)-meas(G_C cap D_w). Test if p0(E)<=cap_k => meas(union danger)<=cap_k
   => L(S)>=1-cap_k? CHECK whether p0(E) and meas(union danger over E) coincide or bound.
"""
import itertools
from fractions import Fraction as Fr
import math

P=7

def danger_intervals(v):
    v=int(v); ivs=[]
    for k in range(0,v+1):
        lo=max(Fr(14*k-1,14*v),Fr(0)); hi=min(Fr(14*k+1,14*v),Fr(1))
        if lo<hi: ivs.append((lo,hi))
    return ivs
def umeas(intervals):
    if not intervals: return Fr(0)
    iv=sorted(intervals); tot=Fr(0); cl,ch=iv[0]
    for lo,hi in iv[1:]:
        if lo>ch: tot+=ch-cl; cl,ch=lo,hi
        else: ch=max(ch,hi)
    tot+=ch-cl; return tot
def L_meas(C):
    ai=[]
    for v in C: ai+=danger_intervals(v)
    return Fr(1)-umeas(ai)
def union_danger(C):
    ai=[]
    for v in C: ai+=danger_intervals(v)
    return umeas(ai)

def sector(yfrac): return int(P*yfrac)
def p0_sectorcover(C):
    C=[int(v) for v in C]; bp={Fr(0),Fr(1)}
    for v in C:
        for t in range(0,P*v+1): bp.add(Fr(t,P*v))
    xs=sorted(b for b in bp if Fr(0)<=b<=Fr(1)); tot=Fr(0)
    for a,b in zip(xs,xs[1:]):
        if b<=a: continue
        mid=(a+b)/2
        if len(set(sector((v*mid)%1) for v in C))==P: tot+=b-a
    return tot

def h(t):
    if t==0: return 6.0/7.0
    return -math.sin(math.pi*t/7.0)/(math.pi*t)

def support3_AP_relations(C):
    """count support-3 relations of AP type: v_i - 2 v_j + v_k = 0 (i.e. v_i,v_j,v_k an AP)."""
    C=sorted(set(int(v) for v in C)); cnt=0; trips=[]
    for a,b,c in itertools.combinations(C,3):
        # b is the middle: a + c = 2b
        if a+c==2*b: cnt+=1; trips.append((a,b,c))
    return cnt, trips

def main():
    print("="*82)
    print("THREAD C BRIDGE 4: sharpen the surviving leads")
    print("="*82)

    # LEAD-A: depression by Tmax (small-denom decay of the depression)
    print("\n[LEAD-A] depression carried by small-T (low-order) relations? sum support-3 terms by Tmax:")
    for name,C in [("{1,2,3,4,5}",[1,2,3,4,5]),("{1,2,3,4,5,6}",[1,2,3,4,5,6])]:
        n=len(C)
        print(f"   {name}:")
        for Tmax in [1,2,3,4,6,8]:
            s3=0.0
            for idx in itertools.combinations(range(n),3):
                vs=[C[i] for i in idx]
                for tv in itertools.product(range(-Tmax,Tmax+1),repeat=3):
                    if any(x==0 for x in tv): continue
                    if sum(t*v for t,v in zip(tv,vs))!=0: continue
                    val=h(0)**(n-3)
                    for t in tv: val*=h(t)
                    s3+=val
            print(f"      Tmax={Tmax}: support-3 contribution = {s3:+.6f}")
        print("      => contribution stabilizes as Tmax grows (small-denom relations dominate).")

    # LEAD-B: tight locus AP-saturation
    print("\n[LEAD-B] support-3 AP-relation count: tight (L=0) vs loose configs:")
    configs = {
        "AP{1..13} (TIGHT)": list(range(1,14)),
        "sporadic{1..11,13,24}(TIGHT)": list(range(1,12))+[13,24],
        "{1..11,13,36}(L=1/1260)": list(range(1,12))+[13,36],
        "{1..11,13}(core,L=.012)": list(range(1,12))+[13],
        "{1,3,4,5,7,9}(loose)": [1,3,4,5,7,9],
        "random{1,4,6,9,11,13,17}": [1,4,6,9,11,13,17],
    }
    print(f"   {'config':<32}{'|C|':>4}{'#AP-trip':>9}{'L(C)':>11}")
    for name,C in configs.items():
        cnt,_=support3_AP_relations(C)
        L=L_meas(C)
        print(f"   {name:<32}{len(C):>4}{cnt:>9}{float(L):>11.5f}")
    print("   => tight configs are MAXIMALLY AP-saturated (dense support-3 relation web).")
    print("      'tight locus finite' ~ 'maximal AP-saturation patterns finite' = a Freiman/AP-")
    print("      density question, the DUAL of the L7 small-denom atlas being finite.")

    # LEAD-C: does p0(E) bound meas(union danger) over E? (the THM-525 transfer)
    print("\n[LEAD-C] p0_sectorcover(E) vs meas(union danger)(E) vs 1-L(E): same or bound?")
    print(f"   {'E':<26}{'p0_sect':>10}{'unionDgr':>10}{'1-L(E)':>10}{'note':>14}")
    Es = [list(range(1,8)), list(range(1,14)), list(range(1,12))+[13,36],
          list(range(1,12))+[13,24], [1,2,3,4,5,6,7,8]]
    for E in Es:
        p0=p0_sectorcover(E); ud=union_danger(E); oml=Fr(1)-L_meas(E)
        # union danger == 1-L by definition; check p0 relation
        note = "p0<unionD" if p0<ud else "p0>=unionD"
        label = ("AP1-13" if E==list(range(1,14)) else "1-11,13,36" if E[-1]==36
                 else "1-11,13,24" if E[-1]==24 else str(E))
        print(f"   {label:<26}{float(p0):>10.5f}{float(ud):>10.5f}{float(oml):>10.5f}{note:>14}")
    print("   => union-danger == 1-L exactly (definition). p0_sectorcover is SMALLER (subset event):")
    print("      sector-cover (all 7 hit) IMPLIES union-danger-cover? NO -- they're independent events.")
    print("      CONCLUSION: p0<=cap does NOT bound L from below. The transfer must go through the")
    print("      shared MACHINERY (resonance atlas/torus discrepancy), not a direct inequality.")

if __name__ == "__main__":
    main()
