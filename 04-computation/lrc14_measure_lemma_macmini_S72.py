#!/usr/bin/env python3
"""mac-mini-S72: the FAR-ELEMENT MEASURE LEMMA (sharpening kps's 1/13-B/L to an exact threshold).
LEMMA: let C be a speed set, G_C = {t in [0,1): ||v t||>=1/13 for all v in C} its level-1/13
safe set (a union of arcs). If G_C has an arc of width > 1/L, then M(C u {L}) >= 1/13.
PROOF: an arc of width>1/L spans a full L-period, so contains a point t0 with ||L t0||=1/2>=1/13
(the midpoint between consecutive L-resonances); t0 in G_C gives ||v t0||>=1/13 for v in C.
=> M(C u {L}) >= min(1/13,1/2)=1/13. EXACT (not asymptotic), NO undershoot.
So M(C u {L}) >= 1/13 whenever L > 1/w_max(G_C), provided M(C)>1/13 (=> w_max>0).
TEST: does the lemma REACH the multi-killer extremals (largest outlier > 1/w_max)?"""
from fractions import Fraction as F
from math import gcd

def safe_arcs_width_max(C, eta=F(1,13)):
    """largest arc of G_C = {t: ||v t||>=eta for all v in C}, exact via breakpoints."""
    # breakpoints where some ||v t|| = eta: v t = m +- eta => t = (m +- eta)/v
    pts=set([F(0),F(1)])
    for v in C:
        for m in range(0,v+1):
            for s in (eta,-eta):
                t=(m+s)/v
                if 0<=t<=1: pts.add(t)
    P=sorted(pts)
    wmax=F(0); 
    for a,b in zip(P,P[1:]):
        mid=(a+b)/2
        if all(min((v*mid)%1,1-((v*mid)%1))>=eta for v in C):
            if b-a>wmax: wmax=b-a
    return wmax

def M_exact(S,Qmax):
    best=F(0)
    for q in range(2,Qmax+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            mind=q
            for v in S:
                r=(a*v)%q; d=r if r<=q-r else q-r
                if d<mind: mind=d
                if d==0: break
            if mind>0 and F(mind,q)>best: best=F(mind,q)
    return best

onethird=F(1,13)
print("FAR-ELEMENT MEASURE LEMMA reach on multi-killer configs:")
print("C = S minus largest outlier L; lemma applies iff L > 1/w_max(G_C(1/13)).\n")
print(f"{'S':34s} | L | M(C) | w_max(G_C) | 1/w_max | L>1/wmax? | M(S)")
print("-"*104)
cases=[
 [*range(1,12),13,84],       # k=11 extremal, M=7/89
 [*range(1,11),11,13,84],    # 
 [*range(1,11),22,13,84],    # k=10, M=2/23
 [*range(1,12),13,168],
 [*range(1,12),26,84],
 [*range(1,10),10,11,13,84], # k=9-ish
 [*range(1,12),14,156],      # cover 14 by 14? 14 is outlier; 156=lcm(12,13)
]
for S in cases:
    S=sorted(set(S))
    if len(S)!=13: continue
    L=max(S); C=[x for x in S if x!=L]
    MC=M_exact(C,150)
    wmax=safe_arcs_width_max(C)
    inv = F(1,1)/wmax if wmax>0 else None
    applies = (wmax>0 and F(L)>inv)
    MS=M_exact(S,min(2*L,300))
    invs = f"{float(inv):.2f}" if inv else "inf"
    print(f"{str(S):34s} | {L:3d} | {float(MC):.4f} | {float(wmax):.5f} | {invs:>7} | {str(applies):5s} | {float(MS):.5f}={MS}")
print()
print("If 'L>1/wmax'=True for the extremals, the measure lemma PROVES M(S)>=1/13 directly")
print("(no finite check, no undershoot) for those. Remaining: M(C)=1/13 (dilated AP, shallow")
print("witness) + any config where L<=1/wmax (small L & small wmax = near-tight C).")
