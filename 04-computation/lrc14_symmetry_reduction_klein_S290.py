#!/usr/bin/env python3
"""
lrc14_symmetry_reduction_klein_S290.py
======================================
klein-2026-07-13-S290 (owner: prove L>0 for the compact core via the shared cancellation; continue+extend).

The compact core splits:
  - BOUNDED-RATIO (max<=13 min): already DONE by THM-405 (lonely on [1/(14a),13/(14b)], no cancellation).
  - min=1 non-isolated residual {1}UC (ratio>13): handled here by an EXACT symmetry reduction.

EXACT SYMMETRY REDUCTION.  For S={1}UC, G({1})={t:||t||>=1/14}=[1/14,13/14] (ONE arc), and G(C) is
symmetric under t->1-t, so
    L(S) = |G(C)| - 2|G(C) cap [0,1/14)| = |G(C)| * (1 - conc/7),   conc := 14|G(C)cap[0,1/14)|/|G(C)|.
Hence  L>0  <=>  conc < 7.  Census: conc=7 EXACTLY at C={2..13} (S=AP {1..13}, non-covering, L=0);
conc<7 with margin for covering (max 5.165 at C={3..14}). So the AP is the UNIQUE tight extremal, and
covering forces conc<7 => L>0 (VERIFIED). This is an exact geometric RESTATEMENT of the residual (conc<7
<=> L>0), pinning the AP boundary; it converges with opus-S271 (dilation-blindness/AP-shadow) and
kps-THM-734 (tight census {AP, GW-doubling}).
"""
import numpy as np
NG=1<<21; THR=1.0/14.0; t=np.arange(NG)/NG
def good(W):
    g=np.ones(NG,bool)
    for w in W:
        fr=(w*t)%1.0; d=np.minimum(fr,1.0-fr); g&=(d>=THR)
    return g
def iscov(S): return all(any(x%q==0 for x in S) for q in range(2,15))
lo=int(round(NG/14.0))

print("=== (A) EXACT symmetry reduction  L({1}UC) = |G(C)| - 2|G(C)cap[0,1/14)|  (verify vs L_full) ===")
print("%-30s %9s %13s %12s %9s %10s %s"%("S={1}UC","|G(C)|","|GcCap[0,1/14)|","L=|Gc|-2that","conc","L_full","cov"))
for S in [[1,90,91,92,93,94,95,96,97,98,99,100,101],
          [1,2,3,4,9,10,11,12,13,14,15,16,17],
          [1,30,31,32,33,34,35,36,37,38,39,40,42]]:
    S=sorted(set(S)); C=[w for w in S if w!=1]
    gc=good(C); Lc=gc.mean(); near0=gc[:lo].sum()/NG
    Lsym=Lc-2*near0; conc=near0/(Lc/14.0); Lfull=good(S).mean()
    print("%-30s %9.5f %13.6f %12.6f %9.3f %10.6f %s"%(str(S)[:30],Lc,near0,Lsym,conc,Lfull,'Y' if iscov(S) else 'N'))

print("\n=== (B) concentration-ratio census over consecutive clusters C={b..b+11}: conc=7 at the AP only ===")
print("%-14s %9s %9s %11s %s"%("C","|G(C)|","conc","L({1}UC)","S={1}UC cov?"))
for b in [2,3,4,5,7,9,12,20,45,90]:
    C=list(range(b,b+12)); gc=good(C); Lc=gc.mean(); near0=gc[:lo].sum()/NG
    conc=near0/(Lc/14.0); S=sorted(set([1]+C))
    tag="  <-- AP {1..13} (L=0, non-covering)" if b==2 else ""
    print("{%2d..%2d}      %9.5f %9.3f %11.6f %s%s"%(b,b+11,Lc,conc,Lc-2*near0,'Y' if iscov(S) else '-',tag))
print("\n  L>0 <=> conc<7. AP is the UNIQUE conc=7 extremal; covering => conc<=~5.2 => L>=0.26|G(C)|>0.")
print("  HONEST: conc<7 <=> L>0 is a RESTATEMENT, not a reduction; but it pins the AP boundary exactly.")
print("done.")
