#!/usr/bin/env python3
"""
lrc14_windowed_overlap_klein_S294.py
====================================
klein-2026-07-13-S294 (owner: prove the windowed overlap on [0,1/14)).

W = |bad_c ∩ bad_c' ∩ [0,1/14)|, the LOCAL version of THM-739 the one-interval bound (S292) needs.

CORRECTION 2026-09-06: the original explanatory resonance formula was false:
it omitted containment clipping, multiple resonant partners and clipping to
the observation window. See the corrected microscopic addendum in THM-739
(pairwise-coprime-bad-overlap slug). The numerical body below is unchanged
and computes literal sampled masks; its finite observations are not an
all-height window theorem. Printed historical conclusions retain that finite
scope and do not supply a proof of near-origin asymptotics.
"""
import numpy as np
from math import gcd
NG=1<<24; THR=1.0/14.0; t=np.arange(NG)/NG
def badmask(c):
    fr=(c*t)%1.0; d=np.minimum(fr,1.0-fr); return d<THR
def B2(x):
    x=x%1.0; return x*x-x+1.0/6.0
def overlap_full(c,cp): return 1.0/49.0 + (1.0/(c*cp))*(B2((cp-c)/14.0)-B2((cp+c)/14.0))
lo=int(round(NG/14.0))
def W_direct(c,cp): return (badmask(c)&badmask(cp))[:lo].sum()/NG
print("W = |bad_c ∩ bad_c' ∩ [0,1/14)|.  bulk=(1/14)*overlap_full ≈ 1/686=%.6f"%(1.0/686.0))
print("%10s %5s %12s %12s %8s  %s"%("(c,c')","|c-c'|","W_direct","bulk","W/bulk","regime"))
for c,cp in [(3,5),(11,13),(99,101),(50,99),(23,45),(90,101),(100,171),(200,313),(3,200),(17,101)]:
    if gcd(c,cp)!=1: continue
    W=W_direct(c,cp); bulk=overlap_full(c,cp)/14.0
    reg="CLOSE (cluster) -> W large" if abs(cp-c)<=5 else ("far -> W~bulk" if abs(cp-c)>50 else "mid")
    print("  (%3d,%3d) %5d %12.8f %12.8f %8.2f  %s"%(c,cp,abs(cp-c),W,bulk,W/bulk,reg))
print()
print("=> W is LARGE for close speeds (clusters): pairwise near-0 overlap does NOT decorrelate;")
print("   the 2-speed one-interval refinement fails for clusters; near-0 equidistribution is MULTI-speed.")
print("done.")
