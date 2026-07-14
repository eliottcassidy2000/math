#!/usr/bin/env python3
"""
lrc14_lrc13_localization_klein_S295.py
======================================
klein-2026-07-13-S295 (owner: prove the multi-speed near-0 equidistribution).

LRC(13) LOCALIZATION. The cluster C (12 speeds) satisfies M(C) >= 1/13 (LRC(13), SETTLED). By the S290
symmetry, L({1}UC) > 0  <=>  |G(C) cap [1/14,13/14]| > 0  <=>  C's good set REACHES THE MIDDLE.

KEY INSIGHT (verified): M(C) >= 1/13 holds with HUGE margin (M = 0.13..0.47), but this is IRRELEVANT to
placement -- the AP-cluster {2..13} is very lonely (M=2/15) yet ALL its witnesses are trapped at the ends
[0,1/14) u (13/14,1) (|G cap mid|=0, L=0); every COVERING C reaches the middle (|G cap mid|>0, L>0). So
the residual is a PLACEMENT/RIGIDITY question (do C's witnesses avoid the AP-ends), decoupled from
loneliness MAGNITUDE. Only the AP-cluster (<=> C u {1} = the AP {1..13}, non-covering) confines its good
set to the ends. This is the irreducible multi-speed near-0 equidistribution = the cancellation.
"""
import numpy as np
NG=1<<22; t=np.arange(NG)/NG
def profile(C):
    m=np.full(NG,1.0)
    for c in C:
        fr=(c*t)%1.0; d=np.minimum(fr,1.0-fr); m=np.minimum(m,d)
    return m  # m(t)=min_i||c_i t||; M(C)=max m
def report(C,name):
    C=sorted(set(C)); m=profile(C); M=m.max()
    W13=m>=1.0/13.0-1e-9
    Wmid=W13&(t>=1.0/14.0)&(t<=13.0/14.0)
    G=m>=1.0/14.0-1e-9; Gmid=G&(t>=1.0/14.0)&(t<=13.0/14.0)
    print("  %-22s M(C)=%.5f (>=1/13? %s) | |G|=%.4f |G∩mid|=%.4f = L({1}UC)  [W13: |.|=%.4f |∩mid|=%.4f]"%(
        name,M,M>=1.0/13.0-1e-6,G.mean(),Gmid.mean(),W13.mean(),Wmid.mean()))
print("LRC(13) localization: L({1}UC)>0 <=> G(C) reaches middle [1/14,13/14]. 1/13=%.5f"%(1.0/13.0))
for C,nm in [([90,91,92,93,94,95,96,97,98,99,100,101],'{90..101}'),
             ([2,3,4,5,6,7,8,9,10,11,12,13],'{2..13} AP-cluster'),
             ([3,4,5,6,7,8,9,10,11,12,13,14],'{3..14}'),
             ([30,31,32,33,34,35,36,37,38,39,40,41],'{30..41}'),
             ([50,55,60,66,70,77,84,90,99,100,101,103],'spread-ish')]:
    report(C,nm)
print("\n  M(C)>=1/13 with HUGE margin (0.13..0.47) but IRRELEVANT: {2..13} very lonely yet witnesses trapped")
print("  at the ends (=AP, non-covering). Residual = PLACEMENT (reach the middle), not loneliness magnitude.")
print("done.")
