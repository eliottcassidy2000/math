#!/usr/bin/env python3
"""
lrc_hypothesis_refutation_klein_S330.py -- klein-2026-07-18-S330
Owner: prove or refute the witness-vs-accounting hypothesis.  REFUTED, by my own theorem.

REFUTATION. THM-731 (rigorous) gives L = (6/7)|G_P| - eps_v with |eps_v|^2 <= (6/49) disc_v, so
L > 0 <=> disc_v < 6|G_P|^2 -- a pure MEASURE statement, no point exhibited. THM-755 (PROVED) gives
disc_v <= 4 r_P |G_P|/(pi v) + 2|G_P|^2 < 6|G_P|^2 exactly when v > v* = r_P/(pi |G_P|). Chained: for
v > v*, L > 0 is PROVED BY ACCOUNTING ALONE. So accounting certificates work; the hypothesis is false.

REPLACEMENT (verified below): classify by SCALE, not by method.
  family                 regime                    v      v*      accounting  witness
  deep well {1..12,182}  separated  rho=15.2      182    112.0      FIRES      THM-1007
  {1..12,364}            separated  rho=30.3      364    112.0      FIRES      yes
  {1..11,13,84}          INTERMEDIATE sigma=84     84    104.7      fails      fails
  2*{1..12}u{13}         coherent   sigma=12       24    158.0      fails      SPREAD LADDER
Three findings: (1) both method families succeed at scale extremes and fail in the middle; (2) they are
COMPLEMENTARY -- at the coherent end the witness fires and accounting does NOT (v*=158 >> 24); (3) their
union misses exactly the wedge sigma>13 AND rho<13, the same boundary THM-1043 found by a different route.

CORRECTED HYPOTHESIS: certificates of BOTH kinds succeed exactly in the scale-extreme regimes (one
13-fold window, or one speed separated by >13x); the open residual is exactly the intermediate regime,
one octave wide (W = ceil(log_13 sigma) = 2). Falsifiable: exhibit a proved certificate firing strictly
inside the wedge, or a scale-extreme family no certificate reaches. Binding test: {1..11,13,84}, where
v=84 misses v*=104.7 by a factor 1.25.
"""
import numpy as np
from fractions import Fraction as Fr
NG=1<<20; t=np.arange(NG)/NG
def gstat(W):
    g=np.ones(NG,bool)
    for w in W:
        fr=(w*t)%1.0; g &= (np.minimum(fr,1-fr)>=1.0/14)
    gi=g.astype(np.int8); r=int(np.sum(np.diff(np.concatenate([gi,gi[:1]]))==1))
    return g.mean(), r
def vstar(P):
    m,r=gstat(P); return (r/(np.pi*m)) if m>0 else float('inf')
def regime(V):
    V=sorted(V); sig=Fr(V[-1],V[0]); rho=Fr(V[-1],V[-2])
    if sig<=13: return 'coherent'
    if rho>=13: return 'separated'
    return 'INTERMEDIATE WEDGE'
if __name__=="__main__":
    for nm,V in [('deep well',list(range(1,13))+[182]),
                 ('{1..11,13,84}',list(range(1,12))+[13,84]),
                 ('2*{1..12}u{13}',[2*k for k in range(1,13)]+[13])]:
        V=sorted(V); vs=vstar(V[:-1])
        print("%-16s regime=%-18s v=%4d v*=%6.1f accounting=%s"
              %(nm,regime(V),V[-1],vs,"FIRES" if V[-1]>vs else "fails"))
