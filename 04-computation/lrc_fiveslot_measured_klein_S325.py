#!/usr/bin/env python3
"""
lrc_fiveslot_measured_klein_S325.py -- klein-2026-07-18-S325
Owner: work the higher-order direction with opus's five-slot ledger (THM-1026).

THM-1026 prices the level-5 truncation B5 = 1 - S1 + S2 - S3 + S4 - S5 at EQUIDISTRIBUTED slot values
S_k = C(13,k)(1/7)^k, where it clears at +0.1221, and lists S3(upper)/S4(lower)/S5(upper) as the three
open slots. I measured the slots on real families. THE PREMISE DOES NOT HOLD.

S_k VALIDATED EXACTLY: full inclusion-exclusion sum_k (-1)^k S_k reproduces the true uncovered measure
to 0.00e+00 on AP / deep well / random covering. So what follows is not a computation artefact.

MEASURED S_k / EQUIDISTRIBUTED  (the factor the open-slot bounds would have to accommodate):
                          S3      S4       S5
    AP {1..13}           5.49x   32.47x   211.02x
    deep well            4.98x   26.44x   156.04x
    GW                   4.92x   28.00x   174.12x
    {1..12,14} artefact  5.61x   32.57x   207.94x
    random covering      1.2x     3.0x     15.2x     <- even the MOST generic family is 15x off at S5
Equidistribution is not a mild idealization at the higher slots; it is off by 1-2 orders of magnitude.

CONSEQUENCE -- B5 IS NEGATIVE EVERYWHERE, INCLUDING WHERE IT SHOULD BE EASIEST:
    AP -9.72 | deep well -6.87 | GW -7.86 | {1..12,14} -9.49 | random covering -0.53 .. -0.58
ODD-TRUNCATION SCAN (odd B_m is a valid LOWER bound on the uncovered measure):
    AP        : B1 -0.86, B3 -3.23, B5 -9.72, B7 -9.20, B9 -2.48, B11 -0.13, B13 +0.0000 (= true 0)
    random cov: B1 -0.86, B3 -1.09, B5 -2.78, B7 -2.37, B9 -0.52, B11 +0.089, B13 +0.1224 (= true)
The first level that clears on a GENERIC covering family is B11, not B5. Level 5 is far too shallow.

MECHANISM (chains to klein-S324): the danger events are positively correlated (Var(X) = +77% over
independence for the AP, +5% for random covering). That correlation COMPOUNDS with order, so S_k blows up
(15x-211x at k=5) and the alternating truncation stays dominated by the odd terms until very high order.

WHAT THIS DOES NOT KILL: non-alternating certificates. The interval-survival tail lemma (THM-1004/1015) is
a direct measure argument with no inclusion-exclusion cancellation, and it does not degrade this way --
that is why it closed the clustered stratum additively.
"""
import numpy as np
from math import comb
NG=1<<20; tt=np.arange(NG)/NG
def slots(S):
    X=np.zeros(NG,dtype=np.int32)
    for v in S:
        fr=(v*tt)%1.0; X += (np.minimum(fr,1-fr)<1.0/14)
    bc=np.bincount(X)/NG
    Sk=[sum(bc[n]*comb(n,k) for n in range(len(bc))) for k in range(len(S)+1)]
    return Sk, bc[0]
def trunc(Sk,m): return sum((-1)**k*Sk[k] for k in range(m+1))
if __name__=="__main__":
    for nm,S in [("AP",list(range(1,14))),("deep well",list(range(1,13))+[182]),
                 ("random cov",[3,5,7,8,9,11,13,16,22,27,35,44,52])]:
        Sk,L=slots(S)
        print("%-11s B5=%+8.4f  B11=%+8.4f  full=%+.6f  true=%.6f"
              %(nm,trunc(Sk,5),trunc(Sk,11),trunc(Sk,13),L))
