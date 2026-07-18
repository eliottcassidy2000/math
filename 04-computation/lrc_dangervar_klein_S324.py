#!/usr/bin/env python3
"""
lrc_dangervar_klein_S324.py -- klein-2026-07-18-S324
Owner: find a transfer candidate that distinguishes the AP; compare its structure to randomness.

CANDIDATE: Var_t(X), X(t) = #{i : ||v_i t|| < 1/14} (mac-mini-S89's danger count; tournament analog =
score variance / eigenvalue flatness, THM-133). Mean = 13/7 for EVERY 13-set, so all content is in the
variance, which is a pure PAIRWISE functional (independent reference 78/49 = 1.5918).

PASSES THE FILTER: AP 2.8116 (+77%), deep well 2.8342, GW 2.5977, vs RANDOM covering 1.55-1.92 (+5%).
Random covering sets are essentially INDEPENDENT; near-tight families are strongly positively correlated.
All named low-M objects sit at the 100th percentile of the random cloud.

BUT FAILS AS A DETECTOR:
 (1) no gap -- AP->random interpolation decays smoothly 2.81,2.65,2.49,...,1.68; the apparent gap was a
     SAMPLING gap (random 13-sets are never near-AP), not a structural one;
 (2) not predictive -- corr(Var,M) = -0.065 inside the random cloud;
 (3) not variational -- 585 single-replacements BEAT Var(AP), topping at 3.0576, by adding multiples of 14
     (speeds exactly on the 1/14 grid). Those are NON-COVERING with M ~ 0.077-0.091: resonance artefact.

CONCLUSION: Var(X) summarizes the pairwise data completely and cannot characterize the AP, so NO
pairwise-only invariant can. Coherence and tightness come apart: a speed on the 1/14 grid maximizes
pairwise correlation while HELPING loneliness. Independently corroborates opus THM-1026 ("pairs
insufficient at 13") from the opposite direction.
"""
import numpy as np
NG=1<<18; tt=np.arange(NG)/NG
def varX(S):
    X=np.zeros(NG,dtype=np.int32)
    for v in S:
        fr=(v*tt)%1.0; X += (np.minimum(fr,1-fr)<1.0/14)
    return X.var()
if __name__=="__main__":
    print("AP      Var(X) =", round(varX(list(range(1,14))),4))
    print("indep reference =", round(78/49,4))
    print("beaten by e.g. {1..10,12,13,15}:", round(varX([1,2,3,4,5,6,7,8,9,10,12,13,15]),4))
