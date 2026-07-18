#!/usr/bin/env python3
"""
lrc14_rigidity_radius4_klein_S315.py
====================================
klein-2026-07-17-S315 (owner: prove radius 4, and more).

HONEST OUTCOME: radius 4 is NOT proved. Recorded here: the sizing, the counting DEAD END, the stress
evidence, and the method's hard ceiling.

SIZING. r=4 tail needs sum 1/w_i < 9L/4 (1-8*delta = 9/25); B has 8 speeds so M(B) >= 1/9 by LRC(9).
  W_joint = 16/(9 L_min) = 148.1 (worst quadruple (3,4,5,12), L_min=0.012) -> box ~200
  -> finite core = 495 * C(188,4) ~ 2.4e10 tuples. Infeasible to enumerate.

THE DEAD END (do not retry). Reframe the bitmask table as a HITTING-SET count: an uncovered r-tuple must
hit every B_c = box \ S_c, so N <= sum_i deg(w_i) <= (sum of top-r degrees). Hence
      (top-r degree sum) < N   ==>   every r-tuple is covered, whole finite core in ONE inequality.
FAILS at every r>=2 and every pool (QMAX=60,120,200,300; N=551..13699). e.g. r=4,QMAX=300: 3361 vs 2369.
WHY: deg(w)/N averages ~0.16 but its MAX reaches ~0.63 (highly composite w kill every candidate with q|w),
so degrees are far too correlated for a union bound -- and a bigger pool makes the ratio WORSE.

EVIDENCE (not proof). Adversarial stress over all 495 quadruples: hardest-26 w by kill-degree taken
exhaustively (C(26,4) each) + random tuples across the box = 7,548,750 4-tuples, ZERO uncovered.

HARD CEILING (proved). The tail lemma needs 1-2 r delta > 0, i.e. r < 1/(2 delta) = 6.25:
  r=1..6 give 1-2r*delta = +.84,+.68,+.52,+.36,+.20,+.04  (OK);  r=7 gives -0.12 (DEAD).
So this method reaches radius <= 6 and provably cannot reach radius 7.
"""
import numpy as np, itertools, random
from math import gcd
AP=list(range(1,13)); num,den=2,25; WMAX=200; LO=13; QMAX=60
def tables(QMAX=QMAX):
    cands=[(p,q) for q in range(2,QMAX+1) for p in range(1,q//2+1) if gcd(p,q)==1]
    P=np.array([p for p,_ in cands]); Q=np.array([q for _,q in cands])
    ws=np.arange(LO,WMAX+1)
    R=(np.outer(P,ws))%Q[:,None];  S=np.minimum(R,Q[:,None]-R)*den>=num*Q[:,None]
    vs=np.arange(1,13)
    Rv=(np.outer(P,vs))%Q[:,None]; Sv=np.minimum(Rv,Q[:,None]-Rv)*den>=num*Q[:,None]
    return S,Sv,ws
if __name__=="__main__":
    S,Sv,ws=tables()
    print("counting criterion (top-r deg sum < N) per radius:")
    for r in (2,3,4,5,6):
        ok=True
        for rem in itertools.combinations(AP,r):
            keep=np.ones(Sv.shape[0],bool)
            for v in AP:
                if v not in rem: keep &= Sv[:,v-1]
            C=S[keep]; N=C.shape[0]
            if N==0: ok=False; break
            if np.sort((~C).sum(axis=0))[-r:].sum()>=N: ok=False; break
        print("   r=%d : %s"%(r,"holds" if ok else "FAILS (dead end)"))
    print("\nhard ceiling: r < 1/(2*delta) = %.2f"%(1/(2*(2/25))))
