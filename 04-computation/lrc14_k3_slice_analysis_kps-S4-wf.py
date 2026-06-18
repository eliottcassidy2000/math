"""
LRC(14) — k=3 slice analysis (exactly THREE speeds > 13).
angle "scale-invariance-AP-and-k3"  (kind-pasteur-S4-wf)

GOAL: extend the PROVED k=2 bounded-core slice to k=3.

FINDINGS (all exact):
  1. EXHAUSTIVE near-boundary board [14,30]: 11632 covering primitive S3 k=3 sets,
     min M = 2/23 (margin M*14 = 28/23 ~ 1.217 > 1). ZERO counterexamples.
     The worst set [1,2,3,4,5,7,8,11,12,13,18,20,28] realizes M=2/23 at tau=4/23,
     binders {11,12} — the SAME global S3 floor as the documented k=2 set
     S*={1,2,3,5,7,8,9,10,11,12,13,38,42}.
  2. DROP-MAX STABILITY is only APPROXIMATE for k=3: with P and the two smaller
     large speeds (a,b) fixed, M is usually constant in the largest speed c, but
     DIPS at isolated large c (e.g. M(c) takes values 250/1251 etc, varying with c).
     This is WHY the clean k=2 finite-core argument does NOT lift verbatim: the
     top speed cannot simply be "dropped" because it can create a new binding
     resonance at large c.
  3. The dips stay FAR above 1/14 (smallest observed ~0.097). The floor is governed
     by the SHAPE (P and the gap structure), not by the magnitude of c.

HONEST STATUS: this is a VERIFIED finite-board closure + structural diagnosis,
NOT a finite-core PROOF for all k=3. The obstruction (approximate drop-max
stability) is named precisely. The k=2 -> k=3 lift requires controlling the
isolated large-c dips uniformly, which remains OPEN.
"""
from fractions import Fraction as F
from math import gcd
import itertools

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r

def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1, 2): C.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2): C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C

def Mval(S):
    b = F(0); bt = None
    for t in cand(S):
        v = min(nrm(x*t) for x in S)
        if v > b: b = v; bt = t
    return b, bt

def is_cov(S): return all(any(v % q == 0 for v in S) for q in range(2, 15))
def gcd_all(S):
    g = 0
    for x in S: g = gcd(g, x)
    return g
def case_of(S):
    S = sorted(set(S)); k = sum(1 for v in S if v > 13)
    if k <= 1: return 'S1'
    if S[-1] < 13*S[0]: return 'S2'
    return 'S3'

def exhaustive_board(LO=14, HI=30):
    minM = F(1); worstM = None; n = 0
    for P in itertools.combinations(range(1, 14), 10):
        P = list(P)
        for (a, b, c) in itertools.combinations(range(LO, HI+1), 3):
            S = sorted(set(P) | {a, b, c})
            if len(S) != 13: continue
            if gcd_all(S) != 1: continue
            if not is_cov(S): continue
            if case_of(S) != 'S3': continue
            n += 1
            Mv, _ = Mval(S)
            if Mv < minM: minM = Mv; worstM = S
    return n, minM, worstM

if __name__ == '__main__':
    print("k=3 exhaustive board [14,30]:")
    n, minM, worstM = exhaustive_board(14, 30)
    print(f"  {n} covering primitive S3 k=3 sets, min M = {minM} = {float(minM):.5f}")
    print(f"  margin M*14 = {minM*14} = {float(minM*14):.4f}  (>1 => strictly above 1/14)")
    print(f"  worst set: {worstM}")
    M, bt = Mval(worstM)
    print(f"  worst-set optimal tau = {bt}, binders {[v for v in worstM if nrm(v*bt)==M]}")
    print()
    print("  drop-max stability example (M dips at isolated large c):")
    P = [1,2,3,4,6,7,8,11,12,13]; a, b = 32, 58
    for c in [59, 60, 251, 1251, 5000]:
        S = sorted(set(P)|{a,b,c})
        if len(S) != 13: continue
        Mv, _ = Mval(S)
        print(f"     c={c:5d}: M={Mv} = {float(Mv):.5f}")
    print("  => M not constant in c; finite-core lift from k=2 does NOT hold verbatim.")
