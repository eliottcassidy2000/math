#!/usr/bin/env python3
"""
[FIX S20]: removed the buggy `if gcd(a,q)!=1: continue` -- it skipped witnesses at
denominators that DIVIDE a pairwise sum/diff (e.g. q=11 via 22=4+18 as the non-coprime 2a/22),
underestimating M.  Now checks ALL a. See MISTAKE log.
mac-mini-2026-07-06-S17 (HYP-4452) -- the LEAVE-ONE-OUT ALIGNMENT lens for the density floor.

NECESSARY CONDITION FOR COVERING: if S covers at beta (M(S)<beta), then for EVERY j,
   Safe(S\\{v_j}, beta)  is a subset of  A_j = {t: ||v_j t|| < beta}
i.e. the dropped runner's danger arcs must contain the ENTIRE hole of the 11-subfamily
(which is nonempty since an 11-family has M >= 1/12 > 2/25).  Covering IS this nesting.

WHY IT IS THE QUANTITATIVE (n-specific) LENS: A_{v_j} is v_j arcs of width 2beta/v_j at
positions a/v_j; Safe(S\\{v_j}) is holes of width ~ gap-width = 1/(n(2n-1)).  Containment
needs each thin hole to nest inside an arc at an ALIGNED position -- a lattice rigidity
the AP's harmonic arcs {a/k} achieve.  The hole width shrinks with n (1/91 -> 1/325), so
at n=13 only the AP-lattice aligns => the floor.  (The second gap is n-SPECIFIC: nonempty
at n=7 {1,5,6,11,16,17} M=5/33, empty at n=13 -- HYP-4442.)
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations

def arcs_danger(v, beta):
    out=[]; j=0
    while F(j-beta,v)<1:
        lo=max(F(j-beta,v),F(0)); hi=min(F(j+beta,v),F(1))
        if lo<hi: out.append((lo,hi))
        j+=1
    return out
def merge(iv):
    iv=sorted(iv); out=[]
    for lo,hi in iv:
        if out and lo<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],hi))
        else: out.append((lo,hi))
    return out
def safe_set(S, beta):
    d=merge([a for v in S for a in arcs_danger(v,beta)])
    safe=[]; prev=F(0)
    for lo,hi in d:
        if lo>prev: safe.append((prev,lo))
        prev=max(prev,hi)
    if prev<1: safe.append((prev,F(1)))
    return safe
def meas(iv): return sum((b-a for a,b in iv), F(0))
def contained(A,B):
    for lo,hi in A:
        p=lo
        for bl,bh in sorted(B):
            if bl<=p<bh: p=bh
            if p>=hi: break
        if p<hi: return False
    return True
def Mfast(S):
    S=sorted(set(S)); Q=set()
    for v in S: Q.add(2*v)
    for a,b in combinations(S,2): Q.add(a+b); Q.add(abs(a-b))
    Q.discard(0); best=F(0)
    for q in Q:
        for a in range(1,q):
            mn=min(min((v*a)%q,q-((v*a)%q)) for v in S)
            v=F(mn,q)
            if v>best: best=v
    return best

def test(name,S,beta):
    S=sorted(S); M=Mfast(S)
    print(f"{name}: S={S}  M={M}  covers(M<beta={beta})? {M<beta}")
    if M>=beta:
        print(f"   (does not cover; safe measure={float(meas(safe_set(S,beta))):.4f})"); return
    allok=True
    for j in range(len(S)):
        sub=S[:j]+S[j+1:]; Ssafe=safe_set(sub,beta); Aj=arcs_danger(S[j],beta)
        ok=contained(Ssafe,Aj); allok &= ok
        print(f"   drop {S[j]:>3}: Safe(sub) meas={float(meas(Ssafe)):.5f}, nests in A_{S[j]}? {ok}")
    print(f"   ALL leave-one-out holes nest: {allok}  (the covering-alignment mechanism)")

if __name__=="__main__":
    test("AP {1..12} (n=13)", list(range(1,13)), F(2,25))
    test("n=7 gap member", [1,5,6,11,16,17], F(2,13))
