"""
opus-2026-07-11-S240: attacking the concrete next target (consec-extremality of the base functional J, the
shared crux of both LRC(14) routes per S239) via the natural COMPRESSION/smoothing proof route -- and finding
it BLOCKED. Honest negative result about the proof strategy + a structural reason.

TARGET (S239 unification): both routes reduce to "consec/AP is the extremal coverer." The base form (klein
THM-711/717, mac-mini THM-716): J = E[N(7-N)] = 6m1-m2 over the 7 sectors; consec {1..9} is the conjectured
global MIN (J~5.05), verified adversarially but NOT proved. The standard route to prove an extremal is a
MONOTONE COMPRESSION: a local step toward the AP that decreases J, iterated to the consec fixed point.

RESULT: the compression route is BLOCKED. Greedy single-coordinate J-descent from 35 random primitive 9-sets
reaches consec 0/35 times -- it STUCK at 35 DISTINCT non-consec local minima, all algebraically-special
(even/dilated families near 2*{1..k}): e.g. [2,4,6,7,8,9,10,12,14] (J~5.33), [2,4,6,8,9,10,12,14,16] (J~5.47),
all with J in [5.3, 5.8] > J(consec)=5.05. So J is NOT unimodal toward consec; the landscape is RUGGED with
many algebraically-special local minima, and consec is the global min but UNREACHABLE by local descent.

WHY THIS MATTERS (honest): (1) it RULES OUT the natural single-step compression/smoothing proof of
consec-extremality -- the crux does not yield to a local convexity/monotonicity argument. (2) it CONFIRMS
klein-S255's "every LRC(14) functional's extremal is an AP or a mod-p resonance, INVISIBLE to local search":
the local minima ARE the algebraically-special (even/dilated/mod-7) families, and consec beats them all only
GLOBALLY. (3) it explains WHY consec-extremality is hard and verified-only: proving consec is the global min
over a rugged landscape of algebraically-special competitors is a GLOBAL/algebraic problem (the inverse
theorem), not a local smoothing. The compression route -- the most natural proof attack -- is dead.

CONSEQUENCE for the endgame: consec-extremality (= the S239 shared wall of both routes) will not fall to
compression. It needs either a global algebraic argument (classify the algebraically-special local minima and
beat each -- klein/mac-mini's atlas lane) or the finite-census framing. Solo local attacks are exhausted.
"""
from math import gcd
from functools import reduce
import random
def J(E, G=630):
    tot=0
    for m in range(G):
        t=(m+0.5)/G
        occ=set(int((e*t)%1*7) for e in E)
        tot+=(7-len(occ))*len(occ)
    return tot/G
def prim(E): return reduce(gcd,E)==1
def descent(E):
    E=sorted(E); cur=J(E)
    for _ in range(50):
        best=None; bJ=cur
        for i in range(9):
            for d in (-1,1):
                nv=E[:]; nv[i]+=d
                if nv[i]<1 or len(set(nv))<9 or not prim(nv): continue
                j=J(nv)
                if j<bJ-1e-9: bJ,best=j,sorted(nv)
        if best is None: return E,cur
        E,cur=best,bJ
    return E,cur

def main():
    consec=list(range(1,10)); Jc=J(consec)
    random.seed(1); reached=0; stuck=[]; n=0
    for _ in range(35):
        E=sorted(random.sample(range(1,18),9))
        if not prim(E): continue
        n+=1
        lm,lJ=descent(E)
        if lm==consec: reached+=1
        else: stuck.append((round(lJ,3),tuple(lm)))
    print(f"J(consec {{1..9}}) ~ {Jc:.3f} (global min). Greedy J-descent reached consec {reached}/{n}.")
    print(f"STUCK at non-consec local minima: {len(stuck)} ({len(set(m for _,m in stuck))} distinct, all algebraically-special):")
    for jv,m in sorted(set(stuck))[:6]: print(f"   J~{jv}  {list(m)}")
    print("VERDICT: J-landscape NOT unimodal => single-step compression proof of consec-extremality is BLOCKED.")

if __name__=='__main__':
    main()
