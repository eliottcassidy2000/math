#!/usr/bin/env python3
"""
lrc14_sector_domination_macmini_0618s7b.py  (mac-mini-2026-06-18-S7, ANGLE B)

THE RIGOROUS PIECE.  POINTWISE DOMINATION LEMMA (proved trivially):
  If E subset {0,1,...,N}, then for EVERY x, {sigma_e(x): e in E} subset {sigma_e(x): e in 0..N}.
  Hence  S7(E) subset S7({0..N}), so  meas(S7(E)) <= meas(S7(AP_{N+1})).

GOAL: leverage this + meas(S7) SCALE-INVARIANCE (meas(S7(dE))=meas(S7(E))) to bound meas(S7(E))
for ARBITRARY primitive |E|=k, k=8..11 (dangerous rows), and see how close to cap_k we get.

OBSERVATION.  Any |E|=k with min 0 has max(E)=N>=k-1.  E subset {0..N} so
  meas(S7(E)) <= meas(S7(AP_{N+1})).
meas(S7(AP_m)) is INCREASING in m (more residues, larger cover): 0,0,0,0,0,.147,.327,.416,.504,
.581,.624,.676,.705 for m=1..14. So a LARGE-span E gets a WEAK bound (AP_{N+1} is big).
=> the bound is only useful for SMALL span. For span exactly k-1 (E=AP_k) it is sharp.
=> the hard cases are MEDIUM/LARGE span E -- exactly the AP-rich residual.

But we can do better using a SECOND structural fact: how many of the AP_{N+1} indices does E
actually occupy?  If E omits index e from {0..N}, then at theta where ONLY index e would have
supplied the 7th residue, E fails to cover. Quantify the loss per omitted index.

ALSO: combine with the dual cap.  We need meas(S7(E)) <= cap_k. The pointwise lemma gives
meas(S7(E)) <= meas(S7(AP_{N+1})). When is meas(S7(AP_{N+1})) <= cap_k?  cap_k is FIXED (the
P-side budget). meas(S7(AP_m)) crosses cap_8=0.3815 around m=8-9. So pointwise domination
CLOSES the cap exactly for E with span N <= N*(k) where meas(S7(AP_{N+1})) <= cap_k.

THIS SCRIPT:
 (1) tabulate meas(S7(AP_m)) and find, for each k, the max span N*(k) s.t. AP_{N+1} <= cap_k.
     => pointwise lemma certifies ALL primitive E of span <= N*(k).  Residual = span > N*(k).
 (2) For the residual (span > N*(k)), is it FINITE mod scale? Use scale-invariance: divide by
     gcd. The residual is {primitive E, |E|=k, span in (N*(k), ???]}. Bounded above by what?
     If unbounded, the pointwise lemma alone does NOT close it -- need the relation-height tail.
 (3) Quantify: what FRACTION of the dangerous-box shapes does pointwise domination certify?
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
sys.stdout.reconfigure(line_buffering=True)

def measS7(E):
    E=sorted(set(E)); Enz=[e for e in E if e!=0]
    bps=set([F(0),F(1)])
    for e in Enz:
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); total=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        res=set(int(7*e*xm)%7 for e in E)
        if len(res)==7: total+=x1-x0
    return total

def measS7_AP(m):
    return measS7(tuple(range(m)))

# cap_k from canon (min_{|P|=13-k} meas(G_P))
cap = {8:F(2243,5880), 9:F(1979,4004), 10:F(55,91), 11:F(66,91), 12:F(6,7), 13:F(1)}

print("="*92)
print("(1) meas(S7(AP_m)) and the max certifiable span N*(k): pointwise lemma closes span<=N*(k)")
print("="*92)
apvals={m:measS7_AP(m) for m in range(1,22)}
print("  m:    " + " ".join(f"{m}" for m in range(7,22)))
print("  AP_m: " + " ".join(f"{float(apvals[m]):.3f}" for m in range(7,22)))
print()
for k in sorted(cap):
    ck=cap[k]
    # E subset {0..N}, |E|=k, so N>=k-1. pointwise: meas(S7(E))<=meas(S7(AP_{N+1})).
    # certified if meas(S7(AP_{N+1}))<=cap_k, i.e. AP_{N+1}<=ck.
    Nstar=None
    for N in range(k-1, 60):
        if apvals.get(N+1, F(2)) <= ck:
            Nstar=N
        else:
            break
    print(f"  k={k}: cap_k={float(ck):.4f}. Pointwise lemma certifies ALL primitive E with span<=N*={Nstar} "
          f"(since meas(S7(AP_{{N*+1}}))={float(apvals.get((Nstar+1) if Nstar else 0,F(0))):.4f}<=cap_k).")

print()
print("="*92)
print("(2) Residual after pointwise lemma: span > N*(k). Is it finite mod scale?")
print("    NO -- span is unbounded (e.g. {0,1,2,...,k-2, BIG}). So pointwise alone leaves an")
print("    infinite residual. BUT that residual has small |E intersect {0..N*}|... quantify.")
print("="*92)
for k in [8]:
    ck=cap[k]
    Nstar=None
    for N in range(k-1,60):
        if apvals.get(N+1,F(2))<=ck: Nstar=N
        else: break
    print(f"  k={k}: N*={Nstar}. A residual shape e.g. E=(0,1,2,3,4,5,6, BIG):")
    for big in [Nstar+1, Nstar+5, 20, 50, 100]:
        E=(0,1,2,3,4,5,6,big)
        if gcd(big,1)!=0: pass
        v=measS7(E)
        print(f"      big={big}: span={big}, meas(S7)={float(v):.5f} <=cap_k({float(ck):.4f})? {v<=ck}")

print()
print("="*92)
print("(3) Fraction of dangerous-box shapes certified by pointwise lemma (span<=N*(k))")
print("="*92)
def gen(k,maxE):
    out=[]
    for rest in itertools.combinations(range(1,maxE+1),k-1):
        E=(0,)+rest; g=0
        for e in E: g=gcd(g,e)
        if g!=1: continue
        out.append(E)
    return out
for k in [8,9,10,11]:
    ck=cap[k]
    Nstar=None
    for N in range(k-1,60):
        if apvals.get(N+1,F(2))<=ck: Nstar=N
        else: break
    maxE = {8:16,9:17,10:16,11:16}[k]  # the dangerous boxes from canon
    shapes=gen(k,maxE)
    cert=sum(1 for E in shapes if (max(E)-min(E))<=Nstar)
    # also count how many of the UNCERTIFIED ones actually exceed cap (should be 0)
    viol=0; uncert_ok=0
    for E in shapes:
        if (max(E)-min(E))>Nstar:
            v=measS7(E)
            if v>ck: viol+=1
            else: uncert_ok+=1
    print(f"  k={k}: N*={Nstar}, box maxE<={maxE}, total shapes={len(shapes)}; "
          f"pointwise-certified (span<=N*)={cert}; uncertified={len(shapes)-cert} "
          f"(of which {uncert_ok} OK by direct compute, {viol} VIOLATE cap)")
print("\nDONE.")
