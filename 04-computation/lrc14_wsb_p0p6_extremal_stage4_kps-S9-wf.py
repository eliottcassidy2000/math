#!/usr/bin/env python3
"""
lrc14_wsb_p0p6_extremal_stage4_kps-S9-wf.py  (kind-pasteur-2026-06-19-S9)

Stage 4.  Stage 3 BREAKTHROUGH (in the PROBABILITY basis, not the moment basis):
  * consec UNIQUELY maximizes p_0=meas(S7), p_6, AND p_0+p_6 over same-k sets (0 beaters in bank).
  * The N=6 set for consec_k is the SINGLE interval [0, 1/(7(k-1))); p_6 = 1/(7(k-1)).
  * This is ORTHOGONAL to the canon's 'non-separability' which was in the MOMENT basis S_r.

CAUTION/SELF-SKEPTICISM: if consec maximizes p_0=meas(S7) DIRECTLY over same-k sets, that is
STRONGER than HYP-2607 (we wouldn't even need L_y!).  THM-536 already shows meas(S7) extremality
is the open crux.  So I must check: is 'consec maximizes p_0 among same-k' actually TRUE over a
WIDE box, or does it fail for large span (where THM-536's subset-dom is the only handle)?  If it
fails, the L_y combination is what rescues it.  Pin this down precisely.

GOALS
 (G1) p_0 = meas(S7) extremality of consec among same-k: test over WIDE boxes (maxE up to 20).
      Does consec stay the max?  Find first counterexample if any.  (This is the real test of
      whether meas(S7) itself, not just L_y, is consec-extremal.)
 (G2) PROVE p_6(consec_k) = 1/(7(k-1)) AND it is the GLOBAL max of P(N=6) over ALL k-sets
      (= the all-in-[0,1/7)-arc event).  Rigorous argument + exhaustive wide check.
 (G3) The all-in-arc event P(N=6)(E)=meas{x: frac(e x) in [0,1/7) for all e in E}.  This is a
      'simultaneous-small-fractional-part' set.  Show meas <= 1/(7*max(E))?? or 1/(7*(k-1))?
      Which is the right universal bound?  (max(E) vs k-1)
 (G4) Decompose L_y(consec)-L_y(E) into [p0(c)-p0(E)] + [p6(c)-p6(E)] + (1/10)[p3(c)-p3(E)].
      Since first two are >=0 (G1/G2, if confirmed), the ONLY possibly-adverse term is the p_3 one.
      Bound (1/10)(p3(E)-p3(c)) <= [p0(c)-p0(E)]+[p6(c)-p6(E)]?  Quantify the slack.
"""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

def law_N(E):
    E = sorted(set(E)); bps = set([F(0), F(1)])
    for e in E:
        if e == 0: continue
        for j in range(7):
            c = F(j, 7)
            for m in range(e): bps.add((c + m) / e)
    bps = sorted(z for z in bps if F(0) <= z < F(1))
    p = [F(0)] * 7
    for i in range(len(bps)):
        x0 = bps[i]; x1 = bps[i + 1] if i + 1 < len(bps) else F(1)
        if x1 <= x0: continue
        xm = (x0 + x1) / 2
        hit = set()
        for e in E:
            fr = (e * xm) % 1
            hit.add((fr.numerator * 7) // fr.denominator)
        nempty = len([j for j in range(1, 7) if j not in hit])
        p[nempty] += x1 - x0
    return p

print("="*78)
print("(G1) is p_0=meas(S7) itself consec-maximal among SAME-k sets?  WIDE box test")
print("="*78)
for k in (8,9):
    pc = law_N(list(range(k)))
    p0c = pc[0]
    print(f"--- k={k}  p0(consec)={float(p0c):.6f} ---")
    for B in (12, 15, 18, 20):
        over=0; worst=(F(0),None)
        cnt=0
        for r in itertools.combinations(range(1,B+1),k-1):
            E=[0]+list(r); cnt+=1
            p0=law_N(E)[0]
            if p0>p0c+F(1,10**12):
                over+=1
                if p0>worst[0]: worst=(p0,E)
        msg = f"shapes maxE<={B}: {cnt}, p0>consec: {over}"
        if worst[1]: msg += f"  WORST p0={float(worst[0]):.6f} E={worst[1]}"
        print(f"    {msg}")

print()
print("="*78)
print("(G2/G3) p_6(consec_k)=1/(7(k-1)); is it the GLOBAL max of the all-in-[0,1/7) event?")
print("="*78)
def measAllInArc(E,a,b):
    E=sorted(set(e for e in E if e!=0)); bps=set([F(0),F(1)])
    for e in E:
        for t in (a,b):
            for m in range(e): bps.add((t+m)/e)
    bps=sorted(z for z in bps if F(0)<=z<F(1)); tot=F(0)
    for i in range(len(bps)):
        x0=bps[i]; x1=bps[i+1] if i+1<len(bps) else F(1)
        if x1<=x0: continue
        xm=(x0+x1)/2
        if all(a<=(e*xm)%1<b for e in E): tot+=x1-x0
    return tot
# test: max P(N=6) over k-sets with maxE<=B; is it 1/(7(k-1)) (consec)? and never exceeded?
for k in (7,8,9):
    cc = measAllInArc(list(range(k)),F(0),F(1,7))
    glob=(F(0),None)
    for B in (k-1, 12, 16):
        if B<k-1: continue
        best=(F(0),None)
        for r in itertools.combinations(range(1,B+1),k-1):
            E=[0]+list(r)
            m=measAllInArc(E,F(0),F(1,7))
            if m>best[0]: best=(m,E)
        if best[0]>glob[0]: glob=best
        print(f"  k={k} maxE<={B}: max P(N=6)={best[0]}={float(best[0]):.6f}  consec=1/(7(k-1))={cc}"
              f"  consec-is-max={best[0]<=cc}")
    # hypothesis: P(N=6)(E) <= 1/(7*max(E))  for any E (the binding largest speed)
    print(f"    note consec maxE={k-1}, so 1/(7*maxE)=1/(7(k-1))=consec exactly")

print()
print("="*78)
print("(G3b) universal bound test: P(N=6)(E) <= 1/(7*max(E))?  for arbitrary E")
print("="*78)
import random
random.seed(1)
viol=0; tested=0
for trial in range(400):
    kk=random.randint(4,9)
    E=[0]+sorted(random.sample(range(1,20),kk-1))
    m=measAllInArc(E,F(0),F(1,7))
    bound=F(1,7*max(E))
    tested+=1
    if m>bound+F(1,10**12):
        viol+=1
        if viol<=5: print(f"  VIOLATION E={E} P(N=6)={float(m):.6f} > 1/(7maxE)={float(bound):.6f}")
print(f"  tested {tested} random E: P(N=6) > 1/(7 maxE) violations: {viol}")
print("  (if 0: P(N=6)<=1/(7 maxE), and consec uniquely achieves equality having min possible maxE=k-1)")

print()
print("="*78)
print("(G4) decompose L_y(consec)-L_y(E) and check the p3-term is dominated  (k=8)")
print("="*78)
k=8; pc=law_N(list(range(k)))
bank=[[0]+list(r) for r in itertools.combinations(range(1,13),k-1)]
worst_p3=(F(-10),None); n_p3adverse=0
for E in bank:
    p=law_N(E)
    d0=pc[0]-p[0]; d6=pc[6]-p[6]; d3=pc[3]-p[3]
    # L_y diff = d0 + (1/10)d3 + d6 ; want >=0.  adverse if d3<0 (E has more p3).
    if d3 < 0:
        n_p3adverse+=1
        # check d0+d6 covers -(1/10)d3
        slack = d0+d6 + F(1,10)*d3   # = L_y diff
        if slack < worst_p3[0] or worst_p3[1] is None:
            pass
        ratio = (d0+d6) / (-F(1,10)*d3) if d3<0 else None
        if worst_p3[1] is None or (ratio is not None and ratio < worst_p3[0]):
            worst_p3=(ratio,E,d0,d3,d6,slack)
print(f"  consec p0={float(pc[0]):.5f} p3={float(pc[3]):.5f} p6={float(pc[6]):.5f}")
print(f"  shapes where p3(E)>p3(consec) (adverse for the 1/10 term): {n_p3adverse}/{len(bank)}")
if worst_p3[1]:
    ratio,E,d0,d3,d6,slack=worst_p3
    print(f"  TIGHTEST adverse: E={E}")
    print(f"    d0=p0c-p0E={float(d0):.5f}  d3=p3c-p3E={float(d3):.5f}  d6=p6c-p6E={float(d6):.5f}")
    print(f"    (d0+d6)/(-d3/10) = {float(ratio):.4f}  (>1 means p0,p6 surplus covers the p3 deficit)")
    print(f"    L_y diff (must be >=0) = {float(slack):.6f}")
print("\nDONE stage 4.")
