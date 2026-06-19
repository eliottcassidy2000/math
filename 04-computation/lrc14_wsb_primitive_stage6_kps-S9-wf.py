#!/usr/bin/env python3
"""
lrc14_wsb_primitive_stage6_kps-S9-wf.py  (kind-pasteur-2026-06-19-S9)

Stage 6.  Established so far (exact):
  * Lemma B (meas(S7) same-k extremality): consec maximizes meas(S7) among same-k sets, k=8,
    maxE<=16 (11440 sets), ZERO beaters.  STRONGER than HYP-2607 (needs no L_y).
  * P(N=6)<=1/(7maxE) holds for PRIMITIVE E; the violations were all NON-primitive (dilations).
  * component-at-0 of the all-in-[0,1/7) set is EXACTLY [0,1/(7 maxE)).

NOTE on the LRC reduction: the cluster speeds in the actual reduction are PRIMITIVE (gcd=1) by
THM-531 scale-invariance (we may divide by gcd). So restricting to primitive E is legitimate.

GOALS
 (G1) PRIMITIVE-restricted Lemma B: consec maximizes meas(S7) over primitive same-k sets.
      Test k=8 (wide), k=9 (moderate).  Also test k=8 over a RANDOM wide-spread bank (maxE up to
      60) to probe the wide regime (where exhaustive is infeasible) -- the wide-spread bound regime.
 (G2) PRIMITIVE P(N=6) bound: P(N=6)(E) <= 1/(7(k-1)) for primitive E, equality iff consec.
      Exhaustive primitive k=4..8; random primitive wide.
 (G3) component structure of meas(S7) near 0: is there an analogue of [0,1/(7maxE))?  i.e. does
      the meas(S7) component at 0 also have a clean form?  (S7 = all 7 sectors hit; near 0 the
      orbit {ex} is a tiny cluster near 0, so only sector 0 is hit -- S7 does NOT contain a nbhd
      of 0.  So meas(S7) has NO 0-component; its structure is different.  Document.)
 (G4) THE WIDE-SPREAD SIDE: confirm large-spread meas(S7) collapses (HYP-2608 (a)).  For primitive
      random E with growing spread, meas(S7) -> small.  Tabulate vs cap_8.  Identify the RESONANT
      wide configs (the canon warns about w==0 mod 7): test arithmetic-progression-like wide sets
      and 7-divisible structured sets specifically.
"""
import sys, itertools, random
from math import gcd
from functools import reduce
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

def gcd_all(E):
    return reduce(gcd, [e for e in E if e != 0], 0)

def measS7(E):
    """meas{ all 7 sectors {0..6} hit by orbit {frac(e x): e in E} }."""
    E=sorted(set(E)); bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for j in range(7):
            c=F(j,7)
            for m in range(e): bps.add((c+m)/e)
    bps=sorted(z for z in bps if F(0)<=z<F(1)); tot=F(0)
    for i in range(len(bps)):
        x0=bps[i]; x1=bps[i+1] if i+1<len(bps) else F(1)
        if x1<=x0: continue
        xm=(x0+x1)/2; hit=set()
        for e in E:
            fr=(e*xm)%1; hit.add((fr.numerator*7)//fr.denominator)
        if len(hit)==7: tot+=x1-x0
    return tot

def measPN6(E):
    """P(N=6)=meas{all frac(ex) in [0,1/7)}."""
    E=sorted(set(e for e in E if e!=0)); a,b=F(0),F(1,7); bps=set([F(0),F(1)])
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

print("="*78)
print("(G1) PRIMITIVE Lemma B: consec maximizes meas(S7) among primitive same-k sets")
print("="*78)
for k,B in [(8,14),(9,12)]:
    pc=measS7(list(range(k)))
    over=0; cnt=0; worst=(F(0),None)
    for r in itertools.combinations(range(1,B+1),k-1):
        E=[0]+list(r)
        if gcd_all(E)!=1: continue
        cnt+=1
        m=measS7(E)
        if m>pc+F(1,10**12):
            over+=1
            if m>worst[0]: worst=(m,E)
    msg=f"  k={k} primitive maxE<={B}: {cnt} sets, meas(S7)(consec)={float(pc):.6f}, beaters={over}"
    if worst[1]: msg+=f" WORST {float(worst[0]):.6f} @ {worst[1]}"
    print(msg)

print()
print("  k=8 RANDOM primitive WIDE-SPREAD bank (maxE up to 60) -- the wide regime:")
random.seed(7); pc8=measS7(list(range(8)))
over=0; n=0; mx=(F(0),None); buckets={}
for _ in range(3000):
    M=random.randint(20,60)
    rest=sorted(random.sample(range(1,M+1),7))
    E=[0]+rest
    if gcd_all(E)!=1: continue
    if max(E)!=M: continue
    n+=1; m=measS7(E)
    sp=max(E)
    if m>pc8: over+=1
    if m>mx[0]: mx=(m,E)
    bk = sp//10*10
    buckets.setdefault(bk,[]).append(float(m))
print(f"    tested {n} random primitive wide E; meas(S7)>consec(0.3272): {over}; "
      f"max meas(S7)={float(mx[0]):.5f} @ {mx[1]}")
for bk in sorted(buckets):
    vals=buckets[bk]; print(f"    spread~[{bk},{bk+10}): n={len(vals)} mean meas(S7)={sum(vals)/len(vals):.4f} max={max(vals):.4f}")

print()
print("="*78)
print("(G2) PRIMITIVE P(N=6) bound: <= 1/(7(k-1)), equality iff consec")
print("="*78)
for k,B in [(4,12),(5,12),(6,12),(7,12),(8,13)]:
    bound=F(1,7*(k-1)); over=0; cnt=0; eq=[]
    for r in itertools.combinations(range(1,B+1),k-1):
        E=[0]+list(r)
        if gcd_all(E)!=1: continue
        cnt+=1; m=measPN6(E)
        if m>bound+F(1,10**12): over+=1
        if m==bound: eq.append(E)
    print(f"  k={k} primitive maxE<={B}: {cnt} sets; P(N=6)>1/(7(k-1))={bound}: {over}; "
          f"#equality={len(eq)} (consec in eq: {list(range(k)) in eq})")

print()
print("="*78)
print("(G4) wide-spread RESONANT configs (w==0 mod 7) -- do they break the collapse?")
print("="*78)
pc8=measS7(list(range(8)))
print(f"  consec_8 meas(S7)={float(pc8):.5f}, cap_8~0.38146")
# structured wide sets: APs {0,d,2d,...}, 7-divisible, etc.
tests = {
 "AP step 1 (consec)": list(range(8)),
 "AP step 2 (nonprim)": [0,2,4,6,8,10,12,14],
 "AP step 3 (nonprim)": [0,3,6,9,12,15,18,21],
 "wide 7-mult {0,7,14..}": [0,7,14,21,28,35,42,49],
 "0..6 + 7": [0,1,2,3,4,5,6,7],
 "0..6 + 14": [0,1,2,3,4,5,6,14],
 "0..6 + 21": [0,1,2,3,4,5,6,21],
 "geometric 2^i": [0,1,2,4,8,16,32,64],
 "primes": [0,2,3,5,7,11,13,17],
 "wide random spread~50": [0,3,11,19,27,34,41,49],
}
for name,E in tests.items():
    g=gcd_all(E); m=measS7(E); flag="OVER cap!" if m>F(2243,5880) else ""
    print(f"  {name:30s} gcd={g} maxE={max(E)}: meas(S7)={float(m):.5f} {flag}")
print("\nDONE stage 6.")
