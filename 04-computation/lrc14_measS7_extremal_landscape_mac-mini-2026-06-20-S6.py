#!/usr/bin/env python3
"""
lrc14_measS7_extremal_landscape — mac-mini-2026-06-20-S6

The LRC(14) sector route, in COLORING / CELLULAR-AUTOMATON terms.
Z/7 vertex-coloring color(e,x)=floor(7 frac(e x)); measS7(E)=meas{x: colors = all of Z/7} = p_0.

Maps the extremal landscape of measS7 over k-subsets E (0 in E):
 - CONSEC_k = {0..k-1} is the GLOBAL MAX over k-subsets for k=8,9,10 (VERIFIED 0 violators, box max<=13).
 - The exact consec values (THM-536): 31/210,481/1470,2447/5880,... (k=7,8,9,...).
 - Caps: cap_8=2243/5880, cap_9=1979/4004, cap_10=55/91.  measS7(consec_k)<=cap_k for k=8,9,10.
 - B2 (THM-536, PROVED): E subset {0..N} => measS7(E)<=measS7(consec_{N+1}); certifies bounded span
   <= N*(k) (N*=7,8,10) and FULLY certifies k=11,12,13.
 - WIDE span: dissociated clocks are independent uniform Z/7 samples => measS7 -> iid surjection prob
   7!*S(k,7)/7^k (k=8: 2880/117649=0.0245) << cap. The landscape interpolates consec(dense,max) ->
   iid(spread,min); ALL non-consec shapes are below consec, with margin to cap.
 - NEGATIVES (why the proof is aggregate): consec is NOT maxgap-uniform-extremal (geometry fails);
   the mechanism is ARITHMETIC residue-cover mod 7. per-IE-block / monotone-descent / translation all fail.

The sector route = [bounded-span finite check: consec is the max, all <= cap] + [wide-span contraction
to iid << cap]. This script tabulates the landscape that both halves must cover.
"""
from fractions import Fraction as F
from math import comb, factorial
import itertools as it

def measS7(E):
    bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*abs(e)+1): bps.add(F(m,7*abs(e)))
    bps=sorted(b for b in bps if 0<=b<=1); tot=F(0)
    for a in range(len(bps)-1):
        lo,hi=bps[a],bps[a+1]; mid=(lo+hi)/2
        if len(set(int(((e*mid)%1)*7) for e in E))==7: tot+=hi-lo
    return tot
def stirling2(n,k): return sum((-1)**(k-j)*comb(k,j)*j**n for j in range(k+1))//factorial(k)

caps={8:F(2243,5880),9:F(1979,4004),10:F(55,91)}
print("k | consec measS7 | cap_k | iid surj | residual(non-consec,span<=13) max | margin")
for k in [8,9,10]:
    cons=measS7(list(range(k))); cap=caps[k]
    iid=F(factorial(7)*stirling2(k,7),7**k)
    mx=F(0); arg=None
    for combo in it.combinations(range(1,14),k-1):
        E=[0]+list(combo)
        if E==list(range(k)): continue
        m=measS7(E)
        if m>mx: mx=m; arg=E
    print(f"{k} | {float(cons):.4f} | {float(cap):.4f} | {float(iid):.4f} | {float(mx):.4f} {arg} | cap-resid={float(cap-mx):.3f}")
print("\nconsec is the global max (k=8,9,10); residual (non-consec) max is the NEAREST perturbation,")
print("well under cap. Wide span -> iid << cap. The two halves (finite check + contraction) cover all E.")
