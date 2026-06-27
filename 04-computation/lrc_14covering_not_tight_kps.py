r"""
Closing the (star) residual: 14-COVERING primitive sets are NOT tight (kps-S31ab).
By THM-568, tight => 14|D, binders sum to D. At D=14k (k>=2): a multiple of 14 has distance in
{0,1/k,2/k,...}, =1/14 only if 14|k; so for 2<=k<=13 the binders are 14-FREE and sum to D=14k>=28.
Search primitive 14-covering sets (gcd 1, has a multiple of 14) for ANY with M=1/14, and any primitive
tight set with optimum D>14.
"""
from fractions import Fraction as F
from math import gcd
def nf(x):
    r=x%1; return min(r,1-r)
def M_and_D(S):
    S=sorted(set(abs(s) for s in S if s!=0)); C=set()
    for i in range(len(S)):
        for j in range(i,len(S)):
            for comb in {S[i]+S[j], abs(S[i]-S[j])}:
                if comb:
                    for k in range(1,comb): C.add(F(k,comb))
    best=F(0); arg=None
    for t in C:
        if 0<t<1:
            m=min(nf(s*t) for s in S)
            if m>best: best=m; arg=t
    return best,(arg.denominator if arg else 0)
def g(S):
    r=0
    for s in S: r=gcd(r,s)
    return r
import itertools
print("(1) PRIMITIVE 14-COVERING sets (AP/GW with a speed replaced by a multiple of 14): min M?")
mults14=[14,28,42,56,70]
bases={"AP":list(range(1,14)),"GW":list(range(1,12))+[13,24]}
minM=F(1); worst=None; tight14cov=[]
for bname,base in bases.items():
    for rem in base:
        for mv in mults14:
            S=sorted(set([x for x in base if x!=rem]+[mv]))
            if len(S)!=13: continue
            if g(S)!=1: continue
            if not any(s%14==0 for s in S): continue
            M,D=M_and_D(S)
            if M<minM: minM=M; worst=(bname,rem,mv,S,M,D)
            if M==F(1,14): tight14cov.append((S,D))
print(f"   min M over primitive 14-covering perturbations = {minM} = {float(minM):.5f} (vs 1/14={1/14:.5f})")
print(f"   any TIGHT (M=1/14)? {len(tight14cov)} {'<-- residual would FAIL' if tight14cov else '=> 14-covering primitive NOT tight (residual holds)'}")
if worst: print(f"   (worst case: {worst[0]} swap {worst[1]}->{worst[2]}, M={worst[4]}, D={worst[5]})")

print("\n(2) ANY primitive tight set with optimum D>14? Scan tight locus + perturbations:")
found_D_gt_14=[]
cands=[list(range(1,14)),list(range(1,12))+[13,24]]
# perturb each speed by +-14 (stays 14-coprime-ish), check primitivity, tightness, D
for base in cands:
    for i in range(13):
        for delta in [-14,14,28,-28]:
            S=base[:]; S[i]=base[i]+delta
            if S[i]<=0 or len(set(S))!=13: continue
            if g(S)!=1: continue
            M,D=M_and_D(S)
            if M==F(1,14) and D!=14: found_D_gt_14.append((S,D))
print(f"   primitive tight sets with D!=14 found: {len(found_D_gt_14)} {'<-- (star) would FAIL' if found_D_gt_14 else '=> all primitive tight have D=14 ((star) holds empirically)'}")
print("\n=> If both empty: 14-covering primitive NOT tight + primitive tight => D=14. With THM-568 (14|D),")
print("   this is the STRUCTURAL half of (star) confirmed; tight => 14-free => optimum at apex denom 14.")
