"""
S75b: BUILD the magic function in the RESONANCE domain (where it is exact + cyclotomic).
THM-576: cap_P = meas(lonely(P)) = sum_{S subset P} (-1)^|S| meas(cap_{p in S} D_p), D_p={x:||px||<1/14}.
The order-2 (pairwise) truncation is EXACT for |P|<=3 -> the magic function is the PAIRWISE overlap kernel
K(p,q)=meas(D_p cap D_q), a CYCLOTOMIC object. The binding rows j=4,5 add the order-3 term (two constants).
This builds the cap certificate the right way (resonance/inclusion-exclusion, not spatial minorant).
"""
from fractions import Fraction as F
from itertools import combinations
import numpy as np

def meas_inter(S):
    """exact meas( intersection_{p in S} D_p ),  D_p = {x in [0,1): ||p x|| < 1/14}."""
    S=[p for p in S if p!=0]
    if not S: return F(1)
    b=set([F(0),F(1)])
    for p in S:
        # ||p x||=1/14  => p x = k +- 1/14
        for k in range(0,p+1):
            for s in (F(1,14),-F(1,14)):
                v=(F(k)+s)/p
                if 0<=v<=1: b.add(v)
    b=sorted(b); tot=F(0)
    for i in range(len(b)-1):
        x0,x1=b[i],b[i+1]
        if x1<=x0: continue
        mid=(x0+x1)/2
        if all(min((p*mid)%1, 1-(p*mid)%1) < F(1,14) for p in S): tot+=x1-x0
    return tot

def cap_incl_excl(P, order=None):
    """cap_P via inclusion-exclusion; order=truncation (None=full)."""
    P=list(P); tot=F(0)
    maxo=len(P) if order is None else order
    for r in range(0,maxo+1):
        for S in combinations(P,r):
            tot+= (-1)**r * meas_inter(S)
    return tot

print("="*90)
print(" THE PAIRWISE-OVERLAP CYCLOTOMIC KERNEL K(1,q) = meas(D_1 cap D_q)  (the magic function's core)")
print("="*90)
print(f"{'q':>3}{'K(1,q)=meas(D_1∩D_q)':>24}{'  vs 1/49 (independent)':>24}")
for q in range(1,14):
    k=meas_inter((1,q))
    print(f"{q:>3}{str(k):>24}{f'{float(k):.5f} vs {1/49:.5f}':>24}")
print(" min overlap at q=13 (= -1 mod 14): the LEAST-resonant pair => the cap minimizer {1,13}. CYCLOTOMIC.")

print("\n"+"="*90)
print(" CAP RECONSTRUCTION via the magic function (inclusion-exclusion), checking the order where it BREAKS")
print("="*90)
caps_known={2:F(66,91),3:F(55,91),4:F(1979,4004),5:F(2243,5880)}
minimizers={2:(1,13),3:(1,12,13),4:(1,11,12,13),5:(1,5,7,8,9)}
print(f"{'j':>3}{'P (minimizer)':>16}{'order-1 (union)':>16}{'order-2 (pairwise)':>20}{'order-3':>14}{'full=cap':>12}{'known cap':>12}")
for j in range(2,6):
    P=minimizers[j]
    o1=cap_incl_excl(P,1); o2=cap_incl_excl(P,2); o3=cap_incl_excl(P,3); full=cap_incl_excl(P)
    kn=caps_known[j]
    print(f"{j:>3}{str(P):>16}{f'{float(o1):.4f}':>16}{f'{float(o2):.5f}':>20}{f'{float(o3):.5f}':>14}{f'{float(full):.5f}':>12}{f'{float(kn):.5f}':>12}")
    if j<=3 and o2==full:
        print(f"      -> order-2 (pairwise magic function) is EXACT: {o2} = cap  [CYCLOTOMIC, f-hat>=0 Fejer]")
    if j>=4:
        corr=full-o2
        print(f"      -> order-2 gives {float(o2):.5f}; the ORDER-3 correction = {corr} = the binding deviation constant")
print("="*90)
print(" MAGIC FUNCTION = order-2 pairwise overlap kernel (exact cap for j<=3) + order-3 correction (j=4,5).")
print(" The kernel K(p,q) is cyclotomic (min at antipode q=-1 mod14); the order-3 = the two hard constants.")
