#!/usr/bin/env python3
"""Symbolic coverage extremality: the inclusion-exclusion ANATOMY of cap_k = min meas(lonely(P)) (S64).

cap_k = min_{|P|=j} meas(lonely(P)),  j = 13-k,  lonely(P) = {x in [0,1): ||p x|| >= 1/14 for all p in P}.
THM-576: cap_k = C(14-j,2)/C(14,2) = C(k+1,2)/91, PROVED only for j=1,2; j=3 verified by search; j=4,5 dips.

This decomposes meas(lonely(minimizer_j)) by INCLUSION-EXCLUSION ORDER:
   meas(lonely(P)) = sum_{T subset P} (-1)^|T| F(T),   F(T) = meas{x: ||p x||<1/14 for all p in T}
                   = sum_{r=0}^{j} (-1)^r [ sum_{|T|=r} F(T) ]   =  sum_r ORDER_r.
ORDER_0 = 1, ORDER_1 = -j/7. The claim cap = C(14-j,2)/91 = 1 - j/7 + C(j,2)/91 means the net of
ORDERS >= 2 must equal +C(j,2)/91. We test: do ORDERS >= 3 VANISH (so cap = pairwise truncation) for
j<=3, and is the DIP at j=4,5 exactly the first nonvanishing high order?  All EXACT rationals.
"""
import sys, itertools
from fractions import Fraction as F
from math import comb
sys.stdout.reconfigure(line_buffering=True)

def norm(x):
    f = x - int(x)
    if f < 0: f += 1
    return min(f, 1 - f)

def joint_forbidden_measure(T, kk1=14):
    """F(T) = meas{x in [0,1): ||p x|| < 1/kk1 for ALL p in T}.  Exact via breakpoints."""
    T = [p for p in T if p != 0]
    if not T: return F(1)
    bp = set([F(0), F(1)])
    for p in T:
        for j in range(0, p + 1):
            for s in (F(-1), F(1)):
                x = F(j, p) + s * F(1, kk1 * p)
                if 0 <= x <= 1: bp.add(x)
    bp = sorted(bp); acc = F(0)
    for i in range(len(bp) - 1):
        a, b = bp[i], bp[i + 1]
        if b <= a: continue
        mid = (a + b) / 2
        if all(norm(p * mid) < F(1, kk1) for p in T): acc += b - a
    return acc

def lonely_orders(P, kk1=14):
    """Return [ORDER_0..ORDER_j] where ORDER_r = (-1)^r sum_{|T|=r} F(T), and the running partial sums."""
    P = list(P); j = len(P); orders = []
    for r in range(j + 1):
        s = F(0)
        for T in itertools.combinations(P, r):
            s += joint_forbidden_measure(T, kk1)
        orders.append((-1)**r * s)
    return orders

# THM-576 minimizers (the top-cluster, then the j=5 BREAK)
MINIMIZER = {
    1: [1],
    2: [1, 13],
    3: [1, 12, 13],
    4: [1, 11, 12, 13],
    5: [1, 5, 7, 8, 9],     # the BREAK (middle-spread, 3-correlated)
}
print("=" * 92)
print(" INCLUSION-EXCLUSION ANATOMY of cap_k = meas(lonely(minimizer_j)),  k = 13-j")
print("=" * 92)
for j in range(1, 6):
    P = MINIMIZER[j]; k = 13 - j
    orders = lonely_orders(P)
    total = sum(orders)
    pascal = F(comb(14 - j, 2), 91)          # the pair-Pascal target C(14-j,2)/91 = C(k+1,2)/91
    pair_trunc = sum(orders[:3])             # ORDER_0 + ORDER_1 + ORDER_2  (Bonferroni even-truncation)
    high = sum(orders[3:])                   # ORDERS >= 3 (the net higher-order correction)
    print(f"\n j={j} (k={k}), minimizer P={P}:")
    print(f"   orders r=0..{j}: " + ", ".join(f"O{r}={str(o)}" for r, o in enumerate(orders)))
    print(f"   ORDER_0+1 = 1 - j/7 = {str(sum(orders[:2]))} = {float(sum(orders[:2])):.5f}")
    print(f"   ORDER_2 (pairwise total, signed +) = {str(orders[2]) if j>=2 else '0'} "
          f"(target C(j,2)/91 = {str(F(comb(j,2),91))})")
    print(f"   pairwise truncation O0+O1+O2 = {str(pair_trunc)} = {float(pair_trunc):.6f}")
    print(f"   ORDERS>=3 (net high-order) = {str(high)} = {float(high):.6f}  "
          f"{'(VANISH)' if high==0 else '(nonzero = the DIP source)'}")
    print(f"   TRUE meas(lonely) = {str(total)} = {float(total):.6f}")
    print(f"   pair-Pascal C(14-j,2)/91 = {str(pascal)} = {float(pascal):.6f}   "
          f"dip = pascal - true = {str(pascal-total)}")

print("\n" + "=" * 92)
print(" READING")
print("=" * 92)
print(" If ORDERS>=3 vanish for j<=3, then cap = O0+O1+O2 = pairwise truncation, and the only task is the")
print(" CLOSED FORM of ORDER_2 (the total pairwise overlap) = C(j,2)/91. Then cap = 1 - j/7 + C(j,2)/91 =")
print(" C(14-j,2)/91 SYMBOLICALLY (k>=10). The DIP at j=4,5 = the first nonvanishing ORDERS>=3 (triple/quad).")
