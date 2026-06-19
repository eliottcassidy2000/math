#!/usr/bin/env python3
"""
lrc14_B7_iid_floor_kps-S6-wf.py
The asymptotic (large-spread / iid) floor of the 7-arc NET MINORANT B_7(E).

B_7(E)=meas{x: at least one of the 7 fixed width-1/7 arcs A_i=[(2i+1)/14,(2i+3)/14)
            is empty of {frac(e x): e in E}}  is a RIGOROUS lower bound on mu_{1/7}(E).
The 7 arcs A_0..A_6 TILE the circle (a partition into 7 equal cells of measure 1/7).

For an EQUIDISTRIBUTED (iid-uniform) configuration of the 8 points, the probability
that all 7 cells are occupied is the classical surjective-occupancy count:
   P_all = sum_{j=0}^{7} (-1)^j C(7,j) (1 - j/7)^8,
so the iid floor of the minorant is  B_7^iid = 1 - P_all.

This is the value B_7(E) approaches as spread(E) -> infinity (equidistribution of
{e x} for the differences).  We compute it EXACTLY and compare to thr_8.

KEY POINT (why the B_7 route looks closeable):  both the bounded-spread MINIMUM of
B_7 (~0.94, exhaustive) and the large-spread iid FLOOR (~0.9755) lie FAR ABOVE
thr_8=0.6185 (margin >= 0.32).  B_7 never approaches thr_8.  A rigorous uniform
B_7(E) >= thr_8 (=> mu_{1/7}(E) >= thr_8 for ALL E, closing binding k=8) therefore
needs only a CRUDE equidistribution/discrepancy estimate (e.g. Erdos-Turan: once the
discrepancy of {e x : e in E-differences} is below an explicit constant, B_7 stays
within the iid value +- small), plus the finite bounded-spread exhaustion.  The
enormous margin is what distinguishes this from the net-upper-bound route (whose
relevant quantity is tiny and offers no such slack for a crude bound).
"""
import sys
from fractions import Fraction as F
from math import comb
sys.stdout.reconfigure(line_buffering=True)

if __name__ == "__main__":
    thr8 = F(3637, 5880)
    # iid floor for k=8 points, 7 cells
    P_all_8 = sum((-1)**j * comb(7, j) * F(7-j, 7)**8 for j in range(8))
    B7_iid_8 = 1 - P_all_8
    print("B_7 iid floor (large-spread limit of the 7-arc net minorant):")
    print(f"  P(all 7 cells occupied by 8 iid pts) = {P_all_8} = {float(P_all_8):.6f}")
    print(f"  B_7^iid(k=8) = {B7_iid_8} = {float(B7_iid_8):.6f}")
    print(f"  thr_8        = {thr8} = {float(thr8):.6f}")
    print(f"  B_7^iid - thr_8 = {float(B7_iid_8 - thr8):.6f}   (>=0: {B7_iid_8 >= thr8})")
    print()
    # general k: m=14*Vmax cluster of size k>=8; here the cluster has |E|=k points,
    # 7 cells still (arcs of width 1/7).  P(all 7 cells occupied by k iid pts):
    print("General k: iid floor B_7^iid(k) vs thr_k (7 cells, k points):")
    THRA = {8:F(3637,5880),9:F(2025,4004),10:F(36,91),11:F(25,91),12:F(1,7),13:F(0)}
    for k in range(8, 14):
        P_all = sum((-1)**j * comb(7, j) * F(7-j, 7)**k for j in range(8))
        B7 = 1 - P_all
        thr = THRA[k]
        print(f"  k={k:2d}: B_7^iid={float(B7):.6f} ({B7})  thr_k={float(thr):.6f}  "
              f"margin={float(B7-thr):+.6f}  {'OK' if B7>=thr else '*** below thr ***'}")
    print("\nNOTE: B_7 is a lower bound on mu; iid floor>=thr_k at every k with margin>=0.30.")
    print("The bounded-spread exhaustive minimum (separate script) is also >> thr_k.")
