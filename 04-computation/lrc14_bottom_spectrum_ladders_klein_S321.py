#!/usr/bin/env python3
"""lrc14_bottom_spectrum_ladders_klein_S321.py -- klein-2026-07-19-S321.

THE BOTTOM SPECTRUM IS TWO LADDERS ROOTED AT THE AP, WITH PLATEAUS.

Exact-M verifications (referee engine = lrc14_subgap_referee_klein_S319.exact_M,
complete-breakpoint method, gate-validated on six known families) for:

(1) THE 12m-LADDER on two different 12-bases:
    L1(m) = {1..11,13} u {12m}      (THM-1230's ladder; L1(1) = {1..13} = the AP)
    L2(m) = {1..9,11,13,20} u {12m} (the deformed base found by the S321 census)
    Law: M(base u {12m}) = max(m/(12m+5), c_base), plateaus where c_base wins:
      L1: 1/14, 1/14, 3/41, 4/53, 1/13(=5/65), 6/77, ...  (c_{L1} = 1/14)
      L2: 2/27, 2/27, 2/27, 4/53, 1/13,        6/77, ...  (c_{L2} = 2/27)
    Mechanism: the far element 12m kills the base's q=12 maximizer instantly
    (12m * a/12 == 0 mod 1); what remains is the two-peak competition between
    the base-internal pinch (c_base; L2's is the pair (7,20) at q=27) and the
    universal far pinch on the active pair (5, 12m) worth m/(12m+5).
    The ladder accumulates at 1/12 = M(base) exactly (both bases).

(2) THE 13s-LADDER (Kravitz rungs): K(s) = {1..12} u {13s}, K(1) = the AP.
    M = s/(13s+1) exactly for s = 2..6 (2/27, 3/40, 4/53, 5/66, 6/79).

(3) Consequences recorded:
    - h_min(2/27) = 20 (via L2-type base at m=1); h_min(3/41) <= 36 < h_min(3/40)
      in (37..39]: the k=2 mediant realizes BEFORE the k=1 Kravitz rung despite
      being the smaller value -- height order inverts value order across strata.
    - Realizer multiplicity = plateau length (+ cross-ladder coincidences):
      2/27 has >= 4 realizers by h=36: K(2) and the L2-plateau m=1,2,3.
    - Both ladders share the root K(1) = L1(1) = {1..13}; 4/53 and 1/13 are
      cross-ladder coincidences (L1(4)=L2(4)=4/53 value-wise with K(4)).

Run: python3 lrc14_bottom_spectrum_ladders_klein_S321.py
     | tee ../05-knowledge/results/lrc14_bottom_spectrum_ladders_klein_S321.out
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction
from lrc14_subgap_referee_klein_S319 import exact_M

def show(name, V, predict=None):
    m = exact_M(sorted(V))
    tag = ""
    if predict is not None:
        tag = f"  predict {predict}  [{'OK' if m == predict else 'MISMATCH'}]"
    print(f"{name:34s} M = {str(m):8s}{tag}")
    return m

print("== (1a) L1(m) = {1..11,13} u {12m} ==")
L1base = list(range(1, 12)) + [13]
for m in range(1, 8):
    pred = max(Fraction(m, 12*m+5), Fraction(1, 14))
    if Fraction(m, 12*m+5) < Fraction(1, 14):
        pred = Fraction(1, 14)
    show(f"L1({m}) = base u {{{12*m}}}", L1base + [12*m] if 12*m not in L1base else L1base + [12*m], pred)

print("== (1b) L2(m) = {1..9,11,13,20} u {12m} ==")
L2base = list(range(1, 10)) + [11, 13, 20]
for m in range(1, 11):
    pred = max(Fraction(m, 12*m+5), Fraction(2, 27))
    show(f"L2({m}) = base u {{{12*m}}}", L2base + [12*m], pred)
show("L2(18)", L2base + [216], Fraction(18, 221))

print("== (1c) base caps (12-sets) ==")
show("M({1..9,11,13,20}) [12-set]", L2base, Fraction(1, 12))

print("== (2) K(s) = {1..12} u {13s} (Kravitz rungs) ==")
K = list(range(1, 13))
for s in range(2, 7):
    show(f"K({s}) = {{1..12,{13*s}}}", K + [13*s], Fraction(s, 13*s+1))

print("== (3) census cross-checks (complete sub-1/13 table, h <= 36) ==")
for V, ex in [
    ([1,2,3,4,5,6,7,8,9,10,11,12,13], Fraction(1,14)),
    ([1,2,3,4,5,6,7,8,9,10,11,13,24], Fraction(1,14)),
    ([1,2,3,4,5,6,7,8,9,11,12,13,20], Fraction(2,27)),
    ([1,2,3,4,5,6,7,8,9,11,13,20,24], Fraction(2,27)),
    ([1,2,3,4,5,6,7,8,9,10,11,12,26], Fraction(2,27)),
    ([1,2,3,4,5,6,7,8,9,11,13,20,36], Fraction(2,27)),
    ([1,2,3,4,5,6,7,8,9,10,11,13,36], Fraction(3,41)),
    ([1,2,3,4,5,6,7,8,9,10,11,12,39], Fraction(3,40)),
    ([1,2,3,4,5,6,7,8,9,10,11,12,52], Fraction(4,53)),
]:
    show(str(V), V, ex)
print("done.")
