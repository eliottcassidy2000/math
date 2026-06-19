#!/usr/bin/env python3
"""
lrc14_angleC_threegap_mechanism_macmini_0618s7.py  (mac-mini-2026-06-18-S7)

REFINED three-gap mechanism for AP-extremality of meas(S7).

From the run [C] we learned: meas(S7) is NOT monotone in additive ENERGY. So the naive
"energy => cover" deployment FAILS. This script isolates the CORRECT three-gap invariant.

The right object (Angle C, literature-grounded):
  meas(S7(E)) = integral over x of 1{orbit P(x) hits all 7 sectors}.
  The orbit P(x) = {frac(e_i x)} is a union of dilated arithmetic progressions whose gap
  structure obeys the (Liang) 3d-distance theorem. The cover event is best understood
  SECTOR by SECTOR with INCLUSION-EXCLUSION (THM-532 form):
     meas(S7) = M7(k) + corr(E),  M7(k) = main (independence) term,
  and corr(E) is a sum over the OFFSET RELATION LATTICE. The relation lattice of the AP is
  the RICHEST possible (a_i - a_j realize EVERY value 1..k-1 with max multiplicity), so corr
  is largest for the AP. THIS is the precise mechanism, and it matches additive structure but
  via the *difference multiset multiplicity*, not the additive energy of sums.

This script:
 (1) Tests whether meas(S7(E)) is monotone in the DIFFERENCE-MULTISET richness:
     r(E) = sum over d of (mult of d as a difference)^2  =  the "difference energy"
     (this is the count of additive QUADRUPLES a-b=c-d, i.e. the genuine additive energy of
     the DIFFERENCE set, which for the AP is maximal). Compare to meas(S7) across many shapes.
 (2) Confirms the AP UNIQUELY maximizes the difference energy among primitive k-sets in a box,
     and checks correlation rank between meas(S7) and difference-energy.
 (3) Tests a CLEAN three-gap SUFFICIENT certificate for the dangerous rows:
     "if E has difference energy < threshold_k then meas(S7(E)) < cap_k" — i.e. whether a
     single scalar (difference energy) certificate separates all shapes from cap.
 (4) Reports the worst (largest meas(S7)) NON-AP shape per row and its difference energy,
     to see how close a non-AP can get and whether the certificate would be tight.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from collections import Counter
sys.stdout.reconfigure(line_buffering=True)
SEV=F(1,7)

def measS7(E):
    E=sorted(set(E))
    bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps)
    total=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        secs=set()
        for e in E:
            secs.add(int(((e*xm)%1)*7))
        if len(secs)==7: total+=x1-x0
    return total

def diff_energy(E):
    """additive energy of the difference structure: sum_d mult(d)^2 where d ranges over all
    ordered differences a-b (a,b in E). For AP {0..k-1} this is maximal among k-sets."""
    c=Counter()
    El=list(E)
    for a in El:
        for b in El:
            c[a-b]+=1
    return sum(v*v for v in c.values())

def normalize_primitive(E):
    E=sorted(set(E)); E=[e-E[0] for e in E]
    g=0
    for e in E: g=gcd(g,e)
    if g>1: E=[e//g for e in E]
    return tuple(E)

# exact caps from canon
CAP={8:F(2243,5880),9:F(1979,4004),10:F(55,91),11:F(66,91),12:F(6,7),13:F(1)}

print("="*92)
print("ANGLE C refined: difference-energy as the three-gap invariant behind AP-extremality")
print("="*92)

# (1)+(2)+(3)+(4): per dangerous row, full small-box census
def census(k,B):
    AP=tuple(range(k))
    ap_s7=measS7(list(AP)); ap_de=diff_energy(AP)
    rows=[]
    seen=set()
    for combo in itertools.combinations(range(1,B+1),k-1):
        En=normalize_primitive((0,)+combo)
        if En in seen: continue
        seen.add(En)
        rows.append((measS7(list(En)), diff_energy(En), En))
    rows.sort(reverse=True)
    return AP,ap_s7,ap_de,rows

for k,B in [(8,12),(9,12),(10,12),(11,12)]:
    AP,ap_s7,ap_de,rows=census(k,B)
    cap=CAP[k]
    # is AP the unique max of meas(S7)? of diff_energy?
    max_s7=rows[0][0]; ap_is_s7max = (rows[0][2]==AP and (len(rows)==1 or rows[1][0]<ap_s7))
    max_de=max(r[1] for r in rows); ap_is_demax=(ap_de==max_de and sum(1 for r in rows if r[1]==max_de)==1)
    # rank correlation between meas(S7) and diff_energy (Spearman, quick)
    import statistics
    s7v=[float(r[0]) for r in rows]; dev=[float(r[1]) for r in rows]
    # pearson on raw
    n=len(rows); ms=sum(s7v)/n; md=sum(dev)/n
    cov=sum((a-ms)*(b-md) for a,b in zip(s7v,dev))
    sds=(sum((a-ms)**2 for a in s7v))**.5; sdd=(sum((b-md)**2 for b in dev))**.5
    pear=cov/(sds*sdd) if sds*sdd>0 else 0
    # best NON-AP shape
    nonap=[r for r in rows if r[2]!=AP]
    best_nonap=nonap[0] if nonap else None
    print(f"\n--- k={k} (box maxE<={B}, {len(rows)} primitive shapes), cap={float(cap):.5f} ---")
    print(f"  meas(S7(AP))      = {float(ap_s7):.6f}   diff_energy(AP) = {ap_de}")
    print(f"  AP is UNIQUE meas(S7) max? {ap_is_s7max}   AP is UNIQUE diff_energy max? {ap_is_demax}")
    print(f"  Pearson corr(meas(S7), diff_energy) over all shapes = {pear:.4f}")
    if best_nonap:
        print(f"  best NON-AP: meas(S7)={float(best_nonap[0]):.6f}  diff_energy={best_nonap[1]}  E={best_nonap[2]}")
    # (3) single-scalar certificate: is there a diff-energy threshold separating cap?
    #     we want: meas(S7) <= cap. Check if meas(S7) is monotone enough that
    #     "all shapes have meas(S7) <= meas(S7(AP)) < cap" already (AP-extremality) -- the
    #     certificate is then just AP-extremality. Report the max meas(S7) vs cap.
    print(f"  max meas(S7) over box = {float(max_s7):.6f}  < cap? {max_s7<cap}  "
          f"(slack {float(cap-max_s7):.6f})")

print("\n[CONCLUSION CHECK] Does diff_energy PERFECTLY rank meas(S7) (so a clean monotone")
print("three-gap certificate 'meas(S7) increasing in diff_energy' would hold)?")
print("If Pearson<1 and AP is still both maxima, the monotone is APPROXIMATE: AP-extremality")
print("holds but is NOT reducible to a single-scalar diff-energy inequality.")
print("\nDONE.")
