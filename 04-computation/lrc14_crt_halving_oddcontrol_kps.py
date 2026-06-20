#!/usr/bin/env python3
"""
lrc14_crt_halving_oddcontrol_kps.py  (kind-pasteur 2026-06-19)

THE n=14=2*7 CRT / HALVING ANGLE.  Exact rational arithmetic.

Question 1 (CRT decoupling of the predicate):
  For the lonely predicate ||v tau|| > 1/14, does CRT on the gap 1/14 = (1/2)(1/7)
  give independent control via the residue-mod-7 sector and a 2-adic/dyadic part?
  We test: does the EVEN sub-core {even speeds} add ANY constraint beyond the ODD
  sub-core, given Glaisher even = 2^a * odd?

Question 2 (odd controls even):
  Each even speed e=2^a*b (b odd) has danger band at k/e of half-width 1/(14 e).
  As e doubles, the band SHRINKS by half and the band CENTERS k/e refine those of b.
  Conjecture: the lonely set L(E) is determined up to controlled error by the ODD
  speeds, because even speeds are "doubling refinements".  We MEASURE the exact
  contribution of even speeds to the gap loss.

EXACT: lonely set = [0,1) minus union of bands [k/v - 1/(14v), k/v + 1/(14v)].
"""
import sys
from fractions import Fraction as F
from math import gcd
sys.stdout.reconfigure(line_buffering=True) if hasattr(sys.stdout,'reconfigure') else None

def lonely_meas_and_free(E):
    """Exact lonely measure on [0,1) for speed set E, target 1/14. Returns (meas, free_intervals)."""
    bands=[]
    for v in E:
        if v==0: continue
        w=F(1,14*v)
        for k in range(0, v+1):
            c=F(k,v)
            blo=max(c-w,F(0)); bhi=min(c+w,F(1))
            if blo<bhi: bands.append((blo,bhi))
        # wraparound near 0 and 1 handled by k=0 and k=v giving [.,.] clamped
    bands.sort()
    merged=[]
    for s,e in bands:
        if merged and s<=merged[-1][1]:
            merged[-1]=(merged[-1][0],max(merged[-1][1],e))
        else: merged.append((s,e))
    covered=sum(e-s for s,e in merged)
    meas=F(1)-covered
    free=[]; prev=F(0)
    for s,e in merged:
        if prev<s: free.append((prev,s))
        prev=max(prev,e)
    if prev<F(1): free.append((prev,F(1)))
    return meas, free

def odd_part(v):
    while v%2==0: v//=2
    return v

print("="*84)
print("CRT/HALVING ANGLE: does the ODD sub-core control the lonely set? (EXACT)")
print("="*84)

# The drop-6 minimizer core (positive form) and consec.
C0=(1,2,3,4,5,7,8,9,10,11,12,13)   # drop-6 core
CONSEC=tuple(range(1,14))
cores={
  "drop-6 core C0":           C0,
  "consec {1..13}":           CONSEC,
  "drop-12 core":             (1,2,3,4,5,6,7,8,9,10,11,13),
}

print("\n(1) FULL lonely meas vs ODD-ONLY sub-core vs EVEN-ONLY sub-core")
print(f"   {'core':<22} {'L(full)':>14} {'L(odd-only)':>14} {'L(even-only)':>14} {'odds':>5} {'evens':>5}")
for name,E in cores.items():
    Lfull,_=lonely_meas_and_free(E)
    odds=sorted(v for v in E if v%2==1)
    evens=sorted(v for v in E if v%2==0)
    Lodd,_=lonely_meas_and_free(odds)
    Leven,_=lonely_meas_and_free(evens)
    print(f"   {name:<22} {float(Lfull):>14.8f} {float(Lodd):>14.8f} {float(Leven):>14.8f} {len(odds):>5} {len(evens):>5}")

print("\n   ==> If L(odd-only) == L(full), evens add NOTHING (odd sub-core controls).")
print("       If L(odd-only) > L(full), evens DO cut extra (doubling refinements bite).")

# (2) Incremental: add even speeds one at a time to the odd sub-core, measure each cut.
print("\n"+"-"*84)
print("(2) INCREMENTAL even-speed cuts (start from odd sub-core, add evens by Glaisher order)")
print("-"*84)
for name,E in cores.items():
    odds=sorted(v for v in E if v%2==1)
    evens=sorted((v for v in E if v%2==0), key=lambda v:(odd_part(v),v))
    cur=list(odds)
    L,_=lonely_meas_and_free(cur)
    print(f"\n   {name}: start odd-only L={float(L):.8f}")
    for e in evens:
        prev=L
        cur=sorted(cur+[e])
        L,_=lonely_meas_and_free(cur)
        b=odd_part(e); a=0; t=e
        while t%2==0: t//=2; a+=1
        print(f"      +{e:>2} (=2^{a}*{b}, doubling of odd {b}): L {float(prev):.8f} -> {float(L):.8f}  cut={float(prev-L):+.8f}")

# (3) THE HALVING REFINEMENT: does the band of 2v sit INSIDE the band-complement of v?
# Danger band of 2v at k/(2v) has half-width 1/(28v). For odd k, k/(2v) is a NEW center
# not coinciding with any j/v. Test whether these new centers land in v's SAFE gaps (cutting)
# or v's danger bands (redundant).
print("\n"+"-"*84)
print("(3) HALVING: where do the bands of 2v land relative to v's bands?")
print("    band(2v) centers = k/(2v). EVEN k -> j/v (refines v's own center, REDUNDANT).")
print("    ODD k -> new midpoint center. Does it fall in v's safe gap (CUTS) or v's band (redundant)?")
print("-"*84)
for v in [1,3,5,7]:
    new_centers=[F(k,2*v) for k in range(1,2*v,2)]  # odd k only
    # v's danger bands: [j/v - 1/(14v), j/v + 1/(14v)]
    in_safe=0; in_band=0
    for c in new_centers:
        # nearest j/v
        nearest_dist=min(abs(c-F(j,v)) for j in range(0,v+1))
        if nearest_dist > F(1,14*v): in_safe+=1
        else: in_band+=1
    print(f"   v={v}: {len(new_centers)} new (odd-k) centers from 2v; {in_safe} land in v's SAFE gap (potential cut), {in_band} inside v's band")
print("   ==> new centers in SAFE gaps are the ONLY way 2v cuts extra lonely measure.")

print("\nDONE.")
