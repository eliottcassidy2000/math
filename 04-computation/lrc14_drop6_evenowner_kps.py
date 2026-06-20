#!/usr/bin/env python3
"""
lrc14_drop6_evenowner_kps.py (kind-pasteur 2026-06-19)
WHY is drop-6 the minimizer through the EVEN/ODD (halving) lens?
- odd base L = 75454/315315 is fixed.
- adding the 5 surviving evens cuts it to the final L; drop-e leaves out one even.
- drop-6 gives the SMALLEST final L (7/858). We show: removing 6 = 2*3 is the LEAST
  damaging deletion because 6's danger bands are the most REDUNDANTLY covered by the
  other surviving even speeds (esp. 12=2*6's halving-tower parent, and odd 3).
- Examine the 4 critical surviving intervals of drop-6 (from HYP-2651): who owns each
  boundary, and is the owner EVEN or ODD?
EXACT.
"""
import sys
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True) if hasattr(sys.stdout,'reconfigure') else None

def lonely_free(E):
    bands=[]
    for v in E:
        if v==0: continue
        w=F(1,14*v)
        for k in range(0,v+1):
            c=F(k,v); blo=max(c-w,F(0)); bhi=min(c+w,F(1))
            if blo<bhi: bands.append((blo,bhi,v,k))
    bands.sort()
    merged=[]
    for s,e,v,k in bands:
        if merged and s<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],e))
        else: merged.append((s,e))
    free=[]; prev=F(0)
    for s,e in merged:
        if prev<s: free.append((prev,s))
        prev=max(prev,e)
    if prev<F(1): free.append((prev,F(1)))
    meas=F(1)-sum(e-s for s,e in merged)
    return meas, free, bands

C0=[1,2,3,4,5,7,8,9,10,11,12,13]
meas,free,bands=lonely_free(C0)
print(f"drop-6 core {C0}: L={meas}={float(meas):.8f}, {len(free)} free intervals")
for s,e in free:
    Lown=[(v,k) for blo,bhi,v,k in bands if bhi==s]
    Rown=[(v,k) for blo,bhi,v,k in bands if blo==e]
    def parity(owns): return ",".join(f"{v}({'E' if v%2==0 else 'O'},k={k})" for v,k in owns)
    print(f"  ({s} , {e}) len {e-s}")
    print(f"     LEFT owner: {parity(Lown)}    RIGHT owner: {parity(Rown)}")

# Now: who owns the boundaries -- count EVEN vs ODD owners across the minimizer.
print("\nBoundary-owner parity tally for drop-6 minimizer:")
owners=[]
for s,e in free:
    for blo,bhi,v,k in bands:
        if bhi==s or blo==e: owners.append(v)
from collections import Counter
cnt=Counter(owners)
ev=sum(c for v,c in cnt.items() if v%2==0); od=sum(c for v,c in cnt.items() if v%2==1)
print(f"  owners multiset: {dict(sorted(cnt.items()))}")
print(f"  EVEN-owner boundary count = {ev}, ODD-owner = {od}")

# Compare: for EACH drop-e (e even), what is the final L and what even is removed?
print("\nDrop-each-even ledger: which even speed, its dyadic type, resulting L:")
def lonely(E):
    return lonely_free(E)[0]
def dyad(v):
    a=0;b=v
    while b%2==0:b//=2;a+=1
    return f"2^{a}*{b}"
rows=[]
for e in [2,4,6,8,10,12]:
    core=[v for v in range(1,14) if v!=e]
    L=lonely(core); rows.append((L,e))
    print(f"  drop {e:>2} (={dyad(e)}): L={float(L):.8f}={L}")
rows.sort()
print(f"\n  MINIMIZER: drop {rows[0][1]} (={dyad(rows[0][1])}) with L={rows[0][0]}")
print(f"  Note 6=2*3 is the unique LEVEL-1 doubling of 3 whose REMOVAL is cheapest.")

# Critical test: 6=2*3. Its halving-parent in the tower is 12=2^2*3=2*6. With 12 present,
# the doubling-tower of odd 3 still has its deepest member -> 6's removal is well-covered.
# Contrast 8=2^3*1: removing 8 (drop-8, L=0.0452, the MAX) is most damaging because 8 is
# the TOP of the odd-1 tower {2,4,8}; nothing else covers its fine bands.
print("\nTOWER-TOP vs TOWER-MIDDLE removal test:")
print("  odd 1 tower in consec: {2,4,8}  -> removing TOP 8 is worst (L=0.0452)")
print("  odd 3 tower in consec: {6,12}   -> removing 6 leaves 12 (deeper); cheapest (L=0.0082)")
print("  odd 5 tower in consec: {10}     -> removing 10 (only member) gives L=0.0241")
