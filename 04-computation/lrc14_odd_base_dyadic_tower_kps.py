#!/usr/bin/env python3
"""
lrc14_odd_base_dyadic_tower_kps.py (kind-pasteur 2026-06-19)
The odd sub-core {1,3,5,7,9,11,13} is the SAME for every [1,13]\{even} core.
Its lonely measure is a CONSTANT base. Even speeds = dyadic doublings refine it.
Investigate: (a) exact value of odd-base L; (b) the dyadic tower 1->2->4->8 etc.;
(c) whether the cap is reached purely by which doublings are present.
EXACT.
"""
import sys
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True) if hasattr(sys.stdout,'reconfigure') else None

def lonely(E):
    bands=[]
    for v in E:
        if v==0: continue
        w=F(1,14*v)
        for k in range(0,v+1):
            c=F(k,v); blo=max(c-w,F(0)); bhi=min(c+w,F(1))
            if blo<bhi: bands.append((blo,bhi))
    bands.sort(); merged=[]
    for s,e in bands:
        if merged and s<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],e))
        else: merged.append((s,e))
    return F(1)-sum(e-s for s,e in merged)

ODD=[1,3,5,7,9,11,13]
Lodd=lonely(ODD)
print(f"odd-base core {ODD}")
print(f"  L(odd-base) = {Lodd} = {float(Lodd):.10f}")
print(f"  as fraction: {Lodd.numerator}/{Lodd.denominator}")

# Is the odd base the SAME for any [1,13]\{even}? The odd speeds are always {1,3,5,7,9,11,13}.
# So yes by construction. The variation is purely in the even/doubling content.
print("\nEVEN content of each [1,13]\{e} core (e even) -- the doubling multiset:")
def dyad(v):
    a=0; b=v
    while b%2==0: b//=2; a+=1
    return (b,a)
for e in [2,4,6,8,10,12]:
    evens=[v for v in range(1,14) if v%2==0 and v!=e]
    print(f"  drop {e:>2}: evens={evens}  dyadic={[ (dyad(v)) for v in evens]}")

# Full lonely for each, and how much each even contributes
print("\nFull L and the even-removal sensitivity (drop each even from consec, measure rise):")
CONSEC=list(range(1,14))
Lfull=lonely(CONSEC)
print(f"  consec full L={float(Lfull):.10f} = {Lfull}")
for e in [2,4,6,8,10,12]:
    core=[v for v in range(1,14) if v!=e]
    Lc=lonely(core)
    print(f"  drop {e:>2}: L={float(Lc):.10f} = {Lc}   (rise from consec = {float(Lc-Lfull):+.8f})")

# THE TOWER: doubling chain 1->2->4->8 and 3->6->12 and 5->10. Build lonely with progressively
# more of the dyadic tower present, to see the "halving cascade" geometry.
print("\nDYADIC TOWER CASCADE (odd base + progressively deeper doublings):")
ODD=[1,3,5,7,9,11,13]
towers={1:[2,4,8],3:[6,12],5:[10]}
# add by 2-adic LEVEL across all odds: level1={2,6,10}, level2={4,12}, level3={8}
levels={1:[2,6,10],2:[4,12],3:[8]}
cur=list(ODD); L=lonely(cur)
print(f"  level 0 (odd base):            L={float(L):.10f}")
for lev in [1,2,3]:
    cur=sorted(cur+levels[lev]); L=lonely(cur)
    print(f"  + level {lev} {str(levels[lev]):<12} -> L={float(L):.10f}")
print("  (consec {1..13} = odd base + all dyadic levels => L=0)")

# Does each dyadic LEVEL roughly HALVE the lonely measure? (the 'halving' in 14=2*7)
print("\nHALVING CHECK: ratio L(after level k)/L(after level k-1):")
cur=list(ODD); prev=lonely(cur); seq=[prev]
for lev in [1,2,3]:
    cur=sorted(cur+levels[lev]); L=lonely(cur); seq.append(L)
for i in range(1,len(seq)):
    print(f"  level {i}: L={float(seq[i]):.8f}  ratio to prev = {float(seq[i]/seq[i-1]):.5f}")
