#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Understanding the FAMILY INFIMUM of the floor R' = m_S/(m_R m_Q)  (klein-S16).

The covering family is INFINITE, so a scan only UPPER-bounds the true inf. Probe its structure:
(1) single-drops R={1..13}\{j}, Q={1,2}: R' vs which speed j is dropped (does apex-7 bind?);
(2) the cusp approach: sequences with m_R -> 0; does R' -> positive limit or 0? is the inf attained?
(3) expand: pair-drops, speeds>13, other Q -- does anything beat the 0.344 scan-inf?
Contrast with the FINITE Z_7-core family (inf rho_j = 4cos^2(3pi/7), PROVED THM-590).
"""
import sys, os, math, itertools
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
M = __import__("lrc14_floor_CV_sheetcount_bound_macmini_20260629")
lonely_set, measure = M.lonely_set, M.measure

def Rp(R,Q):
    R=tuple(sorted(set(x for x in R if x%14!=0))); Q=tuple(sorted(set(Q)))
    mR=measure(lonely_set(R)); mQ=measure(lonely_set(Q))
    if mR==0 or mQ==0: return None,mR,mQ
    S=tuple(sorted(set(R)|set(14*q for q in Q)))
    return measure(lonely_set(S))/(mR*mQ), mR, mQ

print("="*78); print(" (1) single-drops R={1..13}\\{j}, Q={1,2}: R' vs dropped speed j"); print("="*78)
base=list(range(1,14))
rows=[]
for j in base:
    R=[x for x in base if x!=j]
    rp,mR,mQ=Rp(R,[1,2])
    rows.append((float(rp),j,float(mR)))
for rp,j,mR in sorted(rows):
    flag=" <-- APEX-7 (binding)" if j==7 else (" (=14-resonance partner)" if j in (7,) else "")
    print(f"   drop j={j:>2}: R'={rp:.4f}  m_R={mR:.4f}{flag}")
print(f"   => min over single-drops at j={min(rows)[1]} (R'={min(rows)[0]:.4f}); max at j={max(rows)[1]} (R'={max(rows)[0]:.4f})")

print("\n"+"="*78); print(" (2) cusp approach: m_R -> 0. R={1..13}\\{7} then re-add danger toward the boundary"); print("="*78)
# sequence: drop fewer / add speeds to push m_R down toward the {1..13} boundary (m_R=0)
seqs=[
  ("{1..13}\\{7}", [x for x in base if x!=7]),
  ("{1..13}\\{7}+{15}", [x for x in base if x!=7]+[15]),
  ("{1..13}\\{7}+{15,16}", [x for x in base if x!=7]+[15,16]),
  ("{1..13}\\{7}+{15,16,17}", [x for x in base if x!=7]+[15,16,17]),
  ("{1..13}\\{7}+{15..20}", [x for x in base if x!=7]+[15,16,17,18,19,20]),
]
for name,R in seqs:
    rp,mR,mQ=Rp(R,[1,2])
    print(f"   R={name:<26} m_R={float(mR):.5f}  R'={(float(rp) if rp else float('nan')):.5f}")
print("   (watch whether R' -> a positive limit (floor bounded) or -> 0 (floor fails) as m_R shrinks)")

print("\n"+"="*78); print(" (3) does anything beat the 0.344 scan-inf? pair-drops + high speeds + other Q"); print("="*78)
best=(1e9,None,None)
cands=[]
# pair drops, Q={1,2}
for j,k in itertools.combinations(base,2):
    R=[x for x in base if x not in (j,k)]
    cands.append((R,[1,2]))
# single drop + add high speeds, Q={1,2} and {1,2,3}
for j in base:
    for extra in ([15],[15,16],[15,21],[15,16,17,18]):
        cands.append(([x for x in base if x!=j]+extra,[1,2]))
        cands.append(([x for x in base if x!=j]+extra,[1,2,3]))
for R,Q in cands:
    rp,mR,mQ=Rp(R,Q)
    if rp is None: continue
    if float(rp)<best[0]: best=(float(rp),tuple(sorted(set(R))),tuple(Q))
print(f"   min R' over {len(cands)} extra candidates: {best[0]:.5f}  at R={best[1]} Q={best[2]}")
print(f"   vs the {{1..13}}\\{{7}} scan-inf 0.34394: {'LOWER (true inf < scan)' if best[0]<0.3439 else 'NOT lower (0.344 robust here)'}")

print("\n"+"="*78); print(" CONTRAST: finite vs infinite family"); print("="*78)
print(f"   rho_j over Z_7-cores: FINITE family (2^7) -> inf = MIN = 4cos^2(3pi/7) = {4*math.cos(3*math.pi/7)**2:.5f} (PROVED, THM-590, attained at the doublet).")
print(f"   R'   over coverings : INFINITE family    -> scan only UPPER-bounds the true inf; inf is the m_R->0 CUSP LIMIT.")
print( "   The 2-adic descent FINITIZES R' into rho_j (infinite -> finite); that finitization is what makes the infimum PROVABLE.")
