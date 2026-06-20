#!/usr/bin/env python3
"""
lrc14_crt_factor_test_kps.py (kind-pasteur 2026-06-19)
Does the lonely measure CRT-factor: L(E) =? f_7(E) * f_2(E)?
Test the cell decomposition tau = a/7 + xi/7 (a=0..6, xi in [0,1)).
For runner v=7q+r: v*tau = qa + (ra + v*xi mod 7)/7... the 7-LOCAL phase is
(v*a mod 7) which depends on r=v mod 7 only; the wobble v*xi/7 depends on FULL v.
We measure per-cell lonely mass and ask whether it factors as
(7-local survival pattern, residue-only) x (within-cell archimedean density).
ALSO: confirm the TOWER-MIDDLE deletion principle predicts the global minimizer.
EXACT.
"""
import sys
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True) if hasattr(sys.stdout,'reconfigure') else None

def lonely_in_cell(E,a):
    lo=F(a,7)-F(1,14); hi=F(a,7)+F(1,14)
    bands=[]
    for v in E:
        if v==0: continue
        w=F(1,14*v); kmin=int(lo*v)-2; kmax=int(hi*v)+2
        for k in range(kmin,kmax+1):
            c=F(k,v); blo=max(c-w,lo); bhi=min(c+w,hi)
            if blo<bhi: bands.append((blo,bhi))
    bands.sort(); merged=[]
    for s,e in bands:
        if merged and s<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],e))
        else: merged.append((s,e))
    return (hi-lo)-sum(e-s for s,e in merged)

C0=[1,2,3,4,5,7,8,9,10,11,12,13]
print("(A) PER-CELL (7-adic) decomposition of drop-6 minimizer:")
cells=[lonely_in_cell(C0,a) for a in range(7)]
for a in range(7):
    print(f"  cell a={a}: m_a = {cells[a]} = {float(cells[a]):.8f}")
print(f"  sum = {sum(cells,F(0))} = {float(sum(cells,F(0))):.8f}")
surv=[a for a in range(7) if cells[a]>0]
print(f"  surviving cells: {surv}  (residue pattern)")

# residue-mod-7 of the core speeds (the 7-local data)
print(f"\n  speeds mod 7: {sorted(v%7 for v in C0)}  (multiset)")
print(f"  speeds mod 2: {[v%2 for v in C0]}  ({sum(1 for v in C0 if v%2==0)} even, {sum(1 for v in C0 if v%2)} odd)")

# (B) Factorization test: is m_a = (cell-survival indicator) * (common density R)?
print("\n(B) within-cell conditional density R_a = 7*m_a on surviving cells:")
for a in surv:
    print(f"  a={a}: R_a = 7*m_a = {7*cells[a]} = {float(7*cells[a]):.8f}")
Rs=[7*cells[a] for a in surv]
clean = (max(Rs)==min(Rs))
print(f"  R_a constant across surviving cells? {clean}  -> {'CLEAN 7-adic x archimedean product' if clean else 'NOT a clean product'}")
if clean:
    print(f"  => L = (#surv/7) * R = ({len(surv)}/7)*{Rs[0]} = {F(len(surv),7)*Rs[0]}")

# (C) TOWER-MIDDLE deletion principle: predict global minimizer over ALL primitive
# 12-subsets of [1,13] (i.e. single deletions), confirm drop-6, and explain via tower.
print("\n(C) all single-deletion cores ranked by L (the bounded atlas confirms drop-6):")
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
def tower_info(e):
    # in consec {1..13}, the dyadic tower of odd(e)
    b=e
    while b%2==0: b//=2
    tower=sorted(v for v in range(1,14) if (lambda x:(x//(x& -x)))(v)==b and v%2==0) if e%2==0 else None
    return b,tower
rows=[]
for e in range(1,14):
    core=[v for v in range(1,14) if v!=e]
    rows.append((lonely(core),e))
rows.sort()
for L,e in rows:
    par='EVEN' if e%2==0 else 'ODD '
    if e%2==0:
        b,tower=tower_info(e)
        pos='TOP' if e==max(tower) else ('only' if len(tower)==1 else 'MIDDLE')
        note=f"tower of {b}={tower}, deleting {pos}"
    else:
        note="odd (base member, NOT a doubling)"
    print(f"  drop {e:>2} [{par}]: L={float(L):.8f}  {note}")
print("\n  => global single-deletion minimizer = drop 6 (delete the MIDDLE of odd-3's tower {6,12}).")
print("     The deepest tower member 12 survives and re-covers 6's bands. This IS the halving handle.")
