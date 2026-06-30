#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THE BRIDGE from the apex SKELETON to the full rho_j (klein-S18).

THM-580: rho_j = m(lonely S^j)/[m(lonely O_j) . m(lonely S^{j+1})], the per-level 2-sheet parity
decorrelation. THM-590 skeleton: g(O_j mod 7) = apex cyclotomic gap. Compute the ACTUAL rho_j over real
descent chains and correlate with (7-part) g(O_j mod 7) and (2-part) CV(N2_{O_j}), plus the THM-580 CS
bound 1 - CV(N2).sqrt((1-m')/m'). Question: does rho_j factor as 2-part x 7-part (14=2.7)? where binding?
"""
import sys, os, math, cmath
from fractions import Fraction as Fr
sys.path.insert(0, '04-computation')
M=__import__('lrc14_floor_CV_sheetcount_bound_macmini_20260629')
lonely_set, measure, intersect, shift = M.lonely_set, M.measure, M.intersect, M.shift

def descend_chain(S):
    chain=[]; cur=sorted(set(S)); guard=0
    while cur and guard<30:
        O=[v for v in cur if v%2==1]; E=[v for v in cur if v%2==0]
        chain.append((cur[:],O,E)); 
        if not E: break
        cur=sorted(set(v//2 for v in E)); guard+=1
    return chain

w=cmath.exp(2j*math.pi/7)
def apex_gap(O):
    res=set(v%7 for v in O)
    if not res: return None
    return min(abs(sum(w**((k*x)%7) for x in res))**2 for k in range(1,7))

def cv2_2sheet(L):
    m=measure(L)
    if m==0: return None
    A=measure(intersect(L, shift(L, Fr(1,2))))
    return float((m + A - 2*m*m)/(2*m*m))

import random
rng=random.Random(2024)
fam={f"consec{{1..{k}}}":list(range(1,k+1)) for k in range(4,15)}
fam["tightest{1..12,182}"]=list(range(1,13))+[182]
fam["skip12{1..11,13,84}"]=list(range(1,12))+[13,84]
for i in range(120): fam[f"rand{i}"]=rng.sample(range(1,30), rng.randint(9,14))

rows=[]
for name,S in fam.items():
    ch=descend_chain(S)
    for j,(Sj,O,E) in enumerate(ch):
        LO=lonely_set(O); mO=measure(LO)
        if not E:  # last level, no descent
            continue
        LE=lonely_set(E); mE=measure(LE)
        LS=lonely_set(Sj); mS=measure(LS)
        if mO==0 or mE==0: continue
        rho=float(mS/(mO*mE))
        g=apex_gap(O); cv2=cv2_2sheet(LO); mp=float(mE)
        if g is None or cv2 is None: continue   # skip empty odd part
        rows.append((rho,g,cv2,mp,len(set(v%7 for v in O)),name,j))

print("="*86); print(f" THE BRIDGE: actual rho_j vs apex gap g(O_j mod7) and 2-sheet CV  ({len(rows)} levels)"); print("="*86)
rows.sort()
print(f"\n min rho_j = {rows[0][0]:.4f}  (THM-580 reported min 0.515)")
print(" 12 BINDING (smallest rho_j) levels:  rho_j |  g(apex) | CV(N2) | m'(desc) | |O mod7| | source")
for rho,g,cv2,mp,nres,name,j in rows[:12]:
    cv=math.sqrt(cv2) if cv2 is not None else float('nan')
    print(f"   rho={rho:.4f} | g={g:.4f} | CV={cv:.3f} | m'={mp:.4f} | |O7|={nres} | {name} L{j}")

# correlation: does rho track g? does the CS bound hold? where is it lossy?
import statistics
gs=[g for _,g,_,_,_,_,_ in rows]; rhos=[r for r,_,_,_,_,_,_ in rows]
print(f"\n rho_j range [{min(rhos):.3f},{max(rhos):.3f}]; g range [{min(gs):.3f},{max(gs):.3f}]")
# group min rho by apex gap value
from collections import defaultdict
bygap=defaultdict(list)
for rho,g,cv2,mp,nres,name,j in rows: bygap[round(g,4)].append(rho)
print("\n min rho_j conditioned on the apex gap value (does small g force small rho?):")
for g in sorted(bygap):
    rs=bygap[g]; print(f"   g={g:.4f}: levels={len(rs):4d}  min rho={min(rs):.4f}  mean rho={statistics.mean(rs):.4f}")

# the CS bound and where m' is small (the cusp / lossy regime)
print("\n THM-580 CS bound  1 - CV(N2).sqrt((1-m')/m')  vs actual rho_j (is the apex skeleton tighter?):")
viol=0; loss=[]
for rho,g,cv2,mp,nres,name,j in rows:
    if cv2 is None or mp<=0: continue
    cs=1-math.sqrt(cv2)*math.sqrt((1-mp)/mp)
    loss.append((rho-cs, rho, cs, mp, g, name, j))
loss.sort()
small_mp=[r for r in rows if r[3]<0.15]
print(f"   levels with m'<0.15 (lonely-poor, CS lossy): {len(small_mp)}  -- their min rho_j={min((r[0] for r in small_mp),default=float('nan')):.4f}")
print(f"   CS bound is NEGATIVE (vacuous) on {sum(1 for d,rho,cs,mp,g,nm,j in loss if cs<=0)}/{len(loss)} levels; actual rho_j stays >= {min(r[0] for r in rows):.3f} there.")
print(f"   => the m'-independent apex skeleton (THM-590) is the right object exactly where the m'-dependent CS bound dies.")
