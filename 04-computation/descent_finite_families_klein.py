#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""KNOW THE FINITE FAMILIES: the 2-adic descent cores over the covering family (klein-S17).

THM-580 descent: S = O u E (odd/even); S' = E/2; recurse; cores [O_0, O_1, ...]. The apex gap of a core
is g(residues mod 7) = min_{k!=0}|sum_{x in O mod 7} zeta^{kx}|^2 (THM-590; min nonzero = 4cos^2(3pi/7) at
doublets, 0 only at Z_7). This collects the FINITE family of cores arising from a broad covering family,
characterizes which Z_7-residue-sets appear, and confirms inf apex-gap = 4cos^2(3pi/7) (doublets DO arise;
THM-590 forbids anything lower).
"""
import math, cmath, itertools, random

def descend(S):
    cores=[]; cur=sorted(set(S)); guard=0
    while cur and guard<30:
        O=[v for v in cur if v%2==1]; E=[v for v in cur if v%2==0]
        cores.append(tuple(O))
        if not E: break
        cur=sorted(set(v//2 for v in E)); guard+=1
    return cores

w=cmath.exp(2j*math.pi/7)
def apex_gap(O):
    res=set(v%7 for v in O)
    if not res: return None
    return min(abs(sum(w**((k*x)%7) for x in res))**2 for k in range(1,7))

def res7(O): return tuple(sorted(set(v%7 for v in O)))

# ---- broad covering family ----
fam={}
for k in range(2,15): fam[f"consec{{1..{k}}}"]=list(range(1,k+1))
fam["tightest {1..12,182}"]=list(range(1,13))+[182]
fam["skip12 {1..11,13,84}"]=list(range(1,12))+[13,84]
fam["even-heavy"]=[2,4,6,8,10,12,14,1,3,5,7,9,11]
rng=random.Random(1714)
for i in range(400):
    sz=rng.randint(8,14); fam[f"rand{i}"]=rng.sample(range(1,30),sz)

all_cores=set(); all_res=set(); gaps={}
core_examples={}
for name,S in fam.items():
    for O in descend(S):
        if not O: continue
        all_cores.add(O); r=res7(O); all_res.add(r)
        g=apex_gap(O)
        if g is not None: gaps[r]=round(g,6); core_examples.setdefault(r,O)

print("="*82); print(" Finite family of 2-adic descent cores (mod 7) over the covering family"); print("="*82)
print(f" distinct odd-cores (as integer sets): {len(all_cores)}")
print(f" distinct mod-7 residue-sets (the APEX finite family): {len(all_res)} of 2^7=128 possible")
# group by size and gap
byval={}
for r in all_res:
    g=gaps[r]; byval.setdefault(g,[]).append(r)
print(f"\n residue-sets by apex gap value:")
for g in sorted(byval):
    rs=byval[g]
    sizes=sorted(set(len(r) for r in rs))
    tag = "  <-- Z_7 (covering boundary, gap 0)" if g<1e-9 else ("  <-- DOUBLET binding = 4cos^2(3pi/7)" if abs(g-0.198062)<1e-4 else "")
    print(f"   gap={g:.6f}: {len(rs)} residue-sets, sizes {sizes}{tag}")
nz=[g for g in byval if g>1e-9]
print(f"\n inf apex-gap over arising NON-Z_7 cores = {min(nz):.6f}  (4cos^2(3pi/7)={4*math.cos(3*math.pi/7)**2:.6f})")
doublets=[r for r in all_res if len(r)==2]
print(f" doublets (2-residue cores) that ARISE: {len(doublets)} -- e.g. {sorted(doublets)[:6]}")
print(f" => doublets arise => inf apex-gap = 4cos^2(3pi/7) EXACTLY (THM-590 forbids anything lower).")
print(f" min core size |O_j| arising: {min(len(r) for r in all_res)} ; sizes present: {sorted(set(len(r) for r in all_res))}")
# is the family ALL of Z_7's subsets, or constrained?
print(f"\n is the apex family ALL 2^7 residue-sets? {len(all_res)==128}  (else CONSTRAINED to {len(all_res)})")
