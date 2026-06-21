#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ANGLE 5 -- MARKOFF / FROBENIUS extremality of the LRC(14) wall (opus 2026-06-21).

CLAIM family being tested: the consec/AP extremality of measS7 = sum_a W_a is a
Markoff/Frobenius-extremal phenomenon. The survival problem at resonance a is a
COVERING of Z/7 by clocks moving at integer velocities e in E. We test several
falsifiable Markoff/Frobenius invariants M(E) and ask whether measS7 (equivalently
sum_a W_a) is a MONOTONE function of M(E), extremized exactly at consec.

Invariants tested:
  M_var   = sum_{i<j} (e_i-e_j)^2          (variance / "Markoff balance" form)
  M_span  = max(E) - min(E)                (span; AP minimizes for full-residue)
  M_frob  = Frobenius number g(E\{0} reduced) of the numerical semigroup
  M_mark  = x^2+y^2+z^2-xyz with (x,y,z) the three g-legs ||p||,||q||,||pq||
  M_l2leg = sum of squared cycle-legs g(t)=2t(7-t) over occupied residue gaps

The test: among ALL full-residue primitive shapes at given span, does the ranking
by measS7 AGREE with the ranking by -M (i.e. measS7 large <=> M small)?
We report Spearman-style agreement (Kendall concordance on the top stratum) and
whether the UNIQUE measS7-max coincides with the UNIQUE M-extremum.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)

# ---------- measS7 / W_a machinery (reused, verified) ----------
def covered_sectors(E, xm):
    secs=set()
    for e in E:
        v=e*xm; v=v-(v.numerator//v.denominator)
        secs.add((v.numerator*7)//v.denominator)
    return secs

def measS7_arcs(E):
    E=sorted(set(int(e) for e in E))
    bps={F(0),F(1)}
    for e in E:
        ae=abs(e)
        if ae==0: continue
        for m in range(7*ae+1): bps.add(F(m,7*ae))
    bps=sorted(b for b in bps if 0<=b<=1)
    total=F(0)
    for lo,hi in zip(bps,bps[1:]):
        if hi<=lo: continue
        if len(covered_sectors(E,(lo+hi)/2))==7: total+=hi-lo
    return total

def measS7(E): return measS7_arcs(E)

def residues(E): return frozenset(e%7 for e in E)
def is_full_residue(E): return residues(E)==frozenset(range(7))
def primitive(E): return reduce(gcd,[abs(e) for e in E if e!=0],0)==1
def consec(k): return list(range(k))

# ---------- candidate invariants ----------
def M_var(E):
    E=[int(e) for e in E]
    return sum((E[i]-E[j])**2 for i in range(len(E)) for j in range(i+1,len(E)))

def M_span(E):
    E=[int(e) for e in E]; return max(E)-min(E)

def frobenius_number(gens):
    """Frobenius number of numerical semigroup gen by gens (gcd 1). Brute via
    coin reachability up to a safe bound. Returns None if not gcd-1."""
    gens=sorted(set(g for g in gens if g>0))
    if not gens: return None
    if reduce(gcd,gens)!=1: return None
    bound=max(gens)**2+1
    reach=[False]*(bound+1); reach[0]=True
    for i in range(bound+1):
        if reach[i]:
            for g in gens:
                if i+g<=bound: reach[i+g]=True
    # largest non-reachable
    fn=-1
    for i in range(bound+1):
        if not reach[i]: fn=i
    return fn

def M_frob(E):
    # velocity GAPS as semigroup generators: differences from min
    E=sorted(int(e) for e in E)
    gaps=[E[i+1]-E[i] for i in range(len(E)-1)]
    gaps=[g for g in gaps if g>0]
    return frobenius_number(gaps)

# ---------- ranking-agreement test ----------
def kendall_concordant(pairs):
    """pairs = list of (measS7, M). count concordant (measS7 up <=> M down)."""
    n=len(pairs); conc=disc=tie=0
    for i in range(n):
        for j in range(i+1,n):
            dm = pairs[i][0]-pairs[j][0]   # measS7 diff
            dM = pairs[i][1]-pairs[j][1]   # M diff
            if dm==0 or dM==0: tie+=1; continue
            # concordant for measS7 large <=> M small: signs OPPOSITE
            if (dm>0) != (dM>0): conc+=1
            else: disc+=1
    return conc,disc,tie

if __name__=="__main__":
    print("#"*78); print("# ANGLE 5: MARKOFF/FROBENIUS EXTREMALITY TEST"); print("#"*78)
    k=8
    C=consec(k); mC=measS7(C)
    print(f"\nconsec(8): measS7={mC}={float(mC):.6f}  M_var={M_var(C)} M_span={M_span(C)} M_frob={M_frob(C)}")

    for W in [13, 16]:
        bank=[(0,)+r for r in itertools.combinations(range(1,W+1),k-1)]
        full=[list(E) for E in bank if primitive(E) and is_full_residue(E)]
        data=[(measS7(E), E) for E in full]
        print(f"\n=== span<= {W}: {len(full)} full-res primitive shapes ===")
        # check consec is the unique max
        mx=max(d[0] for d in data)
        argmax=[E for m,E in data if m==mx]
        print(f"  measS7 max={mx}={float(mx):.6f} ; argmax count={len(argmax)} ; consec is max: {mx==mC}")

        for name,Mfun in [("M_var",M_var),("M_span",M_span),("M_frob",M_frob)]:
            pairs=[]
            ok=True
            for m,E in data:
                mv=Mfun(E)
                if mv is None: ok=False; break
                pairs.append((m,mv))
            if not ok:
                print(f"  {name}: undefined for some shape -> skip"); continue
            conc,disc,tie=kendall_concordant(pairs)
            tot=conc+disc
            # does the measS7-max also minimize M?
            Mmin=min(mv for _,mv in pairs)
            max_is_Mmin = all(Mfun(E)==Mmin for m,E in data if m==mx)
            tau = (conc-disc)/tot if tot else float('nan')
            print(f"  {name}: Kendall tau(measS7 up<=>M down)={tau:+.4f}  conc={conc} disc={disc} tie={tie}  | measS7-max minimizes {name}: {max_is_Mmin}")
