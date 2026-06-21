#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ANGLE 3 REFINED -- the Freiman functional must live on the TORUS (opus 2026-06-21).

Finding from lrc14_freiman_doubling_vs_measS7: argmax measS7 == argmin INTEGER
doubling == consec, BUT the INTERIOR ranking disagrees (Spearman ~ -0.57): the
near-AP dilations [0,2,4,6,7,8,10,12], [0,1,3,5,7,9,11,13] rank 2nd-4th on measS7
yet have HIGH integer doubling. They are penalized by integer doubling but rewarded
by measS7. Common thread: they share consec's residue multiset {0,0,1,2,3,4,5,6}.

So the additive structure that measS7 actually rewards is NOT integer doubling.
measS7 is a covering measure on R/Z resolved into 7 sectors; the natural additive
ambient is the FREQUENCY side. The covering is controlled by the Fourier/L4 energy
of the spectrum. CANDIDATE FUNCTIONALS that respect the torus geometry:

  (F1) Residue-energy in Z/7: additive energy of the multiset {e mod 7}.
  (F2) Magnitude-weighted doubling: instead of |E+E|, count distinct VALUES but
       weight by smallness -> use the "spread" S(E)=sum of pairwise |gaps|.
  (F3) The TRUE governing quantity from HYP-2745/2750: the cover near a/7 is
       governed by min|e| per residue (the survival LEG). The aggregate functional
       is L(E) = sum over a of (cover width), whose local model is sum over residues
       of survival-leg = some symmetric function of {min|e| per residue}.
  (F4) L4-Fourier energy on Z (the additive energy IS the L4 of the indicator's
       transform) -- already tested, moderate.
  (F5) NEW: a "weighted additive energy" E_w = sum over (a,b,c,d), a+b=c+d, of a
       DECAY weight w(max|.|) -- rewards additive quadruples among SMALL elements.

DECISIVE TEST: which functional gives Spearman closest to -1 (or +1) AND gets the
full ranking right? If a torus-aware functional achieves near-perfect rank
correlation, the Freiman route is ALIVE (max measS7 <=> min torus-doubling => AP).
If NONE beat the integer doubling's -0.57, the additive-rigidity route is a
DEAD-END for the interior (only correct at the extreme point).
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

def occupancy_pi(E):
    E = sorted(set(int(e) for e in E))
    bps = {F(0), F(1)}
    for e in E:
        ae = abs(e)
        if ae == 0: continue
        for m in range(7*ae+1): bps.add(F(m, 7*ae))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    pi = [F(0)]*8
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        xm = (lo+hi)/2
        hit = set()
        for e in E:
            v = e*xm; v = v - (v.numerator // v.denominator)
            hit.add((v.numerator*7)//v.denominator)
        pi[len(hit)] += hi - lo
    return pi
def measS7(E): return occupancy_pi(E)[7]

# ---- functionals ----
def int_doubling(E): return F(len({a+b for a in E for b in E}), len(E))
def int_energy(E):
    c=defaultdict(int)
    for a in E:
        for b in E: c[a+b]+=1
    return sum(v*v for v in c.values())

def res7_energy(E):
    # additive energy of the residue multiset in Z/7 (multiset!)
    c=defaultdict(int)
    for a in E:
        for b in E: c[(a+b)%7]+=1
    return sum(v*v for v in c.values())

def spread(E):
    s=sorted(E); return s[-1]-s[0]   # diameter

def weighted_energy(E, decay=F(1,2)):
    # additive quadruples weighted by decay^max(|elements|): rewards small additive structure
    c=defaultdict(list)
    for a in E:
        for b in E: c[a+b].append((a,b))
    tot=F(0)
    for s,pairs in c.items():
        for (a,b) in pairs:
            for (cc,d) in pairs:
                m=max(abs(a),abs(b),abs(cc),abs(d))
                tot += decay**m
    return tot

def survival_leg_sum(E):
    # F3: sum over residue classes r of min|e| with e mod 7 == r (the slowest rep).
    # SMALLER = better cover (slow drift). Consec: legs {0:0,1:1,2:2,3:3,4:4,5:5,6:6}=21.
    byres=defaultdict(list)
    for e in E: byres[e%7].append(abs(e))
    return sum(min(byres[r]) for r in range(7) if r in byres)

def harmonic_leg(E):
    # the survival width is ~ 1/(slope) = 1/(7e); aggregate ~ sum 1/min|e|? Try sum of
    # reciprocal-of-largest-leg... actually width near a/7 ~ governed by the e that
    # drifts fastest leaving a residue. Use sum over residues of 1/(min nonzero |e|).
    byres=defaultdict(list)
    for e in E: byres[e%7].append(abs(e))
    tot=F(0)
    for r in range(7):
        if r not in byres: continue
        m=min(byres[r])
        tot += F(1, m) if m>0 else F(8)  # residue-0 doubled-zero gets a big bonus
    return tot

def residues(E): return frozenset(e%7 for e in E)
def is_full_residue(E): return residues(E)==frozenset(range(7))
def primitive(E): return reduce(gcd,[abs(e) for e in E if e!=0],0)==1
def consec(k): return list(range(k))

def spearman(xs,ys):
    def ranks(vals):
        idx=sorted(range(len(vals)),key=lambda i:vals[i]); r=[0.0]*len(vals); i=0
        while i<len(idx):
            j=i
            while j+1<len(idx) and vals[idx[j+1]]==vals[idx[i]]: j+=1
            avg=(i+j)/2.0+1
            for k in range(i,j+1): r[idx[k]]=avg
            i=j+1
        return r
    rx,ry=ranks(xs),ranks(ys); n=len(xs)
    mx=sum(rx)/n; my=sum(ry)/n
    cov=sum((rx[i]-mx)*(ry[i]-my) for i in range(n))
    sx=sum((rx[i]-mx)**2 for i in range(n))**.5; sy=sum((ry[i]-my)**2 for i in range(n))**.5
    return cov/(sx*sy) if sx and sy else float('nan')

if __name__=="__main__":
    print("#"*78); print("# REFINED TORUS FREIMAN FUNCTIONALS vs measS7"); print("#"*78)
    for k,W in [(8,12),(8,13)]:
        C=consec(k)
        bank=[(0,)+r for r in itertools.combinations(range(1,W+1),k-1)]
        full=[list(E) for E in bank if primitive(E) and is_full_residue(E)]
        print(f"\n{'='*70}\nk={k} span<= {W}: {len(full)} shapes")
        ms=[float(measS7(E)) for E in full]
        funcs = {
            "int_doubling (min)":      ([float(int_doubling(E)) for E in full], -1),
            "int_energy (max)":        ([float(int_energy(E)) for E in full], +1),
            "res7_energy (max)":       ([float(res7_energy(E)) for E in full], +1),
            "spread/diameter (min)":   ([float(spread(E)) for E in full], -1),
            "survival_leg_sum (min)":  ([float(survival_leg_sum(E)) for E in full], -1),
            "harmonic_leg (max)":      ([float(harmonic_leg(E)) for E in full], +1),
            "weighted_energy d=1/2":   ([float(weighted_energy(E,F(1,2))) for E in full], +1),
            "weighted_energy d=3/4":   ([float(weighted_energy(E,F(3,4))) for E in full], +1),
        }
        print(f"  {'functional':<28}{'Spearman(measS7,.)':>20}  {'argopt==consec?':>16}")
        for name,(vals,direction) in funcs.items():
            rho=spearman(ms,vals)
            if direction>0:
                opt=full[max(range(len(full)),key=lambda i:vals[i])]
            else:
                opt=full[min(range(len(full)),key=lambda i:vals[i])]
            print(f"  {name:<28}{rho:>+20.4f}  {str(tuple(opt)==tuple(C)):>16}")
        # Specifically: does survival_leg_sum or harmonic_leg break ties better than int?
        # Show the two near-AP runner-ups vs a low-doubling weak shape.
        probes=[[0,1,2,3,4,5,6,7],[0,2,4,6,7,8,10,12],[0,1,3,5,7,9,11,13],
                [0,1,2,3,4,5,6,8],[0,2,3,4,5,6,7,8]]
        print("\n  PROBE shapes (the interior-ranking battleground):")
        print(f"  {'E':<28}{'measS7':>9}{'intDoub':>9}{'res7En':>8}{'legSum':>8}{'harmLeg':>9}")
        for E in probes:
            if not (primitive(E) and is_full_residue(E)):
                tag="(out of stratum)"
            else: tag=""
            print(f"  {str(E):<28}{float(measS7(E)):>9.5f}{float(int_doubling(E)):>9.3f}"
                  f"{res7_energy(E):>8}{survival_leg_sum(E):>8}{float(harmonic_leg(E)):>9.3f} {tag}")
