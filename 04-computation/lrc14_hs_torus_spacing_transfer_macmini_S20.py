#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ANGLE 7 -- DECISIVE transfer test in the HUFFER-SHEPP variable.

The HS proof technique (1987) is NOT majorization on speeds; it is a PAIRWISE
TRANSFER on the COVERING-ARC geometry: enlarge one arc, shrink another by e,
coverage does not decrease (Schur-CONVEX in arc lengths => spread-out wins).
The mechanism is the reflection lemma "a fixed-length interval is more likely
covered the closer it is to the center."

LRC ANALOGUE at the binding instant.  At x just off a/7, the 7 clocks paint the
7 sectors; W_a is the window where the painted UNION = all sectors.  As x sweeps,
each clock e moves at speed 7e.  The cell [a/7-1/14, a/7+1/14] of width 1/7 is
covered iff at every x in it the 7 residues {floor(7 e x)} hit all of {0..6}.

THE RIGHT HUFFER-SHEPP VARIABLE: at the central instant x=a/7 the clocks sit at
phases {e a mod 7}.  For FULL coverage to PERSIST as x moves, the sectors must be
covered with SLACK; the slack of sector s = the (signed) distance to the nearest
clock that will sweep INTO s.  The survival width W_a = min over the 6 binding
sectors of (slack / relative speed).  This is EXACTLY a "covering arc length"
min-functional.  The transfer that the HS proof uses = move a clock to enlarge
the smallest covering slack at the expense of a larger one.

FALSIFIABLE HYP-HS3 (the corrected sign):  measS7 = sum_a W_a behaves like
Huffer-Shepp coverage: it is maximized by the SPREAD-OUT covering-arc profile,
and consec achieves the max because the AP speed set {0,7,14,..} makes the clock
DRIFT SPEEDS {0,49,98,..7e} maximally spread (the e=0 clock = a stationary
infinite-slack covering arc; the rest a clean geometric ladder).  Concretely:

  TEST 1 (drift-speed spread): does measS7 INCREASE with the spread (variance /
    majorization) of the drift-speed multiset {7e}?  consec has the canonical
    {0,7,..,42} ladder.  HS predicts spread-out drift => more coverage.

  TEST 2 (the actual HS transfer): take consec, do a single Robin-Hood step that
    EQUALIZES two drift speeds (move e_i toward e_j).  HS predicts coverage DROPS.
    Do this for EVERY adjacent pair and confirm consec sits at a local Schur-max
    in the drift-speed variable (Schur-CONVEX: balanced=min, so equalizing drops).

  TEST 3 (global): over the full-residue stratum, correlate measS7 with the
    spread of the SLACK-NORMALIZED drift speeds; report whether consec is the
    UNIQUE Schur-maximal (most spread) point.
"""
import itertools, statistics, math
from fractions import Fraction as F
from math import gcd
from functools import reduce
from collections import defaultdict

def occupancy_pi(E):
    E=sorted(set(E)); bps={F(0),F(1)}
    for e in E:
        if e==0: continue
        for a in range(7*abs(e)+1): bps.add(F(a,7*abs(e)))
    bps=sorted(b for b in bps if 0<=b<=1); pi=[F(0)]*8
    for lo,hi in zip(bps,bps[1:]):
        if hi<=lo: continue
        xm=(lo+hi)/2; hit=set()
        for e in E:
            v=e*xm; v=v-(v.numerator//v.denominator); hit.add((v.numerator*7)//v.denominator)
        pi[len(hit)]+=hi-lo
    return pi
def measS7(E): return occupancy_pi(E)[7]
def primitive(E): return reduce(gcd,[e for e in E if e!=0],0)==1
def consec(k): return list(range(k))

def majorizes(x,y):
    xs=sorted(x,reverse=True); ys=sorted(y,reverse=True)
    if len(xs)!=len(ys) or sum(xs)!=sum(ys): return False
    px=py=F(0)
    for i in range(len(xs)-1):
        px+=xs[i]; py+=ys[i]
        if px<py: return False
    return True

if __name__=="__main__":
    print("="*78)
    print("TEST 1: correlate measS7 with SPREAD (variance) of drift speeds {7|e|}.")
    print("        HS-convex predicts POSITIVE corr (spread-out drift => more cover).")
    print("="*78)
    W=14
    bank=[(0,)+r for r in itertools.combinations(range(1,W+1),7)]
    bank=[E for E in bank if primitive(E) and set(e%7 for e in E)=={0,1,2,3,4,5,6}]
    cache={}
    def m_(E):
        if E not in cache: cache[E]=measS7(list(E))
        return cache[E]
    rows=[]
    for E in bank:
        drift=[7*abs(e) for e in E]          # speeds 7e ; e=0 -> 0 (stationary)
        rows.append((float(m_(E)), statistics.pvariance(drift), E))
    ms=[r[0] for r in rows]; vs=[r[1] for r in rows]
    n=len(ms); mm=sum(ms)/n; mv=sum(vs)/n
    cov=sum((a-mm)*(b-mv) for a,b in zip(ms,vs))/n
    sdm=math.sqrt(sum((a-mm)**2 for a in ms)/n); sdv=math.sqrt(sum((b-mv)**2 for b in vs)/n)
    corr=cov/(sdm*sdv) if sdm*sdv else 0
    print(f"  corr(measS7, drift-variance) = {corr:.4f}  ({len(bank)} full-res shapes)")
    rows.sort(reverse=True)
    print("  TOP 6 measS7:")
    for r in rows[:6]: print(f"     measS7={r[0]:.5f}  drift-var={r[1]:.1f}  {r[2]}")

    print()
    print("="*78)
    print("TEST 2 (the actual HS transfer): from consec, EQUALIZE two drift speeds")
    print("        (Robin-Hood toward balance). HS-convex predicts measS7 DROPS.")
    print("="*78)
    C=consec(8); mC=measS7(C)
    print(f"  consec drift speeds = {sorted(7*e for e in C)}  measS7={float(mC):.6f}")
    # a Robin-Hood step: pick i<j, move e_i up by t and e_j down by t (toward equal).
    # implement on the speed values themselves (keep distinct, keep set valid).
    drops=0; rises=0; eqs=0; ties=0; tested=0
    for i in range(len(C)):
        for j in range(len(C)):
            if i==j: continue
            ei,ej=C[i],C[j]
            if ei>=ej: continue
            # move ei up by 1, ej down by 1 => more balanced (closer together)
            E2=list(C); E2[i]=ei+1; E2[j]=ej-1
            if len(set(E2))!=len(E2): continue
            if min(E2)<0: continue
            if not primitive(E2): continue
            tested+=1
            d=measS7(E2)-mC
            if d<-F(1,10**15): drops+=1
            elif d>F(1,10**15): rises+=1
            else: eqs+=1
    print(f"  balancing transfers tested = {tested}")
    print(f"     measS7 DROPPED  (HS-convex predicts this) = {drops}")
    print(f"     measS7 ROSE     (violates HS-convex)       = {rises}")
    print(f"     measS7 UNCHANGED                            = {eqs}")
    verdict = "CONFIRMS HS-convex (balancing always drops)" if rises==0 and drops>0 else \
              ("REFUTES HS-convex (balancing sometimes rises)" if rises>0 else "INCONCLUSIVE")
    print(f"  TEST 2 verdict: {verdict}")

    print()
    print("="*78)
    print("TEST 3: is consec the UNIQUE drift-spread (Schur) maximum?  Over same")
    print("        drift-SUM groups, test majorization sign of measS7.")
    print("="*78)
    by_dsum=defaultdict(list)
    for E in bank:
        by_dsum[sum(7*abs(e) for e in E)].append(E)
    conv=0; conc=0; comp=0
    for s,grp in by_dsum.items():
        for Ea,Eb in itertools.combinations(grp,2):
            da=sorted((7*abs(e) for e in Ea),reverse=True)
            db=sorted((7*abs(e) for e in Eb),reverse=True)
            if majorizes(da,db):
                comp+=1; d=m_(Ea)-m_(Eb)
                if d>F(1,10**15): conv+=1
                if d<-F(1,10**15): conc+=1
            elif majorizes(db,da):
                comp+=1; d=m_(Eb)-m_(Ea)
                if d>F(1,10**15): conv+=1
                if d<-F(1,10**15): conc+=1
    print(f"  same-drift-sum majorization-comparable pairs = {comp}")
    print(f"  Schur-CONVEX  evidence (more drift-spread beat balanced) = {conv}")
    print(f"  Schur-CONCAVE evidence (balanced beat more drift-spread) = {conc}")
    if comp>0:
        if conc==0 and conv>0: print("  => measS7 SCHUR-CONVEX in drift-speed vector: MATCHES Huffer-Shepp!")
        elif conv==0 and conc>0: print("  => measS7 SCHUR-CONCAVE in drift-speed vector (opposite HS).")
        else: print(f"  => NEITHER pure (conv={conv}, conc={conc}); not clean Schur in drift vector.")
