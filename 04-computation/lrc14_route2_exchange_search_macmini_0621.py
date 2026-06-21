#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ROUTE 2 -- EXCHANGE LEMMA search. consec = {0,1,...,k-1}. Any other full-residue
shape is obtained by REPLACING some magnitudes by larger ones (keeping residues=Z/7,
one residue doubled at residue (k-1 mod 7) since consec doubles residue 0 with 0,k-... 
wait consec={0,...,k-1}: residues are 0..k-1 mod 7. For k=8: residues {0,1,2,3,4,5,6,0}
-> residue 0 doubled (by 0 and 7). k=9: {0..8 mod7}={0,1,2,3,4,5,6,0,1} -> residues 0,1
doubled. k=10: residues 0,1,2 doubled.

The 'AP refill order' = consec uses the SMALLEST possible magnitudes for each residue.
Any deviation REPLACES a small magnitude by a larger one (in the same residue class, to
stay full-residue) OR reorders which residues get doubled.

EXCHANGE TEST: from consec, do a single 'lift' (increase one element by 7, staying full-
residue & primitive & sorted-distinct). Measure delta WIN. If EVERY single lift strictly
decreases WIN, and WIN is 'lift-monotone' (any shape reachable from consec by lifts, each
decreasing), that proves consec-max via a monotone descent argument.

A 'lift' on element e (residue r=e%7): replace e by e+7 (same residue). To keep DISTINCT
we may need e+7 not already present; if present, this is a different move. We enumerate
all full-residue shapes and check: is the partial order by 'sum of magnitudes' (or by the
lift poset) compatible with WIN being monotone decreasing from consec?
"""
import sys, itertools
from fractions import Fraction as F
from collections import Counter
sys.stdout.reconfigure(line_buffering=True)
HALF=F(1,14)
def sector_of(e,a,y):
    pos=F(e*a)+F(7*e)*y; return (pos.numerator//pos.denominator)%7
def covered(E,a,y): return len({sector_of(e,a,y) for e in E})==7
def breakpoints(E,a):
    bps={F(0),-HALF,HALF}
    for e in E:
        if e==0: continue
        lo_val=F(7*e)*(-HALF)+F(e*a); hi_val=F(7*e)*(HALF)+F(e*a)
        lo_i=min(lo_val,hi_val);hi_i=max(lo_val,hi_val); m=lo_i.numerator//lo_i.denominator
        while m<=hi_i.numerator//hi_i.denominator+1:
            y=F(m-e*a,7*e)
            if -HALF<=y<=HALF: bps.add(y)
            m+=1
    return sorted(bps)
def window_TpTm(E,a):
    bps=breakpoints(E,a); ivals=list(zip(bps,bps[1:])); Tp=F(0)
    for lo,hi in ivals:
        if hi<=0: continue
        lo2=max(lo,F(0))
        if covered(E,a,(lo2+hi)/2): Tp=hi
        else:
            if lo2==F(0): Tp=F(0)
            break
    Tp=min(Tp,HALF); Tm=F(0)
    for lo,hi in reversed(ivals):
        if lo>=0: continue
        hi2=min(hi,F(0))
        if covered(E,a,(lo+hi2)/2): Tm=-lo
        else:
            if hi2==F(0): Tm=F(0)
            break
    Tm=min(Tm,HALF); return Tp,Tm
def WIN(E): return sum(sum(window_TpTm(E,a)) for a in range(1,7))
def is_full_residue(E): return frozenset(e%7 for e in E)==frozenset(range(7))
def consec(k): return list(range(k))

def single_lift_test(k, span):
    """From consec, enumerate all full-residue shapes; check if WIN is monotone under
       the 'cover'/lift partial order. Concretely: does there exist a chain from consec
       to every shape along which WIN strictly decreases at each step?
       Simpler robust test: define poset by componentwise <= on SORTED magnitude vector.
       Is WIN monotone-decreasing wrt this poset (consec is the bottom)?"""
    C=consec(k); Cwin=WIN(C)
    bank=[(0,)+c for c in itertools.combinations(range(1,span+1),k-1)]
    bank=[E for E in bank if is_full_residue(E)]
    # poset: E1 <= E2 if sorted(E1) componentwise <= sorted(E2)
    def le(E1,E2): return all(x<=y for x,y in zip(sorted(E1),sorted(E2)))
    # Check WIN monotone: for all comparable pairs E1<=E2 (E1!=E2), WIN(E1) >= WIN(E2)?
    bankl=[list(E) for E in bank]
    wins={tuple(E):WIN(E) for E in bankl}
    viol=0; checked=0; examples=[]
    for i in range(len(bankl)):
        for j in range(len(bankl)):
            if i==j: continue
            E1,E2=bankl[i],bankl[j]
            if le(E1,E2):
                checked+=1
                if wins[tuple(E1)] < wins[tuple(E2)]:
                    viol+=1
                    if len(examples)<5: examples.append((E1,E2,wins[tuple(E1)],wins[tuple(E2)]))
    print(f"k={k} span<={span}: comparable pairs={checked}, WIN-monotone violations={viol}")
    for E1,E2,w1,w2 in examples:
        print(f"   VIOL {E1}<= {E2} but WIN {float(w1):.6f} < {float(w2):.6f}")
    return viol

for k,span in [(8,14),(9,16),(10,18)]:
    single_lift_test(k,span)
