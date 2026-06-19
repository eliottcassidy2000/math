#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ANGLE A part 4: probe a PROOF that consec maximizes p_0 = meas(S7).

p_0 = meas{x : every sector 1..6 is hit}  (sector 0 always hit by e=0).
Inclusion-exclusion over the 6 inner sectors:
   p_0 = sum_{A subset {1..6}} (-1)^{|A|} J(A,E),
   J(A,E) = meas{x: orbit avoids EVERY sector in A}.
This is exactly L_y with the FULL alternating g (the indicator 1[N=0]).
So p_0 = sum_r (-1)^r S_r (the true generating identity).

PROOF IDEAS TO TEST:
 (a) For a SINGLE point e: the set {x: frac(ex) in sector j} has measure 1/7 and is
     a union of e arcs each length 1/(7e). The orbit COVERS sector j iff some e does.
 (b) Consider the "covering map": as x sweeps [0,1), each e traces the circle e times.
     consec = {0,1,...,k-1}: speeds 0,1,..,k-1. The union of preimages.

We test a concrete COUPLING/REARRANGEMENT claim:
   Claim R: replacing any e_i in E by a value that makes the gaps more "consecutive"
   (a single adjacent-transposition toward consec) does not DECREASE p_0.
   I.e. p_0 is monotone along a path E -> ... -> consec of elementary moves.

Elementary move: pick E with a "gap" (some value v missing, v-1 in E, and a larger
value present); slide the largest element down into the gap. Track p_0.

We also test the SHIFT/COMPRESSION operad: does p_0(E) <= p_0(compress(E)) where
compress pushes elements toward {0,1,2,...}? Test all single-element decrements that
keep the set valid (distinct, contains 0).
kind-pasteur-2026-06-19 ANGLE-A.
"""
import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd

def frac(v): return v-(v.numerator//v.denominator)
def sect(v): return (v.numerator*7)//v.denominator

def p0_meas(E):
    E=sorted(set(E))
    bps=set([Fraction(0),Fraction(1)])
    for e in E:
        if e==0: continue
        for a in range(0,7*e+1):
            bps.add(Fraction(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1)
    tot=Fraction(0)
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        mid=(lo+hi)/2
        hit=set(sect(frac(e*mid)) for e in E)
        if all(jj in hit for jj in range(1,7)):
            tot+=(hi-lo)
    return tot

def single_decrements(E):
    """all valid sets obtained by decreasing one element by 1 (keep distinct, contain 0, sorted)."""
    Es=sorted(E); out=[]
    for i in range(len(Es)):
        v=Es[i]-1
        if v<0: continue
        if v in Es: continue
        new=sorted(Es[:i]+[v]+Es[i+1:])
        if new[0]!=0: continue
        out.append(tuple(new))
    return out

if __name__=="__main__":
    for k in [8,9]:
        print(f"\n{'='*70}\nk={k}: COMPRESSION-MONOTONICITY of p_0\n{'='*70}")
        C=tuple(range(k))
        p0c=p0_meas(C)
        print(f"consec p0={p0c}={float(p0c):.6f}")
        # test: for EVERY non-consec set (small spread), is there a single-element
        # decrement that does NOT decrease p0 and moves strictly toward consec?
        # Stronger: does EVERY valid single decrement never decrease p0?
        bad_any=0; bad_exists_increase_path=0; tested=0
        viol_examples=[]
        all_dec_monotone=True
        dec_viol=[]
        for combo in itertools.combinations(range(1,{8:13,9:12}[k]+1),k-1):
            E=(0,)+combo
            if E==C: continue
            if reduce(gcd,E)!=1: continue
            tested+=1
            p0E=p0_meas(E)
            decs=single_decrements(E)
            # does some decrement not decrease p0?
            vals=[(d,p0_meas(d)) for d in decs]
            # monotonicity: ALL decrements >= p0E ?
            for d,pd in vals:
                if pd<p0E:
                    all_dec_monotone=False
                    if len(dec_viol)<8: dec_viol.append((E,d,float(p0E),float(pd)))
        print(f"  tested {tested} non-consec sets (spread<= {{8:13,9:12}}[k])")
        print(f"  'EVERY single decrement never decreases p0': {'TRUE' if all_dec_monotone else 'FALSE'}")
        if dec_viol:
            print(f"  decrement violations (E -> decremented, p0E, p0dec):")
            for v in dec_viol: print(f"     {v}")
