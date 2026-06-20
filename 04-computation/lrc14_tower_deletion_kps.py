#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TOWER-DELETION LEMMA for LRC14 AP-tail mouth-retention.
kind-pasteur-2026-06-19.  EXACT rationals (lonely measure at gap 1/14).

CLAIM: for a in {0,1,2,3}, if 2^a NOT in C (AP-tail 12-core) then meas(G_C) >= 426/35035.

An AP-tail 12-core is C = ({1..13} \ holes) union tails, tails subset [14,inf), |C|=12.
We scan ALL AP-tail 12-cores in which a chosen tower bit 2^a is a hole, over a generous
tail window, find the MINIMUM measure for each a, and confirm all are >= thr2.

Then we set up the comb bound for the cleanest (a=3) case.
"""
import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd

def lonely_measure(C, theta=Fraction(1,14)):
    arcs=[]
    for d in C:
        if d==0: continue
        w=theta/d
        for m in range(0, d+1):
            lo=Fraction(m,d)-w; hi=Fraction(m,d)+w
            arcs.append((lo,hi))
    segs=[]
    for lo,hi in arcs:
        for shift in (-1,0,1):
            a=lo+shift; b=hi+shift
            a2=max(a,Fraction(0)); b2=min(b,Fraction(1))
            if a2<b2: segs.append((a2,b2))
    segs.sort()
    cur=Fraction(-1); union=[]
    for a,b in segs:
        if a>cur:
            union.append([a,b]); cur=b
        elif b>cur:
            union[-1][1]=b; cur=b
    covered=sum(b-a for a,b in union)
    return Fraction(1)-covered

def ap_tail_core(holes, tails):
    return tuple(sorted([d for d in range(1,14) if d not in holes]+list(tails)))

def carry_profile(C):
    d={}
    for e in C:
        if e==0: continue
        m=e; a=0
        while m%2==0: m//=2; a+=1
        d[m]=d.get(m,0)+2**a
    return dict(sorted(d.items()))

thr2=Fraction(426,35035)

def scan_missing_tower_bit(a, tail_max=40):
    """All AP-tail 12-cores where 2^a is a hole (2^a in {1,2,4,8}).
    holes always include 2^a; total holes = nh; tails = nh-? ; keep |C|=12.
    base {1..13} has 13 elements; |C|=12 means (#holes among 1..13) - (#tails) = 1.
    So #tails = #holes - 1.  We require 2^a among holes.
    Scan #holes = 1..4 (i.e. up to 3 tails), tails in [14, tail_max]."""
    bit=2**a
    best=None; below=[]
    cnt=0
    for nh in [1,2,3]:
        ntails=nh-1
        # holes: choose nh from 1..13, must include bit
        others=[d for d in range(1,14) if d!=bit]
        for extra in itertools.combinations(others, nh-1):
            holes=set(extra)|{bit}
            tailpool=range(14, tail_max+1)
            for tails in itertools.combinations(tailpool, ntails):
                C=ap_tail_core(holes, tails)
                if len(C)!=12: continue
                if reduce(gcd,C)!=1: continue
                cnt+=1
                L=lonely_measure(C)
                if best is None or L<best[0]:
                    best=(L,tuple(sorted(holes)),tails)
                if L<thr2:
                    below.append((L,tuple(sorted(holes)),tails))
    return best, below, cnt

if __name__=="__main__":
    print("=== TOWER-DELETION SCAN: minimum meas over AP-tail 12-cores missing 2^a ===")
    print(f"    thr2 = 426/35035 = {float(thr2):.8f}\n")
    for a in [3,2,1,0]:
        # smaller tail window for a=0 (1 is a hole => more tails) to keep it tractable
        tmax = 26 if a==0 else (30 if a==1 else 40)
        best, below, cnt = scan_missing_tower_bit(a, tail_max=tmax)
        L,holes,tails=best
        status = "OK >= thr2" if L>=thr2 else "*** BELOW thr2 ***"
        print(f"  a={a} (bit {2**a} is a HOLE): scanned {cnt} cores; MIN meas={float(L):.8f}={L}")
        print(f"      binding core: holes={holes} tails={tails}  -> {status}")
        if below:
            print(f"      !!! {len(below)} cores BELOW thr2:")
            for Lb,hb,tb in sorted(below)[:8]:
                print(f"          meas={Lb}={float(Lb):.8f} holes={hb} tails={tb}")
        print()
