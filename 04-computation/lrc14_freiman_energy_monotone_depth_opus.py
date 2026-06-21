#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ANGLE 3 -- HOW DEEP DOES THE ENERGY<->measS7 ORDER GO? (opus 2026-06-21).

Apex coincidence verified k=8,9,10: argmax measS7 == argmax additive_energy ==
consec (unique AP). At k=9,10 even the 2nd place agrees. This script measures the
EXACT DEPTH of the order-isomorphism between additive energy and measS7, and runs
the decisive falsification:

  (R1) MONOTONE-IN-TOP-TIER: among shapes with energy >= threshold, is measS7
       strictly order-preserving? Find the FIRST inversion pair (E1,E2):
       energy(E1) > energy(E2) but measS7(E1) < measS7(E2). Report its energy gap
       and measS7 gap -- a SMALL inversion deep in the tail is harmless; an
       inversion NEAR THE TOP would threaten the bridge.
  (R2) TOP-t agreement: for t=1,2,3,5,10, do the top-t by energy == top-t by measS7
       (as sets)? Report the depth at which they first diverge.
  (R3) CONDITIONAL bridge: restrict to "AP-like" shapes (doubling <= 2.1). Is the
       energy<->measS7 order PERFECT there? (the inversions may all live among
       high-doubling dilated sets, which are exactly the torus-incommensurate ones
       Freiman theory would call "non-rigid").
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
def add_energy(E):
    c=defaultdict(int)
    for a in E:
        for b in E: c[a+b]+=1
    return sum(v*v for v in c.values())
def doubling(E): return F(len({a+b for a in E for b in E}), len(E))
def residues(E): return frozenset(e%7 for e in E)
def is_full_residue(E): return residues(E)==frozenset(range(7))
def primitive(E): return reduce(gcd,[abs(e) for e in E if e!=0],0)==1
def consec(k): return list(range(k))

if __name__=="__main__":
    print("#"*78); print("# ENERGY <-> measS7 ORDER DEPTH + INVERSION HUNT"); print("#"*78)
    for k,W in [(8,13),(9,13),(10,14)]:
        C=consec(k)
        bank=[(0,)+r for r in itertools.combinations(range(1,W+1),k-1)]
        full=[list(E) for E in bank if primitive(E) and is_full_residue(E)]
        if not full:
            print(f"\nk={k} W={W}: empty"); continue
        rows=[{'E':E,'m':measS7(E),'e':add_energy(E),'d':doubling(E)} for E in full]
        print(f"\n{'='*70}\nk={k} span<= {W}: {len(rows)} shapes")
        # R2: top-t agreement
        by_m=sorted(rows,key=lambda r:-r['m']); by_e=sorted(rows,key=lambda r:(-r['e'],-r['m']))
        for t in [1,2,3,5,10]:
            sm={tuple(r['E']) for r in by_m[:t]}; se={tuple(r['E']) for r in by_e[:t]}
            print(f"  top-{t}: |energy-top ∩ measS7-top| = {len(sm&se)}/{t}")
        # R1: first inversion (sorted by energy desc, find earliest measS7 drop-violation)
        # An inversion: a pair where energy strictly orders one way, measS7 the other.
        # Find the inversion involving the HIGHEST-energy shape possible.
        s=sorted(rows,key=lambda r:-r['e'])  # energy desc
        first_inv=None
        for i in range(len(s)):
            for j in range(i+1,len(s)):
                if s[i]['e']>s[j]['e'] and s[i]['m']<s[j]['m']:
                    first_inv=(s[i],s[j]); break
            if first_inv: break
        if first_inv:
            a,b=first_inv
            rank_a=s.index(a); rank_b=s.index(b)
            print(f"  FIRST inversion (highest-energy shape involved, energy-rank {rank_a} vs {rank_b}):")
            print(f"    E1={a['E']} energy={a['e']} measS7={float(a['m']):.6f} doub={float(a['d']):.3f}")
            print(f"    E2={b['E']} energy={b['e']} measS7={float(b['m']):.6f} doub={float(b['d']):.3f}")
            print(f"    (energy higher but measS7 LOWER -> inversion at energy-rank {rank_a})")
        else:
            print("  NO inversion: energy and measS7 are PERFECTLY order-isomorphic!")
        # R3: restrict to low-doubling (AP-like)
        for cap in [F(21,10), F(23,10)]:
            sub=[r for r in rows if r['d']<=cap]
            inv=False; worst=None
            sd=sorted(sub,key=lambda r:-r['e'])
            for i in range(len(sd)):
                for j in range(i+1,len(sd)):
                    if sd[i]['e']>sd[j]['e'] and sd[i]['m']<sd[j]['m']:
                        inv=True; worst=(sd[i],sd[j]); break
                if inv: break
            print(f"  doubling<= {float(cap):.2f}: {len(sub)} shapes, energy<->measS7 order PERFECT? {not inv}")
            if inv:
                a,b=worst
                print(f"      inversion: {a['E']}(e={a['e']},m={float(a['m']):.4f}) vs {b['E']}(e={b['e']},m={float(b['m']):.4f})")
