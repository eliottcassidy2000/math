#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DOES the trace/pairing collapse actually act on the REAL relation lattice sum?

We KNOW (exact): K(n)=w(n) D7(c), w real, and K(n)+K(-n)=2 w(n) Re D7(c).
So the n<->-n PAIRING is a genuine lattice symmetry (lattice is symmetric) and
the correction sum_{n in Lambda} K(n) is REAL with imaginary parts gone.  GOOD.

The Galois DILATION c->a c is NOT a lattice symmetry: it changes residues without
a corresponding integer-lattice map.  So the integer trace collapse does NOT directly
telescope sum_{n in Lambda} K(n).  We TEST this honestly:

  (1) confirm sum_{|n|<=L} K(n) is real (imag ~ 0) on actual relation lattices -- the
      pairing telescoping is REAL and usable.
  (2) measure: per RESIDUE CLASS c, the lattice-weight sum  W_c(L) = sum_{n in Lambda, n=c} w(n).
      Then correction = sum_c W_c(L) Re? no: = sum_c (sum_{n=c} w(n)) D7(c) = sum_c W_c D7(c).
      Since W_c is REAL and D7(c) cyclotomic, the correction = sum_c W_c D7(c).
      If the W_c were EQUAL across a dilation orbit, we'd get sum over orbits of W * Tr(D7).
      Measure how close W_{ac} are to each other within an orbit (the obstruction).
  (3) The honest bound: |correction| <= (max_c |D7(c)|) * sum_c |W_c|.  But sum|W_c| ~ A(L)
      diverges. The REAL gain must come from D7-cancellation across residues sharing weight.
      Quantify: does grouping by dilation orbit REDUCE the partial-sum oscillation envelope?
"""
import sys, itertools, cmath, math
from collections import defaultdict
from fractions import Fraction
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

Z=cmath.exp(2j*math.pi/7)
QR={1,2,4}; NQR={3,5,6}
# numeric D7
Tlist=[T for r in range(7) for T in itertools.combinations(range(1,7),r)]
SGN={T:(-1)**len(T) for T in Tlist}
SIG={(T,m):sum(Z**((-m*t)%7) for t in T) for T in Tlist for m in range(1,7)}
PREF={m:(1-Z**((-m)%7)) for m in range(1,7)}
_Dcache={}
def D7(c):
    v=_Dcache.get(c)
    if v is not None: return v
    pref=1.0+0.0j
    for cj in c: pref*=PREF[cj]
    acc=0.0+0.0j
    for T in Tlist:
        p=1.0+0.0j
        for cj in c: p*=SIG[(T,cj)]
        acc+=SGN[T]*p
    v=pref*acc; _Dcache[c]=v; return v

def relations_support6(E,L):
    k=len(E); out=[]
    for idxs in itertools.combinations(range(k),6):
        es=[E[i] for i in idxs]
        dep=max(range(6),key=lambda i:abs(es[i])); e_dep=es[dep]
        if e_dep==0: continue
        free=[i for i in range(6) if i!=dep]; efree=[es[i] for i in free]
        for vfree in itertools.product(range(-L,L+1),repeat=5):
            if any(c==0 or c%7==0 for c in vfree): continue
            s=sum(c*e for c,e in zip(vfree,efree))
            if s%e_dep!=0: continue
            vd=-s//e_dep
            if vd==0 or vd%7==0 or abs(vd)>L: continue
            combo=[0]*6
            for i,c in zip(free,vfree): combo[i]=c
            combo[dep]=vd
            out.append(tuple(combo))
    return out

def w_real(vals):
    """w(n) = real value of prod 1/(2 pi i n_j) for support 6 (i^6=-1)."""
    p=1.0
    for v in vals: p*= 1.0/v
    return -p/((2*math.pi)**6)  # i^6=-1

if __name__=="__main__":
    print("LRC(14) does QR/trace collapse act on the REAL lattice sum? (kps)")
    for E,lab in [(list(range(8)),"AP8 consec"), ([0,1,3,5,7,9,11,13],"WIDE8 odd-AP")]:
        print(f"\n=== {lab}  E={E} ===")
        for L in [4,6]:
            rels=relations_support6(E,L)
            corr=0.0+0.0j
            Wc=defaultdict(float)   # weight sum per residue class
            for n in rels:
                c=tuple(v%7 for v in n)
                wv=w_real(n)
                corr += wv*D7(c)
                Wc[c]+=wv
            # imag check
            print(f"  L={L}: #rel={len(rels)}  correction Re={corr.real:+.6e} Im={corr.imag:+.2e}  (Im/Re={abs(corr.imag)/max(abs(corr.real),1e-30):.1e})")
            # within-orbit weight spread: for each dilation orbit, how equal are W_{ac}?
            orbits=defaultdict(list)
            for c in Wc:
                key=min(tuple((a*cj)%7 for cj in c) for a in range(1,7))
                orbits[key].append((c,Wc[c]))
            spreads=[]
            orbit_contrib_via_trace=0.0
            orbit_contrib_exact=0.0
            for key,members in orbits.items():
                ws=[w for _,w in members]
                if len(ws)>1 and max(abs(x) for x in ws)>0:
                    spreads.append((max(ws)-min(ws))/max(abs(max(ws)),abs(min(ws)),1e-30))
                # exact contribution of this orbit:
                oc=sum(w*D7(c) for c,w in members)
                orbit_contrib_exact+=oc.real
            avgspread=sum(spreads)/len(spreads) if spreads else float('nan')
            print(f"        #residue classes hit={len(Wc)}  #dilation-orbits={len(orbits)}  mean within-orbit weight spread={avgspread:.3f}")
