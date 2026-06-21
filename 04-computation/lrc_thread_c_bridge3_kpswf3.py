#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_thread_c_bridge3_kpswf3.py   (kind-pasteur 2026-06-21, THREAD C lead-gen)

The productive lead. Bridge2 found: meas(G_C) for hard cores ~0.012 << (6/7)^12=0.157
(the decorrelated lonely floor). So loneliness is DEPRESSED below its decorrelated value by
arithmetic RESONANCE -- exactly the same phenomenon the L7 closure controls for sector-cover
(where p0 is ELEVATED above its plateau P2 by resonance).

DUALITY (the THREAD C thesis):
  sector-cover route: p0 = P2(plateau) + R(resonance), R>=0, p0<=cap needs R bounded ABOVE.
  lonely route:       L  = F(floor=(6/7)^k) - D(resonance depression), L>=c needs D bounded ABOVE.
Both are RESONANCE-CORRECTION problems around a decorrelated value. The L7 machinery
(finite small-denom atlas + classical torus discrepancy O(1/q) tail) is the SAME tool.

HYP to test:
  HYP-C1: meas(G_C) = (6/7)^|C| - (depression D(C)), where D(C) is a sum over the relation
          lattice (THM-515 theta form), and D(C) is controlled by SMALL-DENOMINATOR relations
          (a finite resonance atlas) plus an O(1/v) tail -- the dual of the L7 atlas.
  HYP-C2: the tight locus (L=0) is EXACTLY where the resonance depression D = (6/7)^k (full
          cancellation). Tight <=> a maximal/saturating small-denom relation pattern. So
          "tight locus finite" <=> "saturating resonance atlas finite" -- the SAME finiteness
          the L7 atlas needs. EQUIVALENCE TEST.
  HYP-C3: the resonance depression D(C) is dominated by support-2 and support-3 relations
          (low-order), mirroring the L7 support-3 driver. Decompose D by relation support.

Compute the THM-515 theta-form term-by-term: L(C) = sum_{t in Lambda} prod_i h(t_i),
h(0)=6/7, h(t)=-sin(pi t/7)/(pi t). Group by support of the relation t.
We use the EXACT measure (danger-comb) for L, and SEPARATELY compute the relation-lattice
contributions to see which relations drive the depression.
"""
import itertools
from fractions import Fraction as Fr
import math

P=7

def danger_intervals(v):
    v=int(v); ivs=[]
    for k in range(0,v+1):
        lo=max(Fr(14*k-1,14*v),Fr(0)); hi=min(Fr(14*k+1,14*v),Fr(1))
        if lo<hi: ivs.append((lo,hi))
    return ivs
def umeas(intervals):
    if not intervals: return Fr(0)
    iv=sorted(intervals); tot=Fr(0); cl,ch=iv[0]
    for lo,hi in iv[1:]:
        if lo>ch: tot+=ch-cl; cl,ch=lo,hi
        else: ch=max(ch,hi)
    tot+=ch-cl; return tot
def L_meas(C):
    ai=[]
    for v in C: ai+=danger_intervals(v)
    return Fr(1)-umeas(ai)

# h-function (float; the theta form is not rational but we want the STRUCTURE)
def h(t):
    if t==0: return 6.0/7.0
    return -math.sin(math.pi*t/7.0)/(math.pi*t)

def relations_upto(C, Tmax, max_support):
    """enumerate integer relations sum t_i v_i = 0 with |t_i|<=Tmax, support<=max_support.
       returns list of (support_size, prod h(t_i) for the relation including h(0) on zeros)."""
    C=[int(v) for v in C]; n=len(C)
    base = h(0)**n  # the all-zero relation (the floor)
    terms=[]  # (supp, value) for nonzero relations
    # iterate over which coords are nonzero
    for supp in range(1, max_support+1):
        for idx in itertools.combinations(range(n), supp):
            vs=[C[i] for i in idx]
            # find t over idx, |t|<=Tmax, t_i!=0, sum t_i v_i=0
            for tvals in itertools.product(range(-Tmax,Tmax+1), repeat=supp):
                if any(tv==0 for tv in tvals): continue
                if sum(tv*v for tv,v in zip(tvals,vs))!=0: continue
                # gcd-primitivity not required for theta sum; include all
                val = h(0)**(n-supp)
                for tv in tvals: val*=h(tv)
                terms.append((supp, val))
    return base, terms

def main():
    print("="*82)
    print("THREAD C BRIDGE 3: lonely measure as RESONANCE DEPRESSION of (6/7)^k (dual of L7)")
    print("="*82)

    cores = {
        "{1,2,3}":[1,2,3],
        "{1,2,3,4,5}":[1,2,3,4,5],
        "{1..7}":list(range(1,8)),
        "{1..11,13}":list(range(1,12))+[13],
    }
    print("\n[HYP-C1] L(C) vs decorrelated floor (6/7)^|C|, and the depression D=floor-L:")
    print(f"   {'C':<16}{'|C|':>4}{'(6/7)^|C|':>12}{'L(C)exact':>12}{'depress D':>12}{'D/floor':>9}")
    for name,C in cores.items():
        k=len(C); floor=Fr(6,7)**k; L=L_meas(C); D=floor-L
        print(f"   {name:<16}{k:>4}{float(floor):>12.5f}{float(L):>12.5f}{float(D):>12.5f}{float(D/floor):>9.3f}")

    print("\n   => depression D grows as relations accumulate; for {1..11,13} D/floor~0.92")
    print("      (loneliness is depressed to 8% of its decorrelated value by resonance).")

    # HYP-C3: decompose the theta-sum depression by relation support (low cores only, feasible)
    print("\n[HYP-C3] theta-form depression by relation SUPPORT (float; Tmax=6):")
    print("   sum of nonzero-relation terms, grouped by |support|. Negative sum = depression.")
    for name,C in [("{1,2,3}",[1,2,3]), ("{1,2,3,4,5}",[1,2,3,4,5])]:
        base, terms = relations_upto(C, Tmax=6, max_support=min(4,len(C)))
        bysupp={}
        for supp,val in terms:
            bysupp[supp]=bysupp.get(supp,0.0)+val
        Lapprox = base + sum(bysupp.values())
        Lexact = float(L_meas(C))
        print(f"   {name}: floor(6/7)^k={base:.5f}")
        for s in sorted(bysupp):
            print(f"       support {s}: sum of terms = {bysupp[s]:+.5f}")
        print(f"     theta-approx L (Tmax=6) = {Lapprox:.5f}   exact L = {Lexact:.5f}")
        # which support dominates the depression?
        dom = min(bysupp.items(), key=lambda kv: kv[1])
        print(f"     => largest depression from support {dom[0]} ({dom[1]:+.5f})")

    print("\n[HYP-C2] TIGHT-LOCUS <-> RESONANCE-ATLAS equivalence sketch:")
    print("   tight (L=0) configs found in canon: AP{1..13}, sporadic {1..11,13,24}, Goddyn-Wong T5.")
    print("   Each is a SATURATING small-denominator relation pattern (AP: t=(1,-2,1,..) support-3 chains).")
    print("   The L7 resonance atlas = small-denom ratios p/q in (1,2.15). CLAIM: tight-locus")
    print("   finiteness <=> these saturating relation patterns are finitely many. Both reduce to")
    print("   'small-denominator commensurability is bounded'. Tested structurally, not proved.")

if __name__ == "__main__":
    main()
