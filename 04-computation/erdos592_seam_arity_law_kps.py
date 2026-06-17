#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THE SEAM-ARITY LAW (THM-521): why the Erdős-592 seam is at p=2, and how it tracks arity.
kind-pasteur-2026-06-16-S4. Resolves the orphan T778; answers the snippet handoff "explain
the v2/v3 asymmetry algebraically". Builds on THM-469 (sum-free grading ⟺ p=2).

A "k-term Schur relation" on a level L_0(p)={n≥1: p∤n} (the units): a_1+...+a_k=b, all in L_0.
The Erdős-592 TRIANGLE closure is k=2 (gaps add along a chain: g(x,z)=g(x,y)+g(y,z)).
LAW (both elementary, proved in THM-521; verified here):
  (FULL)     L_0(p) is k-term-sum-free  ⟺  p=2 AND k even.
  (DIAGONAL) the degenerate k·a=b is killed by v_p  ⟺  p|k.
Triangle: k=2 (even) ⟹ seam p=2 (and 2|2, so FULL and DIAGONAL coincide at the triangle).
CHAR-2 framing: p=2 is the unique field where the all-ones 2-term Schur sum 1+1=0 is a non-unit.
"""
import sys; sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None
from itertools import combinations_with_replacement

def vp(n,p):
    k=0
    while n%p==0: n//=p; k+=1
    return k

def full_sumfree(p,k,N=3000):
    units=[n for n in range(1,N+1) if n%p!=0]; U=set(units)
    small=[u for u in units if u<=20][:8]
    for c in combinations_with_replacement(small,k):
        s=sum(c)
        if s<=N and s%p!=0: return False
    return True

print("(FULL) L_0(p) k-term-sum-free? predict TRUE ⟺ p=2 & k even")
for k in (2,3,4,5,6):
    r=[f"p{p}:{'SF' if full_sumfree(p,k) else '-'}" for p in (2,3,5,7)]
    print(f"  k={k} {'even' if k%2==0 else 'odd':<4}: {'  '.join(r)}")
print("(DIAGONAL) k·a escapes L_v ⟺ p|k")
for k in (2,3,4,5,6,12):
    killers=[p for p in (2,3,5,7,11) if all(vp(k*a,p)!=vp(a,p) for a in range(1,60) if a%p!=0)]
    print(f"  k={k:2d}: killers={killers}  divisors-of-k={[p for p in (2,3,5,7,11) if k%p==0]}")
print("\nTRIANGLE (k=2, even): seam p=2; char-2 = unique field with 1+1=0 (all-ones 2-term sum vanishes).")
print("BARRIER COROLLARY: the ONLY fully-sum-free valuation grading for the (even-arity-2) triangle is p=2,")
print("which runs out at the linear wall t=2n+1 (HYP-2396) — so a STRONG Specker witness (Q(n,t) SAT ∀t)")
print("CANNOT come from any valuation/invariant grading; it must use non-invariant (Larson partial-sum) features.")
