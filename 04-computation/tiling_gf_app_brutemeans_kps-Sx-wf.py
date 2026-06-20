#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Optimized BRUTE means over all tilings (exact Fraction), flushed per n.
Confirms per-subset/closed-form means against true enumeration:
  E[c3], E[c5], E[alpha2-tri] for n<=7 ; full E[alpha2] at n=8.
Adjacency stored as bitmask rows (out[i] = bitset of j with i->j) for speed.
"""
import sys, time
from itertools import product, combinations, permutations
from fractions import Fraction
from math import comb
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

def tiles(n):
    return [(a,b) for a in range(3,n+1) for b in range(1,a-1)]

# Precompute, for each subset (as tuple), the list of (ordered cyclic perms) to test.
# We test cycles by: fix smallest vertex, iterate permutations of the rest.
def cyclic_orderings(sub):
    sub=tuple(sorted(sub)); s0=sub[0]; rest=sub[1:]
    return [ (s0,)+p for p in permutations(rest) ]

def count_dcycles(seqs, out):
    c=0
    for seq in seqs:
        ok=True
        k=len(seq)
        for i in range(k):
            u=seq[i]; w=seq[(i+1)%k]
            if not (out[u]>>w)&1: ok=False;break
        if ok: c+=1
    return c

def run(n, do_c5=True, do_full_a2=False):
    T=tiles(n); F=len(T); N=1<<F
    base_out=[0]*(n+1)
    for k in range(n,1,-1):
        base_out[k]|=(1<<(k-1))
    trips=list(combinations(range(1,n+1),3))
    tri_seqs={t:cyclic_orderings(t) for t in trips}
    fives=list(combinations(range(1,n+1),5)) if do_c5 else []
    five_seqs={t:cyclic_orderings(t) for t in fives}
    # disjoint triple pairs
    tri_pairs=[(i,j) for i in range(len(trips)) for j in range(i+1,len(trips))
               if set(trips[i]).isdisjoint(trips[j])]
    # for full alpha2 at n=8 also need 3-5 disjoint pairs (5 cycles vs 3 cycles)
    odd_cycles=[]  # for full a2: list of (idx into a combined structure)
    sum_c3=0; sum_c5=0; sum_a2tri=0; sum_a2full=0
    t0=time.time()
    for bv in product((0,1),repeat=F):
        out=base_out[:]
        for (a,b),bit in zip(T,bv):
            if bit==0: out[a]|=(1<<b)
            else: out[b]|=(1<<a)
        # c3 via scores (fast) ; also need per-triple cycle indicator for a2tri
        tri_ind={}
        c3=0
        for t in trips:
            cc=count_dcycles(tri_seqs[t],out)
            tri_ind[t]=cc; c3+=cc
        sum_c3+=c3
        # a2 triangle-pair
        a2tri=0
        for (i,j) in tri_pairs:
            a2tri += tri_ind[trips[i]]*tri_ind[trips[j]]
        sum_a2tri+=a2tri
        if do_c5 or do_full_a2:
            five_ind={}
            c5=0
            for t in fives:
                cc=count_dcycles(five_seqs[t],out)
                five_ind[t]=cc; c5+=cc
            sum_c5+=c5
            if do_full_a2:
                # full a2 = disjoint pairs among all odd cycles (3 and 5)
                # 3-3 disjoint (a2tri) + 3-5 disjoint + 5-5 disjoint (none fit in n=8)
                a2full=a2tri
                for t3 in trips:
                    if tri_ind[t3]==0: continue
                    s3=set(t3)
                    for t5 in fives:
                        if five_ind[t5]==0: continue
                        if s3.isdisjoint(t5):
                            a2full += tri_ind[t3]*five_ind[t5]
                # 5-5 disjoint needs 10 vertices > 8, skip
                sum_a2full+=a2full
    res={"Ec3":Fraction(sum_c3,N),"Ea2tri":Fraction(sum_a2tri,N)}
    if do_c5 or do_full_a2: res["Ec5"]=Fraction(sum_c5,N)
    if do_full_a2: res["Ea2full"]=Fraction(sum_a2full,N)
    res["t"]=time.time()-t0; res["N"]=N
    return res

def Ec5_formula(n):
    n=Fraction(n)
    return (Fraction(1,160)*n**5 - Fraction(1,16)*n**4 + Fraction(9,32)*n**3
            - Fraction(7,8)*n**2 + Fraction(147,80)*n - Fraction(7,4))

CLAIM_Ec5={5:Fraction(19,16),6:Fraction(49,8),7:Fraction(315,16)}
CLAIM_Ea2={5:Fraction(0),6:Fraction(15,16),7:Fraction(93,16)}

print("BRUTE confirmation of means (exact rationals)", flush=True)
print("n |   Ec3      Ec5(formula)        Ea2tri(claim)    | t(s)", flush=True)
for n in range(5,8):
    r=run(n, do_c5=True, do_full_a2=False)
    f5=Ec5_formula(n); m5=(r["Ec5"]==f5==CLAIM_Ec5[n])
    ma=(r["Ea2tri"]==CLAIM_Ea2[n])
    ec3f=Fraction(comb(n,3)+(n-2),4); m3=(r["Ec3"]==ec3f)
    print(f"{n}: c3={r['Ec3']}({'OK' if m3 else 'X'}) c5={r['Ec5']} vs {f5} ({'OK' if m5 else 'X'})"
          f"  a2tri={r['Ea2tri']} vs {CLAIM_Ea2[n]} ({'OK' if ma else 'X'})  | {r['t']:.1f}s", flush=True)

print("\nn=8 FULL alpha2 (3-3 + 3-5 disjoint), claim 1139/32, tri claim 173/8, cross 447/32", flush=True)
r8=run(8, do_c5=True, do_full_a2=True)
print(f"n=8: Ea2tri={r8['Ea2tri']} (claim 173/8={Fraction(173,8)})", flush=True)
print(f"n=8: Ea2full={r8['Ea2full']} (claim 1139/32={Fraction(1139,32)})", flush=True)
print(f"n=8: cross=full-tri={r8['Ea2full']-r8['Ea2tri']} (claim 447/32={Fraction(447,32)})", flush=True)
print(f"n=8: Ec5={r8['Ec5']} (formula {Ec5_formula(8)}, claim 199/4={Fraction(199,4)})", flush=True)
m_tri = r8['Ea2tri']==Fraction(173,8)
m_full = r8['Ea2full']==Fraction(1139,32)
m_cross = (r8['Ea2full']-r8['Ea2tri'])==Fraction(447,32)
print(f"n=8 verdict: tri={'OK' if m_tri else 'X'} full={'OK' if m_full else 'X'} cross={'OK' if m_cross else 'X'}  ({r8['t']:.1f}s)", flush=True)
