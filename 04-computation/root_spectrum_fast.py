#!/usr/bin/env python3
"""
root_spectrum_fast.py — oracle-2026-05-17-S1

Root spectrum: root gap, ULC, Vieta, SC asymmetry, forbidden cases.
Runs in <60 seconds.
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))

from itertools import permutations
from math import comb
from collections import defaultdict
import numpy as np
import random, time

from tournament_lib import (find_odd_cycles, hamiltonian_path_count,
                             tournament_from_bits, random_tournament)

def ip_coeffs(cycles, n):
    m = len(cycles)
    if m == 0: return [1]
    vsets = [frozenset(c) for c in cycles]
    adj_bits = [0]*m
    for a in range(m):
        for b in range(a+1,m):
            if vsets[a]&vsets[b]: adj_bits[a]|=1<<b; adj_bits[b]|=1<<a
    max_d = n//3
    coeffs = [0]*(max_d+2); coeffs[0]=1; coeffs[1]=m
    pairs = [(a,b) for a in range(m) for b in range(a+1,m) if not(adj_bits[a]>>b&1)]
    coeffs[2] = len(pairs)
    if max_d>=3:
        trips = [(a,b,c) for a,b in pairs for c in range(b+1,m)
                 if not(adj_bits[a]>>c&1) and not(adj_bits[b]>>c&1)]
        coeffs[3] = len(trips)
        if max_d>=4:
            for a,b,c in trips:
                for d in range(c+1,m):
                    if not(adj_bits[a]>>d&1) and not(adj_bits[b]>>d&1) and not(adj_bits[c]>>d&1):
                        coeffs[4]+=1
    while len(coeffs)>1 and coeffs[-1]==0: coeffs.pop()
    return coeffs

def rhos(coeffs):
    if len(coeffs)<=1: return []
    rs = np.roots(list(reversed(coeffs)))
    return sorted([-r.real for r in rs if abs(r.imag)<1e-7 and r.real<-1e-10], reverse=True)

def is_sc(T, n):
    comp = [[1-T[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]
    for p in permutations(range(n)):
        if all(T[p[i]][p[j]]==comp[i][j] for i in range(n) for j in range(n) if i!=j):
            return True
    return False

# ────────────────────────────────────────────────────────────
# Part 1: n=6 exhaustive (15s)
# ────────────────────────────────────────────────────────────

def n6_core():
    n=6; t0=time.time()
    print("="*70)
    print("N=6 EXHAUSTIVE ROOT SPECTRUM")
    print("="*70)

    iso_map = {}
    for bits in range(2**(n*(n-1)//2)):
        T = tournament_from_bits(n, bits)
        H = hamiltonian_path_count(T)
        sc_seq = tuple(sorted(sum(T[i]) for i in range(n)))
        cycles = find_odd_cycles(T)
        coeffs = ip_coeffs(cycles, n)
        I6 = sum(coeffs[k]*6**k for k in range(len(coeffs)))
        d3 = sum(1 for c in cycles if len(c)==3)
        d5 = sum(1 for c in cycles if len(c)==5)
        key = (H, I6, sc_seq, d3, d5)
        if key not in iso_map:
            iso_map[key] = dict(H=H, I6=I6, coeffs=coeffs, sc=sc_seq,
                                a1=coeffs[1] if len(coeffs)>1 else 0,
                                a2=coeffs[2] if len(coeffs)>2 else 0,
                                dc3=d3, dc5=d5, repr=T)

    print(f"  Enumeration: {time.time()-t0:.1f}s → {len(iso_map)} distinct classes (true: 56)")
    print(f"  Key = (H, I6, scores, dc3, dc5); {56-len(iso_map)} classes still merged")

    # SC checks for representatives
    for d in iso_map.values():
        d['is_sc'] = is_sc(d['repr'], n)

    classes = sorted(iso_map.values(), key=lambda d:(d['H'],d['I6'],d['dc3']))

    print(f"\n{'H':>5}{'I6':>6}{'a1':>4}{'a2':>4}{'d3':>4}{'d5':>4}{'SC':>3}",
          f"{'ρ₁':>10}{'ρ₂':>10}{'ratio':>9}", "ULC")
    print("-"*70)

    sep = defaultdict(list)
    sc_r, ns_r = [], []
    ulc_fail = 0; gap_fail = 0

    for d in classes:
        H,I6,a1,a2,deg = d['H'],d['I6'],d['a1'],d['a2'],len(d['coeffs'])-1
        co = d['coeffs']
        rv = rhos(co)
        if len(rv)==2:
            r1,r2=rv; rat=r2/r1 if r1>1e-10 else 0
            r1s,r2s,rs = f"{r1:.5f}",f"{r2:.5f}",f"{rat:.5f}"
            (sc_r if d['is_sc'] else ns_r).append((rat,H,a1,a2))
            for r in rv:
                if 0.25<r<1/3: gap_fail+=1
        elif len(rv)==1:
            r1=rv[0]; r1s,r2s,rs=f"{r1:.5f}","  =ρ₁"," 1.0000"
            if 0.25<r1<1/3: gap_fail+=1
        else:
            r1s,r2s,rs="  —","","  —"

        ulc=True
        for k in range(1,deg):
            ck,ckm,ckp=comb(deg,k),comb(deg,k-1),comb(deg,k+1)
            if ck and ckm and ckp:
                if (co[k]/ck)**2 < (co[k-1]/ckm)*(co[k+1]/ckp)-1e-9:
                    ulc=False; ulc_fail+=1
        sc_s="SC" if d['is_sc'] else "  "
        print(f"{H:>5}{I6:>6}{a1:>4}{a2:>4}{d['dc3']:>4}{d['dc5']:>4}{sc_s:>3}",
              f"{r1s:>10}{r2s:>10}{rs:>9}", "✓" if ulc else "✗")
        sep[(H,I6)].append((a1,a2,d['is_sc']))

    # Summary
    cols = {k:v for k,v in sep.items() if len(v)>1}
    unique_by_H_I6 = len(sep) - sum(len(v)-1 for v in cols.values())
    print(f"\n{'─'*70}")
    print(f"(H,I6) separation: {unique_by_H_I6}/{len(classes)} unique  (→ does NOT separate all classes)")
    print(f"Root gap (-1/3,-1/4): {gap_fail} violations ✓" if gap_fail==0 else f"ROOT GAP: {gap_fail} violations!")
    print(f"Ultra-log-concavity: {'all OK ✓' if ulc_fail==0 else f'{ulc_fail} violations!'}")
    a13 = any(d['a1']==3 and d['a2']==0 for d in classes)
    print(f"Forbidden (a1=3,a2=0): {'NOT found ✓' if not a13 else 'FOUND!'}")

    if sc_r and ns_r:
        sc_min = min(r for r,*_ in sc_r)
        ns_min = min(r for r,*_ in ns_r)
        print(f"\nRoot ratio ρ₂/ρ₁ (degree-2 only):")
        print(f"  SC ({len(sc_r)} classes): min={sc_min:.5f}")
        print(f"  NS ({len(ns_r)} classes): min={ns_min:.5f}")
        mc = min(sc_r, key=lambda x:x[0])
        print(f"  Most asymmetric SC: H={mc[1]}, a1={mc[2]}, a2={mc[3]}")
        print(f"  SC has smaller min: {'✓' if sc_min<ns_min else '✗'}")

    # Vieta
    sv,pv = [],[]
    for d in classes:
        co=d['coeffs']
        if len(co)>=3 and co[2]>0:
            rv=rhos(co)
            if len(rv)==2:
                r1,r2=rv
                sv.append(abs(r1+r2-co[1]/co[2]))
                pv.append(abs(r1*r2-1.0/co[2]))
    if sv:
        print(f"\nVieta: ρ₁+ρ₂=α₁/α₂ err={max(sv):.1e}, ρ₁ρ₂=1/α₂ err={max(pv):.1e}",
              "✓" if max(sv)<1e-6 else "✗")

    return classes

# ────────────────────────────────────────────────────────────
# Part 2: n=5 Vieta (degree-1 formula)
# ────────────────────────────────────────────────────────────

def n5_vieta():
    print("\n"+"="*70)
    print("N=5 VIETA: r = -1/α₁ = -2/(H-1)")
    n=5; errs=[]
    for _ in range(3000):
        T=random_tournament(n)
        H=hamiltonian_path_count(T)
        cycles=find_odd_cycles(T)
        co=ip_coeffs(cycles,n)
        if len(co)>=2 and co[1]>0 and H>1:
            errs.append(abs(-1/co[1]+2/(H-1)))
    print(f"  {len(errs)} cases, max|err|={max(errs):.1e} {'✓' if max(errs)<1e-9 else '✗'}")

# ────────────────────────────────────────────────────────────
# Part 3: n=7,8,9 gap + ULC + forbidden (no SC check)
# ────────────────────────────────────────────────────────────

def large_n(n, k):
    print(f"\n{'─'*70}\nN={n}: {k} samples")
    gap,a13,ulc_f,min_r=0,0,0,float('inf')
    deg_seen = defaultdict(int)
    for _ in range(k):
        T=random_tournament(n)
        cycles=find_odd_cycles(T)
        co=ip_coeffs(cycles,n)
        a1=co[1] if len(co)>1 else 0
        a2=co[2] if len(co)>2 else 0
        if a1==3 and a2==0: a13+=1
        rv=rhos(co)
        for r in rv:
            if 0.25<r<1/3: gap+=1
        if len(rv)==2: min_r=min(min_r,rv[1]/rv[0])
        d=len(co)-1; deg_seen[d]+=1
        for ki in range(1,d):
            ck,ckm,ckp=comb(d,ki),comb(d,ki-1),comb(d,ki+1)
            if ck and ckm and ckp:
                if (co[ki]/ck)**2<(co[ki-1]/ckm)*(co[ki+1]/ckp)-1e-9: ulc_f+=1
    deg_str=" ".join(f"deg{d}:{c}" for d,c in sorted(deg_seen.items()))
    print(f"  {deg_str}")
    print(f"  Gap: {gap}✓, (a1=3,a2=0): {a13}✓, ULC: {ulc_f}✓, min ratio: {min_r:.5f}" if gap==0 and a13==0 and ulc_f==0
          else f"  gap={gap}, a13_a20={a13}, ulc={ulc_f}, min_ratio={min_r:.5f}")

if __name__=='__main__':
    random.seed(42); np.random.seed(42)
    t0=time.time()

    classes = n6_core()
    n5_vieta()
    large_n(7, 2000)
    large_n(8, 300)
    large_n(9, 50)

    print(f"\n[Total: {time.time()-t0:.1f}s]")
