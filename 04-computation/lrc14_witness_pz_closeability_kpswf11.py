#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_witness_pz_closeability_kpswf11.py  (kind-pasteur 2026-06-22, THREAD 3)

FOLLOW-UP to the Paley-Zygmund probe: is the PZ floor mu_{1/7} >= (EN)^2/E[N^2]
CLOSEABLE with the existing sector-moment tools?

PZ needs:  (a) a LOWER bound on E N  (first sector-emptiness moment), and
           (b) an UPPER bound on E[N^2] (= E N + 2 * sum_{j<j'} P(j,j' both empty)).

KEY STRUCTURAL FACTS to test:
 (A) E N = sum_{j=1..6} P(sector j empty).  P(sector j empty) = meas{x: no e_i x
     in [j/7,(j+1)/7)}.  For the DECORRELATED (iid) model each sector is empty w.p.
     (6/7)^(k-1) (k-1 nonzero speeds, 0 sits in sector 0).  So the decorrelated
     E N_dec = 6*(6/7)^(k-1).  Compare to the true E N: is true E N >= E N_dec?
     (If the single-sector-empty prob is >= its decorrelated value, E N has a clean
      closed lower bound -- the FIRST MOMENT is the EASY, decorrelated quantity.)
 (B) Cross term P(j, j' both empty): is it close to / below the decorrelated
     (6/7)^... no -- two disjoint sectors empty for iid is ((5/7))^(k-1) only if the
     k-1 points each avoid BOTH sectors = (5/7)^(k-1).  So E[N^2]_dec = EN_dec +
     2*C(6,2)*(5/7)^(k-1).  Test whether true E[N^2] <= some clean multiple.
 (C) The DECORRELATED Paley-Zygmund value PZ_dec = (EN_dec)^2/E[N^2]_dec -- is it
     itself >= m_P with room?  If yes, and if the true moments are controlled
     relative to decorrelated by the SAME Delsarte/coupon machinery (THM-534/HYP-2708
     death-chain), the floor closes through already-proved tools.

We compute everything EXACTLY and compare true vs decorrelated, k=8..13.
"""
import itertools
from fractions import Fraction as Fr
from math import gcd

m_P = Fr(14249, 252252)

def sector_breakpoints(E):
    pts = set()
    for e in E:
        if e == 0: continue
        for t in range(0, 7*e + 1):
            pts.add(Fr(t, 7*e))
    pts.add(Fr(0)); pts.add(Fr(1))
    return sorted(p for p in pts if Fr(0) <= p <= Fr(1))

def empty_inner(E, x):
    hit = set()
    for e in E:
        j = int((Fr(e)*x % 1) * 7)
        if j == 7: j = 6
        hit.add(j)
    return set(range(1,7)) - hit

def moments(E):
    bps = sector_breakpoints(E)
    EN = Fr(0); EN2 = Fr(0)
    single = {a: Fr(0) for a in range(1,7)}
    pair_off = Fr(0)   # sum over ORDERED j!=j' of P(both empty)
    for a, b in zip(bps, bps[1:]):
        if b <= a: continue
        mid = (a+b)/2; L = b-a
        empt = empty_inner(E, mid); n = len(empt)
        EN += n*L; EN2 += n*n*L
        for s in empt: single[s] += L
        pair_off += n*(n-1)*L
    return EN, EN2, single, pair_off

def decorrelated(k):
    """iid model: k-1 nonzero speeds, each lands uniformly; sector 0 holds e=0."""
    m = k-1
    EN_dec = 6*Fr(6,7)**m
    # two distinct inner sectors both empty: each of m points avoids both -> (5/7)^m
    pair_dec = 6*5*Fr(5,7)**m   # ordered pairs
    EN2_dec = EN_dec + pair_dec
    PZ_dec = (EN_dec*EN_dec)/EN2_dec
    return EN_dec, EN2_dec, PZ_dec

def consec(k): return list(range(k))

def fr(x): return f"{float(x):.6f}"

print("="*92)
print("LRC(14) WITNESS FLOOR -- PALEY-ZYGMUND CLOSEABILITY  [kpswf11, THREAD 3]")
print("="*92)
print(f"m_P = {m_P} = {fr(m_P)}")
print()
print("Comparing TRUE sector-emptiness moments vs DECORRELATED (iid coupon) model.")
print("Decorrelated E N = 6(6/7)^(k-1);  E[N^2]_dec = EN_dec + 30(5/7)^(k-1).")
print()
print(f"{'k':>3} | {'EN_true':>9} {'EN_dec':>9} {'EN>=dec?':>8} | "
      f"{'EN2_true':>9} {'EN2_dec':>9} {'EN2<=dec?':>9} | "
      f"{'PZ_true':>8} {'PZ_dec':>8} {'PZdec>=mP?':>10}")
for k in range(8, 14):
    E = consec(k)
    EN, EN2, single, pair_off = moments(E)
    PZ = (EN*EN)/EN2
    EN_d, EN2_d, PZ_d = decorrelated(k)
    print(f"{k:>3} | {fr(EN):>9} {fr(EN_d):>9} {str(EN>=EN_d):>8} | "
          f"{fr(EN2):>9} {fr(EN2_d):>9} {str(EN2<=EN2_d):>9} | "
          f"{fr(PZ):>8} {fr(PZ_d):>8} {str(PZ_d>=m_P):>10}")

print()
print("-"*92)
print("[B] Is single-sector-empty prob >= decorrelated (6/7)^(k-1) for ALL E?")
print("    (would give a clean closed-form LOWER bound on E N)")
print("-"*92)
print(f"{'k':>3} {'#E':>6} {'min P(sec empty)':>16} {'(6/7)^(k-1)':>12} {'always>=dec?':>12}")
for k in range(8, 13):
    bank = []
    for rest in itertools.combinations(range(1,15), k-1):
        Etup = (0,)+rest
        g = 0
        for e in Etup: g = gcd(g,e)
        if g==1: bank.append(Etup)
    dec = Fr(6,7)**(k-1)
    worst = None; allok = True
    for Etup in bank:
        EN, EN2, single, pair_off = moments(list(Etup))
        mn = min(single.values())
        if worst is None or mn < worst: worst = mn
        if mn < dec: allok = False
    print(f"{k:>3} {len(bank):>6} {fr(worst):>16} {fr(dec):>12} {str(allok):>12}")

print()
print("-"*92)
print("[C] Decorrelated PZ floor across k -- is it a clean uniform constant >> m_P?")
print("-"*92)
for k in range(8, 25):
    EN_d, EN2_d, PZ_d = decorrelated(k)
    print(f"  k={k:>2}: EN_dec={fr(EN_d):>9}  PZ_dec={fr(PZ_d):>9}  PZ_dec/m_P={float(PZ_d/m_P):>6.2f}x  PZ_dec>=m_P:{PZ_d>=m_P}")

print()
print("-"*92)
print("[D] The honest gap: does TRUE PZ >= DECORRELATED PZ (so dec is a valid floor)?")
print("    or does correlation push true PZ BELOW dec? (tested on consec + wide)")
print("-"*92)
for k in range(8,14):
    E = consec(k); EN,EN2,_,_ = moments(E); PZ=(EN*EN)/EN2
    _,_,PZ_d = decorrelated(k)
    tag = "true>=dec (dec is floor)" if PZ>=PZ_d else "true<dec (corr LOWERS PZ!)"
    print(f"  k={k:>2} consec: PZ_true={fr(PZ)}  PZ_dec={fr(PZ_d)}  {tag}")
# wide
for k in (9,10,11):
    base=list(range(k-1))
    for f in (15,40,100):
        E=base+[f]; EN,EN2,_,_=moments(E); PZ=(EN*EN)/EN2; _,_,PZ_d=decorrelated(k)
        tag="true>=dec" if PZ>=PZ_d else "true<dec"
        print(f"  k={k:>2} f={f:>3}: PZ_true={fr(PZ)}  PZ_dec={fr(PZ_d)}  {tag}")
print("="*92)
