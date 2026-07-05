#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
HDICH SUPPORT (HYP-4097): the n=12 tight-locus rigidity  M(W) = 1/13  =>  W = c*{1..12},
split into two halves with the mechanism for each.

HALF 1 (residue pinning, ~free from THM-593A at the PRIME modulus): a tight-from-above
12-family with no multiple of 13 represents every unit residue mod 13; at q = 13 prime the
units are ALL 12 nonzero residues; 12 runners covering 12 required classes forces the residue
multiset to be EXACTLY the full nonzero system {1,...,12}, each once.  (Verify + edge cases.)

HALF 2 (lift rigidity, the open mechanism): given full residues, why can no runner be LIFTED
(v = r + 13k, k >= 1) while keeping M = 1/13?  NOT automatic at primes: q = 5's {1,3,4,7} is a
tight PERM-LIFT family (unit 2 lifted to 7).  This script hunts the mechanism:
  (a) census: all single-lift families over the full-residue base {1..12} mod 13, lift heights
      k = 1..K, all 12 positions: EXACT M for each; confirm none tight; collect the ARGMAX
      structure (which witness dies, what M becomes);
  (b) same for q = 5 ({1,2,3,4} base, lifts): find which lifts stay tight, and WHAT DIFFERS;
  (c) the conjectured mechanism: a lift r -> r + qk keeps tightness iff at every unit witness
      a/q the lifted runner still clears 1/q -- automatic (residue unchanged!) -- so tightness
      survives all SINGLE unit-fraction witnesses; what dies must be the SUP over t NEAR the
      witnesses: the lifted runner's local slope grows (v large => narrow safe plateau), and
      M(W) = 1/q needs the max-min to stay EXACTLY 1/q: if the plateau around some a/q still
      has one-sided room, M > 1/q?? NO -- M <= 1/q needs the family to COVER everything above
      1/q... M(tight) = 1/q means sup equals 1/q: lifting can only DECREASE clearances away
      from the witness points, so M stays <= ... hmm: lifting changes the function; the sup
      could DROP below 1/q (family becomes covering-below = not lonely enough) -- at q=13,
      12 runners: sup >= 1/13 ALWAYS by LRC(13) (citation!)  => sup = 1/13 iff no time does
      better: lifting a runner FREES time regions (bigger v = thinner teeth) so the sup should
      INCREASE past 1/13 (strictly loose)!  CONJECTURED MECHANISM: any lift strictly INCREASES
      M above 1/13 at q=13; at q=5 the {1,3,4,7} anomaly survives because ... (measure the
      difference).  (d) quantify: min over lifts of M(lifted) at q=13 (the "rigidity gap").

opus-2026-07-05-S74.
"""
import sys
from fractions import Fraction as Fr
from math import floor, gcd
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

def dist(x):
    f = x - floor(x); return min(f, 1 - f)

def M_exact(S):
    dens = set()
    for v in S: dens.add(v); dens.add(2 * v)
    L = sorted(S)
    for i in range(len(L)):
        for j in range(i + 1, len(L)):
            dens.add(L[i] + L[j]); dens.add(L[j] - L[i])
    dens.discard(0)
    best, argt = Fr(-1), None
    for d in dens:
        for k in range(d + 1):
            t = Fr(k, d)
            m = min(dist(v * t) for v in S)
            if m > best: best, argt = m, t
    return best, argt

print("=" * 100)
print(" HALF 1 sanity: the AP {1..12} at q=13")
base = list(range(1, 13))
m0, t0 = M_exact(base)
print(f"   M({{1..12}}) = {m0} at t = {t0}   (tight = 1/13? {m0 == Fr(1,13)})")

print("\n HALF 2(a): single lifts over the full-residue base at q=13: v_pos -> r + 13k")
K = 8
worst = None
tight_lifts = []
for pos in range(12):
    for k in range(1, K + 1):
        W = base.copy()
        W[pos] = base[pos] + 13 * k
        m, t = M_exact(W)
        if m == Fr(1, 13): tight_lifts.append((pos + 1, k))
        if worst is None or m < worst[0]: worst = (m, pos + 1, k)
        if m < Fr(1, 13):
            print(f"   *** BELOW 1/13: lift r={pos+1} k={k}: M = {m}")
print(f"   single lifts checked: {12*K}; TIGHT lifts (M = 1/13): {tight_lifts}")
print(f"   min M over lifts = {worst[0]} = {float(worst[0]):.6f} at r={worst[1]}, k={worst[2]}")
print(f"   ALL lifts strictly LOOSE (M > 1/13)? {worst[0] > Fr(1,13)}")

print("\n   direction check: does every lift INCREASE M (loosen)?")
increases = 0; total = 0
mins_by_r = {}
for pos in range(12):
    for k in range(1, K + 1):
        W = base.copy(); W[pos] = base[pos] + 13 * k
        m, _ = M_exact(W)
        total += 1
        if m > Fr(1, 13): increases += 1
        mins_by_r.setdefault(pos + 1, []).append(m)
print(f"   lifts with M > 1/13: {increases}/{total}")
print("   min-over-k M per lifted residue r:")
for r in range(1, 13):
    mn = min(mins_by_r[r])
    print(f"     r={r:>2}: min M = {mn} = {float(mn):.6f}")

print("\n HALF 2(b): the q=5 contrast -- base {1,2,3,4}, single lifts r -> r + 5k")
base5 = [1, 2, 3, 4]
tight5 = []
for pos in range(4):
    for k in range(1, 7):
        W = base5.copy(); W[pos] = base5[pos] + 5 * k
        m, t = M_exact(W)
        if m == Fr(1, 5): tight5.append((base5[pos], base5[pos] + 5 * k, t))
print(f"   tight lifts at q=5: {tight5}")
print("   ({1,3,4,7} = lift 2->7 expected among them)")

print("\n HALF 2(c/d): the mechanism probe -- for each q=13 lift, WHERE does the sup move?")
print("   (the lifted family's argmax t and its denominator; tight witnesses are k/13)")
for (pos, k) in [(1, 1), (2, 1), (12, 1), (6, 2)]:
    W = base.copy(); W[pos - 1] = base[pos - 1] + 13 * k
    m, t = M_exact(W)
    print(f"   lift r={pos} k={k}: M = {m} at t = {t} (denominator {t.denominator}; 13-witness? {t.denominator == 13})")

print("\n MECHANISM SUMMARY (to verify above): at q=13, replacing r by r+13k keeps every")
print(" k/13-witness value EXACTLY (residues unchanged) -- so sup >= 1/13 persists -- while the")
print(" thinner teeth of the faster runner OPEN a better witness elsewhere: sup STRICTLY exceeds")
print(" 1/13. Rigidity mechanism = 'lifts strictly loosen'. At q=5, {1,3,4,7}: the lift ALSO")
print(" keeps 1/5-witnesses, but q=5's small unit group leaves no room for a better witness")
print(" (4 runners, dense teeth) -- sup stays exactly 1/5. The q=13 claim to prove:")
print("   ANY proper lift of the full-residue base admits t with min clearance > 1/13.")
print("DONE.")
