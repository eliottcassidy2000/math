#!/usr/bin/env python3
"""
paley_betti_survey.py — opus-2026-03-13-S71

Compute Paley tournament Betti numbers for p=3,7,11,19,23.
Goal: test chi(P_p)=p conjecture and find Betti pattern.

Known results:
  P_3:  β=[1,1,0], chi=0 (exception: p=3)
  P_7:  β=[1,0,0,0,6,0,0], chi=7
  P_11: β=[1,0,0,0,0,5,15,0,0,0,0], chi=11
  P_19: β=?, chi=? (Ω partial through m=8)
  P_23: β=?, chi=?

Pattern so far:
  P_7:  nonzero at m=(p-3)/2=2... no, at m=4=(p-3)/2+2
  P_11: nonzero at m=5,6
  Hypothesis: β_m=0 for m < (p-3)/2 (except β_0=1)
"""

import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

from circulant_homology import PaleyHomology
import time

# First, compute Omega dims for P_19 (partial — may take a while)
print("="*70)
print("PALEY TOURNAMENT BETTI SURVEY")
print("="*70)

# Use cached values for small primes
for p in [3, 7, 11]:
    h = PaleyHomology(p=p)
    betti = h.betti_numbers()
    omega = h.omega_dims()
    chi = sum((-1)**m * b for m, b in enumerate(betti))
    chi_omega = sum((-1)**m * o for m, o in enumerate(omega))
    print(f"\nP_{p}:")
    print(f"  Ω/p = {[o//p if o%p==0 else o for o in omega]}")
    print(f"  β   = {betti}")
    print(f"  chi(β) = {chi}")
    print(f"  chi(Ω) = {chi_omega}")

    # Find nonzero Betti positions
    nonzero = [(m, b) for m, b in enumerate(betti) if b > 0 and m > 0]
    print(f"  Nonzero β (m>0): {nonzero}")
    if nonzero:
        positions = [m for m, b in nonzero]
        print(f"  First nonzero m = {positions[0]}, (p-3)/2 = {(p-3)//2}")

# Now try P_19 — just Omega dims first
print(f"\n{'='*70}")
print("P_19: OMEGA DIMS (recomputing to verify)")
print("="*70)

h19 = PaleyHomology(p=19)
t0 = time.time()
# Start with just omega dims up to degree 8 (known to be feasible)
omega19 = h19.omega_dims(max_degree=8, use_cache=False, verbose=True)
t1 = time.time()
print(f"\nP_19 Ω (m=0..8): {omega19}")
print(f"P_19 Ω/19: {[o//19 if o%19==0 else o for o in omega19]}")
chi_partial = sum((-1)**m * o for m, o in enumerate(omega19))
print(f"P_19 chi (partial, m=0..8): {chi_partial}")
print(f"Time: {t1-t0:.1f}s")

# Try extending to m=9
print(f"\nExtending to m=9...")
try:
    t0 = time.time()
    omega19_9 = h19.omega_dims(max_degree=9, use_cache=False, verbose=True)
    t1 = time.time()
    print(f"P_19 Ω (m=0..9): {omega19_9}")
    chi_9 = sum((-1)**m * o for m, o in enumerate(omega19_9))
    print(f"P_19 chi (partial, m=0..9): {chi_9}")
    print(f"Time: {t1-t0:.1f}s")
except Exception as e:
    print(f"Failed at m=9: {e}")

print("\nDONE.")
