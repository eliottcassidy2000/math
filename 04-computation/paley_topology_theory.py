#!/usr/bin/env python3
"""
TOPOLOGICAL INTERPRETATION OF PALEY PATH HOMOLOGY

For P_7: β = (1,0,0,0,6,0) → 6 independent 4-dimensional holes.
This is NOT S^4 (which would have β_4=1), nor is it 6 × S^4.
Rather, it's a space with the homology of a WEDGE of 6 copies of S^4:
  P_7 ~ ∨^6 S^4 (homologically)

For P_3: β = (1,1,0) → S^1 (circle)

KEY: The per-eigenspace decomposition gives:
  P_7: each k≠0 contributes β_4=1 → each eigenspace is like S^4
  Total is the "direct sum" of 6 spheres

This suggests a SUSPENSION structure:
  P_3 ~ S^1 (1-sphere)
  P_7 ~ ? (4-sphere × 6?)

Is there a suspension or join operation connecting these?

Let's also investigate: what do non-Paley circulant tournaments look like?
"""
import numpy as np
import sys, time
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

import os
old_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')
from path_homology_fourier_v3 import fourier_betti_v3, enumerate_step_sequences
sys.stdout = old_stdout

def qr(p):
    return set((a*a) % p for a in range(1, p))

# ===== All circulant tournaments at p=7 =====
print("=" * 70)
print("ALL CIRCULANT TOURNAMENTS AT p=7: TOPOLOGY COMPARISON")
print("=" * 70)

from itertools import combinations

p = 7
all_ct = []
for S_comb in combinations(range(1, p), (p-1)//2):
    S_set = set(S_comb)
    neg_S = {(p - s) % p for s in S_set}
    if S_set & neg_S == set():
        # Compute F = illegal merges
        F = set()
        L = set()
        for a in S_set:
            for b in S_set:
                if a != b:
                    m = (a + b) % p
                    if m != 0:
                        if m in S_set:
                            L.add(m)
                        else:
                            F.add(m)

        betti = fourier_betti_v3(S_set, p, max_dim=p-1)
        chi = sum((-1)**k * betti[k] for k in range(len(betti)))
        is_qr = (S_set == qr(p))
        label = " ★PALEY★" if is_qr else ""

        all_ct.append((S_set, betti, chi, len(F), len(L)))
        print(f"  S={sorted(S_set)}: β={betti}, χ={chi}, |F|={len(F)}, |L|={len(L)}{label}")

# ===== The Euler characteristic pattern =====
print(f"\nEuler characteristic values:")
chi_vals = set(chi for _, _, chi, _, _ in all_ct)
for chi_val in sorted(chi_vals):
    sets_with_chi = [(S, b) for S, b, c, _, _ in all_ct if c == chi_val]
    print(f"  χ={chi_val}: {len(sets_with_chi)} tournaments")
    for S, b in sets_with_chi:
        print(f"    S={sorted(S)}: β={b}")

# ===== Geometric classification =====
print(f"\nGeometric classification:")
for S, betti, chi, nF, nL in all_ct:
    nonzero = [(i, betti[i]) for i in range(len(betti)) if betti[i] > 0 and i > 0]
    if len(nonzero) == 0:
        topo = "point (contractible)"
    elif len(nonzero) == 1:
        d, b = nonzero[0]
        if b == 1:
            topo = f"S^{d} ({d}-sphere)"
        else:
            topo = f"∨^{b} S^{d} (wedge of {b} {d}-spheres)"
    else:
        parts = [f"β_{d}={b}" for d, b in nonzero]
        topo = ", ".join(parts)
    print(f"  S={sorted(S)}: {topo}")

# ===== The key question: dimension = p-3? =====
print(f"\n\n{'='*70}")
print("DIMENSION vs p FOR PALEY: WHERE DOES β CONCENTRATE?")
print("="*70)

print("""
P_3: β_1 = 1 (but from k=0, not k≠0)
P_7: β_4 = 6 (from k≠0, dim = p-3 = 4)

If d = p-3 for all Paley with p ≥ 7:
  P_11: β_8 = 10 (pred: 10 × S^8?)
  P_19: β_{16} = 18 (pred: 18 × S^{16}?)
  P_23: β_{20} = 22 (pred: 22 × S^{20}?)

This would mean P_p homologically resembles ∨^{p-1} S^{p-3}.
The Euler characteristic would be:
  χ = 1 + (-1)^{p-3} · (p-1) = 1 + (p-1) = p  (for p-3 even ↔ p ≡ 3 mod 4 ✓!)
  since p-3 ≡ 0 mod 4 when p ≡ 3 mod 4

This is CONSISTENT with χ = p!

More precisely: p-3 mod 2 = (p-3) mod 2.
  p ≡ 3 mod 4 → p-3 ≡ 0 mod 4, so p-3 is even.
  (-1)^{p-3} = 1.
  χ = 1 + 1·(p-1) = p. ✓✓✓
""")

# ===== All circulant tournaments at p=11 =====
print("=" * 70)
print("ALL CIRCULANT TOURNAMENTS AT p=11")
print("=" * 70)

p = 11
all_ct_11 = []
t0 = time.time()
for S_comb in combinations(range(1, p), (p-1)//2):
    S_set = set(S_comb)
    neg_S = {(p - s) % p for s in S_set}
    if S_set & neg_S == set():
        is_qr = (S_set == qr(p))

        # Only compute up to dim 4 for non-Paley (speed)
        if is_qr:
            max_d = 4  # Will get full P_11 later
        else:
            max_d = 4

        betti = fourier_betti_v3(S_set, p, max_dim=max_d)
        chi = sum((-1)**k * betti[k] for k in range(len(betti)))
        label = " ★PALEY★" if is_qr else ""
        all_ct_11.append((S_set, betti, chi, is_qr))
        print(f"  S={sorted(S_set)}: β={betti}, χ≈{chi}{label}", flush=True)

t1 = time.time()
print(f"\n  Total: {len(all_ct_11)} circulant tournaments ({t1-t0:.1f}s)")

# Are there other χ=p tournaments?
chi_p = [(S, b) for S, b, c, _ in all_ct_11 if c == p]
print(f"  χ≈{p}: {len(chi_p)} tournaments (through dim {max_d})")

print("\nDone.")
