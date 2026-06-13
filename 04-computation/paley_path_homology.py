#!/usr/bin/env python3
"""
PATH HOMOLOGY OF PALEY TOURNAMENTS

The Paley tournament P_p (p prime ≡ 3 mod 4) is the circulant tournament
on Z_p with connection set S = QR(p) (quadratic residues mod p).

Properties:
- |S| = (p-1)/2, S ∩ (-S) = ∅ (since -1 is a QNR when p ≡ 3 mod 4)
- P_p is a tournament (exactly one edge between each pair)
- P_p is self-complementary: P_p ≅ P_p^op (from the QR structure)
- P_p is vertex-transitive (circulant)

Key questions:
1. What are the path Betti numbers of Paley tournaments?
2. Do they match the tournament phase pattern (only odd β_p)?
3. How does QR structure relate to β_2 = 0?
4. Does the Fourier decomposition simplify for QR connection sets?
"""
import numpy as np
from itertools import combinations
from collections import Counter
import sys
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import path_betti_numbers, circulant_digraph
from path_homology_fourier_v3 import fourier_betti_v3


def quadratic_residues(p):
    """Compute QR(p) = {a^2 mod p : a ∈ {1,...,p-1}}."""
    qr = set()
    for a in range(1, p):
        qr.add((a * a) % p)
    return qr


def is_paley_prime(p):
    """Check if p is prime and ≡ 3 mod 4 (Paley condition)."""
    if p < 3:
        return False
    if p % 4 != 3:
        return False
    for d in range(2, int(p**0.5) + 1):
        if p % d == 0:
            return False
    return True


# ===== List Paley primes =====
print("=" * 70)
print("PALEY TOURNAMENTS: PATH HOMOLOGY")
print("=" * 70)

paley_primes = [p for p in range(3, 100) if is_paley_prime(p)]
print(f"\nPaley primes (≡ 3 mod 4): {paley_primes}")

# ===== Compute path homology for small Paley tournaments =====
print("\n--- Path Homology (full computation) ---")

for p in [3, 7, 11]:
    S = quadratic_residues(p)
    S_comp = set(range(1, p)) - S

    # Verify tournament property: S and n-S partition [1,p-1]
    neg_S = {(p - s) % p for s in S}
    is_tournament = (S & neg_S == set()) and (S | neg_S == set(range(1, p)))

    print(f"\nP_{p}:")
    print(f"  QR = {sorted(S)}")
    print(f"  QNR = {sorted(S_comp)}")
    print(f"  Is tournament? {is_tournament}")

    # F = illegal merges
    F = set()
    L = set()
    for a in S:
        for b in S:
            if a != b:
                m = (a + b) % p
                if m != 0:
                    if m in S:
                        L.add(m)
                    else:
                        F.add(m)
    print(f"  F (illegal merges) = {sorted(F)}")
    print(f"  L (legal merges) = {sorted(L)}")
    print(f"  |F|={len(F)}, |L|={len(L)}")
    print(f"  F = QNR? {F == S_comp}")

    # Full computation
    A = circulant_digraph(p, sorted(S))
    max_d = min(p - 1, 6)
    betti_full = path_betti_numbers(A, p, max_dim=max_d)
    print(f"  β = {betti_full}")

# ===== Fourier computation for larger Paley primes =====
print("\n\n--- Fourier v3 Computation ---")

for p in paley_primes:
    S = quadratic_residues(p)
    S_comp = set(range(1, p)) - S

    # F structure
    F = set()
    L = set()
    for a in S:
        for b in S:
            if a != b:
                m = (a + b) % p
                if m != 0:
                    if m in S:
                        L.add(m)
                    else:
                        F.add(m)

    max_d = min(p - 1, 6)
    betti = fourier_betti_v3(S, p, max_dim=max_d)
    print(f"\nP_{p}: QR = {sorted(S)}")
    print(f"  |F|={len(F)}, |L|={len(L)}, F=QNR? {F == S_comp}")
    print(f"  β = {betti}")

    # Check odd-only property
    odd_only = all(betti[k] == 0 for k in range(2, len(betti), 2))
    print(f"  β_even = 0? {odd_only}")

    # If not tournament, also compute for QNR
    neg_S = {(p - s) % p for s in S}
    is_tournament = (S & neg_S == set())
    if not is_tournament:
        # p ≡ 1 mod 4: QR is NOT a tournament connection set
        print(f"  NOTE: p ≡ {p%4} mod 4, not a tournament")

# ===== Compare Paley vs random tournaments =====
print("\n\n" + "=" * 70)
print("PALEY vs RANDOM: Is Paley special?")
print("=" * 70)

import random
random.seed(42)

for p in [7, 11]:
    S = quadratic_residues(p)
    A_paley = circulant_digraph(p, sorted(S))
    betti_paley = path_betti_numbers(A_paley, p, max_dim=min(p-1, 5))

    # All circulant tournaments at this prime
    print(f"\np={p}: Paley β = {betti_paley}")
    print(f"All circulant tournaments on Z_{p}:")

    # A circulant tournament needs S with |S| = (p-1)/2 and S ∩ (-S) = ∅
    all_circulant_tournaments = []
    for S_comb in combinations(range(1, p), (p-1)//2):
        S_set = set(S_comb)
        neg_S = {(p - s) % p for s in S_set}
        if S_set & neg_S == set():
            max_d = min(p - 1, 5)
            if p <= 7:
                A = circulant_digraph(p, sorted(S_set))
                betti = path_betti_numbers(A, p, max_dim=max_d)
            else:
                betti = fourier_betti_v3(S_set, p, max_dim=max_d)
            all_circulant_tournaments.append((S_set, betti))
            print(f"  S={sorted(S_set)}: β={betti}", flush=True)

    print(f"  Total circulant tournaments: {len(all_circulant_tournaments)}")

# ===== QR sum structure (Gauss sums connection) =====
print("\n\n" + "=" * 70)
print("GAUSS SUM STRUCTURE: QR sums and eigenvalues")
print("=" * 70)

for p in [7, 11, 19, 23]:
    S = quadratic_residues(p)
    print(f"\np={p}: QR={sorted(S)}")

    # For each eigenvalue λ = ω^k, compute Σ_{s∈S} λ^s (= Gauss sum character)
    omega = np.exp(2j * np.pi / p)
    for k in range(p):
        lam = omega ** k
        gauss_char = sum(lam ** s for s in S)
        # M_1(λ) = [λ^s - 1 for s ∈ S]
        M1_entries = [lam**s - 1 for s in sorted(S)]
        rank1 = 1 if any(abs(x) > 1e-8 for x in M1_entries) else 0
        if k < 5 or abs(gauss_char) > 0.1:
            print(f"  k={k}: Σ λ^s = {gauss_char:.4f}, |Σ|={abs(gauss_char):.4f}, rank(M1)={rank1}")

    # Gauss sum value
    g = sum(omega ** (a*a) for a in range(1, p))
    print(f"  Gauss sum g = {g:.4f}, |g| = {abs(g):.4f}, sqrt(p) = {p**0.5:.4f}")

print("\nDone.")
