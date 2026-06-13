#!/usr/bin/env python3
"""
TANG-YAU SYMBOL MATRIX ANALYSIS

Deep dive into how the Fourier symbol matrices M_p(λ) control path homology.

For C_n^S:
- M_1(λ) = row vector [λ^{s_i} - 1]
- β_0 = #{λ : λ^s = 1 for all s ∈ S} = # {k : k*gcd(S) ≡ 0 (mod n)}
- β_1 depends on ker(M_1) and im(M_2)

KEY FORMULA: dim(ker M_1(λ)) = |S| - rank(M_1(λ))
  = |S| - 1 if any λ^{s_i} ≠ 1
  = |S| if λ^{s_i} = 1 for all i (but then rank = 0)

For 2-generators S={a,b}:
  M_1(λ) = [λ^a - 1, λ^b - 1]
  ker(M_1) has dim |S|-1 = 1 when rank = 1
  So β_1 contribution from eigenspace λ = 1 - rank(M_2(λ))

This means β_1 = Σ_λ max(0, ker_dim - im_dim) where
the im_dim comes from Ω_2 and the M_2 symbol matrix.

ANALYSIS: When does β_2 appear? When can β_3 appear?
What is the role of the "illegal merged steps" F = (S+S)\\S?
"""
import numpy as np
from itertools import combinations
import sys
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_fourier_v3 import fourier_betti_v3, enumerate_step_sequences

def analyze_symbol(S_set, n, verbose=True):
    """Analyze the symbol matrix structure for C_n^S."""
    roots = [np.exp(2j * np.pi * k / n) for k in range(n)]
    S_list = sorted(S_set)
    d = len(S_list)

    # Illegal merged steps
    F = set()
    for a in S_set:
        for b in S_set:
            if a != b:
                merged = (a + b) % n
                if merged != 0 and merged not in S_set:
                    F.add(merged)

    # Legal merged steps (in S)
    L = set()
    for a in S_set:
        for b in S_set:
            if a != b:
                merged = (a + b) % n
                if merged != 0 and merged in S_set:
                    L.add(merged)

    if verbose:
        print(f"\n  C_{n}^{S_set}:")
        print(f"    |S|={d}, F(illegal)={F}, L(legal)={L}")

    # Step sequences
    omega2_seqs = enumerate_step_sequences(S_set, n, 2)
    omega3_seqs = enumerate_step_sequences(S_set, n, 3)

    if verbose:
        print(f"    |A_2 step seqs|={len(omega2_seqs)}, |A_3 step seqs|={len(omega3_seqs)}")

    # Per-eigenvalue analysis
    eigenvalue_data = []
    for k in range(n):
        lam = roots[k]
        # M_1(λ)
        M1 = np.array([[lam**s - 1 for s in S_list]])
        rank1 = np.linalg.matrix_rank(M1, tol=1e-8)
        ker1 = d - rank1
        eigenvalue_data.append({
            'k': k, 'rank1': rank1, 'ker1': ker1,
        })

    if verbose:
        # Show which eigenspaces contribute to β_0 and β_1
        b0_contributions = [e for e in eigenvalue_data if e['rank1'] == 0]
        b1_ker = sum(e['ker1'] for e in eigenvalue_data)
        print(f"    β_0 contributions: {len(b0_contributions)} eigenspaces")
        print(f"    Total ker(∂_1) = {b1_ker}")

    # Compute full Betti
    betti = fourier_betti_v3(S_set, n, max_dim=min(n-1, 6))

    if verbose:
        print(f"    β = {betti}")

    return betti, F, L


# ===== 1. Classification by F (illegal merges) =====
print("=" * 70)
print("SYMBOL MATRIX ANALYSIS: ROLE OF ILLEGAL MERGES F = (S+S)\\S")
print("=" * 70)

print("\n--- n=7 (prime): 2-generator sets ---")
for S in combinations(range(1, 7), 2):
    analyze_symbol(set(S), 7)

print("\n\n--- n=7: 3-generator sets ---")
for S in combinations(range(1, 7), 3):
    analyze_symbol(set(S), 7)

# ===== 2. The F-size vs topology connection =====
print("\n\n" + "=" * 70)
print("|F| (# illegal merges) vs TOPOLOGY")
print("=" * 70)

for n in [5, 7, 11]:
    print(f"\nn={n}:")
    data_by_F = {}
    for size in range(1, min(n, 5)):
        for S in combinations(range(1, n), size):
            S_set = set(S)
            betti, F, L = analyze_symbol(S_set, n, verbose=False)
            key = (len(S_set), len(F))
            if key not in data_by_F:
                data_by_F[key] = []
            data_by_F[key].append((S_set, betti))

    for key in sorted(data_by_F.keys()):
        examples = data_by_F[key]
        betti_types = {}
        for S_set, betti in examples:
            bt = tuple(betti)
            while bt and bt[-1] == 0: bt = bt[:-1]
            betti_types[bt] = betti_types.get(bt, 0) + 1
        print(f"  |S|={key[0]}, |F|={key[1]}: {len(examples)} sets, types={dict(betti_types)}", flush=True)

# ===== 3. When |F|=0 (S is closed under pairwise sums mod n) =====
print("\n\n" + "=" * 70)
print("CLOSED CONNECTION SETS: F = ∅ (S+S ⊆ S ∪ {0})")
print("=" * 70)

for n in [5, 7, 11, 13]:
    print(f"\nn={n}:")
    for size in range(1, min(n, 6)):
        for S in combinations(range(1, n), size):
            S_set = set(S)
            F = set()
            for a in S_set:
                for b in S_set:
                    if a != b:
                        merged = (a + b) % n
                        if merged != 0 and merged not in S_set:
                            F.add(merged)
            if not F:
                betti = fourier_betti_v3(S_set, n, max_dim=min(n-1, 5))
                print(f"  S={S_set}: β={betti}")

# ===== 4. The s+s' ∈ S condition and torus detection =====
print("\n\n" + "=" * 70)
print("TORUS DETECTION: When does β=(1,2,1,...)?")
print("=" * 70)

print("\nTang-Yau: For S={a,b} with gcd(a,b,n)=1 and a+b ∉ S:")
print("  The torus appears when F = {a+b mod n}.")
print("  Equivalently: exactly ONE illegal merge.")

for n in [7, 11, 13, 17]:
    print(f"\nn={n}:")
    torus_count = 0
    total = 0
    for S in combinations(range(1, n), 2):
        S_set = set(S)
        a, b = S
        total += 1
        merged = (a + b) % n
        F_size = 0
        for x in S_set:
            for y in S_set:
                if x != y:
                    m = (x + y) % n
                    if m != 0 and m not in S_set:
                        F_size += 1

        if n <= 13:
            betti = fourier_betti_v3(S_set, n, max_dim=min(n-1, 4))
            is_torus = tuple(betti[:3]) == (1, 2, 1)
        else:
            is_torus = None

        if is_torus and F_size <= 4:
            torus_count += 1

    if n <= 13:
        print(f"  {torus_count}/{total} are tori", flush=True)

print("\nDone.")
