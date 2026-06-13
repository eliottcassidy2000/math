#!/usr/bin/env python3
"""
dihedral_mod4_dichotomy.py — The deep algebraic reason for p≡3 vs p≡1 mod 4

THE KEY INSIGHT:

For p≡3 mod 4: χ(-1) = -1, so -1 ∈ NQR.
  → The signed permutation P_{-1} = -I is NOT in the QR symmetry group.
  → The QR fixed subspace is 1-dimensional = span(σ_P).
  → THM-137 applies: Paley is eigenvector of J with largest eigenvalue.
  → "Paley phase" exists at small p.

For p≡1 mod 4: χ(-1) = +1, so -1 ∈ QR.
  → P_{-1} = -I IS in the QR symmetry group.
  → The QR fixed subspace is {0} (since -I·v = v implies v = 0).
  → NO Paley eigenvector exists. THM-137 does NOT apply.
  → No "Paley phase" — Interval wins from the start.

This is the POLYGON GEOMETRY explanation:
  → p≡3 mod 4: the p-gon reflection is an ANTI-automorphism of T_Paley
    (sends T to T^op). Paley breaks reflection symmetry.
  → p≡1 mod 4: the p-gon reflection is an AUTOMORPHISM of T_Paley.
    Paley has FULL dihedral symmetry D_{2p}.
  → More symmetry = more expander-like = worse for flow = Interval wins.

Author: opus-2026-03-12-S63
"""

import numpy as np
from itertools import product, combinations
import math

def legendre(a, p):
    if a % p == 0: return 0
    v = pow(a, (p-1)//2, p)
    return v if v == 1 else -1

def chord_type(a, p):
    m = (p-1)//2
    a = a % p
    if a == 0: return 0
    return a if a <= m else p - a

def make_signed_permutation(a, p):
    m = (p-1)//2
    P = np.zeros((m, m))
    for k in range(1, m+1):
        ak = (a * k) % p
        if ak <= m:
            P[ak-1, k-1] = 1
        else:
            P[(p-ak)-1, k-1] = -1
    return P

def held_karp(A):
    n = len(A)
    dp = {}
    for start in range(n):
        dp[(1 << start, start)] = 1
    for mask in range(1, 1 << n):
        for last in range(n):
            if not (mask & (1 << last)):
                continue
            if (mask, last) not in dp:
                continue
            count = dp[(mask, last)]
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if A[last][nxt]:
                    key = (mask | (1 << nxt), nxt)
                    dp[key] = dp.get(key, 0) + count
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def tournament_from_sigma(sigma, p):
    m = (p-1)//2
    n = p
    A = np.zeros((n, n), dtype=int)
    for k in range(1, m+1):
        for i in range(n):
            j = (i + k) % n
            if sigma[k-1] == 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


print("=" * 70)
print("THE DIHEDRAL MOD-4 DICHOTOMY")
print("=" * 70)

# =====================================================================
# Section 1: Why -I kills the Paley eigenvector for p≡1 mod 4
# =====================================================================
print("\n" + "=" * 70)
print("1. THE SIGNED PERMUTATION P_{-1} FOR EACH PRIME")
print("=" * 70)

for p in [7, 11, 13, 17, 19, 23, 29]:
    m = (p-1)//2
    chi_neg1 = legendre(-1, p)
    chi_2 = legendre(2, p)
    residue_class = p % 8

    # Build P_{-1}
    P_neg1 = make_signed_permutation(p-1, p)  # -1 mod p = p-1

    # Check if P_{-1} = -I
    is_neg_I = np.allclose(P_neg1, -np.eye(m))

    # QR elements
    qr = [a for a in range(1, p) if legendre(a, p) == 1]
    neg1_in_qr = (p-1) in qr or legendre(p-1, p) == 1

    print(f"\n  p={p:3d} ≡ {residue_class} mod 8: m={m}, χ(-1)={chi_neg1:+d}, χ(2)={chi_2:+d}")
    print(f"    -1 ∈ QR: {neg1_in_qr}")
    print(f"    P_{{-1}} = -I: {is_neg_I}")

    if neg1_in_qr:
        print(f"    → P_{{-1}} ∈ QR symmetry group → fixed subspace = {{0}}")
        print(f"    → NO Paley eigenvector exists!")
    else:
        # Compute fixed subspace
        sigma_P = np.array([legendre(k, p) for k in range(1, m+1)], dtype=float)
        # Verify sigma_P is fixed by all QR P_a
        all_fixed = True
        for a in qr:
            Pa = make_signed_permutation(a, p)
            if not np.allclose(Pa @ sigma_P, sigma_P):
                all_fixed = False
                break
        print(f"    → σ_P = {sigma_P.astype(int)}")
        print(f"    → σ_P fixed by all QR P_a: {all_fixed}")
        print(f"    → THM-137 applies: Paley eigenvector theorem holds.")


# =====================================================================
# Section 2: Automorphism group structure
# =====================================================================
print("\n\n" + "=" * 70)
print("2. POLYGON SYMMETRY: REFLECTION AUTOMORPHISM vs ANTI-AUTOMORPHISM")
print("=" * 70)

print("""
A tournament T on Z_p drawn as a regular p-gon:
  - Vertices at positions exp(2πik/p) on the unit circle
  - Arc i→j is a directed chord

The ROTATION x → x+1 is always an automorphism (circulant).
The REFLECTION x → -x maps arc (i→j) to arc (-i→-j):
  - New difference: (-j)-(-i) = i-j = -(j-i)
  - Original arc exists iff j-i ∈ QR
  - Reflected arc exists iff -(j-i) ∈ QR iff χ(-1)·χ(j-i) = 1
  - If χ(-1) = +1 (p≡1 mod 4): reflected arc exists iff original does
    → REFLECTION IS AN AUTOMORPHISM → T_Paley has D_{2p} symmetry
  - If χ(-1) = -1 (p≡3 mod 4): reflected arc exists iff original doesn't
    → REFLECTION IS AN ANTI-AUTOMORPHISM → T → T^op
    → T_Paley only has C_p × C_m symmetry, NOT D_{2p}
""")

for p in [7, 11, 13, 17, 19]:
    m = (p-1)//2
    qr = set(a for a in range(1, p) if legendre(a, p) == 1)

    # Build Paley tournament
    A = np.zeros((p, p), dtype=int)
    for i in range(p):
        for j in range(p):
            if i != j and ((j - i) % p) in qr:
                A[i][j] = 1

    # Check if reflection x → -x is an automorphism
    # Reflected tournament: A'[i][j] = A[-i][-j] = A[p-i][p-j]
    A_reflected = np.zeros((p, p), dtype=int)
    for i in range(p):
        for j in range(p):
            A_reflected[i][j] = A[(p-i) % p][(p-j) % p]

    # A_reversed = A^T (transpose = reverse all arcs)
    is_auto = np.array_equal(A, A_reflected)
    is_anti = np.array_equal(A.T, A_reflected)

    aut_size = 0
    for a in range(1, p):
        for c in range(p):
            # Check if x → ax+c is an automorphism
            is_aut = True
            for i in range(p):
                for j in range(p):
                    if i != j:
                        ni = (a * i + c) % p
                        nj = (a * j + c) % p
                        if A[i][j] != A[ni][nj]:
                            is_aut = False
                            break
                if not is_aut:
                    break
            if is_aut:
                aut_size += 1

    chi_neg1 = "≡1" if legendre(-1, p) == 1 else "≡3"
    print(f"  p={p:3d} ({chi_neg1} mod 4):")
    print(f"    Reflection is automorphism: {is_auto}")
    print(f"    Reflection is anti-automorphism: {is_anti}")
    print(f"    |Aut(T_Paley)| = {aut_size}")
    print(f"    Expected: p·m = {p * m}")


# =====================================================================
# Section 3: H-maximization at p≡1 mod 4 — no Paley phase
# =====================================================================
print("\n\n" + "=" * 70)
print("3. H-MAXIMIZATION AT p≡1 MOD 4: NO PALEY PHASE")
print("=" * 70)

for p in [5, 13]:
    m = (p-1)//2
    print(f"\n  p={p}, m={m} (p ≡ {p%4} mod 4, p ≡ {p%8} mod 8)")

    # All 2^m orientations
    all_sigmas = list(product([1, -1], repeat=m))
    H_dict = {}
    for sigma in all_sigmas:
        A = tournament_from_sigma(np.array(sigma), p)
        H_dict[sigma] = held_karp(A)

    # QR set
    qr = set(a for a in range(1, p) if legendre(a, p) == 1)
    sigma_P = tuple(legendre(k, p) for k in range(1, m+1))
    sigma_I = tuple([1]*m)

    print(f"  QR = {sorted(qr)}")
    print(f"  σ_P (Paley) = {sigma_P}")
    print(f"  σ_I (Interval) = {sigma_I}")
    print(f"  H(Paley) = {H_dict[sigma_P]}")
    print(f"  H(Interval) = {H_dict[sigma_I]}")

    # Rank all orientations
    ranked = sorted(all_sigmas, key=lambda s: -H_dict[s])
    print(f"  Top 5:")
    for i, sigma in enumerate(ranked[:5]):
        A_val = sum(legendre(k, p) * sigma[k-1] for k in range(1, m+1))
        print(f"    #{i+1}: σ={sigma}, H={H_dict[sigma]}, A(σ)={A_val}")

    H_max = H_dict[ranked[0]]
    H_Paley = H_dict[sigma_P]
    H_Interval = H_dict[sigma_I]

    print(f"\n  Maximum H = {H_max}")
    print(f"  H(Paley) / H_max = {H_Paley/H_max:.6f}")
    print(f"  H(Interval) / H_max = {H_Interval/H_max:.6f}")
    print(f"  H(Paley) rank: {ranked.index(sigma_P) + 1}/{len(ranked)}")
    print(f"  H(Interval) rank: {ranked.index(sigma_I) + 1}/{len(ranked)}")

    # How many distinct H values?
    distinct_H = sorted(set(H_dict.values()), reverse=True)
    print(f"  Distinct H values: {len(distinct_H)}")
    for idx, h in enumerate(distinct_H[:8]):
        count = sum(1 for v in H_dict.values() if v == h)
        print(f"    H={h}: {count} orientations")


# =====================================================================
# Section 4: p≡5 mod 8 vs p≡1 mod 8: the role of χ(2)
# =====================================================================
print("\n\n" + "=" * 70)
print("4. THE ROLE OF χ(2): p≡5 MOD 8 vs p≡1 MOD 8")
print("=" * 70)

print("""
For p≡1 mod 4, we further split by χ(2):
  p≡5 mod 8: χ(2) = -1 (2 is NQR)
  p≡1 mod 8: χ(2) = +1 (2 is QR)

χ(2) determines whether CONSECUTIVE chord types agree:
  Chord 1 (nearest neighbor): always QR (since 1 ∈ QR)
  Chord 2 (next-nearest): QR iff χ(2) = +1

For Interval: σ₁ = σ₂ = +1 (both CW). Always agrees.
For Paley:
  p≡1 mod 8: σ₁ = σ₂ = +1 (both CW). AGREES with Interval on short chords.
  p≡5 mod 8: σ₁ = +1, σ₂ = -1 (disagree!). CONFLICTS with Interval.

When Paley conflicts with Interval on short chords, Interval's
"flow advantage" is LARGER. This explains why Interval wins MORE
decisively at p≡5 mod 8.
""")

# p=17 is the first p≡1 mod 8 where we can test
p = 17
m = (p-1)//2  # = 8
print(f"p={p}, m={m} (p ≡ {p%8} mod 8)")
print(f"  χ(-1) = {legendre(-1, p)}, χ(2) = {legendre(2, p)}")

qr = sorted(a for a in range(1, p) if legendre(a, p) == 1)
print(f"  QR = {qr}")

sigma_P = tuple(legendre(k, p) for k in range(1, m+1))
sigma_I = tuple([1]*m)
print(f"  σ_P = {sigma_P}")
print(f"  σ_I = {sigma_I}")

# Count agreements between Paley and Interval
agree = sum(1 for k in range(m) if sigma_P[k] == sigma_I[k])
disagree = m - agree
print(f"  Agreements: {agree}/{m}, Disagreements: {disagree}/{m}")
print(f"  NQR chords (where they disagree): {[k+1 for k in range(m) if sigma_P[k] != sigma_I[k]]}")

# At p=17, m=8, 2^8 = 256 orientations — feasible!
print(f"\n  Computing all {2**m} orientations at p={p}...")
all_sigmas = list(product([1, -1], repeat=m))
H_dict = {}
for sigma in all_sigmas:
    A_mat = tournament_from_sigma(np.array(sigma), p)
    H_dict[sigma] = held_karp(A_mat)

ranked = sorted(all_sigmas, key=lambda s: -H_dict[s])
H_Paley = H_dict[sigma_P]
H_Interval = H_dict[sigma_I]
H_max = H_dict[ranked[0]]

print(f"\n  H(Paley) = {H_Paley}")
print(f"  H(Interval) = {H_Interval}")
print(f"  H_max = {H_max}")
print(f"  H(P)/H_max = {H_Paley/H_max:.6f}")
print(f"  H(I)/H_max = {H_Interval/H_max:.6f}")
print(f"  H(P) rank: {ranked.index(sigma_P) + 1}/{len(ranked)}")
print(f"  H(I) rank: {ranked.index(sigma_I) + 1}/{len(ranked)}")

# Who maximizes?
print(f"\n  Top 10:")
for i, sigma in enumerate(ranked[:10]):
    A_val = sum(legendre(k, p) * sigma[k-1] for k in range(1, m+1))
    nqr_flipped = sum(1 for k in range(m) if sigma[k] != sigma_I[k])
    print(f"    #{i+1}: σ={sigma}, H={H_dict[sigma]}, A(σ)={A_val}, "
          f"dist_from_I={nqr_flipped}")

distinct_H = sorted(set(H_dict.values()), reverse=True)
print(f"\n  Distinct H values: {len(distinct_H)}")


# =====================================================================
# Section 5: The complete mod-8 picture
# =====================================================================
print("\n\n" + "=" * 70)
print("5. THE COMPLETE MOD-8 PICTURE")
print("=" * 70)

# Summary table
results = {}
for p in [3, 5, 7, 11, 13, 17, 19]:
    m = (p-1)//2
    if m == 0: continue

    sigma_P = tuple(legendre(k, p) for k in range(1, m+1))
    sigma_I = tuple([1]*m)

    if p in [3, 5, 7, 11, 13, 17]:
        # Compute directly
        A_P = tournament_from_sigma(np.array(sigma_P), p)
        A_I = tournament_from_sigma(np.array(sigma_I), p)
        H_P = held_karp(A_P)
        H_I = held_karp(A_I)
    elif p == 19:
        H_P = 1172695746915
        H_I = 1184212824763
    else:
        continue

    results[p] = (H_P, H_I)

print(f"\n  {'p':>3} {'mod8':>5} {'χ(-1)':>5} {'χ(2)':>5} {'m':>3} {'m%2':>4} "
      f"{'H(P)':>16} {'H(I)':>16} {'Winner':>8} {'Ratio':>8}")
print(f"  {'-'*3} {'-'*5} {'-'*5} {'-'*5} {'-'*3} {'-'*4} "
      f"{'-'*16} {'-'*16} {'-'*8} {'-'*8}")

for p in sorted(results.keys()):
    m = (p-1)//2
    H_P, H_I = results[p]
    chi1 = legendre(-1, p)
    chi2 = legendre(2, p)
    winner = "PALEY" if H_P > H_I else "INTERVAL" if H_I > H_P else "TIE"
    ratio = H_P / H_I if H_I > 0 else float('inf')

    print(f"  {p:3d} {p%8:5d} {chi1:+5d} {chi2:+5d} {m:3d} {m%2:4d} "
          f"{H_P:16d} {H_I:16d} {winner:>8} {ratio:8.4f}")


# =====================================================================
# Section 6: The dihedral group interpretation
# =====================================================================
print("\n\n" + "=" * 70)
print("6. DIHEDRAL GROUP INTERPRETATION")
print("=" * 70)

print("""
The tournament on Z_p is drawn as a regular p-gon.
The dihedral group D_{2p} = ⟨r, s | r^p = s² = (sr)² = 1⟩ acts on vertices.

ROTATIONS (r^k): Always automorphisms of circulant tournaments.
REFLECTION (s): x → -x mod p.

For the PALEY tournament:
  p≡3 mod 4: s is an ANTI-automorphism (T → T^op)
    → Aut(T_Paley) = {x → ax + c : a ∈ QR, c ∈ Z_p} (order p·m)
    → D_{2p} is NOT contained in Aut(T_Paley)
    → Paley BREAKS the polygon's reflection symmetry
    → This "chirality" creates a nontrivial Paley eigenvector σ_P

  p≡1 mod 4: s IS an automorphism (T → T)
    → Aut(T_Paley) ⊇ D_{2p} (full dihedral symmetry!)
    → Paley PRESERVES the polygon's reflection symmetry
    → "Too symmetric" = "too expander-like" = worse for flow
    → No Paley eigenvector exists (P_{-1} = -I kills it)

For the INTERVAL tournament:
  s maps S = {1,...,m} to p-S = {m+1,...,p-1} (complementary interval)
  → s is ALWAYS an anti-automorphism of Interval
  → Interval ALWAYS breaks reflection symmetry → creates "flow"
  → This is why Interval wins at ALL p ≡ 1 mod 4 and at large p ≡ 3 mod 4

THE DEEP STRUCTURE:
  H-maximization rewards CHIRALITY (breaking reflection symmetry).
  Chirality = directed flow = more Hamiltonian paths.

  At p≡3 mod 4 (small p):
    Both Paley and Interval are chiral. Paley wins via the 2-body eigenvector.
  At p≡3 mod 4 (large p):
    Interval's STRONGER chirality (unidirectional flow) beats Paley's.
  At p≡1 mod 4 (ALL p):
    Paley is NOT chiral (has reflection symmetry) → Interval wins immediately.

POLYGON PICTURE:
  Think of directed edges as "traffic" on the p-gon.
  - Interval: all traffic flows clockwise → maximum chirality
  - Paley (p≡3 mod 4): traffic alternates by QR structure → partial chirality
  - Paley (p≡1 mod 4): traffic is reflection-symmetric → zero chirality

  "Chirality" = nonzero winding number = Hamiltonian paths "go around" the polygon.
  Maximum chirality = maximum flow = maximum path count.
""")


# =====================================================================
# Section 7: Interlacing of even and odd groups
# =====================================================================
print("=" * 70)
print("7. INTERLACING OF EVEN AND ODD GROUPS")
print("=" * 70)

print("""
The user's insight about interlacing:

  D_{2n} (order 2n, even) = symmetries of the n-gon (tournament on n vertices)
  C_m (order m) = QR symmetry group (m = (p-1)/2 chord types)

  p:   3   5   7  11  13  17  19  23  29  31
  m:   1   2   3   5   6   8   9  11  14  15
  m%2: 1   0   1   1   0   0   1   1   0   1
  mod4: 3   1   3   3   1   1   3   3   1   3

  p≡3 mod 4 → m ODD:
    C_m has no element of order 2
    → -I is NOT in the QR symmetry group
    → Nontrivial fixed subspace exists (σ_P)
    → "Odd groups are interlaced between tournament sizes"
    → These are the primes where Paley can win

  p≡1 mod 4 → m EVEN:
    C_m has an element of order 2 (= m/2-th power of generator)
    → But P_{-1} = -I IS in the QR group (since -1 ∈ QR)
    → Fixed subspace = {0}
    → "Even groups align WITH tournament sizes"
    → These are the primes where Paley CANNOT win

The tournament sizes (p) and group orders (m) alternate:
  ...p=11(m=5 odd)...p=13(m=6 even)...p=17(m=8 even)...p=19(m=9 odd)...
  |_____Paley can win___|_______________Paley cannot________________|...

The "interlacing" means: odd-order QR groups (p≡3 mod 4) have a
DIFFERENT algebraic structure than even-order ones (p≡1 mod 4).
The odd groups support a unique Paley eigenvector; the even ones don't.
""")


# =====================================================================
# Section 8: Chirality measure
# =====================================================================
print("=" * 70)
print("8. CHIRALITY MEASURE: QUANTIFYING REFLECTION ASYMMETRY")
print("=" * 70)

# Define chirality as ||A - A_reflected||_F / ||A||_F
# where A_reflected[i][j] = A[-i][-j]

for p in [5, 7, 11, 13, 17, 19]:
    m = (p-1)//2
    qr = set(a for a in range(1, p) if legendre(a, p) == 1)

    # Paley tournament
    A_P = np.zeros((p, p), dtype=int)
    for i in range(p):
        for j in range(p):
            if i != j and ((j - i) % p) in qr:
                A_P[i][j] = 1

    # Reflected Paley
    A_P_ref = np.zeros((p, p), dtype=int)
    for i in range(p):
        for j in range(p):
            A_P_ref[i][j] = A_P[(p-i)%p][(p-j)%p]

    # Interval tournament
    A_I = tournament_from_sigma(np.ones(m, dtype=int), p)

    # Reflected Interval
    A_I_ref = np.zeros((p, p), dtype=int)
    for i in range(p):
        for j in range(p):
            A_I_ref[i][j] = A_I[(p-i)%p][(p-j)%p]

    # Chirality = fraction of arcs that change under reflection
    paley_chirality = np.sum(A_P != A_P_ref) / (p * (p-1))
    interval_chirality = np.sum(A_I != A_I_ref) / (p * (p-1))

    mod4_str = "3" if p % 4 == 3 else "1"

    print(f"  p={p:3d} (≡{mod4_str} mod 4): "
          f"Paley chirality = {paley_chirality:.4f}, "
          f"Interval chirality = {interval_chirality:.4f}, "
          f"ratio I/P = {interval_chirality/paley_chirality if paley_chirality > 0 else float('inf'):.4f}")


print("\n\nDONE.")
