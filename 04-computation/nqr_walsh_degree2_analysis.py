#!/usr/bin/env python3
"""
nqr_walsh_degree2_analysis.py — Analyzing the degree-2 NQR Walsh sum

The Walsh Difference Formula says:
  H(Int) - H(Pal) = 2 · Σ_{ψ(S)=-1} hat{H}[S]

At degree 2: hat{H}[{i,j}] = J[i,j] = interaction matrix entry.
The degree-2 NQR sum is: Σ_{ψ({i,j})=-1} J[i,j]

ψ({i,j}) = legendre(i, p) * legendre(j, p)
  = +1 iff i,j are BOTH QR or BOTH NQR
  = -1 iff exactly one of i,j is NQR

So the degree-2 NQR sum = Σ over "mixed" pairs {QR, NQR}.

For p=7, m=3: QR∩{1,2,3} = {1,2}, NQR∩{1,2,3} = {3}
  NQR pairs: {1,3}, {2,3}
  QR pairs: {1,2}

The sign of Σ_NQR J[i,j] tells us: do mixed QR-NQR chord pairs
contribute positively or negatively to H?

Author: opus-2026-03-12-S65
"""

import numpy as np
from itertools import product as iprod, combinations

def legendre(a, p):
    if a % p == 0: return 0
    v = pow(a, (p-1)//2, p)
    return v if v == 1 else -1

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

def held_karp(A):
    n = len(A)
    dp = {}
    for start in range(n):
        dp[(1 << start, start)] = 1
    for mask in range(1, 1 << n):
        for last in range(n):
            if not (mask & (1 << last)): continue
            if (mask, last) not in dp: continue
            count = dp[(mask, last)]
            for nxt in range(n):
                if mask & (1 << nxt): continue
                if A[last][nxt]:
                    key = (mask | (1 << nxt), nxt)
                    dp[key] = dp.get(key, 0) + count
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

print("=" * 70)
print("DEGREE-2 NQR WALSH SUM: THE J MATRIX CONNECTION")
print("=" * 70)

for p in [7, 11, 13]:
    m = (p-1)//2
    print(f"\n{'='*50}")
    print(f"p={p}, m={m}")
    print(f"{'='*50}")

    # Classify chords as QR or NQR
    qr_chords = [k for k in range(1, m+1) if legendre(k, p) == 1]
    nqr_chords = [k for k in range(1, m+1) if legendre(k, p) == -1]
    print(f"  QR chords: {qr_chords}")
    print(f"  NQR chords: {nqr_chords}")

    # Compute all H values
    all_sigs = list(iprod([1, -1], repeat=m))
    H_vals = {}
    for sig in all_sigs:
        A = tournament_from_sigma(np.array(sig), p)
        H_vals[sig] = held_karp(A)

    # Compute J matrix (Walsh degree-2 coefficients)
    J = np.zeros((m, m))
    for i in range(m):
        for j in range(i, m):
            # hat{H}[{i,j}] = (1/2^m) Σ_σ H(σ) σ_i σ_j
            coeff = 0
            for sig in all_sigs:
                coeff += H_vals[sig] * sig[i] * sig[j]
            coeff /= (1 << m)
            J[i][j] = coeff
            J[j][i] = coeff

    print(f"\n  J matrix (interaction):")
    for i in range(m):
        row = " ".join(f"{J[i][j]:>8.1f}" for j in range(m))
        chi_i = "QR" if legendre(i+1, p) == 1 else "NQR"
        print(f"    chord {i+1} ({chi_i}): [{row}]")

    # Classify pairs
    print(f"\n  Pair classification:")
    qr_qr_sum = 0
    nqr_nqr_sum = 0
    mixed_sum = 0
    diag_qr = 0
    diag_nqr = 0

    for i in range(m):
        chi_i = legendre(i+1, p)
        # Diagonal
        if chi_i == 1:
            diag_qr += J[i][i]
        else:
            diag_nqr += J[i][i]

        for j in range(i+1, m):
            chi_j = legendre(j+1, p)
            psi = chi_i * chi_j
            if psi == 1:
                if chi_i == 1:
                    qr_qr_sum += J[i][j]
                else:
                    nqr_nqr_sum += J[i][j]
            else:
                mixed_sum += J[i][j]

    print(f"  QR×QR pairs: Σ J = {qr_qr_sum:.2f} ({len(list(combinations(qr_chords, 2)))} pairs)")
    print(f"  NQR×NQR pairs: Σ J = {nqr_nqr_sum:.2f} ({len(list(combinations(nqr_chords, 2)))} pairs)")
    print(f"  Mixed pairs: Σ J = {mixed_sum:.2f} ({len(qr_chords)*len(nqr_chords)} pairs)")
    print(f"  QR diagonal: Σ J_ii = {diag_qr:.2f}")
    print(f"  NQR diagonal: Σ J_ii = {diag_nqr:.2f}")

    # The degree-2 NQR Walsh sum combines diagonal + off-diagonal NQR terms
    # ψ({i}) = legendre(i,p). But we only care about degree 2 here.
    # For degree 2: ψ({i,j}) = legendre(i*j, p)
    # NQR pairs at degree 2: those with chi_i * chi_j = -1

    nqr_deg2_sum = mixed_sum
    print(f"\n  Degree-2 NQR sum = {nqr_deg2_sum:.2f}")
    print(f"  Contribution to H(Int)-H(Pal) = {2*nqr_deg2_sum:.2f}")

    # Compare with Q(sigma_P) = sigma_P^T J sigma_P
    sig_pal = tuple(legendre(k, p) for k in range(1, m+1))
    if p % 4 == 3:
        Q_pal = sum(J[i][j] * sig_pal[i] * sig_pal[j]
                    for i in range(m) for j in range(m))
        Q_int = sum(J[i][j] for i in range(m) for j in range(m))
        print(f"\n  Q(Paley) = σ_P^T J σ_P = {Q_pal:.2f}")
        print(f"  Q(Interval) = 1^T J 1 = {Q_int:.2f}")
        print(f"  Q(Pal) - Q(Int) = {Q_pal - Q_int:.2f}")

        # Relationship: Q(Pal) - Q(Int) relates to the degree-2 contribution
        # Q(σ) = Σ_{i,j} J[i,j] σ_i σ_j
        # Q(Pal) = Σ J[i,j] χ(i) χ(j)
        # Q(Int) = Σ J[i,j]
        # Q(Pal) - Q(Int) = Σ J[i,j] (χ(i)χ(j) - 1)
        #   = -2 Σ_{mixed pairs} J[i,j] - 2 Σ_{NQR diag} J[i,i]
        # Wait, need to be more careful with signs.

        # χ(i)χ(j) = 1 for QR-QR or NQR-NQR pairs, -1 for mixed
        # So: Q(Pal) - Q(Int) = -2 Σ_{mixed} J[i,j] for off-diag
        #   plus diagonal: χ(i)² = 1, so diagonal contributes 0

        # Actually for ALL pairs (i,j) including i=j:
        # If both QR: χ(i)χ(j)=1, contributes 0 to Q(Pal)-Q(Int)
        # If mixed: χ(i)χ(j)=-1, contributes -2J[i,j]
        # If both NQR: χ(i)χ(j)=1, contributes 0
        # Diagonal: χ(i)²=1, contributes 0

        diff_formula = -2 * sum(J[i][j] for i in range(m) for j in range(m)
                                if legendre(i+1, p) * legendre(j+1, p) == -1)
        print(f"\n  Q(Pal)-Q(Int) via formula: {diff_formula:.2f}")
        print(f"  Direct computation:        {Q_pal - Q_int:.2f}")

        # So the degree-2 contribution to H(Pal)-H(Int) = -(Q(Pal)-Q(Int))
        # Wait, let's be precise.
        # H(Int) - H(Pal) at degree 2 = 2 * Σ_{ψ(S)=-1, |S|=2} hat{H}[S]
        # ψ({i,j}) = χ(i)χ(j) = -1 for mixed pairs
        # hat{H}[{i,j}] = J[i,j]
        # So deg-2 contribution = 2 * Σ_{mixed (i,j), i<j} J[i,j]

        print(f"\n  RELATIONSHIP:")
        print(f"  deg-2 H(Int)-H(Pal) = 2 * Σ_mixed J = {2*mixed_sum:.2f}")
        print(f"  Q(Pal) - Q(Int) = -2 * Σ_all_mixed J = {-2*sum(J[i][j] for i in range(m) for j in range(m) if i!=j and legendre(i+1,p)*legendre(j+1,p)==-1):.2f}")

# Now the key question: what determines the sign of the mixed sum?
print("\n" + "=" * 70)
print("THE KEY STRUCTURAL QUESTION")
print("=" * 70)
print("""
At degree 2, the NQR Walsh sum = Σ_{mixed pairs {i,j}} J[i,j]

J[i,j] measures the INTERACTION between chord i and chord j in H:
  J[i,j] > 0: flipping both chords i,j together INCREASES H more than expected
  J[i,j] < 0: flipping both chords together DECREASES H

For Paley to win: mixed pairs must have NEGATIVE interaction (J < 0)
This means: flipping a QR chord and an NQR chord together is BAD for H.
Equivalently: it's GOOD to have QR and NQR chords pointing in OPPOSITE directions.
This is exactly what Paley does (σ_P has QR chords = +1, NQR chords = -1).

For Interval to win: mixed pairs must have POSITIVE interaction
This means: flipping QR and NQR chords together is GOOD for H.
Equivalently: it's GOOD to have ALL chords pointing the SAME direction.
This is exactly what Interval does (all σ = +1).

The TRANSITION happens when the J[i,j] entries for mixed pairs collectively
flip from negative to positive average.
""")

# Check which specific J entries are positive/negative for mixed pairs
for p in [7, 11]:
    m = (p-1)//2
    all_sigs = list(iprod([1, -1], repeat=m))
    H_vals = {}
    for sig in all_sigs:
        A = tournament_from_sigma(np.array(sig), p)
        H_vals[sig] = held_karp(A)

    J = np.zeros((m, m))
    for i in range(m):
        for j in range(m):
            coeff = 0
            for sig in all_sigs:
                coeff += H_vals[sig] * sig[i] * sig[j]
            J[i][j] = coeff / (1 << m)

    print(f"\np={p}: Mixed pair J entries:")
    for i in range(m):
        for j in range(i+1, m):
            chi_i = legendre(i+1, p)
            chi_j = legendre(j+1, p)
            if chi_i * chi_j == -1:
                print(f"  J[{i+1},{j+1}] = {J[i][j]:>8.2f}  "
                      f"(chord {i+1} {'QR' if chi_i==1 else 'NQR'}, "
                      f"chord {j+1} {'QR' if chi_j==1 else 'NQR'})")

print("\nDONE.")
