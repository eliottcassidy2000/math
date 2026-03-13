#!/usr/bin/env python3
"""
walsh_character_excess.py — Why is the QR/NQR product excess EXACTLY -1?

DISCOVERY: Among all nonempty subsets S ⊂ {1,...,m} where m=(p-1)/2:
  #{S : prod(S) ∈ QR} = (2^m - 2)/2 = 2^{m-1} - 1
  #{S : prod(S) ∈ NQR} = (2^m - 2)/2 + 1 = 2^{m-1}
  (The empty set is excluded; prod(∅) = 1 ∈ QR would make it 2^{m-1}.)

Actually including empty set: QR count = 2^{m-1}, NQR = 2^{m-1}.
But excluding it: QR = 2^{m-1} - 1, NQR = 2^{m-1}.

PROOF ATTEMPT: Consider the map ψ: P({1,...,m}) → {±1} defined by
  ψ(S) = legendre(prod(S), p) = prod_{k∈S} legendre(k, p)

This is a GROUP HOMOMORPHISM from (F_2^m, XOR) to ({±1}, ×).

The kernel has index at most 2. If ψ is surjective, kernel has index 2
and exactly half the subsets map to +1. If ψ is trivial (always +1),
all map to +1.

ψ is trivial iff legendre(k, p) = +1 for ALL k = 1,...,m.
This happens iff all of 1,...,m are QR mod p.

For p=7: QR = {1,2,4}, {1,2,3} not all QR (3 is NQR).
For p=11: QR = {1,3,4,5,9}, {1,2,3,4,5} not all QR (2 is NQR).
For p=13: QR = {1,3,4,9,10,12}, {1,2,3,4,5,6} not all QR.

So ψ is surjective for ALL primes p > 3 (since 2 is NQR for p = 3 mod 8,
and for p = 1 mod 8, some k ≤ m is always NQR).

When ψ is surjective: kernel has size 2^{m-1}, image of +1 has 2^{m-1} subsets.
So #{prod(S) ∈ QR} = 2^{m-1} (including S=∅).
Excluding S=∅: #{prod(S) ∈ QR} = 2^{m-1} - 1, #{prod(S) ∈ NQR} = 2^{m-1}.

THIS PROVES THE EXACT -1 EXCESS.

CONSEQUENCE FOR H:
  H(Paley) = Σ_S hat{H}[S] · ψ(S)
           = Σ_{ψ(S)=+1} hat{H}[S] - Σ_{ψ(S)=-1} hat{H}[S]

  H(Interval) = Σ_S hat{H}[S]

  So: H(Int) - H(Pal) = 2 · Σ_{ψ(S)=-1} hat{H}[S]

  The difference is TWICE the sum of Walsh coefficients over NQR-product subsets.

Author: opus-2026-03-12-S65
"""

import numpy as np
from itertools import combinations

def legendre(a, p):
    if a % p == 0: return 0
    v = pow(a, (p-1)//2, p)
    return v if v == 1 else -1

print("=" * 70)
print("WALSH CHARACTER BALANCE: PROOF AND CONSEQUENCES")
print("=" * 70)

print("\n--- Part 1: Verify the group homomorphism structure ---")
for p in [5, 7, 11, 13, 17, 19, 23, 29, 31]:
    m = (p-1)//2
    leg_vec = tuple(legendre(k, p) for k in range(1, m+1))

    # Count QR and NQR products (including empty set)
    qr_count = 0  # empty set has product 1 ∈ QR
    nqr_count = 0
    for mask in range(0, 1 << m):
        prod_val = 1
        for j in range(m):
            if mask & (1 << j):
                prod_val *= leg_vec[j]
        if prod_val == 1:
            qr_count += 1
        else:
            nqr_count += 1

    # ψ surjective?
    surjective = nqr_count > 0

    print(f"  p={p:>2}, m={m}: QR={qr_count}, NQR={nqr_count}, "
          f"total={qr_count+nqr_count}=2^{m}={1<<m}, "
          f"surjective={surjective}, "
          f"QR=2^{m-1}={1<<(m-1)}? {qr_count == 1<<(m-1)}")

print("\n--- Part 2: The Walsh difference formula ---")
print("""
THEOREM: H(Interval) - H(Paley) = 2 · Σ_{S: ψ(S)=-1} hat{H}[S]

where hat{H}[S] is the Walsh-Fourier coefficient of H at subset S.

Since hat{H}[S] depends only on the STRUCTURE of S (number of even-length
path components, total degree), and the NQR subsets are essentially
"half of all subsets", this gives:

  H(Int) - H(Pal) = 2 · (sum of hat{H} over half the subsets)

At small m: the specific WHICH subsets are NQR matters (structured)
At large m: approaches random → H(Int) - H(Pal) ≈ 0 (relative to H)

But the data shows H(Int)/H(Pal) → 1.06... So the difference is O(H).
This means the NQR subsets are NOT random — they are biased toward
subsets with LARGER hat{H}[S]!
""")

# Verify at p=7: compute hat{H}[S] for all S, split by ψ
print("--- Part 3: Verification at p=7 ---")
p = 7
m = 3

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

from itertools import product as iprod

# Compute H for all orientations
all_sigs = list(iprod([1, -1], repeat=m))
H_vals = {}
for sig in all_sigs:
    A = tournament_from_sigma(np.array(sig), p)
    H_vals[sig] = held_karp(A)

# Walsh-Fourier transform
walsh_coeffs = {}
for mask in range(1 << m):
    S = tuple(j for j in range(m) if mask & (1 << j))
    coeff = 0
    for sig in all_sigs:
        chi = 1
        for j in S:
            chi *= sig[j]
        coeff += H_vals[sig] * chi
    walsh_coeffs[S] = coeff / (1 << m)

# Now split by ψ
leg_vec = tuple(legendre(k, p) for k in range(1, m+1))

print(f"\np={p}, m={m}")
print(f"  H(Interval) = {H_vals[(1,1,1)]}")

# Paley orientation
sig_pal = tuple(legendre(k, p) for k in range(1, m+1))
print(f"  H(Paley) = {H_vals[sig_pal]}")
print(f"  Difference = {H_vals[(1,1,1)] - H_vals[sig_pal]}")

nqr_sum = 0
qr_sum = 0
print(f"\n  {'S':>12} {'|S|':>4} {'hat{H}[S]':>12} {'ψ(S)':>6}")
for mask in range(1 << m):
    S = tuple(j for j in range(m) if mask & (1 << j))
    coeff = walsh_coeffs[S]

    psi = 1
    for j in S:
        psi *= leg_vec[j]

    label_S = '{' + ','.join(str(j+1) for j in S) + '}' if S else '∅'
    print(f"  {label_S:>12} {len(S):>4} {coeff:>12.2f} {'+1' if psi == 1 else '-1':>6}")

    if psi == 1:
        qr_sum += coeff
    else:
        nqr_sum += coeff

print(f"\n  Σ hat{{H}}[S] for ψ=+1: {qr_sum:.2f}")
print(f"  Σ hat{{H}}[S] for ψ=-1: {nqr_sum:.2f}")
print(f"  H(Pal) = qr_sum - nqr_sum = {qr_sum - nqr_sum:.2f}")
print(f"  H(Int) = qr_sum + nqr_sum = {qr_sum + nqr_sum:.2f}")
print(f"  Difference = 2 * nqr_sum = {2 * nqr_sum:.2f}")
print(f"  Actual diff = {H_vals[(1,1,1)] - H_vals[sig_pal]}")

# Do the same at p=11
print("\n--- Part 4: Verification at p=11 ---")
p = 11
m = 5

all_sigs = list(iprod([1, -1], repeat=m))
H_vals = {}
for sig in all_sigs:
    A = tournament_from_sigma(np.array(sig), p)
    H_vals[sig] = held_karp(A)

walsh_coeffs = {}
for mask in range(1 << m):
    S = tuple(j for j in range(m) if mask & (1 << j))
    coeff = 0
    for sig in all_sigs:
        chi = 1
        for j in S:
            chi *= sig[j]
        coeff += H_vals[sig] * chi
    walsh_coeffs[S] = coeff / (1 << m)

leg_vec = tuple(legendre(k, p) for k in range(1, m+1))
sig_pal = tuple(legendre(k, p) for k in range(1, m+1))

print(f"\np={p}, m={m}")
print(f"  H(Interval) = {H_vals[tuple(1 for _ in range(m))]}")
print(f"  H(Paley) = {H_vals[sig_pal]}")
print(f"  Difference = {H_vals[tuple(1 for _ in range(m))] - H_vals[sig_pal]}")

# Split by ψ and analyze by Walsh degree
qr_by_deg = {}
nqr_by_deg = {}
for mask in range(1 << m):
    S = tuple(j for j in range(m) if mask & (1 << j))
    coeff = walsh_coeffs[S]
    deg = len(S)

    psi = 1
    for j in S:
        psi *= leg_vec[j]

    if psi == 1:
        qr_by_deg[deg] = qr_by_deg.get(deg, 0) + coeff
    else:
        nqr_by_deg[deg] = nqr_by_deg.get(deg, 0) + coeff

print(f"\n  Walsh difference by degree:")
print(f"  {'deg':>4} {'Σ_QR':>14} {'Σ_NQR':>14} {'2*NQR':>14}")
total_qr = 0
total_nqr = 0
for deg in range(m+1):
    qr_val = qr_by_deg.get(deg, 0)
    nqr_val = nqr_by_deg.get(deg, 0)
    total_qr += qr_val
    total_nqr += nqr_val
    print(f"  {deg:>4} {qr_val:>14.2f} {nqr_val:>14.2f} {2*nqr_val:>14.2f}")

print(f"  {'total':>4} {total_qr:>14.2f} {total_nqr:>14.2f} {2*total_nqr:>14.2f}")
print(f"  Check: H(Pal) = {total_qr - total_nqr:.2f}, H(Int) = {total_qr + total_nqr:.2f}")

# KEY: Which degrees contribute most to the difference?
print(f"\n  Most important Walsh degrees for H(Int)-H(Pal):")
for deg in range(m+1):
    nqr_val = nqr_by_deg.get(deg, 0)
    pct = 2*nqr_val / (H_vals[tuple(1 for _ in range(m))] - H_vals[sig_pal]) * 100 if H_vals[tuple(1 for _ in range(m))] != H_vals[sig_pal] else 0
    print(f"  degree {deg}: 2*Σ_NQR = {2*nqr_val:>10.2f} ({pct:>7.1f}% of diff)")

# NEW QUESTION: Is it always true that H(Pal) > H(Int) iff nqr_sum < 0?
print(f"\n--- Part 5: Sign of NQR sum determines winner ---")
for p in [7, 11]:
    m = (p-1)//2
    leg_vec = tuple(legendre(k, p) for k in range(1, m+1))

    all_sigs = list(iprod([1, -1], repeat=m))
    H_vals = {}
    for sig in all_sigs:
        A = tournament_from_sigma(np.array(sig), p)
        H_vals[sig] = held_karp(A)

    sig_int = tuple(1 for _ in range(m))
    sig_pal = tuple(legendre(k, p) for k in range(1, m+1))

    walsh_coeffs = {}
    for mask in range(1 << m):
        S = tuple(j for j in range(m) if mask & (1 << j))
        coeff = 0
        for sig in all_sigs:
            chi = 1
            for j in S:
                chi *= sig[j]
            coeff += H_vals[sig] * chi
        walsh_coeffs[S] = coeff / (1 << m)

    nqr_sum = sum(walsh_coeffs[S] for mask in range(1 << m)
                  for S in [tuple(j for j in range(m) if mask & (1 << j))]
                  if all(True for _ in [None])  # just compute psi
                  and (lambda s: (1 if not s else
                      (1 if all(True for _ in [None]) else -1)))(None) == 1  # placeholder
                  )

    # Simpler computation
    nqr_sum = 0
    for mask in range(1 << m):
        S = tuple(j for j in range(m) if mask & (1 << j))
        psi = 1
        for j in S:
            psi *= leg_vec[j]
        if psi == -1:
            nqr_sum += walsh_coeffs[S]

    print(f"  p={p}: NQR sum = {nqr_sum:.2f}, "
          f"sign = {'<0' if nqr_sum < 0 else '>0'}, "
          f"H(Int)-H(Pal) = {H_vals[sig_int] - H_vals[sig_pal]}")

print("""
CONCLUSION:
  - At p=7: NQR sum < 0 → H(Pal) > H(Int) (Paley wins)
  - At p=11: NQR sum < 0 → H(Pal) > H(Int) (Paley wins)
  - At large p: NQR sum > 0 (expected) → H(Int) > H(Pal)

The sign of the NQR Walsh sum is the EXACT criterion for the phase transition.
The question reduces to: when does Σ_{ψ(S)=-1} hat{H}[S] change sign?
""")

print("\nDONE.")
