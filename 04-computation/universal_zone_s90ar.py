#!/usr/bin/env python3
"""
universal_zone_s90ar.py — The universal zone and its boundaries
opus-2026-03-16-S90ar

The trace triangle T(k,n) = Tr(M_k^n) has a beautiful structure:
- Universal zone (n <= k): T(k,n) = 2^n - 1 (Mersenne)
- Boundary (n = k): T(k,k) = 2^k - 1 (diagonal = Mersenne)
- Individual zone (n > k): T(k,n) depends on k
- First departure: T(k,k+1) = 2^{k+1} - k - 2

Let's explore the FULL individual zone and find patterns.
"""

from sympy import Matrix, Rational
from math import log, log2

print("""


          THE UNIVERSAL ZONE AND ITS BOUNDARIES


""")

# ========================================================================
# Build the trace triangle with exact arithmetic
# ========================================================================

def knacci_matrix_sym(k):
    M_data = []
    for i in range(k):
        row = []
        for j in range(k):
            if i == 0:
                row.append(1)
            elif j == i - 1:
                row.append(1)
            else:
                row.append(0)
        M_data.append(row)
    return Matrix(M_data)

# Build using recurrence instead of matrix power (much faster)
def trace_sequence(k, length):
    """Compute Tr(M_k^n) for n = 0, 1, ..., length-1 using Newton + recurrence."""
    t = [0] * length
    t[0] = k  # Tr(I) = k

    # Newton range: n = 1, ..., min(k, length-1)
    for n in range(1, min(k, length)):
        # t(n) = sum_{m=1}^{n-1} t(m) + n (Newton with all-(-1) coefficients)
        t[n] = sum(t[m] for m in range(1, n)) + n

    # Boundary: n = k (also Newton but the last one)
    if k < length:
        t[k] = sum(t[m] for m in range(1, k)) + k

    # Recurrence range: n > k
    for n in range(k + 1, length):
        t[n] = sum(t[n - j] for j in range(1, k + 1))

    return t

# ========================================================================
# Part 1: Display the extended triangle
# ========================================================================

print("=" * 100)
print("Part 1: Extended trace triangle T(k,n) = Tr(M_k^n)")
print("=" * 100)

max_k = 10
max_n = 15

print(f"\n{'k':>3} |", end="")
for n in range(max_n):
    print(f"  n={n:2d}", end="")
print()
print("-" * (6 + 6 * max_n))

for k in range(2, max_k + 1):
    seq = trace_sequence(k, max_n)
    print(f"{k:3d} |", end="")
    for n in range(max_n):
        v = seq[n]
        if 1 <= n <= k and v == 2**n - 1:
            # In universal zone - mark with *
            print(f" *{v:>3}", end="")
        else:
            print(f"  {v:>3}", end="")
    print()

print("\n(* marks the universal zone where T(k,n) = 2^n - 1)")

# ========================================================================
# Part 2: The deficit triangle D(k,n) = (2^n - 1) - T(k,n)
# ========================================================================

print("\n" + "=" * 100)
print("Part 2: Deficit D(k,n) = (2^n - 1) - T(k,n) in the individual zone")
print("=" * 100)

print(f"\n{'k':>3} |", end="")
for n in range(max_n):
    print(f"  n={n:2d}", end="")
print()
print("-" * (6 + 6 * max_n))

for k in range(2, max_k + 1):
    seq = trace_sequence(k, max_n)
    print(f"{k:3d} |", end="")
    for n in range(max_n):
        mersenne = 2**n - 1
        deficit = mersenne - seq[n]
        if n == 0:
            print(f"  {'-':>3}", end="")
        elif deficit == 0:
            print(f"    .", end="")
        else:
            print(f"  {deficit:>3}", end="")
    print()

print("""
The deficit is 0 in the universal zone and grows in the individual zone.

First deficit column (n = k+1): D(k,k+1) = k+1.
Second deficit column (n = k+2): ?
""")

# Compute deficit at n = k+2
print("Deficit at n = k+2:")
for k in range(2, 12):
    seq = trace_sequence(k, k + 3)
    d1 = (2**(k+1) - 1) - seq[k+1]
    d2 = (2**(k+2) - 1) - seq[k+2]
    print(f"  k={k:2d}: D(k,k+1) = {d1:4d}, D(k,k+2) = {d2:4d}, ratio D2/D1 = {d2/d1:.4f}")

print("""
D(k,k+1) = k + 1 exactly.
D(k,k+2) grows... let me check the formula.

At n = k+2: t(k+2) = t(k+1) + t(k) + ... + t(2)
  = (2^{k+1} - k - 2) + (2^k - 1) + ... + (2^2 - 1)
  = (2^{k+1} - k - 2) + sum_{m=2}^{k} (2^m - 1)
  = (2^{k+1} - k - 2) + (2^{k+1} - 4 - (k-1))
  = 2^{k+2} - 2k - 5.

So D(k,k+2) = (2^{k+2} - 1) - (2^{k+2} - 2k - 5) = 2k + 4.
""")

# Verify
print("Verification of D(k,k+2) = 2k + 4:")
for k in range(2, 12):
    seq = trace_sequence(k, k + 3)
    d2 = (2**(k+2) - 1) - seq[k+2]
    predicted = 2*k + 4
    print(f"  k={k:2d}: D(k,k+2) = {d2:4d}, predicted = {predicted:4d}, match: {d2 == predicted}")

# ========================================================================
# Part 3: General deficit formula
# ========================================================================

print("\n" + "=" * 100)
print("Part 3: General deficit D(k,k+j) for small j")
print("=" * 100)

# Compute D(k,k+j) for j = 1, 2, 3, 4, 5
for j in range(1, 7):
    print(f"\nj={j}: D(k,k+{j}) for k = 2..10:")
    vals = []
    for k in range(2, 11):
        seq = trace_sequence(k, k + j + 1)
        d = (2**(k+j) - 1) - seq[k+j]
        vals.append(d)
        print(f"  k={k}: {d}", end="")
    print()

    # Check if polynomial in k
    # First differences
    diffs = [vals[i+1] - vals[i] for i in range(len(vals)-1)]
    print(f"  First diffs: {diffs}")
    diffs2 = [diffs[i+1] - diffs[i] for i in range(len(diffs)-1)]
    print(f"  Second diffs: {diffs2}")

# ========================================================================
# Part 4: The growth rate in the individual zone
# ========================================================================

print("\n" + "=" * 100)
print("Part 4: Asymptotic growth rate in the individual zone")
print("=" * 100)

print("""
In the universal zone (n <= k): T(k,n) = 2^n - 1, growth rate = 2.
In the individual zone (n > k): T(k,n) follows k-nacci recurrence,
  growth rate = tau_k < 2.

The growth rate DROPS from 2 to tau_k at the transition n = k.

The ratio T(k,n+1)/T(k,n) should approach tau_k as n -> inf:
""")

for k in [3, 5, 7]:
    seq = trace_sequence(k, 30)
    print(f"\nk={k}: T(k,n+1)/T(k,n):")
    for n in range(1, 25):
        if seq[n] > 0:
            ratio = seq[n+1] / seq[n]
            zone = "U" if n < k else "I"
            print(f"  n={n:2d}: {ratio:.8f} [{zone}]", end="")
            if n == k - 1:
                print(" <-- transition", end="")
            if n == k:
                print(" <-- FIRST individual", end="")
            print()

import numpy as np
for k in [3, 5, 7]:
    M = np.zeros((k, k))
    for j in range(k):
        M[0][j] = 1
    for i in range(1, k):
        M[i][i-1] = 1
    eigs = np.linalg.eigvals(M)
    tau = max(abs(e) for e in eigs)
    print(f"  tau_{k} = {tau:.8f}")

print("""
The transition is SHARP:
  - In the universal zone: ratio = 2^{n+1}/(2^n) = 2.0 (exactly, for large n).
  - At n = k: ratio drops abruptly.
  - For n >> k: ratio approaches tau_k (the k-nacci constant).

The universal zone is a "regime" where the system behaves as if
it has infinite dimension (all k-nacci agree, rate = 2).
The individual zone is where finite dimension matters.

The TRANSITION at n = k is a phase transition:
  - Growth rate drops from 2 to tau_k.
  - The deficit D(k,n) starts growing.
  - The system "discovers" its own dimension.
""")

# ========================================================================
# Part 5: Connection to the Collatz-like iteration
# ========================================================================

print("=" * 100)
print("Part 5: The iteration x -> 2x - 1 vs k-nacci")
print("=" * 100)

print("""
From the augmented polynomial: L^{k+1} = 2*L^k - 1.

So the eigenvalues satisfy the COLLATZ-LIKE iteration:
  x_{n+1} = 2*x_n - 1.

Starting from x_0 = L_i (an eigenvalue of M_k):
  x_1 = 2*L_i - 1
  x_2 = 2*(2*L_i - 1) - 1 = 4*L_i - 3
  x_3 = 8*L_i - 7
  x_n = 2^n * L_i - (2^n - 1) = 2^n * (L_i - 1) + 1.

So: L_i^{k+n} = x_n = 2^n * (L_i - 1) + 1  (if the relation holds repeatedly).

Wait, that's wrong. The relation L^{k+1} = 2*L^k - 1 is an eigenvalue identity,
not a recurrence for iterating on a SINGLE eigenvalue.

For the DOMINANT eigenvalue tau_k:
  tau_k^{k+1} = 2*tau_k^k - 1.

But tau_k^{k+2} != 2*tau_k^{k+1} - 1 in general.
Instead, tau_k^{k+2} = tau_k * tau_k^{k+1} (by definition).

Let me check:
""")

for k in [3, 5, 7]:
    M = np.zeros((k, k))
    for j in range(k):
        M[0][j] = 1
    for i in range(1, k):
        M[i][i-1] = 1
    eigs = np.linalg.eigvals(M)
    tau = max(e.real for e in eigs if abs(e.imag) < 1e-10)

    # Check tau^{k+1} = 2*tau^k - 1
    lhs = tau**(k+1)
    rhs = 2 * tau**k - 1
    print(f"  k={k}: tau^{{k+1}} = {lhs:.8f}, 2*tau^k - 1 = {rhs:.8f}, match: {abs(lhs - rhs) < 1e-6}")

    # Check tau^{k+2} = 2*tau^{k+1} - 1?
    lhs2 = tau**(k+2)
    rhs2 = 2 * tau**(k+1) - 1
    print(f"       tau^{{k+2}} = {lhs2:.8f}, 2*tau^{{k+1}} - 1 = {rhs2:.8f}, match: {abs(lhs2 - rhs2) < 1e-6}")
    # tau satisfies the char poly, so tau^{k+2} = tau*tau^{k+1}
    rhs3 = tau * tau**(k+1)
    print(f"       tau*tau^{{k+1}} = {rhs3:.8f}, match: {abs(lhs2 - rhs3) < 1e-6}")

print("""
Correct: tau^{k+1} = 2*tau^k - 1 (from augmented poly).
But tau^{k+2} != 2*tau^{k+1} - 1 (the relation only holds once).

The augmented polynomial gives a SINGLE-STEP relation
from power k to power k+1, not an iterated one.

However, combining with the k-nacci recurrence:
  Tr(M^{k+j}) follows a linear recurrence (the k-nacci).
  The augmented relation Tr(M^{k+1}) = 2*Tr(M^k) - k
  is the FIRST STEP of this recurrence.
""")

# ========================================================================
# Part 6: The full picture
# ========================================================================

print("=" * 100)
print("Part 6: THE FULL PICTURE")
print("=" * 100)

print("""
The trace triangle T(k,n) has three regimes:

REGIME 1: UNIVERSAL (1 <= n <= k)
  T(k,n) = 2^n - 1 (Mersenne numbers)
  Independent of k
  Growth rate: 2 (each step doubles and subtracts 1)
  Origin: Newton's identities with all-(-1) coefficients

REGIME 2: TRANSITION (n = k to n = k + O(1))
  T(k,k) = 2^k - 1 (last Mersenne)
  T(k,k+1) = 2^{k+1} - k - 2 (first departure, deficit k+1)
  T(k,k+2) = 2^{k+2} - 2k - 5 (deficit 2k+4, linear in k)
  Growth rate: dropping from 2 toward tau_k

REGIME 3: INDIVIDUAL (n >> k)
  T(k,n) ~ C * tau_k^n (dominated by largest eigenvalue)
  Growth rate: tau_k < 2
  Different for each k

THE TRANSITION AT n = k IS A PHASE TRANSITION.

In physics language:
  - Universal zone: the system is "in the universality class of 2"
  - Individual zone: the system shows its own critical exponent tau_k
  - The transition at n = k: the system "discovers its dimension"

In information theory:
  - Universal zone: k bits of behavior are fully determined
  - Individual zone: additional bits require knowledge of k
  - At the boundary: 2^k - 1 is the last universally determined trace value

In tournament theory:
  - The transfer matrix M_k encodes a k-step comparison structure
  - For the first k powers, ALL k-step structures behave identically
  - Beyond k powers, the specific k matters

THE MERSENNE NUMBER 2^k - 1 IS THE AMOUNT OF UNIVERSAL BEHAVIOR
THAT A k-DIMENSIONAL SYSTEM EXHIBITS BEFORE ITS INDIVIDUALITY EMERGES.
""")

# ========================================================================
# Part 7: The phi special case
# ========================================================================

print("=" * 100)
print("Part 7: The golden ratio special case (k=2)")
print("=" * 100)

print("""
At k=2: M_2 = [[1,1],[1,0]] (the Fibonacci matrix!)

tau_2 = phi = (1+sqrt(5))/2 = 1.618...

Tr(M_2^n) = Lucas numbers: 2, 1, 3, 4, 7, 11, 18, 29, ...

Universal zone (n=1,2): Tr = 1, 3 = 2^1-1, 2^2-1. MERSENNE!
Individual zone (n>2): Tr = Lucas numbers (NOT Mersenne).

tau_2^2 = phi^2 = phi + 1 = 2.618 = 1/(2-phi) = 1/0.382.
  = 2^2 - 1 + (phi - 1 - 1) = 3 + (phi - 2) = 3 - 0.382.
  Hmm, not quite. Tr(M_2^2) = 3 = 2^2 - 1 exactly.
  But phi^2 = 2.618 != 3. The trace includes ALL eigenvalues:
  phi^2 + psi^2 = 2.618 + 0.382 = 3. Yes!

The Fibonacci matrix is the k=2 case of the universal trace theorem.
Lucas numbers BEGIN with Mersenne (1, 3) and then diverge.
The first 2 Lucas numbers are the first 2 Mersenne numbers.

For k=2, tau_2^2 = 1/(2-phi). And:
  1/(2-phi) = 1/(3-phi-1) = 1/((3-phi-1)) ... hmm.
  2 - phi = 2 - (1+sqrt(5))/2 = (4-1-sqrt(5))/2 = (3-sqrt(5))/2.
  1/(2-phi) = 2/(3-sqrt(5)) = 2(3+sqrt(5))/((3-sqrt(5))(3+sqrt(5))) = 2(3+sqrt(5))/4
            = (3+sqrt(5))/2 = phi + 1 = phi^2.

So phi^2 = (3+sqrt(5))/2. Known identity.
The identity tau_k^k = 1/(2-tau_k) generalizes phi^2 = phi+1 = 1/(2-phi)!

AT k=2: The GOLDEN RATIO identity phi^2 = phi + 1 IS the k=2 case
of the k-nacci Mersenne identity Tr(M_k^k) = 2^k - 1.

THIS IS BEAUTIFUL: the most famous identity in number theory
(phi^2 = phi + 1) is a special case of our theorem.
""")

print("opus-2026-03-16-S90ar: THE UNIVERSAL ZONE AND ITS BOUNDARIES")
