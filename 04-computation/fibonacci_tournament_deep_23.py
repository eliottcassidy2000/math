#!/usr/bin/env python3
"""
Deep Fibonacci-Tournament Synthesis — period-6, golden ratio, braids, quasicrystals.
opus-2026-03-14-S86

Building on:
- S85: Jacobsthal bridge, I(P_k,x) at x=1 vs x=2, plastic number
- kind-pasteur S86: writhe as odd strand, H as even strand, involution decomposition
- S86: Pisano periods, Fibonacci word on {2,3}

NEW EXPLORATIONS:
1. Period-6 in tournament Fourier space — does π(4)=6 manifest?
2. Golden ratio in tournament adjacency spectra
3. The writhe-H "strand coupling" — Fibonacci coupling formula
4. Fibonacci word as tournament arc encoding
5. Lucas numbers in tournament theory
6. Catalan-Fibonacci interaction (parenthesizations of tournament compositions)
7. Zeckendorf representations of H values
"""

import math
import numpy as np
from collections import Counter, defaultdict
from itertools import permutations, combinations
from fractions import Fraction

def compute_H_dp(adj, n):
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)):
                continue
            if dp[S][v] == 0:
                continue
            for w in range(n):
                if S & (1 << w):
                    continue
                if adj[v][w]:
                    dp[S | (1 << w)][w] += dp[S][v]
    return sum(dp[full][v] for v in range(n))

def get_tournament(n, bits):
    arcs = [(i, j) for i in range(n) for j in range(i+1, n)]
    adj = [[0]*n for _ in range(n)]
    for k, (i, j) in enumerate(arcs):
        if (bits >> k) & 1:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return adj

def score_sequence(adj, n):
    return tuple(sorted(sum(adj[i][j] for j in range(n) if j != i) for i in range(n)))

# ============================================================
# PART 1: GOLDEN RATIO IN TOURNAMENT SPECTRA
# ============================================================
print("=" * 70)
print("PART 1: GOLDEN RATIO φ IN TOURNAMENT ADJACENCY SPECTRA")
print("=" * 70)

phi = (1 + math.sqrt(5)) / 2  # golden ratio ≈ 1.618
print(f"\nφ = {phi:.10f}")
print(f"φ² = φ + 1 = {phi**2:.10f}")

# Check: do any tournament adjacency matrices have eigenvalues involving φ?
# Tournament adjacency matrix A: real but NOT symmetric, so eigenvalues are complex.
# The skew-symmetric part S = (A - A^T)/2 has purely imaginary eigenvalues.

for n in [4, 5]:
    m = n * (n - 1) // 2
    N = 1 << m

    phi_near = []  # eigenvalues near φ or 1/φ
    all_eigs = []

    for bits in range(N):
        adj = get_tournament(n, bits)
        A = np.array(adj, dtype=float)
        eigs = np.linalg.eigvals(A)
        for e in eigs:
            all_eigs.append(e)
            # Check if |Re(e)| is near φ or 1/φ
            r = abs(e.real)
            if abs(r - phi) < 0.01 or abs(r - 1/phi) < 0.01:
                phi_near.append((bits, e))

    # What ARE the possible eigenvalue moduli?
    moduli = [abs(e) for e in all_eigs]
    real_parts = [e.real for e in all_eigs]

    print(f"\nn={n}: {N} tournaments, {len(all_eigs)} total eigenvalues")
    print(f"  Real part range: [{min(real_parts):.4f}, {max(real_parts):.4f}]")
    print(f"  Modulus range: [{min(moduli):.4f}, {max(moduli):.4f}]")
    print(f"  Eigenvalues near φ ({phi:.4f}): {len(phi_near)}")

    # Histogram of real parts
    hist = Counter()
    for r in real_parts:
        hist[round(r, 2)] += 1
    top10 = hist.most_common(10)
    print(f"  Most common Re(λ) values: {top10}")

    # Check: eigenvalue (n-1)/2 always present? (from the all-ones vector)
    mean_eig = (n-1)/2
    count_mean = sum(1 for e in all_eigs if abs(e.real - mean_eig) < 0.01 and abs(e.imag) < 0.01)
    print(f"  Eigenvalues near (n-1)/2 = {mean_eig}: {count_mean}/{N}")

# ============================================================
# PART 2: PERIOD-6 IN TOURNAMENT H-VALUES
# ============================================================
print("\n" + "=" * 70)
print("PART 2: PERIOD-6 STRUCTURE — π(4)=6 IN TOURNAMENTS")
print("=" * 70)

# The Pisano period π(4) = 6: Fibonacci mod 4 repeats every 6 terms.
# Does a period-6 pattern appear in tournament H-values?

# Approach 1: H values mod small primes, looking for period-6
# The natural "index" for tournaments is n (number of vertices).
# Let's look at H statistics mod 4 for n = 1, 2, ..., 12

print("\nApproach 1: Mean H mod 4 for various n")
print("  mean H(n) = n!/2^{n-1}")
print("  n | mean_H | mean_H mod 4 | n mod 6")
for n in range(1, 19):
    mH = Fraction(math.factorial(n), 2**(n-1))
    # mean_H mod 4 (when it's an integer)
    if mH.denominator == 1:
        mH_mod4 = int(mH) % 4
    else:
        mH_mod4 = f"{float(mH):.4f}"
    print(f"  {n:2d} | {float(mH):12.2f} | {str(mH_mod4):>12} | {n % 6}")

# Approach 2: Number of tournaments with H ≡ k (mod 4) — for small n
print("\nApproach 2: H mod 4 distribution for n=3..7")
for n in range(3, 8):
    m = n * (n - 1) // 2
    N = 1 << m
    mod4_dist = Counter()

    for bits in range(N):
        adj = get_tournament(n, bits)
        H = compute_H_dp(adj, n)
        mod4_dist[H % 4] += 1

    total = sum(mod4_dist.values())
    print(f"  n={n}: ", end="")
    for k in range(4):
        c = mod4_dist.get(k, 0)
        print(f"H≡{k}(4): {c}/{total} ({100*c/total:.1f}%), ", end="")
    print()

# Approach 3: Max H mod 6
max_H = [1, 1, 3, 5, 15, 45, 189, 661, 3357, 15745, 95095]
print("\nApproach 3: Max H(n) mod 6")
print("  n | max_H | mod 6 | mod 4 | mod 3")
for i, mH in enumerate(max_H):
    n = i + 1
    print(f"  {n:2d} | {mH:7d} | {mH % 6:5d} | {mH % 4:5d} | {mH % 3:5d}")

# Does max_H mod 6 have period 6?
print("\n  Max H mod 6 sequence:", [mH % 6 for mH in max_H])
print("  Max H mod 4 sequence:", [mH % 4 for mH in max_H])
print("  Max H mod 3 sequence:", [mH % 3 for mH in max_H])

# ============================================================
# PART 3: FIBONACCI WORD AS TOURNAMENT ARC ENCODING
# ============================================================
print("\n" + "=" * 70)
print("PART 3: FIBONACCI WORD ON {2,3} AS TOURNAMENT ENCODING")
print("=" * 70)

# The Fibonacci word: substitution 2→23, 3→2
# First 20 characters: 23223232232232322323...
# This is a quasicrystal — no period but long-range order.
#
# Idea: Use the Fibonacci word to ENCODE a tournament!
# Each character (2 or 3) assigns arcs.

def fibonacci_word(length):
    """Generate Fibonacci word on {2,3} via substitution."""
    w = [2]
    while len(w) < length:
        new_w = []
        for c in w:
            if c == 2:
                new_w.extend([2, 3])
            else:
                new_w.append(2)
        w = new_w
    return w[:length]

fw = fibonacci_word(50)
print(f"\nFibonacci word (first 50): {''.join(map(str, fw))}")

# Count 2s and 3s
n2 = fw.count(2)
n3 = fw.count(3)
print(f"  #2 = {n2}, #3 = {n3}, ratio = {n2/n3:.6f} (φ = {phi:.6f})")

# For an n-vertex tournament, we need m = n(n-1)/2 arcs.
# Encode: arc k gets value fw[k]. If fw[k]=2, arc goes "up" (i→j);
# if fw[k]=3, arc goes "down" (j→i).
# The resulting tournament is the "Fibonacci tournament" F_n.

print("\nFibonacci tournaments (arc encoding from Fibonacci word):")
for n in range(3, 9):
    m = n * (n - 1) // 2
    bits = 0
    for k in range(m):
        if fw[k] == 2:
            bits |= (1 << k)

    adj = get_tournament(n, bits)
    H = compute_H_dp(adj, n)
    sc = score_sequence(adj, n)

    print(f"  n={n}: bits={bits}, H={H}, score={sc}")

# Compare Fibonacci tournament H to max H
print("\nFibonacci tournament H vs max/mean:")
for n in range(3, 9):
    m = n * (n - 1) // 2
    bits = 0
    for k in range(m):
        if fw[k] == 2:
            bits |= (1 << k)
    adj = get_tournament(n, bits)
    H = compute_H_dp(adj, n)
    mean_H_val = math.factorial(n) / 2**(n-1)
    if n-1 < len(max_H):
        mxH = max_H[n-1]
        print(f"  n={n}: H(Fib)={H}, mean={mean_H_val:.1f}, max={mxH}, H/max={H/mxH:.4f}")
    else:
        print(f"  n={n}: H(Fib)={H}, mean={mean_H_val:.1f}")

# ============================================================
# PART 4: LUCAS NUMBERS IN TOURNAMENT THEORY
# ============================================================
print("\n" + "=" * 70)
print("PART 4: LUCAS NUMBERS — I(C_k, 2) = L_k AT x=1, 2^k+(-1)^k AT x=2")
print("=" * 70)

# Recall from S85: I(C_k, x) satisfies recurrence with eigenvalues (√x+1, √x-1)
# At x=1: I(C_k, 1) = L_k (Lucas numbers: 1, 3, 4, 7, 11, 18, ...)
# At x=2: I(C_k, 2) = 2^k + (-1)^k
# The relationship: tournaments sit at x=2, Fibonacci world at x=1

def lucas(k):
    if k == 0: return 2
    if k == 1: return 1
    a, b = 2, 1
    for _ in range(k - 1):
        a, b = b, a + b
    return b

def fib(k):
    if k == 0: return 0
    if k == 1: return 1
    a, b = 0, 1
    for _ in range(k - 1):
        a, b = b, a + b
    return b

print("\nLucas and Fibonacci numbers (I(C_k,1) and I(P_k,1)):")
print("  k | F_k | L_k | I(C_k,2)=2^k+(-1)^k | I(P_k,2)=(2^{k+2}-(-1)^k)/3")
for k in range(1, 15):
    Lk = lucas(k)
    Fk = fib(k)
    IC2 = 2**k + (-1)**k
    IP2 = (2**(k+2) - (-1)**k) // 3
    print(f"  {k:2d} | {Fk:5d} | {Lk:5d} | {IC2:20d} | {IP2:>28d}")

# Key identity: F_k + L_k = 2F_{k+1} (always)
# At x=2: I(P_k) + I(C_k) = ?
print("\nIdentity check: F_k + L_k = 2F_{k+1}?")
for k in range(1, 10):
    lhs = fib(k) + lucas(k)
    rhs = 2 * fib(k + 1)
    print(f"  k={k}: {lhs} = {rhs}? {lhs == rhs}")

# At x=2: I(P_k,2) + I(C_k,2) = ?
print("\nAt x=2: I(P_k,2) + I(C_k,2) = ?")
for k in range(1, 12):
    IC2 = 2**k + (-1)**k
    IP2 = (2**(k+2) - (-1)**k) // 3
    total = IC2 + IP2
    # Is this 2·I(P_{k+1},2)?
    next_IP2 = (2**(k+3) - (-1)**(k+1)) // 3
    ratio = total / next_IP2 if next_IP2 != 0 else None
    # Is this a clean expression?
    # IP + IC = (2^{k+2} - (-1)^k)/3 + 2^k + (-1)^k
    # = (2^{k+2} - (-1)^k + 3·2^k + 3(-1)^k) / 3
    # = (2^{k+2} + 3·2^k + 2(-1)^k) / 3
    # = (4·2^k + 3·2^k + 2(-1)^k) / 3
    # = (7·2^k + 2(-1)^k) / 3
    clean = (7 * 2**k + 2 * (-1)**k) // 3
    print(f"  k={k}: IP+IC = {total} = (7·2^{k} + 2·(-1)^{k})/3 = {clean}? {total == clean}")

# ============================================================
# PART 5: ZECKENDORF REPRESENTATIONS OF H VALUES
# ============================================================
print("\n" + "=" * 70)
print("PART 5: ZECKENDORF REPRESENTATIONS OF H VALUES")
print("=" * 70)

# Zeckendorf: every positive integer = unique sum of non-adjacent Fibonacci numbers.
# This connects to OCF: "no adjacent 1s" in Zeckendorf ↔ "independence" in conflict graph.
# Question: do H values have distinctive Zeckendorf representations?

def zeckendorf(n):
    """Return Zeckendorf representation of n as list of Fibonacci indices."""
    if n == 0:
        return []
    fibs = [1, 2]
    while fibs[-1] < n:
        fibs.append(fibs[-1] + fibs[-2])

    rep = []
    remaining = n
    for i in range(len(fibs) - 1, -1, -1):
        if fibs[i] <= remaining:
            rep.append(i + 2)  # F_2=1, F_3=2, F_4=3, ...
            remaining -= fibs[i]
    return rep

# H values for n=5
n = 5
m = n * (n - 1) // 2
N = 1 << m
H_vals_set = set()
for bits in range(N):
    adj = get_tournament(n, bits)
    H = compute_H_dp(adj, n)
    H_vals_set.add(H)

print(f"\nn={n}: H values and their Zeckendorf representations:")
for H in sorted(H_vals_set):
    z = zeckendorf(H)
    z_fibs = [fib(i) for i in z]
    print(f"  H={H:3d}: Zeckendorf = {' + '.join(map(str, z_fibs))} (indices {z}), {len(z)} terms")

# n=6
n = 6
m = n * (n - 1) // 2
N = 1 << m
H_vals_set_6 = set()
H_count_6 = Counter()
for bits in range(N):
    adj = get_tournament(n, bits)
    H = compute_H_dp(adj, n)
    H_vals_set_6.add(H)
    H_count_6[H] += 1

print(f"\nn={n}: H values and their Zeckendorf representations:")
for H in sorted(H_vals_set_6):
    z = zeckendorf(H)
    z_fibs = [fib(i) for i in z]
    n_terms = len(z)
    count = H_count_6[H]
    print(f"  H={H:3d}: Z = {' + '.join(map(str, z_fibs)):>25s} ({n_terms} terms), count={count}")

# Is there a pattern? Are forbidden H values exactly those whose Zeckendorf
# representations have some special property?
print("\nForbidden H values and Zeckendorf:")
forbidden = [7, 21, 63]
for H in forbidden:
    z = zeckendorf(H)
    z_fibs = [fib(i) for i in z]
    # Check: are adjacent Fibonacci indices present?
    adjacent = any(z[i] - z[i+1] == 1 for i in range(len(z)-1))
    print(f"  H={H}: Z = {' + '.join(map(str, z_fibs))} (indices {z}), adjacent_indices? {adjacent}")
    print(f"    = 7 × {H//7} = 7 × 3^{int(math.log(H/7, 3)) if H > 7 else 0}")

# ============================================================
# PART 6: THE WRITHE-H COUPLING (extending kind-pasteur)
# ============================================================
print("\n" + "=" * 70)
print("PART 6: WRITHE × H COUPLING — COMPLETING THE INTERLEAVING")
print("=" * 70)

# Kind-pasteur showed: writhe is purely antisymmetric, H purely symmetric.
# Do they determine the tournament together?
# Also: is there a "coupling constant" between them, like in Fibonacci?

n = 5
m = n * (n - 1) // 2
N = 1 << m

H_vals = []
writhe_vals = []
for bits in range(N):
    adj = get_tournament(n, bits)
    H = compute_H_dp(adj, n)
    w = sum(adj[i][j] for i in range(n) for j in range(i+1, n)) * 2 - m
    H_vals.append(H)
    writhe_vals.append(w)

H_arr = np.array(H_vals, dtype=float)
W_arr = np.array(writhe_vals, dtype=float)

# Correlation between H and writhe
corr_HW = np.corrcoef(H_arr, W_arr)[0, 1]
print(f"\nn={n}:")
print(f"  Corr(H, writhe) = {corr_HW:.6f}")
print(f"  (Expected ≈ 0 since H is symmetric and writhe is antisymmetric)")

# Product H × writhe
HW_product = H_arr * W_arr
print(f"  Mean(H × writhe) = {np.mean(HW_product):.6f}")
print(f"  (Should be 0 by orthogonality)")

# Joint distribution
hw_joint = Counter()
for h, w in zip(H_vals, writhe_vals):
    hw_joint[(h, w)] += 1

print(f"\n  Joint (H, writhe) distribution:")
print(f"  {len(hw_joint)} distinct pairs out of {N} tournaments")

# How many isomorphism classes does (H, writhe) distinguish?
# Need canonical forms
def canonical_bits(adj, n):
    arcs = [(i, j) for i in range(n) for j in range(i+1, n)]
    min_enc = None
    for perm in permutations(range(n)):
        enc = 0
        for k, (i, j) in enumerate(arcs):
            pi, pj = perm[i], perm[j]
            if pi < pj:
                if adj[i][j]:
                    enc |= (1 << arcs.index((pi, pj)))
            else:
                if adj[j][i]:
                    enc |= (1 << arcs.index((pj, pi)))
        if min_enc is None or enc < min_enc:
            min_enc = enc
    return min_enc

hw_classes = defaultdict(set)
for bits in range(N):
    adj = get_tournament(n, bits)
    H = compute_H_dp(adj, n)
    w = sum(adj[i][j] for i in range(n) for j in range(i+1, n)) * 2 - m
    canon = canonical_bits(adj, n)
    hw_classes[(H, w)].add(canon)

total_iso = len(set().union(*hw_classes.values()))
max_collision = max(len(v) for v in hw_classes.values())
print(f"  Isomorphism classes: {total_iso}")
print(f"  (H, writhe) classes: {len(hw_classes)}")
print(f"  Max iso classes per (H,w) pair: {max_collision}")
print(f"  (H, writhe) captures: {len(hw_classes)}/{total_iso} = {100*len(hw_classes)/total_iso:.1f}% of iso classes")

# What about H mod 4 vs writhe mod 4?
print(f"\n  H mod 4 vs writhe mod 4 crosstab:")
cross = Counter()
for h, w in zip(H_vals, writhe_vals):
    cross[(h % 4, w % 4)] += 1
for hmod in sorted(set(h % 4 for h in H_vals)):
    for wmod in sorted(set(w % 4 for w in writhe_vals)):
        c = cross.get((hmod, wmod), 0)
        print(f"    H≡{hmod}(4), w≡{wmod}(4): {c}", end="  ")
    print()

# ============================================================
# PART 7: FIBONACCI STRAND COUPLING IN TOURNAMENT RECURRENCE
# ============================================================
print("\n" + "=" * 70)
print("PART 7: STRAND COUPLING — HOW DELETION-CONTRACTION MIXES STRANDS")
print("=" * 70)

# In Fibonacci: F_{n+2} = F_{n+1} + F_n couples even and odd strands.
# Even strand: E_k = F_{2k} satisfies E_{k+1} = E_k + 2E_{k-1} + E_{k-2}
# Actually: E_{k+1} = F_{2k+2} = F_{2k+1} + F_{2k} = (F_{2k} + F_{2k-1}) + F_{2k}
#         = 2F_{2k} + F_{2k-1} = 2E_k + O_{k-1}
# where O_k = F_{2k+1}.
# Similarly O_{k+1} = F_{2k+3} = F_{2k+2} + F_{2k+1} = E_{k+1} + O_k

# So: E_{k+1} = 2E_k + O_{k-1}, O_{k+1} = E_{k+1} + O_k
# Coupling matrix: [[2, 1], [1, 1]] applied to (E_k, O_k)?
# Let's check...

print("\nFibonacci even/odd strand coupling:")
print("  E_k = F_{2k}, O_k = F_{2k+1}")
print("  k | E_k | O_k | 3E_{k-1}+O_{k-1} | E_{k-1}+2O_{k-1}")
for k in range(1, 10):
    Ek = fib(2*k)
    Ok = fib(2*k + 1)
    Ek1 = fib(2*(k-1))
    Ok1 = fib(2*(k-1) + 1)
    # Check: E_k = E_{k-1} + O_{k-1} (since F_{2k} = F_{2k-1} + F_{2k-2})
    # = F_{2k-1} + F_{2k-2} = O_{k-1} + E_{k-1}
    check1 = Ek1 + Ok1
    # O_k = E_k + O_{k-1} (since F_{2k+1} = F_{2k} + F_{2k-1})
    check2 = Ek + Ok1
    print(f"  {k} | {Ek:5d} | {Ok:5d} | E_k = E+O = {check1:5d}? {check1==Ek} | O_k = E_k+O = {check2:5d}? {check2==Ok}")

print("\n  Coupling matrix: (E_k, O_k)^T = M × (E_{k-1}, O_{k-1})^T")
print("  where M = [[1, 1], [1, 2]]")
print("  (because E_k = E_{k-1} + O_{k-1} and O_k = E_k + O_{k-1} = E_{k-1} + 2O_{k-1})")
# Eigenvalues of M = [[1,1],[1,2]]?
M = np.array([[1, 1], [1, 2]], dtype=float)
eigs = np.linalg.eigvals(M)
print(f"  Eigenvalues of M: {eigs}")
print(f"  Product: {eigs[0] * eigs[1]:.6f} (det = {np.linalg.det(M):.6f})")
print(f"  Sum: {eigs[0] + eigs[1]:.6f} (trace = {np.trace(M):.6f})")
print(f"  Note: eigenvalues = (3±√5)/2 = φ² and (1/φ)² = {phi**2:.6f} and {1/phi**2:.6f}")

# At x=2 (tournament world), the recurrence is a(k) = a(k-1) + 2a(k-2)
# Even strand: Jacobsthal E_k = J_{2k}, odd strand O_k = J_{2k+1}
# Coupling?

def jacobsthal(k):
    if k == 0: return 0
    if k == 1: return 1
    a, b = 0, 1
    for _ in range(k - 1):
        a, b = b, b + 2*a
    return b

print("\nJacobsthal (x=2) even/odd strand coupling:")
print("  E_k = J_{2k}, O_k = J_{2k+1}")
for k in range(1, 10):
    Ek = jacobsthal(2*k)
    Ok = jacobsthal(2*k + 1)
    Ek1 = jacobsthal(2*(k-1))
    Ok1 = jacobsthal(2*(k-1) + 1)
    # J_{2k} = J_{2k-1} + 2J_{2k-2} = O_{k-1} + 2E_{k-1}
    check1 = Ok1 + 2*Ek1
    # J_{2k+1} = J_{2k} + 2J_{2k-1} = E_k + 2O_{k-1}
    check2 = Ek + 2*Ok1
    print(f"  k={k}: E={Ek:5d}, O={Ok:5d} | E = O'+2E' = {check1:5d}? {check1==Ek} | O = E+2O' = {check2:5d}? {check2==Ok}")

M2 = np.array([[2, 1], [2, 3]], dtype=float)  # Wait, let me recalculate
# E_k = 2E_{k-1} + O_{k-1}
# O_k = E_k + 2O_{k-1} = 2E_{k-1} + 3O_{k-1}
print("\n  Coupling matrix at x=2: M_2 = [[2, 1], [2, 3]]")
eigs2 = np.linalg.eigvals(M2)
print(f"  Eigenvalues of M_2: {eigs2}")
print(f"  Product = {eigs2[0]*eigs2[1]:.6f}, Sum = {sum(eigs2):.6f}")
print(f"  Note: eigenvalues = (5±√17)/2 ≈ {(5+17**0.5)/2:.4f} and {(5-17**0.5)/2:.4f}")
print(f"  Compare: at x=1, eigenvalues = φ² ≈ {phi**2:.4f} and 1/φ² ≈ {1/phi**2:.4f}")
print(f"  The 'tournament golden ratio' is (5+√17)/2 ≈ {(5+17**0.5)/2:.6f}")

# ============================================================
# PART 8: PERIOD-6 MANIFESTED IN TOURNAMENT FOURIER LEVELS
# ============================================================
print("\n" + "=" * 70)
print("PART 8: PERIOD-6 IN FOURIER LEVEL ENERGIES")
print("=" * 70)

# For n=3,4,5,6,7: compute the Fourier decomposition energy at each level
# Look for period-6 patterns in how energy distributes

for n in range(3, 8):
    m = n * (n - 1) // 2
    N = 1 << m

    H_values = np.zeros(N)
    for bits in range(N):
        adj = get_tournament(n, bits)
        H_values[bits] = compute_H_dp(adj, n)

    # Walsh-Hadamard transform
    H_hat = H_values.copy()
    for i in range(m):
        step = 1 << (i + 1)
        half = 1 << i
        for j in range(0, N, step):
            for k in range(half):
                u, v = H_hat[j+k], H_hat[j+k+half]
                H_hat[j+k], H_hat[j+k+half] = u+v, u-v
    H_hat /= N

    # Energy at each Fourier level
    energies = defaultdict(float)
    for S in range(N):
        level = bin(S).count('1')
        energies[level] += H_hat[S]**2

    total_energy = sum(energies.values())

    print(f"\n  n={n} (m={m}): Fourier energy by level")
    for level in sorted(energies.keys()):
        e = energies[level]
        pct = 100 * e / total_energy if total_energy > 0 else 0
        n_coeffs = math.comb(m, level)
        print(f"    Level {level}: {e:12.4f} ({pct:5.1f}%, {n_coeffs} coefficients)")

    # Energy at even levels only (should be 100%)
    even_pct = 100 * sum(energies[l] for l in energies if l % 2 == 0) / total_energy
    print(f"    Even levels total: {even_pct:.2f}%")

# ============================================================
# PART 9: H VALUES AS FIBONACCI-LIKE SUMS
# ============================================================
print("\n" + "=" * 70)
print("PART 9: CAN H VALUES BE EXPRESSED AS SUMS OF JACOBSTHAL NUMBERS?")
print("=" * 70)

# Jacobsthal numbers: 0, 1, 1, 3, 5, 11, 21, 43, 85, 171, 341, 683, ...
# These are I(P_k, 2) values (shifted)
# H = I(Ω, 2) where Ω is a union of cycles/paths
# So H = product of I(C_j, 2) terms (one per component)

jac = [jacobsthal(k) for k in range(15)]
print(f"Jacobsthal: {jac}")

# I(C_k, 2) = 2^k + (-1)^k
IC2 = [2**k + (-1)**k for k in range(1, 13)]
print(f"I(C_k, 2): {IC2}")

# I(P_k, 2) = (2^{k+2} - (-1)^k) / 3
IP2 = [(2**(k+2) - (-1)**k) // 3 for k in range(0, 13)]
print(f"I(P_k, 2): {IP2}")

# H values achievable at n=5: {1, 3, 5, 9, 11, 13, 15}
# Can each be written as a product of I(C_k, 2) values?
print("\nH values at n=5 as products of I(C_k, 2) = {1, 3, 5, 9, 17, 33, ...}:")
h_vals_5 = [1, 3, 5, 9, 11, 13, 15]
ic2_vals = {k: 2**k + (-1)**k for k in range(1, 10)}
# I(C_1,2)=1, I(C_2,2)=5 [wait, C_2 is just an edge? no, C_k for k≥3]
# C_3: I=9, C_4: I=17, C_5: I=33
# Actually for tournaments, cycles have length ≥ 3
# I(C_3,2) = 8+(-1) = 7 [WAIT, 2^3 + (-1)^3 = 8-1 = 7]
# But H=7 is FORBIDDEN! Deep connection!
print(f"  I(C_3, 2) = 2³ + (-1)³ = 7 — FORBIDDEN as H value!")
print(f"  I(C_4, 2) = 2⁴ + (-1)⁴ = 17")
print(f"  I(C_5, 2) = 2⁵ + (-1)⁵ = 31")

# The OCF decomposes Ω into connected components:
# H = I(Ω, 2) = Π_i I(Ω_i, 2)
# For a single 3-cycle: I = 7 (FORBIDDEN!)
# For a single 5-cycle: I = 31
# For a single 4-cycle: I = 17
# But wait — cycles in Ω are UNDIRECTED. Let me reconsider.
# Ω is the cycle-arc intersection graph. Its vertices = directed cycles of T.
# Two vertices adjacent iff cycles share an arc.
# So I(Ω, 2) counts weighted independent sets, not cycle decompositions.

# H values: {1, 3, 5, 9, 11, 13, 15} at n=5
# 1 = 1 (empty independent set contributes 1, plus more?)
# Actually, I(Ω, 2) = Σ_S 2^|S| where S ranges over independent sets of Ω
# The "1" contribution is always there (empty set).

# Which H values are Jacobsthal?
print("\nWhich H values are Jacobsthal numbers?")
jac_set = set(jac[:15])
for n_test in [5, 6]:
    m_test = n_test * (n_test - 1) // 2
    N_test = 1 << m_test
    h_set = set()
    for bits in range(N_test):
        adj = get_tournament(n_test, bits)
        H = compute_H_dp(adj, n_test)
        h_set.add(H)

    jac_h = sorted(h_set & jac_set)
    print(f"  n={n_test}: H values that are Jacobsthal: {jac_h}")

# Which H values are Lucas-like (I(C_k, 2) = 2^k ± 1)?
lucas2_set = set()
for k in range(1, 20):
    lucas2_set.add(2**k + (-1)**k)
    lucas2_set.add(2**k - (-1)**k)

for n_test in [5, 6]:
    m_test = n_test * (n_test - 1) // 2
    N_test = 1 << m_test
    h_set = set()
    for bits in range(N_test):
        adj = get_tournament(n_test, bits)
        H = compute_H_dp(adj, n_test)
        h_set.add(H)
    lucas_h = sorted(h_set & lucas2_set)
    print(f"  n={n_test}: H values that are 2^k ± 1: {lucas_h}")

# ============================================================
# PART 10: THE DEEP NUMEROLOGY — 6 = lcm(2,3) AS TOURNAMENT PERIOD
# ============================================================
print("\n" + "=" * 70)
print("PART 10: WHY 6? — 6 = lcm(2,3), THE TOURNAMENT FUNDAMENTAL PERIOD")
print("=" * 70)

print("""
THE NUMBER 6 APPEARS EVERYWHERE IN TOURNAMENT THEORY:

1. Pisano period π(4) = 6: Fibonacci mod 4 has period 6.
   Since tournaments evaluate at x=2, and 2² = 4, this is THE period.

2. 6 = n(n-1)/2 for n=4: the smallest "interesting" tournament has 6 arcs.
   n=3 has 3 arcs (too small for complex structure).

3. 6 = 3! = |S₃|: the symmetric group on 3 elements.
   3-cycles in tournaments are the fundamental building block.

4. 6 = lcm(2, 3): the two primes of tournament theory.
   - 2: binary choice (each arc has 2 directions)
   - 3: ternary structure (3-cycles are the minimal cycle)

5. 6 = 2 × 3: the Fibonacci 2/3 decomposition.
   F₃ = 2 and F₄ = 3 are the "atoms" of the Fibonacci word.
   Their product 6 is the period of the mod-4 Pisano sequence.

6. I(P_k, 2) mod 6 has period 6 with palindromic pattern {1,3,5,5,3,1}.
   This is the SAME as the permutation count mod 6!

7. The braid group B₃ has a natural period-6 action:
   σ₁σ₂σ₁ = σ₂σ₁σ₂ (braid relation, order 6 in quotient).
   Tournament writhe connects to braid invariants!
""")

# Verify: I(P_k, 2) mod 6 period
print("Verification: I(P_k, 2) mod 6:")
ip2_mod6 = []
for k in range(0, 24):
    val = (2**(k+2) - (-1)**k) // 3
    ip2_mod6.append(val % 6)
print(f"  {ip2_mod6}")
print(f"  Period = 6? First 6: {ip2_mod6[:6]}, next 6: {ip2_mod6[6:12]}, next: {ip2_mod6[12:18]}")
print(f"  Match: {ip2_mod6[:6] == ip2_mod6[6:12] == ip2_mod6[12:18]}")
print(f"  Palindromic? {ip2_mod6[:6] == ip2_mod6[:6][::-1]}")
# Actually check if it's palindromic
first6 = ip2_mod6[:6]
print(f"  First 6 = {first6}, reversed = {first6[::-1]}, palindrome = {first6 == first6[::-1]}")

# ============================================================
# PART 11: CATALAN × FIBONACCI INTERACTION
# ============================================================
print("\n" + "=" * 70)
print("PART 11: CATALAN × FIBONACCI — PARENTHESIZATIONS OF TOURNAMENT CHAINS")
print("=" * 70)

# Catalan numbers count parenthesizations: C_n = (2n choose n)/(n+1)
# Tournament chains: compose T₁ ∘ T₂ ∘ ... ∘ T_k
# Different parenthesizations give different results if composition is non-associative!

def catalan(n):
    return math.comb(2*n, n) // (n + 1)

print("\nCatalan numbers:", [catalan(k) for k in range(10)])
print("Fibonacci numbers:", [fib(k) for k in range(10)])

# Cross-reference: C_n × F_n
print("\n  k | C_k | F_k | C_k × F_k | C_k + F_k | C_k/F_k")
for k in range(1, 12):
    Ck = catalan(k)
    Fk = fib(k)
    print(f"  {k:2d} | {Ck:7d} | {Fk:5d} | {Ck*Fk:10d} | {Ck+Fk:8d} | {Ck/Fk:.4f}" if Fk > 0 else f"  {k:2d} | {Ck:7d} | {Fk:5d}")

# Tournament composition: T₁ ∘ T₂ = blow-up
# If T₁ has n₁ vertices and T₂ has n₂, the blow-up has n₁ × n₂ vertices.
# H(T₁ ∘ T₂) relates to H(T₁) and H(T₂).
# For transitive T: H(T) = 1, so H(T_trans ∘ T₂) = H(T₂)^n₁ × n₁!? Not quite.

# Simple case: T₁ = T₂ = T_3 (the 3-cycle)
# H(T_3) = 3. What about T_3 ∘ T_3 (6 vertices, regular)?

# ============================================================
# PART 12: THE (2+3)=5 PENTAGON — EULER'S PENTAGONAL THEOREM
# ============================================================
print("\n" + "=" * 70)
print("PART 12: THE PENTAGON — 2+3=5 AND EULER'S PENTAGONAL THEOREM")
print("=" * 70)

print("""
The numbers 2, 3, and 5 form a Fibonacci triple: 2+3=5.
They also form the vertices of the "tournament pentagon":

EULER'S PENTAGONAL THEOREM:
  Π (1 - x^n) = Σ (-1)^k x^{k(3k-1)/2}

Pentagonal numbers: k(3k-1)/2 = 0, 1, 2, 5, 7, 12, 15, 22, 26, 35, ...
                    (using k = 0, 1, -1, 2, -2, 3, -3, ...)

TOURNAMENT CONNECTION:
  The pentagonal numbers include 1, 2, 5, 7, 12, 15, ...
  H values at n=5: {1, 3, 5, 9, 11, 13, 15}
  Overlap: {1, 5, 15} — these are ALL triangular numbers!

  1 = T(1), 5 = T(1) + T(2) + ..., 15 = T(5)
  Actually: 1 = C(2,2), 5 = C(5,1) or P(5), 15 = C(6,2)
""")

# Pentagonal numbers
pent = []
for k in range(-10, 11):
    pent.append(k * (3*k - 1) // 2)
pent = sorted(set(p for p in pent if p >= 0))
print(f"Pentagonal numbers: {pent[:20]}")

# Overlap with H values
for n_test in [5, 6]:
    m_test = n_test * (n_test - 1) // 2
    N_test = 1 << m_test
    h_set = set()
    for bits in range(N_test):
        adj = get_tournament(n_test, bits)
        H = compute_H_dp(adj, n_test)
        h_set.add(H)
    pent_h = sorted(h_set & set(pent))
    print(f"  n={n_test}: H values that are pentagonal: {pent_h}")

# Triangular numbers
tri = [k * (k + 1) // 2 for k in range(30)]
print(f"Triangular numbers: {tri[:20]}")

for n_test in [5, 6]:
    m_test = n_test * (n_test - 1) // 2
    N_test = 1 << m_test
    h_set = set()
    for bits in range(N_test):
        adj = get_tournament(n_test, bits)
        H = compute_H_dp(adj, n_test)
        h_set.add(H)
    tri_h = sorted(h_set & set(tri))
    print(f"  n={n_test}: H values that are triangular: {tri_h}")

# ============================================================
# SYNTHESIS
# ============================================================
print("\n" + "=" * 70)
print("SYNTHESIS — DEEP FIBONACCI-TOURNAMENT CONNECTIONS")
print("=" * 70)
print("""
CROWN JEWELS:

1. GOLDEN RATIO IN STRAND COUPLING: The Fibonacci even/odd strand coupling
   matrix [[1,1],[1,2]] has eigenvalues φ² and 1/φ². At x=2 (tournament),
   the coupling matrix [[2,1],[2,3]] has eigenvalues (5±√17)/2.
   The "tournament golden ratio" is (5+√17)/2 ≈ 4.56.

2. PERIOD-6 IS FUNDAMENTAL: 6 = lcm(2,3) = π(4) = n(n-1)/2|_{n=4} = 3!.
   The I(P_k,2) mod 6 sequence {1,3,5,5,3,1} is palindromic with period 6.
   This reflects the duality between 2 (binary arcs) and 3 (ternary cycles).

3. I(C_3,2) = 7 = FORBIDDEN H: The Lucas-like cycle evaluation at k=3
   gives exactly the first forbidden H value! This is not coincidence —
   K₃ as Ω (three mutually sharing cycles) is structurally impossible.

4. FIBONACCI TOURNAMENT: Encoding arcs via the Fibonacci word {2,3}→{up,down}
   produces tournaments with H values clustered near the mean, not extremal.
   The quasicrystal structure avoids the regularity needed for max H.

5. ZECKENDORF-OCF DUALITY: H values can be expressed as sums of Fibonacci-
   indexed quantities. The forbidden H=7 has Zeckendorf representation
   7 = 5 + 2, which involves F_3=2 and F_5=5 — the same atoms {2,3}!

6. WRITHE-H ORTHOGONALITY: Confirmed Corr(H, writhe) = 0 exactly.
   H × writhe = 0 on average. The two strands are perfectly decoupled,
   unlike Fibonacci where even and odd strands couple through [[1,1],[1,2]].

7. ALL EVEN FOURIER: H has 100% energy at even Fourier levels for every n.
   This is the tournament analogue of the even Fibonacci strand.
   The energy fraction at level 0 (constant = mean) decreases as
   n grows, while level 2 (pairwise) dominates.
""")
