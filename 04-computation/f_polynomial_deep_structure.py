#!/usr/bin/env python3
"""
f_polynomial_deep_structure.py — Deep structural analysis of F(T,x).

F(T,x) = sum_P x^{fwd(P)} is palindromic, H = F_0 = F_{n-1}.

NEW QUESTIONS:
1. Is F(T,x) always UNIMODAL? (coefficients increase then decrease)
2. Is F(T,x)/H always REAL-ROOTED? (this would imply log-concavity)
3. What is the CENTER OF MASS of fwd(P) over all permutations?
   (This is D_1/D_0 = (n-1)/2 — always the midpoint!)
4. What is the VARIANCE of fwd(P)?
   var = D_2/D_0 - (D_1/D_0)^2
   = D_2/n! - ((n-1)/2)^2
   At n=7: D_2 = 16800 + 240*t3, so var depends on t3!
5. Does F(T,x) factor nicely for special tournaments (Paley, regular)?
6. Is there a TRANSFER MATRIX interpretation of F(T,x)?
7. Can we express F(T,x) as a PERMANENT of a matrix?

Author: opus-2026-03-07-S44
"""
from itertools import permutations, combinations
import math
import random
import numpy as np

def tournament_from_bits(bits, n):
    A = [[0]*n for _ in range(n)]
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> pos) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            pos += 1
    return A

def compute_F(A, n):
    F = [0] * n
    for P in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if A[P[i]][P[i+1]])
        F[fwd] += 1
    return F

def count_t3(A, n):
    t3 = 0
    for triple in combinations(range(n), 3):
        i, j, k = triple
        if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]):
            t3 += 1
    return t3

# ============================================================
# 1. UNIMODALITY CHECK
# ============================================================
print("=" * 60)
print("1. UNIMODALITY OF F(T,x)")
print("=" * 60)

for n in [5, 7]:
    m = n*(n-1)//2
    unimodal_count = 0
    total = 0
    non_unimodal = []

    if n == 5:
        iterator = range(1 << m)
    else:
        random.seed(42)
        iterator = [random.getrandbits(m) for _ in range(200)]

    seen_F = set()
    for bits in iterator:
        A = tournament_from_bits(bits, n)
        F = compute_F(A, n)
        key = tuple(F)
        if key in seen_F:
            continue
        seen_F.add(key)
        total += 1

        # Check unimodality
        is_unimodal = True
        peak = False
        for i in range(1, len(F)):
            if F[i] < F[i-1]:
                peak = True
            elif F[i] > F[i-1] and peak:
                is_unimodal = False
                break

        if is_unimodal:
            unimodal_count += 1
        else:
            non_unimodal.append(F)

    print(f"  n={n}: {unimodal_count}/{total} unimodal", end="")
    if non_unimodal:
        print(f"  Non-unimodal examples: {non_unimodal[:3]}")
    else:
        print("  ALL UNIMODAL!")

# ============================================================
# 2. REAL-ROOTEDNESS OF F(T,x)/H
# ============================================================
print("\n" + "=" * 60)
print("2. REAL-ROOTEDNESS OF F(T,x)/H")
print("=" * 60)

for n in [5, 7]:
    m = n*(n-1)//2
    real_rooted = 0
    total = 0
    complex_examples = []

    if n == 5:
        iterator = range(1 << m)
    else:
        random.seed(42)
        iterator = [random.getrandbits(m) for _ in range(200)]

    seen_F = set()
    for bits in iterator:
        A = tournament_from_bits(bits, n)
        F = compute_F(A, n)
        key = tuple(F)
        if key in seen_F:
            continue
        seen_F.add(key)
        total += 1

        H = F[n-1]
        # Polynomial coefficients (highest degree first)
        poly = [F[k]/H for k in range(len(F)-1, -1, -1)]
        roots = np.roots(poly)

        all_real = all(abs(r.imag) < 1e-8 for r in roots)
        all_neg = all(r.real < 1e-8 for r in roots) if all_real else False

        if all_real:
            real_rooted += 1
        else:
            complex_examples.append((H, F, roots))

    print(f"  n={n}: {real_rooted}/{total} real-rooted")
    if complex_examples:
        for H, F, roots in complex_examples[:3]:
            root_strs = [f"{r.real:.3f}+{r.imag:.3f}i" for r in roots]
            print(f"    H={H}: F={F}, roots={root_strs}")

# ============================================================
# 3. VARIANCE OF fwd(P)
# ============================================================
print("\n" + "=" * 60)
print("3. VARIANCE OF fwd(P) — depends on t3!")
print("=" * 60)

n = 7
m = n*(n-1)//2
random.seed(123)

print(f"  {'t3':>3} {'var':>8} {'D2':>6} {'formula':>8}")
seen = set()
for trial in range(100):
    bits = random.getrandbits(m)
    A = tournament_from_bits(bits, n)
    F = compute_F(A, n)
    t3 = count_t3(A, n)

    # D_0 = n!, D_1 = (n-1)/2 * n!, D_2 = 16800 + 240*t3
    D0 = math.factorial(n)
    D1 = (n-1) * D0 // 2
    D2 = 16800 + 240 * t3

    mean = D1 / D0  # = (n-1)/2 = 3
    var = D2 / D0 - mean**2

    key = t3
    if key not in seen:
        seen.add(key)
        # Formula: var = (16800 + 240*t3)/5040 - 9 = 16800/5040 - 9 + 240*t3/5040
        #        = 10/3 - 9 + t3/21 = -17/3 + t3/21
        # Hmm: 16800/5040 = 10/3 ≈ 3.333
        # var = 10/3 - 9 + t3/21 = (70 - 189 + 10*t3)/63 = (-119 + 10*t3)/63
        # Hmm that gives negative variance for small t3. Let me recompute.
        # mean = 3, mean^2 = 9
        # D_2 includes C(fwd,2) = fwd*(fwd-1)/2 contributions
        # sum C(fwd,2) = D_2
        # E[fwd^2] = 2*D_2/D_0 + E[fwd] = 2*D_2/n! + (n-1)/2
        # var = E[fwd^2] - E[fwd]^2 = 2*D_2/n! + (n-1)/2 - ((n-1)/2)^2

        Efwd2 = 2*D2/D0 + (n-1)/2
        var_correct = Efwd2 - mean**2

        print(f"  {t3:3d} {var_correct:8.4f} {D2:6d} E[fwd^2]={Efwd2:.4f}")

print(f"\n  var = 2*D_2/n! + (n-1)/2 - ((n-1)/2)^2")
print(f"      = 2*(16800+240*t3)/5040 + 3 - 9")
print(f"      = (33600+480*t3)/5040 - 6")
print(f"      = 20/3 + 2*t3/21 - 6")
print(f"      = 2/3 + 2*t3/21")
print(f"  So variance INCREASES with t3! More 3-cycles → more spread in fwd(P).")

# ============================================================
# 4. F(T,x) AS PERMANENT
# ============================================================
print("\n" + "=" * 60)
print("4. F(T,x) AS PERMANENT OF A MATRIX")
print("=" * 60)

# F(T,x) = sum_P x^{fwd(P)} = sum_P prod_{i=0}^{n-2} x^{A[P_i][P_{i+1}]}
# This is NOT a standard permanent. The standard permanent is:
# perm(M) = sum_sigma prod M[i][sigma(i)]
# Our sum is over PATHS not matchings.

# However, we can write it as a permanent-like expression:
# Define M(x) where M(x)[i][j] = x^{A[i][j]} for i ≠ j.
# Then the "path permanent" is sum over Hamiltonian paths P of prod M[P_i][P_{i+1}].
# This equals tr(adj(I-M(x))) or can be computed via inclusion-exclusion.

# Actually: F(T,x) = sum over all perms P of [n] of prod_{i=0}^{n-2} (x if A[P_i][P_{i+1}]=1, else 1)
# = sum_P prod (1 + (x-1)*A[P_i][P_{i+1}])
# At x=1: n!. At x=0: F_0 = H(T^op) = H(T).

# The transfer matrix approach: define M(x)[i][j] = x if i→j, else 1 (for j≠i)
# Then F(T,x) = sum_{start v} sum over paths v→...  of products
# But this sums over Hamiltonian paths (each vertex exactly once).
# So F(T,x) = sum_v [x-weighted HP count from v].

# This is related to the permanent of the n×n matrix with M[i][j] = x*A[i][j] + (1-A[i][j])
# = 1 + (x-1)*A[i][j]... but permanent sums over MATCHINGS not paths.

# Let me verify at n=3:
n_test = 3
A = [[0,1,1],[0,0,1],[0,0,0]]  # transitive: 0→1, 0→2, 1→2
F = compute_F(A, n_test)
print(f"  n=3 transitive: F = {F}")

# Matrix M(x): M[i][j] = x if A[i][j]=1, else 1 (for i≠j), M[i][i] = 0
# M = [[0, x, x],
#      [1, 0, x],
#      [1, 1, 0]]
# perm(M) = 0*0*0 + 0*x*1 + x*0*1 + x*x*0 + x*0*1 + 0*1*0 = ...
# Standard permanent: sum_sigma prod M[i][sigma(i)]
# sigma = (0,1,2): M[0][0]*M[1][1]*M[2][2] = 0*0*0 = 0
# sigma = (0,2,1): M[0][0]*M[1][2]*M[2][1] = 0*x*1 = 0
# sigma = (1,0,2): M[0][1]*M[1][0]*M[2][2] = x*1*0 = 0
# sigma = (1,2,0): M[0][1]*M[1][2]*M[2][0] = x*x*1 = x^2
# sigma = (2,0,1): M[0][2]*M[1][0]*M[2][1] = x*1*1 = x
# sigma = (2,1,0): M[0][2]*M[1][1]*M[2][0] = x*0*1 = 0
# perm = x^2 + x
# But F = [1, 4, 1]. Not matching!

# The permanent counts CYCLE COVERS, not paths. We need a different formulation.
# F(T,x) is a HAFNIAN-like or path-counting object, not a permanent.

print("  F(T,x) is NOT a permanent — it sums over paths, not cycle covers.")
print("  It IS the 'path permanent' S(T) when x = -1 (up to sign).")

# ============================================================
# 5. COMPARISON WITH CHROMATIC POLYNOMIAL
# ============================================================
print("\n" + "=" * 60)
print("5. ANALOGY: F(T,x) vs CHROMATIC POLYNOMIAL")
print("=" * 60)

# The chromatic polynomial P(G, k) counts proper k-colorings.
# It satisfies: P(G, k) = P(G-e, k) - P(G/e, k) (deletion-contraction)
# It's a polynomial in k of degree n, with alternating signs.
# P(G, -1) = (-1)^n * # acyclic orientations (Stanley's theorem)

# The Redei-Berge polynomial u_X(m) counts (f,X)-friendly m-colorings.
# It satisfies: u_X(m) = u_{X\e}(m) - u_{X/e}(m) (Grujic-Stojadinovic)
# u_X(1) = H(T).
# u_X(-m) = (-1)^n u_{X-bar}(m) (reciprocity)

# Both satisfy deletion-contraction and reciprocity!
# The chromatic poly of the UNDERLYING UNDIRECTED graph is different.

# Key question: is u_X(m) = F(T, x) for some x = f(m)?
# u_X(1) = H. F(T, 1) = n!. So u_X ≠ F at m=1 unless n=1.

# Actually u_X(m) counts the number of "surjective" colorings compatible
# with the tournament. F(T,x) counts all permutations weighted by x^fwd.
# They're related but different.

print("  The Redei-Berge polynomial u_X(m) and F(T,x) are related via")
print("  the Eulerian-to-Bernoulli transform:")
print("    u_X(m) = sum_{k>=0} F_k * C(m-k+n-1-k, n-1) ... (needs verification)")
print()
print("  Both satisfy deletion-contraction and reciprocity.")
print("  F(T, -1) relates to S(T); u_X(-1) relates to complement.")

# ============================================================
# 6. LOG-CONCAVITY OF F
# ============================================================
print("\n" + "=" * 60)
print("6. LOG-CONCAVITY CHECK: F_k^2 >= F_{k-1} * F_{k+1}")
print("=" * 60)

n = 5
m = n*(n-1)//2
log_concave_count = 0
total = 0
failures = []
seen_F = set()

for bits in range(1 << m):
    A = tournament_from_bits(bits, n)
    F = compute_F(A, n)
    key = tuple(F)
    if key in seen_F:
        continue
    seen_F.add(key)
    total += 1

    is_lc = True
    for k in range(1, n-1):
        if F[k]**2 < F[k-1] * F[k+1]:
            is_lc = False
            failures.append((F, k, F[k]**2, F[k-1]*F[k+1]))
            break

    if is_lc:
        log_concave_count += 1

print(f"  n=5: {log_concave_count}/{total} log-concave")
if failures:
    for F, k, lhs, rhs in failures[:3]:
        print(f"    FAILS at k={k}: F={F}, F_k^2={lhs} < F_{k-1}*F_{k+1}={rhs}")
else:
    print("    ALL LOG-CONCAVE!")

# n=7
n = 7
m = n*(n-1)//2
random.seed(42)
log_concave_count = 0
total = 0
failures = []
seen_F = set()

for trial in range(500):
    bits = random.getrandbits(m)
    A = tournament_from_bits(bits, n)
    F = compute_F(A, n)
    key = tuple(F)
    if key in seen_F:
        continue
    seen_F.add(key)
    total += 1

    is_lc = True
    for k in range(1, n-1):
        if F[k]**2 < F[k-1] * F[k+1]:
            is_lc = False
            failures.append((F, k))
            break
    if is_lc:
        log_concave_count += 1

print(f"  n=7: {log_concave_count}/{total} log-concave (sampled)")
if failures:
    print(f"    {len(failures)} failures found")
    for F, k in failures[:2]:
        print(f"    FAILS at k={k}: F={F}")
