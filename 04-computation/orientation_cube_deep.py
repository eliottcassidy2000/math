#!/usr/bin/env python3
"""
orientation_cube_deep.py — Deep analysis of H on the orientation cube {+1,-1}^m

At p=7 (m=3): H = 178.5 + 3.5*s1*s2 - 3.5*s1*s3 - 3.5*s2*s3
  Only degree 0 and 2 terms! Paley maximizes by choosing sigma = (+1,+1,-1).

At p=11 (m=5): 2^5 = 32 orientations, still feasible for exact HP via Held-Karp.
  Can we find the full Walsh expansion? What degree terms appear?

KEY QUESTION: Does the Walsh expansion illuminate WHY Paley wins at p=7,11
but loses at p=19?

Also: compute p=13 interval H (kind-pasteur's S56 open question).

Author: opus-2026-03-12-S62
"""

import numpy as np
import time
from itertools import combinations

def adjacency_matrix(S, p):
    A = np.zeros((p, p), dtype=np.int8)
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in S:
                A[i][j] = 1
    return A

def count_hp_fast(A):
    """Held-Karp DP for Hamiltonian path count."""
    n = len(A)
    full = (1 << n) - 1
    dp = np.zeros((1 << n, n), dtype=np.int64)
    for v in range(n):
        dp[1 << v, v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            c = dp[mask, v]
            if c == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v, u]:
                    dp[mask | (1 << u), u] += c
    return int(np.sum(dp[full]))

def is_qr(a, p):
    if a % p == 0: return False
    return pow(a, (p - 1) // 2, p) == 1

def sigma_to_S(sigma, p):
    """Convert orientation vector to connection set."""
    m = (p - 1) // 2
    S = set()
    for k in range(1, m + 1):
        if sigma[k - 1] == 1:
            S.add(k)
        else:
            S.add(p - k)
    return frozenset(S)

def walsh_transform(H_vals, m):
    """Compute Walsh-Fourier coefficients of H on {+1,-1}^m."""
    coeffs = {}
    for r in range(m + 1):
        for S in combinations(range(m), r):
            coeff = 0
            for bits in range(2**m):
                sigma = tuple(1 if (bits >> i) & 1 else -1 for i in range(m))
                prod = 1
                for k in S:
                    prod *= sigma[k]
                coeff += H_vals[sigma] * prod
            coeff /= 2**m
            if abs(coeff) > 0.001:
                coeffs[S] = coeff
    return coeffs


# ============================================================
# SECTION 1: p=11, m=5, all 32 orientations
# ============================================================

print("=" * 70)
print("ORIENTATION CUBE ANALYSIS")
print("=" * 70)

p = 11
m = 5
S_paley = frozenset(j for j in range(1, p) if is_qr(j, p))

print(f"\np={p}: Computing H for all 2^{m} = {2**m} orientations...")
t0 = time.time()

H_vals_11 = {}
H_list = []

for bits in range(2**m):
    sigma = tuple(1 if (bits >> i) & 1 else -1 for i in range(m))
    S = sigma_to_S(list(sigma), p)
    A = adjacency_matrix(S, p)
    H = count_hp_fast(A)
    H_vals_11[sigma] = H
    H_list.append((sigma, S, H))

elapsed = time.time() - t0
print(f"Done in {elapsed:.1f}s")

# Sort by H
H_list.sort(key=lambda x: -x[2])

# Find which is Paley
sigma_paley = tuple(1 if k in S_paley else -1 for k in range(1, m + 1))
# Wait, need to be more careful:
# sigma_k = +1 if chord type k goes clockwise, i.e., k in S
sigma_paley_11 = tuple(1 if k in S_paley else -1 for k in range(1, m + 1))

print(f"\nPaley QR = {sorted(S_paley)}")
print(f"Paley sigma = {sigma_paley_11}")
print(f"H(Paley) = {H_vals_11[sigma_paley_11]}")

print(f"\nTop 10 by H:")
for rank, (sigma, S, H) in enumerate(H_list[:10], 1):
    note = ""
    if S == S_paley:
        note = " *** PALEY ***"
    elif S == frozenset(range(1, p)) - S_paley:
        note = " *** PALEY compl ***"
    elif S == frozenset(range(1, m + 1)):
        note = " *** INTERVAL ***"
    elif S == frozenset(range(m + 1, p)):
        note = " *** INTERVAL compl ***"
    print(f"  #{rank}: H={H}, sigma={sigma}, S={sorted(S)}{note}")

print(f"\nBottom 5:")
for rank, (sigma, S, H) in enumerate(H_list[-5:], len(H_list) - 4):
    print(f"  #{rank}: H={H}, sigma={sigma}, S={sorted(S)}")

# Unique H values
unique_H = sorted(set(h for _, _, h in H_list), reverse=True)
print(f"\n{len(unique_H)} distinct H values:")
for h in unique_H:
    count = sum(1 for _, _, hh in H_list if hh == h)
    print(f"  H={h}: {count} orientations")

# Walsh-Fourier expansion
print(f"\nWalsh-Fourier expansion (nonzero coefficients):")
coeffs = walsh_transform(H_vals_11, m)
for S in sorted(coeffs.keys(), key=lambda s: (len(s), s)):
    c = coeffs[S]
    print(f"  hat{{H}}{S} = {c:.4f}  (degree {len(S)})")

# Check: are odd-degree coefficients zero?
max_odd = max((abs(c) for S, c in coeffs.items() if len(S) % 2 == 1), default=0)
print(f"\nMax odd-degree coefficient: {max_odd:.6f} (should be ~0)")

# Degree structure
print(f"\nDegree structure:")
for d in range(m + 1):
    deg_coeffs = {S: c for S, c in coeffs.items() if len(S) == d}
    if deg_coeffs:
        total_energy = sum(c**2 for c in deg_coeffs.values())
        print(f"  degree {d}: {len(deg_coeffs)} terms, total L2 energy = {total_energy:.4f}")
        if d <= 2:
            for S, c in sorted(deg_coeffs.items()):
                print(f"    {S}: {c:.4f}")


# ============================================================
# SECTION 2: UNDERSTANDING WHY PALEY MAXIMIZES AT p=7,11
# ============================================================

print(f"\n{'=' * 70}")
print("WHY PALEY MAXIMIZES: INTERACTION MATRIX ANALYSIS")
print("=" * 70)

# From the Walsh expansion, degree-2 coefficients form a matrix J:
#   hat{H}({i,j}) = J[i,j]
# H(sigma) = hat{H}({}) + sum_{i<j} J[i,j] * sigma_i * sigma_j + ...
# To maximize, we want sigma to align with the positive eigenspace of J.

# p=7 first (from previous computation)
print(f"\np=7: Interaction matrix J[i,j] = hat{{H}}({{i,j}}):")
# Recall: hat{H}((0,1)) = 3.5, hat{H}((0,2)) = -3.5, hat{H}((1,2)) = -3.5
J_7 = np.array([[0, 3.5, -3.5],
                 [3.5, 0, -3.5],
                 [-3.5, -3.5, 0]])
print(J_7)
eigvals_7, eigvecs_7 = np.linalg.eigh(J_7)
print(f"Eigenvalues: {eigvals_7}")
print(f"Top eigenvector: {eigvecs_7[:, -1]}")

# The maximizer sigma should align with the top eigenvector
print(f"\nPaley sigma: [1, 1, -1]")
print(f"Top eigenvector sign pattern: {np.sign(eigvecs_7[:, -1])}")

# p=11
print(f"\np=11: Interaction matrix J[i,j] = hat{{H}}({{i,j}}):")
J_11 = np.zeros((m, m))
for S, c in coeffs.items():
    if len(S) == 2:
        i, j = S
        J_11[i, j] = c
        J_11[j, i] = c
print(J_11)
eigvals_11, eigvecs_11 = np.linalg.eigh(J_11)
print(f"Eigenvalues: {eigvals_11}")
print(f"Top eigenvector: {eigvecs_11[:, -1]}")
print(f"Top eigvec sign: {np.sign(eigvecs_11[:, -1])}")
print(f"Paley sigma:     {list(sigma_paley_11)}")

# Does the Paley sigma maximize the quadratic form?
sigma_P = np.array(sigma_paley_11)
quad_P = sigma_P @ J_11 @ sigma_P
print(f"\nQuadratic form at Paley: sigma^T J sigma = {quad_P:.4f}")

# Check the interval
sigma_I = np.ones(m)
quad_I = sigma_I @ J_11 @ sigma_I
print(f"Quadratic form at Interval: sigma^T J sigma = {quad_I:.4f}")

# Check all orientations
print(f"\nQuadratic form values (sigma^T J sigma) for all 2^{m} orientations:")
for bits in range(2**m):
    sigma = np.array([1 if (bits >> i) & 1 else -1 for i in range(m)])
    q = sigma @ J_11 @ sigma
    sigma_t = tuple(sigma.astype(int))
    H = H_vals_11[sigma_t]
    # Is this the best quadratic?
    print(f"  {sigma_t}: Q={q:.2f}, H={H}")


# ============================================================
# SECTION 3: DEGREE-4 AND HIGHER TERMS
# ============================================================

print(f"\n{'=' * 70}")
print("HIGHER-ORDER INTERACTIONS")
print("=" * 70)

# At p=7: only degree 0 and 2, no degree 4+ (m=3 is too small for degree 4)
# Wait, there's no degree 4 since max is 3! But we could have degree 2.
# Actually for m=3, the max subset size is 3, and we showed only degree 0 and 2.

# At p=11: we can have degree 0, 2, 4 (even degrees up to m=5, but m=5 is odd)
# Even subsets of {0,1,2,3,4}: sizes 0, 2, 4
print(f"\np=11: Higher-order term analysis")
print(f"  Total terms: {len(coeffs)}")
for d in [0, 2, 4]:
    terms = [(S, c) for S, c in coeffs.items() if len(S) == d]
    print(f"  Degree {d}: {len(terms)} terms")
    if d == 4:
        for S, c in sorted(terms):
            print(f"    {S}: {c:.4f}")

# The KEY question: does the degree-4 contribution flip sign at p=19?
# At p=7,11: the degree-2 interaction matrix favors Paley.
# At p=19: perhaps the degree-4 terms overwhelm the degree-2?


# ============================================================
# SECTION 4: p=13 — FILL IN THE GAP
# ============================================================

print(f"\n{'=' * 70}")
print("p=13: FILLING THE CROSSOVER GAP")
print("=" * 70)

p = 13
m = 6

# p=13 is 1 mod 4, so no Paley tournament. But we can check interval vs others.
S_interval = frozenset(range(1, m + 1))

print(f"p={p}, m={m} (p mod 4 = {p % 4})")
print(f"No Paley tournament at p=13 (p = 1 mod 4)")
print(f"Interval S = {sorted(S_interval)}")

# But wait — with p=13 (2^13 DP table), this is already very memory-intensive.
# Held-Karp needs 2^13 * 13 entries. That's 8192 * 13 = 106496 entries. Fine.
print(f"Held-Karp table size: 2^{p} * {p} = {2**p * p} entries")

A_int = adjacency_matrix(S_interval, p)
t0 = time.time()
H_int = count_hp_fast(A_int)
elapsed = time.time() - t0
print(f"\nInterval H = {H_int} ({elapsed:.1f}s)")

# Try a few other connection sets at p=13
# Even residues
S_even = frozenset(range(2, p, 2))  # {2,4,6,8,10,12}
valid_even = all((p - j) % p not in S_even for j in S_even)
print(f"\nEven residues S = {sorted(S_even)}, valid? {valid_even}")
if valid_even and len(S_even) == m:
    A_even = adjacency_matrix(S_even, p)
    H_even = count_hp_fast(A_even)
    print(f"H(even) = {H_even}")

# QR set (symmetric since p=1 mod 4, so not a valid tournament)
S_qr = frozenset(j for j in range(1, p) if is_qr(j, p))
valid_qr = all((p - j) % p not in S_qr for j in S_qr)
print(f"\nQR set = {sorted(S_qr)}, valid? {valid_qr}")

# Try some random valid connection sets
import random
random.seed(42)
best_H = H_int
best_S = S_interval
results_13 = [("Interval", S_interval, H_int)]

for trial in range(20):
    elems = list(range(1, p))
    random.shuffle(elems)
    S = set()
    for j in elems:
        if j not in S and (p - j) % p not in S:
            S.add(j)
            if len(S) == m:
                break
    if len(S) != m:
        continue
    S = frozenset(S)
    A = adjacency_matrix(S, p)
    H = count_hp_fast(A)
    results_13.append((f"Random-{trial}", S, H))
    if H > best_H:
        best_H = H
        best_S = S

# Also try interval complement
S_int_comp = frozenset(range(m + 1, p))
valid_comp = all((p - j) % p not in S_int_comp for j in S_int_comp)
if valid_comp and len(S_int_comp) == m:
    A = adjacency_matrix(S_int_comp, p)
    H = count_hp_fast(A)
    results_13.append(("Interval_comp", S_int_comp, H))

# Sort by H
results_13.sort(key=lambda x: -x[2])
print(f"\nResults at p=13:")
for name, S, H in results_13[:10]:
    marker = " ***" if H == best_H else ""
    print(f"  {name:>15}: H={H:>15}, S={sorted(S)}{marker}")

print(f"\nBest H = {best_H} for S = {sorted(best_S)}")


# ============================================================
# SECTION 5: WINDING NUMBER ANALYSIS AT p=11
# ============================================================

print(f"\n{'=' * 70}")
print("WINDING NUMBER AT p=11: PALEY vs INTERVAL")
print("=" * 70)

p = 11
m = 5
S_paley_11 = frozenset(j for j in range(1, p) if is_qr(j, p))
S_interval_11 = frozenset(range(1, m + 1))

# We can't enumerate all HPs at p=11 (too many: ~95000).
# But we can sample or compute the winding distribution from eigenvalues.

# Instead, compute the NET CLOCKWISE STEPS statistic.
# For a circulant T(S): each step in a HP uses some d in S or p-d (for d not in S).
# The "winding" depends on how many steps are "short clockwise" vs "short counterclockwise".

# Actually, for the orientation cube, the winding information is encoded
# in the degree-2 Walsh coefficients.

# Let's instead analyze the MOMENT GENERATING FUNCTION of H.
# The Walsh expansion gives:
#   H(sigma) = sum_S hat{H}(S) * chi_S(sigma)
# This is a POLYNOMIAL in the sigma variables.

# The polynomial at p=11:
print(f"H polynomial on orientation cube at p=11:")
print(f"  H(sigma) = {coeffs.get((), 0):.2f}")
for S, c in sorted(coeffs.items(), key=lambda x: (len(x[0]), x[0])):
    if len(S) == 0:
        continue
    term = " * ".join(f"s{k+1}" for k in S)
    sign = "+" if c > 0 else "-"
    print(f"           {sign} {abs(c):.2f} * {term}")

# Which sigma maximizes this polynomial?
# The degree-2 part is a quadratic form: sigma^T J sigma.
# Paley maximizes the quadratic form (we showed above).
# But do the degree-4 terms change the answer?

print(f"\nDegree-2 form at Paley: {quad_P:.2f}")
print(f"Degree-4 contribution at Paley:", end=" ")
deg4_P = sum(c * np.prod([sigma_paley_11[k] for k in S])
             for S, c in coeffs.items() if len(S) == 4)
print(f"{deg4_P:.2f}")

deg4_I = sum(c for S, c in coeffs.items() if len(S) == 4)
quad_I_11 = sum(c for S, c in coeffs.items() if len(S) == 2)
print(f"Degree-2 form at Interval: {quad_I_11:.2f}")
print(f"Degree-4 contribution at Interval: {deg4_I:.2f}")

# Decomposition: H = H0 + H2 + H4
# where H0 = hat{H}({}), H2 = quad form, H4 = quartic
H0 = coeffs.get((), 0)
H_paley = H_vals_11[sigma_paley_11]
H_interval = H_vals_11[tuple(1 for _ in range(m))]

print(f"\nDecomposition:")
print(f"  H0 (mean) = {H0:.2f}")
print(f"  Paley:    H0 + H2 + H4 = {H0:.2f} + {quad_P:.2f} + {deg4_P:.2f} = {H0 + quad_P + deg4_P:.2f} (actual: {H_paley})")
print(f"  Interval: H0 + H2 + H4 = {H0:.2f} + {quad_I_11:.2f} + {deg4_I:.2f} = {H0 + quad_I_11 + deg4_I:.2f} (actual: {H_interval})")

# The KEY insight: at what point do the degree-4 terms favor the interval?


# ============================================================
# SECTION 6: PALEY AS EIGENVECTOR OF INTERACTION MATRIX
# ============================================================

print(f"\n{'=' * 70}")
print("PALEY SIGMA AS EIGENVECTOR OF J")
print("=" * 70)

print(f"\np=11: Is Paley sigma an eigenvector of J?")
sigma_P_arr = np.array(sigma_paley_11, dtype=float)
Js = J_11 @ sigma_P_arr
print(f"  sigma_P = {sigma_P_arr}")
print(f"  J * sigma_P = {Js}")
ratio = Js / sigma_P_arr
print(f"  J*sigma / sigma = {ratio}")
is_eigvec = np.std(ratio) < 0.01
print(f"  Is eigenvector? {is_eigvec} (std of ratios = {np.std(ratio):.4f})")

# Check: is J circulant or has special structure?
print(f"\n  J matrix structure:")
print(J_11)
# J[i,j] = hat{H}({i,j}) for the orientation cube
# This is NOT the same as the tournament adjacency matrix!
# It's the 2-body interaction between chord orientations.

# What's the relationship between J[i,j] and the chord geometry?
print(f"\n  J[i,j] vs chord distance |i-j|:")
for i in range(m):
    for j in range(i+1, m):
        # Chord types are i+1 and j+1
        # Their "distance" on the orientation cube
        print(f"    chords ({i+1},{j+1}): J = {J_11[i,j]:.4f}, |i-j| = {j-i}")


# ============================================================
# SECTION 7: MULTIPLICATIVE STRUCTURE AND D_{2p}
# ============================================================

print(f"\n{'=' * 70}")
print("MULTIPLICATIVE AUTOMORPHISMS AND ORIENTATION CUBE")
print("=" * 70)

print("""
For Paley T_p, the automorphism group Aut(T_p) acts on the orientation cube.
Multiplying the connection set by a QR element a:
  S -> a*S mod p
maps T_p to an ISOMORPHIC tournament (same H).

This means: if sigma is the orientation of T_p, then a*sigma
(rearranged by the multiplication) also gives the same H.

For p=7: Aut(T_7) has order 7*3 = 21.
The orientation cube has 8 points, so the orbit of Paley sigma
can have at most 21/gcd points. But the cube is tiny (8 points)
so the effective action is small.

KEY: The multiplicative action PERMUTES chord types.
Multiplying by a QR element a maps chord type k to chord type a*k mod p.
This is a PERMUTATION of {1,...,m} (since a*QR = QR).
""")

for p in [7, 11]:
    m = (p - 1) // 2
    S_paley = frozenset(j for j in range(1, p) if is_qr(j, p))
    sigma_P = [1 if k in S_paley else -1 for k in range(1, m + 1)]

    print(f"\np={p}: Multiplicative permutations of chord types")
    # Find a primitive root
    for g in range(2, p):
        if len(set(pow(g, k, p) for k in range(1, p))) == p - 1:
            break
    print(f"  Primitive root: g={g}")

    # QR elements
    qr_elts = sorted(j for j in range(1, p) if is_qr(j, p))
    print(f"  QR elements: {qr_elts}")

    # Action of each QR element on chord types
    for a in qr_elts[:5]:
        perm = [(a * k % p) if (a * k % p) <= m else p - (a * k % p)
                for k in range(1, m + 1)]
        # Also track the sign: if a*k mod p > m, chord orientation flips
        signs = [1 if (a * k % p) <= m else -1 for k in range(1, m + 1)]
        print(f"  a={a}: chord {list(range(1,m+1))} -> {perm}, signs={signs}")
