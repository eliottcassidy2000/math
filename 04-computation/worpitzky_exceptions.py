#!/usr/bin/env python3
"""
INVESTIGATE EXCEPTIONS: F(T,x) with non-real roots or log-concavity failure.

At n=5: 48 tournaments have F with non-real roots, 4 have log-concavity failure.
What's special about them?

Also: explore F(T,x) modulo 2 and its connection to OCF parity.
"""
from itertools import permutations
from math import comb, factorial
import numpy as np
from collections import defaultdict, Counter

def all_tournaments(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def forward_edge_poly(A, n):
    F = [0] * n
    for P in permutations(range(n)):
        if not all(A[P[i]][P[i+1]] for i in range(n-1)):
            continue
        asc = sum(1 for i in range(n-1) if P[i] < P[i+1])
        F[asc] += 1
    return F

def count_3cycles(A, n):
    t3 = 0
    for i in range(n):
        for j in range(i+1,n):
            for k in range(j+1,n):
                if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[k][j] and A[i][k]):
                    t3 += 1
    return t3

# ===== Find the exceptions at n=5 =====
print("=" * 70)
print("EXCEPTIONS AT n=5: NON-REAL ROOTS")
print("=" * 70)

n = 5
d = n - 1
non_real = []
not_lc = []

for A in all_tournaments(n):
    F = forward_edge_poly(A, n)
    H = sum(F)
    t3 = count_3cycles(A, n)
    scores = tuple(sorted([sum(A[i]) for i in range(n)]))

    # Check roots
    coeffs = list(reversed(F))
    if coeffs[0] == 0:
        # Leading coeff is 0, reduce degree
        while coeffs and coeffs[0] == 0:
            coeffs.pop(0)
    if len(coeffs) > 1:
        roots = np.roots(coeffs)
        has_complex = any(abs(r.imag) > 1e-8 for r in roots)
    else:
        has_complex = False

    if has_complex:
        non_real.append((F, H, t3, scores, roots))

    # Check log-concavity
    is_lc = True
    for k in range(1, d):
        if F[k]**2 < F[k-1] * F[k+1]:
            is_lc = False
            break
    if not is_lc:
        not_lc.append((F, H, t3, scores))

print(f"\n{len(non_real)} tournaments with non-real roots:")
# Group by F
by_F = defaultdict(list)
for F, H, t3, scores, roots in non_real:
    by_F[tuple(F)].append((H, t3, scores))

for F_tuple, data in sorted(by_F.items()):
    F = list(F_tuple)
    H = data[0][0]
    # Compute roots for display
    coeffs = list(reversed(F))
    while coeffs and coeffs[0] == 0:
        coeffs.pop(0)
    roots = np.roots(coeffs)
    real_parts = sorted([r.real for r in roots])
    imag_parts = sorted([abs(r.imag) for r in roots], reverse=True)

    t3_vals = sorted(set(d[1] for d in data))
    score_vals = sorted(set(d[2] for d in data))
    print(f"  F={F}, H={H}, t3={t3_vals}, count={len(data)}")
    print(f"    roots: {[f'{r.real:.4f}+{r.imag:.4f}i' for r in roots]}")

print(f"\n\n{len(not_lc)} tournaments failing log-concavity:")
by_F_lc = defaultdict(int)
for F, H, t3, scores in not_lc:
    by_F_lc[tuple(F)] += 1

for F_tuple, count in sorted(by_F_lc.items()):
    F = list(F_tuple)
    print(f"  F={F}, H={sum(F)}, count={count}")
    # Which k fails?
    for k in range(1, d):
        if F[k]**2 < F[k-1] * F[k+1]:
            print(f"    Fails at k={k}: F[{k}]^2={F[k]**2} < F[{k-1}]*F[{k+1}]={F[k-1]*F[k+1]}")

# ===== SECTION 2: F modulo 2 and OCF =====
print("\n\n" + "=" * 70)
print("F(T,x) MODULO 2 AND OCF PARITY")
print("=" * 70)
print("OCF: H(T) ≡ 1 (mod 2) always (Rédei's theorem)")
print("F(T,1) = H ≡ 1 (mod 2)")
print("What are the parities of individual F_k?")

for n in [3, 4, 5]:
    d = n - 1
    print(f"\nn={n}:")

    parity_counts = Counter()
    F_mod2_counts = Counter()

    for A in all_tournaments(n):
        F = forward_edge_poly(A, n)
        H = sum(F)
        F_mod2 = tuple(f % 2 for f in F)
        F_mod2_counts[F_mod2] += 1
        parity_counts[H % 2] += 1

    print(f"  H parity: {dict(parity_counts)} (all odd, as expected by Rédei)")
    print(f"  F mod 2 patterns:")
    for pattern in sorted(F_mod2_counts.keys()):
        print(f"    {pattern}: {F_mod2_counts[pattern]}")

    # Sum of F_k must be odd. So sum(F_k mod 2) must be odd (mod 2).
    # i.e., #{k : F_k is odd} is odd.
    odd_count_dist = Counter()
    for pattern in F_mod2_counts:
        num_odd = sum(pattern)
        odd_count_dist[num_odd] += F_mod2_counts[pattern]
    print(f"  #{'{'}k: F_k odd{'}'} distribution: {dict(sorted(odd_count_dist.items()))}")

# ===== SECTION 3: F(T,2) and OCF =====
print("\n\n" + "=" * 70)
print("F(T,2) vs I(Omega,2) = H")
print("=" * 70)
print("F(T,2) = sum F_k * 2^k. Is there a simple relation to H?")
print("We know F(T,1) = H and F(T,2) = H + sum F_k(2^k - 1)")

for n in [4, 5]:
    d = n - 1
    print(f"\nn={n}:")

    data = []
    for A in all_tournaments(n):
        F = forward_edge_poly(A, n)
        H = sum(F)
        F2 = sum(F[k] * 2**k for k in range(n))
        t3 = count_3cycles(A, n)
        data.append((H, F2, t3, F))

    # F(2) = H + sum_{k>0} F_k(2^k - 1) + F_0(2^0 - 1) = H + sum F_k(2^k - 1)
    # = H + F_1 + 3*F_2 + 7*F_3 + 15*F_4 + ...
    # = 2*H - sum F_k + F_1 + 3*F_2 + 7*F_3 + ... - F_0 - F_1 - ... - F_d
    # Hmm, simpler: F(2) = sum F_k 2^k, H = sum F_k
    # F(2) - H = sum F_k (2^k - 1) = F_1 + 3F_2 + 7F_3 + 15F_4
    # This is determined by the "weighted excess" beyond k=0.

    # Check: is F(2) - H related to anything?
    excess = [(F2 - H, H, t3) for H, F2, t3, F in data]
    by_H_ex = defaultdict(set)
    for ex, H, t3 in excess:
        by_H_ex[H].add(ex)

    print(f"  F(2) - H distribution by H:")
    for H_val in sorted(by_H_ex.keys()):
        vals = sorted(by_H_ex[H_val])
        print(f"    H={H_val}: F(2)-H in {vals[:10]}{'...' if len(vals)>10 else ''}")

    # F(2) modulo small numbers
    F2_mod2 = Counter(F2 % 2 for _, F2, _, _ in data)
    F2_mod4 = Counter(F2 % 4 for _, F2, _, _ in data)
    print(f"\n  F(2) mod 2: {dict(sorted(F2_mod2.items()))}")
    print(f"  F(2) mod 4: {dict(sorted(F2_mod4.items()))}")

# ===== SECTION 4: The "shape" of F beyond H =====
print("\n\n" + "=" * 70)
print("F SHAPE: F/H AS PROBABILITY DISTRIBUTION ON ASCENTS")
print("=" * 70)
print("F_k/H = probability that a random HP has k ascents")

for n in [5]:
    d = n - 1
    print(f"\nn={n}:")

    # For each tournament, compute the normalized shape
    all_shapes = []
    for A in all_tournaments(n):
        F = forward_edge_poly(A, n)
        H = sum(F)
        shape = tuple(round(f/H, 6) for f in F)
        all_shapes.append((shape, H, F))

    # How many distinct shapes?
    shape_counts = Counter(s for s, _, _ in all_shapes)
    print(f"  Distinct shapes: {len(shape_counts)}")
    print(f"  Most common shapes:")
    for shape, count in shape_counts.most_common(10):
        # Find one example
        H = next(h for s, h, _ in all_shapes if s == shape)
        print(f"    {shape} (H={H}, count={count})")

    # Mean of ascent position
    print(f"\n  Mean ascent position E[k] for each tournament:")
    means = []
    for A in all_tournaments(n):
        F = forward_edge_poly(A, n)
        H = sum(F)
        mean_k = sum(k * F[k] for k in range(n)) / H
        means.append(mean_k)

    print(f"    E[k] range: [{min(means):.4f}, {max(means):.4f}]")
    print(f"    Overall mean E[k]: {np.mean(means):.4f}")
    print(f"    Expected for uniform (Eulerian): {sum(k*A_nk for k, A_nk in enumerate([1,26,66,26,1]))/120:.4f}")
    print(f"    = (n-1)/2 = {d/2:.4f}")

    # Variance of ascent position
    vars_k = []
    for A in all_tournaments(n):
        F = forward_edge_poly(A, n)
        H = sum(F)
        mean_k = sum(k * F[k] for k in range(n)) / H
        var_k = sum((k - mean_k)**2 * F[k] for k in range(n)) / H
        vars_k.append(var_k)

    print(f"    Var[k] range: [{min(vars_k):.4f}, {max(vars_k):.4f}]")
    print(f"    Overall mean Var[k]: {np.mean(vars_k):.4f}")

# ===== SECTION 5: F(T,x) and the transfer matrix eigenvalues =====
print("\n\n" + "=" * 70)
print("SPECTRAL CONNECTION: EIGENVALUES OF M_x")
print("=" * 70)
print("""
M_x[v,u] = A[v][u] * x^{[v<u]} is NOT a normal matrix.
Its eigenvalues don't directly give F(T,x).
But we can consider det(I - t*M_x) = "characteristic polynomial".

For the PERMANENT: perm(M_x) would give F(T,x) if M_x were appropriately
defined, but permanent ≠ determinant in general.

However, for small n, let's check if det(I - t*M_x) relates to F(T,x).
""")

for n in [3, 4]:
    print(f"\nn={n}:")
    for idx, A in enumerate(list(all_tournaments(n))[:5]):
        F = forward_edge_poly(A, n)
        H = sum(F)

        # Check permanent of M_x at x=1,2
        # Permanent at x=1 should be... not H (permanent counts perfect matchings,
        # not Hamiltonian paths).

        # Actually, F(T,x) is not the permanent of M_x. It's a path sum.
        # Let me instead compute the "principal minor sum" or the
        # "immanent" or something else.

        # For Hamiltonian PATHS, we sum over all orderings (v1,...,vn) of [n]
        # with A[vi][vi+1]=1 and weight x^{[vi<vi+1]}.
        # This is NOT the permanent of M_x (permanent sums over bijections sigma
        # of product M[i,sigma(i)]).

        # Let me instead look at what happens with the INCLUSION-EXCLUSION
        # version of HP counting.

        # Eigenvalues of the adjacency matrix (unweighted)
        A_mat = np.array(A, dtype=float)
        eigs = np.linalg.eigvals(A_mat)
        eigs_sorted = sorted(eigs, key=lambda x: -abs(x))

        print(f"  T#{idx}: F={F}, H={H}")
        print(f"    Eigenvalues of A: {[f'{e.real:.3f}+{e.imag:.3f}i' for e in eigs_sorted]}")

# ===== SECTION 6: Connection to even-indexed Worpitzky =====
print("\n\n" + "=" * 70)
print("F(T,x) EVALUATED AT x = -1: ALTERNATING HP COUNT")
print("=" * 70)
print("""
F(T,-1) = #{HPs with even ascents} - #{HPs with odd ascents}

By the complement relation: F(T,-1) = (-1)^{n-1} F(T^op, -1)

For n odd: F(T,-1) = F(T^op,-1)
For n even: F(T,-1) = -F(T^op,-1)

Since H(T) = H(T^op) (complement invariance), we get:
  H = F_0 + F_1 + ... + F_{d} (odd and even summed)
  F(-1) = F_0 - F_1 + F_2 - ... (alternating)

So: #{even-ascent HPs} = (H + F(-1))/2
    #{odd-ascent HPs} = (H - F(-1))/2

Is F(-1) always odd (like H)?
""")

for n in [3, 4, 5]:
    d = n - 1
    print(f"\nn={n}:")

    Fm1_parity = Counter()
    for A in all_tournaments(n):
        F = forward_edge_poly(A, n)
        Fm1 = sum(F[k] * (-1)**k for k in range(n))
        Fm1_parity[Fm1 % 2] += 1

    print(f"  F(-1) parity: {dict(Fm1_parity)}")
    # If F(-1) is always odd, then (H + F(-1))/2 is always an integer
    # (since H is always odd)

    # Check (H + F(-1))/2 and (H - F(-1))/2
    for A in list(all_tournaments(n))[:5]:
        F = forward_edge_poly(A, n)
        H = sum(F)
        Fm1 = sum(F[k] * (-1)**k for k in range(n))
        even_asc = (H + Fm1) // 2
        odd_asc = (H - Fm1) // 2
        print(f"  F={F}, H={H}, F(-1)={Fm1}, even_asc={even_asc}, odd_asc={odd_asc}")

print("\n\nDone.")
