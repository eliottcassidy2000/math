#!/usr/bin/env python3
"""
Verify the identity: u_T(m) = m^n * I(Omega(T), 2/m^2)

where u_T(m) = ps_1(U_T)(m) is the principal specialization at m,
and I(Omega(T), x) = sum_k alpha_k * x^k is the independence polynomial.

This means:
  u_T(m) = m^n * sum_k alpha_k * (2/m^2)^k
         = m^n * sum_k alpha_k * 2^k * m^{-2k}
         = sum_k alpha_k * 2^k * m^{n-2k}

Evaluating at m=1: u_T(1) = sum_k alpha_k * 2^k = I(Omega, 2) = H(T)  ✓

Special values:
  u_T(sqrt(2)) = (sqrt(2))^n * I(Omega, 1)
               = 2^{n/2} * I(Omega, 1)
               = 2^{n/2} * (number of independent sets in Omega)

  u_T(1) = H(T) = I(Omega, 2)

  u_T(sqrt(2/3)) = (sqrt(2/3))^n * I(Omega, 3)

  u_T(i) = i^n * I(Omega, -2) where i = sqrt(-1)
  Since n is odd: i^n = i^{2m+1} = i * (-1)^m = i*(-1)^m
  So u_T(i) = i*(-1)^m * I(Omega, -2)

This is interesting because I(Omega, -2) relates to the chromatic polynomial!
For any graph G: I(G, -x) relates to the chromatic polynomial via Whitney's theorem.

opus-2026-03-07-S39
"""
from itertools import permutations, combinations
from collections import defaultdict
import numpy as np


def make_tournament(n, bits):
    A = [[0]*n for _ in range(n)]
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    for k, (i,j) in enumerate(edges):
        if bits & (1 << k):
            A[j][i] = 1
        else:
            A[i][j] = 1
    return A


def compute_uT_coeffs(n, A):
    """u_T(m) coefficients via permutation sum."""
    edge_set = set()
    opp_set = set()
    for i in range(n):
        for j in range(n):
            if i != j:
                if A[i][j]:
                    edge_set.add((i, j))
                else:
                    opp_set.add((i, j))

    coeffs = defaultdict(int)
    for sigma in permutations(range(n)):
        visited = [False] * n
        cycles = []
        for s in range(n):
            if visited[s]:
                continue
            cyc = []
            c = s
            while not visited[c]:
                visited[c] = True
                cyc.append(c)
                c = sigma[c]
            cycles.append(tuple(cyc))

        valid = True
        phi = 0
        has_even = False
        for cyc in cycles:
            if len(cyc) == 1:
                continue
            if len(cyc) % 2 == 0:
                has_even = True
                break
            is_T = all((cyc[i], cyc[(i+1)%len(cyc)]) in edge_set for i in range(len(cyc)))
            is_Top = all((cyc[i], cyc[(i+1)%len(cyc)]) in opp_set for i in range(len(cyc)))
            if not is_T and not is_Top:
                valid = False
                break
            if is_T:
                phi += len(cyc) - 1
        if not valid or has_even:
            continue
        sign = (-1) ** phi
        coeffs[len(cycles)] += sign
    return coeffs


def compute_alpha(n, A):
    """Compute alpha_k (independence numbers of Omega(T))."""
    edge_set = set()
    for i in range(n):
        for j in range(n):
            if i != j and A[i][j]:
                edge_set.add((i, j))

    # Find all directed odd cycles
    all_cycles = []
    for length in range(3, n+1, 2):
        for verts in combinations(range(n), length):
            seen = set()
            for p in permutations(verts):
                if all((p[i], p[(i+1)%length]) in edge_set for i in range(length)):
                    min_idx = list(p).index(min(p))
                    canon = tuple(list(p)[min_idx:] + list(p)[:min_idx])
                    if canon not in seen:
                        seen.add(canon)
                        all_cycles.append((frozenset(verts), canon))

    # Build conflict graph adjacency
    nc = len(all_cycles)
    adj = [[False]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if all_cycles[i][0] & all_cycles[j][0]:  # share vertex
                adj[i][j] = adj[j][i] = True

    # Count independent sets by size
    alpha = defaultdict(int)
    alpha[0] = 1

    # Enumerate all independent sets via BFS
    def count_indep_sets(cycles_adj, nc):
        counts = defaultdict(int)
        counts[0] = 1
        # Bitmask enumeration for small nc
        if nc <= 20:
            for mask in range(1, 1 << nc):
                bits_list = [i for i in range(nc) if mask & (1 << i)]
                is_indep = True
                for a in range(len(bits_list)):
                    for b in range(a+1, len(bits_list)):
                        if cycles_adj[bits_list[a]][bits_list[b]]:
                            is_indep = False
                            break
                    if not is_indep:
                        break
                if is_indep:
                    counts[len(bits_list)] += 1
            return counts
        return counts

    alpha = count_indep_sets(adj, nc)
    return alpha


# === Verify identity at n=5 ===
print("=== Verification: u_T(m) = m^n * I(Omega, 2/m^2) at n=5 ===")
n = 5
test_cases = [0, 42, 100, 341, 682, 1023]
for bits in test_cases:
    A = make_tournament(n, bits)
    coeffs = compute_uT_coeffs(n, A)
    alpha = compute_alpha(n, A)

    # Build u_T(m) from coefficients
    def uT(m):
        return sum(c * m**k for k, c in coeffs.items())

    # Build I(Omega, x) from alpha
    def I_omega(x):
        return sum(a * x**k for k, a in alpha.items())

    # Check identity at several m values
    m_vals = [1.0, 2.0, 0.5, 1.41421356, 3.0]
    match = True
    for m in m_vals:
        lhs = uT(m)
        rhs = m**n * I_omega(2.0 / m**2)
        if abs(lhs - rhs) > 1e-6:
            match = False
            print(f"  MISMATCH bits={bits}, m={m}: u_T={lhs}, m^n*I(2/m^2)={rhs}")

    H = uT(1.0)
    I_at_2 = I_omega(2.0)
    if match:
        print(f"  bits={bits}: H={int(H)}, I(Omega,2)={int(I_at_2)}, identity VERIFIED at 5 m-values")


# === Special value analysis ===
print("\n=== Special values of u_T(m) ===")

# Paley T_7
A = [[0]*7 for _ in range(7)]
QR = {1, 2, 4}
for i in range(7):
    for j in range(7):
        if i != j and (j - i) % 7 in QR:
            A[i][j] = 1

coeffs = compute_uT_coeffs(7, A)
alpha = compute_alpha(7, A)

def uT(m):
    return sum(c * m**k for k, c in coeffs.items())

def I_omega(x):
    return sum(a * x**k for k, a in alpha.items())

print(f"Paley T_7:")
print(f"  alpha = {dict(sorted(alpha.items()))}")
print(f"  u_T coefficients = {dict(sorted(coeffs.items()))}")

print(f"\n  u_T(1) = {uT(1.0):.0f} = H(T) = I(Omega, 2)")
print(f"  u_T(-1) = {uT(-1.0):.0f} = -H(T)")
print(f"  u_T(sqrt(2)) = {uT(2**0.5):.4f} = 2^(7/2) * I(Omega, 1) = {2**3.5 * I_omega(1.0):.4f}")
print(f"  I(Omega, 1) = {I_omega(1.0):.0f} = total independent sets in Omega")
print(f"  I(Omega, -1) = {I_omega(-1.0):.0f} = alternating sum of alpha_k")
print(f"  I(Omega, -2) = {I_omega(-2.0):.0f}")

# u_T(i) where i = sqrt(-1)
import cmath
m = 1j
uT_at_i = sum(c * m**k for k, c in coeffs.items())
m_n = m**7  # i^7 = i^4 * i^3 = 1 * (-i) = -i
I_at_neg2 = I_omega(-2.0)
print(f"\n  u_T(i) = {uT_at_i}")
print(f"  i^7 * I(Omega, -2) = {m_n * I_at_neg2}")
print(f"  Match: {abs(uT_at_i - m_n * I_at_neg2) < 1e-10}")

# u_T(sqrt(2)*i) = (sqrt(2)*i)^7 * I(Omega, -1)
# (sqrt(2)*i)^7 = 2^(7/2) * i^7 = 2^(7/2) * (-i) = -i*2^(7/2)
m2 = 2**0.5 * 1j
uT_at_sqrt2i = sum(c * m2**k for k, c in coeffs.items())
pred = (m2**7) * I_omega(-1.0)
print(f"\n  u_T(sqrt(2)*i) = {uT_at_sqrt2i}")
print(f"  (sqrt(2)*i)^7 * I(Omega,-1) = {pred}")
print(f"  Match: {abs(uT_at_sqrt2i - pred) < 1e-10}")


# === u_T at m=2: what does it count? ===
print(f"\n  u_T(2) = {uT(2.0):.0f} = 2^7 * I(Omega, 1/2)")
print(f"  I(Omega, 1/2) = {I_omega(0.5):.4f}")
print(f"  This is the \"half-fugacity\" partition function of the hard-core model!")

# === Chromatic connection ===
print(f"\n=== Chromatic polynomial connection ===")
print(f"I(Omega, -1) = {I_omega(-1.0):.0f}")
print(f"This is related to the chromatic polynomial: for any graph G,")
print(f"(-1)^|V| * chi(G, -1) counts acyclic orientations (Stanley's theorem).")
print(f"But Omega is the conflict graph, not the tournament itself.")
print(f"For Omega(T_7): I(Omega, -1) = {I_omega(-1.0):.0f}")

# What about sum_k (-1)^k alpha_k ?
alt_sum = sum((-1)**k * alpha.get(k, 0) for k in range(max(alpha.keys())+1))
print(f"Alternating sum = {alt_sum}")
print(f"This counts (even indep sets) - (odd indep sets) = {alt_sum}")
