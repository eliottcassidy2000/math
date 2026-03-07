#!/usr/bin/env python3
"""
Investigation: u_T(m) root structure.

kind-pasteur-S30 discovered that u_T(m) = ps_1(U_T)(m) is an ODD polynomial in m:
  u_T(m) = m * Q_T(m^2)

Questions:
1. Does Q_T(w) always have all real (negative) roots?
2. What is the degree pattern of Q_T?
3. How do the roots depend on the tournament structure?

u_T(m) = sum_{sigma in S(T,T^op), all odd cycles} (-1)^phi(sigma) * m^{parts(sigma)}
       = sum_k S_k * m^k  (only odd k appear)

where S_k = sum_{I: parts(I)=k} 2^k * I(T) is the f-level weighted sum.

Wait -- more precisely, ps_1(p_lambda)(m) = m^{len(lambda)}.
So u_T(m) = sum_{sigma} (-1)^phi(sigma) * m^{#cycles(sigma)}.

Since only all-odd-cycle perms survive (even cancel), and each odd cycle
has 2 directions, the sum simplifies to:
  u_T(m) = sum over independent sets S in Omega(T): 2^|S| * m^(n - sum_{C in S}(|C|-1))
         = m^n + sum_{k>=1} alpha_k * 2^k * m^{n-2k} * ... hmm

Actually: for a perm with t disjoint odd cycles of sizes l_1,...,l_t on
total of l_1+...+l_t vertices, the number of fixed points = n-(l_1+...+l_t),
so total cycle count = t + (n - sum l_i) = t + n - sum l_i.

And m^{#cycles} = m^{t + n - sum l_i}.

Let's just compute it directly for small n.

opus-2026-03-07-S39
"""
from itertools import permutations, combinations
from collections import defaultdict
import numpy as np
from numpy.polynomial import polynomial as P


def compute_uT_polynomial(n, A):
    """Compute u_T(m) = ps_1(U_T)(m) as a polynomial in m.

    For each valid permutation sigma in S(T, T^op) with all-odd nontrivial cycles:
      contribution = (-1)^phi(sigma) * m^{#cycles(sigma)}
    """
    edge_set = set()
    opp_set = set()
    for i in range(n):
        for j in range(n):
            if i != j:
                if A[i][j]:
                    edge_set.add((i, j))
                else:
                    opp_set.add((i, j))

    # coeffs[k] = coefficient of m^k
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

        # Check: all nontrivial cycles must be directed cycles of T or T^op
        valid = True
        phi = 0
        has_even = False
        for cyc in cycles:
            if len(cyc) == 1:
                continue
            if len(cyc) % 2 == 0:
                has_even = True
                break  # Even cycles cancel, skip
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
        num_cycles = len(cycles)  # includes fixed points
        coeffs[num_cycles] += sign

    return coeffs


def poly_from_coeffs(coeffs, max_deg):
    """Convert coefficient dict to numpy polynomial array."""
    p = [0] * (max_deg + 1)
    for k, v in coeffs.items():
        if k <= max_deg:
            p[k] = v
    return np.array(p, dtype=float)


def find_roots_of_Q(coeffs, n):
    """Given u_T(m) = sum c_k m^k (odd k only), extract Q_T(w) where u_T(m) = m*Q_T(m^2).

    If u_T(m) = c_1*m + c_3*m^3 + c_5*m^5 + ...
    Then Q_T(w) = c_1 + c_3*w + c_5*w^2 + ...
    """
    # Extract odd coefficients
    q_coeffs = []
    for k in range(1, n + 1, 2):
        q_coeffs.append(coeffs.get(k, 0))

    if len(q_coeffs) <= 1:
        return np.array([]), q_coeffs

    # Q_T(w) = q_coeffs[0] + q_coeffs[1]*w + q_coeffs[2]*w^2 + ...
    q_array = np.array(q_coeffs, dtype=float)

    # Find roots of Q_T
    if len(q_array) > 1 and q_array[-1] != 0:
        # Use numpy roots (expects highest degree first)
        roots = np.roots(q_array[::-1])
        return roots, q_coeffs
    else:
        return np.array([]), q_coeffs


def make_tournament(n, edges_bits):
    """Create tournament adjacency matrix from bit encoding."""
    A = [[0]*n for _ in range(n)]
    edge_list = [(i,j) for i in range(n) for j in range(i+1,n)]
    for k, (i,j) in enumerate(edge_list):
        if edges_bits & (1 << k):
            A[j][i] = 1
        else:
            A[i][j] = 1
    return A


# === n=3 ===
print("=== n=3 ===")
for bits in range(1 << 3):
    A = make_tournament(3, bits)
    coeffs = compute_uT_polynomial(3, A)
    roots, q = find_roots_of_Q(coeffs, 3)
    print(f"bits={bits}: u_T(m) coeffs={dict(coeffs)}, Q_T={q}, roots={roots}")

# === n=5 exhaustive ===
print("\n=== n=5 exhaustive ===")
n = 5
edge_count = n*(n-1)//2  # 10
all_real_count = 0
complex_count = 0
total = 0
examples_complex = []

for bits in range(1 << edge_count):
    A = make_tournament(n, bits)
    coeffs = compute_uT_polynomial(n, A)
    roots, q = find_roots_of_Q(coeffs, n)
    total += 1

    if len(roots) > 0:
        if all(abs(r.imag) < 1e-10 for r in roots):
            all_real_count += 1
        else:
            complex_count += 1
            if len(examples_complex) < 3:
                examples_complex.append((bits, q, roots))

print(f"Total tournaments: {total}")
print(f"Q_T has all real roots: {all_real_count}")
print(f"Q_T has complex roots: {complex_count}")
if examples_complex:
    for bits, q, roots in examples_complex:
        print(f"  Example: bits={bits}, Q={q}, roots={roots}")

# Print a few examples
print("\nSample u_T polynomials at n=5:")
for bits in [0, 1, 31, 100, 500]:
    A = make_tournament(5, bits)
    coeffs = compute_uT_polynomial(5, A)
    roots, q = find_roots_of_Q(coeffs, 5)

    # Also compute H(T)
    H = sum(v for v in coeffs.values())  # u_T(1) = H(T)... wait, not quite
    # Actually u_T(1) should give H(T)
    uT_at_1 = sum(coeffs.get(k, 0) for k in range(n+1))

    poly_str = " + ".join(f"{coeffs[k]}*m^{k}" for k in sorted(coeffs.keys()) if coeffs[k] != 0)
    print(f"  bits={bits}: u_T(m) = {poly_str}, u_T(1)={uT_at_1}, Q roots={[f'{r:.4f}' for r in roots]}")

# === n=7 sampling ===
print("\n=== n=7 sampling ===")
import random
n = 7
edge_count = 21
all_real = 0
has_complex = 0
sample_size = 5000
all_neg_real = 0

random.seed(42)
for _ in range(sample_size):
    bits = random.randint(0, (1 << edge_count) - 1)
    A = make_tournament(n, bits)
    coeffs = compute_uT_polynomial(n, A)
    roots, q = find_roots_of_Q(coeffs, n)

    if len(roots) > 0:
        if all(abs(r.imag) < 1e-10 for r in roots):
            all_real += 1
            if all(r.real <= 1e-10 for r in roots if abs(r.imag) < 1e-10):
                all_neg_real += 1
        else:
            has_complex += 1

print(f"Sampled {sample_size} tournaments at n=7")
print(f"Q_T has all real roots: {all_real}/{sample_size}")
print(f"Q_T has complex roots: {has_complex}/{sample_size}")
print(f"Q_T has all NEGATIVE real roots: {all_neg_real}/{sample_size}")

# Show some Q_T polynomials for named tournaments
print("\n=== Named tournaments ===")

# Transitive n=5
A = [[0]*5 for _ in range(5)]
for i in range(5):
    for j in range(i+1, 5):
        A[i][j] = 1
coeffs = compute_uT_polynomial(5, A)
roots, q = find_roots_of_Q(coeffs, 5)
print(f"Transitive T_5: Q={q}, roots={roots}")

# C_5 (Paley T_5)
A = [[0]*5 for _ in range(5)]
for i in range(5):
    A[i][(i+1)%5] = 1
    A[i][(i+2)%5] = 1
coeffs = compute_uT_polynomial(5, A)
roots, q = find_roots_of_Q(coeffs, 5)
print(f"C_5 (Paley T_5): Q={q}, roots={roots}")

# Paley T_7
A = [[0]*7 for _ in range(7)]
QR = {1, 2, 4}  # QR mod 7
for i in range(7):
    for j in range(7):
        if i != j and (j - i) % 7 in QR:
            A[i][j] = 1
coeffs = compute_uT_polynomial(7, A)
roots, q = find_roots_of_Q(coeffs, 7)
print(f"Paley T_7: Q={q}, roots={roots}")
uT_at_1 = sum(coeffs.get(k, 0) for k in range(8))
print(f"  u_T(1) = {uT_at_1} (should be H=189)")

# Transitive T_7
A = [[0]*7 for _ in range(7)]
for i in range(7):
    for j in range(i+1, 7):
        A[i][j] = 1
coeffs = compute_uT_polynomial(7, A)
roots, q = find_roots_of_Q(coeffs, 7)
print(f"Transitive T_7: Q={q}, roots={roots}")
uT_at_1 = sum(coeffs.get(k, 0) for k in range(8))
print(f"  u_T(1) = {uT_at_1} (should be H=1)")
