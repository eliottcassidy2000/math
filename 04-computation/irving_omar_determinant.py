#!/usr/bin/env python3
"""
Irving-Omar determinantal approach to OCF and real-rootedness.

The Irving-Omar formula: W_D(z) = det(I + zXA_bar) / det(I - zXA)
where A is the adjacency matrix, A_bar is complement adjacency,
X = diag(x_1,...,x_n).

For tournaments: A + A_bar = J - I (all-ones minus identity).
So A_bar = J - I - A.

Key idea: H(T) = coefficient of x_1*x_2*...*x_n in W_D(z) evaluated at z=1.
The multilinear extraction gives the Hamiltonian path count.

Can we extract the independence polynomial structure from this determinant?

Also investigate: for the SKEW-ADJACENCY matrix S = A - A^T of a tournament,
det(xI - S) is a polynomial with known properties (all eigenvalues purely imaginary).
Is there a connection to I(Omega(T), x)?

Author: opus-2026-03-06-S18
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, hamiltonian_path_count, find_odd_cycles, conflict_graph
import numpy as np
from collections import defaultdict
from itertools import combinations

def indep_poly_coefficients(adj):
    m = len(adj)
    if m == 0:
        return [1]
    nbr = [0] * m
    for i in range(m):
        for j in range(m):
            if adj[i][j]:
                nbr[i] |= 1 << j
    counts = defaultdict(int)
    for mask in range(1 << m):
        ok = True
        seen = 0
        temp = mask
        while temp:
            v = (temp & -temp).bit_length() - 1
            if nbr[v] & seen:
                ok = False
                break
            seen |= 1 << v
            temp &= temp - 1
        if ok:
            counts[bin(mask).count('1')] += 1
    mx = max(counts.keys()) if counts else 0
    return [counts.get(k, 0) for k in range(mx + 1)]

print("=" * 70)
print("IRVING-OMAR DETERMINANTAL ANALYSIS")
print("=" * 70)

# Part 1: Skew-adjacency matrix spectrum vs I(Omega(T), x)
print("\n--- Part 1: Skew-adjacency eigenvalues vs I(Omega(T), x) ---")
print("  Tournament skew-adjacency S = A - A^T has purely imaginary eigenvalues.")
print("  det(xI - S) = char poly. Eigenvalues come in conjugate pairs +/- bi.")
print("  Question: is there a connection between {b_i} and I(Omega(T), x) roots?")

n = 5
m = n * (n - 1) // 2
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    A = np.array(T, dtype=float)
    S = A - A.T  # Skew-adjacency

    eigs_S = np.linalg.eigvals(S)
    # Sort by imaginary part
    eigs_S = sorted(eigs_S, key=lambda x: abs(x.imag), reverse=True)

    cycles = find_odd_cycles(T)
    if not cycles:
        continue

    cg = conflict_graph(cycles)
    coeffs = indep_poly_coefficients(cg)

    if len(coeffs) - 1 >= 2:
        # Compute roots of I(Omega, x)
        p = list(reversed(coeffs))
        roots_I = sorted(np.roots(p), key=lambda x: x.real)

        if bits < 20:
            print(f"\n  bits={bits}, H={hamiltonian_path_count(T)}")
            print(f"    S eigenvalues: {[f'{e.imag:.3f}i' for e in eigs_S if abs(e.imag) > 0.001]}")
            print(f"    I(Omega,x) coeffs: {coeffs}")
            print(f"    I(Omega,x) roots: {[f'{r.real:.4f}' for r in roots_I]}")

# Part 2: Hafnian / Pfaffian connection
print(f"\n\n--- Part 2: Pfaffian of skew-adjacency ---")
print("  For even n, pf(S)^2 = det(S). For tournament: det(S) = product of eigenvalues.")
print("  At odd n, det(S) = 0 (one zero eigenvalue).")

for n in [4, 5, 6]:
    m_bits = n * (n - 1) // 2
    count = min(1 << m_bits, 200)
    det_values = []
    for bits in range(count):
        T = tournament_from_bits(n, bits)
        A = np.array(T, dtype=float)
        S = A - A.T
        det_values.append(round(np.linalg.det(S)))

    unique_dets = sorted(set(det_values))
    print(f"  n={n}: det(S) values = {unique_dets[:20]}")
    if n % 2 == 0:
        # Pfaffian squared = det for skew-symmetric
        pf_values = sorted(set(int(round(np.sqrt(abs(d)))) for d in det_values if d > 0))
        print(f"  n={n}: |pf(S)| values = {pf_values[:20]}")

# Part 3: Connection between det(I + xS) and I(Omega, x)?
print(f"\n\n--- Part 3: det(I + xS) vs I(Omega(T), x) ---")
print("  For skew-symmetric S, det(I + xS) = sum_{k even} pf_k * x^k")
print("  where pf_k involves Pfaffians of principal submatrices.")
print("  Question: any relation to independence polynomial?")

n = 5
m = n * (n - 1) // 2

for bits in [2, 4, 6, 10, 14, 20, 30]:
    if bits >= (1 << m):
        break
    T = tournament_from_bits(n, bits)
    A = np.array(T, dtype=float)
    S = A - A.T
    H = hamiltonian_path_count(T)

    cycles = find_odd_cycles(T)
    if not cycles:
        continue
    cg = conflict_graph(cycles)
    ip_coeffs = indep_poly_coefficients(cg)

    # Compute det(I + xS) as polynomial in x
    # det(I + xS) = sum of all principal minors of S weighted by x^k
    # For skew-symmetric: only even powers survive
    det_coeffs = []
    for k in range(n + 1):
        if k == 0:
            det_coeffs.append(1)
        else:
            # sum of k x k principal minors of S
            total = 0
            for subset in combinations(range(n), k):
                sub = S[np.ix_(list(subset), list(subset))]
                total += np.linalg.det(sub)
            det_coeffs.append(round(total))

    print(f"  bits={bits}, H={H}")
    print(f"    I(Omega,x) = {ip_coeffs}")
    print(f"    det(I+xS) coeffs = {det_coeffs}")

# Part 4: Direct test of H(T) from determinant
print(f"\n\n--- Part 4: H(T) from matrix permanent ---")
print("  H(T) = permanent of A restricted to Hamiltonian paths.")
print("  Actually, H(T) = sum over permutations of product A[pi(i), pi(i+1)].")
print("  This is NOT the permanent. It's a path-permanent.")
print("  Can we express it as a determinant or Pfaffian?")

# At n=4, test whether H(T) relates to any simple matrix invariant
n = 4
m = n * (n - 1) // 2
print(f"\n  n={n}: checking matrix invariants vs H(T)")
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    A = np.array(T, dtype=float)
    S = A - A.T
    H = hamiltonian_path_count(T)

    det_S = round(np.linalg.det(S))
    trace_S2 = round(np.trace(S @ S))  # = -2 * #edges with orientation = -n(n-1)
    trace_S4 = round(np.trace(np.linalg.matrix_power(S, 4)))
    trace_A2 = round(np.trace(A @ A))
    trace_A3 = round(np.trace(np.linalg.matrix_power(A, 3)))  # = 6 * c_3

    if bits < 10:
        c3 = trace_A3 // 6
        scores = sorted([int(sum(T[i])) for i in range(n)])
        print(f"    bits={bits}: H={H}, det(S)={det_S}, tr(A^3)/6={c3}, scores={scores}")

# Part 5: Key insight - trace of A^{2k+1} counts directed (2k+1)-cycles
print(f"\n\n--- Part 5: Directed cycle counts from traces ---")
print("  tr(A^k) counts closed directed walks of length k.")
print("  For tournaments: tr(A^3) = 6*c_3, tr(A^5) = 10*c_5 + ...(shorter cycles)")

n = 7
from tournament_lib import random_tournament
for trial in range(5):
    T = random_tournament(n)
    A = np.array(T, dtype=float)

    traces = {}
    for k in [3, 5, 7]:
        traces[k] = round(np.trace(np.linalg.matrix_power(A, k)))

    cycles = find_odd_cycles(T)
    cycle_by_len = defaultdict(int)
    for c in (cycles or []):
        cycle_by_len[len(c)] += 1

    H = hamiltonian_path_count(T)
    cg = conflict_graph(cycles) if cycles else []
    ip = indep_poly_coefficients(cg) if cycles else [1]

    print(f"\n  trial {trial}: H={H}, I(Omega,x)={ip}")
    print(f"    cycle counts: {dict(cycle_by_len)}")
    print(f"    tr(A^3)={traces[3]} (expect 6*c3={6*cycle_by_len.get(3,0) if cycles else 0})")

# Part 6: The characteristic polynomial of S
print(f"\n\n--- Part 6: char poly of S vs I(Omega, x) ---")
print("  S eigenvalues are purely imaginary: +/- b_i*j.")
print("  char poly = prod(x^2 + b_i^2) * x^(n mod 2)")
print("  = x^e * (x^2 + b_1^2)(x^2 + b_2^2)...")
print("  Can we relate {b_i^2} to I(Omega, x) coefficients?")

n = 5
m = n * (n - 1) // 2
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    A = np.array(T, dtype=float)
    S = A - A.T

    eigs = np.linalg.eigvals(S)
    b_squared = sorted([round(abs(e.imag)**2, 4) for e in eigs if abs(e.imag) > 0.01], reverse=True)

    cycles = find_odd_cycles(T)
    if not cycles:
        continue
    cg = conflict_graph(cycles)
    ip = indep_poly_coefficients(cg)
    H = hamiltonian_path_count(T)

    if len(ip) >= 3 and bits < 50:
        # For degree-2 I(Omega,x) = 1 + a1*x + a2*x^2
        # Roots: x = (-a1 +/- sqrt(a1^2 - 4a2)) / (2*a2)
        a1, a2 = ip[1], ip[2]
        disc = a1**2 - 4*a2
        print(f"  bits={bits}: H={H}, ip={ip}, b^2={b_squared}, disc={disc}")

print(f"\n{'='*70}")
print("DONE")
print("=" * 70)
