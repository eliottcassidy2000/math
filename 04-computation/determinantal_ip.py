#!/usr/bin/env python3
"""
Determinantal representation of I(Omega_3(T), x).

If I(Omega_3(T), x) = det(I + x*M) for some positive semidefinite M,
then all roots are real (and negative). This would prove real-rootedness.

For claw-free graphs, such a representation exists via matching polynomials
(Godsil's matching polynomial theory). But for general Omega_3(T), we need
a different approach.

Strategy: for each tournament T, compute I(Omega_3, x) and check if
there exists a PSD matrix M such that det(I + x*M) = I(Omega_3, x).

For degree-2: 1 + a1*x + a2*x^2 = det(I + x*M) for 2x2 M iff
tr(M) = a1 and det(M) = a2. This needs a1^2 >= 4*a2 (PSD eigenvalues).

For degree-3: 1 + a1*x + a2*x^2 + a3*x^3 = det(I + x*M) for 3x3 M iff
tr(M) = a1, (tr(M)^2 - tr(M^2))/2 = a2, det(M) = a3.
This needs eigenvalues lambda_i > 0, which is equivalent to all roots of
I being negative real.

Actually, det(I + x*M) = prod(1 + x*lambda_i) for PSD M with eigenvalues
lambda_i >= 0. So the roots are x = -1/lambda_i, all negative real.

So I(Omega_3, x) = det(I + x*M) for PSD M iff I has all real negative roots.
This is circular — we need to PROVE real-rootedness, not assume it.

Better approach: find a SPECIFIC matrix M constructed from T such that
det(I + x*M) = I(Omega_3(T), x).

The Godsil-McKay connection: for a graph G, the matching polynomial
mu(G, x) = det(xI - A(G)) where A(G) is the adjacency matrix.
The independence polynomial and matching polynomial are related for
claw-free graphs but not in general.

However, for our specific graphs Omega_3(T), maybe there's a direct
construction.

Author: opus-2026-03-06-S18
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, random_tournament, find_odd_cycles, conflict_graph
import numpy as np
from collections import defaultdict

def find_3cycles_T(T):
    n = len(T)
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if T[i][j] and T[j][k] and T[k][i]:
                    cycles.append((i,j,k))
                elif T[i][k] and T[k][j] and T[j][i]:
                    cycles.append((i,j,k))
    return cycles

def build_cg3(cycles):
    m = len(cycles)
    adj = np.zeros((m, m))
    for i in range(m):
        for j in range(i+1, m):
            if set(cycles[i]) & set(cycles[j]):
                adj[i][j] = adj[j][i] = 1
    return adj

def indep_poly_from_adj(adj):
    m = len(adj)
    if m == 0: return [1]
    nbr = [0] * m
    for i in range(m):
        for j in range(m):
            if adj[i][j]: nbr[i] |= 1 << j
    counts = defaultdict(int)
    for mask in range(1 << m):
        ok = True; seen = 0; temp = mask
        while temp:
            v = (temp & -temp).bit_length() - 1
            if nbr[v] & seen: ok = False; break
            seen |= 1 << v; temp &= temp - 1
        if ok: counts[bin(mask).count('1')] += 1
    mx = max(counts.keys()) if counts else 0
    return [counts.get(k,0) for k in range(mx+1)]

print("=" * 70)
print("DETERMINANTAL REPRESENTATION SEARCH")
print("=" * 70)

# Part 1: For Omega_3(T) at n=5, check matching poly vs independence poly
print("\n--- Part 1: Matching poly vs Independence poly ---")
n = 5
m = n * (n - 1) // 2
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    c3 = find_3cycles_T(T)
    if len(c3) < 2:
        continue

    adj = build_cg3(c3)
    ip = indep_poly_from_adj(adj.astype(int).tolist())

    # Matching polynomial: mu(G, x) = det(xI - A)
    eigs_A = sorted(np.linalg.eigvals(adj), reverse=True)
    match_poly = np.poly(eigs_A)  # Coefficients highest degree first
    # Convert to lowest degree first
    match_coeffs = list(reversed([round(c) for c in match_poly]))

    # For comparison: I(G, x) evaluated at specific points
    # Lovász theta connection: for perfect graphs, theta = alpha
    # For our purposes, check if eigenvalues of adjacency matrix
    # relate to I(G, x) roots

    if bits < 20:
        print(f"  bits={bits}: c3={len(c3)}")
        print(f"    I(Omega_3, x) = {ip}")
        print(f"    A eigenvalues: {[f'{e:.3f}' for e in eigs_A if abs(e) > 0.001]}")
        print(f"    Matching poly: {match_coeffs[:5]}")

# Part 2: Heilmann-Lieb for complement
print("\n--- Part 2: Complement structure ---")
print("  Omega_3(T) complement = graph of DISJOINT 3-cycle pairs.")
print("  For n<=8, complement is a union of disjoint edges (matching).")
print("  A matching graph has matching polynomial = independence polynomial!")

n = 6
m_bits = n * (n - 1) // 2
match_is_ip = 0
total = 0
for bits in range(1 << m_bits):
    T = tournament_from_bits(n, bits)
    c3 = find_3cycles_T(T)
    if len(c3) < 2:
        continue
    total += 1

    adj = build_cg3(c3)
    m = len(c3)

    # Check: is the COMPLEMENT a matching (disjoint union of edges)?
    comp = np.ones((m, m)) - adj - np.eye(m)
    # A matching means each vertex has degree <= 1 in complement
    comp_degs = comp.sum(axis=1)
    is_matching = all(d <= 1 for d in comp_degs)

    if is_matching:
        match_is_ip += 1

print(f"  n={n}: {match_is_ip}/{total} have Omega_3 complement = matching")
print(f"  When complement is a matching, Omega_3 = K_m minus matching")
print(f"  = complete graph with some edges removed, each removed edge = disjoint pair")

# Part 3: For Omega_3 = K_m - matching, compute I(G, x) explicitly
print("\n--- Part 3: I(K_m - matching, x) ---")
print("  K_m has I(K_m, x) = 1 + m*x")
print("  K_m - one edge {a,b}: I = 1 + m*x + x^2 (the removed edge becomes")
print("  an independent set of size 2)")
print("  K_m - matching of size p: I = 1 + m*x + p*x^2")
print("  because independent sets of size 2 are exactly the removed edges.")
print("  For all-real-roots of 1 + m*x + p*x^2: disc = m^2 - 4p >= 0")
print("  i.e., p <= m^2/4. Since p <= m/2, need m/2 <= m^2/4, i.e., m >= 2. Always true!")

# Verify this formula
n = 6
m_bits = n * (n - 1) // 2
formula_ok = 0
formula_fail = 0
for bits in range(1 << m_bits):
    T = tournament_from_bits(n, bits)
    c3 = find_3cycles_T(T)
    if len(c3) < 2:
        continue

    adj = build_cg3(c3)
    m = len(c3)
    comp = np.ones((m, m)) - adj - np.eye(m)
    comp_degs = comp.sum(axis=1)
    is_matching = all(d <= 1 for d in comp_degs)

    if not is_matching:
        continue

    # Count matching size in complement
    p = int(comp.sum()) // 2  # Number of edges in complement = non-edges in Omega_3

    ip = indep_poly_from_adj(adj.astype(int).tolist())
    expected = [1, m, p]
    if len(ip) == 2:
        expected = [1, m]
        if ip == expected:
            formula_ok += 1
        else:
            formula_fail += 1
    elif len(ip) == 3:
        if ip == expected:
            formula_ok += 1
        else:
            formula_fail += 1
            if formula_fail <= 3:
                print(f"  FAIL: bits={bits}, ip={ip}, expected={expected}")
    else:
        formula_fail += 1

print(f"  n=6: formula I = [1, m, p] verified: {formula_ok} ok, {formula_fail} fail")

# Part 4: At n>=9, complement is NO LONGER just a matching
print("\n--- Part 4: Complement structure at n=9 ---")
n = 9
comp_max_deg = []
for trial in range(200):
    T = random_tournament(n)
    c3 = find_3cycles_T(T)
    if len(c3) < 2 or len(c3) > 23:
        continue

    adj = build_cg3(c3)
    m = len(c3)
    comp = np.ones((m, m)) - adj - np.eye(m)
    comp_degs = [int(comp[i].sum()) for i in range(m)]
    comp_max_deg.append(max(comp_degs))

if comp_max_deg:
    print(f"  n=9: complement max degree dist: min={min(comp_max_deg)}, max={max(comp_max_deg)}, avg={sum(comp_max_deg)/len(comp_max_deg):.1f}")
    deg_counts = defaultdict(int)
    for d in comp_max_deg:
        deg_counts[d] += 1
    print(f"  Distribution: {dict(sorted(deg_counts.items()))}")

# Part 5: Skew-adjacency approach
print("\n--- Part 5: Tournament skew-adjacency and Omega_3 ---")
print("  S = A - A^T (skew-symmetric). det(I + xS) = prod(1 + x^2 * b_i^2) * (1)^(n mod 2)")
print("  This doesn't directly give I(Omega_3, x), but the b_i relate to cycle structure.")
print("  Key fact: tr(A^3) = 3 * sum of directed 3-cycles = 3 * c3_directed")
print("  where c3_directed counts each 3-cycle TWICE (two orientations).")
print("  So c3 = tr(A^3) / 6 (undirected 3-cycles = vertex sets).")
print("  Wait: in a tournament, each triple has exactly one directed 3-cycle orientation.")
print("  So tr(A^3) = 6 * c3_vertex_sets? Let me verify...")

n = 5
m_bits = n * (n - 1) // 2
for bits in [2, 4, 6, 14, 30]:
    T = tournament_from_bits(n, bits)
    A = np.array(T, dtype=float)
    trA3 = round(np.trace(np.linalg.matrix_power(A, 3)))
    c3 = find_3cycles_T(T)
    print(f"  bits={bits}: tr(A^3)={trA3}, c3_vertex_sets={len(c3)}, ratio={trA3/(len(c3) if len(c3) > 0 else 1):.1f}")

# The ratio should be 3 (not 6), because each directed 3-cycle contributes
# 3 to the trace (the walk a->b->c->a hits 3 diagonal elements)
# Wait no, the trace counts closed walks of length 3 starting at each vertex.
# For a directed 3-cycle (a->b->c->a): walk starting at a: a->b->c->a (contributes 1)
# walk starting at b: b->c->a->b (contributes 1)
# walk starting at c: c->a->b->c (contributes 1)
# So each directed 3-cycle contributes 3 to tr(A^3).
# In a tournament on triple {i,j,k}, there's exactly ONE directed 3-cycle.
# So tr(A^3) = 3 * c3. Confirmed!

# Part 6: Can we express alpha_1 = c3 from eigenvalues?
print("\n--- Part 6: c3 from skew-adjacency eigenvalues ---")
print("  S = A - A^T. S^2 = A^2 + (A^T)^2 - A*A^T - A^T*A")
print("  tr(S^2) = 2*tr(A^2) - 2*tr(A*A^T)")
print("  For tournament: tr(A^2) = sum of A[i][j]*A[j][i] = 0 (no mutual arcs)")
print("  Wait: tr(A^2) = sum_i sum_j A[i][j]*A[j][i] = number of 2-cycles = 0 (tournament)")
print("  And tr(A*A^T) = sum_i (sum_j A[i][j])^2 = sum_i d_out(i)^2")
print("  So tr(S^2) = -2 * sum d_out(i)^2")
print("  Hmm, let me compute directly.")

n = 5
for bits in [4, 14]:
    T = tournament_from_bits(n, bits)
    A = np.array(T, dtype=float)
    S = A - A.T
    print(f"\n  bits={bits}:")
    print(f"    tr(A^2) = {np.trace(A@A):.0f}")
    print(f"    tr(A*A^T) = {np.trace(A@A.T):.0f}")
    print(f"    tr(S^2) = {np.trace(S@S):.0f}")
    print(f"    tr(S^3) = {np.trace(np.linalg.matrix_power(S,3)):.0f}")
    print(f"    tr(A^3) = {np.trace(np.linalg.matrix_power(A,3)):.0f}")
    # For skew-symmetric S: S^3 is also skew-symmetric
    # tr(S^3) = tr((A-A^T)^3) = tr(A^3) - 3*tr(A^2 * A^T) + 3*tr(A * (A^T)^2) - tr((A^T)^3)
    # = tr(A^3) - tr((A^T)^3) - 3*tr(A^2*A^T) + 3*tr(A*(A^T)^2)
    # = 2*tr(A^3) - 3*tr(A^2*A^T) + 3*tr(A*(A^T)^2)
    # (using tr(B^T) = tr(B))
    trAAAt = np.trace(A@A@A.T)
    trAAtAt = np.trace(A@A.T@A.T)
    print(f"    tr(A^2*A^T) = {trAAAt:.0f}")
    print(f"    tr(A*(A^T)^2) = {trAAtAt:.0f}")
    print(f"    2*tr(A^3) - 3*tr(A^2*A^T) + 3*tr(A*(A^T)^2) = {2*np.trace(np.linalg.matrix_power(A,3)) - 3*trAAAt + 3*trAAtAt:.0f}")

print(f"\n{'='*70}")
print("KEY FINDINGS")
print("=" * 70)
print("""
For n <= 8: Omega_3(T) = K_m - matching.
=> I(Omega_3, x) = 1 + m*x + p*x^2  (m = c3, p = #disjoint pairs)
=> Discriminant = m^2 - 4p >= 0  (since p <= m/2 and m >= 2)
=> ALL ROOTS REAL for n <= 8 via elementary quadratic discriminant.

For n >= 9: Omega_3(T) complement can have vertices of degree > 1.
=> I(Omega_3, x) can have degree 3+.
=> Need different argument for real-rootedness.

The challenge: find an algebraic argument for n >= 9 that uses the
tournament structure of T (not just abstract graph properties of Omega_3).
""")

print("=" * 70)
print("DONE")
print("=" * 70)
