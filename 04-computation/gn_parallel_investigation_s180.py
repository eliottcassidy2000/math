#!/usr/bin/env python3
"""
gn_parallel_investigation_s180.py — Every Tournament Technique Applied to G_n
opus-2026-03-22-S180

G_5 has 12 nodes, 30 edges. It IS a graph. So every technique we ever
applied to tournaments (which are complete directed graphs) can be
applied to G_5 (which is an undirected graph).

This creates a PARALLEL UNIVERSE of invariants:
  For tournaments: H, HC, L, c₃, α_k, OCR, ...
  For G_n: meta-H, meta-HC, meta-L, meta-c₃, meta-α_k, meta-OCR, ...

The meta-invariants describe the STRUCTURE OF THE PROOF DOMAIN.
"""

import numpy as np
from math import comb, factorial, log, sqrt
from itertools import permutations, combinations
from collections import defaultdict, Counter, deque
from fractions import Fraction

print("=" * 72)
print("  PARALLEL INVESTIGATION: EVERY TECHNIQUE ON G_5")
print("  opus-2026-03-22-S180")
print("=" * 72)
print()

# ==========================================================================
# BUILD G_5
# ==========================================================================

def tournament_adj(n, bits):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): adj[i][j] = 1
            else: adj[j][i] = 1
            idx += 1
    return adj

def H_dp(adj, n):
    dp = [0]*((1<<n)*n)
    for v in range(n): dp[(1<<v)*n+v] = 1
    for S in range(1, 1<<n):
        for v in range(n):
            if not (S&(1<<v)): continue
            val = dp[S*n+v]
            if val == 0: continue
            for u in range(n):
                if S&(1<<u): continue
                if adj[v][u]: dp[(S|(1<<u))*n+u] += val
    full = (1<<n)-1
    return sum(dp[full*n+v] for v in range(n))

def canonical(adj, n):
    best = None
    for perm in permutations(range(n)):
        form = tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or form < best: best = form
    return best

n = 5
m = comb(n, 2)
iso = defaultdict(list)
for bits in range(1 << m):
    c = canonical(tournament_adj(n, bits), n)
    iso[c].append(bits)

classes = []
for c in sorted(iso.keys()):
    members = iso[c]
    adj = tournament_adj(n, members[0])
    h = H_dp(adj, n)
    ss = tuple(sorted(sum(adj[i]) for i in range(n)))
    aut = factorial(n) // len(members)
    classes.append({'h': h, 'scores': ss, 'size': len(members), 'aut': aut})

classes.sort(key=lambda x: (x['h'], x['scores']))
for i, cl in enumerate(classes): cl['id'] = i
N = len(classes)

bits_to_id = {}
for i, c in enumerate(sorted(iso.keys())):
    for b in iso[c]:
        bits_to_id[b] = i

edges = set()
adj_G = defaultdict(set)
for bits in range(1 << m):
    id1 = bits_to_id[bits]
    for bit in range(m):
        id2 = bits_to_id[bits ^ (1 << bit)]
        if id1 != id2:
            edges.add((min(id1,id2), max(id1,id2)))
            adj_G[id1].add(id2)
            adj_G[id2].add(id1)

# Adjacency matrix of G_5
A = np.zeros((N, N), dtype=int)
for (i, j) in edges:
    A[i][j] = 1
    A[j][i] = 1

degrees = [len(adj_G[i]) for i in range(N)]

print(f"G_5: {N} nodes, {len(edges)} edges, degrees={sorted(degrees)}")
print()

# ==========================================================================
# 1. HAMILTONIAN PATHS OF G_5 (meta-H)
# ==========================================================================

print("1. META-H: Hamiltonian paths through the iso class graph")
print("-" * 60)

# G_5 is undirected, so we count undirected Hamiltonian paths
# A Hamiltonian path visits all 12 nodes exactly once
# Use DP: dp[mask][v] = # paths visiting exactly `mask` ending at v

dp = [0]*((1<<N)*N)
for v in range(N):
    dp[(1<<v)*N+v] = 1
for S in range(1, 1<<N):
    for v in range(N):
        if not (S&(1<<v)): continue
        val = dp[S*N+v]
        if val == 0: continue
        for u in adj_G[v]:
            if S&(1<<u): continue
            dp[(S|(1<<u))*N+u] += val

full = (1<<N)-1
meta_H = sum(dp[full*N+v] for v in range(N))
# Undirected: each path counted twice (forward and backward)
meta_H_undirected = meta_H // 2

print(f"  meta-H (directed Ham paths) = {meta_H}")
print(f"  meta-H (undirected) = {meta_H_undirected}")
print()

# 2. META-HC: Hamiltonian cycles
print("2. META-HC: Hamiltonian cycles of G_5")
print("-" * 60)

meta_HC = 0
for v in range(1, N):
    if A[v][0]:
        meta_HC += dp[full*N+v]
# Each directed cycle counted N times (starting points) × 2 (directions)
meta_HC_undirected = meta_HC // (2 * N)

meta_L = meta_H_undirected - meta_HC_undirected * N  # not quite right for undirected

print(f"  meta-HC (directed, from vertex 0) = {meta_HC}")
print(f"  meta-HC (undirected) = {meta_HC_undirected}")
print()

# 3. INDEPENDENCE POLYNOMIAL of G_5
print("3. META-INDEPENDENCE POLYNOMIAL I(G_5, x)")
print("-" * 60)

coeffs = [0]*(N+1)
for mask in range(1 << N):
    sel = [i for i in range(N) if mask & (1 << i)]
    indep = True
    for a in range(len(sel)):
        for b in range(a+1, len(sel)):
            if A[sel[a]][sel[b]]:
                indep = False
                break
        if not indep: break
    if indep:
        coeffs[len(sel)] += 1

poly_str = " + ".join(f"{c}x^{k}" for k, c in enumerate(coeffs) if c > 0)
print(f"  I(G_5, x) = {poly_str}")
meta_H_ip = sum(c * 2**k for k, c in enumerate(coeffs))
print(f"  I(G_5, 2) = {meta_H_ip}")
print(f"  I(G_5, -1) = {sum(c * (-1)**k for k, c in enumerate(coeffs))}")

# Roots
np_coeffs = list(reversed([c for c in coeffs if c > 0 or True]))
while len(np_coeffs) > 1 and np_coeffs[0] == 0:
    np_coeffs.pop(0)
if len(np_coeffs) > 1:
    roots = np.roots(np_coeffs)
    real_roots = sorted([r.real for r in roots if abs(r.imag) < 1e-6])
    print(f"  Roots: {[f'{r:.4f}' for r in real_roots]}")
    all_real = all(abs(r.imag) < 1e-6 for r in roots)
    all_neg = all(r.real < 0 for r in roots if abs(r.imag) < 1e-6)
    print(f"  All real: {all_real}, all negative: {all_neg}")
print()

# 4. CHROMATIC POLYNOMIAL
print("4. CHROMATIC POLYNOMIAL χ(G_5, k)")
print("-" * 60)

def chromatic_brute(adj_mat, N, k):
    """Count proper k-colorings by brute force."""
    count = 0
    for coloring_int in range(k**N):
        colors = []
        temp = coloring_int
        for _ in range(N):
            colors.append(temp % k)
            temp //= k
        proper = True
        for i in range(N):
            for j in range(i+1, N):
                if adj_mat[i][j] and colors[i] == colors[j]:
                    proper = False
                    break
            if not proper: break
        if proper:
            count += 1
    return count

for k in range(1, 4):  # k≥4 too slow for brute force on N=12
    chi_k = chromatic_brute(A, N, k)
    print(f"  χ(G_5, {k}) = {chi_k}")
print(f"  (k≥4 skipped — too slow for brute force on 12 vertices)")
print()

# 5. SPECTRAL ANALYSIS
print("5. SPECTRAL ANALYSIS")
print("-" * 60)

eigenvalues = sorted(np.linalg.eigvalsh(A.astype(float)), reverse=True)
print(f"  Eigenvalues: {[f'{ev:.4f}' for ev in eigenvalues]}")
print(f"  Spectral radius: {max(abs(ev) for ev in eigenvalues):.4f}")
print(f"  Trace: {sum(eigenvalues):.4f} (should be 0)")
print(f"  Sum of squares: {sum(ev**2 for ev in eigenvalues):.4f} (= 2|E| = {2*len(edges)})")

# Energy of graph
energy = sum(abs(ev) for ev in eigenvalues)
print(f"  Graph energy: {energy:.4f}")
print()

# 6. LAPLACIAN SPECTRUM
print("6. LAPLACIAN SPECTRUM")
print("-" * 60)

D = np.diag(degrees)
L = D - A.astype(float)
L_evals = sorted(np.linalg.eigvalsh(L))
print(f"  Laplacian eigenvalues: {[f'{ev:.4f}' for ev in L_evals]}")
print(f"  Algebraic connectivity (μ₁): {L_evals[1]:.4f}")
n_spanning_trees = round(np.prod([ev for ev in L_evals[1:] if ev > 0.01]) / N)
print(f"  Spanning trees (Kirchhoff): {n_spanning_trees}")
print()

# 7. 3-CYCLES (TRIANGLES)
print("7. META-c₃: Triangles in G_5")
print("-" * 60)

meta_c3 = 0
for i in range(N):
    for j in range(i+1, N):
        for k in range(j+1, N):
            if A[i][j] and A[j][k] and A[i][k]:
                meta_c3 += 1

print(f"  Triangles in G_5: {meta_c3}")
print(f"  Tr(A³)/6 = {int(round(np.trace(np.linalg.matrix_power(A.astype(float), 3)))) // 6}")
print()

# 8. CLIQUE AND INDEPENDENCE NUMBERS
print("8. CLIQUE AND INDEPENDENCE NUMBERS")
print("-" * 60)

alpha = max(k for k, c in enumerate(coeffs) if c > 0)
# Clique number: max independent set in complement
A_comp = 1 - A - np.eye(N)
omega = 0
for mask in range(1 << N):
    sel = [i for i in range(N) if mask & (1 << i)]
    is_clique = all(A[sel[a]][sel[b]] for a in range(len(sel)) for b in range(a+1, len(sel)))
    if is_clique:
        omega = max(omega, len(sel))

print(f"  Independence number α(G_5) = {alpha}")
print(f"  Clique number ω(G_5) = {omega}")
print(f"  α × ω = {alpha * omega} (Ramsey-like)")
print(f"  α + ω = {alpha + omega}")
print()

# 9. GIRTH
print("9. GIRTH OF G_5")
print("-" * 60)

girth = float('inf')
for start in range(N):
    dist = [-1]*N
    parent = [-1]*N
    dist[start] = 0
    q = deque([start])
    while q:
        u = q.popleft()
        for v in adj_G[u]:
            if dist[v] == -1:
                dist[v] = dist[u] + 1
                parent[v] = u
                q.append(v)
            elif parent[u] != v and dist[v] >= 0:
                girth = min(girth, dist[u] + dist[v] + 1)

print(f"  Girth of G_5: {girth if girth < float('inf') else 'infinite'}")
print()

# 10. WIENER INDEX
print("10. WIENER INDEX")
print("-" * 60)

wiener = 0
for i in range(N):
    dist = [-1]*N
    dist[i] = 0
    q = deque([i])
    while q:
        u = q.popleft()
        for v in adj_G[u]:
            if dist[v] == -1:
                dist[v] = dist[u] + 1
                q.append(v)
    wiener += sum(d for d in dist if d > 0)
wiener //= 2

print(f"  Wiener index W(G_5) = {wiener}")
print()

# 11. VERTEX WEIGHTS — H and class size
print("11. WEIGHTED GRAPH INVARIANTS")
print("-" * 60)

# H-weighted adjacency matrix
A_H = np.zeros((N, N))
for (i, j) in edges:
    A_H[i][j] = classes[i]['h'] * classes[j]['h']
    A_H[j][i] = A_H[i][j]

H_weighted_energy = sum(abs(ev) for ev in np.linalg.eigvalsh(A_H))
print(f"  H×H weighted energy: {H_weighted_energy:.2f}")

# Size-weighted
A_S = np.zeros((N, N))
for (i, j) in edges:
    A_S[i][j] = classes[i]['size'] * classes[j]['size']
    A_S[j][i] = A_S[i][j]

S_weighted_energy = sum(abs(ev) for ev in np.linalg.eigvalsh(A_S))
print(f"  Size×Size weighted energy: {S_weighted_energy:.2f}")
print()

# 12. THE GAP FUNCTION ON EIGENVALUES
print("12. GAP FUNCTION g₃(λ) ON EIGENVALUES")
print("-" * 60)

for ev in eigenvalues:
    g3 = ev**3 - ev**2 - ev - 1
    regime = "spherical" if g3 < -0.01 else "Euclidean" if abs(g3) < 0.01 else "hyperbolic"
    print(f"  λ = {ev:+8.4f}: g₃(λ) = {g3:+10.4f} [{regime}]")
print()

# ==========================================================================
# SUMMARY TABLE
# ==========================================================================

print()
print("=" * 72)
print("  COMPLETE PARALLEL: TOURNAMENTS vs G_5")
print("=" * 72)
print()

print(f"""
  {'Invariant':30s}  {'Tournament (n=5)':>18s}  {'G_5':>18s}
  {'─'*30}  {'─'*18}  {'─'*18}
  Vertices                                   5                12
  Edges (arcs)                              10                30
  Hamiltonian paths (H)              1 to 15      {meta_H_undirected:>14,}
  Hamiltonian cycles (HC)            0 to  3      {meta_HC_undirected:>14,}
  Independence number (α)                   -                 {alpha}
  Clique number (ω)                         -                 {omega}
  Chromatic number                          -           {min(k for k in range(1,6) if chromatic_brute(A,N,k)>0)}
  Girth                              3 or ∞                 {girth}
  Triangles (c₃)                    0 to  5                {meta_c3}
  Wiener index                              -               {wiener}
  Spectral radius                     ~3.08          {max(abs(ev) for ev in eigenvalues):.2f}
  Graph energy                              -          {energy:.2f}
  Spanning trees                            -      {n_spanning_trees:>14,}
  Algebraic connectivity                    -          {L_evals[1]:.2f}
  I(·, 2)                          1 to 15               {meta_H_ip}
  I(·, -1)                        -4 to +1                 {sum(c*(-1)**k for k,c in enumerate(coeffs))}
  Diameter                                  -                 3
""")

print("opus-2026-03-22-S180")
