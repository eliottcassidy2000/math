#!/usr/bin/env python3
"""
parallel_eval_s20ah.py -- kind-pasteur-2026-03-22-S20ah

PARALLEL EVALUATION: Apply every tournament invariant to BOTH
a tournament T AND the iso class graph G_5.

G_5 is an undirected graph on 12 vertices with 30 edges.
Treat it as a DIGRAPH by replacing each undirected edge with two
directed arcs (the "symmetric digraph" or "bidirected graph").

For each invariant, compute it for:
- A specific tournament T_5 (the regular one, H=15)
- The iso class graph G_5 (as a digraph)

Compare: what does each invariant MEAN when applied to G_5?

Author: kind-pasteur-2026-03-22-S20ah
"""
import sys
import numpy as np
from math import comb, log2
from collections import defaultdict
from itertools import permutations
sys.stdout.reconfigure(line_buffering=True)

def count_hp(A, n):
    dp = defaultdict(int)
    for v in range(n): dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[(mask, v)] == 0: continue
            for w in range(n):
                if mask & (1 << w): continue
                if A[v][w]: dp[(mask | (1 << w), w)] += dp[(mask, v)]
    return sum(dp[((1 << n) - 1, v)] for v in range(n))

def count_hc(A, n):
    full = (1 << n) - 1
    dp = defaultdict(int)
    dp[(1, 0)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[(mask, v)] == 0: continue
            for w in range(n):
                if mask & (1 << w): continue
                if A[v][w]: dp[(mask | (1 << w), w)] += dp[(mask, v)]
    return sum(dp[(full, v)] for v in range(n) if A[v][0])

print("=" * 70)
print("  PARALLEL EVALUATION: TOURNAMENT T_5 vs ISO CLASS GRAPH G_5")
print("=" * 70)

# ================================================================
# BUILD G_5 (the iso class graph on 12 vertices)
# ================================================================
n_t = 5  # tournament size
pairs_t = [(i,j) for i in range(n_t) for j in range(i+1, n_t)]
m_t = len(pairs_t)

print(f"\n  Building G_5 from tournament iso classes at n={n_t}...")

# Build iso classes
H_map = {}
canon_map = {}
for bits in range(2**m_t):
    A = np.zeros((n_t,n_t), dtype=np.int8)
    for k, (i,j) in enumerate(pairs_t):
        if (bits >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1
    H_map[bits] = count_hp(A, n_t)
    cf = tuple(A[i][j] for i in range(n_t) for j in range(n_t))
    best = None
    for perm in permutations(range(n_t)):
        form = tuple(A[perm[i]][perm[j]] for i in range(n_t) for j in range(n_t))
        if best is None or form < best: best = form
    if best not in canon_map:
        canon_map[best] = []
    canon_map[best].append(bits)

classes = []
for cf, members in sorted(canon_map.items(), key=lambda x: H_map[x[1][0]]):
    classes.append({'id': len(classes), 'H': H_map[members[0]], 'size': len(members), 'members': members})

N = len(classes)

# Build G_5 adjacency (undirected -> bidirected digraph)
G5 = np.zeros((N, N), dtype=np.int8)
for c in classes:
    for bits in c['members']:
        for k in range(m_t):
            nb = bits ^ (1 << k)
            for c2 in classes:
                if nb in c2['members'] and c['id'] != c2['id']:
                    G5[c['id']][c2['id']] = 1
                    G5[c2['id']][c['id']] = 1
                    break

# Also build the H-oriented version (DAG)
G5_dag = np.zeros((N, N), dtype=np.int8)
for i in range(N):
    for j in range(N):
        if G5[i][j] and classes[i]['H'] < classes[j]['H']:
            G5_dag[i][j] = 1  # lower -> higher

# Choose T_5: the regular tournament (cycle 0->1->2->3->4->0 + 0->2, 1->3, 2->4, 3->0, 4->1)
T5 = np.array([
    [0,1,1,0,0],
    [0,0,1,1,0],
    [0,0,0,1,1],
    [1,0,0,0,1],
    [1,1,0,0,0]
], dtype=np.int8)
H_T5 = count_hp(T5, n_t)

print(f"  G_5: {N} vertices, {G5.sum()//2} undirected edges")
print(f"  T_5 (regular): H = {H_T5}")

# ================================================================
# PARALLEL EVALUATION TABLE
# ================================================================
print(f"\n{'='*70}")
print(f"  PARALLEL EVALUATION TABLE")
print(f"{'='*70}\n")

results = []

def add_result(name, val_T, val_G, interpretation=""):
    results.append((name, val_T, val_G, interpretation))

# --- 1. Basic counts ---
add_result("Vertices (n)", n_t, N, "T has 5 players; G has 12 iso classes")
add_result("Edges (arcs)", comb(n_t, 2), int(G5.sum()//2), "T has 10 arcs; G has 30 undirected edges")
add_result("Density", 1.0, G5.sum() / (N*(N-1)), "T is complete (1.0); G is sparse")

# --- 2. Score sequence ---
T5_scores = tuple(sorted(T5.sum(axis=1).astype(int)))
G5_degrees = tuple(sorted((G5.sum(axis=1)).astype(int)))
add_result("Score/Degree seq", str(T5_scores), str(G5_degrees), "")

# --- 3. H (Hamiltonian paths) ---
H_G5 = count_hp(G5, N) if N <= 12 else "too large"
add_result("H (Ham paths)", H_T5, H_G5, "")

# --- 4. HC (Hamiltonian cycles) ---
HC_T5 = count_hc(T5, n_t)
HC_G5 = count_hc(G5, N) if N <= 12 else "too large"
add_result("HC (Ham cycles)", HC_T5, HC_G5, "")

# --- 5. L = H - n*HC ---
L_T5 = H_T5 - n_t * HC_T5
L_G5 = H_G5 - N * HC_G5 if isinstance(H_G5, int) and isinstance(HC_G5, int) else "?"
add_result("L = H - n*HC", L_T5, L_G5, "")

# --- 6. c3 (3-cycles) ---
c3_T5 = 0
for i in range(n_t):
    for j in range(i+1, n_t):
        for k in range(j+1, n_t):
            if (T5[i][j] and T5[j][k] and T5[k][i]) or (T5[j][i] and T5[i][k] and T5[k][j]):
                c3_T5 += 1

c3_G5 = 0
for i in range(N):
    for j in range(i+1, N):
        for k in range(j+1, N):
            if G5[i][j] and G5[j][k] and G5[k][i]:
                c3_G5 += 1
            if G5[j][i] and G5[i][k] and G5[k][j]:
                c3_G5 += 1

add_result("c3 (3-cycles)", c3_T5, c3_G5, "Directed 3-cycles in digraph")

# --- 7. Strongly connected ---
def is_strongly_connected(A, n):
    visited = {0}; queue = [0]
    while queue:
        v = queue.pop(0)
        for w in range(n):
            if A[v][w] and w not in visited: visited.add(w); queue.append(w)
    if len(visited) < n: return False
    visited = {0}; queue = [0]
    while queue:
        v = queue.pop(0)
        for w in range(n):
            if A[w][v] and w not in visited: visited.add(w); queue.append(w)
    return len(visited) == n

sc_T5 = is_strongly_connected(T5, n_t)
sc_G5 = is_strongly_connected(G5, N)
add_result("Strongly connected", sc_T5, sc_G5, "")

# --- 8. Kings ---
def count_kings(A, n):
    A2 = A @ A
    reach = A + np.clip(A2, 0, 1)
    np.fill_diagonal(reach, 1)
    return sum(1 for v in range(n) if all(reach[v][w] > 0 or v == w for w in range(n)))

kings_T5 = count_kings(T5, n_t)
kings_G5 = count_kings(G5, N)
add_result("Kings", kings_T5, kings_G5, "Vertices reaching all others in <=2 steps")

# --- 9. Arborescences ---
def count_arb(A, n, root=0):
    A_f = A.astype(float)
    D_in = np.diag(A_f.sum(axis=0))
    L = D_in - A_f.T
    idx = [i for i in range(n) if i != root]
    return int(round(np.linalg.det(L[np.ix_(idx, idx)])))

arb_T5 = count_arb(T5, n_t)
arb_G5 = count_arb(G5, N)
add_result("Arborescences (root=0)", arb_T5, arb_G5, "Directed spanning trees")

# --- 10. Spectral radius ---
eig_T5 = np.sort(np.abs(np.linalg.eigvals(T5.astype(float))))[::-1]
eig_G5 = np.sort(np.abs(np.linalg.eigvals(G5.astype(float))))[::-1]
add_result("Spectral radius", f"{eig_T5[0]:.4f}", f"{eig_G5[0]:.4f}", "Largest |eigenvalue|")

# --- 11. Eigenvalues ---
eig_T5_full = sorted(np.linalg.eigvals(T5.astype(float)).real, reverse=True)
eig_G5_full = sorted(np.linalg.eigvals(G5.astype(float)).real, reverse=True)
add_result("Eigenvalues (real)", [f"{e:.2f}" for e in eig_T5_full], [f"{e:.2f}" for e in eig_G5_full[:6]] + ["..."], "")

# --- 12. Independence number ---
def indep_number(A, n):
    best = 0
    for mask in range(2**n):
        verts = [i for i in range(n) if (mask >> i) & 1]
        ok = True
        for a in range(len(verts)):
            for b in range(a+1, len(verts)):
                if A[verts[a]][verts[b]] or A[verts[b]][verts[a]]:
                    ok = False; break
            if not ok: break
        if ok: best = max(best, len(verts))
    return best

alpha_T5 = indep_number(T5, n_t)
alpha_G5 = indep_number(G5, N)
add_result("Independence number", alpha_T5, alpha_G5, "Max set with no arcs between")

# --- 13. Independence polynomial ---
def indep_poly(A, n):
    coeffs = [0] * (n + 1)
    for mask in range(2**n):
        verts = [i for i in range(n) if (mask >> i) & 1]
        ok = True
        for a in range(len(verts)):
            for b in range(a+1, len(verts)):
                if A[verts[a]][verts[b]] or A[verts[b]][verts[a]]:
                    ok = False; break
            if not ok: break
        if ok: coeffs[len(verts)] += 1
    return coeffs

ip_T5 = indep_poly(T5, n_t)
ip_G5 = indep_poly(G5, N)
add_result("I(G, 2)", sum(c * 2**k for k, c in enumerate(ip_T5)),
                        sum(c * 2**k for k, c in enumerate(ip_G5)),
           "Independence polynomial at x=2")

add_result("I(G, x) coeffs", ip_T5[:6], ip_G5[:7], "")

# --- 14. Chromatic number (of underlying undirected) ---
# For T5: underlying graph is K_5, chi = 5
# For G5: need to compute
# (skipping exact computation for speed)
add_result("Chromatic # (undirected)", 5, "<=5 (G_5 has max degree 7)", "")

# --- 15. Diameter ---
def diameter(A, n):
    dist = np.full((n, n), n+1)
    np.fill_diagonal(dist, 0)
    for i in range(n):
        for j in range(n):
            if A[i][j]: dist[i][j] = 1
    for k in range(n):
        for i in range(n):
            for j in range(n):
                if dist[i][k] + dist[k][j] < dist[i][j]:
                    dist[i][j] = dist[i][k] + dist[k][j]
    return int(dist.max()) if dist.max() < n+1 else float('inf')

diam_T5 = diameter(T5, n_t)
diam_G5 = diameter(G5, N)
add_result("Diameter", diam_T5, diam_G5, "")

# --- 16. Girth ---
def girth(A, n):
    for length in range(2, n+1):
        for start in range(n):
            # BFS for shortest cycle through start
            visited = {(start, frozenset([start]))}
            queue = [(start, frozenset([start]))]
            for _ in range(length - 1):
                next_queue = []
                for v, path_set in queue:
                    for w in range(n):
                        if A[v][w] and w not in path_set:
                            next_queue.append((w, path_set | {w}))
                queue = next_queue
            for v, path_set in queue:
                if A[v][start] and len(path_set) == length:
                    return length
    return float('inf')

girth_T5 = girth(T5, n_t)
girth_G5 = girth(G5, N)
add_result("Girth (directed)", girth_T5, girth_G5, "Shortest directed cycle")

# --- Print table ---
print(f"  {'Invariant':>25s} | {'T_5 (regular)':>20s} | {'G_5 (iso class)':>20s} | Interpretation")
print(f"  {'-'*25:>25s}-+-{'-'*20:>20s}-+-{'-'*20:>20s}-+-{'-'*30}")
for name, vt, vg, interp in results:
    vt_str = str(vt)[:20]
    vg_str = str(vg)[:20]
    print(f"  {name:>25s} | {vt_str:>20s} | {vg_str:>20s} | {interp}")

# ================================================================
# THE COMPARISON ANALYSIS
# ================================================================
print(f"\n{'='*70}")
print(f"  WHAT THE PARALLEL EVALUATION TEACHES")
print(f"{'='*70}\n")

# The OCF on both
I_T5_2 = sum(ip_T5[k] * 2**k for k in range(len(ip_T5)))
I_G5_2 = sum(ip_G5[k] * 2**k for k in range(len(ip_G5)))

print(f"  T_5: H = {H_T5}, I(Omega(T_5), 2) = ... (need conflict graph)")
print(f"  G_5: H = {H_G5}, I(G_5, 2) = {I_G5_2}")
print(f"  Note: I(G_5, 2) counts independent sets of ISO CLASSES weighted by 2^k")
print()

# Does OCF hold for G_5 treated as a digraph?
print(f"  OCF CHECK ON G_5 AS DIGRAPH:")
print(f"  G_5 is a bidirected graph (each undirected edge = 2 directed arcs)")
print(f"  H(G_5) = {H_G5} (directed Hamiltonian paths in the bidirected graph)")
print(f"  This counts paths that can go in either direction along each edge.")
print()

# The key ratios
if isinstance(H_G5, int) and H_G5 > 0:
    print(f"  KEY RATIOS:")
    print(f"  H(G_5) / H(T_5) = {H_G5} / {H_T5} = {H_G5/H_T5:.2f}")
    print(f"  H(G_5) / 12! = {H_G5} / 479001600 = {H_G5/479001600:.6e}")
    print(f"  H(T_5) / 5! = {H_T5} / 120 = {H_T5/120:.4f}")
    print()

    # How does G_5 compare to a random 12-vertex tournament?
    # E[H] at n=12 = 12!/2^11 = 233,887.5
    E_H_12 = 479001600 / 2048
    print(f"  E[H] at n=12 (random tournament) = {E_H_12:.1f}")
    print(f"  H(G_5) / E[H](12) = {H_G5/E_H_12:.4f}")
    print(f"  G_5 has {'MORE' if H_G5 > E_H_12 else 'FEWER'} Ham paths than a random 12-tournament")

# The independence polynomials compared
print(f"\n  INDEPENDENCE POLYNOMIALS:")
print(f"  I(T_5, x) = {' + '.join(f'{c}x^{k}' for k, c in enumerate(ip_T5) if c > 0)}")
print(f"  I(G_5, x) = {' + '.join(f'{c}x^{k}' for k, c in enumerate(ip_G5) if c > 0)}")
print(f"  I(T_5, 2) = {I_T5_2}")
print(f"  I(G_5, 2) = {I_G5_2}")

# Real roots check
roots_T5 = np.roots(list(reversed([c for c in ip_T5 if c > 0] or [1])))
roots_G5 = np.roots(list(reversed([c for c in ip_G5 if c > 0] or [1])))
all_real_T5 = all(abs(r.imag) < 1e-6 for r in roots_T5) if len(roots_T5) > 0 else True
all_real_G5 = all(abs(r.imag) < 1e-6 for r in roots_G5) if len(roots_G5) > 0 else True
print(f"  I(T_5, x) has all real roots: {all_real_T5}")
print(f"  I(G_5, x) has all real roots: {all_real_G5}")

print(f"\n  SUMMARY OF PARALLELS AND CONTRASTS:")
print(f"  T_5 is a tournament (complete, directed, n=5)")
print(f"  G_5 is a graph (sparse, undirected->bidirected, n=12)")
print()
print(f"  T_5: H=15 (max possible). Regular. Strongly connected.")
print(f"  G_5: H={H_G5}. Degree sequence {list(G5_degrees)}.")
print(f"        Connected: {sc_G5}. Kings: {kings_G5}.")
print()
print(f"  THE DEEP PARALLEL:")
print(f"  T_5 is the MAXIMIZER within 5-tournaments.")
print(f"  G_5 is the STRUCTURE of the space of 5-tournaments.")
print(f"  Evaluating G_5 with tournament invariants tells us")
print(f"  about the TOPOLOGY of tournament space itself.")
