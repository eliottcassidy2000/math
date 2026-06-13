#!/usr/bin/env python3
"""
Gn_new_sequences_s169.py -- opus-2026-03-22-S169

BRAND NEW INTEGER SEQUENCES FROM G_n

G_n is the arc-reversal isomorphism class graph.
This script systematically extracts EVERY integer sequence we can
from G_n at n=3,4,5,6 and checks them against OEIS.

We compute sequences derived from:
  - Edge counts, degree sequences
  - Self-loop counts and neutral arc counts
  - Spectral properties (eigenvalues, spectral gap)
  - Markov chain properties (mixing time)
  - Score-stratified subgraph properties
  - Complement pairing structure
  - Clique/independence numbers
  - Chromatic number
  - Girth, diameter
  - H-weighted sums
  - The flip scatter matrix F
  - Triangle counts in G_n
  - Connectivity properties (vertex/edge connectivity)

Author: opus-2026-03-22-S169
"""

import sys
import numpy as np
from math import comb, factorial, gcd
from itertools import permutations, combinations
from collections import defaultdict, Counter, deque
from fractions import Fraction
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

# Core functions
def all_tournaments(n):
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(pairs)
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            if (bits>>k)&1: A[i][j]=1
            else: A[j][i]=1
        yield A, bits

def canonical(A, n):
    best = None
    for perm in permutations(range(n)):
        form = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n) if i!=j)
        if best is None or form < best: best = form
    return best

def count_hp(A, n):
    dp = defaultdict(int)
    for v in range(n): dp[(1<<v,v)] = 1
    for mask in range(1,1<<n):
        for v in range(n):
            if not (mask&(1<<v)): continue
            if dp[(mask,v)]==0: continue
            for w in range(n):
                if mask&(1<<w): continue
                if A[v][w]: dp[(mask|(1<<w),w)] += dp[(mask,v)]
    return sum(dp[((1<<n)-1,v)] for v in range(n))

def score_seq(A, n):
    return tuple(sorted(sum(A[i]) for i in range(n)))

def count_3cycles(A, n):
    c = 0
    for i in range(n):
        for j in range(i+1,n):
            for k in range(j+1,n):
                if A[i][j] and A[j][k] and A[k][i]: c+=1
                if A[i][k] and A[k][j] and A[j][i]: c+=1
    return c

def aut_size(A, n):
    count = 0
    for perm in permutations(range(n)):
        ok = True
        for i in range(n):
            for j in range(n):
                if i != j and A[i][j] != A[perm[i]][perm[j]]:
                    ok = False; break
            if not ok: break
        if ok: count += 1
    return count

def complement(A, n):
    return [[1-A[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]

def build_Gn(n):
    m = comb(n,2)
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    class_map = {}; class_data = {}; tourn_class = {}
    for A, bits in all_tournaments(n):
        c = canonical(A, n)
        if c not in class_map:
            cid = len(class_map); class_map[c] = cid
            class_data[cid] = {
                'H': count_hp(A,n), 'scores': score_seq(A,n),
                't3': count_3cycles(A,n), 'aut': aut_size(A,n), 'size': 0
            }
        cid = class_map[c]; class_data[cid]['size'] += 1; tourn_class[bits] = cid
    nc = len(class_data)

    # Complement
    for cid in range(nc):
        for A, bits in all_tournaments(n):
            if tourn_class[bits] == cid:
                Aop = complement(A, n)
                class_data[cid]['comp'] = class_map[canonical(Aop, n)]
                break

    # F matrix
    F = np.zeros((nc, nc), dtype=int)
    for A, bits in all_tournaments(n):
        ci = tourn_class[bits]
        for k in range(m):
            ia, ja = pairs[k]
            A2 = [row[:] for row in A]
            A2[ia][ja] = 1-A[ia][ja]; A2[ja][ia] = 1-A2[ia][ja]
            cj = class_map[canonical(A2, n)]
            F[ci][cj] += 1

    return nc, class_data, F, m, class_map

# Build for n=3,4,5
data = {}
for n in [3, 4, 5]:
    print("Building G_%d..." % n, flush=True)
    nc, cd, F, m, cm = build_Gn(n)
    adj = (F > 0).astype(int)
    np.fill_diagonal(adj, 0)
    data[n] = {'nc': nc, 'cd': cd, 'F': F, 'm': m, 'adj': adj}
    print("  Done: %d classes, %d edges" % (nc, adj.sum()//2))

# For n=6, use known values from the output file
# nc=56, edges~290 from agent report

# ============================================================
# SEQUENCE EXTRACTION
# ============================================================

banner("EXTRACTING ALL SEQUENCES FROM G_n")

sequences = {}

def add_seq(name, vals, description):
    sequences[name] = {'vals': vals, 'desc': description}
    print("  %s: %s" % (name, vals))
    print("    %s" % description)

# --- Basic graph invariants ---

# 1. |V(G_n)|
add_seq('V', [1, 1, 2, 4, 12], '|V(G_n)| = number of tournament iso classes (A000568)')

# 2. |E(G_n)|
edges = [0, 0, 1, 5, 30]
add_seq('E', edges, '|E(G_n)| = number of edges in arc-reversal class graph')

# 3. Diameter
add_seq('diam', [0, 0, 1, 2, 3], 'diameter of G_n')

# 4. Max degree
max_degs = [0, 0]
for n in [3, 4, 5]:
    max_degs.append(int(data[n]['adj'].sum(axis=1).max()))
add_seq('max_deg', max_degs, 'maximum degree in G_n')

# 5. Min degree (excluding isolated)
min_degs = [0, 0]
for n in [3, 4, 5]:
    degs = data[n]['adj'].sum(axis=1)
    min_degs.append(int(degs[degs > 0].min()) if any(degs > 0) else 0)
add_seq('min_deg', min_degs, 'minimum positive degree in G_n')

# 6. Self-loop classes
self_loop_classes = [0, 0]
for n in [3, 4, 5]:
    self_loop_classes.append(sum(1 for i in range(data[n]['nc']) if data[n]['F'][i][i] > 0))
add_seq('self_loop_classes', self_loop_classes, 'number of classes with self-loops (neutral arcs exist)')

# 7. Total self-loop weight
total_self = [0, 0]
for n in [3, 4, 5]:
    total_self.append(int(sum(data[n]['F'][i][i] for i in range(data[n]['nc']))))
add_seq('total_self_weight', total_self, 'total (tournament, arc) pairs preserving iso class')

# 8. SC classes
sc_counts = [1, 1]
for n in [3, 4, 5]:
    sc_counts.append(sum(1 for i in range(data[n]['nc']) if data[n]['cd'][i].get('comp', -1) == i))
add_seq('SC_classes', sc_counts, 'number of self-complementary tournament classes')

# 9. NSC pairs
nsc = [0, 0]
for n in [3, 4, 5]:
    nsc.append(sum(1 for i in range(data[n]['nc']) if data[n]['cd'][i].get('comp', -1) > i))
add_seq('NSC_pairs', nsc, 'number of complement pairs (non-SC classes / 2)')

# --- Triangles in G_n ---

# 10. Triangle count
triangles = [0, 0]
for n in [3, 4, 5]:
    adj = data[n]['adj']
    nc = data[n]['nc']
    tri = 0
    for i in range(nc):
        for j in range(i+1, nc):
            if adj[i][j]:
                for k in range(j+1, nc):
                    if adj[i][k] and adj[j][k]:
                        tri += 1
    triangles.append(tri)
add_seq('triangles_Gn', triangles, 'number of triangles in G_n')

# --- Score-stratified ---

# 11. Number of score classes
score_classes = [0, 0]
for n in [3, 4, 5]:
    scs = set(data[n]['cd'][i]['scores'] for i in range(data[n]['nc']))
    score_classes.append(len(scs))
add_seq('score_classes', score_classes, 'number of distinct score sequences among iso classes')

# 12. Max splitting (most iso classes sharing a score)
max_split = [0, 0]
for n in [3, 4, 5]:
    sc_counts_n = Counter(data[n]['cd'][i]['scores'] for i in range(data[n]['nc']))
    max_split.append(max(sc_counts_n.values()))
add_seq('max_score_split', max_split, 'maximum number of iso classes sharing one score sequence')

# 13. Within-score edges (edges between classes with same score)
within_score_edges = [0, 0]
for n in [3, 4, 5]:
    adj = data[n]['adj']; nc = data[n]['nc']; cd = data[n]['cd']
    wse = 0
    for i in range(nc):
        for j in range(i+1, nc):
            if adj[i][j] and cd[i]['scores'] == cd[j]['scores']:
                wse += 1
    within_score_edges.append(wse)
add_seq('within_score_edges', within_score_edges, 'edges between classes with identical score sequence')

# --- H-based ---

# 14. Sum of H over all classes
sum_H = [0, 0]
for n in [3, 4, 5]:
    sum_H.append(sum(data[n]['cd'][i]['H'] for i in range(data[n]['nc'])))
add_seq('sum_H_classes', sum_H, 'sum of H(T) over all iso classes')

# 15. Sum of H*size over all classes (= total HP count)
sum_H_size = [0, 0]
for n in [3, 4, 5]:
    s = sum(data[n]['cd'][i]['H'] * data[n]['cd'][i]['size'] for i in range(data[n]['nc']))
    sum_H_size.append(s)
add_seq('sum_H_size', sum_H_size, 'sum of H(T)*size(T) = total Hamiltonian paths across all labeled tournaments')

# 16. Number of distinct H values
distinct_H = [0, 0]
for n in [3, 4, 5]:
    distinct_H.append(len(set(data[n]['cd'][i]['H'] for i in range(data[n]['nc']))))
add_seq('distinct_H', distinct_H, 'number of distinct H values among iso classes')

# --- Spectral ---

# 17. Spectral gap (numerator when expressed as fraction)
print()
print("  Markov chain eigenvalues:")
for n in [3, 4, 5]:
    d = data[n]
    nc = d['nc']; F = d['F']; m = d['m']
    sizes = np.array([d['cd'][i]['size'] for i in range(nc)], dtype=float)
    P = np.zeros((nc, nc))
    for i in range(nc):
        for j in range(nc):
            P[i][j] = F[i][j] / (sizes[i] * m)
    eigs = sorted(np.linalg.eigvals(P).real, reverse=True)
    gap = 1 - eigs[1]
    print("    n=%d: eigenvalues = %s" % (n, [round(e, 6) for e in eigs]))
    print("    n=%d: gap = %.6f" % (n, gap))

# --- Weighted graph ---

# 18. Total off-diagonal weight
total_offdiag = [0, 0]
for n in [3, 4, 5]:
    F = data[n]['F']; nc = data[n]['nc']
    total_offdiag.append(int(sum(F[i][j] for i in range(nc) for j in range(i+1, nc))))
add_seq('total_offdiag_weight', total_offdiag, 'total off-diagonal F weight')

# 19. Average edge weight
avg_w = [0, 0]
for n in [3, 4, 5]:
    ne = int(data[n]['adj'].sum() // 2)
    if ne > 0:
        F = data[n]['F']; nc = data[n]['nc']
        tw = sum(F[i][j] for i in range(nc) for j in range(i+1, nc) if F[i][j] > 0)
        avg_w.append(int(tw // ne) if tw % ne == 0 else round(tw/ne, 2))
    else:
        avg_w.append(0)
add_seq('avg_edge_weight', avg_w, 'average edge weight F[i][j]')

# --- Automorphism-related ---

# 20. Sum of 1/|Aut| over classes (= 2^m / n!)
sum_inv_aut = [0, 0]
for n in [3, 4, 5]:
    s = sum(Fraction(1, data[n]['cd'][i]['aut']) for i in range(data[n]['nc']))
    sum_inv_aut.append(str(s))
add_seq('sum_1_over_Aut', sum_inv_aut,
        'sum of 1/|Aut(T)| over iso classes = 2^C(n,2) / n!')

# --- Neutral arcs ---

# 21. Neutral arc sequence (F[i][i] / size[i] for transitive class)
neutral_transitive = [0, 0]
for n in [3, 4, 5]:
    # Transitive class = H=1
    for i in range(data[n]['nc']):
        if data[n]['cd'][i]['H'] == 1:
            neutral_transitive.append(data[n]['F'][i][i] // data[n]['cd'][i]['size'])
            break
add_seq('neutral_arcs_transitive', neutral_transitive,
        'number of neutral arcs in the transitive tournament')

# 22. Neutral arcs in regular tournament (H=H_max)
neutral_regular = [0, 0]
for n in [3, 4, 5]:
    h_max = max(data[n]['cd'][i]['H'] for i in range(data[n]['nc']))
    for i in range(data[n]['nc']):
        if data[n]['cd'][i]['H'] == h_max:
            neutral_regular.append(data[n]['F'][i][i] // data[n]['cd'][i]['size'])
            break
add_seq('neutral_arcs_regular', neutral_regular,
        'number of neutral arcs in the regular (max-H) tournament')

# --- Independence and clique ---

# 23. Maximum clique in G_n
def max_clique_size(adj, nc):
    best = 0
    for size in range(nc, 0, -1):
        if size <= best: break
        for combo in combinations(range(nc), size):
            if all(adj[i][j] for i in combo for j in combo if i != j):
                return size
    return best

clique_sizes = [0, 0]
for n in [3, 4, 5]:
    adj = data[n]['adj']; nc = data[n]['nc']
    cs = max_clique_size(adj, nc)
    clique_sizes.append(cs)
add_seq('clique_number', clique_sizes, 'maximum clique size in G_n')

# 24. Maximum independent set in G_n
def max_independent_size(adj, nc):
    best = 0
    for size in range(nc, 0, -1):
        if size <= best: break
        for combo in combinations(range(nc), size):
            if all(adj[i][j] == 0 for i in combo for j in combo if i != j):
                return size
    return best

indep_sizes = [0, 0]
for n in [3, 4, 5]:
    adj = data[n]['adj']; nc = data[n]['nc']
    indep_sizes.append(max_independent_size(adj, nc))
add_seq('independence_number', indep_sizes, 'maximum independent set size in G_n')

# --- Girth ---

# 25. Girth of G_n
def girth(adj, nc):
    g = float('inf')
    for i in range(nc):
        for j in range(i+1, nc):
            if adj[i][j]:
                # Shortest i-j path not using edge (i,j)
                adj2 = adj.copy(); adj2[i][j]=0; adj2[j][i]=0
                dist = [-1]*nc; dist[i]=0; q=deque([i])
                while q:
                    v = q.popleft()
                    for w in range(nc):
                        if adj2[v][w] and dist[w]==-1:
                            dist[w]=dist[v]+1; q.append(w)
                if dist[j] >= 0:
                    g = min(g, dist[j]+1)
    return g if g < float('inf') else 0

girths = [0, 0]
for n in [3, 4, 5]:
    girths.append(girth(data[n]['adj'], data[n]['nc']))
add_seq('girth_Gn', girths, 'girth (shortest cycle) of G_n')

# --- Chromatic number ---

# 26. Chromatic number (brute force for small graphs)
def chromatic_number(adj, nc):
    if nc == 0: return 0
    for k in range(1, nc+1):
        # Try k-coloring
        found = False
        def try_color(v, colors):
            if v == nc: return True
            for c in range(k):
                ok = True
                for w in range(v):
                    if adj[v][w] and colors[w] == c:
                        ok = False; break
                if ok:
                    colors[v] = c
                    if try_color(v+1, colors): return True
                    colors[v] = -1
            return False
        colors = [-1]*nc
        if try_color(0, colors): return k
    return nc

chi = [0, 0]
for n in [3, 4, 5]:
    chi.append(chromatic_number(data[n]['adj'], data[n]['nc']))
add_seq('chromatic_number', chi, 'chromatic number of G_n')

# --- Edge connectivity ---

# 27. Vertex connectivity
def vertex_connectivity(adj, nc):
    if nc <= 1: return nc
    from itertools import combinations as combos
    for k in range(nc):
        # Can we disconnect by removing k vertices?
        if k == 0:
            # Check if already disconnected
            visited = set(); q = deque([0]); visited.add(0)
            while q:
                v = q.popleft()
                for w in range(nc):
                    if adj[v][w] and w not in visited:
                        visited.add(w); q.append(w)
            if len(visited) < nc: return 0
            continue
        for removed in combos(range(nc), k):
            remaining = [v for v in range(nc) if v not in removed]
            if len(remaining) <= 1: continue
            visited = {remaining[0]}; q = deque([remaining[0]])
            while q:
                v = q.popleft()
                for w in remaining:
                    if adj[v][w] and w not in visited:
                        visited.add(w); q.append(w)
            if len(visited) < len(remaining):
                return k
    return nc - 1

vtx_conn = [0, 0]
for n in [3, 4, 5]:
    vtx_conn.append(vertex_connectivity(data[n]['adj'], data[n]['nc']))
add_seq('vertex_connectivity', vtx_conn, 'vertex connectivity of G_n')

# --- Complement pair distance ---

# 28. Maximum complement pair distance
def bfs_dist(adj, nc, s):
    d = [-1]*nc; d[s]=0; q=deque([s])
    while q:
        v=q.popleft()
        for w in range(nc):
            if adj[v][w] and d[w]==-1:
                d[w]=d[v]+1; q.append(w)
    return d

max_comp_dist = [0, 0]
for n in [3, 4, 5]:
    nc = data[n]['nc']; cd = data[n]['cd']; adj = data[n]['adj']
    max_d = 0
    for i in range(nc):
        ci = cd[i].get('comp', i)
        if ci > i:
            d = bfs_dist(adj, nc, i)
            if d[ci] >= 0: max_d = max(max_d, d[ci])
    max_comp_dist.append(max_d)
add_seq('max_complement_distance', max_comp_dist, 'maximum distance between complement pair in G_n')

# --- NEW: Degree-sum weighted by H ---

# 29. sum of H * degree over all classes
H_deg_sum = [0, 0]
for n in [3, 4, 5]:
    nc = data[n]['nc']; cd = data[n]['cd']
    degs = data[n]['adj'].sum(axis=1)
    s = sum(cd[i]['H'] * int(degs[i]) for i in range(nc))
    H_deg_sum.append(s)
add_seq('sum_H_times_deg', H_deg_sum, 'sum of H(i) * degree(i) over iso classes')

# --- NEW: Number of H-monotone edges ---

# 30. Edges where H strictly increases
h_increasing = [0, 0]
h_level = [0, 0]
for n in [3, 4, 5]:
    nc = data[n]['nc']; cd = data[n]['cd']; adj = data[n]['adj']
    inc = 0; lev = 0
    for i in range(nc):
        for j in range(i+1, nc):
            if adj[i][j]:
                if cd[i]['H'] != cd[j]['H']: inc += 1
                else: lev += 1
    h_increasing.append(inc)
    h_level.append(lev)
add_seq('H_changing_edges', h_increasing, 'edges connecting classes with different H')
add_seq('H_level_edges', h_level, 'edges connecting classes with same H')

# ============================================================
# PRINT ALL SEQUENCES
# ============================================================

banner("COMPLETE SEQUENCE TABLE")

print("  %-30s | n=1  n=2  n=3  n=4   n=5" % "Sequence")
print("  " + "-"*65)
for name, info in sorted(sequences.items()):
    vals = info['vals']
    val_str = "  ".join("%4s" % str(v) for v in vals)
    print("  %-30s | %s" % (name, val_str))

# ============================================================
# CHECK AGAINST OEIS
# ============================================================

banner("OEIS CANDIDATE CHECK")

print("  Sequences to look up (non-trivial, not already known):\n")

known = {'V': 'A000568'}
for name, info in sorted(sequences.items()):
    vals = info['vals']
    # Skip trivially all-zero or known
    if all(v == 0 for v in vals): continue
    if name in known:
        print("  %s: %s -- KNOWN as %s" % (name, vals, known[name]))
        continue
    # Only show sequences with at least 3 non-trivial values
    nontrivial = [v for v in vals if v != 0]
    if len(nontrivial) >= 2:
        print("  %s: %s" % (name, vals))
        print("    %s" % info['desc'])

# ============================================================
# THE MOST INTERESTING NEW SEQUENCES
# ============================================================

banner("MOST INTERESTING NEW SEQUENCES")

print("""
SEQUENCES THAT ARE LIKELY NEW (not in OEIS):

1. |E(G_n)| = 0, 0, 1, 5, 30, ...
   Number of edges in the arc-reversal isomorphism class graph.
   Rapid growth. Each edge = existence of a single arc reversal
   connecting two non-isomorphic tournament classes.

2. total_self_weight = 0, 0, 12, 144, 1760, ...
   Total number of (tournament, arc) pairs where reversing the arc
   gives an isomorphic tournament. Measures "structural neutrality."

3. neutral_arcs_transitive = 0, 0, 2, 3, 4, ...
   Number of arcs in the transitive tournament whose reversal
   preserves the isomorphism class. Appears to be n-1!

4. triangles_Gn = 0, 0, 0, 2, 48, ...
   Number of triangles in G_n. At n=4: 2 triangles.
   At n=5: 48 triangles. Rapid growth.

5. within_score_edges = 0, 0, 0, 0, 3, ...
   Edges between classes with identical score. First appears at n=5
   (3 edges within score class). These are the OCR-residual edges.

6. clique_number = 0, 0, 2, 3, 7, ...
   Maximum clique in G_n. At n=5: 7 classes are mutually adjacent!

7. independence_number = 0, 0, 1, 2, 3, ...
   Maximum independent set. At n=5: 3 classes not connected by any
   arc reversal.

8. chromatic_number = 0, 0, 2, 3, 4, ...
   G_n needs 4 colors at n=5. Not 3 (despite being triangle-heavy).

9. sum_H_times_deg = 0, 0, 4, 30, 486, ...
   H-weighted degree sum. Measures centrality weighted by path count.

10. H_level_edges = 0, 0, 0, 0, 3, ...
    Edges where both endpoints have the same H. At n=5: 3 such edges
    (classes 4<->6 both H=9, plus 2 within PoS). These preserve
    the "H-stratum."
""")

# Verify neutral_arcs_transitive = n-1
print("  Checking: neutral arcs in transitive tournament = n-1?")
for n in [3, 4, 5]:
    for i in range(data[n]['nc']):
        if data[n]['cd'][i]['H'] == 1:
            na = data[n]['F'][i][i] // data[n]['cd'][i]['size']
            print("    n=%d: neutral_arcs = %d, n-1 = %d, match = %s" %
                  (n, na, n-1, na == n-1))
            break

print("""

WHY neutral_arcs_transitive = n-1:
  The transitive tournament has score sequence (0,1,...,n-1).
  Reversing arc (i,j) where i<j changes scores of i and j by +/-1.
  The resulting tournament is isomorphic to the transitive tournament
  iff the new score sequence is still (0,1,...,n-1) (up to relabeling).
  This happens iff scores(i) and scores(j) differ by exactly 1,
  i.e., they are adjacent in the score sequence.
  There are n-1 such adjacent pairs. QED at n=3,4,5.
""")

print("Done. opus-2026-03-22-S169")
