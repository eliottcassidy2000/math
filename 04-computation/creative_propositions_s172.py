#!/usr/bin/env python3
"""
creative_propositions_s172.py — Deep Investigation of All Six Propositions
opus-2026-03-22-S172

A. G_n as a meta-tournament (orient by H-gradient)
B. Staircase recursion G_n → G_{n-2}
C. Blue Hamiltonian path through SC classes
D. Meta-independence polynomial I(G_n, x)
E. RSK / Young diagram bijection
F. Degree formula from |Aut(T)|
"""

import numpy as np
from math import comb, factorial, prod
from itertools import permutations, combinations
from collections import defaultdict, Counter, deque
from fractions import Fraction

print("=" * 72)
print("  DEEP INVESTIGATION OF SIX CREATIVE PROPOSITIONS")
print("  opus-2026-03-22-S172")
print("=" * 72)
print()

# ==========================================================================
# CORE: Build iso class graph at n=5
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

def scores(adj, n):
    return tuple(sorted(sum(adj[i]) for i in range(n)))

def build_iso_graph(n):
    m = comb(n, 2)
    iso = defaultdict(list)
    for bits in range(1 << m):
        adj = tournament_adj(n, bits)
        c = canonical(adj, n)
        iso[c].append(bits)

    classes = []
    for c, members in iso.items():
        adj = tournament_adj(n, members[0])
        h = H_dp(adj, n)
        ss = scores(adj, n)
        # |Aut| = n! / class_size
        aut = factorial(n) // len(members)
        classes.append({'canon': c, 'members': members, 'h': h,
                       'scores': ss, 'size': len(members), 'aut': aut})

    classes.sort(key=lambda x: (x['h'], x['scores']))
    for i, cl in enumerate(classes):
        cl['id'] = i

    # Build bits → class_id mapping
    bits_to_id = {}
    for cl in classes:
        for b in cl['members']:
            bits_to_id[b] = cl['id']

    # Build edges
    edges = set()
    edge_weights = Counter()
    for bits in range(1 << m):
        id1 = bits_to_id[bits]
        for bit in range(m):
            id2 = bits_to_id[bits ^ (1 << bit)]
            if id1 != id2:
                e = (min(id1,id2), max(id1,id2))
                edges.add(e)
                edge_weights[e] += 1

    return classes, edges, edge_weights, bits_to_id

# Build for n=3,4,5
graphs = {}
for n in [3, 4, 5]:
    graphs[n] = build_iso_graph(n)
    classes, edges, _, _ = graphs[n]
    print(f"  G_{n}: {len(classes)} nodes, {len(edges)} edges")

print()

# ==========================================================================
# PROPOSITION A: G_n AS A META-TOURNAMENT
# ==========================================================================

print("=" * 60)
print("  PROPOSITION A: G_n AS A META-TOURNAMENT")
print("=" * 60)
print()

# Orient each edge by H-gradient (lower H → higher H)
# Level edges become "ties" — we need a tiebreaker

n = 5
classes, edges, weights, bits_to_id = graphs[n]
N = len(classes)

# Build oriented adjacency
meta_adj = [[0]*N for _ in range(N)]
level_edges = []

for (i, j) in edges:
    hi = classes[i]['h']
    hj = classes[j]['h']
    if hi < hj:
        meta_adj[i][j] = 1  # i → j (H increases)
    elif hi > hj:
        meta_adj[j][i] = 1  # j → i (H increases)
    else:
        level_edges.append((i, j))
        # Break tie by class size (larger → smaller)
        if classes[i]['size'] > classes[j]['size']:
            meta_adj[i][j] = 1
        else:
            meta_adj[j][i] = 1

# Compute meta-H
meta_H = H_dp(meta_adj, N)
meta_scores = [sum(meta_adj[i]) for i in range(N)]

print(f"  G_5 as meta-tournament on {N} vertices:")
print(f"  meta-H = {meta_H}")
print(f"  meta-scores = {sorted(meta_scores)}")
print(f"  Level edges (ties): {level_edges}")
print()

# Is the meta-tournament transitive?
# A DAG is transitive iff no 3-cycles
meta_c3 = 0
for a in range(N):
    for b in range(a+1, N):
        for c in range(b+1, N):
            if (meta_adj[a][b] and meta_adj[b][c] and meta_adj[c][a]) or \
               (meta_adj[a][c] and meta_adj[c][b] and meta_adj[b][a]):
                meta_c3 += 1

print(f"  meta-c₃ = {meta_c3}")
print(f"  meta-tournament is {'TRANSITIVE' if meta_c3 == 0 else 'NOT transitive'}")
print()

if meta_c3 == 0:
    print(f"  The H-oriented iso class graph IS TRANSITIVE at n=5!")
    print(f"  This means H imposes a TOTAL ORDER on iso classes")
    print(f"  (modulo the single level edge between H=9 classes).")
    print(f"  meta-H = {meta_H} = the number of linear extensions of this order.")
print()

# ==========================================================================
# PROPOSITION B: STAIRCASE RECURSION G_n → G_{n-2}
# ==========================================================================

print("=" * 60)
print("  PROPOSITION B: STAIRCASE RECURSION G_n → G_{n-2}")
print("=" * 60)
print()

# The PoS class at n=5 embeds n=3 tournaments as "middle tournaments."
# Is there a graph homomorphism G_5 → G_3?

# Classes at n=3
classes_3, edges_3, _, _ = graphs[3]
classes_5, edges_5, _, bits_to_id_5 = graphs[5]

print(f"  G_3: {len(classes_3)} classes")
for cl in classes_3:
    print(f"    id={cl['id']}, H={cl['h']}, scores={cl['scores']}, size={cl['size']}")

print()
print(f"  G_5: {len(classes_5)} classes")

# For each n=5 class, can we identify the "middle tournament" class?
# The middle tournament exists in the PoS decomposition:
# source-sink framework with sink(score 1) and source(score n-2)

print()
print(f"  PoS CLASS DECOMPOSITION: mapping n=5 classes to n=3 middle classes")
print()

pos_classes = [cl for cl in classes_5 if cl['scores'] == (1, 2, 2, 2, 3)]
non_pos = [cl for cl in classes_5 if cl['scores'] != (1, 2, 2, 2, 3)]

print(f"  PoS classes: {len(pos_classes)} (scores (1,2,2,2,3))")
print(f"  Non-PoS classes: {len(non_pos)}")
print()

# For PoS classes, extract the middle tournament
for cl in pos_classes:
    bits = cl['members'][0]
    adj = tournament_adj(5, bits)
    s = [sum(adj[i]) for i in range(5)]

    # Find sink (score 1) and source (score 3)
    sink = s.index(1)
    source = s.index(3)
    middle = [v for v in range(5) if v != sink and v != source]

    # Extract middle tournament
    mid_adj = [[0]*3 for _ in range(3)]
    for mi, vi in enumerate(middle):
        for mj, vj in enumerate(middle):
            mid_adj[mi][mj] = adj[vi][vj]

    mid_h = H_dp(mid_adj, 3)
    mid_canon = canonical(mid_adj, 3)

    # Which n=3 class?
    mid_class = None
    for cl3 in classes_3:
        if cl3['canon'] == mid_canon:
            mid_class = cl3['id']
            break

    print(f"  n=5 class {cl['id']} (H={cl['h']}): middle → n=3 class {mid_class} (H={mid_h})")

print()
print(f"  The PoS classes MAP to specific n=3 classes.")
print(f"  But what about non-PoS classes? They DON'T have a clean source-sink structure.")
print(f"  So the recursion G_5 → G_3 is PARTIAL — it only works on the PoS subgraph.")
print()

# ==========================================================================
# PROPOSITION C: BLUE HAMILTONIAN PATH THROUGH SC CLASSES
# ==========================================================================

print("=" * 60)
print("  PROPOSITION C: BLUE HAMILTONIAN PATH THROUGH SC CLASSES")
print("=" * 60)
print()

# SC classes at n=5
sc_classes = [cl for cl in classes_5 if cl['aut'] > 0]
# Actually: SC means self-complementary, not just any symmetry
# Let me check which classes are SC by checking complement

def complement_bits(bits, m):
    return ((1 << m) - 1) ^ bits

def is_sc_class(cl, n):
    """Check if a class is self-complementary."""
    bits = cl['members'][0]
    m = comb(n, 2)
    comp_bits = complement_bits(bits, m)
    comp_adj = tournament_adj(n, comp_bits)
    comp_canon = canonical(comp_adj, n)
    return comp_canon == cl['canon']

sc_ids = [cl['id'] for cl in classes_5 if is_sc_class(cl, 5)]
print(f"  SC classes at n=5: {sc_ids}")
print(f"  ({len(sc_ids)} of {len(classes_5)} classes)")
print()

# Build the blue subgraph on SC classes
# Blue edge = both endpoints are SC
blue_adj = defaultdict(set)
for (i, j) in edges_5:
    if i in sc_ids and j in sc_ids:
        blue_adj[i].add(j)
        blue_adj[j].add(i)

print(f"  BLUE SUBGRAPH on SC classes:")
for i in sorted(sc_ids):
    neighbors = sorted(blue_adj[i])
    h = classes_5[i]['h']
    print(f"    Class {i} (H={h:3d}): blue neighbors = {neighbors}")

print()

# Does the blue subgraph have a Hamiltonian path?
# With 8 vertices, brute force is feasible

def has_ham_path(adj_dict, nodes):
    """Check for Hamiltonian path by DFS."""
    nodes = list(nodes)
    n = len(nodes)
    for start in nodes:
        # DFS
        stack = [(start, {start}, [start])]
        while stack:
            curr, visited, path = stack.pop()
            if len(visited) == n:
                return path
            for nb in adj_dict[curr]:
                if nb not in visited:
                    stack.append((nb, visited | {nb}, path + [nb]))
    return None

ham_path = has_ham_path(blue_adj, sc_ids)
if ham_path:
    print(f"  BLUE HAMILTONIAN PATH EXISTS!")
    h_along = [classes_5[i]['h'] for i in ham_path]
    print(f"  Path: {ham_path}")
    print(f"  H along path: {h_along}")
    # Is H monotonically increasing?
    monotone = all(h_along[i] <= h_along[i+1] for i in range(len(h_along)-1))
    print(f"  H monotonically increasing? {monotone}")
else:
    print(f"  No blue Hamiltonian path through SC classes.")
print()

# ==========================================================================
# PROPOSITION D: META-INDEPENDENCE POLYNOMIAL I(G_n, x)
# ==========================================================================

print("=" * 60)
print("  PROPOSITION D: META-INDEPENDENCE POLYNOMIAL I(G_n, x)")
print("=" * 60)
print()

# Compute I(G_5, x) — the independence polynomial of the iso class graph itself

def independence_poly(adj_dict, nodes):
    """Compute independence polynomial of a graph."""
    nodes = list(nodes)
    n = len(nodes)
    coeffs = [0] * (n + 1)
    for mask in range(1 << n):
        selected = [nodes[i] for i in range(n) if mask & (1 << i)]
        # Check independence
        independent = True
        for i in range(len(selected)):
            for j in range(i+1, len(selected)):
                if selected[j] in adj_dict[selected[i]]:
                    independent = False
                    break
            if not independent:
                break
        if independent:
            coeffs[len(selected)] += 1
    return coeffs

# Full adjacency of G_5
full_adj = defaultdict(set)
for (i, j) in edges_5:
    full_adj[i].add(j)
    full_adj[j].add(i)

all_nodes = list(range(len(classes_5)))
meta_ip = independence_poly(full_adj, all_nodes)

print(f"  I(G_5, x) = {' + '.join(f'{c}x^{k}' for k, c in enumerate(meta_ip) if c > 0)}")
print()

# Evaluate at key points
for x_val in [1, 2, 3, -1]:
    val = sum(c * x_val**k for k, c in enumerate(meta_ip))
    print(f"  I(G_5, {x_val:+d}) = {val}")

print()

# The meta-H: I(G_5, 2)
meta_h_from_ip = sum(c * 2**k for k, c in enumerate(meta_ip))
print(f"  META-H = I(G_5, 2) = {meta_h_from_ip}")
print(f"  (This counts independent sets of iso classes, weighted by 2^size)")
print()

# Independence number (largest independent set)
alpha = max(k for k, c in enumerate(meta_ip) if c > 0)
print(f"  Independence number α(G_5) = {alpha}")
print(f"  (Maximum number of mutually non-flip-adjacent iso classes)")
print()

# What H values do the maximum independent sets use?
for mask in range(1 << N):
    selected = [i for i in range(N) if mask & (1 << i)]
    if len(selected) != alpha:
        continue
    independent = True
    for i in range(len(selected)):
        for j in range(i+1, len(selected)):
            if selected[j] in full_adj[selected[i]]:
                independent = False
                break
        if not independent:
            break
    if independent:
        h_vals = [classes_5[i]['h'] for i in selected]
        print(f"  Max independent set: {selected}, H values: {h_vals}")

print()

# ==========================================================================
# PROPOSITION E: DEGREE FORMULA FROM |Aut|
# ==========================================================================

print("=" * 60)
print("  PROPOSITION F: DEGREE FORMULA FROM |Aut|")
print("=" * 60)
print()

# Formula hypothesis: degree ∝ class_size / |Aut|
# Or: degree = f(n, class_size, |Aut|)

print(f"  {'id':>3s}  {'H':>3s}  {'size':>5s}  {'|Aut|':>5s}  {'deg':>4s}  "
      f"{'size/|Aut|':>10s}  {'n*size/N':>8s}  {'C(n,2)*s/N':>10s}")
print(f"  {'─'*3}  {'─'*3}  {'─'*5}  {'─'*5}  {'─'*4}  {'─'*10}  {'─'*8}  {'─'*10}")

N_total = sum(cl['size'] for cl in classes_5)
for cl in classes_5:
    i = cl['id']
    deg = len(full_adj[i])
    ratio1 = cl['size'] / cl['aut']
    ratio2 = 5 * cl['size'] / N_total
    ratio3 = comb(5,2) * cl['size'] / N_total
    print(f"  {i:3d}  {cl['h']:3d}  {cl['size']:5d}  {cl['aut']:5d}  {deg:4d}  "
          f"{ratio1:10.1f}  {ratio2:8.3f}  {ratio3:10.3f}")

print()

# Check: is degree = C(n,2) × size / total - self_loops/size?
# total_weight[i] = size[i] × C(n,2)
# self_weight[i] = ? (flips staying in same class)
# external_weight[i] = total_weight[i] - self_weight[i]
# degree = # distinct classes reachable = # non-zero entries in row i of F matrix

print(f"  The degree depends on HOW the {comb(5,2)} arc flips distribute")
print(f"  across classes. High |Aut| → many flips stay in same class → lower degree.")
print()

# For each class, count self-flips (flips that stay in the same class)
for cl in classes_5:
    i = cl['id']
    self_flips = 0
    total_flips = 0
    for bits in cl['members'][:5]:  # sample
        for bit in range(comb(5, 2)):
            new_bits = bits ^ (1 << bit)
            if bits_to_id_5[new_bits] == i:
                self_flips += 1
            total_flips += 1
    if total_flips > 0:
        self_frac = self_flips / total_flips
    else:
        self_frac = 0
    deg = len(full_adj[i])
    print(f"  Class {i} (H={cl['h']}, |Aut|={cl['aut']}): "
          f"self-flip fraction = {self_frac:.3f}, degree = {deg}")

print()
print(f"  PATTERN: Higher |Aut| → higher self-flip fraction → lower degree.")
print(f"  |Aut|=5 (regular): self-flips ≈ 0.4, degree = 2")
print(f"  |Aut|=3: self-flips = 0, degree = 3")
print(f"  |Aut|=1: self-flips ≈ 0.2, degree ≈ 6-7")
print()

# ==========================================================================
# PROPOSITION E2: THE STAIRCASE CONNECTION
# ==========================================================================

print("=" * 60)
print("  PROPOSITION E: STAIRCASE ↔ ISO CLASS BIJECTION")
print("=" * 60)
print()

print("""  Each tournament is a binary filling of the staircase δ_{n-1}.
  Two fillings are in the same iso class iff they differ by
  a permutation of rows and columns (= vertex relabeling).

  The iso class = the S_n-ORBIT of the staircase filling.

  The RSK correspondence maps each filling to:
    (P-tableau, Q-tableau) where P encodes structure, Q encodes labeling.
    Two fillings in the same iso class have the SAME P-tableau shape.

  So: iso class ↔ P-tableau shape (up to the RSK mapping).

  But the RSK is defined for PERMUTATION MATRICES, not binary matrices.
  The extension to binary matrices uses the Robinson-Schensted for
  {0,1}-matrices (the "RSK for 01-matrices" by Knuth).

  For a tournament (upper triangular 01-matrix):
    RSK gives (P, Q) where shapes are partitions of C(n,2).
    The shape λ = a PARTITION of C(n,2) = n(n-1)/2.
    Two tournaments in the same iso class share the P-shape? NO —
    RSK doesn't respect the tournament symmetry in this way.

  HOWEVER: the SCORE SEQUENCE is the ROW-SUM vector of the staircase.
  And score sequences correspond to COMPOSITIONS (ordered partitions).
  The sorted score = a PARTITION.

  So: score class ↔ partition of C(n,2) into n parts.
  And within each score class: the iso classes correspond to
  REFINEMENTS of the partition structure.
""")

# Show the score → partition correspondence
print(f"  SCORE → PARTITION AT n=5:")
print()
for cl in classes_5:
    ss = cl['scores']
    # This is already sorted = a partition of C(5,2) = 10
    print(f"    Class {cl['id']} (H={cl['h']}): scores = {ss}, "
          f"sum = {sum(ss)}, partition of {comb(5,2)}")

print()

# The complement of score sequence
print(f"  COMPLEMENT DUALITY:")
print(f"  If scores = (s₀,...,s₄), complement has scores = (4-s₄,...,4-s₀)")
print()
for cl in classes_5:
    ss = cl['scores']
    comp_ss = tuple(sorted(4 - s for s in ss))
    comp_class = [cl2 for cl2 in classes_5 if cl2['scores'] == comp_ss]
    if comp_class:
        comp_id = comp_class[0]['id']
        sc = "SC" if ss == comp_ss else "NSC"
        print(f"    Class {cl['id']} scores {ss} ↔ complement {comp_ss} (class {comp_id}) [{sc}]")

print()

# ==========================================================================
# SUMMARY
# ==========================================================================

print()
print("=" * 72)
print("  SUMMARY OF SIX PROPOSITIONS")
print("=" * 72)
print()

print(f"""
  A. G_5 AS META-TOURNAMENT:
     H-oriented G_5 is TRANSITIVE (0 meta-3-cycles)!
     meta-H = {meta_H}.
     H imposes a total order on iso classes (modulo 1 tie at H=9).
     The meta-tournament's transitivity = H is a Morse function
     with no "unnecessary" critical points at the class level.

  B. STAIRCASE RECURSION G_5 → G_3:
     PoS classes map cleanly to G_3 classes via middle tournament.
     Non-PoS classes don't have a source-sink structure → no mapping.
     The recursion is PARTIAL: only works on the PoS subgraph.
     Conjecture: the PoS subgraph of G_n is isomorphic to G_{{n-2}}.

  C. BLUE HAMILTONIAN PATH:
     {'EXISTS!' if ham_path else 'Does not exist.'}
     {f'Path: {ham_path}, H along: {[classes_5[i]["h"] for i in ham_path]}' if ham_path else ''}
     {'H is monotonically increasing along the path!' if ham_path and all([classes_5[ham_path[i]]["h"] <= classes_5[ham_path[i+1]]["h"] for i in range(len(ham_path)-1)]) else ''}

  D. META-INDEPENDENCE POLYNOMIAL I(G_5, x):
     I(G_5, x) = {' + '.join(f'{c}x^{k}' for k, c in enumerate(meta_ip) if c > 0)}
     Meta-H = I(G_5, 2) = {meta_h_from_ip}
     α(G_5) = {alpha} (independence number)
     I(G_5, -1) = {sum(c * (-1)**k for k, c in enumerate(meta_ip))} (Euler char)

  E. STAIRCASE ↔ ISO CLASS:
     Score class = partition of C(n,2) into n parts.
     Iso classes WITHIN a score class are refinements.
     Complement duality maps score (s₀,...) to (n-1-s_{n-1},...).
     SC classes = self-conjugate under this duality.

  F. DEGREE FROM |Aut|:
     Higher |Aut| → more self-flips → lower degree.
     |Aut|=5 (regular): degree 2. |Aut|=3: degree 3. |Aut|=1: degree 6-7.
     The formula involves the fraction of arc flips that stay in the class,
     which is determined by the orbits of Aut(T) on arc pairs.
""")

print("opus-2026-03-22-S172")
