#!/usr/bin/env python3
"""
omega4_finer.py — opus-2026-03-13-S70

Omega_4 on 5-vertex subtournaments is NOT determined by (score_seq, t3).
What DOES determine it? Investigate finer invariants.

The problematic types:
  scores=(1,1,2,3,3), t3=3: omega4_local ∈ {0, 1}
  scores=(1,2,2,2,3), t3=4: omega4_local ∈ {0, 1, 2}
"""

import numpy as np
from math import comb
from itertools import combinations, permutations
from collections import Counter, defaultdict

def adj_matrix(n, T_bits):
    A = np.zeros((n,n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (T_bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def omega4_on_subset(A, vertices):
    count = 0
    for perm in permutations(vertices):
        v0, v1, v2, v3, v4 = perm
        if (A[v0][v1] and A[v1][v2] and A[v2][v3] and A[v3][v4] and
            A[v0][v2] and A[v1][v3] and A[v2][v4]):
            count += 1
    return count

def count_3cycles(A, vertices):
    count = 0
    vlist = list(vertices)
    for i in range(len(vlist)):
        for j in range(i+1, len(vlist)):
            for k in range(j+1, len(vlist)):
                a, b, c = vlist[i], vlist[j], vlist[k]
                if (A[a][b] and A[b][c] and A[c][a]) or \
                   (A[a][c] and A[c][b] and A[b][a]):
                    count += 1
    return count

def score_sequence(A, vertices):
    scores = []
    for v in vertices:
        s = sum(A[v][w] for w in vertices if w != v)
        scores.append(s)
    return tuple(sorted(scores))

def count_transitive_triples(A, vertices):
    """Count ordered transitive triples (a,b,c) with a→b→c, a→c."""
    count = 0
    for a in vertices:
        for b in vertices:
            if b == a or not A[a][b]: continue
            for c in vertices:
                if c == a or c == b: continue
                if A[b][c] and A[a][c]:
                    count += 1
    return count

def count_4cycles_oriented(A, vertices):
    """Count oriented 4-cycles a→b→c→d→a."""
    count = 0
    vlist = list(vertices)
    for a in vlist:
        for b in vlist:
            if b == a or not A[a][b]: continue
            for c in vlist:
                if c == a or c == b or not A[b][c]: continue
                for d in vlist:
                    if d == a or d == b or d == c: continue
                    if A[c][d] and A[d][a]:
                        count += 1
    return count

def tournament_isomorphism_type(A, vertices):
    """Get a canonical form for the tournament on these vertices."""
    vlist = list(vertices)
    n = len(vlist)
    # Try all permutations, return lexicographically smallest adjacency encoding
    best = None
    for perm in permutations(range(n)):
        encoding = []
        for i in range(n):
            for j in range(i+1, n):
                encoding.append(A[vlist[perm[i]]][vlist[perm[j]]])
        encoding = tuple(encoding)
        if best is None or encoding < best:
            best = encoding
    return best

# ============================================================
# n = 6: Deep analysis of problematic types
# ============================================================
n = 6
num_pairs = n*(n-1)//2
total = 2**num_pairs

print("="*70)
print("OMEGA_4 FINER INVARIANTS: n = 6")
print("="*70)

# Collect data for 5-vertex subtournaments
print("  Collecting data...")
data_by_iso = defaultdict(list)

for bits in range(total):
    A = adj_matrix(n, bits)
    for combo in combinations(range(n), 5):
        val = omega4_on_subset(A, list(combo))
        scores = score_sequence(A, combo)
        t3 = count_3cycles(A, combo)
        iso = tournament_isomorphism_type(A, combo)
        data_by_iso[iso].append(val)
    if bits % 5000 == 0 and bits > 0:
        print(f"    {bits}/{total}...")

print(f"\n  Number of 5-vertex tournament isomorphism types: {len(data_by_iso)}")
print(f"\n  omega4_local by isomorphism type:")

for iso, vals in sorted(data_by_iso.items()):
    unique_vals = sorted(set(vals))
    # Get score sequence and t3 for this type
    # Reconstruct a 5-vertex tournament from iso
    A_tmp = np.zeros((5,5), dtype=int)
    idx = 0
    for i in range(5):
        for j in range(i+1, 5):
            if iso[idx]:
                A_tmp[i][j] = 1
            else:
                A_tmp[j][i] = 1
            idx += 1
    scores = score_sequence(A_tmp, range(5))
    t3 = count_3cycles(A_tmp, range(5))
    tt = count_transitive_triples(A_tmp, range(5))
    c4 = count_4cycles_oriented(A_tmp, range(5))

    print(f"    iso={iso}, scores={scores}, t3={t3}, tt={tt}, c4={c4}: "
          f"omega4_local={unique_vals[0]} (all {len(vals)} agree)")

# Now check: is omega4_local determined by isomorphism type?
print(f"\n  Is omega4_local determined by isomorphism type?")
all_determined = True
for iso, vals in data_by_iso.items():
    if len(set(vals)) > 1:
        all_determined = False
        print(f"    NOT determined: iso={iso}, vals={sorted(set(vals))}")
if all_determined:
    print("    YES — omega4_local is a function of tournament isomorphism type")

# Relationship table: for each iso type, show omega4_local and key invariants
print(f"\n  SUMMARY TABLE:")
print(f"  {'scores':20s} {'t3':>4s} {'tt':>4s} {'c4':>4s} {'o4_local':>8s}")
print(f"  {'-'*20} {'----':>4s} {'----':>4s} {'----':>4s} {'--------':>8s}")

results = []
for iso, vals in sorted(data_by_iso.items()):
    A_tmp = np.zeros((5,5), dtype=int)
    idx = 0
    for i in range(5):
        for j in range(i+1, 5):
            if iso[idx]:
                A_tmp[i][j] = 1
            else:
                A_tmp[j][i] = 1
            idx += 1
    scores = score_sequence(A_tmp, range(5))
    t3 = count_3cycles(A_tmp, range(5))
    tt = count_transitive_triples(A_tmp, range(5))
    c4 = count_4cycles_oriented(A_tmp, range(5))
    o4 = vals[0]  # should be constant
    results.append((scores, t3, tt, c4, o4))
    print(f"  {str(scores):20s} {t3:4d} {tt:4d} {c4:4d} {o4:8d}")

# Check: is omega4_local = C(5,4) - B for some correction B?
# C(5,4) = 5, but omega4_local can be 0,1,2,5
# For the regular tournament: omega4_local = 5 = C(5,4)
# So correction B = 5 - omega4_local
print(f"\n  Correction B = 5 - omega4_local:")
for scores, t3, tt, c4, o4 in results:
    B = 5 - o4
    print(f"    scores={str(scores):20s} t3={t3}, B={B}")

# Check: is omega4_local related to the parity of some count?
print(f"\n  Parity analysis:")
for scores, t3, tt, c4, o4 in results:
    print(f"    scores={str(scores):20s} t3≡{t3%2} c4≡{c4%2} o4={o4}")

# Check if omega4_local is determined by (t3, c4) pair
print(f"\n  Is omega4_local = f(t3, c4)?")
by_t3_c4 = defaultdict(set)
for scores, t3, tt, c4, o4 in results:
    by_t3_c4[(t3, c4)].add(o4)
for (t3, c4), o4_set in sorted(by_t3_c4.items()):
    det = "YES" if len(o4_set) == 1 else "NO"
    print(f"    t3={t3}, c4={c4}: omega4_local ∈ {sorted(o4_set)} [{det}]")

# Check: Omega_4 vs neighborhood cycle counts (like Omega_3 formula)
# For Omega_3: correction = sum_v [c3(N_in) + c3(N_out)]
# Maybe for Omega_4: correction involves c4 of neighborhoods or some higher structure
print(f"\n  Neighborhood analysis for Omega_4 at n=6:")
# For each tournament, compute Omega_4 and neighborhood-based invariants
for bits in [0, 100, 500, 1000, 5000]:
    A = adj_matrix(n, bits)
    o4 = 0
    for start in range(n):
        path = [start]
        def count_reg4(path, depth):
            c = 0
            if depth == 4:
                return 1
            for v in range(n):
                if v not in path and A[path[-1]][v]:
                    if depth >= 1 and not A[path[-2]][v]:
                        continue
                    path.append(v)
                    c += count_reg4(path, depth+1)
                    path.pop()
            return c
        o4 += count_reg4(path, 0)

    # Neighborhood counts
    for v in range(n):
        n_in = [w for w in range(n) if w != v and A[w][v]]
        n_out = [w for w in range(n) if w != v and A[v][w]]

    # Omega_3 of each vertex's in/out neighborhood? No, Omega_3 needs 4 vertices
    # Count c3 in neighborhoods of pairs?
    # Let's try: for each edge (u,v), count c3 in common in-neighborhood, etc.

    print(f"    T={bits}: Omega_4={o4}")

print("\nDONE.")
