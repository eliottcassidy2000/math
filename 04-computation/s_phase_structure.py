#!/usr/bin/env python3
"""
s_phase_structure.py — Characterize S-phase (beta_3>0) tournaments

At n=6 (exhaustive): S-phase has EXACTLY 320 tournaments with:
  - |Pf| in {7, 9}
  - c3 in {2, 8}
  - Eigenvalue triples: (1, 2.646, 2.646) or (sqrt(3), sqrt(3), 3)

At n=7: S-phase has ~8.4% of random tournaments.

Questions:
1. What score sequences appear in S-phase at n=6?
2. What is the cycle structure (c3, c5) for S-phase?
3. Is there an exact combinatorial criterion?
4. What is beta_3 testing — what "3-dimensional hole" exists?

At n=6, S-phase has Omega with 2 or 8 three-cycles.
With c3=2: two 3-cycles that are vertex-disjoint (use all 6 vertices).
With c3=8: eight 3-cycles — maximum possible.

Let's characterize completely.

Author: kind-pasteur-2026-03-08-S40
"""
import sys, os, time
import numpy as np
from itertools import combinations
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w')
from path_homology_v2 import path_betti_numbers
sys.stdout = _saved

def count_3cycles(A, n):
    c3 = 0
    triples = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]:
                    c3 += 1
                    triples.append((i,j,k))
                if A[i][k] and A[k][j] and A[j][i]:
                    c3 += 1
                    triples.append((i,k,j))
    return c3, triples

def count_5cycles(A, n):
    c5 = 0
    for combo in combinations(range(n), 5):
        # Check all directed 5-cycles on these 5 vertices
        from itertools import permutations
        for p in permutations(combo):
            is_cycle = True
            for idx in range(5):
                if A[p[idx]][p[(idx+1)%5]] != 1:
                    is_cycle = False
                    break
            if is_cycle:
                c5 += 1
    # Each 5-cycle counted 5 times (cyclic rotations)
    return c5 // 5

def H_tournament(A, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[mask][v] == 0: continue
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full])

def score_seq(A, n):
    return tuple(sorted([sum(A[i]) for i in range(n)]))

def is_SC(A, n):
    """Check if tournament is self-complementary (iso to complement)."""
    # A tournament is self-comp if there exists sigma s.t. A[sigma(i)][sigma(j)] = 1 - A[i][j]
    # For n=6, check exhaustively over S_6
    from itertools import permutations
    for perm in permutations(range(n)):
        ok = True
        for i in range(n):
            for j in range(n):
                if i == j: continue
                if A[perm[i]][perm[j]] != 1 - A[i][j]:
                    ok = False
                    break
            if not ok: break
        if ok:
            return True
    return False

# ===== n=6 EXHAUSTIVE: characterize S-phase =====
print("=" * 70)
print("n=6 EXHAUSTIVE: S-phase characterization")
print("=" * 70)

n = 6
m = n * (n-1) // 2
total = 1 << m
t0 = time.time()

s_phase_tours = []
c_phase_tours = []

for bits in range(total):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    try:
        beta = path_betti_numbers(A, n, max_dim=4)
    except:
        continue

    b1 = int(beta[1]) if len(beta) > 1 else 0
    b3 = int(beta[3]) if len(beta) > 3 else 0

    if b3 > 0:
        c3, triples = count_3cycles(A, n)
        c5 = count_5cycles(A, n)
        H = H_tournament(A, n)
        sc = score_seq(A, n)
        s_phase_tours.append({
            'bits': bits, 'A': A, 'c3': c3, 'c5': c5,
            'H': H, 'score': sc, 'triples': triples,
            'beta': beta
        })
    elif b1 > 0:
        c3, triples = count_3cycles(A, n)
        H = H_tournament(A, n)
        sc = score_seq(A, n)
        c_phase_tours.append({
            'bits': bits, 'c3': c3, 'H': H, 'score': sc
        })

    if (bits + 1) % 10000 == 0:
        print(f"  {bits+1}/{total} ({time.time()-t0:.1f}s)")

elapsed = time.time() - t0
print(f"\nDone in {elapsed:.1f}s")
print(f"S-phase: {len(s_phase_tours)} tournaments")
print(f"C-phase: {len(c_phase_tours)} tournaments")

# Analyze S-phase
print("\n" + "=" * 70)
print("S-PHASE DETAILED ANALYSIS")
print("=" * 70)

# c3 distribution
c3_counts = Counter(d['c3'] for d in s_phase_tours)
print(f"\nc3 distribution: {dict(sorted(c3_counts.items()))}")

# c5 distribution
c5_counts = Counter(d['c5'] for d in s_phase_tours)
print(f"c5 distribution: {dict(sorted(c5_counts.items()))}")

# H distribution
H_counts = Counter(d['H'] for d in s_phase_tours)
print(f"H distribution: {dict(sorted(H_counts.items()))}")

# Score sequence
score_counts = Counter(d['score'] for d in s_phase_tours)
print(f"\nScore sequences: {dict(sorted(score_counts.items()))}")

# 3-cycle structure
print("\n--- Cycle structure details ---")

# For c3=2: the two 3-cycles must be vertex-disjoint (6 vertices)
c3_2_tours = [d for d in s_phase_tours if d['c3'] == 2]
print(f"\nc3=2: {len(c3_2_tours)} tournaments")
if c3_2_tours:
    for d in c3_2_tours[:5]:
        print(f"  H={d['H']}, score={d['score']}, c5={d['c5']}")
        # Check if 3-cycles are disjoint
        verts_used = set()
        for t in d['triples']:
            verts_used.update(t)
        print(f"    vertices used by 3-cycles: {verts_used}")
        print(f"    vertex-disjoint? {len(verts_used) == 6}")

# For c3=8: check cycle overlap structure
c3_8_tours = [d for d in s_phase_tours if d['c3'] == 8]
print(f"\nc3=8: {len(c3_8_tours)} tournaments")
if c3_8_tours:
    for d in c3_8_tours[:3]:
        print(f"  H={d['H']}, score={d['score']}, c5={d['c5']}")
        # How many vertex-disjoint pairs?
        triples_as_sets = [frozenset(t) for t in d['triples']]
        disjoint_pairs = 0
        for i in range(len(triples_as_sets)):
            for j in range(i+1, len(triples_as_sets)):
                if not (triples_as_sets[i] & triples_as_sets[j]):
                    disjoint_pairs += 1
        print(f"    disjoint 3-cycle pairs: {disjoint_pairs}")

# C-phase analysis
print("\n" + "=" * 70)
print("C-PHASE ANALYSIS (for comparison)")
print("=" * 70)

c3_counts_C = Counter(d['c3'] for d in c_phase_tours)
print(f"\nc3 distribution: {dict(sorted(c3_counts_C.items()))}")

score_counts_C = Counter(d['score'] for d in c_phase_tours)
print(f"Score sequences: {dict(sorted(score_counts_C.items()))}")

H_counts_C = Counter(d['H'] for d in c_phase_tours)
print(f"H distribution: {dict(sorted(H_counts_C.items()))}")

# Compare: which score sequences are EXCLUSIVE to each phase?
s_scores = set(d['score'] for d in s_phase_tours)
c_scores = set(d['score'] for d in c_phase_tours)
print(f"\nS-only scores: {sorted(s_scores - c_scores)}")
print(f"C-only scores: {sorted(c_scores - s_scores)}")
print(f"Common scores: {sorted(s_scores & c_scores)}")

# Check: does the NUMBER of vertex-disjoint 3-cycle pairs determine the phase?
print("\n" + "=" * 70)
print("INDEPENDENCE NUMBER OF Omega_3")
print("=" * 70)

# For S-phase tournaments, compute alpha(Omega_3) = max independent set in conflict graph
for d in s_phase_tours[:10]:
    triples_as_sets = [frozenset(t) for t in d['triples']]
    # Max independent set (no shared vertex)
    best = 0
    for size in range(len(triples_as_sets), 0, -1):
        found = False
        for combo in combinations(range(len(triples_as_sets)), size):
            all_disjoint = True
            for i in range(len(combo)):
                for j in range(i+1, len(combo)):
                    if triples_as_sets[combo[i]] & triples_as_sets[combo[j]]:
                        all_disjoint = False
                        break
                if not all_disjoint:
                    break
            if all_disjoint:
                best = size
                found = True
                break
        if found:
            break
    print(f"  c3={d['c3']}, H={d['H']}, alpha(Omega_3)={best}, score={d['score']}")
