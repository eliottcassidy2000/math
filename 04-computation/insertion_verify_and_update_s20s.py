#!/usr/bin/env python3
"""
insertion_verify_and_update_s20s.py — kind-pasteur-2026-03-22-S20s

Verify the vertex insertion formula at n=4->5 and n=5->6.
Then derive and implement the E-UPDATE formula.

The staircase delta_k is an isosceles right triangle.
Adding one hypotenuse row: k -> k+1 (one vertex, n tiles)
Adding two leg rows: k -> k+2 (two vertices, source+sink)

The isosceles right triangle property:
  hypotenuse = sqrt(2) * leg
  Adding along hypotenuse: adds n-1 new arcs (1 strip)
  Adding along both legs: adds 2n-1 new arcs (2 strips = L-border)
  Ratio: (2n-1)/(n-1) -> 2 as n -> inf
  This is the PSEUDO-DOUBLING: adding two strips approximately doubles
  the interface, just as the hypotenuse is sqrt(2) times the leg.

Author: kind-pasteur-2026-03-22-S20s
"""

import sys
import numpy as np
from collections import defaultdict
import time

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

def compute_M_and_E(A, n):
    """Compute endpoint matrix M and consecutive pair matrix E by HP enumeration."""
    paths = []
    def find_paths(path, visited):
        if len(path) == n:
            paths.append(tuple(path))
            return
        last = path[-1]
        for w in range(n):
            if w not in visited and A[last][w]:
                find_paths(path + [w], visited | {w})
    for start in range(n):
        find_paths([start], {start})

    M = np.zeros((n, n), dtype=np.int64)
    E = np.zeros((n, n), dtype=np.int64)
    for path in paths:
        M[path[0]][path[-1]] += 1
        for j in range(len(path) - 1):
            E[path[j]][path[j+1]] += 1
    return M, E, paths

def insertion_H(M, E, n_old, sig):
    """Compute H(T+v) using the insertion formula.
    sig[i] = 1 if v beats i, 0 if i beats v."""
    start_count = 0
    end_count = 0
    middle_count = 0
    for b in range(n_old):
        if sig[b] == 1:  # v beats b: v can start a path going to b
            start_count += M[:, b].sum()
    for a in range(n_old):
        if sig[a] == 0:  # a beats v: path ending at a can continue to v
            end_count += M[a, :].sum()
    for a in range(n_old):
        for b in range(n_old):
            if sig[a] == 0 and sig[b] == 1:  # a->v->b insertion
                middle_count += E[a, b]
    return int(start_count + end_count + middle_count)

print("=" * 72)
print("  VERTEX INSERTION: VERIFICATION AND E-UPDATE")
print("=" * 72)

# ========================================================================
# PART 1: VERIFY AT n=4->5 (ALL 16 signatures)
# ========================================================================
print("\n  PART 1: VERIFY AT n=4 -> n=5")

# Use a specific n=4 tournament
A4 = np.array([[0,1,1,1],[0,0,1,0],[0,0,0,1],[0,1,0,0]], dtype=np.int8)
n4 = 4
H4 = count_hp(A4, n4)
M4, E4, _ = compute_M_and_E(A4, n4)

print(f"  Base tournament (n=4): H={H4}")

errors = 0
for sig_bits in range(2**n4):
    sig = [(sig_bits >> i) & 1 for i in range(n4)]

    # Build extended tournament
    A5 = np.zeros((5, 5), dtype=np.int8)
    A5[:4, :4] = A4
    for i in range(4):
        if sig[i]:
            A5[4][i] = 1
        else:
            A5[i][4] = 1

    H5_actual = count_hp(A5, 5)
    H5_pred = insertion_H(M4, E4, 4, sig)

    if H5_actual != H5_pred:
        errors += 1
        print(f"    MISMATCH: sig={sig}, actual={H5_actual}, pred={H5_pred}")

print(f"  Verified: {2**n4 - errors}/{2**n4} correct")
print()

# ========================================================================
# PART 2: VERIFY AT n=5->6 (ALL 32 signatures, multiple base tournaments)
# ========================================================================
print("  PART 2: VERIFY AT n=5 -> n=6")

# Test with several n=5 tournaments
import random
random.seed(42)

total_tests = 0
total_errors = 0

for trial in range(10):
    # Random tournament on 5 vertices
    A5 = np.zeros((5, 5), dtype=np.int8)
    for i in range(5):
        for j in range(i+1, 5):
            if random.random() < 0.5:
                A5[i][j] = 1
            else:
                A5[j][i] = 1

    M5, E5, _ = compute_M_and_E(A5, 5)
    H5 = int(M5.sum())

    for sig_bits in range(2**5):
        sig = [(sig_bits >> i) & 1 for i in range(5)]

        A6 = np.zeros((6, 6), dtype=np.int8)
        A6[:5, :5] = A5
        for i in range(5):
            if sig[i]: A6[5][i] = 1
            else: A6[i][5] = 1

        H6_actual = count_hp(A6, 6)
        H6_pred = insertion_H(M5, E5, 5, sig)

        total_tests += 1
        if H6_actual != H6_pred:
            total_errors += 1

print(f"  Verified: {total_tests - total_errors}/{total_tests} correct (10 base tournaments x 32 sigs)")
print()

# ========================================================================
# PART 3: THE E-UPDATE FORMULA
# ========================================================================
print("  PART 3: DERIVING THE E-UPDATE")
print()

print("""
  When vertex v is added to T (n vertices) to get T' (n+1 vertices):

  HP of T' come from three sources:

  SOURCE A: v is the START. HP = v -> b -> (HP of T from b to end).
    For each b with v->b: each HP of T starting at b gives an HP of T'.
    New consecutive pairs: (v, b) for each such HP.
    Old consecutive pairs: all pairs in the suffix HP (unchanged).

  SOURCE B: v is the END. HP = (HP of T from start to a) -> a -> v.
    For each a with a->v: each HP of T ending at a gives an HP of T'.
    New consecutive pairs: (a, v) for each such HP.
    Old consecutive pairs: all pairs in the prefix HP (unchanged).

  SOURCE C: v is in the MIDDLE. HP = (prefix -> a -> v -> b -> suffix).
    For each pair (a,b) with a->v and v->b:
    each HP of T containing consecutive pair (a,b) gives an HP of T'.
    The pair (a,b) is REPLACED by two pairs (a,v) and (v,b).
    All other consecutive pairs are unchanged.

  E-UPDATE FORMULA:

  E'[x, y] for x, y in {0,...,n} (the new vertex v has index n):

  Case 1: x != v AND y != v (old pair):
    E'[x, y] = E[x, y]                          (from sources A, B: unchanged pairs)
              - sum_{a->v, v->y: x precedes a in HP} correction
              Wait: this is wrong. Source C REMOVES the pair (a,b) and adds (a,v), (v,b).
              So the old pair (a,b) is counted E[a,b] times in source C.
              But only the fraction where x precedes a... no.

  Actually: in Source C, for each HP containing (a,b), the pair (a,b) is SPLIT.
  So (a,b) appears in E'... let me think again.

  In Source A (v starts): the HP of T' is v, b, p_2, ..., p_n.
    Consecutive pairs: (v, b), (b, p_2), (p_2, p_3), ..., (p_{n-1}, p_n).
    The pairs (b, p_2), ..., (p_{n-1}, p_n) are the same as in the original HP of T
    starting at b. So these pairs contribute to E' exactly as they do in E, but
    only for HP starting at b.

  In Source B (v ends): the HP of T' is p_1, p_2, ..., p_{n-1}, a, v.
    Pairs: (p_1, p_2), ..., (p_{n-1}, a), (a, v).
    The pairs (p_1, p_2), ..., (p_{n-1}, a) are from the original HP ending at a.

  In Source C (v in middle): the HP of T' is ..., c, a, v, b, d, ...
    where the original HP had ..., c, a, b, d, ...
    Pairs: ..., (c, a), (a, v), (v, b), (b, d), ...
    The pair (a, b) is REMOVED, replaced by (a, v) and (v, b).
    All other pairs (including (c, a) and (b, d)) are UNCHANGED.
""")

# Let me just COMPUTE E' for a specific case and verify.
print("  COMPUTING E' FOR n=3 -> n=4:")
print()

A3 = np.array([[0,1,0],[0,0,1],[1,0,0]], dtype=np.int8)
M3, E3, paths3 = compute_M_and_E(A3, 3)

print(f"  Base tournament (cycle-3): H={int(M3.sum())}")
print(f"  M3 = {M3.tolist()}")
print(f"  E3 = {E3.tolist()}")
print(f"  Paths: {paths3}")
print()

sig = [1, 0, 1]  # v=3 beats 0 and 2, loses to 1
print(f"  Adding v=3 with signature {sig}")

A4_ext = np.zeros((4, 4), dtype=np.int8)
A4_ext[:3, :3] = A3
for i in range(3):
    if sig[i]: A4_ext[3][i] = 1
    else: A4_ext[i][3] = 1

M4_ext, E4_ext, paths4 = compute_M_and_E(A4_ext, 4)

H4_ext = int(M4_ext.sum())
H4_pred = insertion_H(M3, E3, 3, sig)
print(f"  H(T+v) = {H4_ext}, predicted = {H4_pred}, match = {H4_ext == H4_pred}")
print()
print(f"  New paths: {paths4}")
print()
print(f"  E4 (actual):")
print(f"  {E4_ext.tolist()}")
print()

# Now derive E4 from E3 + sig:
# Source A paths (v=3 starts): v beats {0, 2}. For each b in {0, 2}:
#   HP of T starting at b: paths3 starting at b.
paths_A = []
for b in range(3):
    if sig[b] == 1:
        for p in paths3:
            if p[0] == b:
                paths_A.append((3,) + p)

# Source B paths (v=3 ends): v loses to {1}. For each a in {1}:
#   HP of T ending at a: paths3 ending at a.
paths_B = []
for a in range(3):
    if sig[a] == 0:
        for p in paths3:
            if p[-1] == a:
                paths_B.append(p + (3,))

# Source C paths (v=3 in middle): for each (a,b) with a->v (sig[a]=0) and v->b (sig[b]=1):
#   a=1 (sig=0), b in {0, 2} (sig=1).
#   For each HP of T containing consecutive pair (a,b):
#   Insert v between a and b.
paths_C = []
for p in paths3:
    for j in range(len(p) - 1):
        a, b = p[j], p[j+1]
        if sig[a] == 0 and sig[b] == 1:
            new_path = p[:j+1] + (3,) + p[j+1:]
            paths_C.append(new_path)

all_new_paths = sorted(set(paths_A + paths_B + paths_C))
print(f"  Source A paths: {paths_A}")
print(f"  Source B paths: {paths_B}")
print(f"  Source C paths: {paths_C}")
print(f"  Total new: {len(all_new_paths)}, expected: {H4_ext}")
print(f"  Match: {len(all_new_paths) == H4_ext}")
print()

# Compute E4_derived from the path decomposition
E4_derived = np.zeros((4, 4), dtype=np.int64)
for p in all_new_paths:
    for j in range(len(p) - 1):
        E4_derived[p[j]][p[j+1]] += 1

print(f"  E4 (derived from path decomposition):")
print(f"  {E4_derived.tolist()}")
print(f"  Match actual: {np.array_equal(E4_derived, E4_ext)}")
print()

# Now: express E4_derived in terms of E3, M3, and sig.
# E4[x, y] where neither is v=3:
#   From Source A: pairs not involving v. For HP v,b,...: pairs (b,p_2),...
#     = all pairs in HP of T starting at b, for each b with v->b.
#     = sum_{b: sig[b]=1} E3_starting_at_b[x,y]
#   From Source B: pairs not involving v. For HP ...,a,v: pairs ...(p_{n-1},a).
#     = sum_{a: sig[a]=0} E3_ending_at_a[x,y]
#   From Source C: pairs not involving insertion point.
#     For each (a,b) split by v: all pairs except (a,b) are kept.
#     = sum_{a:sig[a]=0, b:sig[b]=1} (E3_through_ab[x,y] - delta[x=a,y=b])
#     Wait: E3_through_ab means "E restricted to HP containing (a,b)".
#     This needs more detailed tracking than just E3.

# This decomposition requires SUB-MATRICES of E: E restricted to HP
# starting at specific vertices, ending at specific vertices, or
# containing specific pairs. This is MORE than just E.

# HOWEVER: we CAN express E4 in terms of simpler quantities:
# For each pair (x,y) in {0,...,3}^2:

print("  E-UPDATE VERIFICATION (term by term):")
print()

# Track contributions from each source
for x in range(4):
    for y in range(4):
        if E4_ext[x][y] == 0: continue

        # Count from each source
        count_A = sum(1 for p in paths_A if any(p[j]==x and p[j+1]==y for j in range(len(p)-1)))
        count_B = sum(1 for p in paths_B if any(p[j]==x and p[j+1]==y for j in range(len(p)-1)))
        count_C = sum(1 for p in paths_C if any(p[j]==x and p[j+1]==y for j in range(len(p)-1)))

        total = count_A + count_B + count_C
        if total != E4_ext[x][y]:
            print(f"  E4[{x},{y}] = {E4_ext[x][y]}: A={count_A}, B={count_B}, C={count_C}, sum={total} MISMATCH")
        else:
            print(f"  E4[{x},{y}] = {E4_ext[x][y]}: A={count_A}, B={count_B}, C={count_C}")
