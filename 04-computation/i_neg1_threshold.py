#!/usr/bin/env python3
"""
i_neg1_threshold.py — opus-2026-03-14-S71e

INVESTIGATION: At what n does I(Omega, -1) become universally non-positive?

Data so far:
  n=3: I(-1)>0 for 75% (transitive: I(-1)=1)
  n=4: I(-1)>0 for 37.5%
  n=5: I(-1)>0 for 11.7%
  n=6: I(-1)>0 for 2.2%
  n=7: I(-1)>0 for 0.25%
  n=8: I(-1) in [-158, -6] for ALL (exhaustive) — ALWAYS NEGATIVE
  n=9: I(-1)=1 for transitive, I(-1)<0 for all 100 random

The transitive tournament ALWAYS has I(-1)=1 (since α_k=0 for all k).
So I(-1) > 0 is possible at ALL n.

But the QUESTION is: aside from transitive (and near-transitive),
does I(-1) ever stay positive?

NEW QUESTION: Is I(-1) ∈ {1, 0, negative} for ALL tournaments?
Does I(-1)=1 characterize the transitive tournament?
Does I(-1)≥0 have a nice characterization?

The near-transitive flip pattern: 1, 0, -1, -3, -7
  Diffs: -1, -1, -2, -4 = 2^0, 2^0, 2^1, 2^2
  Hmm, almost. Let me look more carefully at what each flip does.
"""

import sys
import numpy as np
from itertools import combinations
from collections import Counter
sys.stdout.reconfigure(line_buffering=True)

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms: continue
            for v in range(n):
                if not (mask & (1 << v)): continue
                pm = mask ^ (1 << v)
                t = 0
                for u in range(n):
                    if (pm & (1 << u)) and A[u][v]:
                        t += dp.get((pm, u), 0)
                if t:
                    dp[(mask, v)] = t
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

def count_ham_cycles(A, n):
    if n < 3: return 0
    full_mask = (1 << n) - 1
    dp = {(1 << 0, 0): 1}
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms: continue
            if not (mask & 1): continue
            for v in range(n):
                if not (mask & (1 << v)): continue
                if v == 0 and ms < n: continue
                pm = mask ^ (1 << v)
                if not (pm & 1): continue
                t = 0
                for u in range(n):
                    if (pm & (1 << u)) and A[u][v]:
                        t += dp.get((pm, u), 0)
                if t:
                    dp[(mask, v)] = t
    total = 0
    for v in range(1, n):
        if A[v][0] and (full_mask, v) in dp:
            total += dp[(full_mask, v)]
    return total

def count_directed_k_cycles(A, n, k):
    if k > n: return 0
    total = 0
    for combo in combinations(range(n), k):
        verts = list(combo)
        sub = np.zeros((k, k), dtype=int)
        for i in range(k):
            for j in range(k):
                sub[i][j] = A[verts[i]][verts[j]]
        total += count_ham_cycles(sub, k)
    return total

def compute_alpha(A, n):
    """Compute α₁ (total directed odd cycles)"""
    total = 0
    for k in range(3, n+1, 2):
        total += count_directed_k_cycles(A, n, k)
    return total

# ======================================================================
# PART 1: I(-1) for transitive and near-transitive at various n
# ======================================================================
print("=" * 70)
print("PART 1: TRANSITIVE AND NEAR-TRANSITIVE I(-1)")
print("=" * 70)

for n in [3, 4, 5, 6, 7]:
    tb = n*(n-1)//2
    print(f"\n  n={n}: (transitive = bits 0)")

    # For small n, exhaustively check which tournaments have I(-1) >= 0
    if n <= 6:
        pos_count = 0
        zero_count = 0
        pos_examples = []
        for bits in range(1 << tb):
            A = bits_to_adj(bits, n)
            a1 = compute_alpha(A, n)
            # For n<=7, α₂ needs disjoint pair counting
            # But for I(-1) at n<=5 (α₂=0): I(-1) = 1 - α₁
            if n <= 5:
                a2 = 0
            else:
                # n=6: count disjoint (3,3) pairs
                H = count_ham_paths(A, n)
                a2 = (H - 1 - 2*a1) // 4

            I_neg1 = 1 - a1 + a2
            if I_neg1 > 0:
                pos_count += 1
                if len(pos_examples) < 5:
                    scores = tuple(sorted([sum(A[i]) for i in range(n)]))
                    pos_examples.append((bits, I_neg1, a1, a2, scores))
            elif I_neg1 == 0:
                zero_count += 1

        total = 1 << tb
        print(f"    I(-1) > 0: {pos_count}/{total} = {pos_count/total*100:.2f}%")
        print(f"    I(-1) = 0: {zero_count}/{total} = {zero_count/total*100:.2f}%")
        print(f"    I(-1) < 0: {total-pos_count-zero_count}/{total}")
        if pos_examples:
            print(f"    Examples with I(-1) > 0:")
            for bits, val, a1, a2, scores in pos_examples[:5]:
                print(f"      bits={bits}: I(-1)={val}, α₁={a1}, α₂={a2}, scores={scores}")

# ======================================================================
# PART 2: I(-1) = 1 characterization
# ======================================================================
print("\n" + "=" * 70)
print("PART 2: DOES I(-1)=1 CHARACTERIZE TRANSITIVE?")
print("=" * 70)

for n in [3, 4, 5, 6]:
    tb = n*(n-1)//2
    ones = []
    for bits in range(1 << tb):
        A = bits_to_adj(bits, n)
        a1 = compute_alpha(A, n)
        if n <= 5:
            a2 = 0
        else:
            H = count_ham_paths(A, n)
            a2 = (H - 1 - 2*a1) // 4
        I_neg1 = 1 - a1 + a2
        if I_neg1 == 1:
            scores = tuple(sorted([sum(A[i]) for i in range(n)]))
            ones.append((bits, scores))

    print(f"\n  n={n}: {len(ones)} tournaments with I(-1)=1")
    for bits, scores in ones[:10]:
        print(f"    bits={bits}, scores={scores}")

    # Check if all are isomorphic to transitive
    is_trans = all(s == tuple(range(n)) for _, s in ones)
    print(f"    All transitive? Score check: {is_trans}")

# ======================================================================
# PART 3: Near-transitive edge-flip analysis
# ======================================================================
print("\n" + "=" * 70)
print("PART 3: EDGE-FLIP ANALYSIS FROM TRANSITIVE")
print("=" * 70)

for n in [5, 6, 7]:
    tb = n*(n-1)//2
    print(f"\n  n={n}: flipping edges from transitive (bits=0)")

    # In the transitive tournament with bits=0:
    # Edge (i,j) with i<j: bit index = sum_{r=0}^{i-1}(n-1-r) + (j-i-1)
    # With bits=0: A[j][i]=1 for all i<j (j beats i)
    # Flipping bit k reverses one edge.

    # Map bit index to edge
    idx = 0
    edges = []
    for i in range(n):
        for j in range(i+1, n):
            edges.append((i, j))
            idx += 1

    print(f"    Edge mapping: {edges[:8]}...")

    # Flip each single edge and compute I(-1)
    vals = []
    for flip_idx in range(min(tb, 15)):
        bits = 1 << flip_idx
        A = bits_to_adj(bits, n)
        a1 = compute_alpha(A, n)
        if n <= 5:
            a2 = 0
        elif n <= 7:
            H = count_ham_paths(A, n)
            a2 = (H - 1 - 2*a1) // 4
        I_neg1 = 1 - a1 + a2
        edge = edges[flip_idx]
        vals.append(I_neg1)
        print(f"    flip edge {edge}: α₁={a1:3d}, α₂={a2:3d}, I(-1)={I_neg1}")

# ======================================================================
# PART 4: The 2^k pattern in edge flips
# ======================================================================
print("\n" + "=" * 70)
print("PART 4: HOW MANY NEW CYCLES DOES FLIPPING EDGE (i,j) CREATE?")
print("=" * 70)

# In the transitive tournament 0→1→2→...→(n-1),
# flipping edge (i,j) (so now i→j instead of j→i)
# creates new cycles through vertices i and j.
# A k-cycle through (i,j) needs a path from j to i
# through k-2 intermediate vertices, all from {i+1,...,j-1}.
# (Because in the transitive tournament, all edges point "upward".)
#
# Wait, that's not quite right. Let me think again.
# Transitive: edge u→v iff u < v (in our bits=0 encoding, A[j][i]=1 means j→i...
# Actually let me check.

n = 5
A_trans = bits_to_adj(0, n)
print(f"\n  Transitive n=5 adjacency (bits=0):")
for i in range(n):
    print(f"    row {i}: {list(A_trans[i])}")

# In our encoding with bits=0: for i<j, bit=0 means A[j][i]=1.
# So j→i for all i<j. That means higher-numbered vertices beat lower.
# Score sequence: vertex k has out-degree n-1-k.
# Vertex 0: out-degree n-1, vertex n-1: out-degree 0.
# Wait, let me recheck...
scores_trans = [sum(A_trans[i]) for i in range(n)]
print(f"  Out-degrees: {scores_trans}")
print(f"  So vertex 0 beats nobody, vertex n-1 beats everybody")

# Flipping edge between vertices i and j (i<j):
# Currently j→i (j beats i). After flip: i→j (i beats j).
# Now vertex i can reach j directly.
# To form a directed cycle through i→j, need path from j back to i.
# In the transitive part, all edges go from higher to lower index.
# So j can reach j-1, j-2, ..., 0 in the transitive part.
# i can be reached from i+1, i+2, ..., n-1 in the transitive part.
# A directed cycle must use the flipped edge i→j and then
# return from j to i using only the transitive edges.
# But transitive edges go from higher to lower, so j→(j-1)→...→i works!
# That's a cycle of length j-i+1 using vertices i, i+1, ..., j.
# But wait: j→j-1 is the edge, j-1→j-2 is the edge, ..., i+1→i is the edge.
# And we added i→j. So the cycle is i→j→j-1→...→i+1→i, length j-i+1.
# But wait, this needs i+1→i which is edge from higher to lower = transitive.
# No! Our transitive has HIGHER beats LOWER. So the edge is (i+1)→i?
# Let's check: vertex i+1 has A[i+1][i] = ?

print(f"\n  A[1][0] = {A_trans[1][0]} (does vertex 1 beat vertex 0?)")
print(f"  A[4][3] = {A_trans[4][3]} (does vertex 4 beat vertex 3?)")
print(f"  A[0][1] = {A_trans[0][1]} (does vertex 0 beat vertex 1?)")

# With bits=0: for pair (i,j) with i<j, the bit index maps to A[i][j]=0, A[j][i]=1.
# So A[j][i]=1 means j beats i. Higher index beats lower.
# Edge j→i means A[j][i]=1. ✓

# So transitive has edges j→i for all j>i.
# Score of vertex v: out-degree = number of vertices it beats = v (beats 0,1,...,v-1)
# Wait: A[v][u] = 1 for u < v (v beats u). So out-degree of v = v.
# vertex 0: out-degree 0, vertex n-1: out-degree n-1. Correct!

# Flipping bit for pair (i,j): now A[i][j]=1, A[j][i]=0. So i beats j.
# Previously j beat i. Now i→j is the new edge.
# Can we form a cycle? Need path j → ... → i using transitive edges.
# Transitive: u→v iff u > v. So j can reach anything < j, including i.
# Path: j → j-1 → j-2 → ... → i. Length j-i edges, cycle length j-i+1.
# But this is an (j-i+1)-cycle. It exists iff j-i+1 is odd (odd cycle needed for α₁).
# j-i+1 odd ↔ j-i even ↔ i and j have same parity.

print(f"\n  KEY INSIGHT:")
print(f"  Flipping edge (i,j) in transitive creates a cycle of length j-i+1.")
print(f"  For α₁ count: need j-i+1 odd, i.e., j-i even (same parity).")
print(f"")
print(f"  But it might create MORE than one cycle!")
print(f"  Any subset of {{i, i+1, ..., j}} of odd size ≥ 3 that includes i,j")
print(f"  and forms a Hamiltonian cycle in the subtournament on those vertices")
print(f"  is a new directed odd cycle.")

# Actually, the subtournament on {i, i+1, ..., j} has the edge i→j plus
# all transitive edges u→v for u>v within that set.
# This is a "near-transitive" subtournament with exactly 1 reversed edge.
# How many directed cycles does it have?

for n in [7]:
    print(f"\n  n={n}: edge-flip cycle creation table")
    print(f"  {'edge':>8s} {'gap':>4s} {'dc3':>5s} {'dc5':>5s} {'dc7':>5s} {'alpha1':>7s} {'alpha2':>7s} {'I(-1)':>6s}")

    idx = 0
    edges = []
    for i in range(n):
        for j in range(i+1, n):
            edges.append((i, j))

    for flip_idx in range(len(edges)):
        bits = 1 << flip_idx
        A = bits_to_adj(bits, n)
        dc3 = count_directed_k_cycles(A, n, 3)
        dc5 = count_directed_k_cycles(A, n, 5)
        dc7 = count_ham_cycles(A, n)
        a1 = dc3 + dc5 + dc7
        H = count_ham_paths(A, n)
        a2 = (H - 1 - 2*a1) // 4
        I_neg1 = 1 - a1 + a2
        i, j = edges[flip_idx]
        gap = j - i
        print(f"  ({i},{j}){' ':>4s} {gap:4d} {dc3:5d} {dc5:5d} {dc7:5d} {a1:7d} {a2:7d} {I_neg1:6d}")

# ======================================================================
# PART 5: THE GAP PATTERN
# ======================================================================
print("\n" + "=" * 70)
print("PART 5: CYCLE COUNT vs GAP IN TRANSITIVE EDGE-FLIP")
print("=" * 70)

# For n=7, flipping edge (i,j) with gap g = j-i:
# The key cycle is of length g+1.
# gap=1: 2-cycle (not odd), so dc3=0
# gap=2: 3-cycle on {i, i+1, i+2}, dc3=1
# gap=3: 4-cycle (not odd), but might create a 3-cycle? No, need i→j edge.
# gap=4: 5-cycle on {i,...,i+4}, plus possible 3-cycles on subsets
# gap=5: 6-cycle (not odd for main cycle), but 3-cycles and 5-cycles possible

print("""
  For the transitive tournament with ONE edge (i,j) flipped (i<j, gap g=j-i):

  The subtournament on {i, i+1, ..., j} is "transitive + 1 reversal".
  In this subtournament:
    - All edges k→l for k>l (transitive direction)
    - PLUS the edge i→j (the reversed one)

  The directed odd cycles that USE the edge i→j:
    For each odd-size subset S of {i,i+1,...,j} containing both i and j,
    if S forms a directed cycle using i→j and the transitive edges.

  The key cycle is: i → j → j-1 → ... → i+1 → i (length g+1).
  This exists iff g+1 ≥ 3, i.e., g ≥ 2.

  For g even: g+1 is odd, so the main cycle IS an odd cycle.
    Example: g=2, cycle length 3 (triangle)
    Example: g=4, cycle length 5 (pentagon)
  For g odd: g+1 is even, so the main cycle is NOT odd.
    But there might still be odd subcycles.

  Subcycles of length 2k+1 (odd):
    Choose a subset of size 2k+1 from {i,...,j} including i and j.
    The remaining 2k+1-(2) = 2k-1 vertices from {i+1,...,j-1}.
    The subset has one reversed edge i→j, rest transitive.
    For it to form a directed cycle: i→j→(decreasing chain)→i.
    The chain from j to i must be decreasing = use only consecutive
    vertices in the gap. Actually any subset works because transitive
    means all edges point "down" within the subset.

  CLAIM: For gap g, the number of odd cycles through edge (i,j) is
    sum_{k=0}^{floor((g-2)/2)} C(g-1, 2k) = 2^{g-2}  (for g ≥ 2)

  Because: choose 2k intermediate vertices from g-1 available (indices i+1,...,j-1),
  giving a cycle of length 2k+2. We need 2k+2 odd, so 2k+2 = 2k+2... wait.
  Cycle has i, j, and 2k intermediates = 2k+2 total vertices.
  Odd cycle needs 2k+2 odd → k is odd... that doesn't work.

  Let me just count directly.
""")

# Direct count: for each gap g, how many odd-length directed cycles
# pass through the reversed edge?
print("  Direct cycle count by gap g at n=7:")
n = 7
edges = []
for i in range(n):
    for j in range(i+1, n):
        edges.append((i, j))

gap_data = {}
for flip_idx in range(len(edges)):
    bits = 1 << flip_idx
    A = bits_to_adj(bits, n)
    dc3 = count_directed_k_cycles(A, n, 3)
    dc5 = count_directed_k_cycles(A, n, 5)
    dc7 = count_ham_cycles(A, n)
    i, j = edges[flip_idx]
    g = j - i
    if g not in gap_data:
        gap_data[g] = []
    gap_data[g].append({'dc3': dc3, 'dc5': dc5, 'dc7': dc7,
                        'a1': dc3 + dc5 + dc7, 'edge': (i,j)})

for g in sorted(gap_data.keys()):
    entries = gap_data[g]
    # Check if all entries with same gap have same cycle counts
    dc3s = set(d['dc3'] for d in entries)
    dc5s = set(d['dc5'] for d in entries)
    dc7s = set(d['dc7'] for d in entries)
    print(f"    gap {g}: dc3={dc3s}, dc5={dc5s}, dc7={dc7s}, α₁ = {set(d['a1'] for d in entries)}")
    if len(dc3s) == 1 and len(dc5s) == 1:
        print(f"           UNIVERSAL within same gap!")

# ======================================================================
# PART 6: DO CYCLE COUNTS ONLY DEPEND ON GAP?
# ======================================================================
print("\n" + "=" * 70)
print("PART 6: UNIVERSAL GAP FORMULA")
print("=" * 70)

# For a single edge flip in a transitive tournament,
# the number of k-cycles depends ONLY on the gap, not on position.
# This is because the subtournament on {i,...,j} is always isomorphic
# to the "transitive + 1 reversal" on g+1 vertices.

# Let's verify and find the pattern
for n in [5, 7, 9]:
    print(f"\n  n={n}:")
    edges = []
    for i in range(n):
        for j in range(i+1, n):
            edges.append((i, j))

    gap_cycles = {}
    for flip_idx in range(len(edges)):
        bits = 1 << flip_idx
        A = bits_to_adj(bits, n)
        dc3 = count_directed_k_cycles(A, n, 3)
        dc5 = count_directed_k_cycles(A, n, 5) if n >= 5 else 0
        dc7 = count_directed_k_cycles(A, n, 7) if n >= 7 else 0
        dc9 = count_ham_cycles(A, n) if n >= 9 else 0
        i, j = edges[flip_idx]
        g = j - i
        key = g
        if key not in gap_cycles:
            gap_cycles[key] = {'dc3': set(), 'dc5': set(), 'dc7': set(), 'dc9': set()}
        gap_cycles[key]['dc3'].add(dc3)
        gap_cycles[key]['dc5'].add(dc5)
        gap_cycles[key]['dc7'].add(dc7)
        gap_cycles[key]['dc9'].add(dc9)

    all_universal = True
    for g in sorted(gap_cycles.keys()):
        d = gap_cycles[g]
        u3 = len(d['dc3']) == 1
        u5 = len(d['dc5']) == 1
        u7 = len(d['dc7']) == 1
        u9 = len(d['dc9']) == 1
        universal = u3 and u5 and u7 and u9
        if not universal:
            all_universal = False
        dc3_val = list(d['dc3'])[0] if u3 else d['dc3']
        dc5_val = list(d['dc5'])[0] if u5 else d['dc5']
        dc7_val = list(d['dc7'])[0] if u7 else d['dc7']
        dc9_val = list(d['dc9'])[0] if u9 else d['dc9']
        print(f"    gap {g}: dc3={dc3_val}, dc5={dc5_val}, dc7={dc7_val}, dc9={dc9_val}{'  UNIVERSAL' if universal else '  NOT UNIVERSAL'}")

    print(f"    All gaps universal: {all_universal}")

# Now look at the pattern
print("\n  Pattern for dc3 at gap g (from n=9 data):")
print("    gap 1: dc3=0 (no triangle possible)")
print("    gap 2: dc3=1 (one triangle: {i, i+1, j})")
print("    gap 3: dc3=? (need even-length main cycle)")
print("    gap 4: dc3=?")
print("    gap 5: dc3=?")

# Explicit formula: dc3(gap g) = C(g-1, 1) for g even, C(g-1,1) for g odd?
# No, let me check from the data.

from math import comb
print("\n  dc3 vs gap (from n=9):")
for g in sorted(gap_cycles.keys()):
    dc3_val = list(gap_cycles[g]['dc3'])[0]
    formula_val = comb(g-1, 2) if g >= 2 else 0
    print(f"    gap {g}: dc3={dc3_val}, C(g-1,2)={formula_val}, match={dc3_val==formula_val}")

print("\n  dc5 vs gap (from n=9):")
for g in sorted(gap_cycles.keys()):
    dc5_val = list(gap_cycles[g]['dc5'])[0]
    formula_val = comb(g-1, 4) if g >= 4 else 0
    print(f"    gap {g}: dc5={dc5_val}, C(g-1,4)={formula_val}, match={dc5_val==formula_val}")

print("\n  dc7 vs gap (from n=9):")
for g in sorted(gap_cycles.keys()):
    dc7_val = list(gap_cycles[g]['dc7'])[0]
    formula_val = comb(g-1, 6) if g >= 6 else 0
    print(f"    gap {g}: dc7={dc7_val}, C(g-1,6)={formula_val}, match={dc7_val==formula_val}")

print("\n  dc9 vs gap (from n=9):")
for g in sorted(gap_cycles.keys()):
    dc9_val = list(gap_cycles[g]['dc9'])[0]
    formula_val = comb(g-1, 8) if g >= 8 else 0
    print(f"    gap {g}: dc9={dc9_val}, C(g-1,8)={formula_val}, match={dc9_val==formula_val}")

print("""
  CONJECTURE: For the transitive tournament with one edge (i,j) flipped (gap g=j-i):
    dc_{2k+1} = C(g-1, 2k) for k ≥ 1.

  Explanation: A (2k+1)-cycle through edge (i,j) consists of:
    - The reversed edge i→j
    - A path j→...→i through 2k-1 intermediate vertices from {i+1,...,j-1}
    - There are g-1 intermediates available, choose 2k of them
    Wait, 2k+1 cycle = 2k+1 vertices = i, j, and 2k-1 intermediates.
    Choose 2k-1 from g-1 available: C(g-1, 2k-1)? Let me recheck...

    Actually: the (2k+1)-cycle uses 2k+1 VERTICES. Two are fixed (i and j).
    So we choose 2k+1 - 2 = 2k - 1 intermediates from the g-1 available.
    So dc_{2k+1} should be C(g-1, 2k-1)?

    At g=2, 2k+1=3: C(1, 1) = 1 ✓
    At g=4, 2k+1=3: C(3, 1) = 3. But data says dc3=? Let me check.
    At g=4, 2k+1=5: C(3, 3) = 1.
""")

# Actually I need to recheck. The cycle length is the number of vertices,
# and we need the cycle to be a DIRECTED cycle.
# In the single-flip transitive tournament, any subset containing i,j
# with the "right parity" forms exactly one directed cycle.
# The key is: in a transitive tournament on m vertices with one edge reversed,
# there is exactly ONE directed cycle containing the reversed edge,
# and it visits ALL m vertices.
# Wait, that's not right either. Let me count more carefully.

print("\nDone.")
