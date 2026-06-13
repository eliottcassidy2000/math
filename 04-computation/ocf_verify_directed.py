#!/usr/bin/env python3
"""
ocf_verify_directed.py — opus-2026-03-14-S71d

CRITICAL FIX: Ω(T) vertices are DIRECTED odd cycles, not vertex sets.

A 3-vertex set has at most 1 directed 3-cycle (tournament).
A 5-vertex set can have multiple directed 5-cycles!
  - In a tournament on 5 vertices, the max number of directed 5-cycles is 12.
  - (5-1)!/2 = 12 is the max (all directions consistent).

Previous scripts counted vertex-SETS, not directed cycles.
This caused I(CG,2) ≠ H at n=5.

Let's fix this and re-verify.
"""

import sys, time
import numpy as np
from itertools import combinations, permutations
from collections import defaultdict
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

def find_directed_cycles(A, n, max_length=None):
    """Find all DIRECTED odd cycles as normalized tuples (not vertex sets).

    A directed cycle (v0, v1, ..., v_{k-1}) means v0→v1→...→v_{k-1}→v0.
    Normalize: smallest vertex first, direction chosen so second < last.
    This avoids counting each cycle k times (once per starting vertex).
    """
    if max_length is None:
        max_length = n

    cycles = []
    for length in range(3, max_length + 1, 2):  # odd lengths only
        if length > n:
            break
        for combo in combinations(range(n), length):
            verts = list(combo)
            # Fix first vertex as min, check all orderings of rest
            v0 = verts[0]
            for perm in permutations(verts[1:]):
                order = (v0,) + perm
                if all(A[order[i]][order[(i+1) % length]] for i in range(length)):
                    # Normalize: fix start = min vertex, choose direction
                    # so that the neighbor of v0 in the cycle with smaller index
                    # comes first
                    # Actually: each directed cycle has exactly ONE canonical form
                    # with smallest vertex first.
                    # The cycle v0→perm[0]→...→perm[-1]→v0
                    # vs reverse: v0→perm[-1]→...→perm[0]→v0
                    # Keep the one where perm[0] < perm[-1]
                    if perm[0] < perm[-1]:
                        cycles.append(order)
                    # (the reverse direction would have perm[-1] first)

    return cycles

def conflict_graph_and_ip(cycles, x_vals):
    """Build conflict graph and compute I(Ω, x)."""
    nc = len(cycles)
    # Two directed cycles conflict iff they share a vertex
    cycle_vsets = [frozenset(c) for c in cycles]

    adj = [[False]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if cycle_vsets[i] & cycle_vsets[j]:
                adj[i][j] = True
                adj[j][i] = True

    # Count independent sets
    alpha = [0] * (nc + 1)
    for mask in range(1 << nc):
        bits_list = [i for i in range(nc) if mask & (1 << i)]
        indep = True
        for a in range(len(bits_list)):
            for b in range(a+1, len(bits_list)):
                if adj[bits_list[a]][bits_list[b]]:
                    indep = False
                    break
            if not indep:
                break
        if indep:
            alpha[len(bits_list)] += 1

    results = {}
    for x in x_vals:
        val = sum(alpha[k] * x**k for k in range(nc + 1))
        results[x] = val

    return results, alpha

# ======================================================================
# First, let's understand: how many directed cycles per vertex set?
# ======================================================================
print("=" * 70)
print("DIRECTED vs VERTEX-SET CYCLE COUNTING")
print("=" * 70)

n = 5
A_sample = bits_to_adj(40, n)
print(f"\n  Example: bits=40, n=5")
print(f"  Adjacency matrix:")
for i in range(n):
    row = ' '.join(str(A_sample[i][j]) for j in range(n))
    print(f"    [{row}]")

# Count directed 3-cycles
dc3 = []
for combo in combinations(range(n), 3):
    a, b, c = combo
    # Two possible directed 3-cycles: a→b→c→a or a→c→b→a
    if A_sample[a][b] and A_sample[b][c] and A_sample[c][a]:
        dc3.append((a, b, c))
    if A_sample[a][c] and A_sample[c][b] and A_sample[b][a]:
        dc3.append((a, c, b))

print(f"\n  Directed 3-cycles (raw): {len(dc3)}")
for c in dc3:
    print(f"    {c[0]}→{c[1]}→{c[2]}→{c[0]}")

# In a tournament, each 3-vertex set has exactly 0 or 1 directed 3-cycle
# (not counting direction/starting point). So for 3-cycles, vertex-set
# counting IS correct.

# BUT for 5-cycles, a 5-vertex tournament can have multiple DISTINCT
# directed 5-cycles on the same vertex set!
dc5 = []
for perm in permutations(range(n)):
    if all(A_sample[perm[i]][perm[(i+1) % 5]] for i in range(5)):
        # Normalize: start at min vertex
        min_idx = perm.index(min(perm))
        rotated = perm[min_idx:] + perm[:min_idx]
        # Choose direction: compare rotated[1] vs rotated[-1]
        if rotated[1] < rotated[-1]:
            dc5.append(rotated)
        else:
            reversed_cycle = (rotated[0],) + rotated[1:][::-1]
            dc5.append(reversed_cycle)

dc5_unique = list(set(dc5))
print(f"\n  Directed 5-cycles (normalized): {len(dc5_unique)}")
for c in dc5_unique:
    print(f"    {c[0]}→{c[1]}→{c[2]}→{c[3]}→{c[4]}→{c[0]}")

# ======================================================================
# KEY QUESTION: Does the OCF use directed cycles or vertex sets?
# ======================================================================
print("\n" + "=" * 70)
print("OCF VERIFICATION: DIRECTED CYCLES vs VERTEX SETS")
print("=" * 70)

print(f"\n  The definition says: 'vertices are directed odd cycles of T'")
print(f"  So each DISTINCT directed cycle is a separate vertex of Ω(T).")
print(f"  For 3-cycles: 1 directed cycle per vertex set (in a tournament).")
print(f"  For 5-cycles: MULTIPLE directed cycles possible per vertex set!")

# ======================================================================
# Exhaustive n=5 test: try both interpretations
# ======================================================================
print("\n" + "=" * 70)
print("EXHAUSTIVE n=5: BOTH INTERPRETATIONS")
print("=" * 70)

n = 5
tb = n*(n-1)//2

# Method 1: vertex sets
# Method 2: directed cycles (normalized)

mm_vset = 0
mm_dcyc = 0
for bits in range(1 << tb):
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)

    # Vertex-set method (3-cycles: 1 per set, 5-cycles: 1 per set)
    cyc_vset = []
    for combo in combinations(range(n), 3):
        a, b, c = combo
        if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]):
            cyc_vset.append(frozenset(combo))
    for combo in combinations(range(n), 5):
        verts = list(combo)
        found = False
        for perm in permutations(verts[1:]):
            order = [verts[0]] + list(perm)
            if all(A[order[i]][order[(i+1) % 5]] for i in range(5)):
                found = True
                break
        if found:
            cyc_vset.append(frozenset(combo))

    ix_vs, _ = conflict_graph_and_ip(
        [tuple(sorted(c)) for c in cyc_vset], [2])
    if ix_vs[2] != H:
        mm_vset += 1

    # Directed cycle method
    all_dcycles = []
    # 3-cycles
    for combo in combinations(range(n), 3):
        a, b, c = combo
        if A[a][b] and A[b][c] and A[c][a]:
            all_dcycles.append((a, b, c))
        elif A[a][c] and A[c][b] and A[b][a]:
            all_dcycles.append((a, c, b))

    # 5-cycles (all of them, normalized)
    seen = set()
    for perm in permutations(range(n)):
        if all(A[perm[i]][perm[(i+1) % 5]] for i in range(5)):
            min_idx = list(perm).index(min(perm))
            rotated = perm[min_idx:] + perm[:min_idx]
            if rotated[1] < rotated[-1]:
                canon = rotated
            else:
                canon = (rotated[0],) + rotated[1:][::-1]
            seen.add(canon)
    all_dcycles.extend(list(seen))

    ix_dc, alpha_dc = conflict_graph_and_ip(all_dcycles, [2])
    if ix_dc[2] != H:
        mm_dcyc += 1
        if mm_dcyc <= 5:
            c3 = sum(1 for c in all_dcycles if len(c) == 3)
            c5 = sum(1 for c in all_dcycles if len(c) == 5)
            print(f"  DC MISMATCH bits={bits}: H={H}, I(2)={ix_dc[2]}, c3={c3}, c5={c5}")

print(f"\n  Vertex-set method: {mm_vset}/{1 << tb} mismatches")
print(f"  Directed-cycle method: {mm_dcyc}/{1 << tb} mismatches")

if mm_vset > 0 and mm_dcyc > 0:
    print(f"\n  NEITHER method works! Something deeper is wrong.")
    print(f"  Perhaps α₁ is not #cycles but something else?")
elif mm_vset > 0 and mm_dcyc == 0:
    print(f"\n  DIRECTED CYCLE method is correct!")
    print(f"  Must count distinct directed cycles, not just vertex sets.")
elif mm_vset == 0:
    print(f"\n  VERTEX SET method works (so directed cycles don't matter).")
else:
    print(f"\n  Both methods work (they agree at n=5).")

# Check a specific case
print("\n" + "=" * 70)
print("DETAILED CASE: bits=40")
print("=" * 70)

A = bits_to_adj(40, n)
H = count_ham_paths(A, n)
print(f"  H = {H}")

# Method 1: vertex sets
cyc_vset = []
for combo in combinations(range(n), 3):
    a, b, c = combo
    if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]):
        cyc_vset.append(frozenset(combo))
has_5cycle = False
for perm in permutations(range(n)):
    if all(A[perm[i]][perm[(i+1) % 5]] for i in range(5)):
        has_5cycle = True
        break
if has_5cycle:
    cyc_vset.append(frozenset(range(n)))

print(f"  Vertex sets: {[set(c) for c in cyc_vset]}")
print(f"  α₁ (vset) = {len(cyc_vset)}")

# Method 2: directed cycles
dc3 = []
for combo in combinations(range(n), 3):
    a, b, c = combo
    if A[a][b] and A[b][c] and A[c][a]:
        dc3.append((a, b, c))
    elif A[a][c] and A[c][b] and A[b][a]:
        dc3.append((a, c, b))

dc5 = set()
for perm in permutations(range(n)):
    if all(A[perm[i]][perm[(i+1) % 5]] for i in range(5)):
        min_idx = list(perm).index(min(perm))
        rotated = perm[min_idx:] + perm[:min_idx]
        if rotated[1] < rotated[-1]:
            dc5.add(rotated)
        else:
            dc5.add((rotated[0],) + rotated[1:][::-1])

dc5 = list(dc5)
all_dc = dc3 + dc5
print(f"  Directed 3-cycles: {dc3}")
print(f"  Directed 5-cycles: {dc5}")
print(f"  α₁ (directed) = {len(all_dc)}")

# Now compute I(Ω,2) with directed cycles
ix, alpha = conflict_graph_and_ip(all_dc, [2, 3, 6])
print(f"  α = {alpha}")
print(f"  I(Ω,2) = {ix[2]}, H = {H}, match = {ix[2] == H}")

# If still wrong, try: maybe each directed cycle counted TWICE (both directions)?
# In a tournament 3-cycle: only 1 direction. In a 5-cycle...
# Wait: in a tournament on 5 vertices, the number of distinct Hamiltonian cycles
# (directed) is either 0, 1, or more. Each one is a 5-cycle.
# Normalization: fix start, fix direction → each undirected cycle counted TWICE
# (clockwise and counterclockwise). But in a tournament, only ONE direction works
# for each cycle structure.

# Let me count WITHOUT normalization: just list all starting-point-fixed cycles
dc5_raw = []
for perm in permutations(range(1, n)):
    order = (0,) + perm
    if all(A[order[i]][order[(i+1) % 5]] for i in range(5)):
        dc5_raw.append(order)
print(f"\n  5-cycles starting at 0: {len(dc5_raw)}")
print(f"  Total directed 5-cycles (÷5 for rotation): {len(dc5_raw)}")
# Each directed cycle is counted exactly 5 times (once per starting vertex)
# So unique directed 5-cycles = len(dc5_raw) (since we fixed start=0)
# BUT: we need to check both directions. Each UNDIRECTED 5-cycle gives
# TWO directed ones. In a tournament, 0 or 1 of these directions works.
# Actually in a tournament, for a 5-vertex set with a Hamiltonian cycle,
# both the cycle and its reverse are valid iff... no, in a tournament
# if a→b then NOT b→a, so the reverse of a cycle NEVER works!
# So each undirected 5-cycle gives EXACTLY 1 directed 5-cycle.

# So # directed 5-cycles = # undirected 5-cycles = dc5_raw ÷ 1 (start fixed)
# But we fixed start at 0. Each 5-cycle passes through 0, and with start fixed,
# each cycle is counted once. But WAIT: not all 5-cycles pass through 0!
# At n=5, all vertices are used, so yes, all pass through 0.

print(f"  # directed 5-cycles (exact) = {len(dc5_raw)}")
for c in dc5_raw[:10]:
    print(f"    {'→'.join(str(v) for v in c)}→{c[0]}")

print("\nDone.")
