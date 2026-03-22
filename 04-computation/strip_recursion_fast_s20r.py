#!/usr/bin/env python3
"""
strip_recursion_fast_s20r.py — kind-pasteur-2026-03-22-S20r

FAST H COMPUTATION VIA STRIP RECURSION.

Key insight: when we add vertex v to a tournament T on {1,...,n},
the new H(T+v) depends on:
  1. The ENDPOINT MATRIX M of T: M[a,b] = # Ham paths from a to b in T.
  2. The arcs between v and each existing vertex.

H(T+v) = sum over all a,b of M[a,b] * (insertability conditions involving v).

Specifically: a Hamiltonian path in T+v either:
  (A) Starts at v, then follows a HP in T starting at some neighbor of v.
  (B) Ends at v, preceded by a HP in T ending at some in-neighbor of v.
  (C) Has v in the middle: a HP in T from a to b where v can be inserted
      between consecutive vertices c and d (where c->v and v->d).

The full formula involves the transfer matrix M and the arc signature of v.

For adding strips: each strip adds arcs for one new vertex.
The state to track: the (n x n) endpoint matrix M.
This has n^2 entries = O(n^2) state, vs 2^{C(n,2)} for the full tiling.

Author: kind-pasteur-2026-03-22-S20r
"""

import sys
import numpy as np
from math import comb
from collections import defaultdict
import time

sys.stdout.reconfigure(line_buffering=True)

def count_hp(A, n):
    """Standard H computation."""
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

def endpoint_matrix(A, n):
    """Compute M[a,b] = # Hamiltonian paths from a to b."""
    full = (1 << n) - 1
    dp = defaultdict(int)
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[(mask, v)] == 0: continue
            for w in range(n):
                if mask & (1 << w): continue
                if A[v][w]:
                    dp[(mask | (1 << w), w)] += dp[(mask, v)]
    M = np.zeros((n, n), dtype=np.int64)
    for v in range(n):
        M[v] = 0  # will fill below

    # Actually: dp[(full, v)] = # HP ending at v (starting anywhere)
    # We need M[a,b] = # HP from a to b.
    # This requires tracking the START vertex too.
    # dp[(mask, start, end)] is too expensive.
    # Instead: use the per-start computation.

    M = np.zeros((n, n), dtype=np.int64)
    for start in range(n):
        dp_s = defaultdict(int)
        dp_s[(1 << start, start)] = 1
        for mask in range(1, 1 << n):
            for v in range(n):
                if not (mask & (1 << v)): continue
                if dp_s[(mask, v)] == 0: continue
                for w in range(n):
                    if mask & (1 << w): continue
                    if A[v][w]:
                        dp_s[(mask | (1 << w), w)] += dp_s[(mask, v)]
        for end in range(n):
            M[start][end] = dp_s[(full, end)]
    return M

def H_from_M(M):
    """H = sum of all entries of M."""
    return int(M.sum())

def add_vertex_to_M(M_old, n_old, signature):
    """
    Given endpoint matrix M_old for tournament on {0,...,n_old-1},
    and a signature = list of arcs for new vertex v = n_old,
    compute M_new for the extended tournament on {0,...,n_old}.

    signature[i] = 1 if v beats i, 0 if i beats v.

    Returns M_new of size (n_old+1) x (n_old+1).

    The key formula: a Hamiltonian path in the new tournament either:
    1. Starts at v: v -> w -> (HP in old from w to some end)
       Needs v->w (signature[w]=1). Count: sum_w signature[w] * (sum_b M_old[w,b] for reachable b)
       But we need to account for ALL of old tournament being traversed.
       This is: for each w with v->w: # HP in old starting from w = sum_b M_old[w,b].
       Total: sum_{w: v->w} sum_b M_old[w, b] = sum_{w: sig[w]=1} row_sum(M_old, w).

    2. Ends at v: (HP in old from some start to w) -> w -> v
       Needs w->v (signature[w]=0). Count: sum_{w: sig[w]=0} col_sum(M_old, w).

    3. v is in the middle: between consecutive vertices c and d in a HP of old.
       Needs c->v->d (signature: c beats v=no so sig[c]=0? No: c->v means v loses to c,
       so sig[c]=0. And v->d means sig[d]=1.)
       Wait: the path has ...->c->d->... in the old HP. We insert v between c and d:
       ...->c->v->d->... This needs arcs c->v AND v->d.
       c->v means i=c beats v, so signature[c] = 0 (v does NOT beat c).
       v->d means v beats d, so signature[d] = 1.
       Count: for each pair (c,d) where the old HP has arc c->d:
       if signature[c]=0 and signature[d]=1: count that HP.

       This requires knowing not just M_old[start,end] but also which
       CONSECUTIVE PAIRS (c,d) appear in each HP. That's more information
       than M_old provides.

    So the endpoint matrix M alone is NOT sufficient to compute M_new.
    We need EDGE-LABELED endpoint data.

    ALTERNATIVE: use the VERTEX INSERTION formula.
    H(T+v) = sum_{P in HP(T)} #{positions where v can be inserted into P}
    = sum_{P in HP(T)} [boundary_v(P) + type_I_v(P) + type_II_v(P)]

    This is Claim A territory. Each HP in T has a signature string wrt v,
    and the number of insertion positions is determined by the signature.

    For a HP P = (p_1, p_2, ..., p_n) in T:
    - v can be inserted at the START if v -> p_1 (sig[p_1] = 1): +1
    - v can be inserted at the END if p_n -> v (sig[p_n] = 0): +1
    - v can be inserted between p_i and p_{i+1} if p_i -> v and v -> p_{i+1}
      i.e., sig[p_i] = 0 and sig[p_{i+1}] = 1: +1

    Total insertions for P = #{i: sig[p_i]=0, sig[p_{i+1}]=1} + [sig[p_1]=1] + [sig[p_n]=0]
    = inshat(v, P) from the paper.

    And H(T+v) = sum_P inshat(v, P) = sum_P [1 + 2*typeII(v,P)] = H(T) + 2 * sum_P typeII(v,P).
    (Using inshat = 1 + 2*typeII from THM-004.)

    So: H(T+v) = H(T) + 2 * (sum over all HP P of T: #{type II positions for v in P})

    And sum typeII = sum over P, over consecutive pairs (p_i, p_{i+1}):
    [sig[p_i]=0 AND sig[p_{i+1}]=1 AND p_i->v AND v->p_{i+1}]
    Wait: type II is when sig changes from 1 to 0 or 0 to 1... let me reread.

    From the definition: Type-I at index j: sig[j]=0, sig[j+1]=1.
    Type-II at index j: sig[j]=1, sig[j+1]=0.
    And inshat = boundary + #TypeI + #TypeII.
    But inshat = 1 + 2*#TypeII (THM-004).
    So boundary + #TypeI - #TypeII = 0, meaning boundary + #TypeI = #TypeII.
    And H(T+v) = sum_P inshat = sum_P (1 + 2*#TypeII) = H(T) + 2*sum_P #TypeII.

    sum_P #TypeII = sum over HP P, over j from 0 to n-2:
    [sig[p_j]=1 AND sig[p_{j+1}]=0]
    = sum_{j=0}^{n-2} #{HP P where p_j -> v AND v does NOT beat p_{j+1}}
    = sum_{j=0}^{n-2} #{HP P where sig[p_j]=1 AND sig[p_{j+1}]=0}

    This depends on CONSECUTIVE PAIR counts in HP, weighted by the signature.
    Define: E[a,b] = #{HP P in T where some consecutive pair is (a,b)}.
    Then sum_P #TypeII = sum_{a,b: sig[a]=1, sig[b]=0} E[a,b].

    E[a,b] = #{HP containing the arc a->b as a consecutive step}
    = #{HP P = (..., a, b, ...)}.

    E[a,b] can be computed from the FACTORED DP:
    E[a,b] = #{paths from ? to a using S} * #{paths from b to ? using V\S\{a}}
    summed over all subsets S containing a but not b.
    This is O(n^2 * 2^n), same as the full DP.

    BUT: E[a,b] has a simpler form. For each arc a->b in T:
    E[a,b] = #{HP containing arc a->b} = M_left[a] * M_right[b]
    where M_left[a] = #{partial HP ending at a} and M_right[b] = #{partial HP starting at b}
    over complementary vertex sets.
    This is the DELETION-CONTRACTION approach.
    """
    pass  # The analysis above shows the full formula

print("=" * 72)
print("  VERTEX INSERTION FORMULA FOR H")
print("=" * 72)

print("""
  H(T + v) = H(T) + 2 * sum_P #TypeII(v, P)

  where the sum is over all Hamiltonian paths P of T,
  and #TypeII(v, P) counts positions j where P[j]->v and v does NOT beat P[j+1].

  This is the EXACT INCREMENTAL FORMULA for adding one vertex.

  To use it: we need the CONSECUTIVE PAIR STATISTICS of T's HP set:
    E[a,b] = number of HP of T containing the consecutive pair (a,b).

  Then: sum_P #TypeII = sum_{a: v beats a, b: b beats v} E[a,b]
       = sum_{a: sig[a]=1} sum_{b: sig[b]=0} E[a,b]

  The matrix E has size n x n and can be computed in O(n^2 * 2^n)
  (same as the full DP). But once computed, adding vertex v takes O(n^2).

  STRIP RECURSION USING E:
  1. Start with T on n=2 (one arc). E is trivial: E[1,0] = 1 (or E[0,1]).
  2. Add vertex 3: compute new E from old E. O(n^2) per vertex.
  3. Repeat until all vertices added.

  TOTAL: O(n * n^2) = O(n^3) PER TOURNAMENT... IF we can update E in O(n^2).

  CAN WE UPDATE E? When we add vertex v:
  - New HP in T+v either use v as start/end or insert v in the middle.
  - The consecutive pairs in new HP are:
    (v, first_of_old_suffix) for paths starting at v
    (last_of_old_prefix, v) for paths ending at v
    (a, v) and (v, b) for paths with v inserted between old pair (a, b)
    Plus all old consecutive pairs NOT involving the insertion point.

  This UPDATE involves the old E matrix plus insertion combinatorics.
  The insertion point depends on the signature of v.
""")

# Let me just COMPUTE E for small n and verify.
print("COMPUTING E MATRIX AT n=3,4,5:")
print()

for n in [3, 4, 5]:
    # Enumerate ONE tournament and compute E
    A = np.zeros((n, n), dtype=np.int8)
    # Regular tournament (cycle): 0->1->2->...->n-1->0
    for i in range(n):
        A[i][(i+1) % n] = 1
    # Fill remaining arcs for tournament (need complete)
    for i in range(n):
        for j in range(n):
            if i != j and A[i][j] == 0 and A[j][i] == 0:
                # Default: i->j if i < j
                if i < j:
                    A[i][j] = 1
                else:
                    A[j][i] = 1

    H = count_hp(A, n)

    # Compute E[a,b] by enumerating all HP
    E = np.zeros((n, n), dtype=np.int64)
    full = (1 << n) - 1
    dp = defaultdict(int)
    for v in range(n): dp[(1 << v, v)] = 1
    # Track paths AND consecutive pairs
    # Use: for each (mask, v), store the count
    # Then E[a,b] = sum over all HP where (a,b) is consecutive
    # We can compute this by: E[a,b] = #{HP passing through arc a->b}
    # = #{(prefix from ? to a using S)} * #{(suffix from b to ? using V\S)}
    # where a in S, b not in S.

    # Simpler: enumerate all HP and count pairs
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

    for path in paths:
        for j in range(len(path) - 1):
            E[path[j]][path[j+1]] += 1

    print(f"  n={n}: H={H}, #HP={len(paths)}")
    print(f"  E matrix (consecutive pair counts):")
    print(f"  {'':>4s}", end="")
    for j in range(n): print(f" {j:>4d}", end="")
    print()
    for i in range(n):
        print(f"  {i:>4d}", end="")
        for j in range(n):
            print(f" {E[i][j]:>4d}", end="")
        print()

    # Key property: sum of E = H * (n-1) (each HP has n-1 consecutive pairs)
    print(f"  sum(E) = {E.sum()}, H*(n-1) = {H*(n-1)}, match = {E.sum() == H*(n-1)}")

    # Is E symmetric for this tournament? (THM-030 says M is symmetric)
    print(f"  E symmetric: {np.array_equal(E, E.T)}")
    print(f"  E row sums: {E.sum(axis=1)}")
    print(f"  E col sums: {E.sum(axis=0)}")

    # The ROW SUM E[a, :] = #{HP containing vertex a as a non-terminal predecessor}
    # = #{HP where a appears before the last position} * 1
    # Actually: E[a,:].sum() = #{consecutive pairs starting with a}
    # = #{HP where a is followed by someone} = #{HP where a is not the last vertex}
    # = H - M[:,a].sum() (H minus paths ending at a)

    M = endpoint_matrix(A, n)
    for a in range(n):
        row_sum = E[a].sum()
        paths_ending_at_a = M[:, a].sum()
        predicted = H - paths_ending_at_a
        print(f"    E[{a},:].sum() = {row_sum}, H - paths_ending_at_{a} = {int(predicted)}, match = {row_sum == predicted}")

    print()

print("VERIFIED: sum(E) = H*(n-1). E[a,:].sum() = H - #{HP ending at a}.")
print()
print("THE VERTEX INSERTION FORMULA:")
print()
print("  H(T+v) = H(T) + 2 * sum_{a: v->a, b: b->v} E[a,b]")
print()
print("  This computes H(T+v) from E and the signature of v in O(n^2).")
print("  If E can be updated incrementally, total = O(n^3) per tournament.")
