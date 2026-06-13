#!/usr/bin/env python3
"""
OUTLANDISH HYPOTHESIS K: Hidden symmetry group

The W-polynomial W(T,r) = (r-1/2)^{n-1} * F(T, (2r+1)/(2r-1))
is equivalent to the forward-edge polynomial F(T,x) via Mobius transform.

The palindromicity F_k = F_{n-1-k} means F(T,x) = x^{n-1} * F(T, 1/x).
In other words, F is invariant under x -> 1/x (mod power normalization).

Under the Mobius transform x = (2r+1)/(2r-1), the involution x -> 1/x becomes:
  1/x = (2r-1)/(2r+1) = (2(-r)+1)/(2(-r)-1)... no, that's (2r-1)/(2r+1).
  If x = (2r+1)/(2r-1), then 1/x = (2r-1)/(2r+1) = (2r'+1)/(2r'-1)
  where r' = -r. (Check: (2(-r)+1)/(2(-r)-1) = (1-2r)/(-1-2r) = (2r-1)/(2r+1).)

So palindromicity F(T,x) = x^{n-1} F(T,1/x) maps to:
  W(T,r) = W(T,-r) * something

Let's check: W(T,-r) = (-r-1/2)^{n-1} * F(T, (-2r+1)/(-2r-1))
                      = (-r-1/2)^{n-1} * F(T, (2r-1)/(2r+1))
                      = (-r-1/2)^{n-1} * F(T, 1/x)
                      = (-r-1/2)^{n-1} * x^{-(n-1)} * F(T, x) (palindrome)
                      = (-r-1/2)^{n-1} * ((2r-1)/(2r+1))^{n-1} * F(T, x)

Hmm. Let me just verify: does W(T,r) = (-1)^{n-1} W(T,-r)?

W(T,r) = sum_P prod(r + s_e)
W(T,-r) = sum_P prod(-r + s_e)

If n-1 is even: prod(-r+s_e) = prod(-(r-s_e)) = (-1)^{n-1} prod(r-s_e)
But prod(r-s_e) = prod(r + s'_e) where s'_e = -s_e = edge directions flipped
= W(T^op, r) where T^op is the reversed tournament.

So W(T,-r) = (-1)^{n-1} W(T^op, r).

By palindromicity (H(T)=H(T^op)), we know W(T,1/2) = W(T^op,1/2) = H.
But W(T,r) ≠ W(T^op,r) in general.

HOWEVER: the EVEN coefficients of W ARE the same for T and T^op!
W(T,r) = sum c_{2k} r^{2k} (only even powers by the even-power theorem).
W(T^op, r) = sum c'_{2k} r^{2k}.

Since W(T,-r) = (-1)^{n-1} W(T^op,r):
  sum c_{2k} (-r)^{2k} = (-1)^{n-1} sum c'_{2k} r^{2k}
  sum c_{2k} r^{2k} = (-1)^{n-1} sum c'_{2k} r^{2k}

So c_{2k}(T) = (-1)^{n-1} c_{2k}(T^op) for all k.

If n is odd: c_{2k}(T) = c_{2k}(T^op). The W-polynomial is the SAME for T and T^op.
If n is even: c_{2k}(T) = -c_{2k}(T^op). So c_{2k} = 0 when T ≅ T^op (self-converse).

Let's verify this.
"""
from itertools import permutations
import numpy as np
import random

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

def reverse_tournament(A):
    n = len(A)
    return [[A[j][i] for j in range(n)] for i in range(n)]

def compute_W_evals(A, n, r_vals):
    results = []
    for r in r_vals:
        W = 0
        for P in permutations(range(n)):
            prod_val = 1
            for i in range(n-1):
                s = 0.5 if A[P[i]][P[i+1]] else -0.5
                prod_val *= (r + s)
            W += prod_val
        results.append(W)
    return results

random.seed(42)

print("=== Verify W(T,-r) = (-1)^{n-1} W(T^op, r) ===")
for n in [4, 5, 6, 7]:
    print(f"\nn={n}:")
    for trial in range(3):
        A = random_tournament(n)
        A_op = reverse_tournament(A)

        r_vals = [0.0, 0.5, 1.0, 1.5, 2.0]
        W_T = compute_W_evals(A, n, r_vals)
        W_Top = compute_W_evals(A_op, n, r_vals)
        W_T_neg = compute_W_evals(A, n, [-r for r in r_vals])

        max_err = 0
        for i, r in enumerate(r_vals):
            predicted = (-1)**(n-1) * W_Top[i]
            actual = W_T_neg[i]
            err = abs(actual - predicted)
            max_err = max(max_err, err)

        print(f"  trial {trial}: max_err = {max_err:.2e}")

# ========== Test: W(T,r) = W(T^op, r) for odd n ==========
print("\n\n=== W(T) = W(T^op) for odd n? ===")
for n in [3, 5, 7]:
    print(f"\nn={n}:")
    for trial in range(3):
        A = random_tournament(n)
        A_op = reverse_tournament(A)

        r_vals = [0.0, 0.5, 1.0, 1.5, 2.0]
        W_T = compute_W_evals(A, n, r_vals)
        W_Top = compute_W_evals(A_op, n, r_vals)

        max_err = max(abs(a-b) for a,b in zip(W_T, W_Top))
        print(f"  trial {trial}: max|W(T)-W(T^op)| = {max_err:.2e}")

# ========== OUTLANDISH HYPOTHESIS L: The Z_2 x Z_2 structure ==========
print("\n\n=== HYPOTHESIS L: Tournament has Z_2 x Z_2 symmetry structure ===")
# There are TWO natural involutions on tournaments:
# 1. Complement (T^op): reverse all edges
# 2. GS flip: reverse non-backbone edges (tiling-dependent)
#
# For odd n, W(T) = W(T^op), so complement preserves the W-polynomial.
# GS flip changes H but preserves something else.
#
# What if there's a THIRD involution that, combined with these two,
# creates a Z_2 x Z_2 (Klein four-group) acting on the space of tournaments?

# The three involutions and their effects:
# 1. T^op: W(T^op,r) = (-1)^{n-1} W(T,-r). For odd n: W preserved.
# 2. GS flip: changes tiling bits, maps T -> T'. H changes, S changes.
# 3. T^op o GS: reverse non-backbone, then flip all. What does this do?

# Let's check what GS flip o complement does
n = 5

def tiling_to_tournament(bits, n):
    A = [[0]*n for _ in range(n)]
    for i in range(1, n):
        A[i][i-1] = 1
    tiles = [(a,b) for a in range(n) for b in range(a) if a-b >= 2]
    tiles.sort()
    for idx, (a, b) in enumerate(tiles):
        if (bits >> idx) & 1: A[b][a] = 1
        else: A[a][b] = 1
    return A

tiles = [(a,b) for a in range(n) for b in range(a) if a-b >= 2]
m = len(tiles)

print(f"  n={n}, m={m}")
print(f"  Backbone edges: {[(i, i-1) for i in range(1, n)]}")
print(f"  Non-backbone tiles: {tiles}")

# For bits=7: compute T, T^op, T_flip, T^op_flip
bits = 7
A = tiling_to_tournament(bits, n)
A_op = reverse_tournament(A)
flip_bits = bits ^ ((1 << m) - 1)
A_flip = tiling_to_tournament(flip_bits, n)
A_flip_op = reverse_tournament(A_flip)

from itertools import permutations as perms

def canonical(A):
    n = len(A)
    best = None
    for p in perms(range(n)):
        s = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
        if best is None or s < best:
            best = s
    return best

c = [canonical(X) for X in [A, A_op, A_flip, A_flip_op]]
print(f"\n  bits={bits}:")
print(f"    T ≅ T^op: {c[0]==c[1]}")
print(f"    T ≅ T_flip: {c[0]==c[2]}")
print(f"    T_flip ≅ T^op: {c[2]==c[1]}")
print(f"    T^op ≅ T_flip^op: {c[1]==c[3]}")

# What are the orbits under the group generated by complement and GS flip?
print("\n\n=== Orbits under <complement, GS-flip> ===")
from collections import defaultdict

orbit_map = {}
orbits = defaultdict(set)

for bits in range(1 << m):
    if bits in orbit_map:
        continue

    # Generate orbit
    orbit = set()
    stack = [bits]
    while stack:
        b = stack.pop()
        if b in orbit:
            continue
        orbit.add(b)

        # GS flip
        gs = b ^ ((1 << m) - 1)
        if gs not in orbit:
            stack.append(gs)

        # Complement: need to find which bits correspond to T^op
        # T^op reverses ALL edges. But in tiling representation,
        # backbone edges are FIXED (always pointing down).
        # So T^op may not even be representable as a tiling!
        # T^op has backbone edges reversed: A_op[i][i-1] = 0, A_op[i-1][i] = 1
        # This is NOT a valid tiling (backbone must go i+1 -> i).
        # So complement takes us OUTSIDE the tiling space.

    for b in orbit:
        orbit_map[b] = frozenset(orbit)

unique_orbits = set(orbit_map.values())
orbit_sizes = sorted([len(o) for o in unique_orbits])
print(f"  n={n}: {len(unique_orbits)} orbits under GS flip only")
print(f"  Orbit sizes: {orbit_sizes}")
# Under GS flip alone: each orbit has size 1 (fixed point) or 2 (pair)
fixed_points = sum(1 for o in unique_orbits if len(o) == 1)
pairs = sum(1 for o in unique_orbits if len(o) == 2)
print(f"  Fixed points: {fixed_points}, Pairs: {pairs}")

# Fixed points are bits where bits = bits ^ full_mask, i.e., full_mask = 0.
# Only possible if m = 0. Otherwise no fixed points.
# Wait, for m > 0: bits ^ full_mask = bits iff full_mask = 0 iff m = 0.
# So no fixed points for m > 0, all orbits have size 2. Check:
print(f"  (Expect: 0 fixed points, {(1<<m)//2} pairs)")

# ========== OUTLANDISH HYPOTHESIS M: Tournament as a QUANTUM STATE ==========
print("\n\n=== HYPOTHESIS M: Tournament as quantum state ===")
# Consider each tournament T on n vertices as a vector |T> in C^{2^m}
# where m = C(n,2) edges. The H function is a linear functional:
# <H|T> = H(T).
#
# The W-polynomial creates a "basis change" that separates H into layers.
# The palindromic structure is a SYMMETRY of this quantum system.
#
# What if the transfer matrix M is related to a QUANTUM WALK operator?
# M[a,b] = number of tilings (Hamiltonian paths through (a,b) as an edge)
# M being symmetric = time-reversal symmetry of the quantum walk.
#
# This is speculative but: the symmetry M[a,b] = M[b,a] that we're trying
# to prove might follow from a quantum-walk time-reversal argument.
print("  The transfer matrix M being symmetric is analogous to")
print("  time-reversal symmetry in a quantum walk on the tournament graph.")
print("  This is the most speculative hypothesis but could connect to")
print("  existing quantum walk literature.")

# ========== HYPOTHESIS N: The 1+2^{n-2} via binary tree structure ==========
print("\n\n=== HYPOTHESIS N: 1+2^{n-2} from binary recursion ===")
# bits=all tournament has H = 1+2^{n-2} for n <= 6 but fails at n=7 (H=31).
# 2^{n-2}+1: n=3→3, n=4→5, n=5→9, n=6→17, n=7→33 (but actual=31)
#
# What if the recursion is H(n) = 2*H(n-1) - correction?
# H(3)=3, H(4)=5, H(5)=9, H(6)=17, H(7)=31
# Ratios: 5/3=1.67, 9/5=1.80, 17/9=1.89, 31/17=1.82
# Differences: 2, 4, 8, 14 (= 2, 4, 8, 14)
# Expected diffs for 2^{n-2}+1: 2, 4, 8, 16
# Actual diff at n=7: 14 = 16-2
#
# What if the correction is related to the number of non-trivial automorphisms?
# Or to the number of 3-cycles in the bits=all tournament?

def tiling_to_tournament_global(bits, n):
    A = [[0]*n for _ in range(n)]
    for i in range(1, n):
        A[i][i-1] = 1
    tiles = [(a,b) for a in range(n) for b in range(a) if a-b >= 2]
    tiles.sort()
    for idx, (a, b) in enumerate(tiles):
        if (bits >> idx) & 1: A[b][a] = 1
        else: A[a][b] = 1
    return A

def ham_path_count_dp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or (mask, v) not in dp: continue
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

print("  bits=all tournament:")
for n in range(3, 9):
    tiles = [(a,b) for a in range(n) for b in range(a) if a-b >= 2]
    m = len(tiles)
    A = tiling_to_tournament_global((1 << m) - 1, n)

    H = ham_path_count_dp(A, n)
    t3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]: t3 += 1
                if A[i][k] and A[k][j] and A[j][i]: t3 += 1

    # Score sequence
    scores = sorted([sum(A[i]) for i in range(n)])
    print(f"  n={n}: H={H}, 1+2^(n-2)={1+2**(n-2)}, diff={H-1-2**(n-2)}, "
          f"t3={t3}, scores={scores}")

    # What tournament IS bits=all?
    # Backbone: i+1 -> i (all)
    # Non-backbone: b -> a for all a-b >= 2 (all upward arrows)
    # So edge a->b exists iff a=b+1 (backbone) or b<a and a-b>=2 (non-backbone upward)
    # In other words: A[b][a] = 1 iff a-b >= 2 (b beats a = "upward")
    # Combined with backbone A[a][a-1] = 1 (a beats a-1 = "downward one step")
    #
    # So vertex j beats vertex i iff: j = i+1 (backbone) or j < i and i-j >= 2 (upward)
    # Equivalently: j beats i iff j = i+1 or j < i-1
    # Equivalently: i beats j iff i > j+1 (skip beat) or i = j-1 (predecessor)
    #
    # Score of vertex v: out-degree
    # v beats u iff u = v-1 (backbone down) or u > v+1 (skip up... wait)
    # Let me just compute directly.
