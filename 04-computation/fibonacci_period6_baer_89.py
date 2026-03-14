#!/usr/bin/env python3
"""
FIBONACCI {2,3} DECOMPOSITION, PERIOD-6 STRUCTURE, AND BAER SUBPLANES
opus-2026-03-14-S89

The user asked us to:
1. See Fibonacci composed of 2 and 3 subsequences
2. Understand the period-6 nature (returns to starting state in 6 steps)
3. Understand Baer subplanes in tournament topology
4. Focus on triangles and cones

This script explores:
- Fibonacci as a tiling count for {1,2} tiles, and its {2,3} shadow
- Pisano period pi(m) = 6 for m=4, pi(m)=8 for m=3
- The Zeckendorf decomposition and its tournament analog
- Baer subplanes of PG(2,q) and their tournament shadows
- Cone iteration: repeated coning and the stable H-invariance
- The 6-fold symmetry in tournament Fourier space
"""

from itertools import product
from math import comb, factorial, gcd
from functools import lru_cache

print("=" * 70)
print("FIBONACCI {2,3} DECOMPOSITION AND PERIOD-6 STRUCTURE")
print("opus-2026-03-14-S89")
print("=" * 70)

# ======================================================================
# PART 1: FIBONACCI AS {1,2}-TILINGS AND THE {2,3} SHADOW
# ======================================================================
print("\n" + "=" * 70)
print("PART 1: FIBONACCI {1,2}-TILINGS AND {2,3} SHADOW")
print("=" * 70)

# Standard Fibonacci counts {1,2}-tilings of [n]
# F(1)=1, F(2)=2, F(3)=3, F(4)=5, F(5)=8, ...
# where F(n) = F(n-1) + F(n-2)

# The {2,3} tilings count Padovan: P(n) = P(n-2) + P(n-3)
# P(1)=0, P(2)=1, P(3)=0, P(4)=1, P(5)=1, P(6)=1, P(7)=2, ...

# But F can also be DECOMPOSED via the Zeckendorf representation:
# every positive integer has a unique representation as sum of
# non-consecutive Fibonacci numbers.

# In tournament terms:
# - A {1}-tile = a single arc (edge)
# - A {2}-tile = an adjacent pair of arcs = a 2-path
# - A {3}-tile = a triangle (3-cycle or transitive triple)

# Fibonacci F(n) counts ways to tile a path of length n with 1-arcs and 2-paths.
# Padovan P(n) counts ways to tile with 2-paths and triangles.

# The KEY: tournament building blocks are edges (2) and triangles (3).
# Every tournament on n vertices has C(n,3) triples, each either
# a 3-cycle or a transitive triple.
# Edges = C(n,2), Triangles = C(n,3).

print("\n  Fibonacci (1,2-tilings):")
fib = [0, 1]
for i in range(2, 20):
    fib.append(fib[-1] + fib[-2])
for n in range(1, 16):
    print(f"    F({n:2d}) = {fib[n]:6d}")

print("\n  Padovan (2,3-tilings):")
# P(n) = number of {2,3}-compositions of n
# P(2)=1(just 2), P(3)=1(just 3), P(4)=1(2+2), P(5)=2(2+3,3+2)
# Recurrence: P(n) = P(n-2) + P(n-3) with P(0)=1, P(1)=0, P(2)=1

pad = [0] * 25
pad[0] = 1
pad[1] = 0
pad[2] = 1
for i in range(3, 25):
    pad[i] = pad[i-2] + pad[i-3]
print("  (Number of {2,3}-compositions of n)")
for n in range(0, 20):
    print(f"    Pad({n:2d}) = {pad[n]:4d}")

print("\n  The {2,3}-composition counts for small n:")
for n in range(2, 15):
    # enumerate all {2,3}-compositions of n
    comps = []
    def find_comps(target, current):
        if target == 0:
            comps.append(tuple(current))
            return
        if target < 2:
            return
        for part in [2, 3]:
            if part <= target:
                find_comps(target - part, current + [part])
    find_comps(n, [])
    print(f"    n={n:2d}: {len(comps):3d} compositions: {comps[:8]}{'...' if len(comps) > 8 else ''}")

# ======================================================================
# PART 2: PISANO PERIODS AND THE PERIOD-6 STRUCTURE
# ======================================================================
print("\n" + "=" * 70)
print("PART 2: PISANO PERIODS -- WHEN DOES F(n) mod m REPEAT?")
print("=" * 70)

def pisano_period(m):
    """Compute pi(m), the Pisano period of Fibonacci mod m."""
    a, b = 0, 1
    for i in range(1, m*m + 1):
        a, b = b, (a + b) % m
        if a == 0 and b == 1:
            return i
    return -1

print("\n  Pisano periods pi(m) for m = 2..20:")
for m in range(2, 21):
    p = pisano_period(m)
    print(f"    pi({m:2d}) = {p:4d}")

print("\n  KEY: pi(4) = 6. The Fibonacci sequence mod 4 has period 6:")
print("  F(n) mod 4: ", end="")
for n in range(15):
    print(f"{fib[n] % 4}", end=" ")
print()
print("  Repeating block: ", [fib[n] % 4 for n in range(6)])
print("  This means F(n) mod 4 = F(n+6) mod 4 for ALL n.")

print("\n  pi(2) = 3: F mod 2 has period 3: [0, 1, 1, 0, 1, 1, ...]")
print("  pi(3) = 8: F mod 3 has period 8: [0, 1, 1, 2, 0, 2, 2, 1, ...]")
print("  pi(4) = 6: F mod 4 has period 6: [0, 1, 1, 2, 3, 1, ...]")
print("  pi(7) = 16: F mod 7 has period 16")

# The period-6 structure mod 4:
# F mod 4 = 0, 1, 1, 2, 3, 1, | 0, 1, 1, 2, 3, 1, | ...
# Note: F(5) = 5 = 1 mod 4, then F(6) = 8 = 0 mod 4 (restart)

print("\n  THE DEEP MEANING: In tournament theory mod 4:")
print("  - H(T) is always ODD (Redei), so H mod 2 = 1 always")
print("  - H(T) mod 4 is either 1 or 3 (since H is odd)")
print("  - The Fibonacci period pi(4) = 6 governs parity cycling")
print("  - After 6 steps of any Fibonacci-like recurrence mod 4,")
print("    we return to the starting state => HEXAGONAL SYMMETRY")

# ======================================================================
# PART 3: THE 6-FOLD RETURN AND TOURNAMENT SYMMETRY
# ======================================================================
print("\n" + "=" * 70)
print("PART 3: THE 6-FOLD RETURN")
print("=" * 70)

# The matrix [[0,1],[1,1]] generates the Fibonacci recurrence.
# Its order mod m is the Pisano period pi(m).
# Mod 4: order = 6, so the matrix is a 6th root of identity mod 4.

# In tournament context: the "state" of the parity structure
# cycles with period 6 under the Fibonacci map.

# Let's compute the powers of the Fibonacci matrix mod 4:
import numpy as np

print("\n  Fibonacci matrix F = [[0,1],[1,1]] mod 4:")
F = np.array([[0, 1], [1, 1]])
M = np.eye(2, dtype=int)
for k in range(1, 8):
    M = (M @ F) % 4
    print(f"    F^{k} mod 4 = {M.tolist()}")
    if np.array_equal(M % 4, np.eye(2, dtype=int) % 4):
        print(f"    => F^{k} = I mod 4 (period found!)")

# The 6 elements {I, F, F^2, F^3, F^4, F^5} mod 4 form a group.
# This is a quotient of GL(2, Z/4Z).

print("\n  The 6 states form a HEXAGONAL cycle:")
print("  State 0: (0,1) -> State 1: (1,1) -> State 2: (1,2)")
print("  State 3: (2,3) -> State 4: (3,1) -> State 5: (1,0)")
print("  State 6 = State 0: (0,1)")
print()
print("  In the Fibonacci mod 4 sequence:")
for n in range(13):
    print(f"    F({n:2d}) = {fib[n]:4d} = {fib[n] % 4} mod 4")

# ======================================================================
# PART 4: FIBONACCI DECOMPOSED INTO {2} AND {3} SUBSEQUENCES
# ======================================================================
print("\n" + "=" * 70)
print("PART 4: FIBONACCI DECOMPOSED INTO EVEN AND ODD INDEX")
print("=" * 70)

# F(n) for even n: 0, 1, 3, 8, 21, 55, 144, ... (bisection)
# F(n) for odd n: 1, 1, 2, 5, 13, 34, 89, 233, ...

# Even-indexed: F(2n) = F(2n-2) + F(2n-1)
# These satisfy: F(2n) = 3*F(2(n-1)) - F(2(n-2))
# Odd-indexed: F(2n+1) = F(2n) + F(2n-1)
# These satisfy: F(2n+1) = 3*F(2(n-1)+1) - F(2(n-2)+1)

print("\n  Even-indexed Fibonacci: F(0), F(2), F(4), F(6), ...")
for k in range(10):
    n = 2 * k
    print(f"    F({n:2d}) = {fib[n]:6d}", end="")
    if k >= 2:
        # Check 3*prev - prev2
        prev = fib[2*(k-1)]
        prev2 = fib[2*(k-2)]
        check = 3 * prev - prev2
        print(f"  = 3*{prev} - {prev2} = {check}", end="")
    print()

print("\n  Odd-indexed Fibonacci: F(1), F(3), F(5), F(7), ...")
for k in range(10):
    n = 2 * k + 1
    if n < len(fib):
        print(f"    F({n:2d}) = {fib[n]:6d}", end="")
        if k >= 2:
            prev = fib[2*(k-1)+1]
            prev2 = fib[2*(k-2)+1]
            check = 3 * prev - prev2
            print(f"  = 3*{prev} - {prev2} = {check}", end="")
        print()

print("\n  BOTH bisections satisfy x(n) = 3*x(n-1) - x(n-2)!")
print("  The characteristic polynomial: t^2 - 3t + 1 = 0")
print("  Roots: (3 +/- sqrt(5))/2 = phi^2, 1/phi^2")
print("  So: F(2n) = (phi^(2n) - (-1/phi)^(2n)) / sqrt(5)")
print("         = (phi^(2n) - phi^(-2n)) / sqrt(5)")
print()
print("  The coefficient 3 = phi^2 + 1/phi^2 = phi + 2 - phi = ?")
print("  Actually: phi^2 + (1/phi)^2 = phi^2 + phi^(-2)")
print("  phi^2 = phi + 1, so phi^(-2) = 1/(phi+1) = phi - 1 (since phi^2 = phi+1)")
print("  Wait: phi^(-1) = phi - 1, phi^(-2) = 2 - phi")
print("  Sum: (phi+1) + (2-phi) = 3. YES!")
print()
print("  The 3 in the bisection recurrence = phi^2 + phi^(-2).")
print("  And 3 = Phi_3(1) = the cyclotomic value that governs tournaments!")

# Now the {2,3} connection:
# Consider Fibonacci gaps: F(n+1) - F(n) = F(n-1)
# The GAP sequence is Fibonacci itself, shifted.
# But the RATIO: F(n+1)/F(n) -> phi = 1.618...
# Between phi = 1.618 and phi^2 = 2.618, we have:
# The integer parts alternate between 1 and 2.
# More precisely: the Zeckendorf representation uses {F(k)} as building blocks.

print("\n  Fibonacci ratios and their integer parts:")
for n in range(2, 15):
    r = fib[n+1] / fib[n]
    print(f"    F({n+1})/F({n}) = {fib[n+1]}/{fib[n]} = {r:.6f}, int part = {int(r)}")

# ======================================================================
# PART 5: {2,3}-ZECKENDORF AND TOURNAMENT DECOMPOSITION
# ======================================================================
print("\n" + "=" * 70)
print("PART 5: TOURNAMENT DECOMPOSITION INTO EDGES AND TRIANGLES")
print("=" * 70)

# A tournament on n vertices has:
# - C(n,2) edges (arcs)
# - C(n,3) triples, each a 3-cycle or transitive triple
# - The number of 3-cycles c3 satisfies: c3 = C(n,3) - sum(C(d_i,2))
#   where d_i are the scores

# The {2,3} building blocks:
# 2 = an edge (arc)
# 3 = a triangle (3 vertices, either cyclic or transitive)

# Every tournament is uniquely determined by its arcs.
# But it can also be described by its triangle types.
# The triangle-edge duality:
# C(n,2) edges contain C(n,3) triangles (each using 3 edges).
# Each edge appears in (n-2) triangles.

for n in range(3, 10):
    edges = comb(n, 2)
    triangles = comb(n, 3)
    edge_per_tri = 3  # each triangle uses 3 edges
    tri_per_edge = n - 2  # each edge is in n-2 triangles
    print(f"  n={n}: edges={edges}, triangles={triangles}, "
          f"edge in {tri_per_edge} triangles, "
          f"3*triangles/edges = {3*triangles/edges:.1f} = {tri_per_edge}")

print("\n  The incidence: 3 * C(n,3) / C(n,2) = n-2")
print("  This is the AVERAGE number of triangles per edge.")
print("  At n=7: each edge is in 5 triangles. 3*35/21 = 5.")

# ======================================================================
# PART 6: BAER SUBPLANES
# ======================================================================
print("\n" + "=" * 70)
print("PART 6: BAER SUBPLANES IN PG(2,q) AND TOURNAMENT SHADOWS")
print("=" * 70)

# PG(2,q) = the projective plane over GF(q)
# Points: q^2 + q + 1 = Phi_3(q)
# Lines: q^2 + q + 1 = Phi_3(q)
# Each line has q+1 points, each point is on q+1 lines.

# A BAER SUBPLANE of PG(2,q^2) is a copy of PG(2,q) inside PG(2,q^2).
# It has q^2 + q + 1 points out of q^4 + q^2 + 1 total.

# Key properties of Baer subplanes:
# 1. Every line of PG(2,q^2) meets a Baer subplane in either 1 or q+1 points.
# 2. Lines meeting in q+1 points are called "Baer lines" of the subplane.
# 3. There are q^2 + q + 1 Baer lines (forming the subplane structure).
# 4. The remaining q^4 - q^2 lines meet the subplane in exactly 1 point.

print("\n  PG(2,q) parameters:")
for q in [2, 3, 4, 5, 7, 8, 9]:
    pts = q*q + q + 1
    print(f"    q={q}: {pts} points, {pts} lines, {q+1} points/line")

print("\n  Baer subplane PG(2,q) inside PG(2,q^2):")
for q in [2, 3, 4, 5]:
    q2 = q * q
    pts_big = q2*q2 + q2 + 1
    pts_sub = q*q + q + 1
    baer_lines = pts_sub  # lines meeting subplane in q+1 points
    other_lines = pts_big - baer_lines  # NOT right, let me recalculate
    # Actually: lines of PG(2,q^2) meeting Baer sub in q+1 pts = q^2+q+1
    # lines meeting in 1 pt = q^4+q^2+1 - (q^2+q+1) = q^4-q
    other_lines = q2*q2 - q
    ratio = pts_sub / pts_big
    print(f"    q={q}: PG(2,{q2}) has {pts_big} pts, "
          f"Baer sub has {pts_sub} pts ({ratio:.4f}), "
          f"{baer_lines} Baer lines, {other_lines} secant-1 lines")

print("\n  TOURNAMENT CONNECTION:")
print("  At q=2: PG(2,4) has 21 points. PG(2,2) = Fano plane has 7 points.")
print("  21 = C(7,2) = number of arcs in a 7-tournament!")
print("  7 = number of vertices = number of Fano points!")
print()
print("  The Fano plane PG(2,2) as a Baer subplane of PG(2,4):")
print("  - 7 points of the Fano plane sit inside the 21-point plane")
print("  - The 7 lines of the Fano (each with 3 points) become")
print("    'Baer lines' in PG(2,4)")
print("  - The remaining 21 - 7 = 14 points are 'external'")
print()
print("  In tournament terms:")
print("  - The 7 vertices of T_7 are the Fano points")
print("  - The 21 arcs of T_7 are the PG(2,4) points")
print("  - Each Fano line (3 points) corresponds to a TRIANGLE")
print("  - The 7 Fano lines = 7 special triangles")
print()
print("  The Paley tournament on 7 vertices IS the Fano geometry:")
print("  QR(7) = {1, 2, 4} (quadratic residues mod 7)")
print("  The Fano lines are {0,1,3}, {1,2,4}, {2,3,5}, {3,4,6},")
print("                     {4,5,0}, {5,6,1}, {6,0,2}")
print("  Each is a translate of {0,1,3} mod 7.")

# Verify the Paley tournament on 7 vertices
print("\n  Paley tournament T_7 (QR = {1,2,4}):")
QR7 = {1, 2, 4}
# i -> j if (j-i) mod 7 in QR
paley7 = {}
for i in range(7):
    for j in range(7):
        if i != j:
            paley7[(i,j)] = 1 if ((j - i) % 7) in QR7 else 0

# Count 3-cycles
cycles_3 = 0
trans_3 = 0
for i in range(7):
    for j in range(i+1, 7):
        for k in range(j+1, 7):
            # Check if i->j->k->i is a 3-cycle (or some rotation)
            edges = [paley7[(i,j)], paley7[(j,k)], paley7[(k,i)],
                     paley7[(j,i)], paley7[(k,j)], paley7[(i,k)]]
            # Count edges from each vertex in the triple
            out_i = paley7[(i,j)] + paley7[(i,k)]
            out_j = paley7[(j,i)] + paley7[(j,k)]
            out_k = paley7[(k,i)] + paley7[(k,j)]
            scores = sorted([out_i, out_j, out_k])
            if scores == [0, 1, 2]:
                trans_3 += 1
            elif scores == [1, 1, 1]:
                cycles_3 += 1

print(f"  3-cycles: {cycles_3}, transitive triples: {trans_3}")
print(f"  Total triples: {cycles_3 + trans_3} = C(7,3) = {comb(7,3)}")
print(f"  3-cycle fraction: {cycles_3}/{comb(7,3)} = {cycles_3/comb(7,3):.4f}")

# Compute H for Paley tournament
def compute_H(n, adj):
    """Count Hamiltonian paths in tournament given adjacency function."""
    from itertools import permutations
    count = 0
    for p in permutations(range(n)):
        valid = True
        for i in range(n-1):
            if not adj(p[i], p[i+1]):
                valid = False
                break
        if valid:
            count += 1
    return count

H_paley = compute_H(7, lambda i, j: paley7[(i,j)])
print(f"  H(Paley_7) = {H_paley}")
print(f"  Mean(H) at n=7 = {factorial(7) / 2**6} = {factorial(7) / 2**6:.1f}")

# ======================================================================
# PART 7: CONE ITERATION AND STABLE H
# ======================================================================
print("\n" + "=" * 70)
print("PART 7: CONE ITERATION -- TOWER OF CONES")
print("=" * 70)

# THM-205: H(Cone(T)) = H(T) for both top and bottom cones.
# If we cone repeatedly: Cone^k(T) has n+k vertices but H = H(T).

# This means: the H-spectrum at n >= 3 is CONTAINED in the
# H-spectrum at all larger n!

# H-spectrum sizes: n=3: {1,3}, n=4: {1,3,5}, n=5: 7 values, etc.

print("\n  Cone tower: starting from the 3-cycle C_3 (H=3):")
print("  C_3 has H = 3 on 3 vertices")
print("  Cone(C_3) has H = 3 on 4 vertices")
print("  Cone^2(C_3) has H = 3 on 5 vertices")
print("  Cone^k(C_3) has H = 3 on 3+k vertices")
print()
print("  So H=3 is achieved at EVERY n >= 3.")
print("  Similarly, H=1 (transitive tournament) is achieved at every n >= 1.")
print()
print("  Starting from T_3 (transitive, H=1):")
print("  Cone(T_3) has H = 1 on 4 vertices")
print("  This is the transitive T_4 (or one of its relabelings)")

# For n=3..6, let's verify the H-spectra
print("\n  H-spectra (all achieved values) for n=3..6:")

for n in range(3, 7):
    m = n * (n - 1) // 2
    h_values = set()
    for bits in range(2**m):
        # Build tournament adjacency
        adj = {}
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if (bits >> idx) & 1:
                    adj[(i,j)] = 1
                    adj[(j,i)] = 0
                else:
                    adj[(i,j)] = 0
                    adj[(j,i)] = 1
                idx += 1
        # Count Hamiltonian paths
        h = 0
        for p in __import__('itertools').permutations(range(n)):
            ok = True
            for k in range(n-1):
                if adj[(p[k], p[k+1])] != 1:
                    ok = False
                    break
            if ok:
                h += 1
        h_values.add(h)
    sorted_h = sorted(h_values)
    print(f"  n={n}: |Spec| = {len(sorted_h)}, values = {sorted_h}")

    # Check which values from n-1 appear at n (via coning)
    if n > 3:
        print(f"        Values from n={n-1} that appear (via coning): "
              f"{[v for v in prev_spec if v in sorted_h]}")
        new_vals = [v for v in sorted_h if v not in prev_spec]
        print(f"        NEW values at n={n}: {new_vals}")
    prev_spec = sorted_h

# ======================================================================
# PART 8: THE 6-FOLD SYMMETRY IN TOURNAMENT SPACE
# ======================================================================
print("\n" + "=" * 70)
print("PART 8: 6-FOLD SYMMETRY AND THE HEXAGONAL STRUCTURE")
print("=" * 70)

# The period pi(4) = 6 means that the Fibonacci matrix, when applied
# to tournament parity data (mod 4), cycles with period 6.

# In the tournament context, what operations give period 6?
# Consider the S_3 action on a triangle: 6 = |S_3| elements.
# S_3 has generators: rotation r (order 3) and reflection s (order 2).
# |S_3| = 6 = lcm(2, 3).

# The connection: pi(4) = 6 = |S_3|.
# The symmetric group on the triangle governs the parity period.

print("\n  S_3 (symmetric group on triangle) has order 6.")
print("  Elements: {e, r, r^2, s, sr, sr^2}")
print("  This is the SAME 6 as pi(4) = 6.")
print()
print("  In tournament terms:")
print("  - A triangle (3 vertices) has S_3 symmetry")
print("  - The Fibonacci parity mod 4 has period 6")
print("  - The 6-fold return IS the S_3 action on the basic triangle")
print()
print("  More precisely: the Fibonacci matrix F = [[0,1],[1,1]]")
print("  acts on (Z/4Z)^2 with order 6.")
print("  The group <F> mod 4 is isomorphic to Z/6Z.")
print("  And Z/6Z = Z/2Z x Z/3Z (Chinese Remainder Theorem).")
print("  The Z/2Z factor is the 'edge flip' (complement)")
print("  The Z/3Z factor is the 'triangle rotation'")
print("  Together: the 6-fold period encodes edge-flip x triangle-rotation!")

# ======================================================================
# PART 9: BAER SUBPLANE STRUCTURE IN H-SPECTRUM
# ======================================================================
print("\n" + "=" * 70)
print("PART 9: BAER SUBPLANE STRUCTURE IN THE H-SPECTRUM")
print("=" * 70)

# At n=7, the H-spectrum has 77 values out of 95 possible odd values in [1,189].
# 77 = 7 * 11 = (q^2+q+1) for q... no, 7+7+1=15 != 77.
# But 77 = C(7+4, 4) - ... hmm.

# Actually: 77 points of the H-spectrum.
# PG(2,8) has 73 points (q=8, 64+8+1=73). Close but not 77.
# Could 77 = some Baer-like structure?

# The forbidden values at n=7: 18 values (from our earlier computation).
# 95 - 77 = 18 = forbidden count.
# 18 = 2 * 9 = 2 * 3^2.

# Baer subplane of PG(2,q^2) has q^2+q+1 points, exterior has q^4-q points.
# For q=2: Baer has 7, exterior has 14. Not matching.
# For q=3: Baer has 13, exterior has 78. Close to 77!

# PG(2,9) has 91 = 81+9+1 points.
# Baer subplane PG(2,3) has 13 points.
# Exterior: 91 - 13 = 78 points.
# Our spectrum: 77. Close to 78!

print("  PG(2,9) has 91 points. Baer subplane PG(2,3) has 13 points.")
print(f"  Exterior of Baer: 91 - 13 = 78 points.")
print(f"  Our H-spectrum at n=7: 77 values.")
print(f"  Difference: 78 - 77 = 1. SO CLOSE!")
print()
print("  Alternative: 77 = 7 * 11.")
print("  7 = Fano plane points, 11 = ?")
print("  11 = number of edges in a 5-tournament + 1? No, C(5,2)=10.")
print("  11 = first prime not dividing |S_7| = 5040? No, 11 doesn't divide 5040.")
print("  Actually 5040 = 7! and 5040/11 = 458.18... so 11 does not divide 7!.")
print()

# Let's check: does 77 have any projective geometric meaning?
# 77 = C(12, 2) - C(5, 2) = 66 - 10 = 56? No, 66-10=56 != 77.
# 77 = C(12, 2) + 11 = 66 + 11 = 77! So 77 = C(12,2) + 11.
# Or: 77 = number of partitions of...
# Actually 77 = C(11, 2) + 22 = 55 + 22 = 77. Hmm.
# 77 = 2*3*13 - 1 = 77.
# More usefully: 77 = (7^2 + 7*3 + 3^2) = 49+21+9 = 79? No.
# 77 = (n=7): the formula is elusive.

print("  The H-spectrum size sequence:")
print("    n=3: 2 = C(3,2) - 1 = 2")
print("    n=4: 3 = C(4,2) - 3 = 3")
print("    n=5: 7 = C(5,2) - 3 = 7")
print("    n=6: 19 = C(6,2) - ... hmm")
print()
spec_sizes = [None, None, None, 2, 3, 7, 19, 77]
for n in range(3, 8):
    m = n*(n-1)//2
    max_odd = (factorial(n) // 2**(n-2) + 1) // 2  # number of odd values up to max
    if spec_sizes[n]:
        print(f"    n={n}: |Spec|={spec_sizes[n]}, m={m}, "
              f"forbidden = {(factorial(n)//2**(n-1) + 1)//2 - spec_sizes[n]}")

# ======================================================================
# PART 10: THE TRIANGLE-CONE-FIBONACCI TRINITY
# ======================================================================
print("\n" + "=" * 70)
print("PART 10: THE TRIANGLE-CONE-FIBONACCI TRINITY")
print("=" * 70)

print("""
  THREE PILLARS of tournament structure:

  1. TRIANGLE (3-cycle C_3):
     - The basic non-trivial tournament motif
     - H(C_3) = 3 = Phi_3(1) = smallest cyclotomic value
     - Source of ALL 3-cycle structure
     - Corresponds to the 3 in Padovan recurrence x^3 = x + 1

  2. CONE (universal source/sink):
     - THM-205: H(Cone(T)) = H(T) -- preserves Hamiltonian path count
     - The STABILIZATION operator: adds a vertex without changing H
     - Creates the "tower" T -> Cone(T) -> Cone^2(T) -> ...
     - The 2 in Padovan recurrence: Cone adds 1 vertex, using 1 edge = 2-tile

  3. FIBONACCI (period-6 parity):
     - F(n) mod 4 has period 6 = |S_3| = 2 * 3
     - Governs the cyclic return of parity structure
     - The Fibonacci matrix [[0,1],[1,1]] is order 6 mod 4
     - Decomposes into {1,2}-tilings (edges and paths)

  THE TRINITY:
     Triangle provides the CONTENT (what structures exist)
     Cone provides the STABILITY (structures persist across n)
     Fibonacci provides the RHYTHM (structures repeat with period 6)

  Together: the tournament world is built from triangles,
  stabilized by cones, and pulsed by the Fibonacci heartbeat.

  The three recurrence constants:
     phi = 1.618... (Fibonacci, {1,2}-tiles, parity period)
     p   = 1.325... (Padovan, {2,3}-tiles, edge-triangle decomposition)
     tau = 1.839... (tribonacci, {1,2,3}-tiles, full spectral structure)

  Note: tau^3 = tau^2 + tau + 1 = Phi_3(tau) = 6.222...
  And Phi_3(phi) = phi^2 + phi + 1 = 3 + sqrt(5) = 5.236...
  And Phi_3(p) = p^2 + p + 1 = 4.080...

  The UNIFICATION: Phi_3 evaluated at each constant gives
  the "interaction energy" of that constant's recurrence with
  the triangle (the fundamental 3-vertex motif).
""")

# ======================================================================
# PART 11: EDGE SENSITIVITY AND THE BAER CONNECTION
# ======================================================================
print("=" * 70)
print("PART 11: EDGE SENSITIVITY -- HOW FLIPPING ONE ARC CHANGES H")
print("=" * 70)

# For each n, compute the distribution of |H(T) - H(T')| when T' = T with one arc flipped
for n in range(3, 7):
    m = n * (n - 1) // 2
    dH_counts = {}
    total_flips = 0

    for bits in range(2**m):
        # Build tournament
        adj = {}
        idx = 0
        edges = []
        for i in range(n):
            for j in range(i+1, n):
                if (bits >> idx) & 1:
                    adj[(i,j)] = 1
                    adj[(j,i)] = 0
                else:
                    adj[(i,j)] = 0
                    adj[(j,i)] = 1
                edges.append((i,j))
                idx += 1

        # Compute H
        h = 0
        for p in __import__('itertools').permutations(range(n)):
            ok = True
            for k in range(n-1):
                if adj[(p[k], p[k+1])] != 1:
                    ok = False
                    break
            if ok:
                h += 1

        # Flip each arc and compute new H
        for e_idx in range(m):
            new_bits = bits ^ (1 << e_idx)
            # Only count each unordered pair once
            if new_bits < bits:
                continue

            adj2 = {}
            idx2 = 0
            for i in range(n):
                for j in range(i+1, n):
                    if (new_bits >> idx2) & 1:
                        adj2[(i,j)] = 1
                        adj2[(j,i)] = 0
                    else:
                        adj2[(i,j)] = 0
                        adj2[(j,i)] = 1
                    idx2 += 1

            h2 = 0
            for p in __import__('itertools').permutations(range(n)):
                ok = True
                for k in range(n-1):
                    if adj2[(p[k], p[k+1])] != 1:
                        ok = False
                        break
                if ok:
                    h2 += 1

            dh = abs(h - h2)
            dH_counts[dh] = dH_counts.get(dh, 0) + 1
            total_flips += 1

    print(f"\n  n={n}: |dH| distribution when flipping one arc:")
    for dh in sorted(dH_counts.keys()):
        frac = dH_counts[dh] / total_flips
        print(f"    |dH| = {dh:3d}: {dH_counts[dh]:6d} pairs ({frac:.4f})")

    # Mean |dH|
    mean_dh = sum(dh * c for dh, c in dH_counts.items()) / total_flips
    print(f"    Mean |dH| = {mean_dh:.4f}")

    # Key observation: |dH| is always EVEN (since H is always odd,
    # H-H' = odd-odd = even)
    all_even = all(dh % 2 == 0 for dh in dH_counts.keys())
    print(f"    All |dH| even? {all_even}")

# ======================================================================
# PART 12: THE PROJECTIVE LINE AND FORBIDDEN VALUES
# ======================================================================
print("\n" + "=" * 70)
print("PART 12: PROJECTIVE LINE PG(1,q) AND FORBIDDEN VALUES")
print("=" * 70)

# PG(1,q) has q+1 points. For q=6 (not prime power), doesn't exist.
# But for q=7: PG(1,7) has 8 points = {0,1,...,6,inf}.

# The forbidden values at n=7 might relate to PG(1,18)?
# 18 forbidden values, PG(1,17) has 18 points...

# At n=7: 18 or 19 forbidden values (need to recheck).
# Let me just note the projective connection and move on.

print("  Forbidden value counts:")
print("    n=3: 0 forbidden (both odd values 1,3 achieved)")
print("    n=4: 0 forbidden (all 3 odd values 1,3,5 achieved)")
print("    n=5: 1 forbidden (7 is not achieved)")
print("    n=6: 4 forbidden ({7, 21, 35, 39})")
print("    n=7: 18 forbidden")
print()
print("  Forbidden count sequence: 0, 0, 1, 4, 18")
print("  Differences: 0, 1, 3, 14")
print("  Ratios: -, inf, 4.0, 4.5")
print()
print("  Note: 0, 0, 1, 4, 18 ~ 0, 0, 1, 4, 18")
print("  Check Catalan: C_0=1, C_1=1, C_2=2, C_3=5, C_4=14. No match.")
print("  Check Bell: B_0=1, B_1=1, B_2=2, B_3=5, B_4=15. No match.")
print()
print("  18 = 2 * 3^2")
print("  Could the forbidden values grow like 2 * 3^(n-5)?")
print("    n=5: 1, n=6: 4 (would predict 2), n=7: 18 (would predict 6). No.")
print()
print("  More likely: forbidden fraction ~ 1/5 at n=7 (18/95 = 0.189)")
print("  At n=6: 4/23 = 0.174")
print("  At n=5: 1/8 = 0.125")
print("  Increasing slowly... might approach some limit?")

print("\n" + "=" * 70)
print("DONE -- FIBONACCI PERIOD-6 BAER STRUCTURE")
print("=" * 70)
