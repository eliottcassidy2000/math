#!/usr/bin/env python3
"""
fast_graph_sequences_s169.py — Fast Formulas for Graph-Structure Sequences
opus-2026-03-22-S169

The goal: sequences that capture the UNDERLYING GRAPH STRUCTURE
of tournaments, computable for ALL n without exhaustive enumeration.

The key insight: many quantities averaged over ALL tournaments
have CLOSED FORMS because the average tournament is "random"
(each arc independently 50/50).

For a random tournament on n vertices:
  E[s_v] = (n-1)/2 for all v
  E[s_v²] = (n-1)(2n-1)/12 + ((n-1)/2)²
  E[c₃] = C(n,3)/4 (each triple has prob 1/4 of being a cycle)
  E[H] = n!/2^{n-1}

These give EXACT formulas that work for ALL n.
"""

from math import comb, factorial, gcd, prod
from fractions import Fraction
from collections import Counter
import time

print("=" * 72)
print("  FAST FORMULAS FOR GRAPH-STRUCTURE SEQUENCES")
print("  opus-2026-03-22-S169")
print("=" * 72)
print()

# ==========================================================================
# SEQUENCE A: E[c₃] = C(n,3)/4 — AVERAGE 3-CYCLE COUNT
# ==========================================================================

print("SEQUENCE A: E[c₃] = C(n,3)/4 — Average 3-cycles")
print("-" * 60)
print()

print(f"  Each triple of vertices forms a 3-cycle with probability 2/8 = 1/4.")
print(f"  (2 cyclic orientations out of 8 possible.)")
print(f"  E[c₃] = C(n,3) × 1/4 = n(n-1)(n-2)/24")
print()

seq_ec3 = []
for n in range(2, 20):
    ec3 = Fraction(comb(n, 3), 4)
    seq_ec3.append(ec3)
    if n <= 12:
        print(f"  n={n:2d}: E[c₃] = {ec3} = {float(ec3):.2f}")

print(f"\n  Numerator sequence: {[s.numerator for s in seq_ec3[:15] if s.denominator in [1,2,3,4,6,8,12,24]]}")
print(f"  Integer values at n ≡ 0,1 (mod 4): ", end="")
print([int(Fraction(comb(n,3),4)) for n in range(2,20) if comb(n,3) % 4 == 0])
print()

# ==========================================================================
# SEQUENCE B: E[H] = n!/2^{n-1} — AVERAGE HAMILTONIAN PATHS
# ==========================================================================

print()
print("SEQUENCE B: E[H] = n!/2^{n-1}")
print("-" * 60)
print()

seq_eh = []
for n in range(1, 16):
    eh = Fraction(factorial(n), 1 << (n-1))
    seq_eh.append(eh)
    print(f"  n={n:2d}: E[H] = {eh} = {float(eh):.4f}")

print(f"\n  As integers × 2^k: {[factorial(n) for n in range(1,10)]}")
print(f"  This is n! / 2^{{n-1}} = A001700(n-1) related")
print()

# ==========================================================================
# SEQUENCE C: E[c₅] — AVERAGE 5-CYCLE COUNT
# ==========================================================================

print()
print("SEQUENCE C: E[c₅] — Average directed 5-cycles")
print("-" * 60)
print()

# A directed 5-cycle on 5 specific vertices:
# Choose 5 vertices: C(n,5) ways.
# Choose a cyclic ordering: (5-1)!/2 = 12 directed cycles per vertex set.
# Each has probability (1/2)^5 = 1/32 of being realized.
# So E[c₅] = C(n,5) × 12 × 1/32 = C(n,5) × 3/8

print(f"  For 5 specific vertices: 12 directed 5-cycles, each prob 1/32.")
print(f"  E[c₅] = C(n,5) × 12/32 = C(n,5) × 3/8")
print()

seq_ec5 = []
for n in range(2, 16):
    ec5 = Fraction(comb(n, 5) * 3, 8)
    seq_ec5.append(ec5)
    if n <= 12:
        print(f"  n={n:2d}: E[c₅] = {ec5} = {float(ec5):.2f}")

print()

# ==========================================================================
# SEQUENCE D: E[c_{2k+1}] — AVERAGE DIRECTED ODD CYCLES OF LENGTH 2k+1
# ==========================================================================

print()
print("SEQUENCE D: E[c_{2k+1}] — General odd cycle average")
print("-" * 60)
print()

# For a directed cycle of length m on m specific vertices:
# Number of directed cycles = (m-1)!/2 (for odd m) or (m-1)! (but halve for direction)
# Wait: directed cycles on m labeled vertices = (m-1)! (fix one vertex, permute rest, one direction)
# Actually: the number of DIRECTED Hamiltonian cycles on K_m = (m-1)!
# Each directed cycle uses m arcs, each independently correct with prob 1/2.
# So P(specific directed cycle exists) = (1/2)^m
# Number of directed cycles on m labeled vertices = (m-1)!
# E[c_m on specific m vertices] = (m-1)! × (1/2)^m
# E[c_m total] = C(n,m) × (m-1)! / 2^m = C(n,m) × (m-1)! / 2^m

print(f"  E[c_m] = C(n,m) × (m-1)! / 2^m  for directed m-cycles")
print()

for m in [3, 5, 7, 9, 11]:
    coeff = Fraction(factorial(m-1), 2**m)
    print(f"  m={m:2d}: coefficient = (m-1)!/2^m = {coeff} = {float(coeff):.6f}")
    seq = [Fraction(comb(n, m) * factorial(m-1), 2**m) for n in range(m, min(m+8, 20))]
    print(f"    E[c_{m}](n={m}..{m+7}): {[float(s) for s in seq]}")
    print()

# ==========================================================================
# SEQUENCE E: Var(c₃) — VARIANCE OF 3-CYCLE COUNT
# ==========================================================================

print()
print("SEQUENCE E: Var(c₃)")
print("-" * 60)
print()

# Var(c₃) = E[c₃²] - E[c₃]²
# c₃ = Σ_{triples} X_t where X_t = 1 if triple t is a 3-cycle
# E[X_t] = 1/4, Var(X_t) = 3/16
# Cov(X_t, X_s) depends on how many vertices t and s share:
#   0 shared: independent → Cov = 0
#   1 shared: Cov = E[X_t X_s] - 1/16
#   2 shared: dependent (share an edge)

# E[X_t X_s] when 1 vertex shared:
# Two triples sharing 1 vertex involve 5 distinct vertices.
# Prob(both cycles) = P(cycle on {a,b,c} AND cycle on {a,d,e})
# These share vertex a. The two triples are independent conditional on
# a's arcs to {b,c} and {d,e}. But the arcs {b,c,d,e} don't overlap.
# Actually: the 5 arcs of the first triple and 5 arcs of the second
# share exactly 0 arcs (different pairs). So they're independent!
# Wait: triple {a,b,c} uses arcs ab, bc, ca. Triple {a,d,e} uses ad, de, ea.
# No arc overlap (all different pairs). Independent!
# So Cov(X_t, X_s) = 0 when 1 shared vertex.

# E[X_t X_s] when 2 shared vertices:
# Triples sharing 2 vertices: {a,b,c} and {a,b,d}.
# Arcs: first uses ab,bc,ca. Second uses ab,bd,da.
# They share arc ab! So NOT independent.
# P(both cycles) = P(a→b→c→a AND a→b→d→a) + similar
# This requires careful counting.

# Let me just compute exactly for small n
print(f"  COMPUTING Var(c₃) EXACTLY:")

def c3_count(adj, n):
    count = 0
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                if (adj[a][b] and adj[b][c] and adj[c][a]) or \
                   (adj[a][c] and adj[c][b] and adj[b][a]):
                    count += 1
    return count

def tournament_adj(n, bits):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            idx += 1
    return adj

seq_varc3 = []
for n in range(3, 8):
    if n > 6:
        print(f"  n={n}: (skipping, too large)")
        continue
    t0 = time.time()
    c3_vals = []
    for bits in range(1 << comb(n, 2)):
        c3_vals.append(c3_count(tournament_adj(n, bits), n))
    total = len(c3_vals)
    ec3 = Fraction(sum(c3_vals), total)
    ec3_sq = Fraction(sum(v*v for v in c3_vals), total)
    var_c3 = ec3_sq - ec3 * ec3
    seq_varc3.append(var_c3)
    elapsed = time.time() - t0
    print(f"  n={n}: E[c₃]={ec3}, E[c₃²]={ec3_sq}, Var(c₃)={var_c3} = {float(var_c3):.4f} ({elapsed:.1f}s)")

print()

# Can we get a formula for Var(c₃)?
# Var(c₃) = Σ_{triples t} Var(X_t) + 2 Σ_{t<s, |t∩s|=2} Cov(X_t, X_s)
# = C(n,3) × 3/16 + 2 × (# pairs sharing 2 vertices) × Cov(X_t, X_s | 2 shared)

# # pairs of triples sharing exactly 2 vertices:
# Choose the 2 shared: C(n,2). Choose the 2 non-shared: C(n-2, 2).
# But each pair is counted once. Actually: choose 2 shared {a,b}, then
# choose c ≠ a,b and d ≠ a,b,c. That gives n(n-1)(n-2)(n-3)/...
# More carefully: # pairs of triples sharing exactly 2 vertices =
# C(n,2) × C(n-2, 1) × C(n-3, 1) / 2 = ... this is getting complicated.
# Let me just derive from data.

print(f"  FORMULA DERIVATION from data:")
for i, var in enumerate(seq_varc3):
    n = i + 3
    # Var(c₃) should be a polynomial in n
    print(f"    n={n}: Var(c₃) = {var}")

# Check: is Var(c₃) = a·n⁴ + b·n³ + ... ?
# n=3: 3/4
# n=4: 3
# n=5: 285/16  ... wait that's Var(H) not Var(c₃). Let me check.
# Actually at n≤4, H = 1+2c₃, so Var(H) = 4·Var(c₃)
# n=3: Var(H) = 3/4, so Var(c₃) = 3/16? No...
# At n=3: c₃ ∈ {0, 1}. 6 transitive (c₃=0), 2 cyclic (c₃=1).
# E[c₃] = 2/8 = 1/4. E[c₃²] = 2/8 = 1/4. Var = 1/4 - 1/16 = 3/16.

print()

# ==========================================================================
# SEQUENCE F: E[total_arb] — AVERAGE TOTAL ARBORESCENCES
# ==========================================================================

print()
print("SEQUENCE F: E[total_arb] = AVERAGE TOTAL ARBORESCENCES")
print("-" * 60)
print()

# By the Matrix Tree Theorem:
# arb_r(T) = det(L_r) where L is the Laplacian.
# For a random tournament: E[L[i][i]] = E[in-degree(i)] = (n-1)/2
# E[L[i][j]] = -E[adj[j][i]] = -1/2 for i≠j
# So E[L] = (n-1)/2 × I - 1/2 × (J - I) = n/2 × I - 1/2 × J

# E[det(L_r)] is NOT simply det(E[L_r]) because det is nonlinear.
# But we can compute it directly for small n.

import numpy as np

def total_arb(adj, n):
    """Sum of arborescences over all roots."""
    total = 0
    for root in range(n):
        L = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                if i != j:
                    L[i][j] = -adj[j][i]
                    L[i][i] += adj[j][i]
        indices = [i for i in range(n) if i != root]
        minor = L[np.ix_(indices, indices)]
        total += int(round(np.linalg.det(minor)))
    return total

seq_etarb = []
for n in range(2, 7):
    t0 = time.time()
    arb_vals = []
    for bits in range(1 << comb(n, 2)):
        arb_vals.append(total_arb(tournament_adj(n, bits), n))
    total_t = len(arb_vals)
    e_arb = Fraction(sum(arb_vals), total_t)
    seq_etarb.append(e_arb)
    elapsed = time.time() - t0
    print(f"  n={n}: E[total_arb] = {e_arb} = {float(e_arb):.4f} ({elapsed:.1f}s)")

print()
print(f"  E[total_arb] sequence: {[float(s) for s in seq_etarb]}")
print()

# Is there a pattern?
# n=2: E[total_arb] = ?
# n=3: ?
# n=4: ?
# Let me see if E[total_arb] = n^{n-2} (the number of labeled trees on n vertices)
# That would be the case if all orientations of each tree are equally likely.
print(f"  COMPARISON WITH n^{{n-2}} (labeled trees):")
for i, e in enumerate(seq_etarb):
    n = i + 2
    labeled_trees = n**(n-2) if n >= 2 else 1
    print(f"    n={n}: E[total_arb] = {float(e):.2f}, n^{{n-2}} = {labeled_trees}, "
          f"ratio = {float(e)/labeled_trees:.4f}")

print()

# ==========================================================================
# SEQUENCE G: E[kings] — AVERAGE NUMBER OF KING VERTICES
# ==========================================================================

print()
print("SEQUENCE G: E[kings]")
print("-" * 60)
print()

# A king = vertex that reaches all others in ≤2 steps.
# In a random tournament: what is P(vertex v is a king)?
# v is NOT a king if ∃u such that adj[v][u]=0 AND no w with adj[v][w]=adj[w][u]=1
# This means: u beats v, AND all vertices that v beats are also beaten by u.

def count_kings(adj, n):
    kings = 0
    for v in range(n):
        is_king = True
        for u in range(n):
            if u == v: continue
            if adj[v][u]: continue
            if any(adj[v][w] and adj[w][u] for w in range(n) if w != v and w != u):
                continue
            is_king = False
            break
        if is_king: kings += 1
    return kings

seq_eking = []
for n in range(2, 7):
    t0 = time.time()
    king_vals = []
    for bits in range(1 << comb(n, 2)):
        king_vals.append(count_kings(tournament_adj(n, bits), n))
    total_t = len(king_vals)
    e_kings = Fraction(sum(king_vals), total_t)
    seq_eking.append(e_kings)
    elapsed = time.time() - t0
    print(f"  n={n}: E[kings] = {e_kings} = {float(e_kings):.4f} ({elapsed:.1f}s)")

print()
print(f"  E[kings] sequence: {[float(s) for s in seq_eking]}")

# Is E[kings] ≈ n for large n? (In random tournaments, most vertices are kings)
for i, e in enumerate(seq_eking):
    n = i + 2
    print(f"    n={n}: E[kings]/n = {float(e)/n:.4f}")

print()

# ==========================================================================
# SEQUENCE H: E[HC] — AVERAGE HAMILTONIAN CYCLES
# ==========================================================================

print()
print("SEQUENCE H: E[HC] = (n-1)!/2^n — Average Hamiltonian cycles")
print("-" * 60)
print()

# A directed Hamiltonian cycle visits all n vertices and returns.
# Fix start vertex: (n-1)! orderings of the rest.
# Each uses n arcs, each independently correct with prob 1/2.
# E[HC] = (n-1)!/2^n

seq_ehc = []
for n in range(2, 16):
    ehc = Fraction(factorial(n-1), 2**n)
    seq_ehc.append(ehc)
    print(f"  n={n:2d}: E[HC] = {ehc} = {float(ehc):.6f}")

print()
print(f"  E[HC] / E[H] = E[HC] × 2^{{n-1}} / n! = (n-1)! / (2 × n!) = 1/(2n)")
print(f"  So E[L] = E[H] - n × E[HC] = n!/2^{{n-1}} - n × (n-1)!/2^n")
print(f"          = n!/2^{{n-1}} - n!/2^n = n!/2^n = E[H]/2")
print()
print(f"  E[L] = E[H]/2 = n!/(2^n)  — EXACT for all n!")
print()

# Verify
for n in range(2, 7):
    eh = Fraction(factorial(n), 1 << (n-1))
    ehc = Fraction(factorial(n-1), 2**n)
    el = eh - n * ehc
    el_formula = Fraction(factorial(n), 2**n)
    print(f"  n={n}: E[H]={eh}, n×E[HC]={n*ehc}, E[L]={el}, formula={el_formula}, "
          f"match={'✓' if el == el_formula else '✗'}")

print()

# ==========================================================================
# SUMMARY: ALL FAST SEQUENCES
# ==========================================================================

print()
print("=" * 72)
print("  COMPLETE CATALOG OF FAST SEQUENCES (valid for ALL n)")
print("=" * 72)
print()

print(f"""
  CLOSED-FORM SEQUENCES (O(1) computation):

  A. E[c₃](n) = C(n,3)/4 = n(n-1)(n-2)/24
     Average 3-cycle count. Grows as n³/24.

  B. E[H](n) = n!/2^{{n-1}}
     Average Hamiltonian path count. Known.

  C. E[c₅](n) = C(n,5) × 3/8
     Average directed 5-cycle count.

  D. E[c_m](n) = C(n,m) × (m-1)!/2^m
     Average directed m-cycle count (any odd m).

  E. E[HC](n) = (n-1)!/2^n
     Average Hamiltonian cycle count.

  F. E[L](n) = n!/2^n = E[H]/2
     Average linear residual. EXACTLY HALF of E[H].
     This means: on average, half the Hamiltonian paths are closeable.

  G. Base(n) = 1 + n(n-1)(2n-1)/6
     H of the regular tournament (at n≤4 exact).
     = 1 + Σ k² = pyramidal number + 1.

  H. H_extreme(n) = (2^{{(n-1)/2}} + 1) × H_extreme(n-2)
     H for the T→S source-sink tournament (recursive).

  COMPUTED SEQUENCES (require enumeration at each n):

  I. Var(c₃)(n) = 3/16, ?, ... (no known closed form yet)

  J. E[total_arb](n) = {[float(s) for s in seq_etarb]}
     (no known closed form)

  K. E[kings](n) = {[float(s) for s in seq_eking]}
     (approaches n for large n?)

  THE KEY DISCOVERY: E[L] = E[H]/2 = n!/2^n

  This means: on average, EXACTLY HALF of all Hamiltonian paths
  in a random tournament can be extended to Hamiltonian cycles.
  The other half are "dead ends" — they reach the last vertex
  but can't return to the start.

  This 50/50 split is the BINARY STRUCTURE at the deepest level:
  each path either closes (probability approaching 1/2) or doesn't.
  The residual L counts the non-closeable paths.
  E[L] = E[H]/2 is the EXACT expression of this binary balance.
""")

print("opus-2026-03-22-S169")
