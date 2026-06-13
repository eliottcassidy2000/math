#!/usr/bin/env python3
"""
Algebraic investigation: WHY is β_m^{orb} = (m-3)/2?

Strategy: compute β_m from the alternating sum formula.
β_m = Budget_m - R_{m+1}. Budget_m = Σ_{d=1}^m (-1)^{m-d} Ω_d^{orb}.

We know formulas for low-degree Ω:
  Ω_0 = 1
  Ω_1 = 1
  Ω_2 = m - 1
  Ω_3 = (m-1)(2m-3)/2

Questions:
1. Is there a GENERATING FUNCTION for the orbit Omega sequence?
2. Can we compute Budget_m in closed form from the Omega generating function?
3. Is there a TOPOLOGICAL explanation (e.g., link to chessboard/matching complex)?

First approach: look for the generating function of Omega_orb.

opus-2026-03-13-S71b
"""

from fractions import Fraction
import math

# Known orbit Omega sequences
omega_orb = {
    3: [1, 1, 2, 3, 3, 2, 1],
    5: [1, 1, 4, 14, 41, 92, 140, 138, 90, 36, 6],
}

# Total orbit counts (|A_d|/m for d >= 1)
A_orb = {
    3: [1, 1, 3, 7, 13, 15, 9],
    5: [1, 1, 5, 22, 86, 286, 794, 1747, 2879, 3149, 1729],
}

# =====================================================
# Approach 1: Generating function analysis
# =====================================================

print("=== GENERATING FUNCTION ANALYSIS ===\n")

for m in [3, 5]:
    O = omega_orb[m]
    print(f"m={m}: Omega_orb = {O}")

    # EGF: Σ Ω_d * x^d / d!
    egf_coeffs = [Fraction(O[d], math.factorial(d)) for d in range(len(O))]
    print(f"  EGF coefficients: {[float(c) for c in egf_coeffs]}")

    # OGF: Σ Ω_d * x^d
    # Check if this factors nicely

    # Polynomial: P(x) = Σ Ω_d * x^d
    # Is P(1) = Σ Ω_d interesting?
    P1 = sum(O)
    print(f"  P(1) = Σ Ω_d = {P1}")

    # P(-1) = Σ (-1)^d Ω_d = chi = 1
    P_neg1 = sum((-1)**d * O[d] for d in range(len(O)))
    print(f"  P(-1) = chi = {P_neg1}")

    # Derivative P'(1)?
    P_prime_1 = sum(d * O[d] for d in range(len(O)))
    print(f"  P'(1) = Σ d*Ω_d = {P_prime_1}")

    # Budget_m = Σ_{d=1}^m (-1)^{m-d} Ω_d
    budget = sum((-1)**(m-d) * O[d] for d in range(1, m+1))
    print(f"  Budget_m = {budget}")
    print(f"  β_m = (m-3)/2 = {(m-3)//2}")

    # R_{m+1} = Budget - β_m
    R_m1 = budget - (m-3)//2
    print(f"  R_{m+1} = {R_m1}")

    print()

# =====================================================
# Approach 2: Can we express Ω_d^{orb} via known sequences?
# =====================================================

print("=== SEQUENCE ANALYSIS ===\n")

# P_7 (m=3): [1, 1, 2, 3, 3, 2, 1]
# This looks like it could be related to ballot numbers or Narayana numbers
# or the h-vector of some polytope.

# Check: is it the h-vector of a simplicial complex?
# For a Cohen-Macaulay complex, h-vector is nonneg and unimodal. ✓
# For the order polytope or something related to the Birkhoff polytope?

# Narayana numbers N(n,k):
# N(3,1)=1, N(3,2)=3, N(3,3)=1
# N(4,k): 1, 6, 6, 1
# Not matching.

# Catalan numbers: 1, 1, 2, 5, 14, 42...
# Matches Ω_0,Ω_1,Ω_2 for m=3 but Ω_3=3≠5.

# Let's check if Ω_d^{orb} satisfies a recurrence
print("m=5 orbit Omega: [1, 1, 4, 14, 41, 92, 140, 138, 90, 36, 6]")
print("Ratios:")
O = omega_orb[5]
for d in range(1, len(O)):
    if O[d-1] != 0:
        print(f"  Ω_{d}/Ω_{d-1} = {O[d]}/{O[d-1]} = {Fraction(O[d], O[d-1])} = {O[d]/O[d-1]:.4f}")

# Check second differences
print("\nDifferences:")
diff1 = [O[d+1] - O[d] for d in range(len(O)-1)]
diff2 = [diff1[d+1] - diff1[d] for d in range(len(diff1)-1)]
print(f"  First:  {diff1}")
print(f"  Second: {diff2}")

# =====================================================
# Approach 3: Factor the Omega polynomial
# =====================================================

print("\n=== POLYNOMIAL FACTORING ===\n")

# For P_7 (m=3): P(x) = 1 + x + 2x^2 + 3x^3 + 3x^4 + 2x^5 + x^6
# Let's check: P(x) = (1 + x + x^2)^2 ?
# (1 + x + x^2)^2 = 1 + 2x + 3x^2 + 2x^3 + x^4. That's degree 4, not 6.

# Try: P(x) = (1 + x)^a * (1 + x + x^2)^b
# (1+x)(1+x+x^2) = 1 + 2x + 2x^2 + x^3. Nope.
# (1+x)^2 (1+x+x^2) = 1 + 3x + 4x^2 + 3x^3 + x^4. Nope.

# Direct approach: treat P(x) for m=3 as a polynomial and find roots
import numpy as np
coeffs_m3 = [1, 1, 2, 3, 3, 2, 1]
roots_m3 = np.roots(coeffs_m3[::-1])  # numpy wants high to low
print("P_7 polynomial roots:")
for r in sorted(roots_m3, key=lambda x: (abs(x.imag), x.real)):
    if abs(r.imag) < 1e-10:
        print(f"  {r.real:.6f}")
    else:
        print(f"  {r.real:.6f} + {r.imag:.6f}i")

coeffs_m5 = [1, 1, 4, 14, 41, 92, 140, 138, 90, 36, 6]
roots_m5 = np.roots(coeffs_m5[::-1])
print("\nP_11 polynomial roots:")
for r in sorted(roots_m5, key=lambda x: (abs(x.imag), x.real)):
    if abs(r.imag) < 1e-10:
        print(f"  {r.real:.6f}")
    else:
        print(f"  {r.real:.6f} + {r.imag:.6f}i")

# =====================================================
# Approach 4: Connection to chi = p and the formula
# =====================================================

print("\n=== BUDGET FORMULA ANALYSIS ===\n")

# Budget_m = Σ_{d=1}^m (-1)^{m-d} Ω_d^{orb}
# = Ω_m - Ω_{m-1} + Ω_{m-2} - ... ± Ω_1
#
# For m=3: 3 - 2 + 1 = 2
# For m=5: 92 - 41 + 14 - 4 + 1 = 62
#
# Budget_m determines β_m: β_m = Budget - R_{m+1}
# And from top recursion: R_{m+2} = Σ_{d=m+2}^{2m} (-1)^{d-m-2} Ω_d
# R_{m+1} = Budget_{m+1} - β_{m+1} = Budget - β_m (since Budget_m = Budget_{m+1} and β_m = β_{m+1})
#
# So: β_m = Budget - (Budget - β_m) → 2β_m = Budget + Budget - 2R_{m+1}
# That's circular.
#
# The NON-CIRCULAR approach: we need the FULL Omega sequence to compute
# R_{m+1} from the top recursion.
# R_{m+1} = Budget_{m+1} - β_{m+1}
# Budget_{m+1} = Ω_{m+1} - R_{m+2}
# R_{m+2} is determined from top recursion.
#
# So: R_{m+1} = Ω_{m+1} - R_{m+2} - β_{m+1}
#   = Ω_{m+1} - R_{m+2} - β_m (since β_m = β_{m+1})
#   = Budget_{m+1} - β_m
#
# And β_m = Budget_m - R_{m+1} = Budget - (Budget - β_m) ← still circular.
#
# The resolution: R_{m+1} is determined by the ACTUAL BOUNDARY MAP, not just
# the Omega sequence. The Omega sequence determines Budget but not β_m.

# However, there may be additional constraints from representation theory.

# Key observation: β_m = m(m-3)/2 for the FULL k=0 eigenspace.
# β_m / m = (m-3)/2 = orbit Betti.
# Is (m-3)/2 the Betti number of some KNOWN complex?

# Recall: (m-3)/2 = # diagonal types of m-gon = # non-trivial gap types
# A gap type is g ∈ {2, 3, ..., (m-1)/2} (since gap g ≡ gap m-g)
# So (m-3)/2 = |(m-1)/2 - 2 + 1| = (m-1)/2 - 1 = (m-3)/2. ✓

# This is also:
# - β_1 of the cycle graph C_m (which has β_0=1, β_1=1)... no, that's 1 not (m-3)/2
# - dimension of the space of "chord diagrams" on m points... not quite
# - rank of some matrix related to the m-gon

# Let me check: is it the number of 2-element subsets of {0,...,m-1}
# minus the m "adjacent" pairs?
# C(m,2) - m = m(m-1)/2 - m = m(m-3)/2. YES! ← this is the TOTAL β_m

# And the ORBIT β_m = (m-3)/2 = (C(m,2) - m) / m = (m-1)/2 - 1

# So β_m = C(m,2) - m = edges of K_m minus edges of C_m = edges of complement of C_m
# = diagonals of convex m-gon

# The complement graph K_m \ C_m has m(m-3)/2 edges and β_1(K_m\C_m) = m(m-3)/2 - m + 1...
# no, that's the circuit rank.

# =====================================================
# Approach 5: Look at this from the FULL complex perspective
# =====================================================

print("=== FULL COMPLEX ANALYSIS ===\n")

# Full Betti: β = [1, 0, ..., 0, m(m-3)/2, m(m+1)/2, 0, ..., 0]
# Positions m and m+1 in a complex of total dimension 2m.
#
# This is the homology of a space homotopy equivalent to:
#   S^0 ∨ (∨^{m(m-3)/2} S^m) ∨ (∨^{m(m+1)/2} S^{m+1})
# (Wedge of spheres, if the space has no torsion.)
#
# Total cells: Σ |A_d| = ?
# Euler char: p = 2m+1.
#
# β_{m+1} = C(m+1,2) = triangular number.
# β_m = m(m-3)/2 = diagonals of m-gon.
#
# Remarkably: β_{m+1} = C(m+1,2) is ALSO the number of edges of K_{m+1}.
# And β_m = C(m,2) - m = edges of K_m minus edges of C_m.
#
# Sum: β_m + β_{m+1} = m(m-3)/2 + m(m+1)/2 = m(m-1) = 2·C(m,2).
# Difference: β_{m+1} - β_m = 2m = p-1.
#
# β_m = C(m,2) - m
# β_{m+1} = C(m+1,2) = C(m,2) + m
# So β_{m+1} = β_m + 2m. The "extra" 2m comes from 2m = p-1 nonzero eigenspaces.
# Each contributes 1 to β_{m+1} but 0 to β_m.

print("Pattern verification:")
for m in range(3, 20, 2):
    p = 2*m + 1
    beta_m = m*(m-3)//2
    beta_m1 = m*(m+1)//2
    chi = 1 - beta_m + beta_m1  # (-1)^m = -1 for odd m, (-1)^{m+1} = 1
    print(f"  p={p:3d} (m={m:2d}): β_m={beta_m:5d} = C({m},2)-{m}={m*(m-1)//2}-{m}, "
          f"β_{{m+1}}={beta_m1:5d} = C({m+1},2), "
          f"χ = {chi:3d} {'✓' if chi == p else '✗'}")

# =====================================================
# Approach 6: Representation-theoretic constraints
# =====================================================

print("\n=== REPRESENTATION THEORY ===\n")

# The Z_p action gives p eigenspaces. Within each:
# - Ω dims are equal (THM-125)
# - Acyclicity at d ≠ m, m+1 (verified)
# - β_0^{(k)} = δ_{k,0}
#
# From χ^{(k)}: for k=0, χ = 1. For k≠0, χ = 0.
# This gives β_m^{(k)} = β_{m+1}^{(k)} for each k.
#
# From rank shift at d=1: R_1^{(k)} = 1 for k≠0, R_1^{(0)} = 0.
# This propagates: R_d^{(k)} - R_d^{(0)} = (-1)^{d+1} for 1 ≤ d ≤ m.
# At d=m+1: shift = β_m^{(0)} - 1.
# At d ≥ m+2: shift = 0.
#
# Result: β_m^{(k)} = 0 for k≠0, β_{m+1}^{(k)} = 1 for k≠0.
# β_m and β_{m+1} for k=0 are the free parameters.
#
# From β_m^{(0)} = β_{m+1}^{(0)} = B (budget equality):
# β_m = B (entirely from k=0)
# β_{m+1} = B + (p-1) (B from k=0, 1 each from k≠0)
#
# So β_{m+1} - β_m = p - 1 = 2m. ✓
# And β_m = B remains to be determined.
#
# The ONLY unknown is B = β_m^{(0)} = m · β_m^{orb}.
# Once B is known, everything follows.

# Is B determined by the TRACE of the Z_p action on homology?
# Lefschetz number: L(g) = Σ (-1)^d tr(g* | H_d) where g is a generator of Z_p.
# For a cellular map g with no fixed cells at d > 0:
# L(g) = 1 (the single fixed vertex) + Σ_{d≥1} (-1)^d · 0 = 1.
# Wait, the Lefschetz number depends on the number of fixed CELLS.
# For Z_p acting on P_p: g permutes the p vertices cyclically.
# Fixed cells at d=0: the single vertex in the orbit complex... but in the full complex,
# there's 1 vertex (the empty partial sum) which IS fixed. So fixed cells at d=0: 1.
# For d ≥ 1: g permutes the d-cells. A d-cell is fixed iff g·σ = σ.
# Since g acts by cyclic shift on labels, a diff-seq (s_1,...,s_d) is fixed iff
# all s_i are multiplied by some constant q... but that's the Z_m action, not Z_p.
#
# Actually, Z_p acts by v ↦ v+1 (addition), not multiplication.
# Z_m = QR acts by multiplication.
# The eigenspace decomposition is under Z_p (not Z_m).
# A diff-seq σ = (s_1,...,s_d) is "shifted" by Z_p: σ ↦ σ (unchanged, since diff-seqs
# are differences, invariant under vertex translation).
# So ALL diff-seqs are fixed by Z_p! That means Z_p acts TRIVIALLY on the chain complex!
#
# Wait, that can't be right. Let me reconsider.
# The Z_p action on the Paley tournament permutes vertex labels: v ↦ v + 1 mod p.
# A Hamiltonian path visits vertices v_0, v_1, ..., v_{p-1}.
# Under Z_p, the path is relabeled. The diff-seq is based on the connection set:
# s_i = v_i - v_{i-1} mod p. Under v ↦ v+1: s_i = (v_i+1) - (v_{i-1}+1) = s_i.
# So the diff-seq IS invariant under Z_p.
#
# But the PARTIAL SUMS change: P_0 = v_0, P_1 = v_0 + s_1 = v_1, etc.
# Under Z_p: P_k ↦ P_k + 1. So partial sum 0 becomes 1 (no longer 0!).
# The "allowed" condition on diff-seqs includes: all partial sums are distinct and nonzero.
# After the shift, partial sums are all distinct and non-1 (instead of non-0).
#
# So Z_p acts on the SET of diff-seqs by changing which starting vertex is used.
# The diff-seq algebra doesn't change, but the CONSTRAINT (partial sums ≠ 0) becomes
# (partial sums ≠ 1), which changes the allowed set!
#
# In the standard definition: a diff-seq σ is in A_d iff its partial sums are
# all distinct and nonzero. Under v ↦ v+a: the constraint becomes
# partial sums ≠ a (and still distinct among themselves).
# So the "translated" diff-seq has a DIFFERENT partial sum constraint.
#
# This means Z_p does NOT act on A_d by simply relabeling!
# The eigenspace decomposition must be more subtle.

print("INSIGHT: Z_p acts on diff-seqs by changing the 'forbidden' partial sum value.")
print("The diff-seq (s_1,...,s_d) with partial sums P_k starting from vertex v_0 = 0")
print("becomes the SAME diff-seq starting from v_0 = a, with partial sums P_k + a.")
print("The 'allowed' condition changes: P_k ∉ {0, a} → this is a different complex!")
print()
print("So the correct action is on the PALEY GRAPH (the tournament itself),")
print("and the eigenspace decomposition of the BOUNDARY MAPS is what matters.")
print()

# The eigenspace decomposition is of the full boundary matrix, not just A_d.
# Z_p acts on the space of all possible chains (labeled paths), and the
# constraint matrix C_d and boundary B_d transform under this action.
# The eigenspace dims Ω^{(k)} are the dimensions of the constrained chain
# spaces in each eigenspace.

# So the "k" eigenspace corresponds to character χ_k: Z_p → C*, g ↦ ω^k.
# The chain complex in eigenspace k has the SAME combinatorial structure
# but with different phases. THM-125 shows the Ω dims are the same.

# The Lefschetz number for the generator g of Z_p:
# L(g) = Σ_d (-1)^d tr(g* | C_d)  — on the chain level
# = Σ_d (-1)^d (number of d-cells fixed by g)
#
# If g permutes d-cells freely (orbits of size p), the trace on C_d = 0 for d ≥ 1.
# The only fixed cell is the empty cell (d=0), which contributes 1.
# So L(g) = 1.
#
# On homology: L(g) = Σ_d (-1)^d tr(g* | H_d) = 1.
# tr(g* | H_d) = Σ_k β_d^{(k)} ω^k  (eigenvalue decomposition)
#
# For g generating Z_p: tr(g* | H_m) = Σ_k β_m^{(k)} ω^k
# Since β_m^{(k)} = 0 for k≠0: tr(g* | H_m) = β_m^{(0)} = m(m-3)/2.
# Similarly tr(g* | H_{m+1}) = β_{m+1}^{(0)} + Σ_{k≠0} ω^k
# = m(m-3)/2 + (-1) = m(m-3)/2 - 1 (since Σ_{k≠0} ω^k = -1).
#
# L(g) = 1 + (-1)^m · m(m-3)/2 + (-1)^{m+1} · (m(m-3)/2 - 1)
# = 1 - m(m-3)/2 + m(m-3)/2 - 1 = 0.
#
# But we said L(g) = 1! Contradiction!
#
# Wait, L(g) = 1 assumes g acts freely on cells of degree ≥ 1.
# But g acts on the PATH complex (chains of QR elements with allowed partial sums).
# The vertex (d=0) is fixed (the empty path). Are there fixed higher cells?
#
# Actually, I was confused. Let me reconsider what "cells" are.

# In GLMY homology, the chain groups C_d = Z^{A_d} where A_d = allowed diff-seqs.
# The diff-seqs themselves are sequences of QR elements.
# Z_p acts by vertex translation: if the path starts at vertex 0 and visits
# 0, s_1, s_1+s_2, ..., then translating by a gives path starting at a.
# The diff-seq (s_1,...,s_d) doesn't change, but the STARTING VERTEX does.
# The path (0, s_1, s_1+s_2, ...) becomes (a, a+s_1, a+s_1+s_2, ...).
#
# In the diff-seq formulation with partial sums starting from 0:
# The element σ = (s_1,...,s_d) ∈ A_d has partial sums 0, P_1, ..., P_d.
# The condition is: P_i ≠ P_j for all i≠j, and P_i ≠ 0 for i≥1.
# (The last condition "P_i ≠ 0" means the path doesn't revisit the start.)
#
# Under Z_p: the vertex set is relabeled v ↦ v+a. This maps the path to a
# new path with the SAME diff-seq but starting from a instead of 0.
# The new partial sums (relative to the new start) are the same: 0, P_1, ..., P_d.
# So the diff-seq is UNCHANGED, and A_d is INVARIANT under Z_p.
#
# This means Z_p acts TRIVIALLY on A_d (and hence on C_d).
# Every cell is fixed. tr(g* | C_d) = |A_d| for all g ∈ Z_p.
#
# Then L(g) = Σ_d (-1)^d |A_d| = χ = p. But we also need L(g) = Σ_d (-1)^d tr(g*|H_d).
# Since g acts trivially: tr(g*|H_d) = β_d for all g.
# So L(g) = Σ_d (-1)^d β_d = χ = p for ALL g.
#
# The Lefschetz number is p for every group element, which is consistent
# but gives NO additional information beyond χ = p.

print("Z_p acts TRIVIALLY on the diff-seq complex (diff-seqs are translation-invariant).")
print("So the 'eigenspace decomposition' must be with respect to a DIFFERENT action.")
print()
print("The action used in the eigenspace decomposition is Z_m (multiplication by QR).")
print("Z_m has m eigenspaces. This is the correct decomposition.")
print()

# So the eigenspace decomposition is by Z_m, not Z_p!
# There are m eigenspaces.
# k=0 eigenspace: dim Ω_d^{(0)} = Ω_d^{orb}
# β_d^{(0)} = β_d^{orb}
#
# Wait, but the eigenspace_complete.out says p-1 = 10 nonzero eigenspaces for P_11.
# If Z_m has m=5 eigenspaces (k=0,1,2,3,4), then k≠0 gives 4 eigenspaces.
# But the file says "10 copies" of k≠0.
#
# Let me recheck: maybe it's Z_p acting on something?
# Or Z_{p-1} = Z_{2m}?
#
# The QR subgroup has index 2 in Z_p^*. Z_p^* = Z_{p-1} = Z_{2m}.
# QR = squares mod p = the unique subgroup of index 2 in Z_p^*.
# QR has order m = (p-1)/2.
#
# The action of Z_p^* on diff-seqs: multiplication by any nonzero element.
# But non-QR multiplication REVERSES the tournament (since it maps QR to NR).
# So only QR acts as automorphisms.
#
# Z_m = QR acts on A_d. This gives m eigenspaces (characters of Z_m).
# But THM-125 says eigenspace Ω_d^{(k)} is the same for ALL k.
# And the boundary maps in each eigenspace are related by the character.
#
# The file says "10 copies" = p-1 = 2m. This suggests the eigenspace decomposition
# is over Z_{p-1}, not Z_m. Or perhaps over Z_p.
# But we showed Z_p acts trivially.
#
# Let me look at this more carefully.

# Actually, the eigenspace decomposition could be:
# Z_p acts on VERTICES, and the chain complex decomposes under this action.
# Even though the diff-seqs are invariant, the CHAIN COMPLEX (as a function
# of vertices) can decompose into eigenspaces.
#
# The chain complex C_d has basis A_d × (starting vertex).
# Wait, in GLMY homology, chains are formal sums of PATHS, not just diff-seqs.
# A path is a specific sequence of vertices, not just a diff-seq.
#
# For n vertices: the d-chains are formal sums of (d+1)-tuples (v_0,...,v_d)
# where v_i → v_{i+1} in the tournament (i.e., v_{i+1} - v_i ∈ QR).
# The allowed condition: all vertices are distinct.
#
# Z_p acts: (v_0,...,v_d) ↦ (v_0+a,...,v_d+a).
# This maps allowed paths to allowed paths (since differences are preserved).
# The action is NOT trivial on individual paths (it changes the vertices),
# even though it preserves the diff-seq.
#
# So Z_p acts freely on d-paths for d ≥ 1 (since if v_0+a = v_0 then a=0).
# Each orbit has exactly p paths (one for each starting vertex).
# The number of orbits = |A_d (paths)| / p.
#
# But wait: |A_d (paths)| = p · |A_d (diff-seqs)| for d ≥ 1.
# Because each diff-seq can start at any of p vertices, giving p distinct paths.
# (Well, we need partial sums to avoid 0 relative to start... actually the
# partial sums relative to the starting vertex must be distinct and nonzero.
# For any starting vertex v_0, the partial sums are v_0, v_0+s_1, v_0+s_1+s_2, ...
# Distinctness: v_0+P_i = v_0+P_j iff P_i = P_j (always fails since P_i ≠ P_j).
# Avoiding the start: v_0+P_k ≠ v_0 iff P_k ≠ 0 (always true for k ≥ 1).
# So ANY diff-seq from A_d gives a valid path from ANY starting vertex.)

print("CORRECTED ANALYSIS:")
print("The chain complex C_d (paths) has |C_d| = p · |A_d| for d ≥ 1.")
print("Z_p acts freely on C_d (d ≥ 1), giving p orbits = |A_d| orbits.")
print("The orbit chain complex IS the diff-seq complex A_d.")
print()
print("So the eigenspace decomposition under Z_p gives p eigenspaces:")
print("  k=0: invariant chains (sums over Z_p orbits)")
print("  k=1,...,p-1: nontrivial characters")
print()
print("Each eigenspace has dim = |A_d| (the number of diff-seq orbits).")
print("THM-125 says Ω_d is the same for each eigenspace → all p eigenspaces")
print("have Ω_d^{(k)} = Ω_d (the diff-seq Omega from orbit computation).")
print()
print("There are p eigenspaces, not m. So:")
print("  Full Ω_d (path-based) = p · Ω_d (diff-seq-based) for d ≥ 1")
print(f"  P_11: Full Ω_5 = 11 · {omega_orb[5][5]} = {11 * omega_orb[5][5]}")
print(f"  But the file says Full Ω_5 = 460")
print(f"  11 · 92 = {11*92}, 5 · 92 = {5*92}")
print()
# 5 * 92 = 460. So it's m eigenspaces, not p.
# p * 92 = 1012 ≠ 460.
# So the decomposition IS over Z_m (= QR), not Z_p.
#
# This makes sense: Z_p acts freely on paths, giving orbits = diff-seqs.
# But then the diff-seq complex is the ORBIT complex under Z_p.
# WITHIN the diff-seq complex, Z_m = QR acts (by multiplication).
# This second action decomposes the diff-seq complex into m eigenspaces.
#
# So the two-step decomposition:
# 1. Z_p on paths → orbit = diff-seqs (dimension divided by p)
# 2. Z_m on diff-seqs → eigenspaces (dimension divided by m)
#
# Total: path dim / (p · m) = path dim / (p · (p-1)/2)

print("RESOLUTION: Two-level decomposition")
print("  Level 1: Z_p orbits on paths → diff-seq complex (÷ p)")
print("  Level 2: Z_m eigenspaces on diff-seqs → orbit complex (÷ m)")
print("  Total: path complex dim = p · m · orbit complex dim")
print()

# But then where do the p-1 "nonzero eigenspaces" come from?
# The file says k≠0 has "10 copies" = p-1, not m-1.
# This suggests eigenspace decomposition under Z_p on paths.
# But Z_p gives p eigenspaces of equal dim on C_d, each of dim |A_d| = diff-seqs.
# The k=0 eigenspace is the "diff-seq" complex.
# The other p-1 eigenspaces are non-invariant.
# Each has Ω_d^{(k)} = Ω_d (diff-seq Omega).
#
# Full Betti: β_d = Σ_{k=0}^{p-1} β_d^{(k)}
# k=0 eigenspace = diff-seq complex. Its Betti is what we've been computing.
# But the diff-seq complex ITSELF has a Z_m action and can be further decomposed.
#
# The k=0 (diff-seq) complex has β_m = m(m-3)/2, β_{m+1} = m(m-3)/2.
# Wait, the file says β_5^{(0)} = 5 and β_6^{(0)} = 5.
# But the orbit complex has β_5^{orb} = 1 and β_6^{orb} = 1.
# And Ω_5^{(0)} = 460, not 92.
# Hmm, 460 = 5 · 92 = m · orbit_Ω.
# So Ω^{(0)} = m · Ω^{orb}? But Ω^{(0)} should equal Ω^{orb} for a free action!

# I think there are TWO different "k=0":
# 1. k=0 under Z_p: this is the diff-seq complex, dim = |A_d|
# 2. k=0 under Z_m (within the diff-seq complex): this is the orbit complex, dim = |A_d|/m

# The file labels the first as "k=0". So k=0 has dim |A_d|, not |A_d|/m.
# And k≠0 under Z_p has dim |A_d| each (but there are p-1 of them).
# Total: p · |A_d| = path dimension. ✓

# So the eigenspace decomposition in the file is under Z_p (not Z_m)!
# k=0 = diff-seq complex (all diff-seqs, not just orbit reps)
# k≠0 = character-twisted diff-seq complexes

# k=0 β_5 = 5 = m · β_5^{orb} = m · (m-3)/2 because the diff-seq complex
# has a Z_m action with (m-3)/2 orbit Betti, and each appears m times.

# k≠0 β_5 = 0 (no contribution to β_m from nontrivial Z_p characters)
# k≠0 β_6 = 1 (each nontrivial Z_p character contributes 1 to β_{m+1})

# So:
# Full β_m = β_m^{Z_p,k=0} + Σ β_m^{Z_p,k≠0} = m(m-3)/2 + 0 = m(m-3)/2
# Full β_{m+1} = m(m-3)/2 + (p-1)·1 = m(m-3)/2 + 2m = m(m+1)/2

print("FINAL RESOLUTION:")
print(f"  Eigenspace decomposition is under Z_p (vertex translation)")
print(f"  k=0 = diff-seq complex (dim |A_d|)")
print(f"  k=0 has its own Z_m decomposition (orbit complex, dim |A_d|/m)")
print()
print(f"  β_m comes ENTIRELY from Z_p k=0 eigenspace:")
print(f"    β_m^{{k=0}} = m · β_m^{{orb}} = m · (m-3)/2")
print(f"    β_m^{{k≠0}} = 0 (each)")
print()
print(f"  β_{{m+1}} gets contributions from ALL eigenspaces:")
print(f"    β_{{m+1}}^{{k=0}} = m · (m-3)/2")
print(f"    β_{{m+1}}^{{k≠0}} = 1 (each, p-1 = 2m eigenspaces)")
print(f"    Total: m(m-3)/2 + 2m = m(m+1)/2 = C(m+1,2)")
print()
print(f"  The single β_{{m+1}}^{{k≠0}} = 1 for each k ≠ 0 comes from")
print(f"  the rank shift theorem (HYP-710): the face-0 phase creates")
print(f"  a rank-1 defect that propagates to a single β at degree m+1.")
