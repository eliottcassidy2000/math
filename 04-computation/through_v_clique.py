#!/usr/bin/env python3
"""
THEOREM: Through-v cycles ALWAYS form a clique in Omega(T).

Proof: Any two cycles C1, C2 through v share at least vertex v.
Therefore they are adjacent in Omega(T) (adjacency = sharing a vertex).
QED.

This immediately explains why higher-order IE terms vanish:
- The inclusion-exclusion over through-v cycles has terms for
  subsets {C1,...,Cr} of through-v cycles that are MUTUALLY non-adjacent.
- But since through-v cycles form a clique, no pair is non-adjacent.
- Therefore ALL terms with r >= 2 are zero.
- The IE collapses to just the r=1 term: sum_C f(C).

Combined with f(C) = 2*mu(C) (THM-069), this gives:
  H(T) - H(T-v) = sum_C f(C) = 2 * sum_C mu(C)
which is exactly Claim A!

Wait - this means we have an independent proof of Claim A from OCF!
Let me verify this reasoning is correct...

Actually, let me re-examine. The IE formula is:
  sum_{S indep, S ni some C ni v} 2^|S|
  = sum_r (-1)^{r+1} sum_{C1<...<Cr through v, mutually non-adj}
    2^r * I(Omega - N[C1] - ... - N[Cr], 2)

But ALL through-v cycles ARE adjacent (share v). So the "mutually non-adjacent"
condition with r >= 2 is NEVER satisfied. The IE collapses to r=1:
  = sum_{C through v} f(C) = sum_C 2*mu(C)

THIS IS THE PROOF!

opus-2026-03-07-S36
"""

print("=" * 60)
print("THEOREM: Through-v cycles form a clique in Omega(T)")
print("=" * 60)
print()
print("Proof:")
print("  Let C1, C2 be any two directed odd cycles of T, both containing v.")
print("  Then v in V(C1) ∩ V(C2), so V(C1) ∩ V(C2) ≠ ∅.")
print("  By definition of Omega(T), C1 and C2 are adjacent.")
print("  Since this holds for any pair, through-v cycles form a clique.")
print("  QED.")
print()
print("=" * 60)
print("COROLLARY: Higher-order IE terms vanish")
print("=" * 60)
print()
print("The inclusion-exclusion formula for")
print("  sum_{S indep, S ∋ some C ∋ v} 2^|S|")
print("expands as:")
print("  sum_{r=1}^{nv} (-1)^{r+1} sum_{C1<...<Cr through v, mut. non-adj}")
print("    2^r * I(Omega \\ (N[C1] ∪ ... ∪ N[Cr]), 2)")
print()
print("Since through-v cycles form a clique, NO pair is non-adjacent.")
print("Therefore all terms with r ≥ 2 have an EMPTY sum.")
print("The formula collapses to r=1:")
print("  = sum_{C ∋ v} 2 * I(Omega \\ N[C], 2)")
print("  = sum_{C ∋ v} f(C)")
print("  = sum_{C ∋ v} 2 * mu(C)   [by THM-069: f(C) = 2*mu(C)]")
print()
print("=" * 60)
print("CONSEQUENCE: Clean proof of Claim A from OCF")
print("=" * 60)
print()
print("Given: H(T) = I(Omega(T), 2) [Grinberg-Stanley]")
print("       I(Omega, 2) - I(Omega(T-v), 2) = sum_{S indep, S ∋ some C ∋ v} 2^|S|")
print("         [Claim B, proved]")
print()
print("Step 1: Through-v cycles form a clique [PROVED above]")
print("Step 2: IE collapses: RHS = sum_C f(C) [from Step 1]")
print("Step 3: Graph equality: f(C) = 2*mu(C) [THM-069]")
print("Step 4: Combining: H(T) - H(T-v) = 2 * sum_{C∋v} mu(C)")
print()
print("This is a self-contained proof that H(T) = I(Omega(T), 2)")
print("implies Claim A, without ANY additional machinery!")
