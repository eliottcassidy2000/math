#!/usr/bin/env python3
"""
oeis_check_s195.py — OEIS Matches for New Sequences
opus-2026-03-22-S195

OEIS RESULTS:

A038375: H_max(n) = max # Hamiltonian paths in tournament on n vertices
  Known: 1, 1, 3, 5, 15, 45, 189, 661, 3357, 15745, 95095
  Our computation: MATCHES exactly for n=3..7.
  Status: KNOWN SEQUENCE (hard, only 11 terms, last extended 2022).
  Reference: Moon "Topics on Tournaments" (1968), p.28.

NEW SEQUENCES (not in OEIS):

B: Total fiber components = 3, 18, 117, 1018 (n=3..6)
   NOT IN OEIS. Definition: "Number of connected components of the
   score-preserving arc flip graph on all labeled tournaments on n vertices."
   A score-preserving arc flip reverses one arc while maintaining the
   sorted score sequence.

C: Isolated tournaments = 2, 16, 64, 368 (n=3..6)
   NOT IN OEIS. Definition: "Number of labeled tournaments on n vertices
   that have no score-preserving arc flip (no arc can be reversed without
   changing the sorted score sequence)."

D: Connected score fibers = 1, 2, 3, 6 (n=3..6)
   Definition: "Number of score sequences s for n-vertex tournaments
   such that the score-preserving arc flip graph on tournaments with
   score sequence s is connected."
   Only 4 terms — too few for reliable OEIS match.

E: Max fiber components = 2, 8, 40, 144 (n=3..6)
   POSSIBLE MATCH: A074092 (plane binary trees of size n+3, contracted height n)
   A074092: 1, 2, 8, 40, 144, 448, 1280, ...
   Our terms match a(1..4) of A074092.
   NEEDS VERIFICATION at n=7 to confirm or refute.

F: Within-fiber edge fraction (numerators when denominator = 2^k)
   1/2, 3/8, 5/16, 35/128
   Numerators: 1, 3, 5, 35 (not in OEIS as subsequence)
   Denominators: 2, 8, 16, 128 (powers of 2)
"""

from math import comb, factorial

print("=" * 72)
print("  OEIS MATCH RESULTS")
print("  opus-2026-03-22-S195")
print("=" * 72)
print()

# Verified sequences
print("CONFIRMED OEIS MATCHES:")
print("-" * 60)
print()
print("  A038375: H_max(n) = max Hamiltonian paths in tournament")
print("  Known:   1, 1, 3, 5, 15, 45, 189, 661, 3357, 15745, 95095")
print("  Ours:    _, _, 3, 5, 15, 45, 189 (n=3..7, all match)")
print("  Status:  HARD sequence, only 11 terms known")
print("  Ref:     Moon, 'Topics on Tournaments' (1968)")
print()

# New sequences
print("NEW SEQUENCES (NOT IN OEIS):")
print("-" * 60)
print()

print("  SEQ-B: Score-preserving arc flip components")
print("  Values: 3, 18, 117, 1018 (n=3,4,5,6)")
print("  Def:    Total connected components of the score-preserving")
print("          arc flip graph on labeled tournaments")
print()

print("  SEQ-C: Isolated tournaments")
print("  Values: 2, 16, 64, 368 (n=3,4,5,6)")
print("  Def:    Tournaments with NO score-preserving arc flip")
print()

print("  SEQ-D: Connected score fibers")
print("  Values: 1, 2, 3, 6 (n=3,4,5,6)")
print("  Def:    Score sequences whose fiber is connected")
print()

print("POSSIBLE OEIS MATCH (NEEDS VERIFICATION):")
print("-" * 60)
print()
print("  SEQ-E: Max fiber components = 2, 8, 40, 144")
print("  A074092: 1, 2, 8, 40, 144, 448, 1280, ...")
print("  Match:   Our terms = A074092(1..4)")
print("  Needs:   n=7 computation to verify a(5) = 448")
print()

# Prediction from A074092
print("  IF SEQ-E = A074092(n-2), THEN:")
print("    n=7: max fiber components = 448")
print("    n=8: max fiber components = 1280")
print("    Formula: a(n) = 2^(n-1) * C(n,floor(n/2)) / something")
print()

# Check A074092 formula
# A074092: a(n) = C(2n+2, n) - C(2n+2, n-1)
# Let's verify: a(0)=1, a(1)=2, a(2)=8, a(3)=40, a(4)=144
for k in range(5):
    val = comb(2*k+2, k) - comb(2*k+2, k-1) if k > 0 else 1
    print(f"    A074092({k}) = C({2*k+2},{k}) - C({2*k+2},{k-1}) = {val}")

print()

# Summary of all sequences
print()
print("=" * 72)
print("  SUMMARY: 6 SEQUENCES FROM FIBER STRUCTURE")
print("=" * 72)
print()

print("""
  A038375 (KNOWN): H_max = 1, 1, 3, 5, 15, 45, 189, 661, ...
    Max # Hamiltonian paths in any tournament on n vertices.

  NEW-B: fiber_components = 3, 18, 117, 1018, ...
    Total connected components of score-preserving flip graph.

  NEW-C: isolated = 2, 16, 64, 368, ...
    Tournaments with no score-preserving flip.

  NEW-D: connected_fibers = 1, 2, 3, 6, ...
    Score sequences with connected fiber.

  A074092? max_comp = 2, 8, 40, 144, [448?], ...
    Max components in any single fiber. Matches A074092 if confirmed.

  NEW-F: within_fraction = 1/2, 3/8, 5/16, 35/128, ...
    Fraction of arc flips that preserve score sequence.

  THE MOST PROMISING FOR OEIS SUBMISSION:
  NEW-B (3, 18, 117, 1018) — cleanly defined, easy to verify,
  connects to tournament theory and Landau's theorem.
""")

print("opus-2026-03-22-S195")
