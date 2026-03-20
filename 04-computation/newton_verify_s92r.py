#!/usr/bin/env python3
"""
newton_verify_s92r.py — Rigorous verification of Newton's identity claims
opus-2026-03-20-S92r

SEPARATING PROVEN from SPECULATIVE.
Every claim verified independently from definitions.
"""

from math import factorial, comb, gcd
from fractions import Fraction
import numpy as np

print("=" * 70)
print("RIGOROUS VERIFICATION OF NEWTON'S IDENTITY CLAIMS")
print("=" * 70)
print()

# =============================================================================
# CLAIM 1: e_1 = 1, e_2 = -1, e_3 = +1 for all k-nacci polynomials (k ≥ 3)
# =============================================================================

print("CLAIM 1: Elementary symmetric functions of k-nacci roots")
print("-" * 60)
print()
print("  The k-nacci polynomial: f(x) = x^k - x^{k-1} - x^{k-2} - ... - x - 1")
print()
print("  For a monic polynomial x^k + c_{k-1}x^{k-1} + ... + c_0:")
print("    e_j = (-1)^j × c_{k-j}  (Vieta's formulas)")
print()
print("  Our polynomial: c_{k-j} = -1 for j = 1, 2, ..., k.")
print("    e_1 = (-1)^1 × (-1) = +1")
print("    e_2 = (-1)^2 × (-1) = -1")
print("    e_3 = (-1)^3 × (-1) = +1  (requires k ≥ 3)")
print("    e_4 = (-1)^4 × (-1) = -1  (requires k ≥ 4)")
print("    General: e_j = (-1)^{j+1}  for j ≤ k")
print()

# Verify by computing eigenvalues of companion matrices
print("  Verification via eigenvalues:")
for k in [3, 4, 5, 6, 7]:
    # Build companion matrix
    M = np.zeros((k, k))
    for j in range(k): M[0, j] = 1
    for i in range(1, k): M[i, i-1] = 1

    evals = np.linalg.eigvals(M)
    # Compute e_1, e_2, e_3 from eigenvalues
    e1 = sum(evals).real
    e2 = sum(evals[i]*evals[j] for i in range(k) for j in range(i+1, k)).real
    e3 = 0
    if k >= 3:
        for i in range(k):
            for j in range(i+1, k):
                for l in range(j+1, k):
                    e3 += (evals[i]*evals[j]*evals[l]).real

    print(f"    k={k}: e_1 = {e1:+.8f} (expect +1), "
          f"e_2 = {e2:+.8f} (expect -1), "
          f"e_3 = {e3:+.8f} (expect +1)")

print()
print("  VERDICT: CLAIM 1 is PROVED (algebraically by Vieta + verified numerically).")
print()

# =============================================================================
# CLAIM 2: p_3 = 7 for all k-nacci (k ≥ 3) via Newton's identity
# =============================================================================

print("CLAIM 2: p_3 = 7 universally for k ≥ 3")
print("-" * 60)
print()
print("  Newton's identity for p_3 (with n = 3 ≤ k):")
print("    p_3 = e_1·p_2 - e_2·p_1 + 3·e_3")
print()
print("  Prerequisites:")
print("    p_1 = e_1 = 1")
print("    p_2 = e_1·p_1 - 2·e_2 = 1×1 - 2×(-1) = 1 + 2 = 3")
print("    p_3 = e_1·p_2 - e_2·p_1 + 3·e_3 = 1×3 - (-1)×1 + 3×1 = 3 + 1 + 3 = 7")
print()
print("  Each step uses only e_1, e_2, e_3 = (1, -1, 1).")
print("  These are FIXED for all k ≥ 3 (Claim 1).")
print("  Therefore p_3 = 7 for ALL k-nacci with k ≥ 3.")
print()

# Verify directly
print("  Direct verification via Tr(M_k^3):")
for k in range(2, 12):
    M = np.zeros((k, k))
    for j in range(k): M[0, j] = 1
    for i in range(1, k): M[i, i-1] = 1
    M3 = M @ M @ M
    tr = round(np.trace(M3).real)
    status = "✓" if (k >= 3 and tr == 7) else ("(k<3)" if k < 3 else "FAIL!")
    print(f"    k={k:2d}: Tr(M^3) = {tr:4d}  {status}")

print()
print("  VERDICT: CLAIM 2 is PROVED.")
print("  Proof: algebraic (Newton + Vieta) + verified computationally k=2..11.")
print("  Note: for k=2 (Fibonacci), Tr(M^3) = 4, NOT 7. Claim holds for k ≥ 3 only.")
print()

# =============================================================================
# CLAIM 3: p_5 = 21 for k = 3 (tribonacci) specifically
# =============================================================================

print("CLAIM 3: p_5 = 21 for k = 3 (tribonacci)")
print("-" * 60)
print()
print("  For k = 3: the recurrence for n > 3 is:")
print("    p_n = e_1·p_{n-1} - e_2·p_{n-2} + e_3·p_{n-3}")
print("         = p_{n-1} + p_{n-2} + p_{n-3}")
print()
print("  p_4 = p_3 + p_2 + p_1 = 7 + 3 + 1 = 11")
print("  p_5 = p_4 + p_3 + p_2 = 11 + 7 + 3 = 21")
print()

# Verify
M3 = np.array([[1,1,1],[1,0,0],[0,1,0]], dtype=float)
M3_5 = np.linalg.matrix_power(M3, 5)
tr5 = round(np.trace(M3_5).real)
print(f"  Direct: Tr(M_3^5) = {tr5}")
print(f"  Match: {tr5 == 21}")
print()

# For other k values
print("  p_5 for other k-nacci values:")
for k in range(2, 10):
    M = np.zeros((k, k))
    for j in range(k): M[0, j] = 1
    for i in range(1, k): M[i, i-1] = 1
    M5 = np.linalg.matrix_power(M, 5)
    tr = round(np.trace(M5).real)
    is21 = tr == 21
    print(f"    k={k}: p_5 = {tr:4d}  {'= 21 ✓' if is21 else ''}")

print()
print("  VERDICT: CLAIM 3 is PROVED.")
print("  p_5 = 21 ONLY for k = 3. For k ≥ 4, p_5 > 21.")
print()

# =============================================================================
# CLAIM 4: p_5 = p_3 × p_2 for k = 3
# =============================================================================

print("CLAIM 4: p_5 = p_3 × p_2 for k = 3")
print("-" * 60)
print()
print("  p_3 = 7, p_2 = 3, p_5 = 21.")
print("  7 × 3 = 21 = p_5. ✓")
print()
print("  ALGEBRAIC PROOF:")
print("    p_5 = 2p_3 + 2p_2 + p_1  (from expanding the recurrence)")
print("    p_3 × p_2 = p_3 × p_2")
print("    These are equal iff 2p_3 + 2p_2 + p_1 = p_3 × p_2")
print("    Substituting p_1=1: 2p_3 + 2p_2 + 1 = p_3 × p_2")
print("    Substituting p_3 = p_2 + 4 (from p_3 = p_2 + p_1 + 3e_3 = 3+1+3 = p_2+4):")
print("    2(p_2+4) + 2p_2 + 1 = (p_2+4) × p_2")
print("    4p_2 + 9 = p_2² + 4p_2")
print("    9 = p_2²")
print("    p_2 = 3. ✓ (This is ALWAYS true.)")
print()
print("  So the identity p_5 = p_3 × p_2 is EQUIVALENT to p_2 = 3.")
print("  p_2 = 3 holds for ALL k-nacci (k ≥ 2) since p_2 = e_1²-2e_2 = 1+2 = 3.")
print("  But p_5 = 2p_3 + 2p_2 + p_1 is specific to the k=3 recurrence.")
print("  For k ≥ 4: p_5 includes additional terms, breaking the identity.")
print()
print("  VERDICT: CLAIM 4 is PROVED.")
print("  It is a consequence of p_2 = 3 (universal) combined with")
print("  the tribonacci recurrence length (specific to k=3).")
print()

# =============================================================================
# CLAIM 5: The forbidden values 7 and 21 ARE p_3 and p_5 of tribonacci
# =============================================================================

print("CLAIM 5: Forbidden H-values = tribonacci power sums at positions 3, 5")
print("-" * 60)
print()
print("  VERIFIED FACTS:")
print("    Tr(M_3^3) = 7 = the first forbidden H-value. ✓")
print("    Tr(M_3^5) = 21 = the second forbidden H-value. ✓")
print()
print("  STATUS: This is a VERIFIED CORRELATION, not a proved causal link.")
print()
print("  The claims that ARE proved:")
print("    - H = 7 is forbidden for all tournaments (THM-029, via K_3 obstruction)")
print("    - H = 21 is forbidden for all tournaments (THM-079/115)")
print("    - p_3 = Tr(M_3^3) = 7 (Newton's identity, proved above)")
print("    - p_5 = Tr(M_3^5) = 21 (recurrence computation, proved above)")
print()
print("  What is NOT proved:")
print("    - That p_3 = 7 is the REASON H = 7 is forbidden")
print("    - That the tribonacci polynomial is the MECHANISM of forbiddenness")
print("    - That there's a causal arrow from Newton's identity to the K_3 obstruction")
print()
print("  The connection might be:")
print("    (a) A genuine deep relationship (the tribonacci encodes the cycle structure)")
print("    (b) A numerical coincidence (7 is small and appears in many contexts)")
print("    (c) An indirect consequence of both being controlled by the same (e_1,e_2,e_3)")
print()

# Test option (b): how many other sequences also hit 7 at position 3?
print("  Testing option (b): Is it remarkable that p_3 = 7?")
print("  Other degree-3 polynomials with integer coefficients:")
print()

count_p3_equals_7 = 0
count_total = 0
for a in range(-3, 4):
    for b in range(-3, 4):
        for c in range(-3, 4):
            # Polynomial x^3 + a*x^2 + b*x + c
            e1 = -a
            e2 = b
            e3 = -c
            p1 = e1
            p2 = e1*p1 - 2*e2
            p3 = e1*p2 - e2*p1 + 3*e3
            count_total += 1
            if p3 == 7:
                count_p3_equals_7 += 1

print(f"  Among {count_total} degree-3 polynomials with |coefficients| ≤ 3:")
print(f"  {count_p3_equals_7} have p_3 = 7 ({count_p3_equals_7/count_total*100:.1f}%)")
print(f"  So p_3 = 7 is NOT extremely rare but also not generic.")
print()

# =============================================================================
# CLAIM 6: The three terms of p_3 = 1 + 3 + 3 have meaning
# =============================================================================

print("CLAIM 6: Decomposition p_3 = e_1^3 - 3e_1e_2 + 3e_3 = 1 + 3 + 3")
print("-" * 60)
print()
print("  p_3 = e_1^3 - 3·e_1·e_2 + 3·e_3")
print("       = 1^3 - 3·1·(-1) + 3·1")
print("       = 1 + 3 + 3")
print()
print("  The three terms:")
print("    Term 1: e_1^3 = 1   (cube of the sum of roots)")
print("    Term 2: -3e_1e_2 = 3 (correction from pairwise products)")
print("    Term 3: 3e_3 = 3    (contribution from triple product)")
print()
print("  STATUS: The decomposition is ALGEBRAICALLY CORRECT.")
print("  The 'virtual 2-cycle' interpretation is SPECULATIVE.")
print()
print("  The terms are just the standard Newton's identity coefficients.")
print("  Calling term 2 a 'virtual 2-cycle bonus' is poetic but not rigorous.")
print("  What IS rigorous: e_2 = -1 encodes the coefficient of x^{k-2} in the")
print("  k-nacci polynomial, which corresponds to the sum of products of")
print("  pairs of roots. For tournament theory, this relates to arc-pair structure.")
print()
print("  VERDICT: The decomposition is CORRECT. The interpretation is SPECULATIVE.")
print()

# =============================================================================
# CLAIM 7: Are p_7 = 71, p_9 = 241 achievable as H-values?
# =============================================================================

print("CLAIM 7: p_7 = 71 and p_9 = 241 should be achievable H-values")
print("-" * 60)
print()

# p values
p = [0]*20
p[1], p[2], p[3] = 1, 3, 7
for n in range(4, 20):
    p[n] = p[n-1] + p[n-2] + p[n-3]

print(f"  Tribonacci power sums: {[p[n] for n in range(1, 12)]}")
print()
print("  If ONLY p_3 and p_5 are permanently forbidden (HYP-1352),")
print("  then p_7 = 71 and p_9 = 241 must be achievable at some n.")
print()
print("  p_7 = 71:  achievable at n = 7? (max H at n=7 is 189 > 71. Possible.)")
print("  p_9 = 241: achievable at n = 8? (max H at n=8 is ~6880 >> 241. Possible.)")
print()
print("  STATUS: Not verified computationally (would need n=7 exhaustive check).")
print("  HYP-1352 says only {7,21} are permanent gaps, verified through n=7 sampling.")
print("  So p_7 = 71 is EXPECTED to be achievable but NOT YET CONFIRMED.")
print()

# =============================================================================
# SUMMARY TABLE
# =============================================================================

print("=" * 70)
print("SUMMARY: WHAT IS PROVED vs SPECULATIVE")
print("=" * 70)
print()
print(f"{'Claim':>5} | {'Statement':>50} | {'Status':>12}")
print("-" * 75)
claims = [
    ("1", "e_1=1, e_2=-1, e_3=+1 for k-nacci k≥3", "PROVED"),
    ("2", "p_3 = 7 for all k-nacci k≥3 (Newton)", "PROVED"),
    ("3", "p_5 = 21 for k=3 (tribonacci) only", "PROVED"),
    ("4", "p_5 = p_3 × p_2 = 21 for k=3 only", "PROVED"),
    ("5", "Forbidden {7,21} = tribonacci {p_3,p_5}", "CORRELATION"),
    ("6", "p_3 = 1+3+3 = identity+2cyc+3cyc", "SPECULATIVE"),
    ("7", "p_7=71 is achievable as H-value", "EXPECTED"),
]
for num, stmt, status in claims:
    print(f"{num:>5} | {stmt:>50} | {status:>12}")

print(f"""

DEPENDENCIES:
  Claim 2 depends on Claim 1 (needs e_1, e_2, e_3 values).
  Claim 3 depends on Claim 2 (needs p_3 = 7 as initial condition).
  Claim 4 depends on Claims 2+3 (algebraic identity).
  Claim 5 depends on Claims 2+3 AND the forbidden value theorems (THM-029, THM-079).
  Claim 6 is an interpretation of the algebraic decomposition in Claim 2.
  Claim 7 depends on HYP-1352 (only 7 and 21 forbidden).

WHAT WOULD MAKE CLAIM 5 A THEOREM:
  Need to show that the K_3 obstruction (which proves H≠7) is
  ALGEBRAICALLY EQUIVALENT to the Newton identity p_3 = 7.
  Specifically: need to prove that H(T) = p_3 requires exactly 3
  mutually overlapping cycles, and that the number of such cycles
  is controlled by the elementary symmetric functions of the
  tribonacci polynomial.

  This would require:
  (a) A connection between the transfer matrix M_3 and the conflict graph Omega
  (b) Showing that Tr(M_3^3) counts something related to K_3 components
  (c) Proving that the K_3 impossibility follows from the trace identity

  NONE OF THESE are currently proved. The connection is observational.
""")

print("opus-2026-03-20-S92r: RIGOROUS VERIFICATION")
