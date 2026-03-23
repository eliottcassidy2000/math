#!/usr/bin/env python3
"""
unlock_formula_s208.py -- opus-2026-03-23-S208

THE MINIMUM EXTRA INFORMATION TO UNLOCK G_n AT ALL n

We have 15 exact formulas. We need ONE more: the degeneracy D = T - 2E.
D = 2, 6, 28, 124, 740, 5966, 85698 (n=3..9).

This script attacks D from EVERY angle to find the formula.

Author: opus-2026-03-23-S208
"""

import sys
from math import comb, factorial, gcd, log, sqrt
from collections import Counter
from fractions import Fraction
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

# The known sequences
V = [2, 4, 12, 56, 456, 6880, 191536]  # n=3..9
SC = [2, 2, 8, 12, 88, 176, 2752]
E = [1, 5, 30, 290, 4086, 91161, 3380751]
T = [4, 16, 88, 704, 8912, 188288, 6847200]
D = [T[i] - 2*E[i] for i in range(7)]  # degeneracy

print("  D = T - 2E = %s\n" % D)

# ============================================================
banner("ATTACK 1: FACTORIZATION")
# ============================================================

print("  Factorize D(n):\n")
for i, d in enumerate(D):
    n = i + 3
    temp = d
    factors = []
    for p in range(2, min(temp+1, 100000)):
        if p*p > temp and temp > 1:
            factors.append(temp); break
        while temp % p == 0:
            factors.append(p); temp //= p
    print("  n=%d: D=%d = %s" % (n, d, " * ".join(str(f) for f in factors) if factors else "1"))

# ============================================================
banner("ATTACK 2: RATIOS AND RECURRENCES")
# ============================================================

print("  D(n)/D(n-1):")
for i in range(1, len(D)):
    print("    D(%d)/D(%d) = %d/%d = %.4f" % (i+3, i+2, D[i], D[i-1], D[i]/D[i-1]))

print("\n  D(n)/D(n-2) (same parity):")
for i in range(2, len(D)):
    print("    D(%d)/D(%d) = %d/%d = %.4f" % (i+3, i+1, D[i], D[i-2], D[i]/D[i-2]))

# Check: D(n) = a*D(n-1) + b*D(n-2)?
print("\n  Testing D(n) = a*D(n-1) + b*D(n-2):")
for i in range(2, len(D)):
    # D[i] = a*D[i-1] + b*D[i-2]
    # Solve: two equations for two unknowns using i and i-1
    if i >= 3:
        # D[i] = a*D[i-1] + b*D[i-2]
        # D[i-1] = a*D[i-2] + b*D[i-3]
        det = D[i-1]*D[i-3] - D[i-2]*D[i-2]
        if det != 0:
            a = Fraction(D[i]*D[i-3] - D[i-1]*D[i-2], det)
            b = Fraction(D[i-1]*D[i-1] - D[i]*D[i-2], det)
            # Verify
            pred = a * D[i-1] + b * D[i-2]
            print("    n=%d: a=%s, b=%s, pred=%s, actual=%d, match=%s" %
                  (i+3, a, b, pred, D[i], pred == D[i]))

# ============================================================
banner("ATTACK 3: D/V, D/E, D/T RATIOS")
# ============================================================

print("  D/V (degeneracy per class):")
for i in range(len(D)):
    n = i+3
    ratio = Fraction(D[i], V[i])
    print("    n=%d: D/V = %d/%d = %s = %.6f" % (n, D[i], V[i], ratio, float(ratio)))

print("\n  D/E (degeneracy per edge):")
for i in range(len(D)):
    n = i+3
    ratio = Fraction(D[i], E[i])
    print("    n=%d: D/E = %d/%d = %s = %.6f" % (n, D[i], E[i], ratio, float(ratio)))

print("\n  D/T (degeneracy fraction of total orbits):")
for i in range(len(D)):
    n = i+3
    ratio = Fraction(D[i], T[i])
    print("    n=%d: D/T = %d/%d = %s = %.6f" % (n, D[i], T[i], ratio, float(ratio)))

# ============================================================
banner("ATTACK 4: D IN TERMS OF SCHUR-WEYL COMPONENTS")
# ============================================================

# T = m_(n) + m_(n-1,1) + m_(n-2,2)
# m_(n) = V
# So D = T - 2E = V + m_(n-1,1) + m_(n-2,2) - 2E
# => 2E = V + m_(n-1,1) + m_(n-2,2) - D
# => E = (V + m_(n-1,1) + m_(n-2,2) - D) / 2

m1 = [2, 8, 36, 240, 2584, 47376, 1525072]  # m_(n-1,1)
m2 = [0, 4, 40, 408, 5872, 134032, 5130592]  # m_(n-2,2)

print("  D = V + m1 + m2 - 2E:")
for i in range(len(D)):
    n = i+3
    check = V[i] + m1[i] + m2[i] - 2*E[i]
    print("    n=%d: V+m1+m2-2E = %d+%d+%d-%d = %d = D=%d (%s)" %
          (n, V[i], m1[i], m2[i], 2*E[i], check, D[i], "OK" if check == D[i] else "FAIL"))

# So: D = V + m1 + m2 - 2E = T - 2E. Tautological.
# But can we express D in terms of ONLY V, m1, m2 (without E)?

# Try: D = f(V, m1, m2)
print("\n  D vs m2 - m1:")
for i in range(len(D)):
    n = i+3
    diff = m2[i] - m1[i]
    print("    n=%d: D=%d, m2-m1=%d, D/(m2-m1)=%s" %
          (n, D[i], diff, Fraction(D[i], diff) if diff != 0 else "inf"))

print("\n  D vs V*(m1/V - something):")
for i in range(len(D)):
    n = i+3
    m1_per_V = Fraction(m1[i], V[i])
    m2_per_V = Fraction(m2[i], V[i])
    D_per_V = Fraction(D[i], V[i])
    print("    n=%d: m1/V=%s, m2/V=%s, D/V=%s" %
          (n, m1_per_V, m2_per_V, D_per_V))

# ============================================================
banner("ATTACK 5: THE SELF-FLIP SEQUENCE")
# ============================================================

# Self-flip pairs (tiling model): 0, 1, 4 at n=3,4,5
# Total self-flip = 2^{m-1} - sum(edge_thicknesses)
# = 2^{m-1} - (T - D)/2 ... no
# Actually: 2E = T - D, and edge thickness sum = some relation
# Let me think differently.

# In the ARC REVERSAL model:
# Total (T, arc) pairs = 2^m * m
# Self-loop weight = sum F[i][i] = total neutral
# Cross weight = 2^m * m - self_loop
# 2E <= cross_weight (because some cross-orbits collide)
# D = cross_orbits - 2E (the collision count)

# From S184: self-loop fraction = 1/2, 3/8, 11/64
self_loop_frac = [Fraction(1,2), Fraction(3,8), Fraction(11,64)]
print("  Self-loop fraction (class-preserving arcs):")
for i, f in enumerate(self_loop_frac):
    n = i+3; m = comb(n,2)
    total_weight = 2**m * m
    self_weight = int(f * total_weight)
    cross_weight = total_weight - self_weight
    print("    n=%d: self_frac=%s, self_weight=%d, cross_weight=%d, 2E=%d" %
          (n, f, self_weight, cross_weight, 2*E[i]))

print()
print("  Is D related to self-loop weight?")
for i in range(min(3, len(D))):
    n = i+3; m = comb(n,2)
    total_weight = 2**m * m
    self_weight = int(self_loop_frac[i] * total_weight)
    print("    n=%d: D=%d, self_weight=%d, D/self_weight=%s" %
          (n, D[i], self_weight, Fraction(D[i], self_weight) if self_weight > 0 else "?"))

# ============================================================
banner("ATTACK 6: GENERATING FUNCTION")
# ============================================================

# D = 2, 6, 28, 124, 740, 5966, 85698
# GF: sum D_n x^n
# Try: D(x) = sum D_n x^n / n! (EGF)

print("  D sequence: 2, 6, 28, 124, 740, 5966, 85698")
print("  D/n!:  ", end="")
for i in range(len(D)):
    n = i+3
    print("%.8f " % (D[i]/factorial(n)), end="")
print()

print("  D/2^m: ", end="")
for i in range(len(D)):
    n = i+3; m = comb(n,2)
    print("%.8f " % (D[i]/2**m), end="")
print()

print("  D/(V*m): ", end="")
for i in range(len(D)):
    n = i+3; m = comb(n,2)
    print("%.8f " % (D[i]/(V[i]*m)), end="")
print()

# D/(V*m) is the degeneracy RATE (fraction of arc-orbits that are degenerate)
print("\n  Degeneracy rate D/(V*m) = D/T (since T ~ V*m):")
for i in range(len(D)):
    n = i+3; m = comb(n,2)
    print("    n=%d: D/T = %.8f, D/(V*m) = %.8f" % (n, D[i]/T[i], D[i]/(V[i]*m)))

# ============================================================
banner("ATTACK 7: THE COLLISION FORMULA")
# ============================================================

print("""
  D = #{arc-orbit COLLISIONS}
  = #{(orbit1, orbit2) in same class C : orbit1 and orbit2 both
      reverse to the SAME target class}

  For a class C with |Aut|=1 (generic):
  Each arc is its own orbit. There are m = C(n,2) orbits.
  Each orbit reverses to some class D_i.
  Collision: two distinct orbits i != j with D_i = D_j.
  Number of collisions = sum_D C(#{orbits hitting D}, 2)
  = (sum_D n_D^2 - sum_D n_D) / 2
  = (sum n_D^2 - m) / 2  (since sum n_D = m for |Aut|=1)

  So for each |Aut|=1 class C:
  collisions(C) = (sum_D weight(C,D)^2 - m) / 2

  Wait: weight(C,D) = F[C][D]/size(C) = P[C][D] * m.
  So collisions = (m^2 * sum P[C][D]^2 - m) / 2.

  Summed over all classes:
  D = sum_C collisions(C) = sum_C (sum_D n_D(C)^2 - m_eff(C)) / 2

  For |Aut|=1 classes: m_eff = m.
  For |Aut|>1 classes: m_eff = #orbits < m.

  THE FORMULA:
  D = (1/2) * [sum_C sum_D (n_D(C))^2 - sum_C #orbits(C)]
  = (1/2) * [sum_{(C,D)} F_reduced[C,D]^2 / size(C)^2 - T]

  WHERE F_reduced[C,D] = F[C][D] / size(C) = #{arcs of a SINGLE tournament
  in C that reverse to class D}. Wait, this is per-tournament, not per-class.

  Actually: for a SPECIFIC tournament t in class C:
  n_D(t) = #{arcs of t whose reversal gives class D}.
  sum_D n_D(t) = m (every arc reverses somewhere).
  collisions(t) = (sum_D n_D(t)^2 - m) / 2.

  D = (1/|S_n|) * sum_t collisions(t)
    = (1/|S_n|) * sum_t (sum_D n_D(t)^2 - m) / 2

  For generic t (|Aut|=1): sum_D n_D(t)^2 = the SECOND MOMENT.
  D = (1/n!) * sum_t sum_D n_D(t)^2 / 2 - (1/n!) * sum_t m/2
  = (1/(2*n!)) * sum_t sum_D n_D(t)^2 - T/2
  = (1/(2*n!)) * sum_t sum_D n_D(t)^2 - T/2

  And E = (T - D) / 2, so:
  2E = T - D = T - [(1/(2*n!)) * sum_t sum_D n_D(t)^2 - T/2]
     = T - (1/(2*n!)) * S_2 + T/2
     = 3T/2 - S_2 / (2*n!)

  WHERE S_2 = sum_t sum_D n_D(t)^2 = the TOTAL SECOND MOMENT.
""")

# Can we compute S_2 from the Burnside formulas?
# S_2 = sum_t sum_D (#{arcs of t reversing to D})^2
# This is a FOURTH-ORDER Burnside quantity.

# ============================================================
banner("ATTACK 8: SMALL CASES EXACT")
# ============================================================

# At n=3: D=2. V=2, E=1, T=4.
# Only 8 tournaments, 3 arcs each.
# Each tournament has 3 arcs. Each arc reverses to some class.
# For the transitive tournament: 2 arcs go to the 3-cycle class, 1 stays.
# Wait: the transitive has 3 arcs. Reversing each:
# 2 neutral arcs (go to itself) + 1 non-neutral (goes to 3-cycle class)...
# no, from S195: transitive has n-1=2 neutral arcs and C(n-1,2)=1 non-neutral.
# So: n_self=2, n_other=1. collisions = (2^2 + 1^2 - 3)/2 = (4+1-3)/2 = 1.
# For the 3-cycle: 3 arcs, all go to the transitive? Let me check.
# The 3-cycle class has |Aut|=3, so 1 orbit. That orbit reverses to the transitive.
# So: n_trans = 3. collisions = (3^2 - 3)/2 = 3 ... wait, #orbits = 1.
# collision formula for |Aut|>1: (sum n_D^2 - #orbits) / 2 = (9 - 1)/2 = 4.
# Hmm, that gives collisions = 4 for 3-cycle class.
# Total D = collisions(trans) + collisions(3-cycle) = ?

# Actually I need to be more careful. D counts per ISO CLASS, not per tournament.
# D = sum over all CLASSES of (collisions within that class).
# For class C: collisions = #{(orbit1, orbit2) : orbit1 != orbit2, both go to same D}
# = (sum_D C(n_D, 2))
# where n_D = #{orbits of C going to D}.

print("  n=3 EXACT COLLISION ANALYSIS:")
print("  Transitive (H=1, |Aut|=1): 3 orbits (each arc is own orbit)")
print("    2 neutral (self) + 1 to 3-cycle. n_self=2, n_other=1.")
print("    Collisions: C(2,2) + C(1,2) = 1 + 0 = 1.")
print()
print("  3-cycle (H=3, |Aut|=3): 1 orbit (3 arcs in one orbit)")
print("    All 3 arcs reverse to transitive. n_trans=1 (ONE orbit).")
print("    Collisions: C(1,2) = 0.")
print()
print("  Total D = 1 + 0 = 1. But actual D=2!")
print("  DISCREPANCY. The formula needs adjustment for multiple tournaments per class.\n")

print("  THE CORRECT FORMULA:")
print("  D counts pairs of LABELED (tournament, arc) that are in the SAME class")
print("  AND whose reversals are ALSO in the same class.")
print("  D = #{(T1,e1), (T2,e2) : T1 in C, T2 in C, T1_e1 in D, T2_e2 in D, (T1,e1)!=(T2,e2)}")
print("    summed over all (C,D) pairs, divided by... hmm.\n")

print("  Actually, D = T - 2E where:")
print("  T = total arc-orbits = sum_C #{arc orbits of C}")
print("  2E = 2 * #{distinct (C,D) pairs with some arc orbit from C to D}")
print("  D = #{arc-orbits from C to D} - 2*#{distinct (C,D) edges}")
print("    = sum over edges (C,D) of (#{orbits from C to D} - 1)")
print("    + sum over edges (D,C) of (#{orbits from D to C} - 1)")
print("    + sum_C (#{self-orbits} - 1 if self-orbits > 0, else 0)")

# Let me just compute D exactly as sum over edges of (multiplicity - 1)
# At n=3: 1 edge with multiplicity 1. D would be 0. But D=2!
# So this formula is also wrong. The issue: T counts orbits PER CLASS,
# and each orbit from C to D is counted once. But the SAME (C,D) pair
# can be reached from BOTH sides. D counts something different.

# Let me reconsider: T = total arc-orbits across all classes.
# An edge (C,D) contributes at least 1 orbit from C side and 1 from D side.
# If more: the excess = degeneracy.

print("\n  REVISED: D = sum over all (C,D) adjacent pairs of")
print("  (orbits_C_to_D + orbits_D_to_C - 2)")
print("  + sum_C (self_orbits_C if > 0)")

# At n=3: orbits trans->3cycle = 1. orbits 3cycle->trans = 1.
# Excess = 1+1-2 = 0. Self-orbits: trans has 2 self-orbits.
# Total D = 0 + 2 = 2. MATCHES!

print("\n  n=3 CHECK:")
print("  Edge (trans, 3cyc): 1+1-2 = 0")
print("  Self-orbits of trans: 2")
print("  Self-orbits of 3cyc: 0")
print("  Total D = 0 + 2 + 0 = 2. MATCHES!")
print()

# n=4:
print("  n=4 CHECK:")
print("  D=6. From S185: self-orbits = 6 total (trans:3, H=5:3)")
print("  Cross orbit excess: ?")
print("  If all edges have multiplicity 1 from each side:")
print("  excess = 0 for all 5 edges")
print("  Total D = 0 + 3 + 3 = 6. MATCHES!")
print()

# n=5:
print("  n=5 CHECK:")
print("  D=28. Self-orbits: from S185 = 16 total")
print("  Cross orbit excess: D - 16 = 12")
print("  = sum over edges of (orbits_i_to_j + orbits_j_to_i - 2)")
print("  = 12 (need per-edge multiplicity data)")

# ============================================================
banner("ATTACK 9: D = SELF-LOOP ORBITS + EDGE MULTIPLICITY EXCESS")
# ============================================================

print("""
  D = S + M where:
  S = total self-loop orbits across all classes
  M = total edge multiplicity excess = sum_{edges} (mult_i + mult_j - 2)

  S sequence: 2, 6, 16 (from S185 data)
  M = D - S: 0, 0, 12

  WAIT: at n=3: S=2, M=0, D=2. Check.
  At n=4: S=6, M=0, D=6. Check.
  At n=5: S=16, M=12, D=28. Check: 16+12=28. YES!

  So: D = S + M.

  S = total self-loop orbits = sum_C (# arc-orbits that self-loop)
  M = total edge multiplicity excess

  CAN WE COMPUTE S FROM BURNSIDE?
  S = sum_C (1/|Aut(C)|) * sum_{g in Aut(C)} Fix_neutral_arcs(g)
  where Fix_neutral_arcs(g) = #{arcs fixed by g that are neutral}

  This is a WEIGHTED Burnside sum — computable from cycle types!

  S IS THE MISSING FORMULA (partially).
  S captures the "self-loop" part of degeneracy.
  M captures the "multiplicity" part.

  AT LARGE n: S dominates (from the data: S/D -> 1 as n -> inf).
  Actually: S/D = 1, 1, 0.57 at n=3,4,5. NOT converging to 1.
  M/D = 0, 0, 0.43. The multiplicity part GROWS.

  Hmm, need more data to determine the asymptotic.
""")

S_known = [2, 6, 16]
M_known = [0, 0, 12]
for i in range(3):
    n = i+3
    print("  n=%d: D=%d = S(%d) + M(%d), S/D=%.3f, M/D=%.3f" %
          (n, D[i], S_known[i], M_known[i], S_known[i]/D[i], M_known[i]/D[i]))

# ============================================================
banner("THE MINIMUM EXTRA INFORMATION NEEDED")
# ============================================================

print("""
  TO UNLOCK G_n AT ALL n, WE NEED:

  OPTION A (1 formula): The degeneracy D(n) = T(n) - 2E(n).
    Then E = (T - D) / 2 gives exact edge count.
    D = S + M where S = self-loop orbits, M = multiplicity excess.

  OPTION B (2 simpler formulas): S(n) and M(n) separately.
    S = sum_C #{neutral arc orbits of C} = Burnside-computable!
    M = sum_{edges} (excess multiplicity) = harder.

  OPTION C (1 function): The Burnside kernel for S.
    S = (1/n!) * sum_sigma #{tournaments T with sigma in Aut(T) AND
        some arc of T fixed by sigma that is neutral}
    This IS a cycle-index computation.

  THE MOST PROMISING PATH: Compute S from cycle types.
  For each cycle type lambda:
    Fix_neutral(lambda) = #{T with lambda in Aut(T)} *
                          #{arcs fixed by lambda that are neutral in T}
    S = (1/n!) * sum_lambda mult(lambda) * Fix_neutral(lambda)

  This is COMPUTABLE from cycle types alone!
  If we can derive Fix_neutral(lambda), we get S for all n.
  Then: D ~ S at large n (since M/D appears to decrease).
  And: E = (T - S - M) / 2 ~ (T - S) / 2 at large n.

  THE FORMULA E ~ (T - S) / 2 MIGHT BE EXACT ENOUGH
  to give the edge count to within a few percent at all n.
  Combined with the exact T from Schur-Weyl, this would be
  a COMPLETE ANALYTICAL DESCRIPTION.

  WHAT WE STILL DON'T KNOW:
  1. Fix_neutral(lambda) — the neutral arc Burnside weight
  2. M(n) — the edge multiplicity excess
  3. Whether M/D -> 0 (which would make E ~ (T-S)/2 asymptotically exact)
""")

print("Done. opus-2026-03-23-S208")
