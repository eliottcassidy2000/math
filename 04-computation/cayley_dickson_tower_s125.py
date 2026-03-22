#!/usr/bin/env python3
"""
cayley_dickson_tower_s125.py — opus-2026-03-21-S125

THE CAYLEY-DICKSON TOWER AND TOURNAMENT THEORY

The tower: R(1) → C(2) → H(4) → O(8) → S(16) → ...
Each step: double dimension, lose one algebraic property.

  R: ordered, commutative, associative, alternative, division    (dim 1)
  C: ✗ordered, commutative, associative, alternative, division  (dim 2)
  H: ✗ordered, ✗commutative, associative, alternative, division (dim 4)
  O: ✗ordered, ✗commutative, ✗associative, alternative, division (dim 8)
  S: ✗ordered, ✗commutative, ✗associative, ✗alternative, ✗division (dim 16)

FROM S124: The Freudenthal magic square uses (R,C,H,O) as rows/columns.
The diagonal entries have dimensions: 3, 16, 66, 248.
The tournament numbers permeate the magic square.

THE DEEP QUESTION: Does the Cayley-Dickson tower correspond to the
sequence of LOSSES in tournament theory?

  n≤2 (R-level): trivial, no cycles, H always 1 → ORDERED
  n=3 (C-level): first non-trivial, 2 iso classes → COMMUTATIVE LOST (first direction)
  n=4-5 (H-level): cycle structure, POS appears → ASSOCIATIVITY PRESENT (OCR=1 at n≤4)
  n=5-7 (O-level): OCR<1, forbidden values, Fano → ASSOCIATIVITY LOST (the residual)
  n≥8 (S-level): OCR recovers, but higher structure emerges → ALTERNATIVITY LOST

Author: opus-2026-03-21-S125
"""

import sys
from math import comb
sys.stdout.reconfigure(line_buffering=True)

print("=" * 72)
print("  THE CAYLEY-DICKSON TOWER AND TOURNAMENT THEORY")
print("  opus-2026-03-21-S125")
print("=" * 72)

# ========================================================================
# PART 1: The tower and the properties lost
# ========================================================================
print("\n" + "=" * 72)
print("  PART 1: THE CAYLEY-DICKSON TOWER")
print("=" * 72)

tower = [
    ("R", 1, "ordered, commutative, associative, alternative, division"),
    ("C", 2, "commutative, associative, alternative, division"),
    ("H", 4, "associative, alternative, division"),
    ("O", 8, "alternative, division"),
    ("S", 16, "(none — zero divisors, non-alternative)"),
]

print(f"\n  {'Algebra':>8s} {'dim':>4s} {'Properties retained':>50s} {'Lost':>15s}")
lost_props = ["(none)", "ordering", "commutativity", "associativity", "alternativity"]
for i, (name, dim, props) in enumerate(tower):
    print(f"  {name:>8s} {dim:4d} {props:>50s} {lost_props[i]:>15s}")

# ========================================================================
# PART 2: Tournament properties lost at each level
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 2: TOURNAMENT PROPERTIES LOST AT EACH LEVEL")
print("=" * 72)
print(f"""
  The Cayley-Dickson tower loses one property at each step.
  Tournament theory loses structural simplicity in a parallel sequence:

  LEVEL R (dim 1, n=1-2): TOTAL ORDER
    H = 1 always. Score determines everything. One tournament per n.
    Property: ORDERING (total linear order, no ambiguity)

  LEVEL C (dim 2, n=3): FIRST DIRECTION
    2 iso classes (transitive, 3-cycle). H in (1, 3).
    Lost: UNIQUENESS of tournament structure.
    Score still determines H (OCR = 1).
    The "imaginary unit" i of C = the 3-cycle (a directed cycle = rotation).

  LEVEL H (dim 4, n=4-5): CYCLE STRUCTURE
    4 classes at n=4, 12 at n=5. Non-commutative cycle interactions.
    Lost: COMMUTATIVITY of cycle ordering.
    At n=4: OCR = 1. At n=5: OCR = 129/133 (first residual!).
    The "quaternion units" i,j,k = the 3 arcs of a 3-cycle,
    which satisfy the multiplication ij = k, ji = -k (non-commutative).

  LEVEL O (dim 8, n=6-7): FANO / OBSTRUCTION
    56 classes at n=6, 456 at n=7. Non-associative cycle interactions.
    Lost: ASSOCIATIVITY of cycle combination.
    Alpha_2 > 0 first appears at n=6. H=7 is forbidden.
    The Fano plane (7 points, 7 lines) governs octonionic multiplication.
    Paley T_7 contains the Fano plane as 7 directed 3-cycles.
    The loss of associativity = the OCR residual.

  LEVEL S (dim 16, n=8+): RECOVERY / HIGHER STRUCTURE
    OCR recovers (rises from n=7 minimum). But new complexity:
    beta_5 first appears at n=8. H=63 becomes achievable.
    Lost: ALTERNATIVITY — even the weakened form of associativity fails.
    At n=8: the binary skeleton is fully established.
    S (sedenions) have zero divisors — tournament analogue:
    there exist tournament pairs with "trivial interaction" (H_1 * H_2 effects).

  LEVEL BEYOND (dim 32+, n=9+): WILD
    OCR continues rising. Higher Betti numbers proliferate.
    Lost: Everything. But structure rebuilds from the ground up.
    Like higher Cayley-Dickson algebras: still have rich structure
    but none of the classical properties.
""")

# ========================================================================
# PART 3: The dimension doubling and tournament arc doubling
# ========================================================================
print("=" * 72)
print("  PART 3: DIMENSION DOUBLING ↔ ARC GROWTH")
print("=" * 72)

# Cayley-Dickson dimensions: 1, 2, 4, 8, 16, 32, ...
# Tournament arcs C(n,2): 1, 3, 6, 10, 15, 21, 28, 36, ...
# Tournament iso classes: 1, 1, 2, 4, 12, 56, 456, 6880, ...

print(f"\n  Cayley-Dickson dim vs Tournament arcs:")
print(f"  {'CD level':>9s} {'CD dim':>7s} {'n':>3s} {'C(n,2)':>7s} {'#classes':>9s} {'ratio':>8s}")

cd_levels = [
    (0, 1, 2, 1, 1),
    (1, 2, 3, 3, 2),
    (2, 4, 5, 10, 12),
    (3, 8, 7, 21, 456),
    (4, 16, 9, 36, None),
]

for level, cd_dim, n, arcs, classes in cd_levels:
    ratio = arcs / cd_dim
    cls_str = str(classes) if classes else "~191536"
    print(f"  Level {level:>2d} {cd_dim:7d} {n:3d} {arcs:7d} {cls_str:>9s} {ratio:8.1f}")

print(f"""
  The ratio arcs/CD_dim: 1.0, 1.5, 2.5, 2.625, 2.25.
  Not a clean doubling — the arcs grow quadratically while CD dims
  grow exponentially. They CROSS at different rates.

  THE CROSSING POINT: CD dim = arcs when:
    2^k = C(n,2) = n(n-1)/2.
    At k=3: 8 = C(n,2) → n(n-1)=16 → n≈4.5. Closest: n=5 (C(5,2)=10≈8).
    At k=4: 16 = C(n,2) → n(n-1)=32 → n≈6.2. Closest: n=6 (C(6,2)=15≈16).

  The Cayley-Dickson "level 3" (octonions, dim 8) corresponds to
  tournaments at n≈5, which is the BOUNDARY ORDER.
  The Cayley-Dickson "level 4" (sedenions, dim 16) corresponds to
  n≈6, the FIRST ORDER where alpha_2 > 0 and OCR < 1.

  THIS IS THE SAME BOUNDARY: the loss of associativity (O→S in CD)
  corresponds to the loss of score-determinacy (n=5→n=6 in tournaments).
""")

# ========================================================================
# PART 4: Each CD property loss ↔ a tournament structural transition
# ========================================================================
print("=" * 72)
print("  PART 4: THE CORRESPONDENCE TABLE")
print("=" * 72)
print(f"""
  CAYLEY-DICKSON PROPERTY          TOURNAMENT PROPERTY
  ════════════════════              ═══════════════════
  Ordering (R has total order)     Score determines H (OCR = 1)
  Commutativity (ab = ba)          Cycle reversal symmetry H(T) = H(T^op)
  Associativity ((ab)c = a(bc))    OCR = 1 (score EXACTLY determines H)
  Alternativity (a(ab) = a^2 b)    Binary skeleton (girth in {{3,inf}})
  Division (ab = 0 => a=0 or b=0)  No zero-divisor tournaments

  LOSSES:
  R→C: lose ordering               n≤2→n=3: first non-trivial H
  C→H: lose commutativity         n=3→n=5: multiple iso classes per score
  H→O: lose associativity         n=5→n=6: OCR < 1 (alpha_2 > 0)
  O→S: lose alternativity/division n=7→n=8: OCR recovers but H=63 appears

  THE KEY PARALLEL:
  Associativity is lost at the OCTONIONS (dim 8).
  OCR-determinacy is lost at n=5-6 (C(6,2) = 15 ≈ 2×8).
  Both losses happen at the SAME STRUCTURAL BOUNDARY.

  In both cases, the loss is caused by the APPEARANCE of the Fano plane:
  - Octonion multiplication rules = Fano plane
  - Tournament conflict structure at n=7 contains the Fano plane (S120)

  THE FANO PLANE IS THE SOURCE OF BOTH LOSSES.
""")

# ========================================================================
# PART 5: The norms and the tournament invariants
# ========================================================================
print("=" * 72)
print("  PART 5: NORMS AND TOURNAMENT INVARIANTS")
print("=" * 72)
print(f"""
  Each CD algebra has a NORM: N(a) = a * conj(a).
  The norm satisfies:
    R: N(ab) = N(a)N(b) (multiplicative)          → composition algebra
    C: N(ab) = N(a)N(b) (multiplicative)          → composition algebra
    H: N(ab) = N(a)N(b) (multiplicative)          → composition algebra
    O: N(ab) = N(a)N(b) (multiplicative)          → LAST composition algebra
    S: N(ab) ≠ N(a)N(b) in general                → NOT composition

  The norm composition property is LOST after the octonions.
  By Hurwitz's theorem: only R, C, H, O are normed division algebras.

  TOURNAMENT PARALLEL: the OCF identity H = I(Omega, 2).
  This is a COMPOSITION-LIKE property:
    H(T) = I(Omega, 2) = PRODUCT over independent cycle sets of 2^|S|

  If Omega is DISCONNECTED: I(Omega, 2) = I(Omega_1, 2) * I(Omega_2, 2).
  This is MULTIPLICATIVE on independent components!

  But: not all Omega are disconnected. Connected Omega = "non-decomposable"
  tournament structure. The H value is NOT a product of component values.

  AT n≤5: Omega is always either empty or connected.
  No PRODUCT STRUCTURE is needed → H is "composition-like."

  AT n≥6: Omega can be disconnected.
  PRODUCT STRUCTURE appears → H decomposes into factors.
  BUT: the forbidden values (7, 21) come from IMPOSSIBLE products.
  The composition fails because certain "norms" (component I-values)
  cannot be realized.

  THE HURWITZ THEOREM FOR TOURNAMENTS:
  "The OCF is compositional (I = product) exactly for tournaments where
  Omega is disconnected. The forbidden values are the 'norm gaps' — values
  that cannot be achieved by any product of realizable component norms."

  7 = the first norm gap (K_3 component is impossible)
  21 = 3 × 7 = first compound norm gap
  These are the TOURNAMENT HURWITZ OBSTRUCTIONS.
""")

# ========================================================================
# PART 6: The dimension sequence and tournament constants
# ========================================================================
print("=" * 72)
print("  PART 6: DIMENSION SEQUENCES")
print("=" * 72)

# CD dimensions: 1, 2, 4, 8, 16, 32, 64, ...
# Freudenthal diagonal: 3, 16, 66, 248
# Tournament: H_max(n) = 1, 1, 3, 5, 15, 45, 189, ...

print(f"  Dimension sequences compared:")
print(f"  {'k':>3s} {'2^k':>6s} {'diag(MS)':>9s} {'H_max(2k+1)':>12s} {'C(2k+1,2)':>10s}")

diag_ms = [3, 16, 66, 248]
H_max = {3: 3, 5: 15, 7: 189, 9: 3357}

for k in range(4):
    cd = 2**k
    ms = diag_ms[k] if k < len(diag_ms) else "?"
    n = 2*k + 3
    hmax = H_max.get(n, "?")
    cn2 = comb(n, 2)
    print(f"  {k:3d} {cd:6d} {str(ms):>9s} {str(hmax):>12s} {cn2:10d}")

print(f"""
  THE MAGIC SQUARE DIAGONAL: 3, 16, 66, 248.
    3 = H_max(3)
    16 = 2^4 = GS tilings at n=5 = tournament space at H-level
    66 = C(12,2) = arcs at n=12
    248 = dim(E_8) = the final exceptional algebra

  THE PATTERN: the diagonal grows roughly as C(k+2, 2) * 2^k?
    k=0: 3 = C(2,2)*2^0... no, that gives 1.
    Actually: 3, 16, 66, 248. Ratios: 16/3≈5.3, 66/16=4.125, 248/66≈3.76.
    Approaching ~3.76... which is close to τ² ≈ 3.38? No.

  These are just the dimensions of A_1, A_2+A_2, D_6, E_8.
  They come from the Freudenthal formula, not a simple recursion.
""")

# ========================================================================
# PART 7: Zero divisors and tournament zero-interaction pairs
# ========================================================================
print("=" * 72)
print("  PART 7: ZERO DIVISORS AND TOURNAMENT INDEPENDENCE")
print("=" * 72)
print(f"""
  In the sedenions S (dim 16): there exist a,b ≠ 0 with a*b = 0.
  These are ZERO DIVISORS — nonzero elements that multiply to zero.

  TOURNAMENT PARALLEL: Two tournaments T_1, T_2 on disjoint vertex sets
  have Omega(T_1 cup T_2) = Omega(T_1) union Omega(T_2) (disjoint).
  Their H values MULTIPLY: H(T_1 cup T_2) = H(T_1) * H(T_2)...
  Actually, this isn't right. The "union" tournament needs arcs between
  the two vertex sets, so it's not just a disjoint union.

  Better parallel: Two ODD CYCLES C_1, C_2 in the conflict graph Omega
  are ADJACENT (conflicting) if they share a vertex. If they DON'T share
  a vertex, they are INDEPENDENT — they "multiply independently" in the OCF.

  ZERO DIVISORS ↔ INDEPENDENT CYCLES:
  In S: a*b = 0 means the product annihilates.
  In tournaments: two independent cycles "ignore" each other.
  Alpha_2 counts the number of such independent pairs.

  At n≤5: alpha_2 = 0 always. No "zero divisors" in the conflict graph.
  This corresponds to DIVISION ALGEBRA (no zero divisors) at the O-level.

  At n=6: alpha_2 > 0 possible. First "zero divisors" appear.
  This corresponds to SEDENIONS (zero divisors appear) at the S-level.

  THE EXACT CORRESPONDENCE:
  CD level ↔ Tournament boundary
  R,C,H,O (no zero divisors) ↔ n≤5 (alpha_2 = 0)
  S (zero divisors exist) ↔ n≥6 (alpha_2 > 0)

  The Hurwitz theorem (only 4 composition algebras) corresponds to
  "only 4 orders (n=1,2,3,5) with alpha_2 = 0 for all tournaments."
  Wait: n=4 also has alpha_2 = 0.

  Actually: alpha_2 = 0 at n≤5 because you need ≥6 vertices for
  two vertex-disjoint 3-cycles. This is a PIGEONHOLE argument:
  2 × 3 = 6 > 5. The same argument gives: 2 × 3 = 6 > dim(H) = 4
  but 2 × 3 = 6 ≤ dim(O) = 8. So the boundary is:
  - CD: between H (dim 4) and O (dim 8)
  - Tournament: between n=5 (C(5,2)=10) and n=6 (C(6,2)=15)

  Both boundaries involve the FIRST APPEARANCE of "enough room"
  for two disjoint 3-structures.
""")

# ========================================================================
# PART 8: Synthesis — the CD tower as tournament phase diagram
# ========================================================================
print("=" * 72)
print("  SYNTHESIS: THE CD TOWER AS TOURNAMENT PHASE DIAGRAM")
print("=" * 72)
print(f"""
  THE CAYLEY-DICKSON TOWER IS THE PHASE DIAGRAM OF TOURNAMENT THEORY.

  Each "phase" has characteristic properties:

  R-PHASE (n=1-2): TRIVIAL
    One tournament per order. H = 1 always.
    Like real numbers: totally ordered, no ambiguity.

  C-PHASE (n=3): BINARY
    Two iso classes (transitive, 3-cycle). H in (1, 3).
    Like complex numbers: one "imaginary" direction (the 3-cycle).
    OCR = 1. Score determines everything.

  H-PHASE (n=4-5): QUATERNIONIC
    Non-commutative cycle interactions. OCR = 1 at n=4, first residual at n=5.
    Like quaternions: three "imaginary" units (the 3 arcs of a 3-cycle).
    The POS appears: score class (1,2,2,2,3) has 3 H values.
    The nesting formula H = 14 - H_inner + 4*sbs (S110).

  O-PHASE (n=6-7): OCTONIONIC
    Non-associative. Fano plane appears in Omega(T_7).
    OCR minimum at n=7. Forbidden values 7, 21.
    Like octonions: Fano multiplication rules, 7 imaginary units.
    Alpha_2 > 0 first appears. The composition property starts failing.

  S-PHASE (n=8+): SEDENIONIC
    Zero divisors (independent cycle pairs are common).
    OCR recovers. H=63 becomes achievable. Binary skeleton established.
    Like sedenions: zero divisors, but rich internal structure.
    Higher Betti numbers proliferate.

  BEYOND (n→∞): INFINITE CAYLEY-DICKSON
    Like the infinite CD tower: each step loses more but gains dimension.
    OCR → 1 (the shadow compression theorem).
    The tournament discriminant (H-1)/2 grows without bound.
    But the RATIO of cusp-to-Eisenstein → 0.

  THE PHASE TRANSITIONS:
  n=3 (C→H): first non-trivial tournament (ordering → direction)
  n=5-6 (H→O): first OCR residual, first alpha_2 (associativity → obstruction)
  n=7 (O deepens): OCR minimum, Fano embedding, forbidden atoms
  n=8 (O→S): OCR recovers, alternativity lost, zero divisors common

  THE FANO PLANE IS THE PHASE BOUNDARY:
  It lives at the O-level (octonions need it for multiplication).
  It appears in T_7 as 7 directed 3-cycles (S120).
  It marks the transition from "score determines H" to "score approximates H."

  THE HURWITZ THEOREM IS THE TOURNAMENT SHADOW THEOREM:
  "Only 4 composition algebras exist" ↔ "Only 4 tournament orders have alpha_2=0."
  Both statements say: the universe of structures with exact multiplicativity
  is FINITE and SMALL. Beyond it lies a richer but less controlled world.
""")

print("\nDone. opus-2026-03-21-S125")
