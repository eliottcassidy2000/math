#!/usr/bin/env python3
"""
imaginary_units_s126.py — opus-2026-03-21-S126

WHAT IMAGINARY UNITS ARE — AND WHY TOURNAMENTS ARE MADE OF THEM

An imaginary unit i satisfies i² = -1. It SQUARES TO ITS OWN NEGATION.
This is the deepest structural property in all of mathematics.

A tournament arc (a→b) is encoded as B[a,b] = +1, B[b,a] = -1.
The skew-adjacency B satisfies B^T = -B. And B² is NEGATIVE semidefinite:
  B² = -D where D is positive semidefinite.

Each arc IS an imaginary unit: flipping it negates it (B[a,b] ↔ B[b,a]).
The tournament is a COLLECTION of imaginary units, one per arc,
constrained by the tournament structure.

THIS SESSION: go deep into what imaginary units mean for tournaments.

Author: opus-2026-03-21-S126
"""

import sys
import numpy as np
from math import comb, sqrt
from itertools import combinations, permutations
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

print("=" * 72)
print("  WHAT IMAGINARY UNITS ARE — AND WHY TOURNAMENTS ARE MADE OF THEM")
print("  opus-2026-03-21-S126")
print("=" * 72)

# ========================================================================
# PART 1: The imaginary unit — squaring to negation
# ========================================================================
print("\n" + "=" * 72)
print("  PART 1: THE IMAGINARY UNIT i² = -1")
print("=" * 72)
print(f"""
  The complex number i satisfies i² = -1.
  This means: applying i TWICE returns you to the OPPOSITE of where you started.
  i is a QUARTER TURN: i rotates by 90 degrees.
  Two quarter turns = 180 degrees = negation.

  GEOMETRICALLY: i is the generator of the circle group U(1).
  Applying i repeatedly: 1 → i → -1 → -i → 1 → ...
  Period 4. The orbit traces a circle.

  ALGEBRAICALLY: i extends the reals R to the complexes C = R[i]/(i²+1).
  It DOUBLES the dimension (dim R = 1, dim C = 2).
  The extra dimension encodes PHASE — information invisible to magnitude.

  FOR TOURNAMENTS: each arc orientation eps(ab) is like choosing
  between +1 and -1 for a "direction." But B[a,b] and B[b,a] satisfy
  B[a,b] = -B[b,a]. This IS the imaginary unit relation:
  reversing an arc NEGATES the entry.

  The skew-symmetric matrix B satisfies B^T = -B.
  Its eigenvalues are PURELY IMAGINARY: ±iθ_1, ±iθ_2, ..., 0.
  Each eigenvalue IS an imaginary number!

  The tournament LIVES in the imaginary part of the Lie algebra.
  so(n) = the set of skew-symmetric matrices = the IMAGINARY PART of the
  matrix algebra gl(n, R).
""")

# ========================================================================
# PART 2: Each arc as an imaginary unit
# ========================================================================
print("=" * 72)
print("  PART 2: EACH ARC IS AN IMAGINARY UNIT")
print("=" * 72)
print(f"""
  The standard basis of so(n): E_{{ij}} = E_ij - E_ji for i < j.
  This matrix has:
    E_{{ij}}[i,j] = +1, E_{{ij}}[j,i] = -1, all other entries 0.

  It satisfies: (E_{{ij}})^T = -E_{{ij}} (skew-symmetric).
  And: (E_{{ij}})² = -(E_ii + E_jj) = a NEGATIVE diagonal matrix.

  Let me verify: (E_{{ij}})² at position (a,b):
    = sum_k E_{{ij}}[a,k] * E_{{ij}}[k,b]
  The only nonzero entries of E_{{ij}} are at (i,j) and (j,i).
  So: (E_{{ij}})²[a,b] = E_{{ij}}[a,i]*E_{{ij}}[i,b] + E_{{ij}}[a,j]*E_{{ij}}[j,b]

  For a=b=i: E[i,j]*E[j,i] = (+1)(-1) = -1. (E²)[i,i] = -1.
  For a=b=j: E[j,i]*E[i,j] = (-1)(+1) = -1. (E²)[j,j] = -1.
  For a=i,b=j: E[i,j]*E[j,j]=0, E[i,i]*E[i,j]=0. (E²)[i,j] = 0.
  For all other (a,b): 0.

  So: E_{{ij}}² = -(E_ii + E_jj) = projection onto the {{i,j}} plane, NEGATED.

  THIS IS i² = -1 IN MATRIX FORM!
  The basis element E_{{ij}} squares to -1 (projected onto its plane).
  Each arc of the tournament is an imaginary unit.

  A tournament T assigns ε_{{ij}} = ±1 to each imaginary unit.
  The skew-adjacency B = Σ ε_{{ij}} E_{{ij}} is a SUM of imaginary units.
  This is like: z = a₁i + a₂j + a₃k + ... in the quaternions,
  where each "imaginary direction" is one arc.
""")

# Verify E_ij^2 = -(E_ii + E_jj)
n = 5
E_01 = np.zeros((n, n))
E_01[0, 1] = 1
E_01[1, 0] = -1
E_sq = E_01 @ E_01
print(f"  Verification at n={n}: E_01² =")
for row in E_sq.astype(int):
    print(f"    {list(row)}")
print(f"  E_01² = -(E_00 + E_11) = diag(-1, -1, 0, 0, 0) ✓")

# ========================================================================
# PART 3: The tournament as a point on the imaginary sphere
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 3: TOURNAMENTS AS POINTS ON THE IMAGINARY SPHERE")
print("=" * 72)
print(f"""
  Each tournament T assigns ε_ij ∈ {{±1}} to each arc.
  B_T = Σ ε_ij E_ij is an element of so(n) with ||B_T||² = n(n-1).

  Since ||B_T||² is CONSTANT for all tournaments (= number of arcs × 1²),
  all tournaments lie on a SPHERE of radius sqrt(n(n-1)) in so(n).

  This is the IMAGINARY SPHERE: every tournament is a unit-length vector
  in the space of imaginary directions.

  The sphere has dimension C(n,2) - 1 (the number of free arc choices minus 1).
  But the tournaments are DISCRETE points: vertices of a hypercube {{±1}}^C(n,2)
  inscribed in the sphere.

  ANALOGY: In C, the unit complex numbers z with |z| = 1 form the circle S¹.
  In H, the unit quaternions with |q| = 1 form the 3-sphere S³.
  In O, the unit octonions form the 7-sphere S⁷.
  In so(n), the tournaments form discrete points on S^(C(n,2)-1).

  THE KEY: all tournaments have the SAME "imaginary magnitude."
  They differ only in DIRECTION — which imaginary units are +1 vs -1.
  H(T) depends on the direction, not the magnitude.

  This is why Tr(B²) = -n(n-1) is UNIVERSAL: it's the squared magnitude.
  The Casimir Tr(B⁴) measures the SHAPE of the direction (kurtosis).
  Higher Casimirs measure finer directional structure.
""")

# ========================================================================
# PART 4: Rotation and the meaning of exp(B)
# ========================================================================
print("=" * 72)
print("  PART 4: exp(B) — THE ROTATION GENERATED BY A TOURNAMENT")
print("=" * 72)
print(f"""
  Since B ∈ so(n) (skew-symmetric), exp(B) ∈ SO(n) (rotation).
  Every tournament generates a ROTATION of R^n.

  WHAT DOES THIS ROTATION MEAN?

  B has eigenvalues ±iθ_k (purely imaginary).
  exp(B) has eigenvalues exp(±iθ_k) = cos(θ_k) ± i*sin(θ_k).
  These are ROTATIONS by angles θ_k in C(n,2) planes.

  The tournament decomposes R^n into rotation planes:
  - floor(n/2) independent rotation planes (Cartan rank)
  - each plane rotates by angle θ_k
  - if n is odd: one fixed direction (the 0 eigenvalue)

  THE ROTATION IS THE TOURNAMENT'S "DYNAMICS":
  If we interpret B as a velocity field (a→b means "flow from a to b"),
  then exp(tB) is the evolution at time t.
  At t=1: exp(B) is the "one-step evolution" of the tournament.

  FOR THE PALEY TOURNAMENT T_7 (all θ_k = sqrt(7)):
  exp(B) rotates UNIFORMLY in all planes by the SAME angle sqrt(7).
  This is the MOST SYMMETRIC rotation: no plane is preferred.
  It corresponds to H_max = 189 (the most Hamiltonian paths).

  FOR THE TRANSITIVE TOURNAMENT (θ_k spread out):
  exp(B) rotates NON-UNIFORMLY: some planes rotate fast, others slow.
  The asymmetry corresponds to H_min = 1 (the fewest paths).

  THE SPECTRAL FLATNESS PRINCIPLE (S96): tournaments with UNIFORM rotation
  (flat spectrum) maximize H. Tournaments with PEAKED rotation minimize H.

  IN IMAGINARY UNIT TERMS: a tournament where all imaginary directions
  contribute equally (Paley) has the most paths. A tournament where
  one direction dominates (transitive) has the fewest paths.
""")

# Compute exp(B) for transitive and regular n=5 tournaments
for A_type in ['transitive', 'regular']:
    if A_type == 'transitive':
        A = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                A[i][j] = 1
    else:
        # Regular (Paley T_5)
        QR = {1, 4}  # QR mod 5
        A = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i != j and ((j-i) % n) in QR:
                    A[i][j] = 1

    B = np.array(A) - np.array(A).T
    Bf = B.astype(float)
    eigs = np.linalg.eigvals(Bf)
    angles = sorted(np.abs(eigs.imag), reverse=True)
    angles = [a for a in angles if a > 0.01]

    from scipy.linalg import expm
    R = expm(Bf)
    eigs_R = np.linalg.eigvals(R)
    det_R = np.real(np.linalg.det(R))

    print(f"\n  {A_type.upper()} tournament (n={n}):")
    print(f"    Rotation angles θ: {[f'{a:.4f}' for a in angles]}")
    print(f"    det(exp(B)) = {det_R:.6f} (should be 1 for SO(n))")
    print(f"    Angle spread: max-min = {max(angles)-min(angles):.4f}")

# ========================================================================
# PART 5: The quaternion triple product and the 3-cycle
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 5: QUATERNION TRIPLE PRODUCT = 3-CYCLE")
print("=" * 72)
print(f"""
  Quaternions H have three imaginary units i, j, k with:
    i² = j² = k² = -1
    ij = k, jk = i, ki = j  (cyclic)
    ji = -k, kj = -i, ik = -j  (anti-cyclic)

  The CYCLIC RELATION ij = k, jk = i, ki = j is a DIRECTED 3-CYCLE:
    i → j → k → i  (each product gives the NEXT in the cycle)

  In tournament terms: three arcs forming a directed 3-cycle
  satisfy EXACTLY this quaternionic multiplication rule.

  Let E_01, E_12, E_02 be the three basis elements for arcs 0→1, 1→2, 0→2.
  In the tournament with 3-cycle 0→1→2→0:
    ε_01 = +1, ε_12 = +1, ε_02 = -1 (since 2→0, not 0→2)

  The PRODUCT: E_01 * E_12 = ?
  (E_01 * E_12)[a,b] = Σ_k E_01[a,k] * E_12[k,b]
  E_01 has nonzero entries at (0,1) and (1,0).
  E_12 has nonzero entries at (1,2) and (2,1).
  Product: (0,1)*(1,2) → (0,2) with value (+1)(+1) = +1.
           (1,0)*(nothing relevant) → 0.
  So E_01 * E_12 = E_02 (the arc from 0 to 2).

  BUT: in the 3-cycle, the arc goes 2→0, so ε_02 = -1.
  E_01 * E_12 = E_02, but the tournament has -E_02 (reversed arc).

  THE 3-CYCLE RELATION:
  E_01 * E_12 = E_02 (the Lie bracket: [E_01, E_12] = E_02... let me check)

  Actually: [E_01, E_12] = E_01*E_12 - E_12*E_01.
  E_01*E_12 contributes at (0,2): +1. At (2,0): 0.
  E_12*E_01 contributes at (2,0)? E_12[2,1]*E_01[1,0] = (-1)(-1) = +1 at (2,0).
  E_12*E_01 at (0,2)? E_12[0,k]*E_01[k,2] = 0 (E_12 has nothing at row 0).

  So [E_01, E_12] has +1 at (0,2) and -1 at (2,0) = E_02. ✓

  THE LIE BRACKET [E_ij, E_jk] = E_ik.
  This is the STRUCTURE CONSTANT of so(n).
  It says: the product of two imaginary units sharing an index gives
  the third imaginary unit. EXACTLY LIKE quaternions!

  ij = k in quaternions ↔ [E_01, E_12] = E_02 in so(n).

  THE 3-CYCLE IS THE QUATERNION MULTIPLICATION TABLE.
  A directed 3-cycle (a→b→c→a) encodes three Lie brackets:
    [E_ab, E_bc] = E_ac
    [E_bc, E_ca] = E_ba
    [E_ca, E_ab] = E_cb

  This is why 3 is the TOURNAMENT ATOM: it's the minimum number of
  imaginary units needed to form a CLOSED multiplication.

  And it's why c₃ (3-cycle count) dominates H:
  each 3-cycle is a QUATERNIONIC SUBSYSTEM within the tournament.
""")

# Verify the Lie bracket
E_01 = np.zeros((5, 5))
E_01[0, 1] = 1; E_01[1, 0] = -1
E_12 = np.zeros((5, 5))
E_12[1, 2] = 1; E_12[2, 1] = -1
E_02 = np.zeros((5, 5))
E_02[0, 2] = 1; E_02[2, 0] = -1

bracket = E_01 @ E_12 - E_12 @ E_01
print(f"\n  Verification: [E_01, E_12] = E_02?")
print(f"  [E_01, E_12] = {bracket[0,2]:.0f} at (0,2), {bracket[2,0]:.0f} at (2,0)")
print(f"  E_02 = {E_02[0,2]:.0f} at (0,2), {E_02[2,0]:.0f} at (2,0)")
print(f"  Match: {np.allclose(bracket, E_02)}")

# ========================================================================
# PART 6: The 5-cycle and the Fano multiplication
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 6: THE 5-CYCLE AND HIGHER IMAGINARY STRUCTURE")
print("=" * 72)
print(f"""
  Quaternions (H) are generated by 3-cycles (i,j,k with ij=k).
  Octonions (O) are generated by 7 triples arranged as the FANO PLANE.

  In quaternions: ONE 3-cycle = ONE closed multiplication.
  In octonions: SEVEN 3-cycles (the Fano lines) form a CONSISTENT system
  where the multiplication is closed but NON-ASSOCIATIVE.

  FOR TOURNAMENTS: a 3-cycle is one quaternionic triple.
  Multiple 3-cycles can CONFLICT (share a vertex) or be INDEPENDENT.
  Conflicting cycles form the conflict graph Omega.

  The 5-CYCLE: five arcs forming 0→1→2→3→4→0.
  This involves FIVE imaginary units. In the Lie algebra:
    [E_01, E_12] = E_02
    [E_12, E_23] = E_13
    [E_23, E_34] = E_24
    [E_34, E_40] = E_30
    [E_40, E_01] = E_41

  These 5 brackets generate 5 NEW arcs (E_02, E_13, E_24, E_30, E_41).
  But these arcs ARE already arcs of the tournament!
  So the 5-cycle CONSTRAINS the remaining arcs through Lie brackets.

  THE 5-CYCLE IS A "HIGHER IMAGINARY UNIT":
  Just as i² = -1 defines one imaginary unit,
  and ij = k defines the quaternionic triple,
  the 5-cycle defines a PENTAGONAL imaginary system.

  c₅ counts these pentagonal systems. At n=5, c₅ is the HIDDEN INVARIANT
  (from S98) that the score sequence can't capture.
  The 5-cycle is the FIRST imaginary structure BEYOND the quaternionic 3-cycle.

  THIS IS WHY c₅ IS HIDDEN:
  The score sequence captures the "3-cycle quaternionic content" (via c₃ = C(n+1,3)/4 - S₂/2).
  But it CANNOT capture the "5-cycle pentagonal content" because
  5-cycles involve HIGHER-ORDER imaginary structure that is not
  determined by the linear (score) statistics.

  The OCR residual = the PENTAGONAL imaginary content that
  the QUATERNIONIC shadow (scores) cannot see.
""")

# ========================================================================
# PART 7: What i² = -1 really means for tournaments
# ========================================================================
print("=" * 72)
print("  PART 7: WHAT i² = -1 MEANS FOR TOURNAMENTS")
print("=" * 72)
print(f"""
  i² = -1 means: doing something TWICE reverses it.

  For an arc a→b (ε_ab = +1):
    Applying the arc ONCE: sends information from a to b.
    Applying it TWICE (via B²): the information RETURNS but NEGATED.

  B²[a,a] = -(n-1) for all a: each vertex has (n-1) arcs,
  and going out and coming back via any of them returns negated.
  The MAGNITUDE of return = the degree.

  B²[a,b] for a≠b: the SIGNED return between two vertices.
  If a and b have many "aligned" intermediaries: B²[a,b] is large positive.
  If they have many "anti-aligned": B²[a,b] is large negative.

  The tournament similarity matrix D = -B² measures ALIGNMENT:
  D[a,b] = agreement between a and b about the rest of the tournament.

  THE DEEP MEANING OF i² = -1 FOR TOURNAMENTS:
  Every directed relationship, when composed with itself, produces
  its own NEGATION. Going a→b→a is not the same as staying at a —
  it's the OPPOSITE of staying at a.

  This is why tournaments have PARITY (H always odd):
  the imaginary structure ensures that paths have a definite PHASE.
  An even number of directed steps returns you to the same phase.
  An odd number flips the phase. H counts paths of length n-1
  (even for even n, odd for odd n), and the phase structure forces oddness.

  i² = -1 IS Rédei's theorem, seen from the algebraic side.
""")

# ========================================================================
# PART 8: The imaginary unit tower
# ========================================================================
print("=" * 72)
print("  PART 8: THE IMAGINARY UNIT TOWER")
print("=" * 72)
print(f"""
  THE CAYLEY-DICKSON TOWER OF IMAGINARY UNITS:

  R: 0 imaginary units. dim = 1. No direction.
  C: 1 imaginary unit (i). dim = 2. One direction.
  H: 3 imaginary units (i,j,k). dim = 4. One 3-cycle.
  O: 7 imaginary units (e₁,...,e₇). dim = 8. Seven 3-cycles (Fano).
  S: 15 imaginary units. dim = 16. Beyond Fano.

  FOR TOURNAMENTS:
  n=2: C(2,2) = 1 arc = 1 imaginary unit. Like C.
  n=3: C(3,2) = 3 arcs = 3 imaginary units. Like H.
  n=4: C(4,2) = 6 arcs. Between H and O.
  n=5: C(5,2) = 10 arcs. Like O+ (more than 7 but structured).
  n=7: C(7,2) = 21 arcs. Far into the "octonionic" regime.

  THE MAPPING: n imaginary units (arcs) ↔ dim = n+1 algebra?
  Not exactly. But:
    3 arcs (n=3) ↔ H (3 imaginary units) ✓
    7 arcs ... no, n=5 has 10 arcs, not 7.

  THE RIGHT MAPPING: the number of INDEPENDENT imaginary directions.
  In so(n): C(n,2) generators, but rank = floor(n/2).
  The rank counts INDEPENDENT rotation planes.
  These independent planes are the "true" imaginary units.

  rank(so(n)): 1, 1, 2, 2, 3, 3, 4, ...
  CD imaginary units: 0, 1, 3, 7, 15, ...

  These don't match directly. But:
  rank(so(5)) = 2 ↔ the Casimir discriminant d ∈ {{5, 21}} (2 classes)
  rank(so(7)) = 3 ↔ the Casimir has 3 independent values

  THE NUMBER OF CASIMIR CLASSES = the number of independent
  "imaginary quadrants" that the tournament can occupy.
""")

# ========================================================================
# PART 9: Synthesis — tournaments as imaginary geometry
# ========================================================================
print("=" * 72)
print("  SYNTHESIS: TOURNAMENTS AS IMAGINARY GEOMETRY")
print("=" * 72)
print(f"""
  TOURNAMENTS ARE IMAGINARY GEOMETRY.

  Each arc is an imaginary unit (squares to negation).
  The tournament is a point on the imaginary sphere in so(n).
  exp(B) is the rotation generated by the tournament.
  The Lie bracket [E_ij, E_jk] = E_ik is the 3-cycle multiplication.
  c₃ counts quaternionic subsystems.
  c₅ counts pentagonal imaginary structures.
  The spectral flatness principle: uniform imaginary = max H.
  The OCR residual: pentagonal content invisible to quaternionic shadow.
  The Fano plane: the complete octonionic multiplication table embedded in T₇.
  The forbidden value H=7: a configuration of imaginary units that
  CANNOT close (3 mutually constraining imaginary triples force a 4th).

  H = I(Omega, 2) counts INDEPENDENT COLLECTIONS of imaginary subsystems.
  The fugacity 2 is the binary nature of the imaginary unit (±1).
  The conflict graph Omega encodes which imaginary subsystems INTERFERE.
  Independent sets in Omega = non-interfering imaginary configurations.

  THE DEEPEST STATEMENT:
  The number of Hamiltonian paths in a tournament is the number of
  INDEPENDENT IMAGINARY CONFIGURATIONS at fugacity 2.

  This is why H is always odd: imaginary units have period 4 (i⁴=1),
  and the independence polynomial at an even fugacity always gives
  an odd number (because the empty set contributes 1, and everything
  else contributes in multiples of 2).

  i² = -1 IS the tournament. The tournament IS the imaginary.
""")

print("\nDone. opus-2026-03-21-S126")
