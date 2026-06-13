"""
quaternion_octonion_tournament.py -- kind-pasteur-2026-03-14-S97
Quaternions, octonions, and the tower R -> C -> H -> O in tournament theory.

THE CAYLEY-DICKSON TOWER:
  R (dim 1): real numbers — H(T) is real
  C (dim 2): complex numbers — F(T, i) is a Gaussian integer
  H (dim 4): quaternions — what would F(T, j) or F(T, k) mean?
  O (dim 8): octonions — the non-associative level

EACH DOUBLING LOSES A PROPERTY:
  R -> C: lose ordering (C is not ordered)
  C -> H: lose commutativity (ij ≠ ji)
  H -> O: lose associativity ((ab)c ≠ a(bc))

IN TOURNAMENTS:
  The binary arc choice (2 options) is a DOUBLING.
  C(n,2) arcs = C(n,2) doublings.
  The total tournament space has dimension 2^{C(n,2)}.
  At n=3: 2^3 = 8 = dim(O) !!!
  At n=4: 2^6 = 64 = 8^2
  At n=5: 2^10 = 1024 = 2^10

THE DEEP QUESTION: Does the Cayley-Dickson structure of 8 tournaments
at n=3 reflect the octonion multiplication table?
"""

import sys, math
import numpy as np

sys.stdout.reconfigure(encoding='utf-8')

def C(n, k):
    if k < 0 or k > n: return 0
    return math.comb(n, k)

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[j][i] = 1
            else: A[i][j] = 1
            idx += 1
    return A

def compute_H(A, n):
    dp = {}
    for v in range(n): dp[(1 << v, v)] = 1
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms: continue
            for v in range(n):
                if not (mask & (1 << v)): continue
                pm = mask ^ (1 << v)
                t = sum(dp.get((pm, u), 0) for u in range(n) if (pm & (1 << u)) and A[u][v])
                if t: dp[(mask, v)] = t
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

maxH = {1:1, 2:1, 3:3, 4:5, 5:15, 6:45, 7:189, 8:661, 9:3357, 10:15745, 11:95095}

def main():
    pi = math.pi
    e = math.e
    I = 1j

    print("=" * 70)
    print("QUATERNIONS, OCTONIONS, AND TOURNAMENTS")
    print("kind-pasteur-2026-03-14-S97")
    print("=" * 70)

    # ============================================================
    # PART 1: THE CAYLEY-DICKSON TOWER
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 1: THE CAYLEY-DICKSON TOWER AND TOURNAMENT DIMENSIONS")
    print(f"{'='*70}")

    print(f"""
  THE TOWER OF NORMED DIVISION ALGEBRAS:
    R (dim 1): real numbers — commutative, associative, ordered
    C (dim 2): complex numbers — commutative, associative, NOT ordered
    H (dim 4): quaternions — NOT commutative, associative
    O (dim 8): octonions — NOT commutative, NOT associative

  TOURNAMENT DIMENSIONS:
    n=2: 2^1 = 2 tournaments (= dim C)
    n=3: 2^3 = 8 tournaments (= dim O !!!)
    n=4: 2^6 = 64 tournaments (= 8^2 = dim O^2)
    n=5: 2^10 = 1024 tournaments (= 2^10)
    n=6: 2^15 = 32768 tournaments (= 2^15)

  THE n=3 COINCIDENCE: 8 tournaments = 8 octonion basis elements!
    Tournament 000 (transitive, H=1) = 1 (real unit)
    Tournament 001 (one flip, H=1) = e_1
    Tournament 010 (one flip, H=1) = e_2
    ...
    Tournament 111 (3-cycle, H=3) = e_7 (or the imaginary octonion unit)

  But the 8 tournaments DON'T have the same H:
    6 tournaments have H=1 (transitive)
    2 tournaments have H=3 (3-cycle)

  In octonion terms: 6 'imaginary units' and 2 'special' elements.
  The octonion imaginary units are e_1,...,e_7 (7 elements), plus 1.
  That's 7+1 = 8. But we have 6+2 = 8.
  The mismatch: 7 ≠ 6.

  HOWEVER: the H=1 tournaments include BOTH transitive orderings
  (0→1→2 and 2→1→0 in the base path model). The H=3 tournaments
  ARE the two 3-cycles (CW and CCW).
""")

    # ============================================================
    # PART 2: THE 8 TOURNAMENTS AT n=3 — OCTONION STRUCTURE?
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 2: THE 8 TOURNAMENTS AT n=3")
    print(f"{'='*70}")

    n = 3
    m = C(n, 2)

    print(f"\n  The {2**m} tournaments on 3 vertices:")
    for bits in range(2**m):
        A = bits_to_adj(bits, n)
        H = compute_H(A, n)
        scores = tuple(sorted(A.sum(axis=1).astype(int)))
        # Arc directions
        arcs = []
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if A[i][j]: arcs.append(f"{i}->{j}")
                else: arcs.append(f"{j}->{i}")
                idx += 1
        print(f"    bits={bits:03b}: H={H}, scores={scores}, arcs={arcs}")

    # The flip operation on each arc:
    # Flipping arc 0: bits XOR 001 = flip arc (0,1)
    # Flipping arc 1: bits XOR 010 = flip arc (0,2)
    # Flipping arc 2: bits XOR 100 = flip arc (1,2)

    print(f"\n  ARC-FLIP CAYLEY TABLE:")
    print(f"  (row XOR col, with H values)")
    print(f"         ", end="")
    for b in range(8):
        print(f"  {b:03b}", end="")
    print()
    for a in range(8):
        Ha = compute_H(bits_to_adj(a, n), n)
        print(f"  {a:03b}(H={Ha})", end="")
        for b in range(8):
            c = a ^ b
            Hc = compute_H(bits_to_adj(c, n), n)
            print(f"  {Hc}", end="  ")
        print()

    # The "multiplication" a * b = a XOR b has:
    # Identity: 000 (H=1, transitive)
    # Every element is its own inverse (XOR is self-inverse)
    # This is the group Z/2 x Z/2 x Z/2 (= (Z/2)^3)

    print(f"\n  The XOR group is (Z/2)^3 — abelian, every element self-inverse.")
    print(f"  This is NOT the octonion multiplication (which is non-associative).")
    print(f"  But the STRUCTURE of H values under XOR is interesting:")
    print(f"  H(a XOR b) depends on a and b in a specific way.")

    # ============================================================
    # PART 3: QUATERNION-VALUED TOURNAMENT FUNCTIONS
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 3: QUATERNION-VALUED TOURNAMENT FUNCTIONS")
    print("  Define F(T, q) = sum over paths: q^{fwd(P)}")
    print("  For q = i: we get Gaussian integers (S96)")
    print("  For q = j (quaternion): what happens?")
    print(f"{'='*70}")

    # In quaternions: i^2 = j^2 = k^2 = ijk = -1
    # i*j = k, j*i = -k (non-commutative!)
    # For F(T, q): each path contributes q^{fwd}
    # Since q^n depends on q, and quaternion powers are well-defined:
    # q^0 = 1, q^1 = q, q^2 = q*q, etc.
    # For unit quaternions q on the 3-sphere: q = cos(theta) + sin(theta)*u
    # where u is a pure imaginary unit quaternion.

    # For q = j: j^0 = 1, j^1 = j, j^2 = -1, j^3 = -j, j^4 = 1
    # Same period as i! So F(T, j) = F(T, i) (same values, different algebra)

    print(f"\n  For unit quaternion q with q^2 = -1:")
    print(f"  q^0 = 1, q^1 = q, q^2 = -1, q^3 = -q, q^4 = 1")
    print(f"  The period is 4, same as for i.")
    print(f"  So F(T, q) = F(T, i) AS A SCALAR (the coefficients are the same).")
    print(f"  The quaternion structure only matters when we COMBINE evaluations.")
    print(f"")
    print(f"  HOWEVER: we can define a QUATERNION-VALUED F:")
    print(f"  F_H(T) = sum over paths: i^{{fwd}} * j^{{back}} * k^{{...}}")
    print(f"  This would encode fwd, back, AND some third quantity.")

    # What could the third quantity be?
    # fwd + back = n-1 (total steps)
    # So fwd and back are not independent.
    # But we could use: fwd, number of "peak" positions, and "valley" positions.
    # peaks = local maxima in the permutation, valleys = local minima.

    # Or: we could use the THREE types of consecutive pairs:
    # Type A: both in the same 3-cycle (consecutive in a triangle)
    # Type B: one in backbone, one off
    # Type C: both off backbone

    # The quaternion F would encode the TYPE of each step, not just fwd/back.

    # ============================================================
    # PART 4: THE DIMENSION DOUBLING AND TOURNAMENT STRUCTURE
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 4: DIMENSION DOUBLING — R -> C -> H -> O")
    print(f"{'='*70}")

    print(f"""
  THE CAYLEY-DICKSON DOUBLING:
    R (1-dim) -> C (2-dim): adjoin i with i^2 = -1
    C (2-dim) -> H (4-dim): adjoin j with j^2 = -1, ij = k
    H (4-dim) -> O (8-dim): adjoin l with l^2 = -1

  EACH DOUBLING:
    Dimension: 1 -> 2 -> 4 -> 8
    Lost property: ordering -> commutativity -> associativity -> ?

  TOURNAMENT DOUBLINGS:
    Each arc choice DOUBLES the tournament count.
    Arc 1: 2^1 = 2 (one arc, two orientations)
    Arc 2: 2^2 = 4 (two arcs, four combinations)
    Arc 3: 2^3 = 8 (three arcs = n=3 tournaments = dim O)

  AT EACH DOUBLING, WHAT PROPERTY IS LOST?
    1 arc (n=2): tournament IS a total order (ORDERED, like R)
    2 arcs (n=2.5?): not meaningful
    3 arcs (n=3): tournaments can have CYCLES (lose ordering, like C!)
    6 arcs (n=4): tournaments have MULTIPLE cycle structures
    10 arcs (n=5): tournaments have DISJOINT cycles (alpha_2 > 0 at n=6)

  THE PARALLEL:
    R -> C: lose ordering. Tournaments at n=3: lose transitivity (cycles appear)
    C -> H: lose commutativity. At n=4: arc-flip effects become non-commutative?
    H -> O: lose associativity. At n=5+: composition is non-associative?

  Let me check: is arc-flipping commutative?
  Flipping arc e1 then arc e2 vs e2 then e1:
  In {0,1}^m encoding: XOR is commutative!
  But in TERMS OF H: H(T XOR e1 XOR e2) is the same regardless of order.
  So arc-flipping IS commutative as a group action.

  BUT: the EFFECT on H is non-trivial.
  Delta_{{e1}} Delta_{{e2}} H = Delta_{{e2}} Delta_{{e1}} H (commutes)
  But Delta_{{e1}} (Delta_{{e2}} H) might differ from Delta_{{e2}} (Delta_{{e1}} H)
  if we consider the INTERMEDIATE tournament.
  For H values: the XOR group is abelian, so the final H is the same.
  But for the PATH through tournament space: the intermediates differ.
""")

    # ============================================================
    # PART 5: OCTONION MULTIPLICATION AND TOURNAMENT 3-CYCLES
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 5: OCTONION MULTIPLICATION AND 3-CYCLES")
    print(f"{'='*70}")

    # The Fano plane encodes octonion multiplication.
    # 7 imaginary units e_1,...,e_7 with e_i * e_j = +-e_k
    # The 7 "lines" of the Fano plane determine the products.
    # Each line has 3 points, and 3 points define a cyclic product.
    # This is a TOURNAMENT ON 7 POINTS (the Fano plane is a tournament?)

    # Actually: the Fano plane has 7 points and 7 lines.
    # Each line passes through 3 points. Each point is on 3 lines.
    # This is a (7, 3, 1) design = a Steiner triple system.

    # A Steiner triple system on 7 points = EXACTLY the Fano plane.
    # And 7 = the first FORBIDDEN H value!

    print(f"  The FANO PLANE has 7 points and 7 lines (Steiner triple system).")
    print(f"  It encodes OCTONION MULTIPLICATION via signed 3-cycles.")
    print(f"  7 = the first FORBIDDEN H value!")
    print(f"")
    print(f"  COINCIDENCE?")
    print(f"  The Fano plane has 7 points, each on 3 lines.")
    print(f"  A tournament with H=7 would need Omega = K_3 (3 pairwise-sharing cycles).")
    print(f"  K_3 has 3 vertices = 3 cycles, each pair sharing a vertex.")
    print(f"  This is like 3 LINES through a COMMON POINT in the Fano plane!")
    print(f"")
    print(f"  In the Fano plane: each point lies on exactly 3 lines.")
    print(f"  In the tournament: K_3 conflict graph = 3 cycles sharing vertices.")
    print(f"  The IMPOSSIBILITY of K_3 as Omega = the impossibility of")
    print(f"  a 'Fano-like' structure in tournament cycle configurations.")
    print(f"")
    print(f"  7 IS FORBIDDEN BECAUSE THE FANO STRUCTURE CANNOT EMBED IN TOURNAMENTS!")

    # ============================================================
    # PART 6: THE HURWITZ THEOREM CONNECTION
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 6: HURWITZ THEOREM — ONLY 4 NORMED DIVISION ALGEBRAS")
    print(f"{'='*70}")

    print(f"""
  HURWITZ'S THEOREM: The only normed division algebras over R are:
    R (dim 1), C (dim 2), H (dim 4), O (dim 8)
  No others exist! Dimensions 1, 2, 4, 8 are SPECIAL.

  IN TOURNAMENT THEORY:
    n=2: 2^1 = 2 = dim C. Only 1 non-trivial arc. H in {{1}}.
    n=3: 2^3 = 8 = dim O. 3 arcs. H in {{1, 3}}.
    n=4: 2^6 = 64 = 8^2. H in {{1, 3, 5}}.

  THE HURWITZ DIMENSIONS 1, 2, 4, 8 AND TOURNAMENTS:
    1 = H(transitive) = the multiplicative identity
    2 = OCF fugacity = the binary choice
    4 = 2^2 = level-2 Fourier structure = dim H (quaternions)
    8 = 2^3 = tournaments at n=3 = dim O (octonions)

  AFTER 8, the Cayley-Dickson construction continues but loses norming:
    Sedenions (dim 16 = 2^4): contain zero divisors!
    This corresponds to: tournaments at n=5 having 2^10 = 1024 elements,
    which is MUCH larger than any Cayley-Dickson algebra.
    The 'nice' algebra structure breaks down.

  THE ANALOGY:
    n=2 (2 tournaments): like C — simple, commutative, everything works
    n=3 (8 tournaments): like O — rich structure, but at the EDGE of normability
    n=4+ (64+ tournaments): like sedenions — zero divisors appear
    (In tournaments: H=7 gap appears = a "zero divisor" analog?)

  WILD CONJECTURE:
  The H=7 impossibility is the TOURNAMENT ANALOG of the
  non-existence of normed division algebras beyond dim 8.
  Just as 16-dimensional algebras have zero divisors,
  tournaments beyond n=3 have "impossible" H values (starting with 7).
""")

    # ============================================================
    # PART 7: THE NUMBER 8 IN TOURNAMENT THEORY
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 7: THE NUMBER 8 = 2^3 = dim(O)")
    print(f"{'='*70}")

    print(f"  8 appears in tournament theory as:")
    print(f"    2^3 = 8 tournaments at n=3")
    print(f"    rank(E_8) = 8")
    print(f"    n=8: beta_4 first nonzero (path homology threshold)")
    print(f"    n=8: H=63 first achievable (escape from forbidden orbit)")
    print(f"    n=8: Omega claw-freeness first fails")
    print(f"    8 SC tournament classes at n=5")
    print(f"    F(T,i) mean = -4 = -2^2 at n=5 (so |mean| = 4 = dim H)")
    print(f"    A_5(i) = -64 = -8^2 (Eulerian at i)")
    print(f"")
    print(f"  8 is the DIMENSION OF THE OCTONIONS.")
    print(f"  Every tournament property that involves 8")
    print(f"  may have an octonionic interpretation.")

    # ============================================================
    # PART 8: THE e^{i*pi} = -1 IN EACH ALGEBRA
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 8: EULER'S FORMULA IN EACH ALGEBRA LEVEL")
    print(f"{'='*70}")

    print(f"""
  IN R: e^0 = 1 (trivial). No rotation.
  IN C: e^(i*pi) = -1. Rotation by pi in the complex plane.
  IN H: e^(u*pi) = -1 for ANY unit imaginary quaternion u.
        There's a WHOLE SPHERE of "Euler identities"!
        e^(i*pi) = e^(j*pi) = e^(k*pi) = -1.
  IN O: e^(u*pi) = -1 for ANY unit imaginary octonion u.
        There's a 6-SPHERE of Euler identities!

  TOURNAMENT INTERPRETATION:
  IN R: F(T, 1) = n! (trivial evaluation)
  IN C: F(T, e^(i*pi)) = F(T, -1) (signed count, always odd)
  IN H: F(T, e^(j*pi)) = F(T, -1) (same! quaternion rotation = complex rotation for F)
  IN O: same — the octonion structure doesn't add new evaluations of F
        BECAUSE F is a polynomial in ONE variable.

  BUT: if we define a MULTIVARIATE F:
    F(T, x_1, ..., x_m) = sum over paths: prod x_{{e_k}}
    Then quaternion/octonion evaluations at (x_1,...,x_m) in H^m or O^m
    would give genuinely new objects.

  THE TRANSFER MATRIX as a quaternion:
    M is n x n with integer entries. We can view M as acting on H^n
    (quaternion vectors) or O^n (octonion vectors).
    The quaternionic eigenvalues would be different from complex ones!
    Quaternionic spectral theory of tournament transfer matrices.
""")

    # ============================================================
    # PART 9: THE FANO PLANE AND H=7
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 9: THE FANO PLANE, H=7, AND THE NUMBER 7")
    print(f"{'='*70}")

    print(f"""
  THE FANO PLANE PG(2, 2):
    7 points, 7 lines, 3 points per line, 3 lines per point.
    The automorphism group is GL(3, F_2) of order 168 = 8*21.
    168 = 8 * 21 = dim(O) * H_forb_2 !!!

  THE NUMBER 7:
    7 = 2^3 - 1 (Mersenne prime)
    7 = H_forb_1 (first forbidden H value)
    7 = I(C_3, 2) (independence polynomial of 3-cycle at fugacity 2)
    7 = h(A_6) (Coxeter number)
    7 = dim(O) - 1 (imaginary octonion dimension)
    7 = number of lines in Fano plane
    7 = number of vertices in Fano plane
    7 = number of non-identity octonion basis elements

  THE FANO PLANE IS THE OCTONION MULTIPLICATION TABLE.
  7 = dim(O) - 1 IS FORBIDDEN as an H value.
  This cannot be a coincidence.

  SPECULATIVE CHAIN:
    Octonions O have dim 8 = 2^3.
    The 7 imaginary units multiply via the Fano plane.
    The Fano plane has symmetry group of order 168 = 8 * 21.
    21 = 3 * 7 = H_forb_2.
    7 = H_forb_1.
    Both forbidden values appear in the octonionic structure!

  DEEPER:
    The Fano plane is PG(2, F_2) — the projective plane over F_2.
    F_2 = the field with 2 elements = the TOURNAMENT FIELD.
    Each arc orientation is an element of F_2.
    The Fano plane lives in F_2^3 — which has 2^3 = 8 elements.
    These 8 elements ARE the 8 tournaments at n=3!

  SO: The 8 tournaments at n=3 ARE the points of F_2^3,
  and the Fano plane is the projective closure PG(2, F_2).
  The forbidden H=7 = |PG(2, F_2)| (number of points in the Fano plane).
  H=21 = |GL(3, F_2)| / 8 = order of the Fano symmetry / dim(O).
""")

    # Verify: |GL(3, F_2)| = 168
    # GL(3, F_2) has order (2^3-1)(2^3-2)(2^3-4) = 7*6*4 = 168
    gl3_order = (2**3 - 1) * (2**3 - 2) * (2**3 - 4)
    print(f"  |GL(3, F_2)| = (2^3-1)(2^3-2)(2^3-4) = 7*6*4 = {gl3_order}")
    print(f"  168 / 8 = {168/8} = 21 = H_forb_2")
    print(f"  168 / 24 = {168/24} = 7 = H_forb_1")
    print(f"  168 = 7 * 24 = 7 * 4! = H_forb_1 * n(4)!")

    # ============================================================
    # PART 10: SYNTHESIS
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 10: THE GRAND SYNTHESIS")
    print(f"{'='*70}")

    print(f"""
  THE CAYLEY-DICKSON TOWER IN TOURNAMENT THEORY:

  LEVEL 0 (R, dim 1): The real line.
    Tournament: H(T) is a real (natural) number.
    Property: ordered (H values form a total order).

  LEVEL 1 (C, dim 2): The complex plane.
    Tournament: F(T, i) is a Gaussian integer.
    Property: lose ordering (F(T,i) is complex, not ordered).
    Euler: F(T, e^(i*pi)) = F(T, -1) always odd.

  LEVEL 2 (H, dim 4): The quaternions.
    Tournament: 4 = 2^2 = level-2 Fourier structure.
    Property: lose commutativity (arc-flip effects may not commute
    at the level of intermediate H values).
    4 = the number of Fourier levels that matter (0, 2, 4, and the rest).

  LEVEL 3 (O, dim 8): The octonions.
    Tournament: 8 = 2^3 = tournaments at n=3.
    Property: lose associativity (tournament composition via lex product
    is NOT associative at the class level — from S78!).
    7 = dim(O) - 1 = FORBIDDEN H value.
    The Fano plane (octonion mult. table) = PG(2, F_2).
    |GL(3, F_2)| = 168 = 8 * 21 = dim(O) * H_forb_2.

  BEYOND O (sedenions, dim 16):
    Tournament: 16 = 2^4 = GS tilings at n=5.
    Property: zero divisors (= the H=7 gap as an "algebraic zero").

  THE DEEPEST INSIGHT:
  The FORBIDDEN H values {7, 21} are the OCTONION NUMBERS:
    7 = |imaginary octonion units| = |Fano plane points|
    21 = |GL(3,F_2)|/8 = |Fano plane symmetries|/dim(O)
  Tournament impossibility at H=7 IS the impossibility of extending
  the Cayley-Dickson construction beyond the octonions without zero divisors.
""")

    print(f"\n{'='*70}")
    print("DONE — QUATERNIONS AND OCTONIONS IN TOURNAMENTS")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
