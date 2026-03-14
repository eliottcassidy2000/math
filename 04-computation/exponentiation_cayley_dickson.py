"""
exponentiation_cayley_dickson.py — Deep exploration of 2^3=8, 3^2=9,
Hilbert's 3rd problem, omega powers, Cayley-Dickson, and n=8 as inflection point.
kind-pasteur-2026-03-14-S65

PART 1:  Exponentiation asymmetry: 2^3=8 vs 3^2=9
PART 2:  The gap 3^2 - 2^3 = 1 and Catalan's conjecture (Mihailescu)
PART 3:  Omega (cube root of unity) — powers and Z[omega] arithmetic
PART 4:  omega^2 = -omega - 1 and the "2/3 omega - 1" relation
PART 5:  Powers omega^1 through omega^8 — pigeonhole at dimension 2
PART 6:  Hilbert's 3rd problem: Dehn invariant and irrationality insertion
PART 7:  Alternating sums and non-negativity in tournament context
PART 8:  n=8 as tournament inflection point — the catalog of failures
PART 9:  Cayley-Dickson tower: R -> C -> H -> O -> S (1,2,4,8,16)
PART 10: Property loss at each doubling — tournament parallels
PART 11: 2^3 in the independence polynomial — I(8) and the octonionic threshold
PART 12: Exponentiation in H(T) — what is H^{1/3}? H^{1/2}?
PART 13: The 8-9 pair: consecutive powers and near-misses
PART 14: Synthesis — the trinity of operations on {2,3}
"""

import numpy as np
from itertools import combinations, permutations
from collections import defaultdict
import math

def random_tournament(n, rng):
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def count_ham_paths(A):
    n = len(A)
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask_size in range(2, n + 1):
        for mask in range(1 << n):
            if bin(mask).count('1') != mask_size:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                prev_mask = mask ^ (1 << v)
                total = 0
                for u in range(n):
                    if (prev_mask & (1 << u)) and A[u][v]:
                        total += dp.get((prev_mask, u), 0)
                if total > 0:
                    dp[(mask, v)] = total
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def count_directed_ham_cycles_sub(A_sub, m):
    if m < 3:
        return 0
    dp = {}
    start = 0
    dp[(1 << start, start)] = 1
    for mask_size in range(2, m + 1):
        for mask in range(1 << m):
            if bin(mask).count('1') != mask_size:
                continue
            if not (mask & 1):
                continue
            for v in range(m):
                if not (mask & (1 << v)):
                    continue
                if v == start and mask_size < m:
                    continue
                prev_mask = mask ^ (1 << v)
                if not (prev_mask & 1):
                    continue
                total = 0
                for u in range(m):
                    if (prev_mask & (1 << u)) and A_sub[u][v]:
                        total += dp.get((prev_mask, u), 0)
                if total > 0:
                    dp[(mask, v)] = dp.get((mask, v), 0) + total
    full = (1 << m) - 1
    count = 0
    for v in range(m):
        if A_sub[v][start]:
            count += dp.get((full, v), 0)
    return count

def get_all_odd_cycles(A):
    """Get all directed odd cycles as frozensets (vertex sets), with multiplicity."""
    n = len(A)
    all_cycles = []
    for length in range(3, n+1, 2):
        for verts in combinations(range(n), length):
            sub = [[A[verts[i]][verts[j]] for j in range(length)] for i in range(length)]
            hc = count_directed_ham_cycles_sub(sub, length)
            for _ in range(hc):
                all_cycles.append(frozenset(verts))
    return all_cycles

def get_alpha_1_2(A):
    """Get alpha_1 (number of odd cycles) and alpha_2 (number of disjoint pairs)."""
    cycles = get_all_odd_cycles(A)
    a1 = len(cycles)
    a2 = 0
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            if cycles[i].isdisjoint(cycles[j]):
                a2 += 1
    return a1, a2

def I_poly(a1, a2, x):
    return 1 + a1*x + a2*x**2

# ============================================================
# PART 1: Exponentiation asymmetry
# ============================================================
def part1():
    print("=" * 70)
    print("PART 1: EXPONENTIATION ASYMMETRY — 2^3=8 vs 3^2=9")
    print("=" * 70)

    print("\nThe three binary operations on {2,3}:")
    print(f"  Addition:        2+3 = 5,   3+2 = 5   (commutative)")
    print(f"  Multiplication:  2*3 = 6,   3*2 = 6   (commutative)")
    print(f"  Exponentiation:  2^3 = 8,   3^2 = 9   (NON-commutative!)")

    print(f"\nExponentiation BREAKS symmetry between 2 and 3.")
    print(f"  2^3 = 8,  3^2 = 9,  difference = {9-8} = 1")
    print(f"  This is Catalan's conjecture (proved by Mihailescu 2002):")
    print(f"  The ONLY consecutive perfect powers are 8 and 9.")

    print(f"\nIn tournament theory:")
    print(f"  2 = evaluation point (H = I(Omega, 2))")
    print(f"  3 = cycle length (3-cycles are fundamental)")
    print(f"  8 = 2^3 = first n where seesaw breaks, beta_3>1")
    print(f"  9 = 3^2 = 'square of triangles' — 9 = C(3+1,2)+3")

    print(f"\nThe operation hierarchy:")
    print(f"  +: 2+3=5 (additive son, prime)")
    print(f"  *: 2*3=6 (multiplicative son, composite = 2*3)")
    print(f"  ^: 2^3=8 (exponential son, 2^3), 3^2=9 (exponential son, 3^2)")
    print(f"  Three sons: 5, 6, 8 (or 5, 6, 9 depending on direction)")

    # Which direction does tournament theory prefer?
    print(f"\nTournament theory preference:")
    print(f"  n=8 is the inflection point (2^3 wins)")
    print(f"  n=9 is where real-roots first fail (3^2 = the response)")
    print(f"  The ASYMMETRY 2^3 < 3^2 mirrors the asymmetry in tournament theory:")
    print(f"  structural breaks happen at 8 BEFORE 9")

# ============================================================
# PART 2: Catalan's gap and near-misses
# ============================================================
def part2():
    print("\n" + "=" * 70)
    print("PART 2: THE GAP 3^2 - 2^3 = 1 AND CATALAN'S THEOREM")
    print("=" * 70)

    print("\nPerfect powers up to 100:")
    powers = set()
    for base in range(2, 11):
        for exp in range(2, 8):
            val = base ** exp
            if val <= 100:
                powers.add((val, base, exp))
    for val, base, exp in sorted(powers):
        print(f"  {base}^{exp} = {val}")

    print(f"\nConsecutive perfect powers: only (8, 9)")
    print(f"  Mihailescu (2002): x^p - y^q = 1 has unique solution 3^2 - 2^3 = 1")

    # Connection to tournament theory
    print(f"\nIn our context:")
    print(f"  H values are always odd: H = 1 + 2*(a1 + 2*a2)")
    print(f"  H=9 is achievable (e.g., a1=4, a2=0: H=9)")
    print(f"  H=8 is NEVER achievable (H always odd)")
    print(f"  So the Catalan pair {8,9} is split: 8 is structurally forbidden, 9 is allowed")
    print(f"  The gap between them (=1) is invisible from H's perspective")

    print(f"\n  But dimensionally: n=8 and n=9 are consecutive orders where")
    print(f"  different structural properties fail:")
    print(f"    n=8: beta_3>1, seesaw breaks, i_*-injectivity fails")
    print(f"    n=9: real roots of I(Omega, x) first fail, claw-free fails")

# ============================================================
# PART 3: Omega and Z[omega]
# ============================================================
def part3():
    print("\n" + "=" * 70)
    print("PART 3: OMEGA (CUBE ROOT OF UNITY) AND Z[omega]")
    print("=" * 70)

    omega = complex(-0.5, math.sqrt(3)/2)
    print(f"\nomega = e^(2*pi*i/3) = -1/2 + sqrt(3)/2 * i")
    print(f"  = {omega}")

    print(f"\nPowers of omega:")
    for k in range(7):
        w = omega ** k
        print(f"  omega^{k} = {w.real:8.4f} + {w.imag:8.4f}i  (|omega^{k}| = {abs(w):.4f})")

    print(f"\nKey relation: omega^2 + omega + 1 = 0")
    print(f"  => omega^2 = -omega - 1")
    val = omega**2 + omega + 1
    print(f"  Check: omega^2 + omega + 1 = {val.real:.6e} + {val.imag:.6e}i")

    print(f"\nThe user's relation 'omega^2 = 2/3 omega - 1':")
    user_val = (2/3)*omega - 1
    print(f"  2/3*omega - 1 = {user_val}")
    print(f"  omega^2 = {omega**2}")
    print(f"  These differ by: {abs(omega**2 - user_val):.6f}")
    print(f"  Standard: omega^2 = -omega - 1 = -(-1/2 + sqrt(3)/2 i) - 1")
    print(f"                    = 1/2 - sqrt(3)/2 i - 1 = -1/2 - sqrt(3)/2 i")

    # Z[omega] = {a + b*omega : a, b in Z}
    print(f"\nZ[omega] = {{a + b*omega : a, b in Z}}")
    print(f"  This is a Euclidean domain (Eisenstein integers)")
    print(f"  Norm: N(a + b*omega) = a^2 - ab + b^2")
    print(f"  Units: {{+-1, +-omega, +-omega^2}} (6 units)")

    print(f"\nI(omega) for a tournament with (a1, a2):")
    print(f"  I(omega) = 1 + a1*omega + a2*omega^2")
    print(f"           = 1 + a1*omega + a2*(-omega-1)")
    print(f"           = (1-a2) + (a1-a2)*omega")
    print(f"  Norm N(I(omega)) = (1-a2)^2 - (1-a2)(a1-a2) + (a1-a2)^2")
    print(f"                   = a1^2 - a1*a2 + a2^2 - a1 - a2 + 1  (Eisenstein norm)")

    # Some example evaluations
    print(f"\nEisenstein norm for small (a1, a2):")
    for a1 in range(0, 12):
        for a2 in range(0, 6):
            N = a1**2 - a1*a2 + a2**2 - a1 - a2 + 1
            if N <= 25:
                r = N % 3
                print(f"  (a1={a1:2d}, a2={a2:2d}): N = {N:3d}  (mod 3 = {r})", end="")
                if r == 2:
                    print("  <-- IMPOSSIBLE (never 2 mod 3)")
                else:
                    print()

# ============================================================
# PART 4: omega^2 = -omega - 1 and the minimal polynomial
# ============================================================
def part4():
    print("\n" + "=" * 70)
    print("PART 4: omega^2 = -omega - 1 AND THE MINIMAL POLYNOMIAL")
    print("=" * 70)

    omega = complex(-0.5, math.sqrt(3)/2)

    print(f"\nThe minimal polynomial of omega over Q: x^2 + x + 1 = 0")
    print(f"  This is Phi_3(x), the 3rd cyclotomic polynomial")
    print(f"  Phi_3(2) = 4 + 2 + 1 = 7")
    print(f"  Phi_3(3) = 9 + 3 + 1 = 13")
    print(f"  Phi_3(5) = 25 + 5 + 1 = 31")
    print(f"  Phi_3(6) = 36 + 6 + 1 = 43")

    print(f"\nThe '2/3 omega - 1' interpretation:")
    print(f"  If we write omega^2 = c*omega + d, then c = -1, d = -1")
    print(f"  But in F_3: omega^2 = -omega - 1 = 2*omega + 2 = 2*omega - 1 (mod 3)")
    print(f"  So '2/3 omega - 1' might mean '2 omega - 1 (mod 3)'")
    print(f"  In F_3[omega]: omega^2 = 2*omega + 2")

    # Verify in F_3
    print(f"\nIn F_3 = Z/3Z:")
    print(f"  omega^2 + omega + 1 = 0 mod 3")
    print(f"  omega^2 = -omega - 1 = 2*omega + 2 (mod 3)")
    print(f"  Note: 2 = -1 mod 3, so omega^2 = -omega + (-1) = -(omega+1)")

    print(f"\nConnection to I(omega) mod 3:")
    print(f"  I(omega) = 1 + a1*omega + a2*omega^2")
    print(f"           = 1 + a1*omega + a2*(2*omega + 2) mod 3")
    print(f"           = (1 + 2*a2) + (a1 + 2*a2)*omega mod 3")

    print(f"\n  Since a1 = total odd cycles, a2 = disjoint pairs:")
    print(f"  The real part (1 + 2*a2) mod 3 depends only on a2")
    print(f"  The omega part (a1 + 2*a2) mod 3 combines both")

    # What does this tell us about H mod 3?
    print(f"\nH mod 3 analysis:")
    print(f"  H = 1 + 2*a1 + 4*a2 = 1 + 2*a1 + a2 (mod 3)")
    print(f"  H mod 3 = (1 + 2*a1 + a2) mod 3")
    for a1_mod in range(3):
        for a2_mod in range(3):
            h_mod = (1 + 2*a1_mod + a2_mod) % 3
            print(f"  a1={a1_mod} mod 3, a2={a2_mod} mod 3 => H={h_mod} mod 3")

# ============================================================
# PART 5: omega^1 through omega^8 — pigeonhole
# ============================================================
def part5():
    print("\n" + "=" * 70)
    print("PART 5: OMEGA POWERS 1 THROUGH 8 — PIGEONHOLE AT DIM 2")
    print("=" * 70)

    omega = complex(-0.5, math.sqrt(3)/2)

    print(f"\nPowers of omega (period 3):")
    for k in range(9):
        w = omega ** k
        k_mod = k % 3
        label = ["1", "omega", "omega^2"][k_mod]
        print(f"  omega^{k} = {label:8s} = ({w.real:8.5f}, {w.imag:8.5f})")

    print(f"\nPigeonhole principle for Z[omega]:")
    print(f"  Z[omega] is a rank-2 Z-module: every element is a + b*omega")
    print(f"  Only 3 distinct powers: 1, omega, omega^2")
    print(f"  After 3 steps, powers cycle back")

    print(f"\nBut the user mentions 'omega^8' as where pigeonhole forces things")
    print(f"  omega^8 = omega^(3*2+2) = omega^2 = -omega - 1")
    print(f"  In terms of INDEX: 8 = 2*3 + 2, so omega^8 = omega^2")

    print(f"\nThe KEY insight: 8 = 2^3 = 2*3 + 2")
    print(f"  The exponent 8 encodes BOTH the exponential (2^3) and modular (2 mod 3) structure")
    print(f"  omega^{8} = omega^{8%3} = omega^2 = conjugate of omega")

    print(f"\nGeneral nth roots of unity and pigeonhole:")
    print(f"  For primitive nth root zeta_n, Z[zeta_n] has rank phi(n) as Z-module")
    print(f"  phi(3) = 2, phi(4) = 2, phi(5) = 4, phi(6) = 2, phi(7) = 6, phi(8) = 4")
    print(f"  At n=8: zeta_8 gives rank phi(8)=4 module")
    print(f"  8 powers to fill 4 dimensions => pigeonhole gives first FORCED dependence")
    print(f"  after 4+1=5 powers (5 = 2+3!)")

    print(f"\n  For omega (n=3): rank 2, period 3, pigeonhole at 2+1=3 powers")
    print(f"  For zeta_8 (n=8): rank 4, period 8, pigeonhole at 4+1=5 powers")
    print(f"  The pigeonhole bound for nth root is phi(n)+1 distinct powers before dependence")

    print(f"\n  phi(n)+1 for small n:")
    for n in range(2, 13):
        phi_n = sum(1 for k in range(1, n) if math.gcd(k, n) == 1)
        print(f"    n={n:2d}: phi({n})={phi_n}, pigeonhole at {phi_n+1} powers")

# ============================================================
# PART 6: Hilbert's 3rd problem
# ============================================================
def part6():
    print("\n" + "=" * 70)
    print("PART 6: HILBERT'S 3RD PROBLEM — DEHN INVARIANT & IRRATIONALITY")
    print("=" * 70)

    print(f"\nHilbert's 3rd problem (1900): Can any two polyhedra of equal volume")
    print(f"be cut into finitely many pieces and reassembled into each other?")
    print(f"\nAnswer: NO (Dehn, 1901)")
    print(f"The Dehn invariant is the obstruction.")

    print(f"\nDehn invariant: D(P) = sum over edges e of length(e) (x) dihedral_angle(e)")
    print(f"where (x) is tensor product in R (x)_Q (R/Q*pi)")
    print(f"Two polyhedra are scissors-congruent iff same volume AND same Dehn invariant.")

    print(f"\nThe IRRATIONALITY insertion:")
    print(f"  Regular tetrahedron: dihedral angle = arccos(1/3)")
    print(f"  arccos(1/3)/pi is IRRATIONAL")
    print(f"  Cube: all dihedral angles = pi/2 (rational multiple of pi)")
    print(f"  => Dehn(cube) = 0 but Dehn(regular tet) != 0")
    print(f"  => They are NOT scissors-congruent despite having equal volume")

    print(f"\nThe proof by contradiction pattern:")
    print(f"  1. Assume scissors-congruence: P = union of pieces = Q")
    print(f"  2. Dehn invariant is additive: D(P) = sum D(piece_i) = D(Q)")
    print(f"  3. But D(cube) = 0 and D(tet) != 0")
    print(f"  4. Contradiction.")

    print(f"\nConnection to tournaments:")
    print(f"  H(T) = I(Omega, 2) uses the EVALUATION at x=2")
    print(f"  I(Omega, omega) uses evaluation at x=omega (irrational)")
    print(f"  The Eisenstein norm N(I(omega)) is an INTEGER invariant")
    print(f"  Like the Dehn invariant, it provides an OBSTRUCTION:")
    print(f"  If two tournaments have different N(I(omega)), they are 'non-congruent'")

    dihedral_tet = math.acos(1/3)
    print(f"\n  arccos(1/3) = {dihedral_tet:.6f} rad = {math.degrees(dihedral_tet):.4f} deg")
    print(f"  arccos(1/3)/pi = {dihedral_tet/math.pi:.6f}... (irrational)")

    print(f"\nAlternating sum non-negativity:")
    print(f"  In the Dehn proof: the key is that contributions CANNOT cancel")
    print(f"  because the angle is irrational over Q*pi")
    print(f"  Similarly: I(Omega, x) has alternating-sign coefficients?")
    print(f"  I(x) = 1 + a1*x + a2*x^2 with a1, a2 >= 0 (no alternation needed)")
    print(f"  But I(x) at x=-1: I(-1) = 1 - a1 + a2 = alternating sum of alpha basis")
    print(f"  I(-1) >= 0 iff a2 >= a1 - 1")

    # Check I(-1) for random tournaments
    rng = np.random.default_rng(42)
    neg_count = 0
    total = 1000
    for _ in range(total):
        A = random_tournament(7, rng)
        a1, a2 = get_alpha_1_2(A)
        val = 1 - a1 + a2
        if val < 0:
            neg_count += 1
    print(f"\n  At n=7: I(-1) < 0 for {neg_count}/{total} tournaments ({100*neg_count/total:.1f}%)")
    print(f"  So alternating sum non-negativity FAILS (a1 can be large with a2=0)")

# ============================================================
# PART 7: Alternating sums in tournament context
# ============================================================
def part7():
    print("\n" + "=" * 70)
    print("PART 7: ALTERNATING SUMS AND NON-NEGATIVITY IN TOURNAMENTS")
    print("=" * 70)

    rng = np.random.default_rng(2026_0314)
    n = 7
    N = 500

    # Collect I(-1), I(-2), I(-3)
    neg_vals = defaultdict(list)
    for _ in range(N):
        A = random_tournament(n, rng)
        a1, a2 = get_alpha_1_2(A)
        for x_val in [-1, -2, -3]:
            v = I_poly(a1, a2, x_val)
            neg_vals[x_val].append(v)

    print(f"\nI(x) at negative integers (n=7, {N} samples):")
    for x_val in [-1, -2, -3]:
        vals = neg_vals[x_val]
        mn, mx = min(vals), max(vals)
        print(f"  I({x_val}): range [{mn}, {mx}], always >=0: {mn >= 0}")

    print(f"\nEuler characteristic perspective:")
    print(f"  I(-1) = chi(Omega) = Euler characteristic of independence complex")
    print(f"  I(-1) = 1 - a1 + a2 (alternating sum of face counts)")
    print(f"  For Omega(T): vertices = odd cycles, edges = disjoint pairs")
    print(f"  chi(Omega) = 1 - |V| + |E| = 1 - a1 + a2")

    print(f"\nWhen is chi(Omega) = 0?")
    print(f"  a2 = a1 - 1 => chi = 0")
    print(f"  This means the 'topological charge' of Omega vanishes")

    # Check what chi values occur
    chi_vals = defaultdict(int)
    for _ in range(2000):
        A = random_tournament(n, rng)
        a1, a2 = get_alpha_1_2(A)
        chi = 1 - a1 + a2
        chi_vals[chi] += 1

    print(f"\nchi(Omega(T)) distribution at n=7:")
    for c in sorted(chi_vals.keys()):
        print(f"  chi = {c:4d}: {chi_vals[c]:4d} ({100*chi_vals[c]/2000:.1f}%)")

# ============================================================
# PART 8: n=8 as tournament inflection point
# ============================================================
def part8():
    print("\n" + "=" * 70)
    print("PART 8: n=8 AS TOURNAMENT INFLECTION POINT (= 2^3)")
    print("=" * 70)

    print(f"\nThe number 8 = 2^3 is where multiple structural properties fail:")
    print()
    print(f"  HOMOLOGICAL FAILURES AT n=8:")
    print(f"  1. beta_3 > 1 first occurs (beta_3 = 2, freq 0.08%)")
    print(f"  2. Seesaw mechanism breaks (beta_3 and beta_4 coexist)")
    print(f"  3. i_*-injectivity fails (even with beta_4 = 0)")
    print(f"  4. beta_4 first appears (n=8 is onset)")
    print(f"  5. beta_5 also first appears at n=8")
    print()
    print(f"  GRAPH-THEORETIC FAILURES:")
    print(f"  6. Quasi-line property of Omega(T) fails at n=8 (49%)")
    print(f"  7. (Claw-free still holds at n=8, fails at n=9)")
    print()
    print(f"  NUMBER-THEORETIC:")
    print(f"  8. n=8 is the first EVEN n with self-complementary maximizer")
    print(f"     (H(SC max) = 661 at n=8)")
    print()
    print(f"  TOTAL: At least 7 independent structural breaks at n=8")

    print(f"\nWhy 8 = 2^3?")
    print(f"  The path homology chain complex has length n-1")
    print(f"  At n=8: chain complex has 7 = n-1 terms")
    print(f"  Betti numbers beta_0 through beta_7 are possible")
    print(f"  beta_2 = 0 always (proved), so effective dimension = 6")
    print(f"  With 6 slots and the seesaw constraint, n=8 is the first n")
    print(f"  where there are ENOUGH degrees of freedom for multiple nonzero betas")

    print(f"\n  Combinatorial: C(8,3) = 56 three-vertex subsets")
    print(f"  Each can be a 3-cycle (2 of 8 orientations)")
    print(f"  Max c3 at n=8: C(8,3)/4 = 14 (regular tournament)")
    print(f"  For comparison: C(7,3) = 35, max c3 = 35/5 = 7")

    print(f"\n  The ratio c3_max / C(n,3) = 1/(n-2):")
    for nn in range(3, 10):
        cmax = math.comb(nn, 3) // (nn - 2)
        ratio = 1 / (nn - 2)
        print(f"    n={nn}: c3_max = {cmax}, ratio = 1/{nn-2} = {ratio:.4f}")

# ============================================================
# PART 9: Cayley-Dickson tower
# ============================================================
def part9():
    print("\n" + "=" * 70)
    print("PART 9: CAYLEY-DICKSON TOWER: R -> C -> H -> O -> S")
    print("=" * 70)

    print(f"\nThe Cayley-Dickson construction doubles dimension at each step:")
    print(f"  R (dim 1):  Real numbers — all properties")
    print(f"  C (dim 2):  Complex numbers — lose ordering")
    print(f"  H (dim 4):  Quaternions — lose COMMUTATIVITY")
    print(f"  O (dim 8):  Octonions — lose ASSOCIATIVITY")
    print(f"  S (dim 16): Sedenions — lose ALTERNATIVITY (and division)")

    print(f"\nDimensions: 1, 2, 4, 8, 16 = 2^0, 2^1, 2^2, 2^3, 2^4")

    print(f"\nProperties lost at each step:")
    tower = [
        ("R", 1, "2^0", "ordered, commutative, associative, alternative, division, composition"),
        ("C", 2, "2^1", "commutative, associative, alternative, division, composition"),
        ("H", 4, "2^2", "associative, alternative, division, composition"),
        ("O", 8, "2^3", "alternative, division, composition"),
        ("S", 16, "2^4", "power-associative, flexible (but NOT division, NOT composition)"),
    ]
    for name, dim, power, props in tower:
        print(f"  {name} (dim {dim:2d} = {power}): {props}")

    print(f"\nProperty loss at dim 2^k:")
    print(f"  2^0 -> 2^1: lose total ordering        (at dim 2)")
    print(f"  2^1 -> 2^2: lose commutativity          (at dim 4)")
    print(f"  2^2 -> 2^3: lose associativity           (at dim 8)")
    print(f"  2^3 -> 2^4: lose alternativity & division (at dim 16)")

    print(f"\nHurwitz theorem: normed division algebras exist ONLY at dim 1, 2, 4, 8")
    print(f"  Equivalently: composition algebras over R = R, C, H, O")
    print(f"  Sum of squares identity: |xy|^2 = |x|^2 * |y|^2")
    print(f"  This corresponds to: 1, 2, 4, 8 squares identities")
    print(f"  (Hurwitz 1, 2, 4, 8 squares theorem)")

    print(f"\n  Topological manifestation: parallelizable spheres only S^0, S^1, S^3, S^7")
    print(f"  Dimensions 0, 1, 3, 7 = (dim of algebra) - 1")
    print(f"  These are exactly 2^k - 1 for k = 0, 1, 2, 3")

# ============================================================
# PART 10: Property loss parallels in tournaments
# ============================================================
def part10():
    print("\n" + "=" * 70)
    print("PART 10: CAYLEY-DICKSON PROPERTY LOSS PARALLELS IN TOURNAMENTS")
    print("=" * 70)

    print(f"\nTournament property losses at powers of 2:")
    print(f"  n=1 (2^0): Trivial tournament. All properties hold vacuously.")
    print(f"  n=2 (2^1): Unique tournament. H=1. Only beta_0=1. No cycles.")
    print(f"  n=4 (2^2): First n with non-trivial structure. C(4,3)=4 triples.")
    print(f"              First n where score sequence matters.")
    print(f"              Omega(T) can be empty or non-empty.")
    print(f"  n=8 (2^3): THE INFLECTION POINT")
    print(f"              - beta_3 can exceed 1 (first time)")
    print(f"              - Seesaw mechanism fails")
    print(f"              - i_*-injectivity fails")
    print(f"              - Quasi-line fails for Omega(T)")
    print(f"  n=16(2^4): Predicted: more properties fail")
    print(f"              - Independence polynomial coefficients more complex")
    print(f"              - Higher Betti numbers emerge")

    print(f"\nThe parallel structure:")
    print(f"  Cayley-Dickson at dim 2^k: lose k-th algebraic property")
    print(f"  Tournament at n=2^k: lose k-th homological property")
    print(f"  k=0: trivial")
    print(f"  k=1: non-trivial but simple")
    print(f"  k=2: first interesting structure")
    print(f"  k=3: MAJOR phase transition (octonions / seesaw break)")
    print(f"  k=4: further degradation (sedenions / ???)")

    print(f"\nSpecifically for k=3 (the CRITICAL step):")
    print(f"  Octonions (dim 8): lose ASSOCIATIVITY of multiplication")
    print(f"  Tournaments (n=8): lose 'ASSOCIATIVITY' of homological seesaw")
    print(f"    The seesaw beta_1*beta_3=0 is a kind of associativity:")
    print(f"    the boundary maps d_2 and d_3 cannot simultaneously have")
    print(f"    nontrivial kernel AND nontrivial cokernel")
    print(f"    At n=8, this constraint breaks: both can be nontrivial")

    print(f"\nAlternative connection via Hurwitz:")
    print(f"  Hurwitz: sum-of-squares identity in n variables exists iff n in {{1,2,4,8}}")
    print(f"  Tournament: exactly at n=8, the 'sum-of-squares' structure of")
    print(f"  the Eisenstein norm breaks in a new way (higher Betti numbers")
    print(f"  allow more complex norm values)")

# ============================================================
# PART 11: I(8) — the octonionic threshold
# ============================================================
def part11():
    print("\n" + "=" * 70)
    print("PART 11: I(8) = 2^3 — THE OCTONIONIC EVALUATION")
    print("=" * 70)

    print(f"\nI(x) = 1 + a1*x + a2*x^2 evaluated at x = 8 = 2^3:")
    print(f"  I(8) = 1 + 8*a1 + 64*a2")
    print(f"  I(2) = 1 + 2*a1 + 4*a2   = H(T)")
    print(f"  I(3) = 1 + 3*a1 + 9*a2")

    print(f"\nRatio I(8)/I(2) for random n=7 tournaments:")
    rng = np.random.default_rng(2026)
    ratios = []
    for _ in range(500):
        A = random_tournament(7, rng)
        a1, a2 = get_alpha_1_2(A)
        h = I_poly(a1, a2, 2)
        i8 = I_poly(a1, a2, 8)
        if h > 0:
            ratios.append(i8 / h)

    print(f"  Mean I(8)/H = {np.mean(ratios):.2f}")
    print(f"  Range: [{min(ratios):.2f}, {max(ratios):.2f}]")

    print(f"\nI(8) = 1 + 8*a1 + 64*a2")
    print(f"     = I(2) + 6*a1 + 60*a2")
    print(f"     = H + 6*(a1 + 10*a2)")
    print(f"  So I(8) - H = 6*(a1 + 10*a2)")
    print(f"  The gap is always divisible by 6 = 2*3!")

    print(f"\nI(8) mod small numbers:")
    print(f"  I(8) mod 2 = 1 (always odd, like H)")
    print(f"  I(8) mod 3 = 1 + 2*a1 + a2 (mod 3) = H mod 3")
    print(f"  I(8) mod 7 = 1 + a1 + a2 (mod 7) = I(1) mod 7 = (1+a1+a2) mod 7")

    print(f"\n  I(8) mod 7 = I(1) mod 7  (since 8 = 1 mod 7)")
    print(f"  I(1) = 1 + a1 + a2 = number of independent sets + 1")
    print(f"  This is the TOTAL count of independent sets in Omega(T)!")

    print(f"\nI(8) mod 8:")
    print(f"  I(8) = 1 + 8*a1 + 64*a2 = 1 mod 8")
    print(f"  I(8) mod 8 = 1 ALWAYS (universal!)")
    print(f"  More generally: I(b) mod b = 1 for all b (CRT tower, HYP-963)")

    # Verify I(8) = 1 mod 8
    print(f"\nVerification I(8) mod 8 = 1:")
    rng = np.random.default_rng(42)
    all_one = True
    for _ in range(200):
        A = random_tournament(7, rng)
        a1, a2 = get_alpha_1_2(A)
        i8 = I_poly(a1, a2, 8)
        if i8 % 8 != 1:
            all_one = False
            print(f"  COUNTEREXAMPLE: a1={a1}, a2={a2}, I(8)={i8}, mod 8 = {i8%8}")
    if all_one:
        print(f"  CONFIRMED: I(8) = 1 mod 8 for all 200 tested tournaments")

# ============================================================
# PART 12: Exponentiation in H(T)
# ============================================================
def part12():
    print("\n" + "=" * 70)
    print("PART 12: EXPONENTIATION IN H(T) — FRACTIONAL POWERS")
    print("=" * 70)

    print(f"\nH = I(2) = 1 + 2*a1 + 4*a2")
    print(f"Can we make sense of 'H to a power'?")

    rng = np.random.default_rng(999)
    print(f"\nH^(1/2) and H^(1/3) for n=7 tournaments:")
    for trial in range(10):
        A = random_tournament(7, rng)
        a1, a2 = get_alpha_1_2(A)
        h = I_poly(a1, a2, 2)
        print(f"  H={h:4d}, H^(1/2)={h**0.5:8.4f}, H^(1/3)={h**(1/3):8.4f}, "
              f"a1={a1:3d}, a2={a2:3d}")

    print(f"\nH^(1/3) interpretation via Claim B:")
    print(f"  H = sum over odd-cycle collections 2^|C|")
    print(f"  = sum_{{k=0}}^{{alpha(Omega)}} f_k * 2^k  where f_k = independent sets of size k")
    print(f"  For n=7: H = 1 + 2*a1 + 4*a2 (only k=0,1,2)")
    print(f"  H^(1/3) would mix the binary layers in a non-integer way")

    print(f"\nMore natural: I(2^(1/3)) = 1 + 2^(1/3)*a1 + 2^(2/3)*a2")
    print(f"  This is NOT an integer, but lives in Z[2^(1/3)]")
    print(f"  The ring Z[2^(1/3)] is a degree-3 extension of Z")
    print(f"  Minimal polynomial: x^3 - 2 = 0 (Eisenstein at p=2)")

    # I(cube_root_2)
    cbrt2 = 2 ** (1/3)
    print(f"\n  2^(1/3) = {cbrt2:.6f}")
    for trial in range(5):
        A = random_tournament(7, rng)
        a1, a2 = get_alpha_1_2(A)
        val = 1 + cbrt2 * a1 + cbrt2**2 * a2
        print(f"  a1={a1:3d}, a2={a2:3d}: I(2^(1/3)) = {val:.4f}")

    print(f"\nThe NUMBER FIELD tower:")
    print(f"  I(2) lives in Z")
    print(f"  I(omega) lives in Z[omega] = Z[(-1+sqrt(-3))/2]")
    print(f"  I(2^(1/3)) lives in Z[2^(1/3)]")
    print(f"  These are three distinct quadratic/cubic extensions of Q")
    print(f"  Union: I lives in Q(omega, 2^(1/3)) which has degree 6 = 2*3 over Q!")
    print(f"  The splitting field of x^3 - 2 is exactly Q(omega, 2^(1/3))")
    print(f"  Galois group: S_3 (order 6 = 2*3 = |S_3|)")

# ============================================================
# PART 13: The 8-9 pair
# ============================================================
def part13():
    print("\n" + "=" * 70)
    print("PART 13: THE 8-9 PAIR — CONSECUTIVE POWERS AND NEAR-MISSES")
    print("=" * 70)

    print(f"\n8 = 2^3 and 9 = 3^2: the ONLY consecutive perfect powers (Catalan/Mihailescu)")

    print(f"\nIn tournament theory:")
    print(f"  n=8: homological inflection (beta_3>1, seesaw breaks)")
    print(f"  n=9: I.P. real-roots first fail (THM-025 counterexample)")
    print(f"  The failures at 8 and 9 are DIFFERENT IN KIND:")
    print(f"    n=8: topological (path homology)")
    print(f"    n=9: algebraic (independence polynomial roots)")

    print(f"\nH values at n=8 and n=9:")
    rng = np.random.default_rng(2026_0314)

    # Small sample at n=7 (feasible)
    print(f"\n  Sample H values at n=7:")
    h_vals_7 = []
    for _ in range(200):
        A = random_tournament(7, rng)
        H = count_ham_paths(A)
        h_vals_7.append(H)
    print(f"    Min H = {min(h_vals_7)}, Max H = {max(h_vals_7)}, Mean = {np.mean(h_vals_7):.1f}")

    print(f"\n  H=8 is impossible (H always odd)")
    print(f"  H=9 = 1 + 2*4 + 4*0 (a1=4, a2=0) — achievable at n>=5")

    print(f"\nThe Catalan equation x^p - y^q = 1:")
    print(f"  3^2 - 2^3 = 1")
    print(f"  Tournament version: I(3) - I(2) = step = a1 + 5*a2")
    print(f"  When does I(3) - I(2) = 1?")
    print(f"    a1 + 5*a2 = 1 => a1=1, a2=0")
    print(f"    This means exactly 1 odd cycle, 0 disjoint pairs")
    print(f"    H = I(2) = 1 + 2 = 3, I(3) = 1 + 3 = 4")
    print(f"    So the 'Catalan' step I(3)-I(2)=1 forces H=3 (smallest non-trivial)")

    print(f"\n  When does I(3)^2 - I(2)^3 = 1? (actual Catalan structure)")
    for a1 in range(0, 40):
        for a2 in range(0, 15):
            i2 = 1 + 2*a1 + 4*a2
            i3 = 1 + 3*a1 + 9*a2
            if i3**2 - i2**3 == 1:
                print(f"    SOLUTION: a1={a1}, a2={a2}: I(2)={i2}, I(3)={i3}")
                print(f"      I(3)^2 = {i3**2}, I(2)^3 = {i2**3}, diff = 1")

# ============================================================
# PART 14: Synthesis
# ============================================================
def part14():
    print("\n" + "=" * 70)
    print("PART 14: SYNTHESIS — THE TRINITY OF OPERATIONS ON {2,3}")
    print("=" * 70)

    print(f"\nThe three binary operations on the pair (2, 3):")
    print(f"")
    print(f"  ADDITION:        2 + 3 = 5    (prime, commutative)")
    print(f"  MULTIPLICATION:  2 * 3 = 6    (composite = 2*3, commutative)")
    print(f"  EXPONENTIATION:  2 ^ 3 = 8    (= 2^3, NON-commutative)")
    print(f"                   3 ^ 2 = 9    (= 3^2)")
    print(f"")
    print(f"  Outputs: {{5, 6, 8, 9}} (four numbers from three operations)")

    print(f"\nTournament theory meaning:")
    print(f"  5 = 2+3: first n with non-trivial Omega(T); I(Omega,x) quadratic")
    print(f"  6 = 2*3: first n with beta_3 > 0; multiplicative structure emerges")
    print(f"  8 = 2^3: inflection point; homological properties fail")
    print(f"  9 = 3^2: algebraic properties fail (real roots); claw-free fails")

    print(f"\nThe ORDERED structure:")
    print(f"  5 < 6 < 8 < 9")
    print(f"  + < * < ^(2,3) < ^(3,2)")
    print(f"  Each operation produces a LARGER result on (2,3)")
    print(f"  The gap pattern: 6-5=1, 8-6=2, 9-8=1")
    print(f"  Symmetric: gaps are 1, 2, 1 (palindromic!)")

    print(f"\nCyclotomic connections:")
    print(f"  Phi_1(2) = 1    (2-1)")
    print(f"  Phi_2(2) = 3    (2+1)")
    print(f"  Phi_3(2) = 7    (2^2+2+1)")
    print(f"  Phi_4(2) = 5    (2^2+1)")
    print(f"  Phi_5(2) = 31   (2^4+2^3+2^2+2+1)")
    print(f"  Phi_6(2) = 3    (2^2-2+1)")
    print(f"  Phi_8(2) = 17   (2^4+1)")
    print(f"  Phi_9(2) = 73   (2^6+2^3+1)")

    print(f"\n  Products:")
    print(f"  2^1 - 1 = 1  = Phi_1(2)")
    print(f"  2^2 - 1 = 3  = Phi_1(2)*Phi_2(2) = 1*3")
    print(f"  2^3 - 1 = 7  = Phi_1(2)*Phi_3(2) = 1*7")
    print(f"  2^6 - 1 = 63 = Phi_1*Phi_2*Phi_3*Phi_6 = 1*3*7*3 = 63")
    print(f"  Note: 63 = 2^6 - 1 is in our forbidden set!")
    print(f"  2^6 = 64, and 63 = forbidden H, and 6 = 2*3")
    print(f"  So 2^(2*3) - 1 is forbidden!")

    print(f"\nThe NUMBER FIELD synthesis:")
    print(f"  Q(omega, 2^(1/3)) is the splitting field of x^3 - 2")
    print(f"  Degree [Q(omega, 2^(1/3)) : Q] = 6 = 2*3")
    print(f"  Galois group = S_3, order 6 = 3! = 2*3")
    print(f"  The tournament independence polynomial I(x) = 1 + a1*x + a2*x^2")
    print(f"  has degree 2 (at n=7), so its splitting field is at most degree 2 over Q")
    print(f"  But evaluating I at 2^(1/3) yields elements of the DEGREE 6 field")
    print(f"  Combining all tournament evaluations requires the full S_3 Galois group")

    print(f"\nFinal observation: FIVE numbers summarize everything")
    print(f"  2, 3, 5=2+3, 6=2*3, 7=Phi_3(2)=2^2+2+1")
    print(f"  These are the first 5 primes: 2, 3, 5, 7 (plus 6=2*3)")
    print(f"  And n=7 is where our quadratic I.P. theory is most complete")
    print(f"  At n=7: alpha(Omega) = 2 always, so I(x) is quadratic")
    print(f"  The theory 'lives' in the quadratic world of 2 and 3")

    print(f"\n  Beyond n=7: I(x) can be cubic or higher")
    print(f"  The quadratic era ends, the cubic era begins")
    print(f"  Exponentiation (2^3=8) marks this transition")
    print(f"  Just as octonions (dim 2^3=8) end the division algebra era")

    print(f"\n{'='*70}")
    print(f"END OF EXPLORATION")
    print(f"{'='*70}")

def main():
    part1()
    part2()
    part3()
    part4()
    part5()
    part6()
    part7()
    part8()
    part9()
    part10()
    part11()
    part12()
    part13()
    part14()

if __name__ == "__main__":
    main()
