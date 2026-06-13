"""
hilbert3_dehn_tournament.py — Deep dive into Hilbert's 3rd problem analogy
and irrationality insertion in tournament parity theory.
kind-pasteur-2026-03-14-S65

The Dehn invariant proof pattern:
  1. Define an additive invariant that lives in a tensor product
  2. Show it separates objects that share a simpler invariant (volume)
  3. The KEY is irrationality: arccos(1/3)/pi is irrational

Tournament parallel:
  1. I(omega) lives in Z[omega], a richer ring than Z
  2. The Eisenstein norm N(I(omega)) separates tournaments with same H
  3. The key is that omega is irrational over Q

PART 1:  Dehn invariant structure: R tensor_Q (R/Q*pi)
PART 2:  Tournament tensor: I(x) as element of Z[x]/(x^2+x+1) tensor Z
PART 3:  Irrationality insertion: s = smallest with omega^s = 1 (s=3)
PART 4:  The proof-by-contradiction template applied to H=21
PART 5:  (omega^a - omega^b)^c in Z[omega] — algebraic closure structure
PART 6:  Alternating sum non-negativity revisited: Dehn vs tournament
PART 7:  The 3-periodicity of omega and the 3-cycle periodicity of tournaments
PART 8:  Norm filtration: tournaments by Eisenstein norm level
PART 9:  Dehn-like obstruction for H achievability
PART 10: Connection to Hilbert's 3rd: scissors congruence of tournaments
"""

import numpy as np
from itertools import combinations
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
    cycles = get_all_odd_cycles(A)
    a1 = len(cycles)
    a2 = 0
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            if cycles[i].isdisjoint(cycles[j]):
                a2 += 1
    return a1, a2

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

def eisenstein_norm(a1, a2):
    """N(I(omega)) = a1^2 - a1*a2 + a2^2 - a1 - a2 + 1"""
    return a1**2 - a1*a2 + a2**2 - a1 - a2 + 1

# ============================================================
# PART 1: Dehn invariant structure
# ============================================================
def part1():
    print("=" * 70)
    print("PART 1: DEHN INVARIANT STRUCTURE")
    print("=" * 70)

    print(f"\nThe Dehn invariant lives in:")
    print(f"  R tensor_Q (R / Q*pi)")
    print(f"  = R tensor_Q V   where V = R / Q*pi")
    print(f"")
    print(f"This is a Q-vector space (uncountable dimension).")
    print(f"Key property: theta in V is ZERO iff theta/pi is RATIONAL.")
    print(f"")
    print(f"For a polyhedron P:")
    print(f"  D(P) = sum_edges  length(e) tensor angle(e)  in  R tensor V")
    print(f"")
    print(f"The PROOF that cube != regular tetrahedron:")
    print(f"  D(cube) = 12 * L * (pi/2 mod Q*pi) = 12 * L * 0 = 0")
    print(f"  D(tet) = 6 * L' * (arccos(1/3) mod Q*pi) != 0")
    print(f"  because arccos(1/3)/pi is IRRATIONAL")
    print(f"")
    print(f"The irrationality is what makes the invariant NONTRIVIAL.")
    print(f"Without it, D would be identically zero.")

    # Niven's theorem: arccos(r) / pi is rational iff r in {0, +-1/2, +-1}
    print(f"\nNiven's theorem: arccos(r)/pi is rational iff r in {{0, +-1/2, +-1}}")
    print(f"  arccos(0) = pi/2       rational")
    print(f"  arccos(1/2) = pi/3     rational")
    print(f"  arccos(1/3) = ?        IRRATIONAL")
    print(f"  arccos(1/4) = ?        IRRATIONAL")
    print(f"")
    print(f"So 1/3 is the FIRST reciprocal integer giving an irrational angle.")
    print(f"And 3 is the NUMBER controlling this irrationality!")

# ============================================================
# PART 2: Tournament tensor structure
# ============================================================
def part2():
    print("\n" + "=" * 70)
    print("PART 2: TOURNAMENT TENSOR STRUCTURE")
    print("=" * 70)

    print(f"\nI(x) = 1 + a1*x + a2*x^2 lives in Z[x].")
    print(f"Evaluating at omega: I(omega) in Z[omega] = Z[x]/(x^2+x+1)")
    print(f"")
    print(f"The 'tournament Dehn invariant' is the omega-evaluation:")
    print(f"  D_T = I(omega) = (1-a2) + (a1-a2)*omega  in Z[omega]")
    print(f"")
    print(f"Two tournaments T1, T2 with same H but different D_T are")
    print(f"'not scissors-congruent' in the omega sense.")

    # Find examples at n=7
    rng = np.random.default_rng(2026_0314)
    h_to_dt = defaultdict(set)
    n = 7
    for _ in range(2000):
        A = random_tournament(n, rng)
        a1, a2 = get_alpha_1_2(A)
        H = 1 + 2*a1 + 4*a2
        real_part = 1 - a2
        omega_part = a1 - a2
        h_to_dt[H].add((real_part, omega_part))

    print(f"\nH values with multiple D_T (first 10):")
    count = 0
    for H in sorted(h_to_dt.keys()):
        dts = h_to_dt[H]
        if len(dts) > 1:
            count += 1
            if count <= 10:
                norms = [r**2 - r*w + w**2 for r,w in dts]
                print(f"  H={H:3d}: D_T in {{{', '.join(f'({r}+{w}w)' for r,w in sorted(dts))}}}")
                print(f"          Norms: {sorted(norms)}")

    print(f"\n  Total H values with multiple omega-types: {count}")
    unique_h = sum(1 for dts in h_to_dt.values() if len(dts) == 1)
    print(f"  H values with unique omega-type: {unique_h}")

    print(f"\nThe omega-evaluation REFINES the H-classification:")
    print(f"  H = I(2) loses information about how a1 and a2 combine")
    print(f"  D_T = I(omega) retains (a1 mod 3, a2 mod 3) information")
    print(f"  Together, (H, D_T) almost determines (a1, a2)")

# ============================================================
# PART 3: Irrationality insertion — s = 3
# ============================================================
def part3():
    print("\n" + "=" * 70)
    print("PART 3: IRRATIONALITY INSERTION — s = 3 (order of omega)")
    print("=" * 70)

    print(f"\nLet s be the smallest positive integer with omega^s = 1.")
    print(f"For omega = e^(2pi*i/3): s = 3.")
    print(f"")
    print(f"This is the 'irrationality insertion' step:")
    print(f"  omega is NOT rational (omega = -1/2 + sqrt(3)/2 * i)")
    print(f"  But omega^3 = 1 (it has FINITE ORDER)")
    print(f"  So omega generates a CYCLOTOMIC field Q(omega)/Q of degree phi(3)=2")
    print(f"")
    print(f"The Dehn parallel:")
    print(f"  Dehn uses angles theta with theta/pi IRRATIONAL")
    print(f"  => theta generates an infinite-dimensional quotient in R/Q*pi")
    print(f"  Tournament uses omega with omega^3 = 1 (RATIONAL relation)")
    print(f"  => omega generates a FINITE (degree 2) extension of Q")
    print(f"")
    print(f"KEY DIFFERENCE: Dehn's irrationality is 'harder' (transcendence degree)")
    print(f"Tournament's omega is algebraic (degree 2)")
    print(f"But BOTH provide non-trivial invariants because the EVALUATION POINT")
    print(f"is not an integer.")

    print(f"\nWhat if we used zeta_s for other s?")
    for s in [2, 3, 4, 5, 6, 7, 8]:
        if s == 1:
            continue
        zeta = complex(math.cos(2*math.pi/s), math.sin(2*math.pi/s))
        phi_s = sum(1 for k in range(1, s) if math.gcd(k, s) == 1)
        # Evaluate Phi_s(2)
        # Phi_s(2) = product_{k coprime to s} (2 - zeta_s^k)
        phi_s_at_2 = 1
        for k in range(1, s):
            if math.gcd(k, s) == 1:
                root = complex(math.cos(2*math.pi*k/s), math.sin(2*math.pi*k/s))
                phi_s_at_2 *= (2 - root)
        phi_val = round(phi_s_at_2.real)
        print(f"  s={s}: zeta_{s}, phi({s})={phi_s}, Phi_{s}(2) = {phi_val}")

    print(f"\nThe s=3 choice is special because:")
    print(f"  1. phi(3)=2 = minimal non-trivial extension (matches I.P. degree at n<=7)")
    print(f"  2. Phi_3(2)=7 = n for our complete theory (at n=7, alpha(Omega)<=2 always)")
    print(f"  3. 3 is the length of the fundamental cycle (3-cycles)")
    print(f"  4. omega^2+omega+1=0 mirrors I(x) = 1+a1*x+a2*x^2 at a1=a2=1")
    print(f"  5. The Eisenstein integers Z[omega] are a PID (Euclidean domain)")

# ============================================================
# PART 4: Proof-by-contradiction for H=21
# ============================================================
def part4():
    print("\n" + "=" * 70)
    print("PART 4: PROOF-BY-CONTRADICTION TEMPLATE FOR H=21")
    print("=" * 70)

    print(f"\nH = 21 = 3 * 7 = Phi_2(2) * Phi_3(2)")
    print(f"Claim: H=21 is impossible at n=7.")
    print(f"")
    print(f"Step 1: If H=21, then 1+2*a1+4*a2 = 21, so a1+2*a2 = 10.")
    print(f"  Possible (a1, a2): (10,0), (8,1), (6,2), (4,3), (2,4), (0,5)")

    print(f"\nStep 2: Compute Eisenstein norm for each candidate:")
    for a1, a2 in [(10,0), (8,1), (6,2), (4,3), (2,4), (0,5)]:
        N = eisenstein_norm(a1, a2)
        print(f"  (a1={a1:2d}, a2={a2}): N = {N:3d}, N mod 3 = {N%3}")

    print(f"\nStep 3: Check which (a1, a2) are achievable at n=7.")
    print(f"  From forbidden_boundary.py analysis:")
    print(f"  - (10, 0): a1=10 with a2=0 is IMPOSSIBLE (gap in achievable a1)")
    print(f"  - (8, 1): need 8 odd cycles with exactly 1 disjoint pair")
    print(f"  - (6, 2): need 6 odd cycles with exactly 2 disjoint pairs")
    print(f"  - (4, 3): need 4 odd cycles with 3 disjoint pairs")
    print(f"  - (2, 4): need 2 odd cycles with 4 disjoint pairs (impossible: C(2,2)=1)")
    print(f"  - (0, 5): need 0 odd cycles with 5 pairs (impossible: 0 cycles => 0 pairs)")

    print(f"\nStep 4: Eliminate remaining candidates")
    print(f"  (2, 4): IMPOSSIBLE — only 2 cycles, max C(2,2)=1 disjoint pair")
    print(f"  (0, 5): IMPOSSIBLE — 0 cycles means 0 pairs")
    print(f"  (4, 3): IMPOSSIBLE — 4 cycles, need C(4,2)=6 pairs, only 3 can be disjoint")
    print(f"          Actually: 3 disjoint pairs from 4 cycles requires each pair to use")
    print(f"          disjoint vertex sets, but 3 pairs = 6 cycles (contradiction: only 4)")

    # More careful: 3 disjoint pairs from 4 cycles
    # A pair is 2 cycles. 3 pairs use 6 cycle-slots from 4 cycles.
    # Each cycle can appear in at most... well, pairs are just edges in Omega graph
    # 4 vertices, 3 edges possible from K_4 minus matching, etc.
    print(f"\n  Wait: a2 = number of edges in Omega(T). With 4 vertices, max edges = C(4,2) = 6.")
    print(f"  So a2=3 with a1=4 IS possible graph-theoretically (triangle on 4 vertices).")
    print(f"  But do these 4 cycles actually exist with 3 disjoint pairs at n=7?")

    # Check computationally
    rng = np.random.default_rng(2026_0314)
    target_found = {(10,0): 0, (8,1): 0, (6,2): 0, (4,3): 0}
    total = 5000
    for _ in range(total):
        A = random_tournament(7, rng)
        a1, a2 = get_alpha_1_2(A)
        key = (a1, a2)
        if key in target_found:
            target_found[key] += 1

    print(f"\n  Computational search ({total} samples at n=7):")
    for key in [(10,0), (8,1), (6,2), (4,3)]:
        print(f"    (a1={key[0]}, a2={key[1]}): found {target_found[key]} times")

    print(f"\nStep 5: The Dehn-style conclusion")
    print(f"  If NONE of the (a1, a2) pairs achieving H=21 are realizable,")
    print(f"  then H=21 is impossible by EXHAUSTION OF THE FIBER.")
    print(f"  This is analogous to Dehn's proof: we show all 'decompositions'")
    print(f"  of H=21 into (a1, a2) components are individually impossible.")

# ============================================================
# PART 5: (omega^a - omega^b)^c in Z[omega]
# ============================================================
def part5():
    print("\n" + "=" * 70)
    print("PART 5: (omega^a - omega^b)^c IN Z[omega]")
    print("=" * 70)

    omega = complex(-0.5, math.sqrt(3)/2)

    print(f"\nAll differences omega^a - omega^b for a,b in {{0,1,2}}:")
    for a in range(3):
        for b in range(3):
            if a == b:
                continue
            diff = omega**a - omega**b
            # Express in Z[omega] basis: find r, s such that diff = r + s*omega
            # omega^0 = 1, omega^1 = omega, omega^2 = -1-omega
            vals = {0: (1, 0), 1: (0, 1), 2: (-1, -1)}
            ra, sa = vals[a]
            rb, sb = vals[b]
            r = ra - rb
            s = sa - sb
            norm = r**2 - r*s + s**2
            print(f"  omega^{a} - omega^{b} = {r} + {s}*omega, "
                  f"norm = {norm}, = sqrt(3)*e^(i*theta)")

    print(f"\nPowers (omega^a - omega^b)^c:")
    for a, b in [(1, 0), (0, 1), (2, 0)]:
        vals = {0: (1, 0), 1: (0, 1), 2: (-1, -1)}
        ra, sa = vals[a]
        rb, sb = vals[b]
        r = ra - rb
        s = sa - sb

        print(f"\n  Base: omega^{a} - omega^{b} = {r} + {s}*omega")
        # Compute powers by repeated multiplication in Z[omega]
        pr, ps = 1, 0  # start with 1
        for c in range(1, 9):
            # (pr + ps*omega) * (r + s*omega)
            # = pr*r + (pr*s + ps*r)*omega + ps*s*omega^2
            # = pr*r - ps*s + (pr*s + ps*r - ps*s)*omega
            new_r = pr*r - ps*s
            new_s = pr*s + ps*r - ps*s
            pr, ps = new_r, new_s
            norm = pr**2 - pr*ps + ps**2
            print(f"    c={c}: ({r}+{s}w)^{c} = {pr} + {ps}*omega, norm = {norm} = {r**2-r*s+s**2}^{c}")

    print(f"\nKey observation: |omega^a - omega^b|^2 = 3 for all a != b (mod 3)")
    print(f"  So norm of (omega^a - omega^b) = 3")
    print(f"  And norm of (omega^a - omega^b)^c = 3^c")
    print(f"")
    print(f"  For c=3: (omega^a - omega^b)^3 has norm 27 = 3^3")
    print(f"  These elements 'saturate' the norm-27 level of Z[omega]")
    print(f"  In tournament terms: I(omega) can never have norm 3^c for certain c")
    print(f"  because the (a1, a2) constraints prevent it")

# ============================================================
# PART 6: Alternating sum non-negativity revisited
# ============================================================
def part6():
    print("\n" + "=" * 70)
    print("PART 6: ALTERNATING SUM NON-NEGATIVITY — DEHN VS TOURNAMENT")
    print("=" * 70)

    print(f"\nDehn invariant: sum length(e) tensor angle(e)")
    print(f"  Non-negativity: each term is a non-negative length times an angle")
    print(f"  The ADDITIVITY is key: cutting a polyhedron, new internal edges")
    print(f"  have paired dihedral angles summing to pi, which is 0 in R/Q*pi")
    print(f"  So internal contributions cancel, boundary is invariant")

    print(f"\nTournament I(x): 1 + a1*x + a2*x^2")
    print(f"  All coefficients (1, a1, a2) are NON-NEGATIVE")
    print(f"  For x >= 0: I(x) >= 1 (always positive)")
    print(f"  For x = -1: I(-1) = 1 - a1 + a2 (alternating, can be negative)")

    print(f"\nThe Dehn proof uses non-negativity of LENGTHS, not angles.")
    print(f"The tournament version:")
    print(f"  'Lengths' = face counts f_k = (1, a1, a2)")
    print(f"  'Angles' = powers of evaluation point x^k = (1, x, x^2)")
    print(f"  At x=2: all 'angles' positive, so I(2) >= 1 (trivially)")
    print(f"  At x=omega: 'angles' are complex, but norm is always >=0")

    print(f"\nThe analog of Dehn's invariant being well-defined:")
    print(f"  Adding/removing a 'cut' in tournament = arc flip (single edge reversal)")
    print(f"  Arc flip changes H by H(T/e) - H(T'/e') (deletion-contraction)")
    print(f"  This delta is NOT zero in general (tournaments are not 'scissors-congruent')")
    print(f"  Each arc flip is an IRREVERSIBLE modification (unlike cutting a polyhedron)")

    # What's the analog of "internal edges cancel"?
    print(f"\nThe deeper analog: COMPLEMENTATION")
    print(f"  For self-complementary tournaments: T and T^c have same score sequence")
    print(f"  The 'complement' operation is like 'reflecting' a polyhedron")
    print(f"  H(T) = H(T^c) for SC tournaments (by definition)")
    print(f"  But I(omega) might differ: I_T(omega) vs I_T^c(omega)")

    # Check I(omega) for T vs T^c
    rng = np.random.default_rng(42)
    print(f"\n  Checking I(omega) for T vs complement T^c at n=5:")
    for trial in range(5):
        A = random_tournament(5, rng)
        Ac = 1 - A
        np.fill_diagonal(Ac, 0)
        a1, a2 = get_alpha_1_2(A)
        a1c, a2c = get_alpha_1_2(Ac)
        h = 1 + 2*a1 + 4*a2
        hc = 1 + 2*a1c + 4*a2c
        r = 1 - a2
        rc = 1 - a2c
        w = a1 - a2
        wc = a1c - a2c
        print(f"    T: (a1={a1:2d},a2={a2}), H={h:3d}, I(w)=({r}+{w}w) | "
              f"T^c: (a1c={a1c:2d},a2c={a2c}), Hc={hc:3d}, I(w)=({rc}+{wc}w)")

# ============================================================
# PART 7: 3-periodicity of omega and 3-cycle periodicity
# ============================================================
def part7():
    print("\n" + "=" * 70)
    print("PART 7: 3-PERIODICITY OF OMEGA AND 3-CYCLE PERIODICITY")
    print("=" * 70)

    print(f"\nomega^3 = 1: the period is 3")
    print(f"3-cycles are the fundamental building blocks of tournaments")
    print(f"")
    print(f"Connection: the cyclotomic polynomial Phi_3(x) = x^2 + x + 1")
    print(f"has the SAME FORM as I(x) = 1 + a1*x + a2*x^2 when a1=a2=1.")
    print(f"")
    print(f"When does I(x) = Phi_3(x)?")
    print(f"  a1 = 1, a2 = 1: exactly 1 odd cycle and 1 disjoint pair")
    print(f"  But 1 cycle and 1 disjoint pair is IMPOSSIBLE (need 2 cycles for a pair)")
    print(f"  So I(x) = Phi_3(x) is NEVER achievable!")
    print(f"  Interesting: the polynomial defining omega is itself unreachable.")

    # What about multiples?
    print(f"\n  I(x) = k * Phi_3(x) requires a1 = k, a2 = k")
    print(f"  Need k odd cycles with k disjoint pairs")
    print(f"  k disjoint pairs from k cycles: need C(k,2) >= k, so k >= 3")

    # Check if a1=a2 ever occurs
    rng = np.random.default_rng(2026)
    equal_count = 0
    equal_examples = []
    for _ in range(3000):
        A = random_tournament(7, rng)
        a1, a2 = get_alpha_1_2(A)
        if a1 == a2:
            equal_count += 1
            if len(equal_examples) < 5:
                equal_examples.append((a1, a2))

    print(f"\n  a1 = a2 frequency at n=7: {equal_count}/3000 ({100*equal_count/3000:.1f}%)")
    if equal_examples:
        for a1, a2 in equal_examples:
            H = 1 + 2*a1 + 4*a2
            N = eisenstein_norm(a1, a2)
            print(f"    a1=a2={a1}: H={H}, N={N}")
            print(f"    I(x) = 1 + {a1}x + {a1}x^2 = {a1}*(x^2+x+1) + (1-{a1})")
            print(f"         = {a1}*Phi_3(x) + {1-a1}")

    print(f"\n  When a1 = a2 = a:")
    print(f"    I(x) = 1 + a*x + a*x^2 = a*(x^2+x+1) + (1-a) = a*Phi_3(x) + (1-a)")
    print(f"    I(omega) = a*Phi_3(omega) + (1-a) = a*0 + (1-a) = 1-a")
    print(f"    So I(omega) is a REAL INTEGER when a1=a2!")
    print(f"    Norm: N(1-a) = (1-a)^2")
    print(f"    These are 'Dehn-trivial' tournaments (real I(omega))")

# ============================================================
# PART 8: Norm filtration
# ============================================================
def part8():
    print("\n" + "=" * 70)
    print("PART 8: NORM FILTRATION — TOURNAMENTS BY EISENSTEIN NORM")
    print("=" * 70)

    rng = np.random.default_rng(2026_0314)
    n = 7
    N_samples = 3000

    norm_to_h = defaultdict(set)
    norm_to_a = defaultdict(set)
    norm_count = defaultdict(int)

    for _ in range(N_samples):
        A = random_tournament(n, rng)
        a1, a2 = get_alpha_1_2(A)
        H = 1 + 2*a1 + 4*a2
        N = eisenstein_norm(a1, a2)
        norm_to_h[N].add(H)
        norm_to_a[N].add((a1, a2))
        norm_count[N] += 1

    print(f"\nEisenstein norm distribution at n=7 ({N_samples} samples):")
    print(f"{'Norm':>6} {'Count':>6} {'%':>6} {'H values':>30} {'mod3':>5}")
    for N in sorted(norm_count.keys())[:30]:
        h_list = sorted(norm_to_h[N])
        h_str = str(h_list) if len(h_list) <= 6 else f"[{h_list[0]}..{h_list[-1]}] ({len(h_list)} vals)"
        pct = 100 * norm_count[N] / N_samples
        print(f"  {N:4d} {norm_count[N]:6d} {pct:5.1f}% {h_str:>35s}  {N%3}")

    # Is norm a good separator?
    print(f"\nNorm uniquely determines H? ", end="")
    unique = sum(1 for N in norm_to_h if len(norm_to_h[N]) == 1)
    total = len(norm_to_h)
    print(f"No: {unique}/{total} norms have unique H ({100*unique/total:.0f}%)")

    print(f"\nH uniquely determines norm? Check by collecting norm per H:")
    h_to_norms = defaultdict(set)
    for N, h_set in norm_to_h.items():
        for h in h_set:
            h_to_norms[h].add(N)
    multi = sum(1 for h in h_to_norms if len(h_to_norms[h]) > 1)
    print(f"  {multi} H values map to multiple norms")
    for h in sorted(h_to_norms.keys())[:15]:
        norms = sorted(h_to_norms[h])
        print(f"    H={h:3d}: norms = {norms}")

# ============================================================
# PART 9: Dehn-like obstruction for H achievability
# ============================================================
def part9():
    print("\n" + "=" * 70)
    print("PART 9: DEHN-LIKE OBSTRUCTION FOR H ACHIEVABILITY")
    print("=" * 70)

    print(f"\nIdea: Use the Eisenstein norm as an OBSTRUCTION for H values.")
    print(f"If H=h, then a1+2*a2 = (h-1)/2 = k.")
    print(f"The achievable (a1, a2) on the line a1 = k-2*a2 must satisfy:")
    print(f"  1. a1, a2 >= 0")
    print(f"  2. (a1, a2) is realizable by some tournament at n=7")
    print(f"  3. The Eisenstein norm N(a1, a2) must be achievable")

    # For forbidden H values, check what norms would be needed
    forbidden = [7, 21, 63, 107, 119, 149]
    for h in forbidden:
        k = (h - 1) // 2
        print(f"\n  H = {h} (k = {k}):")
        any_possible = False
        for a2 in range(0, k // 2 + 1):
            a1 = k - 2 * a2
            if a1 < 0:
                break
            N = eisenstein_norm(a1, a2)
            # Check feasibility: a2 <= C(a1, 2) and a1, a2 within observed ranges
            max_a2 = a1 * (a1 - 1) // 2 if a1 >= 2 else 0
            feasible = "possible" if a2 <= max_a2 else "IMPOSSIBLE (a2 > C(a1,2))"
            if a2 > max_a2:
                feasible = "IMPOSSIBLE"
            print(f"    (a1={a1:3d}, a2={a2:3d}): N={N:5d}, N mod 3 = {N%3}, {feasible}")
            if a2 <= max_a2:
                any_possible = True
        if not any_possible:
            print(f"    ALL decompositions impossible by graph bound!")

    print(f"\nKey insight: for small k (like k=3 giving H=7), ALL (a1,a2) are")
    print(f"individually ruled out by combinatorial constraints.")
    print(f"For larger k (like k=10 giving H=21), the Eisenstein norm")
    print(f"doesn't directly forbid, but tournament realizability does.")

# ============================================================
# PART 10: Scissors congruence of tournaments
# ============================================================
def part10():
    print("\n" + "=" * 70)
    print("PART 10: SCISSORS CONGRUENCE OF TOURNAMENTS")
    print("=" * 70)

    print(f"\nDefine: Two tournaments T1, T2 are 'scissors-congruent' if")
    print(f"  there exists a sequence of arc flips transforming T1 to T2")
    print(f"  such that at each step, some 'measure' is preserved.")
    print(f"")
    print(f"Volume analog: H(T) = total Hamiltonian path count")
    print(f"Dehn analog: I(omega) = independence polynomial at cube root of unity")
    print(f"")
    print(f"The classification:")
    print(f"  Same H AND same I(omega) => 'scissors-congruent'")
    print(f"  Same H, different I(omega) => 'Dehn-inequivalent'")
    print(f"  Different H => not even 'volume-equivalent'")

    print(f"\nIn precise terms: (a1, a2) determines BOTH H and I(omega)")
    print(f"So the classification reduces to the (a1, a2) lattice.")
    print(f"Two tournaments are in the same 'congruence class' iff")
    print(f"they have the same (a1, a2).")

    print(f"\nThe achievable (a1, a2) lattice at n=7:")
    rng = np.random.default_rng(42)
    pairs = defaultdict(int)
    for _ in range(5000):
        A = random_tournament(7, rng)
        a1, a2 = get_alpha_1_2(A)
        pairs[(a1, a2)] += 1

    # Show the lattice
    a1_max = max(a for a, _ in pairs)
    a2_max = max(b for _, b in pairs)
    print(f"  a1 range: [0, {a1_max}], a2 range: [0, {a2_max}]")
    print(f"  Total distinct (a1,a2) pairs: {len(pairs)}")

    # Show the achievable H from each pair
    h_set = set()
    for (a1, a2), count in sorted(pairs.items()):
        H = 1 + 2*a1 + 4*a2
        h_set.add(H)

    all_possible_h = set(range(1, max(h_set)+1, 2))
    gaps = sorted(all_possible_h - h_set)
    print(f"\n  Achievable H values: {len(h_set)}")
    print(f"  Gaps in odd numbers up to {max(h_set)}: {gaps[:20]}")

    print(f"\nThe Dehn invariant story for tournaments:")
    print(f"  Hilbert: 'Can equal-volume polyhedra always be cut into equal pieces?'")
    print(f"  Tournament: 'Can equal-H tournaments always be related by arc flips")
    print(f"              preserving H at every step?'")
    print(f"  Answer (like Dehn): NO, because I(omega) provides additional obstruction.")
    print(f"  Two tournaments with same H but different (a1,a2) decomposition")
    print(f"  are 'Dehn-inequivalent' — they represent the same H via")
    print(f"  fundamentally different cycle structures.")

    print(f"\n  Example: H=45 can be achieved by:")
    h45_pairs = [(a1, a2) for (a1, a2) in pairs if 1+2*a1+4*a2 == 45]
    for a1, a2 in sorted(h45_pairs):
        N = eisenstein_norm(a1, a2)
        print(f"    (a1={a1:2d}, a2={a2}): N={N}, I(omega) = ({1-a2}) + ({a1-a2})*omega")

    print(f"\n{'='*70}")
    print(f"END OF HILBERT'S 3RD EXPLORATION")
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

if __name__ == "__main__":
    main()
