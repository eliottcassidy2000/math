"""
z_omega_tower.py — The Z[omega] / Z[zeta_8] / Z[2^(1/3)] tower
and connections to Cayley-Dickson and tournament inflection.
kind-pasteur-2026-03-14-S65

PART 1:  Z[omega] norm spectrum — which norms are achievable?
PART 2:  I(omega) mod 3 — the reduced tournament invariant
PART 3:  I(zeta_4) = I(i) — what does evaluation at i give us?
PART 4:  I(zeta_8) — the 8th root of unity evaluation
PART 5:  Norm comparison: N_3(I(omega)) vs N_4(I(i)) vs N_8(I(zeta_8))
PART 6:  Quaternionic structure: I(i) = (1-a2) + a1*i (Gaussian integer)
PART 7:  The four evaluations: I(1), I(2), I(omega), I(i) — a quaternion?
PART 8:  Cayley-Dickson doubling and tournament invariant doubling
PART 9:  The discriminant 4*a2 - a1^2 and the phase transition
PART 10: H mod 8 structure (tournament octonionic residue)
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

def eisenstein_norm(a1, a2):
    return a1**2 - a1*a2 + a2**2 - a1 - a2 + 1

def gaussian_norm(a1, a2):
    """N(I(i)) where I(i) = (1-a2) + a1*i"""
    return (1-a2)**2 + a1**2

# ============================================================
# PART 1: Z[omega] norm spectrum
# ============================================================
def part1():
    print("=" * 70)
    print("PART 1: Z[omega] NORM SPECTRUM — ACHIEVABLE NORMS")
    print("=" * 70)

    print(f"\nN(a + b*omega) = a^2 - a*b + b^2")
    print(f"This is a positive definite quadratic form (discriminant -3)")
    print(f"")
    print(f"Which positive integers are representable as a^2 - ab + b^2?")
    print(f"Answer: n is representable iff all prime factors p = 2 mod 3 appear")
    print(f"with EVEN exponent. Equivalently: n has no prime factor = 2 mod 3")
    print(f"with odd exponent.")

    # List representable norms up to 50
    representable = set()
    for a in range(-50, 51):
        for b in range(-50, 51):
            n = a**2 - a*b + b**2
            if 0 < n <= 100:
                representable.add(n)

    print(f"\nRepresentable norms up to 100: {sorted(representable)}")

    # Non-representable (these are exactly those with odd power of prime = 2 mod 3)
    non_rep = set(range(1, 101)) - representable
    print(f"Non-representable: {sorted(non_rep)}")

    print(f"\nPrimes = 2 mod 3 up to 100: ", end="")
    primes_2mod3 = [p for p in range(2, 101) if all(p % d != 0 for d in range(2, p)) and p % 3 == 2]
    print(primes_2mod3)

    print(f"\nNote: the Eisenstein norm for tournaments is NOT the standard norm.")
    print(f"Tournament norm: N(I(omega)) = a1^2 - a1*a2 + a2^2 - a1 - a2 + 1")
    print(f"                             = N_std(a1 - (a2+1)/2, ...) - shifted")
    print(f"Let's check which TOURNAMENT norms occur:")

    rng = np.random.default_rng(2026_0314)
    tour_norms = set()
    for _ in range(5000):
        A = random_tournament(7, rng)
        a1, a2 = get_alpha_1_2(A)
        N = eisenstein_norm(a1, a2)
        tour_norms.add(N)

    tour_norms_list = sorted(tour_norms)
    print(f"  Achieved tournament norms: {tour_norms_list[:30]}...")
    print(f"  Total distinct norms: {len(tour_norms_list)}")

    # Are there gaps?
    all_norms = set(range(min(tour_norms_list), max(tour_norms_list)+1))
    gaps = sorted(all_norms - tour_norms)
    print(f"  Gap norms (not achieved): {gaps[:30]}...")

    # Check which norms mod 3 occur
    mod3_dist = defaultdict(int)
    for n in tour_norms:
        mod3_dist[n % 3] += 1
    print(f"\n  Tournament norms mod 3: {dict(sorted(mod3_dist.items()))}")
    print(f"  NO norm = 2 mod 3? {2 not in mod3_dist or mod3_dist[2] == 0}")

# ============================================================
# PART 2: I(omega) mod 3
# ============================================================
def part2():
    print("\n" + "=" * 70)
    print("PART 2: I(omega) MOD 3 — THE REDUCED INVARIANT")
    print("=" * 70)

    print(f"\nI(omega) = (1-a2) + (a1-a2)*omega in Z[omega]")
    print(f"Reducing mod 3:")
    print(f"  I(omega) mod 3 = ((1-a2) mod 3) + ((a1-a2) mod 3)*omega")
    print(f"  Lives in F_9 = F_3[omega] = (Z/3Z)[x]/(x^2+x+1)")

    rng = np.random.default_rng(2026)
    n = 7
    reduced = defaultdict(int)
    for _ in range(5000):
        A = random_tournament(n, rng)
        a1, a2 = get_alpha_1_2(A)
        r = (1 - a2) % 3
        s = (a1 - a2) % 3
        reduced[(r, s)] += 1

    print(f"\nI(omega) mod 3 distribution at n=7 (5000 samples):")
    for (r, s) in sorted(reduced.keys()):
        count = reduced[(r, s)]
        pct = 100 * count / 5000
        label = f"{r} + {s}*omega"
        norm_mod3 = (r**2 - r*s + s**2) % 3
        print(f"  {label:12s}: {count:4d} ({pct:5.1f}%), norm mod 3 = {norm_mod3}")

    print(f"\nF_9 has 9 elements, but only some appear as I(omega) mod 3.")
    print(f"The non-appearing elements correspond to 'forbidden residues'.")

    # The 9 elements of F_9 and their norms
    print(f"\nAll 9 elements of F_9 = F_3[omega]:")
    for r in range(3):
        for s in range(3):
            norm = (r**2 - r*s + s**2) % 3
            count = reduced.get((r, s), 0)
            print(f"  {r} + {s}*omega: norm = {norm}, occurs {count} times", end="")
            if count == 0:
                print(" <-- NEVER")
            else:
                print()

# ============================================================
# PART 3: I(i) — Gaussian integer evaluation
# ============================================================
def part3():
    print("\n" + "=" * 70)
    print("PART 3: I(i) — EVALUATION AT 4TH ROOT OF UNITY (GAUSSIAN)")
    print("=" * 70)

    print(f"\nI(i) = 1 + a1*i + a2*i^2 = 1 + a1*i - a2 = (1-a2) + a1*i")
    print(f"This is a Gaussian integer in Z[i]!")
    print(f"Gaussian norm: N(I(i)) = (1-a2)^2 + a1^2")

    rng = np.random.default_rng(2026_0314)
    n = 7
    gauss_data = defaultdict(int)
    gauss_norms = defaultdict(int)

    for _ in range(3000):
        A = random_tournament(n, rng)
        a1, a2 = get_alpha_1_2(A)
        real = 1 - a2
        imag = a1
        gn = real**2 + imag**2
        gauss_data[(real, imag)] += 1
        gauss_norms[gn] += 1

    print(f"\nGaussian norm distribution at n=7 (3000 samples):")
    for gn in sorted(gauss_norms.keys())[:20]:
        count = gauss_norms[gn]
        pct = 100 * count / 3000
        print(f"  N_gauss = {gn:6d}: {count:4d} ({pct:5.1f}%)")

    print(f"\nComparison: Eisenstein vs Gaussian norm")
    rng2 = np.random.default_rng(42)
    print(f"  {'a1':>4} {'a2':>4} {'H':>5} {'N_eis':>6} {'N_gau':>6} {'N_eis/N_gau':>10}")
    for _ in range(15):
        A = random_tournament(n, rng2)
        a1, a2 = get_alpha_1_2(A)
        H = 1 + 2*a1 + 4*a2
        ne = eisenstein_norm(a1, a2)
        ng = gaussian_norm(a1, a2)
        ratio = ne/ng if ng > 0 else float('inf')
        print(f"  {a1:4d} {a2:4d} {H:5d} {ne:6d} {ng:6d} {ratio:10.4f}")

    print(f"\nRelation between the two norms:")
    print(f"  N_eis = a1^2 - a1*a2 + a2^2 - a1 - a2 + 1")
    print(f"  N_gau = a1^2 + (1-a2)^2 = a1^2 + a2^2 - 2*a2 + 1")
    print(f"  Difference: N_eis - N_gau = -a1*a2 + 2*a2^2 - a1 + a2")
    print(f"            = a2*(2*a2 - a1) - (a1 - a2)")
    print(f"            = (a2 - 1)*(2*a2 - a1) + (a2 - a1)")

    # When are they equal?
    print(f"\n  N_eis = N_gau iff -a1*a2 + 2*a2^2 - a1 + a2 = 0")
    print(f"                iff a2*(2*a2 - a1 + 1) = a1")
    print(f"  For a2=0: a1=0 (trivial)")
    print(f"  For a2=1: 1*(2 - a1 + 1) = a1 => 3 - a1 = a1 => a1 = 3/2 (no integer)")
    print(f"  So they are equal ONLY at (0,0) among non-negative integers!")

# ============================================================
# PART 4: I(zeta_8) — 8th root of unity
# ============================================================
def part4():
    print("\n" + "=" * 70)
    print("PART 4: I(zeta_8) — 8TH ROOT OF UNITY EVALUATION")
    print("=" * 70)

    zeta8 = complex(math.cos(2*math.pi/8), math.sin(2*math.pi/8))
    print(f"\nzeta_8 = e^(2*pi*i/8) = cos(pi/4) + i*sin(pi/4) = (1+i)/sqrt(2)")
    print(f"       = {zeta8}")
    print(f"  |zeta_8| = 1")
    print(f"  zeta_8^2 = i")
    print(f"  zeta_8^4 = -1")
    print(f"  zeta_8^8 = 1")

    print(f"\nI(zeta_8) = 1 + a1*zeta_8 + a2*zeta_8^2 = 1 + a1*zeta_8 + a2*i")
    print(f"  Real part: 1 + a1*cos(pi/4) = 1 + a1/sqrt(2)")
    print(f"  Imag part: a1*sin(pi/4) + a2 = a1/sqrt(2) + a2")

    print(f"\n  Z[zeta_8] has rank phi(8) = 4 as Z-module")
    print(f"  Basis: 1, zeta_8, zeta_8^2=i, zeta_8^3 = i*zeta_8")
    print(f"  I(zeta_8) = 1 + a1*zeta_8 + a2*zeta_8^2")
    print(f"  Only uses 3 of the 4 basis elements (no zeta_8^3 component)")

    print(f"\nNorm in Z[zeta_8]: N(alpha) = |alpha|^2 * |sigma(alpha)|^2")
    print(f"  where sigma is the Galois conjugation zeta_8 -> zeta_8^3")
    print(f"  sigma(I(zeta_8)) = 1 + a1*zeta_8^3 + a2*zeta_8^6")
    print(f"                   = 1 + a1*zeta_8^3 + a2*(-i)")

    # Compute |I(zeta_8)|^2 for some tournaments
    rng = np.random.default_rng(2026)
    print(f"\n  Sample |I(zeta_8)|^2 values at n=7:")
    for _ in range(10):
        A = random_tournament(7, rng)
        a1, a2 = get_alpha_1_2(A)
        val = 1 + a1*zeta8 + a2*zeta8**2
        mod_sq = abs(val)**2
        # Also sigma
        val_sigma = 1 + a1*zeta8**3 + a2*zeta8**6
        full_norm = abs(val)**2 * abs(val_sigma)**2
        print(f"    a1={a1:3d}, a2={a2:2d}: |I(zeta_8)|^2 = {mod_sq:.2f}, "
              f"full norm = {full_norm:.0f}")

    print(f"\n  |I(zeta_8)|^2 = (1 + a1/sqrt(2))^2 + (a1/sqrt(2) + a2)^2")
    print(f"  = 1 + a1*sqrt(2) + a1^2/2 + a1^2/2 + a1*a2*sqrt(2) + a2^2")
    print(f"  = 1 + a1^2 + a2^2 + (a1 + a1*a2)*sqrt(2)")
    print(f"  NOTE: this has an IRRATIONAL part (involving sqrt(2))!")
    print(f"  So |I(zeta_8)|^2 is NOT in Z unless a1*(1+a2) = 0")
    print(f"  i.e., a1=0 or a2=-1 (impossible since a2>=0)")
    print(f"  So only for trivial a1=0, |I(zeta_8)|^2 is rational!")

# ============================================================
# PART 5: Norm comparison
# ============================================================
def part5():
    print("\n" + "=" * 70)
    print("PART 5: NORM COMPARISON ACROSS ROOTS OF UNITY")
    print("=" * 70)

    omega = complex(-0.5, math.sqrt(3)/2)
    i_unit = complex(0, 1)

    rng = np.random.default_rng(2026_0314)
    n = 7

    print(f"\nComparing |I(x)|^2 at x = omega, i, -1:")
    print(f"  {'a1':>4} {'a2':>4} {'H':>5} {'|I(w)|^2':>10} {'|I(i)|^2':>10} {'I(-1)':>6} {'I(-1)^2':>8}")
    for _ in range(15):
        A = random_tournament(n, rng)
        a1, a2 = get_alpha_1_2(A)
        H = 1 + 2*a1 + 4*a2
        iw = abs(1 + a1*omega + a2*omega**2)**2
        ii = abs(1 + a1*i_unit + a2*i_unit**2)**2
        im1 = 1 - a1 + a2
        print(f"  {a1:4d} {a2:4d} {H:5d} {iw:10.2f} {ii:10.2f} {im1:6d} {im1**2:8d}")

    print(f"\nFormulas:")
    print(f"  |I(omega)|^2 = N_eis = a1^2 - a1*a2 + a2^2 - a1 - a2 + 1")
    print(f"               (integer, always >= 0, never 2 mod 3)")
    print(f"  |I(i)|^2     = N_gau = a1^2 + (1-a2)^2")
    print(f"               (integer, sum of two squares)")
    print(f"  I(-1)^2      = (1-a1+a2)^2")
    print(f"               (perfect square)")

    print(f"\nProduct of norms: |I(omega)|^2 * |I(i)|^2 * I(-1)^2")
    for _ in range(10):
        A = random_tournament(n, rng)
        a1, a2 = get_alpha_1_2(A)
        ne = eisenstein_norm(a1, a2)
        ng = gaussian_norm(a1, a2)
        nm1 = (1 - a1 + a2)**2
        product = ne * ng * nm1
        # Also: product mod 12
        print(f"  a1={a1:3d}, a2={a2:2d}: N_e={ne:5d}, N_g={ng:5d}, "
              f"N_-1={nm1:5d}, product={product:10d}")

# ============================================================
# PART 6: Quaternionic structure of I(i)
# ============================================================
def part6():
    print("\n" + "=" * 70)
    print("PART 6: QUATERNIONIC STRUCTURE OF I(i)")
    print("=" * 70)

    print(f"\nI(i) = (1-a2) + a1*i = Gaussian integer")
    print(f"Gaussian integers Z[i] are the 'complex part' of quaternions H")
    print(f"")
    print(f"Can we define a 'quaternionic tournament invariant'?")
    print(f"  Q(T) = I(1) + I(i)*j = (1+a1+a2) + ((1-a2) + a1*i)*j")
    print(f"  = (1+a1+a2) + (1-a2)*j + a1*ij")
    print(f"  = (1+a1+a2) + (1-a2)*j + a1*k")
    print(f"  (using quaternion convention ij = k)")

    rng = np.random.default_rng(2026)
    print(f"\nQuaternionic tournament invariant Q(T) = w + xi + yj + zk:")
    for _ in range(10):
        A = random_tournament(7, rng)
        a1, a2 = get_alpha_1_2(A)
        w = 1 + a1 + a2  # I(1)
        x = 0  # no i component from I(1)
        y = 1 - a2       # j component from I(i)
        z = a1            # k component
        H = 1 + 2*a1 + 4*a2
        qnorm = w**2 + x**2 + y**2 + z**2
        print(f"  a1={a1:3d}, a2={a2:2d}: Q = {w} + {x}i + {y}j + {z}k, "
              f"|Q|^2 = {qnorm}, H = {H}")

    print(f"\nQuaternionic norm |Q|^2 = (1+a1+a2)^2 + (1-a2)^2 + a1^2")
    print(f"  = I(1)^2 + N_gaussian(I(i))")
    print(f"  = (total indep sets + 1)^2 + Gaussian norm")

    print(f"\nRelation to H:")
    print(f"  H = I(2) = 1 + 2*a1 + 4*a2")
    print(f"  I(1) = 1 + a1 + a2")
    print(f"  N_gau = a1^2 + (1-a2)^2")
    print(f"  |Q|^2 = I(1)^2 + N_gau = (1+a1+a2)^2 + a1^2 + (1-a2)^2")

    print(f"\n  Can we recover H from |Q|^2? Probably not uniquely,")
    print(f"  but |Q|^2 encodes STRICTLY MORE information than H.")

    # What about octonionic?
    print(f"\nOctonionic extension (dim 8 = 2^3):")
    print(f"  Would need evaluations at 4 independent points to fill 8 components")
    print(f"  e.g., I(1), I(i), I(omega), I(omega*i)")
    print(f"  But I(x) is degree 2, so 3 points determine it completely")
    print(f"  An octonionic invariant would be REDUNDANT for quadratic I(x)")
    print(f"  This is why n=7 (quadratic I.P.) doesn't 'see' the octonionic structure")
    print(f"  At n>=8, I(x) can be cubic => need 4 points => octonionic NOT redundant!")
    print(f"  THIS is the Cayley-Dickson connection: at n=8, the invariant")
    print(f"  needs more than quaternionic (dim 4) structure!")

# ============================================================
# PART 7: Four evaluations as quaternion
# ============================================================
def part7():
    print("\n" + "=" * 70)
    print("PART 7: THE FOUR EVALUATIONS I(1), I(2), I(omega), I(i)")
    print("=" * 70)

    omega = complex(-0.5, math.sqrt(3)/2)

    print(f"\nI(x) = 1 + a1*x + a2*x^2 (degree 2)")
    print(f"Three points determine I(x) uniquely.")
    print(f"With FOUR points, we have one redundancy = one CHECK.")
    print(f"")
    print(f"The four evaluations:")
    print(f"  I(1) = 1 + a1 + a2         (total independent sets + empty set)")
    print(f"  I(2) = 1 + 2*a1 + 4*a2     (= H, Hamiltonian paths)")
    print(f"  I(omega) = (1-a2)+(a1-a2)w  (Eisenstein integer)")
    print(f"  I(i) = (1-a2) + a1*i        (Gaussian integer)")

    print(f"\nThe CHECK relation:")
    print(f"  I(2) = 2*I(1) + 2*a2 - 1")
    print(f"  H = 2*(1+a1+a2) + 2*a2 - 1 = 1 + 2*a1 + 4*a2  (correct)")
    print(f"")
    print(f"  More interesting: I(2) from I(omega) and I(i)?")
    print(f"  I(omega): real = 1-a2, omega-coeff = a1-a2")
    print(f"    => a2 = 1-real, a1 = omega-coeff + (1-real)")
    print(f"  I(i): real = 1-a2, imag = a1")
    print(f"    => a2 = 1-Re(I(i)), a1 = Im(I(i))")
    print(f"  Then H = 1 + 2*Im(I(i)) + 4*(1 - Re(I(i)))")
    print(f"         = 5 - 4*Re(I(i)) + 2*Im(I(i))")

    rng = np.random.default_rng(42)
    print(f"\nVerification at n=7:")
    for _ in range(8):
        A = random_tournament(7, rng)
        a1, a2 = get_alpha_1_2(A)
        H = 1 + 2*a1 + 4*a2
        i1 = 1 + a1 + a2
        ii = complex(1-a2, a1)
        iw = complex(1-a2, 0) + complex(a1-a2, 0) * omega
        h_check = 5 - 4*(1-a2) + 2*a1
        print(f"  a1={a1:3d}, a2={a2:2d}: H={H}, I(1)={i1}, "
              f"Re(I(i))={1-a2}, Im(I(i))={a1}, "
              f"5-4Re+2Im={int(h_check)}")

    print(f"\nThe formula H = 5 - 4*Re(I(i)) + 2*Im(I(i)) shows:")
    print(f"  H increases with a1 (more cycles -> more paths)")
    print(f"  H increases with a2 (more disjoint pairs -> more paths)")
    print(f"  The coefficients 4 and 2 are POWERS OF 2")
    print(f"  And the constant 5 = 2+3!")

# ============================================================
# PART 8: Cayley-Dickson doubling and I.P. degree
# ============================================================
def part8():
    print("\n" + "=" * 70)
    print("PART 8: CAYLEY-DICKSON DOUBLING AND I.P. DEGREE TRANSITION")
    print("=" * 70)

    print(f"\nI(x) = 1 + a1*x + a2*x^2 + a3*x^3 + ...")
    print(f"The maximum independence number alpha(Omega(T)) determines the degree.")
    print(f"")
    print(f"  n <= 4: alpha(Omega) = 0 or 1, degree <= 1 (real = dim 1)")
    print(f"  n = 5,6,7: alpha(Omega) = 0,1, or 2, degree <= 2 (complex = dim 2)")
    print(f"  n = 8: alpha(Omega) can reach 3, degree <= 3")
    print(f"  n = 9+: alpha(Omega) can reach 4+, degree <= 4+")

    print(f"\nCayley-Dickson parallel:")
    print(f"  Degree 1 (linear I.P.): R (dim 1) — fully ordered")
    print(f"  Degree 2 (quadratic I.P.): C (dim 2) — lose ordering")
    print(f"  Degree 3 (cubic I.P.): between C and H")
    print(f"  Degree 4 (quartic I.P.): H (dim 4) — lose commutativity")

    print(f"\nBut wait: Cayley-Dickson doubles, while I.P. degree increments by 1.")
    print(f"The connection is through EVALUATION POINTS needed:")
    print(f"  Degree d polynomial needs d+1 points to determine it")
    print(f"  d=1: 2 points = 2^1 (complex)")
    print(f"  d=2: 3 points (not a power of 2!)")
    print(f"  d=3: 4 points = 2^2 (quaternionic)")
    print(f"  d=7: 8 points = 2^3 (octonionic)")

    print(f"\n  The Cayley-Dickson step at 2^k evaluations corresponds to")
    print(f"  I.P. degree 2^k - 1:")
    print(f"    2^1 - 1 = 1: linear I.P. (n <= 4)")
    print(f"    2^2 - 1 = 3: cubic I.P. (n = 8)")
    print(f"    2^3 - 1 = 7: degree-7 I.P. (n = ?)")

    print(f"\n  For degree 3 (n=8): the 4 evaluations I(1), I(2), I(omega), I(i)")
    print(f"  form a COMPLETE basis for recovering I(x).")
    print(f"  This is the 'quaternionic' level of information.")
    print(f"  The fact that n=8 = 2^3 matches the Cayley-Dickson dim for")
    print(f"  quaternion -> octonion transition is suggestive but not proven.")

    # Check what alpha(Omega) values occur at n=8 (just estimate)
    print(f"\nEstimating max alpha(Omega) at n=7 vs n=8:")
    print(f"  n=7: alpha(Omega) <= 2 always (proved, Omega has <= 14+21+1=36 vertices)")
    print(f"  n=8: alpha(Omega) can reach 3 (3 mutually disjoint odd cycles)")
    print(f"  Example: three 3-cycles on disjoint vertex triples")
    print(f"    Vertices {{0,1,2}}, {{3,4,5}}, {{6,7}}: only 2 triples of 3 from 8 vertices")
    print(f"    Wait: C(8,3)=56, can we find 3 DISJOINT 3-cycles?")
    print(f"    Need 9 vertices for 3 disjoint triples, but n=8 has only 8!")
    print(f"    So max disjoint 3-cycles at n=8 is FLOOR(8/3) = 2")
    print(f"    But: a 3-cycle and a 5-cycle on disjoint vertices: 3+5=8 = n!")
    print(f"    So alpha(Omega) >= 2 at n=8, and a1=3 = 1 disjoint 3-5 pair + 1 cycle")
    print(f"    Can we get 3 mutually disjoint odd cycles? No, because 3+3+3=9 > 8")
    print(f"    But 3+5 = 8, and we need a third cycle disjoint from both.")
    print(f"    That's impossible: 8 vertices are all used up!")
    print(f"    So alpha(Omega) = 2 at n=8 as well! (for 3-cycles and 5-cycles)")
    print(f"")
    print(f"    What about at n=9? 3+3+3=9: three disjoint 3-cycles possible!")
    print(f"    So alpha(Omega) can reach 3 at n=9 = 3^2")
    print(f"    The CUBIC transition happens at n=9, not n=8!")
    print(f"    This matches: 9 = 3^2, and cubic I.P. = degree 3")

# ============================================================
# PART 9: Discriminant and phase transition
# ============================================================
def part9():
    print("\n" + "=" * 70)
    print("PART 9: DISCRIMINANT 4*a2 - a1^2 AND PHASE TRANSITION")
    print("=" * 70)

    print(f"\nI(x) = 1 + a1*x + a2*x^2")
    print(f"Roots: x = (-a1 +- sqrt(a1^2 - 4*a2)) / (2*a2)")
    print(f"Discriminant: Delta = a1^2 - 4*a2")
    print(f"  Delta > 0: two real roots (hyperbolic)")
    print(f"  Delta = 0: double root (parabolic)")
    print(f"  Delta < 0: two complex roots (elliptic)")

    rng = np.random.default_rng(2026_0314)
    n = 7
    disc_dist = defaultdict(int)
    type_count = {"hyperbolic": 0, "parabolic": 0, "elliptic": 0}

    for _ in range(5000):
        A = random_tournament(n, rng)
        a1, a2 = get_alpha_1_2(A)
        if a2 == 0:
            disc_dist["a2=0 (linear)"] = disc_dist.get("a2=0 (linear)", 0) + 1
            continue
        delta = a1**2 - 4*a2
        if delta > 0:
            type_count["hyperbolic"] += 1
        elif delta == 0:
            type_count["parabolic"] += 1
        else:
            type_count["elliptic"] += 1

    print(f"\nRoot type distribution at n=7 (5000 samples):")
    for t, c in sorted(type_count.items()):
        print(f"  {t:12s}: {c:4d} ({100*c/5000:.1f}%)")
    print(f"  linear (a2=0): {disc_dist.get('a2=0 (linear)', 0)}")

    print(f"\nThe phase transition at Delta = 0:")
    print(f"  a1^2 = 4*a2 => a2 = a1^2/4")
    print(f"  For integer a2: a1 must be even, say a1 = 2m, then a2 = m^2")
    print(f"  Parabolic tournaments: (a1, a2) = (2m, m^2)")

    for m in range(0, 8):
        a1 = 2*m
        a2 = m**2
        H = 1 + 2*a1 + 4*a2
        ne = eisenstein_norm(a1, a2)
        dbl_root = f"{-a1/(2*a2):.4f}" if a2 > 0 else "inf"
        print(f"  m={m}: (a1={a1:2d}, a2={a2:2d}), H={H:4d}, N_eis={ne:4d}, "
              f"double root at x = {dbl_root}")

    # Are any parabolic tournaments achievable?
    print(f"\nAchievable parabolic points at n=7:")
    rng2 = np.random.default_rng(2026)
    parabolic_found = defaultdict(int)
    for _ in range(10000):
        A = random_tournament(7, rng2)
        a1, a2 = get_alpha_1_2(A)
        if a2 > 0 and a1**2 == 4*a2:
            parabolic_found[(a1, a2)] += 1

    if parabolic_found:
        for (a1, a2), count in sorted(parabolic_found.items()):
            print(f"  (a1={a1}, a2={a2}): found {count} times")
    else:
        print(f"  NONE found in 10000 samples!")
        print(f"  Parabolic tournaments may be very rare or non-existent")

    print(f"\n  Note: for (a1, a2) = (2, 1), Delta = 4-4 = 0")
    print(f"  This requires exactly 2 odd cycles with 1 disjoint pair")
    print(f"  Two disjoint odd cycles: possible at n>=6 (e.g., 3+3=6)")

    # Check at n=6
    rng3 = np.random.default_rng(42)
    para_n6 = 0
    for _ in range(5000):
        A = random_tournament(6, rng3)
        a1, a2 = get_alpha_1_2(A)
        if a1 == 2 and a2 == 1:
            para_n6 += 1
    print(f"  At n=6: (a1=2, a2=1) found {para_n6}/5000 times ({100*para_n6/5000:.1f}%)")

# ============================================================
# PART 10: H mod 8 structure
# ============================================================
def part10():
    print("\n" + "=" * 70)
    print("PART 10: H MOD 8 — THE OCTONIONIC RESIDUE")
    print("=" * 70)

    print(f"\nH = 1 + 2*a1 + 4*a2")
    print(f"H mod 2 = 1 always")
    print(f"H mod 4 = 1 + 2*(a1 mod 2)")
    print(f"H mod 8 = 1 + 2*a1 + 4*(a2 mod 2)  (since 8|4*2)")
    print(f"        = 1 + 2*(a1 mod 4) + 4*(a2 mod 2)")

    rng = np.random.default_rng(2026_0314)
    n = 7
    h_mod8 = defaultdict(int)
    h_mod4 = defaultdict(int)

    for _ in range(5000):
        A = random_tournament(n, rng)
        a1, a2 = get_alpha_1_2(A)
        H = 1 + 2*a1 + 4*a2
        h_mod8[H % 8] += 1
        h_mod4[H % 4] += 1

    print(f"\nH mod 4 distribution at n=7:")
    for r in range(4):
        count = h_mod4.get(r, 0)
        print(f"  H = {r} mod 4: {count:4d} ({100*count/5000:.1f}%)")

    print(f"\nH mod 8 distribution at n=7:")
    for r in range(8):
        count = h_mod8.get(r, 0)
        print(f"  H = {r} mod 8: {count:4d} ({100*count/5000:.1f}%)")

    print(f"\nExplanation:")
    print(f"  H mod 8 = (1 + 2*a1 + 4*a2) mod 8")
    print(f"  H is always odd, so only residues 1, 3, 5, 7 appear")
    print(f"  H mod 4: only 1 and 3 (since H = 1 + 2*a1 mod 4)")
    print(f"    H = 1 mod 4 iff a1 even")
    print(f"    H = 3 mod 4 iff a1 odd")

    print(f"\n  H mod 8:")
    print(f"    1 mod 8: a1 even, a2 even")
    print(f"    3 mod 8: a1 odd, a2 even")
    print(f"    5 mod 8: a1 even, a2 odd")
    print(f"    7 mod 8: a1 odd, a2 odd")

    # This gives the parity structure of (a1, a2) from H mod 8
    print(f"\n  H mod 8 determines (a1 mod 2, a2 mod 2)!")
    print(f"  This is a 2-dimensional mod-2 invariant")
    print(f"  = point in F_2 x F_2 = (Z/2Z)^2")
    print(f"  The group (Z/2Z)^2 has order 4 = 2^2 = dim(quaternions)")

    print(f"\n  Extended: H mod 16 determines (a1 mod 4, a2 mod 2)")
    print(f"  and H mod 32 determines (a1 mod 8, a2 mod 2)")
    print(f"  In general: H mod 2^k determines (a1 mod 2^(k-1), a2 mod 2^(k-2))")

    # The 2-adic valuation tower
    print(f"\n  2-adic structure of H-1:")
    print(f"  H - 1 = 2*(a1 + 2*a2)")
    print(f"  v_2(H-1) >= 1 always")
    print(f"  v_2(H-1) >= 2 iff a1 + 2*a2 even iff a1 even")
    print(f"  v_2(H-1) >= 3 iff a1 + 2*a2 = 0 mod 4 iff a1 = 0 mod 4")

    # Check v_2 distribution
    v2_dist = defaultdict(int)
    rng2 = np.random.default_rng(42)
    for _ in range(5000):
        A = random_tournament(n, rng2)
        a1, a2 = get_alpha_1_2(A)
        H = 1 + 2*a1 + 4*a2
        h_minus_1 = H - 1
        v2 = 0
        while h_minus_1 > 0 and h_minus_1 % 2 == 0:
            v2 += 1
            h_minus_1 //= 2
        v2_dist[v2] += 1

    print(f"\n  v_2(H-1) distribution at n=7:")
    for v in sorted(v2_dist.keys()):
        count = v2_dist[v]
        print(f"    v_2(H-1) = {v}: {count:4d} ({100*count/5000:.1f}%)")

    print(f"\n{'='*70}")
    print(f"END OF Z[omega] TOWER EXPLORATION")
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
