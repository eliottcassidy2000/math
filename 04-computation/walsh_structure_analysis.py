#!/usr/bin/env python3
"""
walsh_structure_analysis.py -- Deep analysis of Walsh-Fourier structure of H

At p=19, the Walsh coefficients J[i,j] take only 4 distinct absolute values.
This is because J[i,j] = f(|i-j| mod m) for some function f (circulant symmetry
acts on the orientation cube variables).

Key question: can we PROVE that the all-zeros vector (= Interval) maximizes
H(sigma) = H_0 + sigma^T J sigma + higher terms?

Strategy:
1. Show J[i,j] = g(i-j mod m) where g is a function on Z_m (J is circulant)
2. Show sigma_Int = (0,...,0) maps to s = (+1,...,+1) in {+1,-1} convention
3. sigma^T J sigma = sum_{i,j} J[i,j] s_i s_j
   For s = (+1,...,+1): Q_Int = sum_{i,j} J[i,j] = m * (sum of g values)
4. This is maximized when all J[i,j] have the same sign (cooperative Ising)
5. Higher-order terms: similar analysis for degree-4,6 Walsh coefficients

Also: analyze the orbit structure at p=13 and p=19 more carefully.
The 18 maximizers at p=19 = |Z_19^*| = 18. Same for p=13 (12 = |Z_13^*|).

Author: kind-pasteur-2026-03-12-S58
"""

import math


def analyze_walsh_circulant():
    """Check if J[i,j] = g(i-j mod m) (circulant structure)."""
    print("=" * 70)
    print("WALSH-FOURIER CIRCULANT STRUCTURE")
    print("=" * 70)

    # From orientation_cube_p19.out: Walsh coefficients J[i,j]
    # Variables 0..8 correspond to pairs (1,18), (2,17), ..., (9,10)
    m = 9
    # Read from output: degree-2 coefficients
    J_data = {
        (0,1): 180376253.45,
        (0,2): 207724893.82,
        (1,2): -6284752.76,
        (0,3): 122158124.55,
        (1,3): 180376253.45,
        (2,3): -122158124.55,
        (0,4): 122158124.55,
        (1,4): 6284752.76,
        (2,4): 6284752.76,
        (3,4): -207724893.82,
        (0,5): -207724893.82,
        (1,5): 207724893.82,
        (2,5): 180376253.45,
        (3,5): -6284752.76,
        (4,5): 122158124.55,
        (0,6): -6284752.76,
        (1,6): 207724893.82,
        (2,6): -122158124.55,
        (3,6): -207724893.82,
        (4,6): -180376253.45,
        (5,6): -180376253.45,
        (0,7): 6284752.76,
        (1,7): 122158124.55,
        (2,7): -180376253.45,
        (3,7): 180376253.45,
        (4,7): 207724893.82,
        (5,7): -122158124.55,
        (6,7): 6284752.76,
        (0,8): -180376253.45,
        (1,8): -122158124.55,
        (2,8): 207724893.82,
        (3,8): -6284752.76,
        (4,8): -180376253.45,
        (5,8): -6284752.76,
        (6,8): 122158124.55,
        (7,8): 207724893.82,
    }

    # Build full symmetric J matrix
    J = [[0.0]*m for _ in range(m)]
    for (i,j), v in J_data.items():
        J[i][j] = v
        J[j][i] = v

    # Check circulant: J[i,j] should depend only on (j-i) mod m
    print("\n  Checking if J[i,j] = g(j-i mod m):")
    g_values = {}
    for d in range(m):
        vals = []
        for i in range(m):
            j = (i + d) % m
            if i != j:
                vals.append(J[i][j])
        if vals:
            g_values[d] = vals
            all_same = all(abs(v - vals[0]) < 1.0 for v in vals)
            print(f"    d={d}: values = {[f'{v:.2f}' for v in vals[:5]]}{'...' if len(vals) > 5 else ''} "
                  f"{'CONSTANT' if all_same else 'NOT CONSTANT'}")

    # What IS the pattern?
    print("\n  J[i,j] values by (j-i) mod 9:")
    for d in range(1, m):
        vals = [J[i][(i+d) % m] for i in range(m)]
        distinct = sorted(set(round(v, 0) for v in vals))
        print(f"    d={d}: {len(distinct)} distinct values: {distinct}")

    # Instead, check if J[i,j] depends on a * (j-i) mod p for some multiplier
    # The orientation cube variable k corresponds to the pair (k+1, p-k-1)
    # In the Fourier world, the eigenvalues of the circulant adjacency matrix
    # at character chi_a give the trace contribution.
    p = 19

    # Actually, the key observation: variable k controls pair (k+1, p-k-1).
    # The Paley character chi(k) = Legendre(k+1/p) maps each variable to +/-1.
    # Walsh coefficient J[i,j] measures how flipping variables i and j together
    # affects H. By circulant symmetry of the TOURNAMENT (not the cube),
    # J[i,j] depends on the "difference" of the pairs, which is NOT just i-j.

    # Let's check: J[i,j] vs (i+1)*(j+1) mod p relationship
    print("\n  J[i,j] vs multiplicative structure:")
    for i in range(m):
        for j in range(i+1, m):
            a, b = i+1, j+1
            diff_add = (b - a) % p
            prod = (a * b) % p
            print(f"    J[{i},{j}] = {J[i][j]:>15.2f}  pair ({a},{b}), diff={diff_add:2d}, prod={prod:2d}")

    # Group by absolute value
    print("\n  Distinct |J| values:")
    abs_vals = sorted(set(round(abs(J[i][j]), 0) for i in range(m) for j in range(i+1, m)))
    for av in abs_vals:
        pairs = [(i,j) for i in range(m) for j in range(i+1, m) if abs(abs(J[i][j]) - av) < 1.0]
        print(f"    |J| = {av:.0f}: {len(pairs)} pairs")
        for i, j in pairs:
            a, b = i+1, j+1
            prod = (a * b) % p
            print(f"      ({i},{j}) = pair ({a},{b}), prod mod {p} = {prod}, Legendre = {pow(prod, (p-1)//2, p)}")

    # Q(Interval) = s_Int^T J s_Int where s_Int = (+1,...,+1) (sigma=0 => s=+1)
    Q_int = sum(J[i][j] for i in range(m) for j in range(m))
    print(f"\n  Q(Interval) = sum J[i,j] = {Q_int:.2f}")

    # Q(Paley): sigma_QR = (0, 1, 1, 0, 0, 0, 0, 1, 0), s_k = 1 - 2*sigma_k
    sigma_qr = [0, 1, 1, 0, 0, 0, 0, 1, 0]
    s_qr = [1 - 2*s for s in sigma_qr]
    Q_qr = sum(J[i][j] * s_qr[i] * s_qr[j] for i in range(m) for j in range(m))
    print(f"  Q(Paley) = {Q_qr:.2f}")
    print(f"  s_Paley = {s_qr}")

    # Optimal sigma for quadratic form
    # For sigma^T J sigma, the maximum over {+1,-1}^m is achieved by the
    # sign pattern of the dominant eigenvector
    print(f"\n  Key insight:")
    print(f"  Q(Int) = {Q_int:.2f} > Q(Paley) = {Q_qr:.2f}")
    print(f"  But SDP bound uses max_eig, not the sum.")
    print(f"  Sum of all J[i,j] = Q(Int) for s=(+1,...,+1)")


def orbit_structure():
    """Analyze why maximizers form orbits under Z_p^*."""
    print(f"\n{'='*70}")
    print("ORBIT STRUCTURE OF MAXIMIZERS")
    print("=" * 70)

    for p in [13, 19]:
        m = (p - 1) // 2
        S_int = set(range(1, m+1))

        print(f"\n  p={p}, m={m}, |Z_p^*| = {p-1}")

        # Orbit of S_int under multiplication by Z_p^*
        orbit = []
        for a in range(1, p):
            scaled = frozenset((a * s) % p for s in S_int)
            if scaled not in orbit:
                orbit.append(scaled)

        print(f"  |Orbit(Int)| = {len(orbit)}")
        print(f"  Expected maximizers = {p-1}")

        # Stabilizer: which a fix S_int?
        stab = [a for a in range(1, p) if frozenset((a*s) % p for s in S_int) == frozenset(S_int)]
        print(f"  Stabilizer = {stab}")
        print(f"  |Stab| = {len(stab)}")
        print(f"  |Orbit| = |Z_p^*|/|Stab| = {(p-1)//len(stab)}")

        # For p=19: orbit should have 18 elements (18 maximizers in output)
        # For p=13: orbit should have 12 elements (12 maximizers in output)

        # Key: S_int = {1,...,m}. Scaling by a gives {a, 2a, ..., ma} mod p.
        # S_int is fixed by a iff {a,...,ma} mod p = {1,...,m} (as sets).
        # This means a*{1,...,m} = {1,...,m} mod p.
        # i.e., a is in the stabilizer of {1,...,m} under Z_p^* action.

        # Equivalently: a*m mod p <= m for all elements (stays in lower half)
        # OR: the map x -> ax permutes {1,...,m}
        # This is very restrictive. For p prime, the only way is a=1
        # (except when -1 also stabilizes, which happens when p ≡ 1 mod 4
        #  and -S = {p-1,...,p-m} = {m+1,...,p-1} ≠ S)

        # For p=19: -1 mod 19 = 18. Does {18, 17, ..., 10} = {1,...,9}? No.
        # So stabilizer is {1}.
        # For p=13: -1 mod 13 = 12. Does {12, 11, 10, 9, 8, 7} = {1,...,6}? No.
        # So stabilizer is {1}.

        # Check which elements of orbit are ALSO valid tournament connection sets
        # (i.e., S ∩ (-S) = ∅)
        valid = 0
        for S_scaled in orbit:
            neg_S = frozenset(p - s for s in S_scaled)
            if not (S_scaled & neg_S):
                valid += 1
        print(f"  # valid tournament connection sets in orbit: {valid}/{len(orbit)}")

        if p % 4 == 3:
            # All should be valid since for p ≡ 3 mod 4, -1 is QNR
            # and S = a*{1,...,m} with a ∈ Z_p^* maps to -S = (-a)*{1,...,m}
            # which is a different orbit element (since -a ≠ a mod p when a ≠ 0)
            print(f"  (p ≡ 3 mod 4: all orbit elements are valid connection sets)")


def additive_energy_orbit():
    """All connection sets in the orbit have the same additive energy."""
    print(f"\n{'='*70}")
    print("ADDITIVE ENERGY IS ORBIT-INVARIANT")
    print("=" * 70)

    for p in [7, 11, 13, 19]:
        m = (p - 1) // 2
        S_int = list(range(1, m+1))

        def add_energy(S, p):
            S_set = set(S)
            count = 0
            for a in S:
                for b in S:
                    for c in S:
                        d = (a + b - c) % p
                        if d in S_set:
                            count += 1
            return count

        E_int = add_energy(S_int, p)

        # Check a few orbit elements
        energies = set()
        for a in range(1, p):
            S_scaled = sorted((a * s) % p for s in S_int)
            E = add_energy(S_scaled, p)
            energies.add(E)

        print(f"  p={p}: E(Int) = {E_int}, orbit energies = {energies}")
        if len(energies) == 1:
            print(f"  CONFIRMED: additive energy is constant on orbit")
        else:
            print(f"  SURPRISE: additive energy varies in orbit!")

    print("""
  This is expected: E(S) = #{(a,b,c,d) in S^4 : a+b=c+d}.
  Under scaling S -> aS: (a*x + a*y = a*z + a*w) iff (x+y = z+w),
  so E(aS) = E(S). Additive energy is scale-invariant.

  But NOT translation-invariant! E(S+t) can differ from E(S).
  The orbit maximizers are related by SCALING, not translation.
""")


def h_as_function_of_additive_energy():
    """Is H monotone in additive energy across ALL circulant tournaments?"""
    print(f"\n{'='*70}")
    print("H vs ADDITIVE ENERGY ACROSS ALL CIRCULANT TOURNAMENTS")
    print("=" * 70)

    # Use p=7 and p=11 (already have full data)
    for p in [7, 11]:
        m = (p - 1) // 2
        pairs = [(k+1, p-(k+1)) for k in range(m)]

        print(f"\n  p={p}, m={m}:")

        # Build all 2^m circulants and compute E(S) and H
        data = []
        for bits in range(2**m):
            S = []
            for k in range(m):
                if bits & (1 << k):
                    S.append(pairs[k][1])
                else:
                    S.append(pairs[k][0])
            S.sort()

            # Additive energy
            S_set = set(S)
            E = 0
            for a in S:
                for b in S:
                    for c in S:
                        d = (a + b - c) % p
                        if d in S_set:
                            E += 1

            data.append((E, S, bits))

        # Sort by E
        data.sort(key=lambda x: x[0], reverse=True)

        print(f"  {'E(S)':>6s}  {'S':>20s}  {'marker':>10s}")
        for E, S, bits in data:
            marker = ""
            S_int = list(range(1, m+1))
            S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
            if S == S_int:
                marker = "INTERVAL"
            elif S == S_qr and p % 4 == 3:
                marker = "PALEY"
            print(f"  {E:>6d}  {str(S):>20s}  {marker:>10s}")

        # Correlation between E and H?
        # We don't have H for all tournaments here, but we know:
        # - At p=7: Interval has highest E but NOT highest H (Paley wins)
        # - At p=19: Interval has highest E AND highest H
        # So additive energy alone doesn't determine the winner at small p


def main():
    analyze_walsh_circulant()
    orbit_structure()
    additive_energy_orbit()
    h_as_function_of_additive_energy()


if __name__ == '__main__':
    main()
