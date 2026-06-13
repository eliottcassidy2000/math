#!/usr/bin/env python3
"""
f_poly_zeros_leeyang.py — Zeros of F(T,x) in the complex plane.

LEE-YANG THEORY: For a partition function Z(x), the distribution of zeros
in the complex x-plane reveals phase transitions. The Lee-Yang theorem
says that for ferromagnetic Ising models, all zeros lie on |x| = 1.

F(T,x) = sum_P x^{fwd(P)} is a partition function with "activity" x
weighting forward edges. Its zeros might reveal:
1. Whether F has a Lee-Yang-like property (zeros on a circle?)
2. The angular distribution of zeros
3. Connections to the hard-core lattice gas (I(Omega,x) zeros)

ALSO: Compare F(T,x) zeros with I(Omega(T),x) zeros.
I(Omega,x) has all real negative zeros for n<=8 (Chudnovsky-Seymour).
F(T,x) is NOT real-rooted in general. Where do the complex zeros go?

PALINDROME: F(T,x) = x^{n-1} F(T, 1/x), so if r is a root, so is 1/r.
The zeros come in reciprocal pairs! This means they lie on a "reciprocal curve."
For zeros on the unit circle |r| = 1, the reciprocal is the conjugate,
so unit-circle zeros come in conjugate pairs (automatic for real coefficients).

Author: opus-2026-03-07-S44
"""
import numpy as np
from itertools import permutations, combinations
import math

def tournament_from_bits(bits, n):
    A = [[0]*n for _ in range(n)]
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> pos) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            pos += 1
    return A

def compute_F(A, n):
    F = [0] * n
    for P in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if A[P[i]][P[i+1]])
        F[fwd] += 1
    return F

def find_odd_cycles(A, n):
    """Find all odd directed cycles."""
    cycles = []
    for k in [3, 5, 7]:
        if k > n:
            break
        for combo in combinations(range(n), k):
            for perm in permutations(combo):
                if all(A[perm[i]][perm[(i+1) % k]] for i in range(k)):
                    canon = min(perm[i:] + perm[:i] for i in range(k))
                    if tuple(canon) not in [tuple(c) for c in cycles]:
                        cycles.append(list(canon))
                    break
    return cycles

def independence_poly(cycles, n_cycles):
    """Compute I(Omega, x) coefficients."""
    if n_cycles == 0:
        return [1]

    # Build adjacency (share vertex = edge)
    adj = [[False]*n_cycles for _ in range(n_cycles)]
    for i in range(n_cycles):
        for j in range(i+1, n_cycles):
            if set(cycles[i]) & set(cycles[j]):
                adj[i][j] = adj[j][i] = True

    # Count independent sets by size
    alpha = [0] * (n_cycles + 1)
    for size in range(n_cycles + 1):
        for S in combinations(range(n_cycles), size):
            independent = True
            for a, b in combinations(S, 2):
                if adj[a][b]:
                    independent = False
                    break
            if independent:
                alpha[size] += 1

    return alpha

# ============================================================
# ZEROS OF F(T,x) IN THE COMPLEX PLANE
# ============================================================
print("=" * 60)
print("ZEROS OF F(T,x) — LEE-YANG ANALYSIS")
print("=" * 60)

for n in [5, 7]:
    m = n*(n-1)//2
    print(f"\n--- n={n} ---")

    seen_F = set()
    all_zeros_data = []

    if n == 5:
        iterator = range(1 << m)
    else:
        import random
        random.seed(42)
        iterator = [random.getrandbits(m) for _ in range(200)]

    for bits in iterator:
        A = tournament_from_bits(bits, n)
        F = compute_F(A, n)
        key = tuple(F)
        if key in seen_F:
            continue
        seen_F.add(key)

        H = F[n-1]
        # Polynomial coefficients (highest degree first for np.roots)
        poly = [F[k] for k in range(len(F)-1, -1, -1)]
        roots = np.roots(poly)

        # Classify roots
        on_unit_circle = sum(1 for r in roots if abs(abs(r) - 1) < 1e-6)
        real_negative = sum(1 for r in roots if abs(r.imag) < 1e-6 and r.real < -1e-6)
        all_real = all(abs(r.imag) < 1e-6 for r in roots)
        moduli = sorted([abs(r) for r in roots])

        all_zeros_data.append({
            'H': H, 'F': F, 'roots': roots,
            'on_unit': on_unit_circle,
            'real_neg': real_negative,
            'all_real': all_real,
            'moduli': moduli
        })

    # Summary statistics
    unit_circle_counts = {}
    for d in all_zeros_data:
        k = d['on_unit']
        unit_circle_counts[k] = unit_circle_counts.get(k, 0) + 1

    print(f"  Total distinct F polynomials: {len(all_zeros_data)}")
    print(f"  Zeros on unit circle distribution: {dict(sorted(unit_circle_counts.items()))}")
    print(f"  All-real-rooted: {sum(1 for d in all_zeros_data if d['all_real'])}/{len(all_zeros_data)}")

    # Check: are ALL zeros on the unit circle for any tournament?
    all_unit = [d for d in all_zeros_data if d['on_unit'] == n-1]
    print(f"  ALL zeros on unit circle: {len(all_unit)}/{len(all_zeros_data)}")

    # Check reciprocal pairing: if r is root, is 1/r also root?
    recip_ok = 0
    for d in all_zeros_data:
        roots = d['roots']
        paired = True
        for r in roots:
            if abs(r) < 1e-10:
                paired = False
                break
            # Check if 1/r is also a root
            dists = [abs(r2 - 1/r) for r2 in roots]
            if min(dists) > 1e-4:
                paired = False
                break
        if paired:
            recip_ok += 1
    print(f"  Reciprocal pairing verified: {recip_ok}/{len(all_zeros_data)}")

    # Show a few examples
    print(f"\n  Examples:")
    for d in all_zeros_data[:5]:
        root_strs = [f"{r.real:+.4f}{r.imag:+.4f}i" for r in sorted(d['roots'], key=lambda z: (abs(z), z.real))]
        moduli_str = [f"{abs(r):.4f}" for r in sorted(d['roots'], key=lambda z: abs(z))]
        print(f"    H={d['H']:3d} F={d['F']} unit={d['on_unit']}")
        print(f"      roots: {root_strs}")
        print(f"      |roots|: {moduli_str}")

# ============================================================
# ANGULAR DISTRIBUTION ON UNIT CIRCLE
# ============================================================
print("\n" + "=" * 60)
print("ANGULAR DISTRIBUTION OF ZEROS")
print("=" * 60)

n = 7
m = n*(n-1)//2
random.seed(42)
seen_F = set()
all_angles = []

for trial in range(200):
    bits = random.getrandbits(m)
    A = tournament_from_bits(bits, n)
    F = compute_F(A, n)
    key = tuple(F)
    if key in seen_F:
        continue
    seen_F.add(key)

    poly = [F[k] for k in range(len(F)-1, -1, -1)]
    roots = np.roots(poly)

    for r in roots:
        angle = np.angle(r) / np.pi  # in units of pi
        all_angles.append((abs(r), angle))

# Histogram of moduli
moduli = [m for m, a in all_angles]
print(f"  Total zeros: {len(all_angles)}")
print(f"  Modulus statistics: min={min(moduli):.4f}, max={max(moduli):.4f}, mean={np.mean(moduli):.4f}")

# How many are near unit circle?
near_unit = sum(1 for m in moduli if abs(m - 1) < 0.1)
print(f"  Near unit circle (|r| in [0.9, 1.1]): {near_unit}/{len(all_angles)} = {near_unit/len(all_angles)*100:.1f}%")

# Angle histogram (bins of pi/6)
angle_vals = [a for m, a in all_angles if abs(m - 1) < 0.2]
if angle_vals:
    hist, edges = np.histogram(angle_vals, bins=12, range=(-1, 1))
    print(f"  Angle histogram (near unit circle, bins of pi/6):")
    for i in range(len(hist)):
        bar = '#' * (hist[i] // 2)
        print(f"    [{edges[i]:+.2f}pi, {edges[i+1]:+.2f}pi): {hist[i]:3d} {bar}")

# ============================================================
# COMPARE F(T,x) ZEROS WITH I(Omega,x) ZEROS
# ============================================================
print("\n" + "=" * 60)
print("F(T,x) ZEROS vs I(Omega(T),x) ZEROS")
print("=" * 60)

n = 5
m = n*(n-1)//2
seen = set()

for bits in range(1 << m):
    A = tournament_from_bits(bits, n)
    F = compute_F(A, n)
    key = tuple(F)
    if key in seen:
        continue
    seen.add(key)

    H = F[n-1]

    # F(T,x) zeros
    poly_F = [F[k] for k in range(len(F)-1, -1, -1)]
    roots_F = np.roots(poly_F)

    # I(Omega, x) zeros
    cycles = find_odd_cycles(A, n)
    alpha = independence_poly(cycles, len(cycles))
    # I(Omega, x) = sum alpha_k x^k
    if len(alpha) > 1:
        poly_I = list(reversed(alpha))
        roots_I = np.roots(poly_I)
        roots_I_str = [f"{r.real:.3f}" for r in sorted(roots_I, key=lambda z: z.real)]
    else:
        roots_I_str = ["(no zeros)"]

    roots_F_str = [f"{r.real:.3f}{r.imag:+.3f}i" for r in sorted(roots_F, key=lambda z: z.real)]

    print(f"  H={H}: F zeros={roots_F_str}")
    print(f"         I zeros={roots_I_str}")

# ============================================================
# KEY OBSERVATION: F(T, -1) = S(T) * (-1)^{n-1}
# So x = -1 is special. The behavior of F near x = -1 matters.
# ============================================================
print("\n" + "=" * 60)
print("F(T,x) NEAR x = -1 (SIGNED PERMANENT CONNECTION)")
print("=" * 60)

n = 5
m = n*(n-1)//2
seen = set()
print(f"  n={n}: F(-1) = (-1)^{n-1} * S(T)")

for bits in range(1 << m):
    A = tournament_from_bits(bits, n)
    F = compute_F(A, n)
    key = tuple(F)
    if key in seen:
        continue
    seen.add(key)

    H = F[n-1]
    F_at_neg1 = sum(F[k] * (-1)**k for k in range(n))
    S = (-1)**(n-1) * F_at_neg1

    # F'(-1) — derivative at x=-1
    F_prime = [k * F[k] * (-1)**(k-1) for k in range(1, n)]
    F_deriv_at_neg1 = sum(F_prime)

    print(f"  H={H:3d}: F(-1)={F_at_neg1:6d}, S={S:6d}, F'(-1)={F_deriv_at_neg1:6d}")

# ============================================================
# NOVEL: MAHLER MEASURE OF F(T,x)
# ============================================================
print("\n" + "=" * 60)
print("MAHLER MEASURE OF F(T,x)")
print("=" * 60)
# M(f) = |a_d| * prod max(1, |r_i|) where r_i are roots
# For palindromic polynomials with reciprocal roots, M(f) has special properties.
# Lehmer's conjecture: smallest Mahler measure > 1 is Lehmer's number ~1.17628

for n in [5, 7]:
    m = n*(n-1)//2
    print(f"\n  n={n}:")
    seen = set()

    if n == 5:
        iterator = range(1 << m)
    else:
        random.seed(42)
        iterator = [random.getrandbits(m) for _ in range(200)]

    mahler_vals = []
    for bits in iterator:
        A = tournament_from_bits(bits, n)
        F = compute_F(A, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)

        H = F[n-1]
        poly = [F[k] for k in range(len(F)-1, -1, -1)]
        roots = np.roots(poly)

        # Mahler measure: |leading coeff| * prod max(1, |r|)
        leading = F[n-1]
        mahler = abs(leading) * np.prod([max(1, abs(r)) for r in roots])

        # Normalized: divide by H to get "per-path" Mahler measure
        mahler_norm = mahler / H

        mahler_vals.append((H, mahler, mahler_norm))

    mahler_vals.sort(key=lambda x: x[1])
    for H, mah, mah_n in mahler_vals[:8]:
        print(f"    H={H:3d}: Mahler={mah:12.2f}, Mahler/H={mah_n:10.4f}")
    if len(mahler_vals) > 8:
        print(f"    ... ({len(mahler_vals)} total)")

# ============================================================
# NOVEL: RESULTANT OF F(T,x) AND I(Omega,x)
# ============================================================
print("\n" + "=" * 60)
print("RESULTANT: Res(F(T,x), I(Omega(T),x))")
print("=" * 60)
# The resultant measures whether two polynomials share a common root.
# If Res = 0, they share a root. This would mean F and I have a common zero —
# a point where both the partition function and the independence polynomial vanish.

n = 5
m = n*(n-1)//2
seen = set()

for bits in range(1 << m):
    A = tournament_from_bits(bits, n)
    F = compute_F(A, n)
    key = tuple(F)
    if key in seen:
        continue
    seen.add(key)

    H = F[n-1]
    cycles = find_odd_cycles(A, n)
    alpha = independence_poly(cycles, len(cycles))

    if len(alpha) <= 1:
        print(f"  H={H}: I(Omega,x) = 1 (no cycles), no shared roots possible")
        continue

    # Compute resultant via numpy
    poly_F = np.array([F[k] for k in range(len(F)-1, -1, -1)], dtype=float)
    poly_I = np.array(list(reversed(alpha)), dtype=float)

    # Resultant = product of F evaluated at roots of I (times leading coeff powers)
    roots_I = np.roots(poly_I)
    F_at_I_roots = [np.polyval(poly_F, r) for r in roots_I]

    # Product
    res_approx = abs(np.prod(F_at_I_roots)) * abs(poly_F[0])**len(roots_I)

    print(f"  H={H}: |Res(F,I)| ≈ {res_approx:.2f}, I_zeros={[f'{r:.3f}' for r in roots_I]}")
    print(f"         F at I-zeros: {[f'{abs(v):.2f}' for v in F_at_I_roots]}")
