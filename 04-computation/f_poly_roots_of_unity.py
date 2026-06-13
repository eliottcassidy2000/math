#!/usr/bin/env python3
"""
f_poly_roots_of_unity.py — F(T,x) evaluated at roots of unity.

DISCOVERY: F(T,x) zeros cluster near angles ±2pi/3 on the unit circle.
This suggests evaluating F at:
  omega = e^{2pi*i/3} (cube root of unity)
  i = e^{pi*i/2} (4th root)
  zeta_k = e^{2pi*i/k} (k-th root)

F(T,1) = n! (always)
F(T,-1) = S(T) * (-1)^{n-1} (signed HP permanent)

What is F(T, omega)? F(T, i)?

For PALINDROMIC polynomials of degree d:
  F(x) = x^d F(1/x)
At x = omega (|omega|=1): F(omega) = omega^d * F(omega-bar) = omega^d * conj(F(omega))
So |F(omega)|^2 = F(omega) * conj(F(omega)) = F(omega) * omega^{-d} F(omega)
=> |F(omega)| = |F(omega)| (tautology)

But: F(omega) * omega^{-d} = conj(F(omega))
=> F(omega) = omega^d * conj(F(omega))
=> arg(F(omega)) = d * arg(omega) - arg(F(omega))  (mod 2pi)
=> 2*arg(F(omega)) = d * 2pi/3  (for omega = e^{2pi*i/3})

For d = n-1:
  If n-1 = 0 mod 3: arg(F(omega)) = 0, so F(omega) is REAL
  If n-1 = 1 mod 3: arg(F(omega)) = pi/3
  If n-1 = 2 mod 3: arg(F(omega)) = 2pi/3

Author: opus-2026-03-07-S44
"""
import numpy as np
from itertools import permutations
import math
import random

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

def eval_F_at(F, x):
    return sum(F[k] * x**k for k in range(len(F)))

# ============================================================
# F(T,x) AT ROOTS OF UNITY
# ============================================================
print("=" * 60)
print("F(T,x) AT ROOTS OF UNITY")
print("=" * 60)

omega = np.exp(2j * np.pi / 3)  # cube root
zeta4 = 1j  # 4th root
zeta5 = np.exp(2j * np.pi / 5)
zeta6 = np.exp(2j * np.pi / 6)

for n in [5, 7]:
    m = n*(n-1)//2
    print(f"\n--- n={n} (degree d={n-1}) ---")
    d = n - 1

    # Phase prediction for palindromic poly at omega
    predicted_phase = (d * 2 * np.pi / 3) / 2
    print(f"  Predicted arg(F(omega)) = {predicted_phase/np.pi:.4f}*pi")

    seen_F = set()

    if n == 5:
        iterator = range(1 << m)
    else:
        random.seed(42)
        iterator = [random.getrandbits(m) for _ in range(200)]

    results = []
    for bits in iterator:
        A = tournament_from_bits(bits, n)
        F = compute_F(A, n)
        key = tuple(F)
        if key in seen_F:
            continue
        seen_F.add(key)

        H = F[n-1]
        F_omega = eval_F_at(F, omega)
        F_i = eval_F_at(F, zeta4)
        F_neg1 = eval_F_at(F, -1)
        F_zeta5 = eval_F_at(F, zeta5)
        F_zeta6 = eval_F_at(F, zeta6)

        results.append({
            'H': H, 'F': F,
            'F_omega': F_omega, 'F_i': F_i,
            'F_neg1': F_neg1,
            'F_zeta5': F_zeta5, 'F_zeta6': F_zeta6
        })

    print(f"\n  {'H':>3} {'|F(w)|':>10} {'arg/pi':>8} {'|F(i)|':>10} {'arg/pi':>8} {'F(-1)':>8} {'|F(z5)|':>10} {'|F(z6)|':>10}")
    for r in sorted(results, key=lambda x: x['H']):
        fw = r['F_omega']
        fi = r['F_i']
        fz5 = r['F_zeta5']
        fz6 = r['F_zeta6']
        print(f"  {r['H']:3d} {abs(fw):10.2f} {np.angle(fw)/np.pi:+8.4f} "
              f"{abs(fi):10.2f} {np.angle(fi)/np.pi:+8.4f} "
              f"{r['F_neg1'].real:8.0f} {abs(fz5):10.2f} {abs(fz6):10.2f}")

# ============================================================
# F(T, omega) MODULAR ARITHMETIC
# ============================================================
print("\n" + "=" * 60)
print("F(T, omega) AS GAUSSIAN INTEGER")
print("=" * 60)
# omega = (-1 + i*sqrt(3))/2
# omega^0 = 1, omega^1 = omega, omega^2 = omega-bar
# F(T, omega) = sum F_k omega^k
# Since omega^3 = 1: F(T, omega) = sum_{j=0}^{2} omega^j * (sum_{k=j mod 3} F_k)

for n in [5, 7]:
    m = n*(n-1)//2
    print(f"\n--- n={n} ---")

    seen_F = set()
    if n == 5:
        iterator = range(1 << m)
    else:
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
        # Decompose into residues mod 3
        S0 = sum(F[k] for k in range(n) if k % 3 == 0)
        S1 = sum(F[k] for k in range(n) if k % 3 == 1)
        S2 = sum(F[k] for k in range(n) if k % 3 == 2)

        # F(omega) = S0 + S1*omega + S2*omega^2
        # = S0 + S1*(-1+i*sqrt(3))/2 + S2*(-1-i*sqrt(3))/2
        # = S0 - (S1+S2)/2 + i*sqrt(3)*(S1-S2)/2
        # Real part = S0 - (S1+S2)/2 = (2*S0 - S1 - S2)/2
        # Since S0 + S1 + S2 = n!: S1 + S2 = n! - S0
        # Real part = (2*S0 - (n! - S0))/2 = (3*S0 - n!)/2

        real_part = (3*S0 - math.factorial(n)) / 2
        imag_part = np.sqrt(3) * (S1 - S2) / 2

        # S0 - S1 + S2 - ... = F(-1) = S(T)*(-1)^{n-1}
        # S0 + S1 + S2 = n!
        # S0 - S1 + S2 = ? (relates to F(-1) differently)

        # Check: is |F(omega)|^2 always divisible by something nice?
        mod_sq = real_part**2 + imag_part**2

        # Also: F(omega) mod 3
        F_omega_mod3_real = int(round(real_part)) % 3
        F_omega_mod3_imag = int(round(imag_part / np.sqrt(3))) % 3  # imag/sqrt(3) should be integer

        print(f"  H={H:3d}: S0={S0:4d} S1={S1:4d} S2={S2:4d} "
              f"Re={real_part:8.1f} Im/√3={imag_part/np.sqrt(3):8.1f} "
              f"|F(w)|²={mod_sq:10.1f}")

# ============================================================
# NOVEL: F(T, q) FOR q = 2, 3 — COUNTING INTERPRETATIONS
# ============================================================
print("\n" + "=" * 60)
print("F(T,q) FOR SMALL INTEGERS q")
print("=" * 60)
# F(T,2) = sum_P 2^{fwd(P)}
# This weights paths by 2^{forward edges}.
# Since each edge is forward or backward:
# F(T,2) = sum_P 2^{fwd(P)} * 1^{bwd(P)}
# = sum_P prod (2 if forward, 1 if backward)
#
# F(T,q) = sum_P q^{fwd(P)} = sum_P prod (q if forward, 1 if backward)
# = sum_P prod_{i} (1 + (q-1)*A[P_i][P_{i+1}])

# At q=2: F(T,2) = sum_P prod (1 + A[P_i][P_{i+1}])
# Since A[P_i][P_{i+1}] in {0,1}:
# prod (1 + A) = prod 2^{A} = 2^{fwd(P)}

# F(T, 1+y) = sum_P (1+y)^{fwd(P)} = sum_P sum_k C(fwd,k) y^k
# = sum_k D_k y^k

# So F(T, 2) = sum D_k = sum_k sum_P C(fwd,k) = sum_P 2^{fwd(P)}

for n in [5, 7]:
    m = n*(n-1)//2
    print(f"\n--- n={n} ---")

    seen_F = set()
    if n == 5:
        iterator = range(1 << m)
    else:
        random.seed(42)
        iterator = [random.getrandbits(m) for _ in range(100)]

    vals = []
    for bits in iterator:
        A = tournament_from_bits(bits, n)
        F = compute_F(A, n)
        key = tuple(F)
        if key in seen_F:
            continue
        seen_F.add(key)

        H = F[n-1]
        F2 = sum(F[k] * 2**k for k in range(n))
        F3 = sum(F[k] * 3**k for k in range(n))
        Fneg2 = sum(F[k] * (-2)**k for k in range(n))
        vals.append((H, F2, F3, Fneg2))

    vals.sort()
    print(f"  {'H':>3} {'F(2)':>10} {'F(3)':>10} {'F(-2)':>10} {'F(2)/2^(n-1)':>14} {'F(3)/3^(n-1)':>14}")
    for H, F2, F3, Fn2 in vals[:15]:
        print(f"  {H:3d} {F2:10d} {F3:10d} {Fn2:10d} {F2/2**(n-1):14.4f} {F3/3**(n-1):14.4f}")

# ============================================================
# NOVEL: F(T, q) / q^{(n-1)/2} AS "CENTERED" VALUE
# ============================================================
print("\n" + "=" * 60)
print("F(T,q) / q^{(n-1)/2} — CENTERED NORMALIZATION")
print("=" * 60)
# Since F is palindromic of degree n-1, the "natural" normalization is
# G(T,q) = F(T,q) / q^{(n-1)/2}
# G(T,1) = n!
# G(T,-1) = F(-1)/(-1)^{(n-1)/2}
# G is self-reciprocal: G(q) = G(1/q)

for n in [5, 7]:
    print(f"\n  n={n}:")
    m = n*(n-1)//2
    seen_F = set()
    if n == 5:
        iterator = range(1 << m)
    else:
        random.seed(42)
        iterator = [random.getrandbits(m) for _ in range(100)]

    for bits in iterator:
        A = tournament_from_bits(bits, n)
        F = compute_F(A, n)
        key = tuple(F)
        if key in seen_F:
            continue
        seen_F.add(key)

        H = F[n-1]
        center = (n-1) / 2
        # G(q) = q^{-center} F(q) = sum F_k q^{k - center}
        # G(2) = 2^{-center} F(2)
        G2 = sum(F[k] * 2**(k - center) for k in range(n))
        G3 = sum(F[k] * 3**(k - center) for k in range(n))

        print(f"    H={H:3d}: G(2)={G2:10.2f}, G(3)={G3:10.2f}, G(2)/n!={G2/math.factorial(n):.6f}")
