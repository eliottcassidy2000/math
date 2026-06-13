#!/usr/bin/env python3
"""
f_poly_mod3_structure.py — Deep analysis of F(T,omega) structure.

DISCOVERY: At n=7, F(T,omega) is ALWAYS REAL (S1 = S2).
Also F(T,omega) is ALWAYS divisible by 9.

Why is S1 = S2? The palindrome F_k = F_{n-1-k} gives:
  S_j = sum_{k=j mod 3} F_k

For n=7 (degree 6):
  S0 = F_0 + F_3 + F_6
  S1 = F_1 + F_4
  S2 = F_2 + F_5

Palindrome: F_0 = F_6, F_1 = F_5, F_2 = F_4, F_3 = F_3.
So S0 = 2*F_0 + F_3, S1 = F_1 + F_4 = F_5 + F_2 = S2. YES!

The palindrome FORCES S1 = S2 when n-1 ≡ 0 mod 3 (i.e., n ≡ 1 mod 3).
For n=7: n-1=6, pairs (1,5) and (2,4) are swapped by palindrome,
and they fall into residue classes 1 and 2 mod 3 respectively... wait:
  k=1 → class 1, k=5 → class 2. But F_1 = F_5.
  k=2 → class 2, k=4 → class 1. But F_2 = F_4.
So S1 = F_1 + F_4 = F_5 + F_2 = S2. ✓

For n=5: n-1=4, pairs (0,4) and (1,3), F_2 unpaired.
  k=0 → class 0, k=4 → class 1. F_0 = F_4.
  k=1 → class 1, k=3 → class 0. F_1 = F_3.
  S0 = F_0 + F_3 = F_4 + F_1 = S1. So S0 = S1 at n=5!
  S2 = F_2 (alone).
  This means Re = (3S0 - n!)/2 = (3(F_0+F_3) - 120)/2
  Im/√3 = (S1-S2)/2 = (F_0+F_3-F_2)/2 = (S0-F_2)/2

GENERAL RULE: The palindrome k ↦ n-1-k induces a permutation on
residue classes mod 3. If n-1 ≡ 0 mod 3: classes {0,1,2} map to {0,2,1}.
So S1 = S2 always. F(omega) is real.

If n-1 ≡ 1 mod 3: k → n-1-k sends class j to class (n-1-j) mod 3 = (1-j) mod 3.
  0 → 1, 1 → 0, 2 → 2. So S0 = S1, S2 fixed. F(omega) has Im/√3 = (S0-S2)/2 ≠ 0 generally.

If n-1 ≡ 2 mod 3: k → n-1-k sends class j to (2-j) mod 3.
  0 → 2, 1 → 1, 2 → 0. So S0 = S2, S1 fixed. F(omega) has Im/√3 = (S1-S0)/2 ≠ 0 generally.
  Wait: actually Im/√3 = (S1-S2)/2 = (S1-S0)/2 ≠ 0.

So F(T,omega) is REAL iff n ≡ 1 mod 3 (i.e., n-1 ≡ 0 mod 3).
F(T,omega) has S0=S1 iff n ≡ 2 mod 3 (n-1 ≡ 1 mod 3).
F(T,omega) has S0=S2 iff n ≡ 0 mod 3 (n-1 ≡ 2 mod 3).

This is a THEOREM about palindromic polynomials at roots of unity!

Now: is F(T,omega) always divisible by 9 at n=7?
F(omega) = (3*S0 - 5040)/2. We need 3*S0 ≡ 5040 mod 18.
5040 = 7! = 2^4 * 3^2 * 5 * 7. 5040/2 = 2520, 5040 mod 18 = 0.
So F(omega) = (3*S0)/2 - 2520. Need 3*S0/2 ≡ 0 mod 9, i.e., S0 ≡ 0 mod 6.

Actually let me just check the divisibility pattern more carefully.

Author: opus-2026-03-07-S44
"""
from itertools import permutations, combinations
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

# ============================================================
# THEOREM: PALINDROME + ROOTS OF UNITY
# ============================================================
print("=" * 60)
print("THEOREM: F(T, omega) REALITY FROM PALINDROME")
print("=" * 60)

for n in [3, 4, 5, 6, 7]:
    d = n - 1
    res = d % 3
    if res == 0:
        prediction = "F(omega) ALWAYS REAL (S1=S2)"
    elif res == 1:
        prediction = "S0=S1, F(omega) generally complex"
    else:
        prediction = "S0=S2, F(omega) generally complex"
    print(f"  n={n}: d={d}, d mod 3 = {res}: {prediction}")

# ============================================================
# DIVISIBILITY OF F(T, omega) BY 9
# ============================================================
print("\n" + "=" * 60)
print("F(T, omega) DIVISIBILITY AT n=7")
print("=" * 60)

n = 7
m = n*(n-1)//2
random.seed(42)
seen_F = set()
f_omega_vals = set()

for trial in range(200):
    bits = random.getrandbits(m)
    A = tournament_from_bits(bits, n)
    F = compute_F(A, n)
    key = tuple(F)
    if key in seen_F:
        continue
    seen_F.add(key)

    S0 = sum(F[k] for k in range(n) if k % 3 == 0)
    S1 = sum(F[k] for k in range(n) if k % 3 == 1)
    S2 = sum(F[k] for k in range(n) if k % 3 == 2)

    # F(omega) = (3*S0 - n!)/2 (real)
    f_omega = (3 * S0 - math.factorial(n)) // 2
    f_omega_vals.add(f_omega)

# Check GCD of all F(omega) values
from math import gcd
from functools import reduce
all_vals = list(f_omega_vals)
g = reduce(gcd, [abs(v) for v in all_vals if v != 0])
print(f"  n=7: {len(all_vals)} distinct F(omega) values")
print(f"  GCD of |F(omega)|: {g}")
print(f"  All F(omega) mod 9: {sorted(set(v % 9 for v in all_vals))}")
print(f"  All F(omega) mod 3: {sorted(set(v % 3 for v in all_vals))}")
print(f"  All F(omega) / 9: {sorted(set(v // 9 for v in all_vals))[:20]}...")

# ============================================================
# F(T, i) STRUCTURE
# ============================================================
print("\n" + "=" * 60)
print("F(T, i) STRUCTURE — EVALUATION AT 4TH ROOT OF UNITY")
print("=" * 60)
# F(T, i) = sum F_k i^k = S0 - S2 + i(S1 - S3) where S_j = sum_{k=j mod 4} F_k
# For n=7: d=6
#   S0(mod4) = F_0 + F_4
#   S1(mod4) = F_1 + F_5
#   S2(mod4) = F_2 + F_6
#   S3(mod4) = F_3
# Palindrome: F_0=F_6, F_1=F_5, F_2=F_4
# So: S0 = F_0 + F_4 = F_6 + F_2 = S2. ← palindrome forces S0 = S2!
# And: S1 = F_1 + F_5 = 2*F_1. S3 = F_3.
# So F(i) = (S0 - S2) + i*(S1 - S3) = 0 + i*(2*F_1 - F_3) = i*(2*F_1 - F_3)
# F(i) is PURE IMAGINARY at n=7!

for n in [5, 7]:
    d = n - 1
    print(f"\n  n={n}, d={d}:")
    # Which S_j's are equal?
    # k ↦ d-k sends class j to (d-j) mod 4
    for j in range(4):
        partner = (d - j) % 4
        if j < partner:
            print(f"    S_{j} = S_{partner} (palindrome pairing)")
        elif j == partner:
            print(f"    S_{j} is self-paired")

    m = n*(n-1)//2
    seen = set()
    if n == 5:
        iterator = range(1 << m)
    else:
        random.seed(42)
        iterator = [random.getrandbits(m) for _ in range(200)]

    fi_vals = []
    for bits in iterator:
        A = tournament_from_bits(bits, n)
        F = compute_F(A, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)

        H = F[n-1]
        # Compute S_j mod 4
        S = [sum(F[k] for k in range(n) if k % 4 == j) for j in range(4)]
        Re_fi = S[0] - S[2]
        Im_fi = S[1] - S[3]

        fi_vals.append((H, Re_fi, Im_fi, S))

    # Check predictions
    fi_re_set = set(v[1] for v in fi_vals)
    fi_im_set = set(v[2] for v in fi_vals)
    print(f"    Re(F(i)) values: {sorted(fi_re_set)[:10]}")
    print(f"    Im(F(i)) values: {sorted(fi_im_set)[:10]}...")

    if len(fi_re_set) == 1 and 0 in fi_re_set:
        print(f"    *** Re(F(i)) = 0 ALWAYS! F(i) is PURE IMAGINARY! ***")
    if len(fi_im_set) == 1 and 0 in fi_im_set:
        print(f"    *** Im(F(i)) = 0 ALWAYS! F(i) is REAL! ***")

    # GCD of Im(F(i))
    if fi_im_set != {0}:
        g = reduce(gcd, [abs(v) for v in fi_im_set if v != 0])
        print(f"    GCD of |Im(F(i))|: {g}")
        print(f"    Im(F(i)) / {g}: {sorted(set(v // g for v in fi_im_set))[:20]}")

# ============================================================
# GENERAL: F(T, zeta_k) phase structure for all k
# ============================================================
print("\n" + "=" * 60)
print("PALINDROME PHASE THEOREM — GENERAL k-TH ROOT")
print("=" * 60)
print("For palindromic polynomial F of degree d:")
print("F(zeta_k) = zeta_k^d * conj(F(zeta_k))")
print("")
print("This means: arg(F(zeta_k)) = d*pi/k + n*pi for some integer n")
print("So F(zeta_k) lies on a FIXED RAY in the complex plane!")
print("")

for n in [5, 7]:
    d = n - 1
    print(f"n={n}, d={d}:")
    for k in [2, 3, 4, 5, 6, 8]:
        # zeta_k^d = e^{2*pi*i*d/k}
        phase = (2 * d / k) % 2  # in units of pi
        # F(zeta_k) = e^{i*pi*phase} * conj(F(zeta_k))
        # Let F = r * e^{i*theta}. Then r*e^{i*theta} = e^{i*pi*phase} * r*e^{-i*theta}
        # => e^{2*i*theta} = e^{i*pi*phase} => theta = pi*phase/2 (mod pi)
        ray_angle = phase / 2  # in units of pi
        if abs(ray_angle) < 1e-10 or abs(ray_angle - 1) < 1e-10:
            label = "REAL"
        elif abs(ray_angle - 0.5) < 1e-10 or abs(ray_angle - 1.5) < 1e-10:
            label = "PURE IMAGINARY"
        else:
            label = f"ray at {ray_angle:.4f}*pi"
        print(f"  zeta_{k}: F lies on {label}")

# ============================================================
# F(T, omega) mod p for small primes
# ============================================================
print("\n" + "=" * 60)
print("F(T, omega) MOD SMALL PRIMES (n=7)")
print("=" * 60)

n = 7
m = n*(n-1)//2
random.seed(42)
seen = set()

f_omega_list = []
for trial in range(200):
    bits = random.getrandbits(m)
    A = tournament_from_bits(bits, n)
    F = compute_F(A, n)
    key = tuple(F)
    if key in seen:
        continue
    seen.add(key)
    S0 = sum(F[k] for k in range(n) if k % 3 == 0)
    f_omega = (3 * S0 - math.factorial(n)) // 2
    f_omega_list.append(f_omega)

for p in [2, 3, 5, 7, 9, 27]:
    residues = sorted(set(v % p for v in f_omega_list))
    print(f"  F(omega) mod {p:2d}: {residues}")

# ============================================================
# CONNECTION: F(omega) and 3-CYCLES
# ============================================================
print("\n" + "=" * 60)
print("F(T, omega) vs NUMBER OF 3-CYCLES t3")
print("=" * 60)

def count_t3(A, n):
    t3 = 0
    for triple in combinations(range(n), 3):
        i, j, k = triple
        if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]):
            t3 += 1
    return t3

n = 7
m = n*(n-1)//2
random.seed(42)
seen = set()

t3_fomega = {}
for trial in range(200):
    bits = random.getrandbits(m)
    A = tournament_from_bits(bits, n)
    F = compute_F(A, n)
    key = tuple(F)
    if key in seen:
        continue
    seen.add(key)

    t3 = count_t3(A, n)
    S0 = sum(F[k] for k in range(n) if k % 3 == 0)
    f_omega = (3 * S0 - math.factorial(n)) // 2

    if t3 not in t3_fomega:
        t3_fomega[t3] = set()
    t3_fomega[t3].add(f_omega)

print(f"  t3 → F(omega) values:")
for t3 in sorted(t3_fomega.keys()):
    vals = sorted(t3_fomega[t3])
    print(f"    t3={t3:2d}: {vals[:10]}{'...' if len(vals) > 10 else ''} ({len(vals)} distinct)")

# Is F(omega) determined by t3? No, since multiple values appear.
# But is F(omega) mod 9 determined by t3?
print(f"\n  t3 → F(omega) mod 9:")
for t3 in sorted(t3_fomega.keys()):
    vals = sorted(t3_fomega[t3])
    mod9 = sorted(set(v % 9 for v in vals))
    print(f"    t3={t3:2d}: F(omega) mod 9 = {mod9}")

# Is F(omega) / 9 related to anything?
print(f"\n  F(omega) / 9 = (S0 - 1680) / 6")
print(f"  (Since F(omega) = (3*S0 - 5040)/2 = 9*(S0/6 - 280))")
print(f"  S0 = F_0 + F_3 + F_6 = 2*H + F_3 (palindrome)")
print(f"  So F(omega)/9 = (2*H + F_3 - 1680)/6")
