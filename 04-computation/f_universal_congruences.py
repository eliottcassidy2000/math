#!/usr/bin/env python3
"""
f_universal_congruences.py — Systematic search for universal congruences of F(T,x).

KNOWN: F(T,1) = n! (trivially universal)
       F(T,-1) = S(T)*(-1)^{n-1}, S mod 2^{n-1} universal
       F(T,omega) ≡ 0 mod 9 at n=7

QUESTION: For which evaluation points x and moduli m is F(T,x) mod m universal?

The D_k = [y^k] F(T, 1+y) have D_k mod 2^{n-1-k} universal.
What about F(T, x) for x = 2, 3, -2, etc.?

F(T, 2) = sum F_k * 2^k = sum D_k (by binomial expansion of (1+1)^k)
F(T, 3) = sum F_k * 3^k = sum D_k * 2^k (since 3 = 1 + 2)

So F(T, 1+y) = sum D_k y^k where D_k mod 2^{n-1-k} is universal.
F(T, 1+y) mod 2^{n-1} = D_0 mod 2^{n-1} = n! mod 2^{n-1} (universal)

More interesting: what modulus works for y = 1 (x = 2)?
F(T, 2) = sum D_k = sum D_k
D_0 = n!, D_1 = C(n,2)*(n-1)!, both universal mod huge powers of 2.
D_2 = 16800 + 240*t3 at n=7, so D_2 mod 240 not universal.
But D_2 mod 16 = (16800 + 240*t3) mod 16 = 0 + 0 = 0 (universal mod 16).

F(2) = sum D_k. The universality of F(2) mod M requires each D_k to be
universal mod M. Since D_k is universal mod 2^{n-1-k}, we need
2^{n-1-k} | M for each k. The weakest constraint is k = n-2: D_{n-2} mod 2.
So F(2) mod 2 should be universal. Is F(2) mod higher powers universal?

Author: opus-2026-03-07-S44
"""
from itertools import permutations, combinations
import math
import random
from functools import reduce
from math import gcd

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
# SYSTEMATIC CONGRUENCE SEARCH
# ============================================================
print("=" * 60)
print("UNIVERSAL CONGRUENCES OF F(T,x) — SYSTEMATIC SEARCH")
print("=" * 60)

for n in [5, 7]:
    m = n*(n-1)//2
    print(f"\n=== n={n} ===")

    # Collect all distinct F polynomials
    seen = set()
    all_F = []

    if n == 5:
        iterator = range(1 << m)
    else:
        random.seed(42)
        iterator = [random.getrandbits(m) for _ in range(300)]

    for bits in iterator:
        A = tournament_from_bits(bits, n)
        F = compute_F(A, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)
        all_F.append(F)

    print(f"  {len(all_F)} distinct F polynomials")

    # For each evaluation point x, find max universal modulus
    for x_val in [-3, -2, -1, 0, 1, 2, 3, 4, 5]:
        vals = [sum(F[k] * x_val**k for k in range(n)) for F in all_F]

        # Find GCD of all pairwise differences
        diffs = set()
        for i in range(len(vals)):
            for j in range(i+1, len(vals)):
                diffs.add(abs(vals[i] - vals[j]))
        diffs.discard(0)

        if not diffs:
            print(f"  F({x_val:+d}): ALL IDENTICAL = {vals[0]}")
        else:
            g = reduce(gcd, diffs)
            residue = vals[0] % g
            print(f"  F({x_val:+d}): universal mod {g}, residue = {residue}, "
                  f"range [{min(vals)}, {max(vals)}]")

    # Also check at omega (cube root of unity)
    # F(omega) = (3*S0 - n!)/2 where S0 = sum F_k for k ≡ 0 mod 3
    omega_vals = []
    for F in all_F:
        S0 = sum(F[k] for k in range(n) if k % 3 == 0)
        f_omega = 3 * S0 - math.factorial(n)  # = 2 * F(omega)
        omega_vals.append(f_omega)

    diffs = set()
    for i in range(len(omega_vals)):
        for j in range(i+1, len(omega_vals)):
            diffs.add(abs(omega_vals[i] - omega_vals[j]))
    diffs.discard(0)
    if diffs:
        g = reduce(gcd, diffs)
        print(f"  2*F(omega): universal mod {g}, residue = {omega_vals[0] % g}")
        # So F(omega) is universal mod g/2 (if g even) or g (if odd)
        if g % 2 == 0:
            print(f"    => F(omega) universal mod {g//2}")

    # F(i) = Re + i*Im
    # Re = S0(mod4) - S2(mod4), Im = S1(mod4) - S3(mod4)
    i_re_vals = []
    i_im_vals = []
    for F in all_F:
        S = [sum(F[k] for k in range(n) if k % 4 == j) for j in range(4)]
        i_re_vals.append(S[0] - S[2])
        i_im_vals.append(S[1] - S[3])

    for label, vals in [("Re(F(i))", i_re_vals), ("Im(F(i))", i_im_vals)]:
        diffs = set()
        for i in range(len(vals)):
            for j in range(i+1, len(vals)):
                diffs.add(abs(vals[i] - vals[j]))
        diffs.discard(0)
        if not diffs:
            print(f"  {label}: ALL IDENTICAL = {vals[0]}")
        else:
            g = reduce(gcd, diffs)
            print(f"  {label}: universal mod {g}, residue = {vals[0] % g}")

# ============================================================
# D_k UNIVERSALITY TABLE
# ============================================================
print("\n" + "=" * 60)
print("D_k UNIVERSALITY TABLE")
print("=" * 60)

for n in [5, 7]:
    m = n*(n-1)//2
    print(f"\n  n={n}:")

    seen = set()
    all_F = []
    if n == 5:
        iterator = range(1 << m)
    else:
        random.seed(42)
        iterator = [random.getrandbits(m) for _ in range(300)]

    for bits in iterator:
        A = tournament_from_bits(bits, n)
        F = compute_F(A, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)
        all_F.append(F)

    for k in range(n):
        # D_k = sum_P C(fwd(P), k)
        dk_vals = []
        for F in all_F:
            dk = sum(F[j] * math.comb(j, k) for j in range(n))
            dk_vals.append(dk)

        diffs = set()
        for i in range(len(dk_vals)):
            for j in range(i+1, len(dk_vals)):
                diffs.add(abs(dk_vals[i] - dk_vals[j]))
        diffs.discard(0)

        if not diffs:
            print(f"    D_{k}: UNIVERSAL = {dk_vals[0]}")
        else:
            g = reduce(gcd, diffs)
            v2 = 0
            gg = g
            while gg % 2 == 0:
                v2 += 1
                gg //= 2
            predicted = 2**(n-1-k)
            print(f"    D_{k}: mod {g} (2^{v2}*{gg}), predicted 2^{n-1-k}={predicted}, "
                  f"ratio={g/predicted:.1f}")

# ============================================================
# NOVEL: F(T, x) mod (x-1)^k — TAYLOR COEFFICIENTS AT x=1
# ============================================================
print("\n" + "=" * 60)
print("TAYLOR EXPANSION OF F(T,x) AT x=1")
print("=" * 60)
# F(T, 1+y) = sum D_k y^k
# Already done via D_k. But let's think about this differently.
# The coefficients a_k of the expansion F(T, 1+y) = sum a_k y^k
# are EXACTLY the D_k. And D_k is the k-th "moment" C(fwd, k).

# NEW: What about expansion at x = -1?
# F(T, -1+y) = sum E_k y^k where E_k = sum F_j * C(j, k) * (-1)^{j-k}
# E_0 = F(-1) = S * (-1)^{n-1}
# E_1 = F'(-1) = sum j * F_j * (-1)^{j-1}

print("Taylor expansion at x = -1: F(T, -1+y) = sum E_k y^k")
for n in [5, 7]:
    m = n*(n-1)//2
    print(f"\n  n={n}:")

    seen = set()
    all_F = []
    if n == 5:
        iterator = range(1 << m)
    else:
        random.seed(42)
        iterator = [random.getrandbits(m) for _ in range(300)]

    for bits in iterator:
        A = tournament_from_bits(bits, n)
        F = compute_F(A, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)
        all_F.append(F)

    for k in range(min(n, 5)):
        # E_k = sum_j F_j * C(j, k) * (-1)^{j-k} / k!  ... actually
        # F(T, -1+y) = sum_j F_j (-1+y)^j = sum_j F_j sum_m C(j,m) (-1)^{j-m} y^m
        # So E_k = sum_j F_j C(j, k) (-1)^{j-k}
        ek_vals = []
        for F in all_F:
            ek = sum(F[j] * math.comb(j, k) * ((-1)**(j-k)) for j in range(n))
            ek_vals.append(int(ek))

        diffs = set()
        for i in range(len(ek_vals)):
            for j in range(i+1, len(ek_vals)):
                diffs.add(abs(ek_vals[i] - ek_vals[j]))
        diffs.discard(0)

        if not diffs:
            print(f"    E_{k} (at x=-1): UNIVERSAL = {ek_vals[0]}")
        else:
            g = reduce(gcd, diffs)
            v2 = 0
            gg = g
            while gg % 2 == 0:
                v2 += 1
                gg //= 2
            print(f"    E_{k} (at x=-1): mod {g} (2^{v2}*{gg}), "
                  f"range [{min(ek_vals)}, {max(ek_vals)}]")

# ============================================================
# NOVEL: EXPANSION AT x = 0
# ============================================================
print("\n" + "=" * 60)
print("TAYLOR AT x=0: F(T, y) = sum F_k y^k (just the F_k themselves)")
print("=" * 60)

for n in [5, 7]:
    m = n*(n-1)//2
    print(f"\n  n={n}:")

    seen = set()
    all_F = []
    if n == 5:
        iterator = range(1 << m)
    else:
        random.seed(42)
        iterator = [random.getrandbits(m) for _ in range(300)]

    for bits in iterator:
        A = tournament_from_bits(bits, n)
        F = compute_F(A, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)
        all_F.append(F)

    for k in range(n):
        fk_vals = [F[k] for F in all_F]
        diffs = set()
        for i in range(len(fk_vals)):
            for j in range(i+1, len(fk_vals)):
                diffs.add(abs(fk_vals[i] - fk_vals[j]))
        diffs.discard(0)

        if not diffs:
            print(f"    F_{k}: UNIVERSAL = {fk_vals[0]}")
        else:
            g = reduce(gcd, diffs)
            v2 = 0
            gg = g
            while gg % 2 == 0:
                v2 += 1
                gg //= 2
            print(f"    F_{k}: mod {g} (2^{v2}*{gg}), range [{min(fk_vals)}, {max(fk_vals)}]")
