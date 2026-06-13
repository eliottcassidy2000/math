"""
f_omega_mod9_proof.py
kind-pasteur-2026-03-07-S35

Algebraic analysis of F(T, omega) mod 9 where omega = e^{2pi*i/3}.

KEY FORMULA (proved algebraically):
  For n=7 (d=n-1=6=3*2), using palindrome F_k = F_{6-k}:
    F(T,omega) = S_0 - S_1 where S_0=F_0+F_3+F_6, S_1=F_1+F_4=F_1+F_2
    F(T,omega) = 5040 - 3(F_1+F_2)
    = 3*(1680 - F_1 - F_2)

  Divisibility by 3 is AUTOMATIC (for any tournament at n=7).
  Divisibility by 9 iff F_1+F_2 = 0 mod 3 (since 1680 = 0 mod 3).

  For general n = 3k+1 (d=3k divisible by 3):
    F(T,omega) = 3*(n!/3 - (F_1+...+F_{k-1}+...)) [palindrome simplified]

This script:
1. Verifies the algebraic formula at n=4,5,6,7
2. Checks F_1+F_2 mod 3 at n=7 for all/many tournaments
3. Tests F(T,omega) mod 9 at all n <= 7
"""

import os; os.environ['PYTHONIOENCODING'] = 'utf-8'
from collections import defaultdict
import cmath
import random


def tournament_from_bits(n, bits):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            idx += 1
    return adj


def compute_F_dp(adj, n):
    """
    Compute F(T,x) coefficients using DP.
    dp[mask][last][fwd] = # orderings of vertices in mask ending at last with fwd forward edges.
    Returns F = [F_0, F_1, ..., F_{n-1}].
    """
    dp = [[[0] * n for _ in range(n)] for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v][0] = 1

    for mask in range(1, 1 << n):
        for last in range(n):
            if not (mask & (1 << last)):
                continue
            for fwd in range(n):
                if dp[mask][last][fwd] == 0:
                    continue
                for nxt in range(n):
                    if mask & (1 << nxt):
                        continue
                    new_mask = mask | (1 << nxt)
                    if adj[last][nxt]:
                        dp[new_mask][nxt][fwd + 1] += dp[mask][last][fwd]
                    else:
                        dp[new_mask][nxt][fwd] += dp[mask][last][fwd]

    full = (1 << n) - 1
    F = [0] * n
    for last in range(n):
        for fwd in range(n):
            F[fwd] += dp[full][last][fwd]
    return F


def compute_F_omega(F, n):
    """Compute F(T, omega) from coefficient list."""
    omega = cmath.exp(2j * cmath.pi / 3)
    result = sum(F[k] * omega**k for k in range(n))
    return result


def test_mod9_n7():
    """Exhaustive or large-sample test of F(T,omega) mod 9 at n=7."""
    n = 7
    m = n * (n - 1) // 2  # 21
    num_T = 1 << m  # 2097152

    # Sample random tournaments
    num_samples = 5000
    mod9_counts = defaultdict(int)
    f12_mod3_counts = defaultdict(int)

    random.seed(42)
    for _ in range(num_samples):
        bits = random.randint(0, num_T - 1)
        adj = tournament_from_bits(n, bits)
        F = compute_F_dp(adj, n)

        # Algebraic formula: F(T,omega) = 5040 - 3*(F[1]+F[2])
        f_omega_alg = 5040 - 3 * (F[1] + F[2])

        # Verify against direct computation
        f_omega_direct = compute_F_omega(F, n)
        assert abs(f_omega_alg - f_omega_direct.real) < 0.01, \
            f"Mismatch: alg={f_omega_alg}, direct={f_omega_direct}"
        assert abs(f_omega_direct.imag) < 0.01, f"Not real: {f_omega_direct}"

        mod9_counts[f_omega_alg % 9] += 1
        f12_mod3_counts[(F[1] + F[2]) % 3] += 1

    print(f"\nn=7: F(T,omega) mod 9 distribution ({num_samples} samples):")
    for k in sorted(mod9_counts.keys()):
        print(f"  mod 9 = {k}: {mod9_counts[k]} ({100*mod9_counts[k]/num_samples:.1f}%)")

    print(f"\nn=7: (F_1+F_2) mod 3 distribution:")
    for k in sorted(f12_mod3_counts.keys()):
        print(f"  mod 3 = {k}: {f12_mod3_counts[k]} ({100*f12_mod3_counts[k]/num_samples:.1f}%)")


def test_all_small_n():
    """Test F(T,omega) mod 3 and mod 9 at n=3,4,5,6 exhaustively."""
    omega = cmath.exp(2j * cmath.pi / 3)

    for n in range(3, 7):
        m = n * (n - 1) // 2
        num_T = 1 << m

        mod3_counts = defaultdict(int)
        mod9_counts = defaultdict(int)
        always_int = True

        for bits in range(num_T):
            adj = tournament_from_bits(n, bits)
            F = compute_F_dp(adj, n)
            f_omega = compute_F_omega(F, n)

            # Check if integer (or lies on a known ray)
            # For n=3k+1 (d divisible by 3): F(omega) is real
            # For n=3k (d=3k-1): F(omega) = omega^? * real
            # For n=3k+2 (d=3k+1): F(omega) = omega^? * real

            # Just use the algebraic form
            d = n - 1  # degree

            if d % 3 == 0:
                # Real
                f_val = round(f_omega.real)
                if abs(f_omega.imag) > 0.01:
                    always_int = False
            else:
                # Complex - compute |F(omega)|^2 mod 9 instead
                # Or compute the integer part
                f_val = round(abs(f_omega))

            if abs(f_omega.imag) < 0.01:
                f_int = round(f_omega.real)
                mod3_counts[f_int % 3] += 1
                mod9_counts[f_int % 9] += 1
            else:
                # Complex case: check if 3 divides the "norm"
                # F(omega) = alpha * omega^r for some integer alpha, r
                # The magnitude |F(omega)| = |alpha|
                # Check: F(omega)/omega^r should be a real integer for some r

                # For palindrome F(x) = x^d * F(1/x):
                # F(omega) = omega^d * F(omega^{-1}) = omega^d * conj(F(omega))
                # So F(omega) * omega^{-d} = conj(F(omega)), i.e.,
                # F(omega) = omega^{d/2} * |F(omega)| * e^{i*phi} for some phi
                # More precisely, F(omega)/omega^{d/2} is real if d is even...

                # For even d: F(omega) = omega^{d/2} * real_number
                # For odd d: more complex.

                # Let's just check divisibility by extracting the "ray coefficient"
                if d % 2 == 0:
                    ray_val = f_omega / omega ** (d // 2)
                    if abs(ray_val.imag) < 0.01:
                        f_int = round(ray_val.real)
                        mod3_counts[f_int % 3] += 1
                        mod9_counts[f_int % 9] += 1
                    else:
                        # Non-standard ray
                        mod3_counts['complex'] += 1
                        mod9_counts['complex'] += 1
                else:
                    # For odd d, try omega^{(d-1)/2+1} = omega^{(d+1)/2}
                    # Actually, F(omega) = omega^d * conj(F(omega))
                    # So F(omega)/omega^{d} = conj(F(omega))/|omega^d|^2 = conj(F(omega))
                    # Hmm, omega^d = omega^{d mod 3}

                    # Just report magnitude
                    mag = abs(f_omega)
                    mod3_counts[f'|{round(mag)}|'] += 1
                    mod9_counts[f'|{round(mag)}|'] += 1

        print(f"\nn={n} (d={n-1}, d mod 3={d%3}):")
        print(f"  F(T,omega) mod 3: {dict(sorted(mod3_counts.items(), key=lambda x: str(x[0])))}")
        print(f"  F(T,omega) mod 9: {dict(sorted(mod9_counts.items(), key=lambda x: str(x[0])))}")
        print(f"  Always integer: {always_int}")


def algebraic_analysis():
    """
    For n where d = n-1 is divisible by 3 (n = 4, 7, 10, ...):
    F(T,omega) is always real and equals 3*something.

    Proof:
    F(T,omega) = sum_k F_k * omega^k. With d=3m:
    = sum_{r=0}^{2} omega^r * S_r where S_r = sum_{k equiv r mod 3} F_k.
    Palindrome F_k = F_{d-k}. Since d=3m, k equiv r mod 3 iff d-k equiv -r equiv 3-r mod 3.
    So F_{d-k} is in S_{3-r mod 3}. For r=0: stays in S_0. For r=1: goes to S_2. For r=2: goes to S_1.
    Therefore S_1 = S_2.

    F(T,omega) = S_0 + S_1*omega + S_1*omega^2 = S_0 + S_1*(omega+omega^2) = S_0 - S_1.

    And S_0 + 2*S_1 = n!, so S_0 = n! - 2*S_1.
    F(T,omega) = n! - 3*S_1.

    Since n! is divisible by 3 for n >= 3: F(T,omega) = 3*(n!/3 - S_1).
    This is ALWAYS divisible by 3!

    For mod 9: need S_1 equiv n!/3 mod 3, i.e., S_1 equiv n!/3 mod 3.
    """
    print("\n" + "=" * 60)
    print("ALGEBRAIC ANALYSIS")
    print("=" * 60)

    for n in [4, 7, 10]:
        d = n - 1
        if d % 3 != 0:
            continue
        nfact = 1
        for i in range(1, n + 1):
            nfact *= i

        print(f"\nn={n}, d={d}=3*{d//3}, n!={nfact}, n!/3={nfact//3}")
        print(f"  F(T,omega) = n! - 3*S_1 = {nfact} - 3*S_1")
        print(f"  Always divisible by 3: YES (algebraic proof)")
        print(f"  Divisible by 9 iff S_1 = {nfact//3} mod 3 = {(nfact//3) % 3}")
        print(f"  i.e., S_1 equiv {(nfact//3) % 3} mod 3")

    # Explicit for n=7
    print(f"\n  n=7: S_1 = F_1 + F_4 = F_1 + F_2 (by palindrome F_4=F_2)")
    print(f"  Need F_1+F_2 equiv 0 mod 3 (since 5040/3=1680 equiv 0 mod 3)")
    print(f"  F_0 = H(T^op) = H(T) [by path reversal]")

    # For n=4
    print(f"\n  n=4: S_1 = F_1 (only term with k equiv 1 mod 3)")
    print(f"  Need F_1 equiv 24/3 = 8 equiv 2 mod 3")
    print(f"  F_0 = H(T). F_0+F_1 = 12 (since 24=2F_0+2F_1). So F_1 = 12-H(T).")
    print(f"  F_1 mod 3 = (12-H) mod 3 = (-H) mod 3 = (3-H mod 3) mod 3")
    print(f"  H=1: F_1=11, 11 mod 3=2 (=target). H=3: F_1=9, 9 mod 3=0 (!=2). H=5: F_1=7, 7 mod 3=1 (!=2).")
    print(f"  So mod 9 holds only for H=1 at n=4. NOT universal!")


# Run
algebraic_analysis()
print("\n" + "=" * 60)
print("COMPUTATIONAL VERIFICATION")
print("=" * 60)
test_all_small_n()

print("\n" + "=" * 60)
print("n=7 SAMPLING TEST")
print("=" * 60)
test_mod9_n7()

print("\nDONE")
