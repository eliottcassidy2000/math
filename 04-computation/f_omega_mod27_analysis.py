"""
f_omega_mod27_analysis.py
kind-pasteur-2026-03-07-S36

Deep analysis of F(T, omega) modular arithmetic.

KNOWN (S35):
  - F(T,omega) divisible by 3 for ALL n >= 3 (algebraic proof)
  - F(T,omega) divisible by 9 for ALL n >= 6 (verified: n=6 exhaustive, n=7 sampled)
  - NOT divisible by 9 at n=4 (H=3,5 give non-zero mod 9) or n=5

THIS SCRIPT:
  1. Exhaustive mod 27 check at n=3,4,5,6
  2. Sampled mod 27 check at n=7
  3. Investigate WHY S_1 = 0 mod 3 at n >= 6
  4. Check mod 81 (3^4) to find the exact 3-adic valuation pattern

KEY FORMULAS (from S35 algebraic analysis):
  When d = n-1 = 0 mod 3 (n = 4, 7, 10, ...):
    F(T,omega) = n! - 3*S_1, where S_1 = sum_{k = 1 mod 3} F_k
    v_3(F(T,omega)) >= 1 always (since 3 | n! for n >= 3)
    v_3 >= 2 iff S_1 = n!/3 mod 3

  When d = 1 mod 3 (n = 3, 6, 9, ...):
    Palindrome gives S_0 = S_2, S_1 = S_1
    F(T,omega) = omega * (S_1 - S_0) = omega * (3*S_1/2 - n!/2)
    Wait, let me rederive carefully.

  When d = 2 mod 3 (n = 5, 8, 11, ...):
    Different palindrome pairings.
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
    """Compute F(T,x) coefficients using DP."""
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


def compute_F_omega_exact(F, n):
    """Compute F(T, omega) exactly using integer arithmetic.

    omega = e^{2pi*i/3} satisfies omega^2 + omega + 1 = 0.
    F(T, omega) = sum_k F_k * omega^k.

    Write result as a + b*omega where a, b are integers.
    Then F(T,omega) = a + b*omega = (a - b/2) + b*sqrt(3)/2 * i

    Actually: 1, omega, omega^2 with omega^2 = -1-omega.
    So F(T,omega) = S_0 + S_1*omega + S_2*omega^2
                  = S_0 + S_1*omega + S_2*(-1-omega)
                  = (S_0 - S_2) + (S_1 - S_2)*omega

    where S_r = sum_{k = r mod 3} F_k.
    """
    S = [0, 0, 0]
    for k in range(n):
        S[k % 3] += F[k]

    # F(T,omega) = (S_0 - S_2) + (S_1 - S_2)*omega
    a = S[0] - S[2]
    b = S[1] - S[2]
    return a, b, S


def analyze_mod_powers(n, exhaustive=True, num_samples=10000):
    """Analyze F(T,omega) mod 3^k for various k."""
    m = n * (n - 1) // 2
    num_T = 1 << m

    if not exhaustive:
        random.seed(42)
        indices = [random.randint(0, num_T - 1) for _ in range(num_samples)]
    else:
        indices = range(num_T)
        num_samples = num_T

    # Track mod 3, 9, 27, 81
    # F(T,omega) = a + b*omega. For this to be divisible by 3^k,
    # we need a = 0 mod 3^k AND b = 0 mod 3^k.
    results = {3: 0, 9: 0, 27: 0, 81: 0, 243: 0}
    total = 0

    # Also track the distribution of (a mod 27, b mod 27)
    ab_mod27 = defaultdict(int)
    # Track S_r mod 3 distribution
    S_mod3 = defaultdict(int)

    for bits in indices:
        adj = tournament_from_bits(n, bits)
        F = compute_F_dp(adj, n)
        a, b, S = compute_F_omega_exact(F, n)

        total += 1

        for p in [3, 9, 27, 81, 243]:
            if a % p == 0 and b % p == 0:
                results[p] += 1

        ab_mod27[(a % 27, b % 27)] += 1
        S_mod3[(S[0] % 3, S[1] % 3, S[2] % 3)] += 1

    print(f"\n{'='*60}")
    print(f"n={n} ({'exhaustive' if exhaustive else f'{num_samples} samples'})")
    print(f"{'='*60}")
    print(f"d = n-1 = {n-1}, d mod 3 = {(n-1) % 3}")
    print(f"n! = {factorial(n)}, v_3(n!) = {v3(factorial(n))}")

    for p in [3, 9, 27, 81, 243]:
        pct = 100 * results[p] / total
        print(f"  3^{p.bit_length()-1:d} = {p:4d} divides F(T,omega): {results[p]}/{total} ({pct:.1f}%)")
        if results[p] == total:
            print(f"    ^^ UNIVERSAL!")

    print(f"\n  (S_0 mod 3, S_1 mod 3, S_2 mod 3) distribution:")
    for key in sorted(S_mod3.keys()):
        pct = 100 * S_mod3[key] / total
        print(f"    {key}: {S_mod3[key]} ({pct:.1f}%)")

    # For palindrome analysis
    d = n - 1
    if d % 3 == 0:
        # S_1 = S_2 (palindrome). F(omega) = S_0 - S_1 = n! - 3*S_1
        print(f"\n  Palindrome: d = 0 mod 3 => S_1 = S_2")
        print(f"  F(T,omega) = n! - 3*S_1 = {factorial(n)} - 3*S_1")
        print(f"  v_3 >= 2 iff S_1 = {factorial(n)//3 % 3} mod 3")
    elif d % 3 == 1:
        print(f"\n  Palindrome: d = 1 mod 3 => S_0 = S_1")
        print(f"  F(T,omega) = omega^2 * (n! - 3*S_0)")
    elif d % 3 == 2:
        print(f"\n  Palindrome: d = 2 mod 3 => S_0 = S_2")
        print(f"  F(T,omega) = omega * (n! - 3*S_0)")

    return results, total


def factorial(n):
    r = 1
    for i in range(1, n+1):
        r *= i
    return r


def v3(n):
    """3-adic valuation of n."""
    if n == 0:
        return float('inf')
    v = 0
    while n % 3 == 0:
        v += 1
        n //= 3
    return v


def study_S1_structure(n, exhaustive=True, num_samples=5000):
    """Deep dive into why S_1 = 0 mod 3 at n >= 6."""
    m = n * (n - 1) // 2
    num_T = 1 << m
    d = n - 1

    if not exhaustive:
        random.seed(42)
        indices = [random.randint(0, num_T - 1) for _ in range(num_samples)]
    else:
        indices = range(num_T)
        num_samples = num_T

    print(f"\n{'='*60}")
    print(f"S_1 structure at n={n}")
    print(f"{'='*60}")

    # The key quantity depends on d mod 3:
    # d = 0: F(omega) = n! - 3*S_1. Mod 9 iff S_1 = n!/3 mod 3.
    # d = 1: F(omega) = omega^2*(n! - 3*S_0). Mod 9 iff S_0 = n!/3 mod 3.
    # d = 2: F(omega) = omega*(n! - 3*S_0). Mod 9 iff S_0 = n!/3 mod 3.

    if d % 3 == 0:
        target_mod = (factorial(n) // 3) % 3
        print(f"  Need S_1 = {target_mod} mod 3 (since n!/3 = {factorial(n)//3} = {target_mod} mod 3)")
    elif d % 3 == 1:
        target_mod = (factorial(n) // 3) % 3
        print(f"  Need S_0 = {target_mod} mod 3 (since n!/3 = {factorial(n)//3} = {target_mod} mod 3)")
    else:
        target_mod = (factorial(n) // 3) % 3
        print(f"  Need S_0 = {target_mod} mod 3 (since n!/3 = {factorial(n)//3} = {target_mod} mod 3)")

    # Track the key quantity mod 3
    key_mod3 = defaultdict(int)

    # Track F_k mod 3 individually
    Fk_mod3_dist = [defaultdict(int) for _ in range(n)]

    # Track relationship with H(T)
    H_vs_key = defaultdict(list)

    for bits in indices:
        adj = tournament_from_bits(n, bits)
        F = compute_F_dp(adj, n)
        a, b, S = compute_F_omega_exact(F, n)
        H = F[n-1]

        if d % 3 == 0:
            key = S[1]
        else:
            key = S[0]

        key_mod3[key % 3] += 1
        H_vs_key[H].append(key % 3)

        for k in range(n):
            Fk_mod3_dist[k][F[k] % 3] += 1

    print(f"\n  Key quantity mod 3:")
    for r in sorted(key_mod3.keys()):
        pct = 100 * key_mod3[r] / num_samples
        print(f"    = {r} mod 3: {key_mod3[r]} ({pct:.1f}%)")

    print(f"\n  Individual F_k mod 3:")
    for k in range(n):
        dist = Fk_mod3_dist[k]
        parts = [f"={r}: {dist[r]}" for r in range(3)]
        print(f"    F_{k}: {', '.join(parts)}")

    # Check: is the congruence forced by some simpler invariant?
    print(f"\n  Key quantity mod 3 by H(T) value:")
    for H in sorted(H_vs_key.keys())[:15]:
        vals = H_vs_key[H]
        counts = [vals.count(r) for r in range(3)]
        if any(c > 0 for c in counts):
            print(f"    H={H}: mod 3 = {counts}")


def study_individual_Fk_mod3(n, exhaustive=True, num_samples=5000):
    """Check if individual F_k are divisible by 3."""
    m = n * (n - 1) // 2
    num_T = 1 << m

    if not exhaustive:
        random.seed(42)
        indices = [random.randint(0, num_T - 1) for _ in range(num_samples)]
    else:
        indices = range(num_T)
        num_samples = num_T

    print(f"\n{'='*60}")
    print(f"Individual F_k mod 3 at n={n}")
    print(f"{'='*60}")

    # For each k, track how often F_k = 0 mod 3
    always_div3 = [True] * n
    counts = [[0, 0, 0] for _ in range(n)]

    for bits in indices:
        adj = tournament_from_bits(n, bits)
        F = compute_F_dp(adj, n)

        for k in range(n):
            r = F[k] % 3
            counts[k][r] += 1
            if r != 0:
                always_div3[k] = False

    for k in range(n):
        div3_pct = 100 * counts[k][0] / num_samples
        status = "ALWAYS = 0" if always_div3[k] else ""
        print(f"  F_{k}: =0: {counts[k][0]} ({div3_pct:.1f}%), =1: {counts[k][1]}, =2: {counts[k][2]}  {status}")


def study_F1_mod3_algebraic(n):
    """
    At n=7 (d=6 = 0 mod 3): need S_1 = F_1 + F_4 = 2*F_1 (palindrome F_4=F_2, wait...)

    Actually at n=7: F_k for k=0,...,6. Palindrome: F_k = F_{6-k}.
    S_1 = sum_{k = 1 mod 3} F_k = F_1 + F_4 = F_1 + F_2 (since F_4 = F_{6-4} = F_2).

    Hmm, F_4 = F_{6-4} = F_2. So S_1 = F_1 + F_2.
    But F_1 and F_2 are different in general.

    At n=6 (d=5 = 2 mod 3): need S_0 = F_0 + F_3 = F_0 + F_2 (palindrome F_3=F_2).
    S_0 = H(T) + F_2. Need this = n!/3 mod 3 = 240 mod 3 = 0 mod 3.
    So H(T) + F_2 = 0 mod 3.

    Can we prove H(T) + F_2 = 0 mod 3 at n=6?
    """
    m = n * (n - 1) // 2
    num_T = 1 << m
    d = n - 1

    print(f"\n{'='*60}")
    print(f"Algebraic analysis of key congruence at n={n}")
    print(f"{'='*60}")

    # Compute S_r decomposition explicitly
    for bits in range(min(16, num_T)):
        adj = tournament_from_bits(n, bits)
        F = compute_F_dp(adj, n)
        H = F[n-1]
        a, b, S = compute_F_omega_exact(F, n)

        print(f"  bits={bits}: F={F}, H={H}")
        print(f"    S = {S}, S mod 3 = {[s%3 for s in S]}")
        print(f"    a={a}, b={b}, a mod 3={a%3}, b mod 3={b%3}")
        f_omega = a + b * cmath.exp(2j * cmath.pi / 3)
        v = v3(round(abs(a))) if a != 0 else 'inf'
        print(f"    F(omega) ~ {f_omega:.1f}, v_3(a)={'inf' if a==0 else v3(abs(a))}, v_3(b)={'inf' if b==0 else v3(abs(b))}")


# Run analysis
for n in [3, 4, 5, 6]:
    analyze_mod_powers(n, exhaustive=True)

# n=7: sampled
analyze_mod_powers(7, exhaustive=False, num_samples=5000)

# Deep S_1 structure
for n in [4, 5, 6]:
    study_S1_structure(n, exhaustive=True)
study_S1_structure(7, exhaustive=False, num_samples=5000)

# Individual F_k mod 3
for n in [3, 4, 5, 6]:
    study_individual_Fk_mod3(n, exhaustive=True)
study_individual_Fk_mod3(7, exhaustive=False, num_samples=5000)

# Algebraic deep dive at n=6
study_F1_mod3_algebraic(6)

print("\nDONE")
