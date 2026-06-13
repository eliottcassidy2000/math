"""
thm086_verify.py
kind-pasteur-2026-03-07-S37

Verify THM-086: c_j(T) = 0 mod 3 for j < 2*floor((n-1)/2) for ALL tournaments.

Also verify the COROLLARY: 3|A(n,k) => 3|F_k(T) for all T
(the Eulerian conjecture).

Extended sampling at n=9,10 to increase confidence.
"""

import os; os.environ['PYTHONIOENCODING'] = 'utf-8'
import random
from math import comb


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


def eulerian_numbers(n):
    if n == 1:
        return [1]
    prev = eulerian_numbers(n - 1)
    A = [0] * n
    for k in range(n):
        A[k] = (k + 1) * prev[k] if k < len(prev) else 0
        if k > 0:
            A[k] += (n - k) * prev[k - 1]
    return A


def test_thm086(n, num_samples, seed=42):
    """Test THM-086 and Eulerian conjecture at given n."""
    m = n * (n - 1) // 2
    num_T = 1 << m
    val = 2 * ((n - 1) // 2)  # predicted universal zero count

    A = eulerian_numbers(n)
    eulerian_zeros = set(k for k in range(n) if A[k] % 3 == 0)

    random.seed(seed)

    # Track failures
    cj_failures = [0] * n  # c_j != 0 mod 3
    fk_failures = [0] * n  # F_k != 0 mod 3 when A(n,k)=0 mod 3

    for trial in range(num_samples):
        bits = random.randint(0, num_T - 1)
        adj = tournament_from_bits(n, bits)
        F = compute_F_dp(adj, n)

        # Compute c_j
        c = [0] * n
        for j in range(n):
            for k in range(j, n):
                c[j] += comb(k, j) * F[k]

        for j in range(n):
            if c[j] % 3 != 0:
                cj_failures[j] += 1

        for k in eulerian_zeros:
            if F[k] % 3 != 0:
                fk_failures[k] += 1

    # Report
    print(f"\n  n={n} ({num_samples} samples), val={val}:")

    # THM-086 check
    thm086_ok = True
    for j in range(val):
        if cj_failures[j] > 0:
            thm086_ok = False
            print(f"    ** THM-086 FAILS: c_{j} has {cj_failures[j]} failures **")

    if thm086_ok:
        print(f"    THM-086: c_j=0 mod 3 for j<{val}: PASS (0 failures)")

    # Check that c_{val} is NOT always 0
    if cj_failures[val] > 0:
        pct = 100 * (num_samples - cj_failures[val]) / num_samples
        print(f"    c_{val}: {pct:.1f}% divisible by 3 (NOT universal, as expected)")
    else:
        print(f"    c_{val}: 100% divisible by 3 (unexpected? may need more samples)")

    # Eulerian conjecture check
    eul_ok = True
    for k in sorted(eulerian_zeros):
        if fk_failures[k] > 0:
            eul_ok = False
            print(f"    ** EULERIAN FAILS: F_{k} has {fk_failures[k]} failures **")

    if eul_ok:
        if eulerian_zeros:
            print(f"    Eulerian conjecture: 3|A(n,k) => 3|F_k(T) at k={sorted(eulerian_zeros)}: PASS")
        else:
            print(f"    Eulerian conjecture: no zeros to check (all A(n,k) nonzero mod 3)")

    return thm086_ok and eul_ok


# Run tests
print("=" * 70)
print("THM-086 and Eulerian Conjecture Verification")
print("=" * 70)

all_pass = True
for n, ns in [(5, 'exhaustive'), (6, 'exhaustive'), (7, 10000), (8, 5000), (9, 3000), (10, 500)]:
    if ns == 'exhaustive':
        m = n * (n - 1) // 2
        ns_int = 1 << m
        # Run exhaustive
        m_bits = n * (n - 1) // 2
        num_T = 1 << m_bits
        val = 2 * ((n - 1) // 2)
        A = eulerian_numbers(n)
        eulerian_zeros = set(k for k in range(n) if A[k] % 3 == 0)

        cj_failures = [0] * n
        fk_failures = [0] * n
        for bits in range(num_T):
            adj = tournament_from_bits(n, bits)
            F = compute_F_dp(adj, n)
            c = [0] * n
            for j in range(n):
                for k in range(j, n):
                    c[j] += comb(k, j) * F[k]
            for j in range(n):
                if c[j] % 3 != 0:
                    cj_failures[j] += 1
            for k in eulerian_zeros:
                if F[k] % 3 != 0:
                    fk_failures[k] += 1

        print(f"\n  n={n} (exhaustive, {num_T} tournaments), val={val}:")
        ok = True
        for j in range(val):
            if cj_failures[j] > 0:
                ok = False
                print(f"    ** THM-086 FAILS: c_{j} has {cj_failures[j]} failures **")
        if ok:
            print(f"    THM-086: c_j=0 mod 3 for j<{val}: PASS")
        if cj_failures[val] > 0:
            pct = 100 * (num_T - cj_failures[val]) / num_T
            print(f"    c_{val}: {pct:.1f}% div by 3 (NOT universal)")
        for k in sorted(eulerian_zeros):
            if fk_failures[k] > 0:
                ok = False
                print(f"    ** EULERIAN FAILS at k={k} **")
        if ok and eulerian_zeros:
            print(f"    Eulerian conjecture: PASS at k={sorted(eulerian_zeros)}")
        all_pass = all_pass and ok
    else:
        ok = test_thm086(n, ns)
        all_pass = all_pass and ok

print("\n" + "=" * 70)
if all_pass:
    print("ALL TESTS PASSED")
else:
    print("SOME TESTS FAILED")
print("=" * 70)

print("\nDONE")
