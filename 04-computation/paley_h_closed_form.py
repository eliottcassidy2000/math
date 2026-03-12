"""
paley_h_closed_form.py — Search for a closed-form expression for H(T_p)

Known values:
  H(T_3) = 3
  H(T_7) = 189
  H(T_11) = 95095

Candidate forms:
  H = n! * f(n) for some rational f(n)?
  H = (2n-1)!! / g(n)?
  H as product of linear factors?

Check: 3 = 3, 189 = 27*7, 95095 = ?

Also check: H/n (orbit count under Z_n rotation)
  H(T_3)/3 = 1
  H(T_7)/7 = 27
  H(T_11)/11 = 8645

Factor: 27 = 3^3, 8645 = 5 * 1729 = 5 * 7 * 13 * 19... hmm
Actually 8645 = 5 * 1729 and 1729 = 7 * 247 = 7 * 13 * 19. So 8645 = 5*7*13*19.

These are (p-2)(p-4)...(2)(1) products? Let me check:
For p=7: 5*3*1 = 15, 27 ≠ 15.
For p=11: 9*7*5*3*1 = 945, 8645 ≠ 945.

Double factorial? (p-2)!! = (p-2)(p-4)...3*1
For p=7: 5!! = 15, 27 ≠ 15.

189 = 3^3 * 7 = 7 * 27
95095 = 5 * 7 * 11 * 13 * 19

Hmm: 95095 = 5*7*11*13*19? Let me verify: 5*7=35, 35*11=385, 385*13=5005, 5005*19=95095. YES!

So H(T_11) = 5 * 7 * 11 * 13 * 19 = 11 * (5*7*13*19) = 11 * 8645.

What about H(T_7) = 189 = 7 * 27 = 7 * 3^3.
And 3^3... that doesn't fit an obvious pattern with 5*7*13*19.

Let me look at this differently.
  H(T_3) = 3 = 3
  H(T_7) = 189 = 7 * 27
  H(T_11) = 95095 = 11 * 8645

For n=3: H/n = 1
For n=7: H/n = 27
For n=11: H/n = 8645

Is there a pattern in 1, 27, 8645?
27 = 3^3
8645 = 5 * 7 * 13 * 19

Hmm, let me try the permanent formula.

Also: OEIS search for 1, 27, 8645...

Actually, let's compute H(T_p) directly at p=19 and p=23 if we can.

Author: opus-2026-03-12-S60
"""
import sys
import math
import numpy as np
from collections import defaultdict
from fractions import Fraction
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)


def hamiltonian_paths_dp(A, n):
    dp = defaultdict(lambda: defaultdict(int))
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    for mask in range(1, full + 1):
        if not dp[mask]:
            continue
        for v in dp[mask]:
            if dp[mask][v] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    dp[mask | (1 << w)][w] += dp[mask][v]
    return sum(dp[full][v] for v in range(n))


def paley_adj(p):
    """Adjacency matrix of Paley tournament T_p (p ≡ 3 mod 4)."""
    qr = set(pow(x, 2, p) for x in range(1, p))
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in qr:
                A[i][j] = 1
    return A


def factorize(n):
    if n <= 1:
        return []
    factors = []
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors.append(d)
            n //= d
        d += 1
    if n > 1:
        factors.append(n)
    return factors


def main():
    print("SEARCH FOR CLOSED-FORM H(T_p)")
    print("=" * 75)

    # Compute H for Paley primes
    paley_primes = [3, 7, 11]
    results = {}

    for p in paley_primes:
        A = paley_adj(p)
        H = hamiltonian_paths_dp(A, p)
        results[p] = H
        print(f"  H(T_{p}) = {H}")
        print(f"    Factorization: {factorize(H)}")
        print(f"    H/p = {H // p}, factors of H/p: {factorize(H // p)}")
        print(f"    H/p! = {Fraction(H, 1) / Fraction(math.factorial(p), 1)}")
        print(f"    (p-1)! = {math.factorial(p-1)}")
        print(f"    H/(p-1)! = {Fraction(H, math.factorial(p-1))}")

    # Try p=13 (but this is p ≡ 1 mod 4, no Paley tournament)
    # Try p=19 if feasible
    p = 19
    print(f"\n  Computing H(T_{p})...")
    A = paley_adj(p)
    import time
    t0 = time.time()
    H = hamiltonian_paths_dp(A, p)
    elapsed = time.time() - t0
    results[p] = H
    print(f"  H(T_{p}) = {H} ({elapsed:.1f}s)")
    print(f"    Factorization: {factorize(H)}")
    print(f"    H/p = {H // p}, factors of H/p: {factorize(H // p)}")

    # Look for patterns
    print(f"\n{'='*75}")
    print("PATTERN SEARCH")
    print(f"{'='*75}")

    for p in sorted(results.keys()):
        H = results[p]
        m = (p - 1) // 2  # degree of tournament regularity
        print(f"\n  p={p}: H={H}")
        # Various ratios
        print(f"    H / p = {H // p}")
        print(f"    H / p! = {H / math.factorial(p):.10f}")
        print(f"    H / (p-1)! = {H / math.factorial(p-1):.10f}")
        print(f"    H / (p-1)!! = ", end="")
        dfact = 1
        for k in range(p-1, 0, -2):
            dfact *= k
        print(f"{H / dfact:.10f}")
        print(f"    H / ((p-1)/2)!² = {H / (math.factorial(m)**2):.10f}")
        # Central binomial coefficient
        cbc = math.factorial(p-1) // (math.factorial(m) ** 2)
        print(f"    C(p-1, (p-1)/2) = {cbc}")
        print(f"    H / C(p-1, m) = {H / cbc:.10f}")
        print(f"    H / (C(p-1,m) * p) = {H / (cbc * p):.10f}")

    # Check: is 95095 = C(10,5) * 11 * something?
    # C(10,5) = 252. 95095 / 252 = 377.36... not integer
    # C(10,5) = 252, 95095 / 11 = 8645, 8645 / 252 = 34.3... no

    # Check OEIS for 3, 189, 95095
    print(f"\n  Check OEIS: 3, 189, 95095")
    print(f"  3 = 3!/2")
    print(f"  189 = 7*27 = 7*3^3")
    print(f"  95095 = 5*7*11*13*19")

    # H(T_3)/3! = 1/2
    # H(T_7)/7! = 189/5040 = 3/80
    # H(T_11)/11! = 95095/39916800 = ?
    print(f"  H/n! ratios: {[Fraction(results[p], math.factorial(p)) for p in sorted(results.keys()) if p <= 11]}")

    # Alternative: H vs the number of Hamiltonian paths in the complete graph K_n
    # K_n has n!/2 directed Hamiltonian paths
    print(f"\n  H / (n!/2):")
    for p in sorted(results.keys()):
        ratio = Fraction(results[p] * 2, math.factorial(p))
        print(f"    p={p}: {ratio} = {float(ratio):.6f}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
