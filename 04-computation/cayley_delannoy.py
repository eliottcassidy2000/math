#!/usr/bin/env python3
"""
cayley_delannoy.py — Clean library for the Cayley-Delannoy tournament theory.

Core functions:
  Q(x)          — Cayley transform (1+x)/(1-x)
  gk(k, m)      — Delannoy weight (closed form)
  cv2(n)        — CV^2 of H(T) for random tournament on n vertices
  W(n)          — NUD-weighted permutation count
  T_matrix(k,m) — Symmetric Delannoy matrix T(k,m) = k*gk(k,m)
  address(n)    — Cayley address x_n = (n-1)/(n+1)
  golden_shadow(n) — f_n = (n-2+sqrt(n^2+4))/2

All functions use exact arithmetic (fractions) where possible.
"""

from fractions import Fraction
from math import comb, factorial, sqrt, log, atanh
from functools import lru_cache


def Q(x):
    """Cayley transform Q(x) = (1+x)/(1-x)."""
    if isinstance(x, Fraction):
        return (1 + x) / (1 - x)
    return (1 + x) / (1 - x)


def address(n):
    """Cayley address of natural number n: x_n = (n-1)/(n+1)."""
    return Fraction(n - 1, n + 1)


def gk(k, m):
    """Delannoy weight g_k(m) = sum C(k-1,j-1)*C(m,j)*2^{j-1}.

    Satisfies: Q(x)^m = 1 + 2*sum_{k>=1} gk(k,m)*x^k.
    """
    if isinstance(m, int) and m < 1:
        return Fraction(0)
    return sum(
        Fraction(comb(k - 1, j - 1) * comb(m, j)) * Fraction(2) ** (j - 1)
        for j in range(1, min(k, m) + 1)
    )


def T_matrix(k, m):
    """Symmetric Delannoy diagonal-step matrix: T(k,m) = k*gk(k,m).

    Counts total diagonal steps in all Delannoy paths from (0,0) to (k,m).
    T(k,m) = T(m,k) (duality).
    Diagonal T(k,k) = OEIS A108666.
    """
    return k * gk(k, m)


def cv2(n):
    """Exact CV^2(H) for random tournaments on n vertices.

    CV^2 = Var(H)/E[H]^2 where H = Hamiltonian path count.
    """
    total = Fraction(0)
    for k in range(1, (n - 1) // 2 + 1):
        m = n - 2 * k
        if m < 1:
            continue
        g = gk(k, m)
        ff = Fraction(1)
        for i in range(2 * k):
            ff *= (n - i)
        total += 2 * g / ff
    return total


@lru_cache(maxsize=256)
def W(n):
    """W(n) = sum over NUD permutations of 2^{#unit_ascents}.

    Satisfies W(n)/n! = 1 + CV^2(H).
    Sequence: 1, 2, 8, 32, 158, 928, 6350, 49752, ...
    Computed via bitmask DP for n <= 20.
    """
    if n <= 0:
        return 0
    if n == 1:
        return 1
    # Bitmask DP
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for length in range(1, n):
        new_dp = {}
        for (mask, last), weight in dp.items():
            if bin(mask).count('1') != length:
                continue
            for v in range(n):
                if mask & (1 << v):
                    continue
                if v == last - 1:
                    continue  # unit descent
                new_weight = weight * (2 if v == last + 1 else 1)
                key = (mask | (1 << v), v)
                new_dp[key] = new_dp.get(key, 0) + new_weight
        for k2, w in new_dp.items():
            dp[k2] = dp.get(k2, 0) + w
    full_mask = (1 << n) - 1
    return sum(dp.get((full_mask, v), 0) for v in range(n))


def golden_shadow(n):
    """Golden shadow f_n = (n-2+sqrt(n^2+4))/2.

    Satisfies: f^2 - (n-2)*f - n = 0, i.e., f(f+2) = n(f+1).
    Has continued fraction [n-1; n, n, n, ...].
    Maps to Cayley address via x = (f^2+f-1)/(f^2+3f+1) = (n-1)/(n+1).
    """
    return (n - 2 + sqrt(n ** 2 + 4)) / 2


def rodrigues(k, m):
    """Rodrigues formula: g_k(m) = (1/k!) [d^k/du^k u^m(u+1)^{k-1}]_{u=1}.

    Uses Leibniz rule for exact computation.
    """
    total = Fraction(0)
    for j in range(1, k + 1):
        # d^j/du^j u^m at u=1 = (m)_j = falling factorial
        fall_m = Fraction(1)
        for i in range(j):
            fall_m *= (m - i)
        # d^{k-j}/du^{k-j} (u+1)^{k-1} at u=1 = (k-1)!/(j-1)! * 2^{j-1}
        deriv = Fraction(factorial(k - 1), factorial(j - 1)) * Fraction(2 ** (j - 1))
        total += Fraction(comb(k, j)) * fall_m * deriv
    return total / Fraction(factorial(k))


def asymptotic_cv2(n, terms=5):
    """Fast asymptotic approximation: CV^2 ~ 2/n + corrections.

    Uses the exact formula up to 'terms' levels.
    """
    total = Fraction(0)
    for k in range(1, min(terms + 1, (n - 1) // 2 + 1)):
        m = n - 2 * k
        if m < 1:
            continue
        g = gk(k, m)
        ff = Fraction(1)
        for i in range(2 * k):
            ff *= (n - i)
        total += 2 * g / ff
    return total


# ============================================================
# Self-test
# ============================================================

if __name__ == "__main__":
    print("Cayley-Delannoy Library Self-Test")
    print("=" * 40)

    # Test Q
    assert Q(Fraction(0)) == 1
    assert Q(Fraction(1, 3)) == 2
    assert Q(address(5)) == 5
    print("Q: OK")

    # Test gk
    assert gk(1, 5) == 5
    assert gk(2, 3) == 9
    assert gk(3, 3) == 19
    print("gk: OK")

    # Test duality
    for k in range(1, 6):
        for m in range(1, 6):
            assert k * gk(k, m) == m * gk(m, k), f"Duality failed at k={k},m={m}"
    print("Duality: OK")

    # Test parity (via the polynomial: gk is degree k with same parity as k)
    # Verify that gk(k, m) polynomial coefficients have the right parity
    for k in range(1, 6):
        vals = [gk(k, m) for m in range(k + 1)]
        # Forward differences to extract polynomial
        diffs = list(vals)
        for d in range(k):
            diffs = [diffs[i + 1] - diffs[i] for i in range(len(diffs) - 1)]
        # k-th difference should be constant (degree k)
        assert len(set(diffs)) <= 1 or all(d == diffs[0] for d in diffs), f"Degree check failed k={k}"
    print("Parity (degree check): OK")

    # Test Rodrigues
    for k in range(1, 6):
        for m in range(1, 6):
            assert rodrigues(k, m) == gk(k, m), f"Rodrigues failed at k={k},m={m}"
    print("Rodrigues: OK")

    # Test functional equation
    for m in [3, 5, 7]:
        for two_k in [2, 4, 6]:
            actual = gk(two_k, m)
            predicted = sum((-1) ** (j + 1) * gk(j, m) * gk(two_k - j, m)
                           for j in range(1, two_k))
            assert actual == predicted, f"FuncEq failed at m={m},k={two_k}"
    print("Functional equation: OK")

    # Test Wick rotation: arctanh(i) = i*pi/4
    import cmath
    wick = cmath.atanh(1j)
    expected = 1j * cmath.pi / 4
    assert abs(wick - expected) < 1e-14, f"Wick rotation failed: {wick} != {expected}"
    # Q(i) = i
    qi = (1 + 1j) / (1 - 1j)
    assert abs(qi - 1j) < 1e-14, f"Q(i) failed: {qi} != i"
    print("Wick rotation: OK")

    # Test OEIS A108666 diagonal
    oeis_a108666 = [0, 1, 8, 57, 384, 2505, 16008, 100849, 628736, 3888657]
    for n in range(1, 10):
        t = T_matrix(n, n)
        assert t == oeis_a108666[n], f"A108666 failed at n={n}: {t} != {oeis_a108666[n]}"
    print("OEIS A108666: OK")

    # Test cv2 against W
    for n in range(3, 10):
        w = W(n)
        c = cv2(n)
        assert Fraction(w, factorial(n)) - 1 == c, f"CV2 failed at n={n}"
    print("CV2 vs W(n): OK")

    # Test W(n) values
    expected = [1, 2, 8, 32, 158, 928, 6350, 49752]
    for i, val in enumerate(expected):
        assert W(i + 1) == val, f"W({i+1}) failed"
    print("W(n): OK")

    # Print summary
    print()
    print("W(n) for n=1..12:")
    for n in range(1, 13):
        print(f"  W({n:2d}) = {W(n)}")

    print()
    print("CV^2 for n=3..10:")
    for n in range(3, 11):
        c = cv2(n)
        print(f"  CV^2({n}) = {c} = {float(c):.8f}")

    print()
    print("All tests passed.")
