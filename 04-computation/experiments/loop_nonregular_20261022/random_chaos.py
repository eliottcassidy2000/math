"""The forest identity is the Walsh (chaos) expansion of kappa for uniformly random tournaments.  Session opus-2026-10-22-S54.

kappa = sum_L N_L/(n)_L, and N_L = sum over undirected linear forests F of K_n with L edges whose paths all have even
length of 2^p(F) * prod_(e in F) s_e (s_e = +-1 the arc signs; an odd path cancels against its reversal).  Distinct
forests are distinct edge sets, so under independent fair arcs the monomials are orthonormal and
  E kappa = 1,   E kappa^2 = sum_L sum_F 4^p(F) / ((n)_L)^2 = sum_lambda 2^p n! / ((n - L - p)! prod m_l! ((n)_L)^2),
the last sum over multisets lambda of even path lengths (L = |lambda|, p = number of parts, m_l multiplicities).
The L = 2 chaos is N_2 = n(n - 1) - ||S 1||^2 (the score energy), with variance share 2(n - 2)/(n(n - 1)).
Prints the exact Var(kappa), the L = 2 share, and the empirical values on the random samples.
Usage: python random_chaos.py [FILE n ...]"""
import sys, os, math
from fractions import Fraction
from functools import lru_cache
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def falling(n, L):
    return math.prod(range(n - L + 1, n + 1)) if L else 1


def second_moment(n):
    """exact E kappa^2 as a Fraction, by chaos level L."""
    by_L = {}

    def parts(maxpart, budget):          # multisets of even parts >= 2, nonincreasing, with L + p <= n
        yield []
        for l in range(min(maxpart, budget), 1, -1):
            if l % 2:
                continue
            for rest in parts(l, budget - l - 1):
                yield [l] + rest

    for lam in parts(n, n - 1):
        L = sum(lam); p = len(lam)
        if L + p > n:
            continue
        mult = math.prod(math.factorial(lam.count(l)) for l in set(lam))
        term = Fraction(2 ** p * math.factorial(n), math.factorial(n - L - p) * mult * falling(n, L) ** 2)
        by_L[L] = by_L.get(L, Fraction(0)) + term
    return by_L


def brute_second_moment(n):
    """E kappa^2 over all 2^C(n,2) tournaments, exact (Fraction), Hamiltonian paths by subset DP."""
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    tot = 0
    for mask in range(1 << len(pairs)):
        A = [[0] * n for _ in range(n)]
        for b, (i, j) in enumerate(pairs):
            if mask >> b & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        dp = [[0] * n for _ in range(1 << n)]
        for v in range(n):
            dp[1 << v][v] = 1
        for S in range(1, 1 << n):
            for v in range(n):
                if dp[S][v]:
                    for w in range(n):
                        if not S >> w & 1 and A[v][w]:
                            dp[S | 1 << w][w] += dp[S][v]
        H = sum(dp[(1 << n) - 1])
        tot += H * H
    return Fraction(tot * 4 ** (n - 1), (1 << len(pairs)) * math.factorial(n) ** 2)


if __name__ == '__main__':
    for n in (3, 4, 5, 6):
        ex = sum(second_moment(n).values())
        bf = brute_second_moment(n)
        print(f'n={n}: E kappa^2 formula {ex} = {float(ex):.6f}; brute force over all tournaments {bf}; equal: {ex == bf}')
    for n in (3, 4, 5, 6, 8, 10, 12, 14, 16, 20, 30, 50):
        by_L = second_moment(n)
        var = sum(v for L, v in by_L.items() if L > 0)
        share2 = by_L.get(2, Fraction(0)) / var
        print(f'n={n}: Var(kappa) = {float(var):.6f} (n Var = {n * float(var):.4f}); L=2 (score) share {float(share2):.4f}; '
              f'by level: ' + ', '.join(f'L{L}: {float(v):.3e}' for L, v in sorted(by_L.items()) if L > 0 and L <= 8))
    args = sys.argv[1:]
    for i in range(0, len(args), 2):
        path, n = args[i], int(args[i + 1])
        iu = np.triu_indices(n, 1)
        K, N2 = [], []
        for line in open(path):
            code, H = line.split()[0], int(line.split()[1])
            A = np.zeros((n, n)); b = np.array([1 if ch == '1' else 0 for ch in code])
            A[iu[0], iu[1]] = b; A[iu[1], iu[0]] = 1 - b
            S = A - A.T; d = S @ np.ones(n)
            K.append(H * 2 ** (n - 1) / math.factorial(n)); N2.append(n * (n - 1) - d @ d)
        K = np.array(K); N2 = np.array(N2)
        proj = 1 + N2 / (n * (n - 1))
        print(f'{os.path.basename(path)} (n={n}, {len(K)} samples): mean kappa {K.mean():.4f}, var {K.var():.5f} '
              f'(exact {float(sum(v for L, v in second_moment(n).items() if L > 0)):.5f}); corr(kappa, 1 + N2/(n)_2) {np.corrcoef(K, proj)[0, 1]:.4f}, '
              f'var of the L=2 projection {proj.var():.5f} (exact {2 * (n - 2) / (n * (n - 1)):.5f})')
