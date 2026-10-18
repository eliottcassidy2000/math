"""S40 audit: THM-4666 statement 6 (spectral ceiling).

 (i) eigenvalues of the n x n sign matrix T_n (+1 above, -1 below the diagonal) against +-i cot((2k-1) pi/(2n));
 (ii) ||S||_op = sqrt(p) for the skew adjacency matrix of P_p;
 (iii) the trace-norm bound B(p) = (sqrt p / 2) sum_k |cot((2k-1) pi/(2p))| against (p^1.5/pi)(log p + gamma + log(4/pi));
 (iv) ratio exact optimum / B(p) at p = 7, 11, 19, 23, 31, 43 (exact optima 7, 15, 43, 69, 105, 183);
 (v) comparison with Alon's explicit Corollary 3.2 bound p^1.5 ceil(log2 p) and the best known witnesses.
Usage: python audit_spectral_bound.py
"""
import math
import numpy as np
import mpmath as mp


def T(n):
    return np.triu(np.ones((n, n)), 1) - np.tril(np.ones((n, n)), -1)


def main():
    for n in (2, 3, 5, 8, 13, 31, 64):
        ev = np.sort(np.linalg.eigvals(T(n)).imag)
        pred = np.sort(np.array([1 / math.tan((2 * k - 1) * math.pi / (2 * n)) if abs(math.cos((2 * k - 1) * math.pi / (2 * n))) > 1e-15 else 0.0
                                 for k in range(1, n + 1)]))
        print('n=%d max |eig - cot| = %.2e' % (n, np.max(np.abs(ev - pred))))
    for p in (7, 11, 43, 103):
        leg = [0] + [1 if pow(a, (p - 1) // 2, p) == 1 else -1 for a in range(1, p)]
        S = np.array([[leg[(v - u) % p] for v in range(p)] for u in range(p)], dtype=float)
        print('p=%d ||S||_op = %.10f, sqrt(p) = %.10f' % (p, np.linalg.norm(S, 2), math.sqrt(p)))
    g = float(mp.euler)
    for n in (10 ** 3, 10 ** 5):
        k = np.arange(1, n + 1)
        s = np.sum(np.abs(1 / np.tan((2 * k - 1) * np.pi / (2 * n))))
        a = (2 * n / math.pi) * (math.log(n) + g + math.log(4 / math.pi))
        print('n=%d sum|cot| = %.10f, asymptotic = %.10f, rel diff %.2e' % (n, s, a, s / a - 1))
    exact = {7: 7, 11: 15, 19: 43, 23: 69, 31: 105, 43: 183}
    for p, e in exact.items():
        k = np.arange(1, p + 1)
        B = math.sqrt(p) / 2 * np.sum(np.abs(1 / np.tan((2 * k - 1) * np.pi / (2 * p))))
        alon = p ** 1.5 * math.ceil(math.log2(p))
        print('p=%d exact optimum %d, ceiling B(p) = %.1f, ratio %.3f; Alon Cor 3.2 bound %d (ratio %.3f)' % (p, e, B, e / B, alon, e / alon))
    # ratio at the large witnesses found in the audit
    for p, e in {263: 4089, 383: 7519, 463: 9553, 479: 11975}.items():
        k = np.arange(1, p + 1)
        B = math.sqrt(p) / 2 * np.sum(np.abs(1 / np.tan((2 * k - 1) * np.pi / (2 * p))))
        print('p=%d best known excess %d, ceiling %.1f, ratio %.3f' % (p, e, B, e / B))


if __name__ == '__main__':
    main()
