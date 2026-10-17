#!/usr/bin/env python3
"""The Fourier transform of THM-1300's Keller map modulo p (S39; row J of local_systems_dictionary_20261017.md).

F = (u^3 z + y^2 u (4 + 3xy), y + 3x u^2 z + 3x y^2 (4 + 3xy), 2x - 3x^2 y - x^3 z), u = 1 + xy, det JF = -2.
For an odd prime p: the fibre-count function N(w) = #{v in F_p^3 : F(v) = w}, its exact moments (THM-1310: N in {0,1,3},
Chebotarev proportions 1/3, 1/2, 1/6 for the S_3 monodromy), and S(a) = sum_v e(a.F(v)/p) = sum_w N(w) e(a.w/p)
(a 3-dimensional DFT of N). Parseval: sum_(a != 0) |S(a)|^2 = p^3 #{(v, v') : v != v', F(v) = F(v')} (checked exactly in
integers on the collision side, in floating point on the Fourier side). An automorphism would have S(a) = 0 for all
a != 0; here the defect V = F_* 1 - 1 (the standard representation of S_3 on the generic fibre) has Fourier energy
p^3 <V, V> ~ p^6.
Usage: python keller_thm1300_fourier_modp_20261017.py
Session opus-2026-10-17-S39.
"""
import numpy as np
from numba import njit


@njit(cache=True)
def fibre_counts(p):
    N = np.zeros((p, p, p), dtype=np.int64)
    for x in range(p):
        for y in range(p):
            u = (1 + x * y) % p
            for z in range(p):
                f1 = (u * u % p * u % p * z + y * y % p * u % p * ((4 + 3 * x * y) % p)) % p
                f2 = (y + 3 * x % p * u % p * u % p * z + 3 * x % p * y % p * y % p * ((4 + 3 * x * y) % p)) % p
                f3 = (2 * x - 3 * x * x % p * y - x * x % p * x % p * z) % p
                N[f1, f2, f3] += 1
    return N


if __name__ == '__main__':
    print(' p   frac N=0   N=1   N=2   N=3  other   collisions/p^3   Parseval rel.err   E|S|^2/p^3   E|S|^4/(E|S|^2)^2   max|S|/p^1.5')
    for p in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41]:
        N = fibre_counts(p)
        tot = p ** 3
        fr = [np.count_nonzero(N == k) / tot for k in range(4)]
        other = np.count_nonzero(N > 3) / tot
        coll = int((N * (N - 1)).sum())
        S = np.fft.fftn(N.astype(np.float64))
        S2 = np.abs(S) ** 2
        S2[0, 0, 0] = 0.0
        pars = abs(S2.sum() - tot * coll) / (tot * coll)
        m2 = S2.sum() / (tot - 1)
        m4 = (S2 ** 2).sum() / (tot - 1)
        print('%3d  %.4f  %.4f  %.4f  %.4f  %.4f   %.4f          %.1e           %.4f        %.4f            %.3f' %
              (p, fr[0], fr[1], fr[2], fr[3], other, coll / tot, pars, m2 / tot, m4 / m2 ** 2, np.sqrt(S2.max()) / p ** 1.5))
        assert pars < 1e-9
    # where the Fourier energy lives: the cone over the projective dual of the Jelonek surface's curve at infinity.
    # Jelonek set (THM-4621): 27a^2c^2 - 18abc + 16a + b^3c - b^2 = 0; top degree c (27 a^2 c + b^3); the cuspidal cubic
    # 27 a^2 c + b^3 = 0 is parametrized by (1 : -3t : t^3), its tangent lines are (2t^3 : t^2 : 1): the dual cubic
    # alpha^2 gamma = 4 beta^3 (it contains (0 : 0 : 1), the dual of the line c = 0).
    print('\n p   #{|S| > p^2/2}  = #{alpha^2 gamma = 4 beta^3} \\ 0   sets equal   energy share on cone   max off-cone |S|/p^1.5   min on-cone |S|/p^2')
    for p in [11, 13, 17, 19, 23, 29, 31, 37, 41]:
        N = fibre_counts(p)
        S = np.abs(np.fft.fftn(N.astype(np.float64)))
        S[0, 0, 0] = 0.0
        a = np.arange(p)
        al, be, ga = np.meshgrid(a, a, a, indexing='ij')
        cone = ((al * al % p) * ga - 4 * (be ** 3 % p)) % p == 0
        cone[0, 0, 0] = False
        big = S > p ** 2 / 2
        eq = bool((big == cone).all())
        share = (S[cone] ** 2).sum() / (S ** 2).sum()
        print('%3d  %8d        %8d                  %s        %.4f                 %.3f                  %.4f' %
              (p, int(big.sum()), int(cone.sum()), eq, share, S[~cone].max() / p ** 1.5, S[cone].min() / p ** 2))
        assert eq and cone.sum() == p * p - 1

    # exact decomposition for alpha != 0: S = M - E + S_Gamma, M = p e(-4 beta^2/alpha) (p [alpha^2 gamma = 4 beta^3] - 1),
    # E = sum_(c != 0, w) e(alpha q(w)/c^2 + beta p(w)/c + gamma c), q = -w^3 + w^2, p(w) = -3w^2 + 4w (the excluded critical b),
    # S_Gamma = sum over the surface 2 - 3xy - x^2 z = 0, x != 0 (where F_3 = 0) of e(alpha F_1 + beta F_2).
    print('\nexact decomposition S = M - E + S_Gamma (alpha != 0), every frequency, and the sizes of the two remainders')
    for p in [5, 7, 11, 13]:
        N = fibre_counts(p)
        S = np.fft.fftn(N.astype(np.float64))         # S(a) = sum_w N(w) e(-a.w/p) (numpy sign convention)
        ep = lambda t: np.exp(-2j * np.pi * (t % p) / p)
        inv = [0] + [pow(c, -1, p) for c in range(1, p)]
        worst = 0.0; maxE = maxG = 0.0
        for al in range(1, p):
            for be in range(p):
                for ga in range(p):
                    onc = (al * al * ga - 4 * be ** 3) % p == 0
                    M = p * ep(-4 * be * be * inv[al]) * (p * onc - 1)
                    E = 0j
                    for c in range(1, p):
                        ic = inv[c]
                        for w in range(p):
                            qw = (-w ** 3 + w * w) % p; pw = (-3 * w * w + 4 * w) % p
                            E += ep(al * qw * ic * ic + be * pw * ic + ga * c)
                    G = 0j
                    for x in range(1, p):
                        ix = inv[x]
                        for y in range(p):
                            z = (2 - 3 * x * y) * ix * ix % p
                            u = (1 + x * y) % p
                            f1 = (u ** 3 * z + y * y * u * (4 + 3 * x * y)) % p
                            f2 = (y + 3 * x * u * u * z + 3 * x * y * y * (4 + 3 * x * y)) % p
                            G += ep(al * f1 + be * f2)
                    worst = max(worst, abs(S[al, be, ga] - (M - E + G)))
                    maxE = max(maxE, abs(E)); maxG = max(maxG, abs(G))
        print('p = %2d: max |S - (M - E + S_Gamma)| = %.1e;  max|E|/p = %.3f, max|S_Gamma|/p = %.3f' % (p, worst, maxE / p, maxG / p))
        assert worst < 1e-6 * p ** 3
    print('ALL CHECKS PASSED')
