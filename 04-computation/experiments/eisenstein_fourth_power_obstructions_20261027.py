"""Local obstructions and local densities for sums of four Eisenstein fourth powers (HYP-9429; opus-2026-10-27-S66).

Companion to 05-knowledge/results/two_grids_ring_or_ideal_class_13_37_20261027.md, section 5.

Elements of Z[omega]/M are pairs (x, y) for x + y omega.  For each modulus M (prime powers 2^4..2^8, 3^3..3^5, 5^2, 5^3,
7^2) the set of fourth-power residues is computed exactly, and the rational residues n (mod M) that are sums of four of
them (zeros allowed) are found by iterated set addition.  By the Chinese remainder theorem composite moduli add
nothing.  Then the 2-adic and 3-adic local densities sigma_p(n) = #{4-tuples mod M with sum = n} / M^6 are computed by
four-fold convolution (FFT on the M x M torus; exact counts rounded), and the persistent misses of HYP-9429 are ranked.
"""
import time

import numpy as np

MISSES = [6, 54, 60, 75, 78, 124, 140]      # no representation with bases <= 120 (6: <= 150)
FILLED = [102, 186]                          # represented with bases <= 120 but not <= 60


def fourth_residues(M):
    a = np.arange(M, dtype=np.int64)
    A, B = np.meshgrid(a, a, indexing='ij')
    A, B = A.ravel(), B.ravel()

    def mul(x1, y1, x2, y2):
        return (x1 * x2 - y1 * y2) % M, (x1 * y2 + y1 * x2 - y1 * y2) % M
    s1, s2 = mul(A, B, A, B)
    return mul(s1, s2, s1, s2)


def rational_residues(M):
    f1, f2 = fourth_residues(M)
    S = np.unique(f1 * M + f2)
    reach = np.zeros(M * M, dtype=bool)
    reach[0] = True
    Sx, Sy = S // M, S % M
    for _ in range(4):
        cur = np.flatnonzero(reach)
        cx, cy = cur // M, cur % M
        nx = (cx[:, None] + Sx[None, :]) % M
        ny = (cy[:, None] + Sy[None, :]) % M
        reach[(nx * M + ny).ravel()] = True
    return {int(n) for n in range(M) if reach[n * M]}, len(S)


def density(M):
    f1, f2 = fourth_residues(M)
    D = np.zeros((M, M))
    np.add.at(D, (f1, f2), 1)
    C = np.rint(np.real(np.fft.ifft2(np.fft.fft2(D) ** 4)))
    assert C.sum() == float(M) ** 8
    return C / float(M) ** 6


def main():
    t0 = time.time()
    f1, f2 = fourth_residues(16)
    res16 = sorted(set(zip(f1.tolist(), f2.tolist())))
    assert res16 == [(0, 0), (0, 1), (0, 9), (1, 0), (7, 7), (9, 0), (15, 15)]
    print('fourth powers mod 16:', res16, '= {0} u mu_3 u 9 mu_3')
    for M in (16, 32, 64, 128, 256, 27, 81, 243, 25, 125, 49):
        rr, ns = rational_residues(M)
        bad = [n for n in range(201) if n % M not in rr]
        if M % 2 == 0:
            assert bad == [n for n in range(201) if n % 8 == 5], M
        else:
            assert bad == [], M
        print(f'  M = {M:3d}: {ns:5d} fourth-power residues; obstructed n in 0..200: '
              f'{"exactly n = 5 mod 8" if bad else "none"}')
    s2, s3 = density(128), density(243)
    # Hensel: a tuple with a unit term lifts once the precision exceeds 2 v(4 x^3) = 4 (at 2) or 0 (at 1 - omega);
    # all-non-unit tuples only reach n = 0 (mod 16), resp. n = 0 (mod 9).  So the density of n is stable from 2^5 for
    # n != 0 (mod 16), and from 3^3 for n != 0 (mod 9); deeper n need deeper moduli (54 = 2 * 27 is checked to 3^5).
    d32, d64, d27, d81 = density(32), density(64), density(27), density(81)
    for n in range(1, 201):
        if n % 16:
            assert abs(d32[n % 32, 0] - s2[n % 128, 0]) < 1e-9 and abs(d64[n % 64, 0] - s2[n % 128, 0]) < 1e-9, n
        if n % 9:
            assert abs(d27[n % 27, 0] - s3[n % 243, 0]) < 1e-9 and abs(d81[n % 81, 0] - s3[n % 243, 0]) < 1e-9, n
    for n in MISSES + FILLED:
        print(f'  n = {n:3d}: 2-adic density mod 32/64/128 = {d32[n % 32, 0]:.4f}/{d64[n % 64, 0]:.4f}/{s2[n % 128, 0]:.4f}; '
              f'3-adic mod 27/81/243 = {d27[n % 27, 0]:.4f}/{d81[n % 81, 0]:.4f}/{s3[n % 243, 0]:.4f}')
    allowed = [n for n in range(1, 201) if n % 8 != 5]
    prod = {n: s2[n % 128, 0] * s3[n % 243, 0] for n in allowed}
    order = sorted(allowed, key=lambda n: (prod[n], n))
    rank = {n: i + 1 for i, n in enumerate(order)}
    med = float(np.median([prod[n] for n in allowed]))
    print(f'local densities (2-adic mod 128 x 3-adic mod 243), n = 1..200, n != 5 mod 8: median {med:.4f}, '
          f'minimum {prod[order[0]]:.4f}')
    for n in MISSES + FILLED:
        tag = 'miss' if n in MISSES else 'filled at B <= 120'
        print(f'  n = {n:3d} ({tag}): density {prod[n]:.4f}, rank {rank[n]} of {len(allowed)}')
    low = [n for n in allowed if abs(prod[n] - prod[order[0]]) < 1e-9]
    print(f'  the lowest density class ({prod[order[0]]:.4f}): {low}')
    print(f'ALL CHECKS PASSED  [{time.time() - t0:.0f} s]')


if __name__ == '__main__':
    main()
