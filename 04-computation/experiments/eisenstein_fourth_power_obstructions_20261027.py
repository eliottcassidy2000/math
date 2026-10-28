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

MISSES = [60, 75, 78, 124, 140]             # no representation with bases <= 650 (independent audit)
LARGE = [6, 54]                              # representable, but only with bases beyond 540 (audit: 6 at 541, 54 = 9 * 6)
FILLED = [102, 186]                          # represented with bases <= 120 but not <= 60


def fourth_residues(M):
    a = np.arange(M, dtype=np.int64)
    A, B = np.meshgrid(a, a, indexing='ij')
    A, B = A.ravel(), B.ravel()

    def mul(x1, y1, x2, y2):
        return (x1 * x2 - y1 * y2) % M, (x1 * y2 + y1 * x2 - y1 * y2) % M
    s1, s2 = mul(A, B, A, B)
    return mul(s1, s2, s1, s2)


def torus_sum(A, B):
    """Boolean sumset of two subsets of (Z/M)^2 given as M x M indicator arrays (FFT convolution; counts are exact
    integers far below 2^53, so rounding at 0.5 is safe).  Memory O(M^2) instead of O(M^2 |S|)."""
    C = np.real(np.fft.ifft2(np.fft.fft2(A.astype(float)) * np.fft.fft2(B.astype(float))))
    return C > 0.5


def rational_residues(M):
    f1, f2 = fourth_residues(M)
    S = np.zeros((M, M), dtype=bool)
    S[f1, f2] = True                                  # includes 0 = 0^4, so "zeros allowed"
    reach = S.copy()
    for _ in range(3):
        reach = torus_sum(reach, S)
    return {int(n) for n in range(M) if reach[n, 0]}, int(S.sum())


def local_criterion_checks():
    """The finite inputs of THM-4826 (7).
    (a) At 2: for every residue n mod 32 with n != 5 (mod 8) there is a solution mod 2^5 of x1^4 + ... + x4^4 = n with a
        unit coordinate (Hensel then lifts it, since v(4 x^3) = 2 for a unit x).
    (b) At the residue fields F_3 (prime 1 - omega), F_7, F_13, F_19 (split primes) and F_25 (the inert prime 5): for every
        residue r of n in the prime field there is a solution of sum x_i^4 = r with some x_i != 0 (the smooth points that
        Hensel lifts away from 2).  Larger residue fields are covered by Weil's bound (q >= 19 if p does not divide n,
        q > 27 if it does)."""
    M = 32
    a = np.arange(M)
    A, B = np.meshgrid(a, a, indexing='ij')
    A, B = A.ravel(), B.ravel()
    f1, f2 = fourth_residues(M)
    unit = (A % 2 == 1) | (B % 2 == 1)               # a + b omega is a unit mod 2 iff not both a, b even
    S_all = np.zeros((M, M), dtype=bool)
    S_all[f1, f2] = True
    S_unit = np.zeros((M, M), dtype=bool)
    S_unit[f1[unit], f2[unit]] = True
    R3 = torus_sum(torus_sum(S_all, S_all), S_all)
    R = torus_sum(S_unit, R3)
    two_adic = all(bool(R[n, 0]) for n in range(M) if n % 8 != 5)
    two_adic_fails = [n for n in range(M) if not R[n, 0]]
    fields = {}
    # split primes p = 1 mod 3: residue field F_p; the fourth powers of F_p
    for p in (7, 13, 19):
        F = {pow(x, 4, p) for x in range(p)}
        nz = {pow(x, 4, p) for x in range(1, p)}
        ok = True
        for r in range(p):
            # some coordinate nonzero: one term from nz, three from F
            S3 = {(u + v + w) % p for u in F for v in F for w in F}
            ok &= any((t + s) % p == r for t in nz for s in S3)
        fields[f'F_{p}'] = ok
    # the prime 1 - omega: residue field F_3
    F3 = {pow(x, 4, 3) for x in range(3)}
    nz3 = {pow(x, 4, 3) for x in range(1, 3)}
    fields['F_3'] = all(any((t + u + v + w) % 3 == r for t in nz3 for u in F3 for v in F3 for w in F3) for r in range(3))
    # the inert prime 5: F_25 = F_5[s]/(s^2 - 2) (2 is a non-square mod 5); rational residues r in F_5
    def mul25(x, y):
        return ((x[0] * y[0] + 2 * x[1] * y[1]) % 5, (x[0] * y[1] + x[1] * y[0]) % 5)
    def p4_25(x):
        s = mul25(x, x)
        return mul25(s, s)
    els = [(u, v) for u in range(5) for v in range(5)]
    F25 = {p4_25(x) for x in els}
    nz25 = {p4_25(x) for x in els if x != (0, 0)}
    S3 = {((u[0] + v[0] + w[0]) % 5, (u[1] + v[1] + w[1]) % 5) for u in F25 for v in F25 for w in F25}
    fields['F_25'] = all(any(((t[0] + s[0]) % 5, (t[1] + s[1]) % 5) == (r, 0) for t in nz25 for s in S3) for r in range(5))
    return two_adic, two_adic_fails, fields


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
    two_adic, fails, fields = local_criterion_checks()
    assert two_adic and fails == [5, 13, 21, 29] and all(fields.values())
    print('local criterion THM-4826 (7): a unit-coordinate solution mod 2^5 exists for every n != 5 (mod 8) '
          f'(failing residues mod 32: {fails}); solutions with a nonzero coordinate for every rational residue over '
          f'{", ".join(k for k in fields)}: {all(fields.values())}')
    s2, s3 = density(256), density(243)      # 2^8: exact 2-adic factor for every n <= 200 (sigma_2(16m) needs it)
    # Hensel: a tuple with a unit term lifts once the precision exceeds 2 v(4 x^3) = 4 (at 2) or 0 (at 1 - omega);
    # all-non-unit tuples only reach n = 0 (mod 16), resp. n = 0 (mod 9).  So the density of n is stable from 2^5 for
    # n != 0 (mod 16), and from 3^3 for n != 0 (mod 9); deeper n need deeper moduli (54 = 2 * 27 is checked to 3^5).
    d32, d64, d27, d81, d128 = density(32), density(64), density(27), density(81), density(128)
    deep = [n for n in range(1, 201) if abs(d128[n % 128, 0] - s2[n % 256, 0]) > 1e-9]
    assert deep == [16, 32, 48, 64, 128, 144, 160, 176, 192], deep   # the mod-128 values are not limits only here
    print(f'  2-adic densities: mod 128 is not yet the limit exactly at n = {deep} (multiples of 16); mod 256 is used')
    for n in range(1, 201):
        if n % 16:
            assert abs(d32[n % 32, 0] - s2[n % 256, 0]) < 1e-9 and abs(d64[n % 64, 0] - s2[n % 256, 0]) < 1e-9, n
        if n % 9:
            assert abs(d27[n % 27, 0] - s3[n % 243, 0]) < 1e-9 and abs(d81[n % 81, 0] - s3[n % 243, 0]) < 1e-9, n
    for n in MISSES + LARGE + FILLED:
        print(f'  n = {n:3d}: 2-adic density mod 32/64/256 = {d32[n % 32, 0]:.4f}/{d64[n % 64, 0]:.4f}/{s2[n % 256, 0]:.4f}; '
              f'3-adic mod 27/81/243 = {d27[n % 27, 0]:.4f}/{d81[n % 81, 0]:.4f}/{s3[n % 243, 0]:.4f}')
    allowed = [n for n in range(1, 201) if n % 8 != 5]
    prod = {n: s2[n % 256, 0] * s3[n % 243, 0] for n in allowed}
    order = sorted(allowed, key=lambda n: (prod[n], n))
    rank = {n: i + 1 for i, n in enumerate(order)}
    med = float(np.median([prod[n] for n in allowed]))
    print(f'local densities (2-adic mod 256 x 3-adic mod 243), n = 1..200, n != 5 mod 8: median {med:.4f}, '
          f'minimum {prod[order[0]]:.4f}')
    for n in MISSES + LARGE + FILLED:
        tag = ('persistent miss (none with bases <= 650)' if n in MISSES else
               'representable only with large bases' if n in LARGE else 'filled at B <= 120')
        print(f'  n = {n:3d} ({tag}): density {prod[n]:.4f}, rank {rank[n]} of {len(allowed)}')
    low = [n for n in allowed if abs(prod[n] - prod[order[0]]) < 1e-9]
    print(f'  the lowest density class ({prod[order[0]]:.4f}): {low}')
    print(f'ALL CHECKS PASSED  [{time.time() - t0:.0f} s]')


if __name__ == '__main__':
    main()
