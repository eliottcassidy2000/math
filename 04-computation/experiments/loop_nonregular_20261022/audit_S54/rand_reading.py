"""Auditor: THM-4761 section 2 table and the sample check, recomputed (my code) on (a) the lane's rand12/rand14
files and (b) my own uniform random tournaments (random12_s2026.rows: 20000 samples, random14_s2027.rows: 3000),
with exact H from forest_audit.exe.  Also checks the identity W = n(n^2 - 1)/3 - 8 c_3 (c_3 = cyclic triples), i.e.
kappa_2 = 8 c_3/(n(n-1)) - (n-2)/3, and the empirical corr(kappa, 1 + kappa_2) against sqrt(share).
Usage: python rand_reading.py > rand_reading.out"""
import os, sys, math
from fractions import Fraction as Fr
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tourn
HERE = os.path.dirname(os.path.abspath(__file__))
UP = os.path.join(HERE, '..')


def analyse(label, codes, Hs, n):
    lk, lr, T1, W, k2, kap = [], [], [], [], [], []
    for code, H in zip(codes, Hs):
        A = tourn.from_code(code); S = (A - A.T).astype(float); M = S / n
        x = np.linalg.eigvals(M).imag
        logrho = 0.5 * np.sum(np.log((1 + x ** 2) / (1 - x ** 2)))
        delta = M @ np.ones(n)
        t1 = -delta @ np.linalg.solve(np.eye(n) + M @ M, delta)
        w = float(np.sum((S.sum(1)) ** 2))
        s = A.sum(1)
        c3 = math.comb(n, 3) - sum(math.comb(int(v), 2) for v in s)
        assert abs(w - (n * (n * n - 1) / 3 - 8 * c3)) < 1e-9
        kp = H * 2 ** (n - 1) / math.factorial(n)
        lk.append(math.log(kp)); lr.append(logrho); T1.append(t1); W.append(w); kap.append(kp)
        k2.append(1 - w / (n * (n - 1)))
    lk, lr, T1, W, k2, kap = map(np.array, (lk, lr, T1, W, k2, kap))
    c = lambda a, b: np.corrcoef(a, b)[0, 1]
    # exact share
    V = sum(Fr(sum(2 ** p * math.prod(range(n - L - p + 1, n + 1)) * math.comb(L // 2 - 1, p - 1) // math.factorial(p)
                   for p in range(1, L // 2 + 1) if L + p <= n), math.prod(range(n - L + 1, n + 1)) ** 2)
            for L in range(2, n, 2))
    share = Fr(2 * (n - 2), n * (n - 1)) / V
    print(f'{label} (n={n}, {len(lk)} samples): corr(log k, log rho + T(1)) {c(lk, lr + T1):.4f}; corr(log k, -W) {c(lk, -W):.4f}; '
          f'corr(log k, log rho) {c(lk, lr):+.4f}; mean log(k/rho) - T(1) {np.mean(lk - lr - T1):+.4f} (sd {np.std(lk - lr - T1):.4f}); '
          f'mean kappa {kap.mean():.4f}, var {kap.var():.5f} (exact {float(V):.5f}); corr(kappa, 1 + kappa_2) {c(kap, k2):.4f} '
          f'vs sqrt(share) {math.sqrt(share):.4f}; identity W = n(n^2-1)/3 - 8 c3 holds on all samples')


def main():
    for fn, n in (('rand12_samples.txt', 12), ('rand14_samples.txt', 14)):
        rows = [l.split() for l in open(os.path.join(UP, fn))]
        analyse(f'lane file {fn}', [r[0] for r in rows], [int(r[1]) for r in rows], n)
    for fn, n in (('random12_s2026.rows', 12), ('random14_s2027.rows', 14)):
        rows = [l.split() for l in open(os.path.join(HERE, fn))]
        analyse(f'my samples {fn}', [r[0] for r in rows], [int(r[1]) for r in rows], n)


if __name__ == '__main__':
    main()
