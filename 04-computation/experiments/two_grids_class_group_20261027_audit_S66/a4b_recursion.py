"""Check the critical-case recursion of the local densities and the stability of the 3-adic values at 3^5 for 81, 162.

sigma_p(n) = sigma_p^unit(n) + sigma_p(n / pi^4-part):  at 2, all-even tuples give sigma_2(16 m) = sigma_2^u(16 m) + sigma_2(m);
at 3, all tuples in (1 - w) give sigma_3(9 m) = sigma_3^u(9 m) + sigma_3(m)  (the scaling factors cancel exactly: 8 real
coordinates, two equations of degree 4 -- the critical dimension).  sigma^u counts tuples with at least one unit.
"""
import numpy as np


def counts(M, nonunit_only, p):
    a = np.arange(M, dtype=np.int64)
    A, B = np.meshgrid(a, a, indexing='ij')
    A, B = A.ravel(), B.ravel()
    if nonunit_only:
        keep = ((A + B) % 3 == 0) if p == 3 else ((A % 2 == 0) & (B % 2 == 0))
        A, B = A[keep], B[keep]
    s, t = (A * A - B * B) % M, (2 * A * B - B * B) % M
    x, y = (s * s - t * t) % M, (2 * s * t - t * t) % M
    D = np.bincount(x * M + y, minlength=M * M).astype(np.int64)
    idx = np.flatnonzero(D)
    w = D[idx]
    sx, sy = idx // M, idx % M
    C = np.zeros(M * M, dtype=np.float64)
    step = max(1, 4_000_000 // len(idx))
    for i0 in range(0, len(idx), step):
        X = (sx[i0:i0 + step, None] + sx[None, :]) % M
        Y = (sy[i0:i0 + step, None] + sy[None, :]) % M
        W = (w[i0:i0 + step, None] * w[None, :]).astype(np.float64)
        C += np.bincount((X * M + Y).ravel(), weights=W.ravel(), minlength=M * M)
    C = C.round().astype(np.int64).reshape(M, M)
    return C


def sigma(C, M, n):
    Qx = (n - np.arange(M)) % M
    Qy = (-np.arange(M)) % M
    S = C[np.ix_(Qx, Qy)]
    tot = sum(int(u) * int(v) for u, v in zip(C.ravel().tolist(), S.ravel().tolist()) if u and v)
    return tot / float(M) ** 6


for p, Ms, ns in ((3, (81, 243), (1, 2, 9, 18, 81, 162, 27, 54)), (2, (128, 256, 512), (1, 2, 3, 4, 16, 32, 64, 128))):
    for M in Ms:
        Call = counts(M, False, p)
        Cnu = counts(M, True, p)
        out = []
        for n in ns:
            s_all = sigma(Call, M, n)
            s_nu = sigma(Cnu, M, n)
            out.append(f'{n}: total {s_all:.6f}, unit {s_all - s_nu:.6f}, all-nonunit {s_nu:.6f}')
        print(f'p = {p}, modulus {M}:')
        for o in out:
            print('   ', o)
