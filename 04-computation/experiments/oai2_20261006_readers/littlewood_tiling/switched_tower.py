#!/usr/bin/env python3
"""(1) Closed form of the tower rows: r_i(y) = w_i(y) - 2 sum_{l: i_l=1, y = i mod 2^l} (-1)^{sum_{m>=l} i_m y_m}.
(2) Rudin-Shapiro switching D H D (D = diag(RS)): skew Hadamard; rows' merit factors and sup norms,
    compared with the PROVED bound sup <= (1+sqrt2)^2 sqrt(2N).
(3) Golay complementary pairs among rows/columns of the tower (small k, exhaustive)."""
import numpy as np, sys, time
from tower_merit import tower, autocorr_int, supnorm

def popcount(a):
    a = np.asarray(a, dtype=np.int64)
    c = np.zeros_like(a)
    while np.any(a):
        c += a & 1
        a >>= 1
    return c

def closed_form(k):
    N = 1 << k
    y = np.arange(N)
    H = np.zeros((N, N), dtype=np.int64)
    for i in range(N):
        row = (-1) ** popcount(i & y)
        for l in range(k):
            if (i >> l) & 1:
                mask = (y % (1 << l)) == (i % (1 << l))
                hi = (i >> l) << l
                row = row - 2 * mask * ((-1) ** popcount(hi & y))
        H[i] = row
    return H

def rs_seq(k):
    N = 1 << k
    y = np.arange(N, dtype=np.int64)
    return ((-1) ** popcount(y & (y >> 1))).astype(np.int8)

def merit_rows(X, bs):
    out = []
    for s in range(0, X.shape[0], bs):
        C = autocorr_int(X[s:s + bs])
        out.append((C[:, 1:] ** 2).sum(axis=1))
    s2 = np.concatenate(out)
    N = X.shape[1]
    return N * N / (2.0 * s2)

if __name__ == "__main__":
    # (1) closed form check
    for k in range(1, 9):
        assert np.array_equal(closed_form(k), tower(k).astype(np.int64)), k
    print("closed form r_i(y) = w_i(y) - 2 sum ... verified k = 1..8")
    # RS sequence check (Golay partner exists): via recursion
    P = np.array([1]); Q = np.array([1])
    for k in range(1, 12):
        P, Q = np.concatenate([P, Q]), np.concatenate([P, -Q])
        assert np.array_equal(P, rs_seq(k)), k
    print("RS(y) = (-1)^{sum y_m y_{m+1}} verified k <= 11")
    # (3) Golay pairs among rows and columns
    for k in range(1, 8):
        H = tower(k); N = H.shape[0]
        R = np.concatenate([H, H.T]).astype(np.int8)
        C = autocorr_int(R)[:, 1:]
        # pairs with C_a + C_b == 0 for all u: hash on C and -C
        d = {}
        for idx in range(R.shape[0]):
            d.setdefault(C[idx].tobytes(), []).append(idx)
        pairs = 0; ex = None
        for idx in range(R.shape[0]):
            key = (-C[idx]).tobytes()
            for j in d.get(key, []):
                if j > idx:
                    pairs += 1; ex = (idx, j)
        print(f"k={k} N={N}: Golay complementary pairs among the 2N rows+columns: {pairs}  example {ex}")
    # (2) switched tower
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    print("k N | switched rows: Fmin Fmean Fmax | sup/sqrtN max over rows (grid, ub) | bound (1+sqrt2)^2 sqrt2 = %.4f | min|P|/sqrtN min over rows" % ((1 + 2 ** 0.5) ** 2 * 2 ** 0.5))
    for k in range(1, kmax + 1):
        t = time.time()
        H = tower(k); N = H.shape[0]
        a = rs_seq(k)
        Hs = (a[:, None] * H * a[None, :]).astype(np.int8)
        if N <= 2048:
            assert np.array_equal(Hs.astype(int) + Hs.T.astype(int), 2 * np.eye(N, dtype=int))
            assert np.array_equal(Hs.astype(np.int64) @ Hs.T.astype(np.int64), N * np.eye(N, dtype=np.int64))
        bs = max(1, (1 << 22) // N)
        F = merit_rows(Hs, bs)
        sups = []; mins = []
        over = 16 if N <= 2048 else 8
        for s in range(0, N, max(1, bs // over)):
            ub, mx, mn = supnorm(Hs[s:s + max(1, bs // over)], over)
            sups.append(ub); mins.append(mn)
        ub = np.concatenate(sups); mn = np.concatenate(mins)
        rt = np.sqrt(N)
        print(f"{k:2d} {N:6d} | {F.min():.4f} {F.mean():.4f} {F.max():.4f} | {ub.max()/rt:.4f} | min-mod (grid) {mn.min()/rt:.2e} median-row-min {np.median(mn)/rt:.3f} [{time.time()-t:.1f}s]", flush=True)
