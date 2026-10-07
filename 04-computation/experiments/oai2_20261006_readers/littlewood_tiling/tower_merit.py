#!/usr/bin/env python3
"""Merit factors and sup norms of the rows/columns of the Sierpinski tower
H_2 = [[1,1],[-1,1]], H_{2n} = [[H,H],[-H^T,H^T]]  (skew Hadamard, order 2^k).

Exact: aperiodic autocorrelations via FFT (rounded to integers, checked),
merit factor F = N^2 / (2 sum_{u>=1} C_u^2).
Numerical: sup norm on the circle from an oversampled FFT (M = 16N) with the
Bernstein-type correction  ||P|| <= max_grid / cos(pi*N/(2M)).
"""
import numpy as np, sys, time

def tower(k):
    H = np.array([[1, 1], [-1, 1]], dtype=np.int8)
    for _ in range(k - 1):
        H = np.block([[H, H], [-H.T, H.T]]).astype(np.int8)
    return H

def autocorr_int(X):
    """X: (m, N) +-1 int8.  Returns (m, N) int64 aperiodic autocorrelations C_0..C_{N-1}."""
    m, N = X.shape
    L = 1 << int(np.ceil(np.log2(2 * N)))
    F = np.fft.rfft(X.astype(np.float64), n=L, axis=1)
    C = np.fft.irfft(np.abs(F) ** 2, n=L, axis=1)[:, :N]
    Ci = np.rint(C).astype(np.int64)
    err = np.max(np.abs(C - Ci))
    assert err < 1e-3, err
    return Ci

def merit(X):
    C = autocorr_int(X)
    N = X.shape[1]
    s = (C[:, 1:] ** 2).sum(axis=1)
    return N * N / (2.0 * s), s

def supnorm(X, over=16):
    m, N = X.shape
    M = over * N
    F = np.fft.fft(X.astype(np.float64), n=M, axis=1)
    a = np.abs(F)
    mx = a.max(axis=1)
    mn = a.min(axis=1)
    return mx / np.cos(np.pi * N / (2 * M)), mx, mn

def batched(fun, X, bs):
    outs = None
    for s in range(0, X.shape[0], bs):
        o = fun(X[s:s + bs])
        if not isinstance(o, tuple):
            o = (o,)
        if outs is None:
            outs = [[x] for x in o]
        else:
            for l, x in zip(outs, o):
                l.append(x)
    return [np.concatenate(l) for l in outs]

if __name__ == "__main__":
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    print("k N | rows: Fmax Fmean Fmin  argmax | cols: Fmax Fmean | rows: sup/sqrtN min..max (upper bd) | min|P| over grid (max over rows)")
    for k in range(1, kmax + 1):
        t = time.time()
        H = tower(k)
        N = H.shape[0]
        # sanity: Hadamard and skew
        if N <= 4096:
            G = H.astype(np.int64) @ H.T.astype(np.int64)
            assert np.array_equal(G, N * np.eye(N, dtype=np.int64))
            assert np.array_equal(H.astype(int) + H.T.astype(int), 2 * np.eye(N, dtype=int))
        bs = max(1, (1 << 22) // N)
        Fr, sr = batched(merit, H, bs)
        Fc, sc = batched(merit, np.ascontiguousarray(H.T), bs)
        if N >= 2:
            ub, mx, mn = batched(lambda X: supnorm(X, 16 if N <= 4096 else 8), H, max(1, bs // 16))
        else:
            ub = mx = mn = np.array([1.0])
        rt = np.sqrt(N)
        print(f"{k:2d} {N:6d} | {Fr.max():.4f} {Fr.mean():.4f} {Fr.min():.6f} {int(Fr.argmax()):6d} | "
              f"{Fc.max():.4f} {Fc.mean():.4f} | sup/sqrtN in [{(mx/rt).min():.4f}, {(mx/rt).max():.4f}] (ub {(ub/rt).min():.4f}) | "
              f"grid-min max {mn.max()/rt:.2e}  [{time.time()-t:.1f}s]", flush=True)
