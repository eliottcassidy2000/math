#!/usr/bin/env python3
"""Exact signed block moments of the Collatz parity through N = 2^LOGN (default 37), in the format of
05-knowledge/results/collatz_boundary_moments_20261013c.npz (renewal_shapes_oct11, N = 2^34), so that the conditional
Mercator-pole certificate (Theorem B3 of collatz_amplitude_boundary_20261013c.md) can be evaluated at N = 2^37.

moments[j] = (sum t(n), sum (n - jB) t(n)) over n in [jB, (j+1)B), B = 2^17, with t(0) := 0 and t(n) = (-1)^tau(n).
Also S(N) = sum_(n<=N) t(n) and H(N) = sum_(n<=N) (3[n = 2 mod 3] - 1) t(n), inclusive of n = N (tau(2^k) = k).
Cross-checks: the first 2^17 rows reproduce the stored 2^34 moments exactly; S(2^34) = 24953018, H(2^34) = 26869.
Parity bits come from collatz_census_extension_2p37_20261015.py (uint64 descent with exact big-integer fallback).
Session opus-2026-10-15-S37.  Usage: python collatz_boundary_moments_2p37_20261015.py LOGN out.npz out.json [old34.npz]
"""
import sys, os, json, time
import importlib.util
import numpy as np
import numba as nb

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('ext', os.path.join(HERE, 'collatz_census_extension_2p37_20261015.py'))
ext = importlib.util.module_from_spec(spec)
_argv = sys.argv
sys.argv = ['x', '30', 'unused.npz']
spec.loader.exec_module(ext)
sys.argv = _argv
LOGN = int(sys.argv[1]) if len(sys.argv) > 1 else 37
OUTN = sys.argv[2] if len(sys.argv) > 2 else 'moments.npz'
OUTJ = sys.argv[3] if len(sys.argv) > 3 else 'moments.json'
OLD = sys.argv[4] if len(sys.argv) > 4 else None
B = 1 << 17


@nb.njit(parallel=True)
def moments(P, nblocks, Bs):
    v = np.zeros((nblocks, 2), dtype=np.int64)
    c3 = np.zeros((nblocks, 3), dtype=np.int64)
    for j in nb.prange(nblocks):
        lo = j * Bs
        s0 = 0
        s1 = 0
        a = np.zeros(3, dtype=np.int64)
        for n in range(max(lo, 1), lo + Bs):
            t = 1 - 2 * int((P[n >> 6] >> np.uint64(n & 63)) & np.uint64(1))
            s0 += t
            s1 += (n - lo) * t
            a[n % 3] += t
        v[j, 0] = s0
        v[j, 1] = s1
        for r in range(3):
            c3[j, r] = a[r]
    return v, c3


if __name__ == '__main__':
    t0 = time.time()
    P, patched = ext.parity_bits(LOGN)
    print(f'parity bits to 2^{LOGN}: {time.time() - t0:.1f} s; {len(patched)} words recomputed exactly', flush=True)
    N = 1 << LOGN
    v, c3 = moments(P, N // B, B)
    tN = -1 if LOGN % 2 else 1
    S = int(v[:, 0].sum()) + tN
    C3 = c3.sum(axis=0).astype(object)
    C3[N % 3] += tN
    H = 3 * int(C3[2]) - S
    k34 = (1 << 34) // B
    S34 = int(v[:k34, 0].sum()) + 1
    C34 = c3[:k34].sum(axis=0)
    H34 = 3 * int(C34[2]) - S34      # n = 2^34 is 1 mod 3 with t = +1, so it enters S34 only
    print(f'S(2^{LOGN}) = {S}, H(2^{LOGN}) = {H}; S(2^34) = {S34} (24953018), H(2^34) = {H34} (26869)', flush=True)
    ok_old = None
    if OLD:
        old = np.load(OLD)['moments']
        ok_old = bool(np.array_equal(old, v[:len(old)]))
        print(f'first {len(old)} rows equal to the stored 2^34 moments: {ok_old}', flush=True)
    np.savez_compressed(OUTN, moments=v, N=np.int64(N), B=np.int64(B))
    json.dump(dict(N=N, logN=LOGN, t_N=tN, S=S, H=H, S_2p34=S34, H_2p34=H34, rows_equal_2p34_moments=ok_old,
                   patched_words=len(patched), seconds=time.time() - t0), open(OUTJ, 'w'), indent=1)
    print(f'wrote {OUTN}, {OUTJ}; total {time.time() - t0:.1f} s', flush=True)
