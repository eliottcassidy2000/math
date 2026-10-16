#!/usr/bin/env python3
"""Streaming extension of the exact Collatz parity census to 2^LOGMAX (default 40) with a stored table below 2^LOW.

T(n) = n/2 or (3n+1)/2, tau(n) = number of T-steps from n to 1, t(n) = (-1)^tau(n), H(N) = sum_(n<=N) w_n t(n) with
w_n = 3[n = 2 mod 3] - 1.  The parity bits below 2^LOW (default 37) come from collatz_census_extension_2p37_20261015.py
(uint64 descent with exact big-integer fallback, a 16 GB table at LOW = 37).  For 2^LOW <= n < 2^LOGMAX the parity is
computed on the fly: iterate T until the value drops below 2^LOW, then look up the table.  Values that would leave
uint64 continue in exact two-limb (128-bit) arithmetic; the script aborts if a value would leave 128 bits.

Pass A (parallel over blocks of 2^20): class counts of t modulo 243 per block, and the signed moments
           (sum t, sum (n - lo) t) of every sub-block of 2^17 (the format of collatz_boundary_moments_20261013c.npz).
Pass B (parallel over blocks): exact walk from the prefix class sums; at every grid point N_j = round(2^(10 + j/64))
           < 2^LOGMAX the class sums modulo 243, H(N_j) and within-block prefix sums of H and H^2; block totals.
Outputs: census NPZ (same keys as the 2^37 extension), moments NPZ (N = 2^LOGMAX, B = 2^17), JSON summary with exact
octave sums of H and H^2, S(N), H(N) (inclusive of n = N, tau(2^k) = k).
Checks: all grid points below 2^37 must reproduce the 2^37 extension exactly; 3000 random n in [2^LOW, 2^LOGMAX) are
compared with direct Python-integer trajectories.
Session opus-2026-10-15-S37.  Usage: python collatz_census_extension_2p40_20261015.py LOGMAX LOW out_census.npz
       out_moments.npz out.json [ext37.npz]
"""
import sys, os, json, time, random
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
LOGMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 40
LOW = int(sys.argv[2]) if len(sys.argv) > 2 else 37
OUTC = sys.argv[3] if len(sys.argv) > 3 else 'census40.npz'
OUTM = sys.argv[4] if len(sys.argv) > 4 else 'moments40.npz'
OUTJ = sys.argv[5] if len(sys.argv) > 5 else 'census40.json'
OLD = sys.argv[6] if len(sys.argv) > 6 else None
M = 243
BLOCK = 1 << 20
SUB = 1 << 17
GUARD = np.uint64((2**64 - 2) * 2 // 3 - 1)


@nb.njit
def tval(P, n, low, guard):
    # returns t(n) in {+1,-1}; 0 signals a 128-bit overflow (abort)
    if n < low:
        return 1 - 2 * int((P[n >> 6] >> np.uint64(n & 63)) & np.uint64(1))
    x = np.uint64(n)
    lo64 = np.uint64(low)
    s = 0
    while x >= lo64:
        if x & np.uint64(1):
            if x > guard:
                # two-limb continuation: value = hi * 2^64 + lo
                hi = np.uint64(0)
                lo = x
                while hi != np.uint64(0) or lo >= lo64:
                    if lo & np.uint64(1):
                        if hi >= np.uint64(1) << np.uint64(62):
                            return 0
                        # y = (x + 1) >> 1
                        l1 = lo + np.uint64(1)
                        h1 = hi + (np.uint64(1) if l1 == np.uint64(0) else np.uint64(0))
                        yl = (l1 >> np.uint64(1)) | ((h1 & np.uint64(1)) << np.uint64(63))
                        yh = h1 >> np.uint64(1)
                        # x = x + y
                        nl = lo + yl
                        carry = np.uint64(1) if nl < lo else np.uint64(0)
                        hi = hi + yh + carry
                        lo = nl
                    else:
                        lo = (lo >> np.uint64(1)) | ((hi & np.uint64(1)) << np.uint64(63))
                        hi = hi >> np.uint64(1)
                    s += 1
                x = lo
                break
            x = x + (x + np.uint64(1)) // np.uint64(2)
        else:
            x = x >> np.uint64(1)
        s += 1
    xi = int(x)
    par = (s + int((P[xi >> 6] >> np.uint64(xi & 63)) & np.uint64(1))) & 1
    return 1 - 2 * par


@nb.njit(parallel=True)
def pass_a(P, nblocks, B, Sb, Mod, low, guard):
    cnt = np.zeros((nblocks, Mod), dtype=np.int64)
    nsub = B // Sb
    mom = np.zeros((nblocks * nsub, 2), dtype=np.int64)
    bad = np.zeros(nblocks, dtype=np.int64)
    for b in nb.prange(nblocks):
        lo = b * B
        n0 = lo if lo > 0 else 1
        r = n0 % Mod
        for n in range(n0, lo + B):
            t = tval(P, n, low, guard)
            if t == 0:
                bad[b] = n
                t = 1
            cnt[b, r] += t
            k = (n - lo) // Sb
            mom[b * nsub + k, 0] += t
            mom[b * nsub + k, 1] += (n - lo - k * Sb) * t
            r += 1
            if r == Mod:
                r = 0
    return cnt, mom, bad


@nb.njit(parallel=True)
def pass_b(P, nblocks, B, Mod, start, grid, glo, ghi, low, guard):
    G = grid.shape[0]
    snapC = np.zeros((G, Mod), dtype=np.int64)
    snapH1 = np.zeros(G, dtype=np.int64)
    snapH2 = np.zeros(G, dtype=np.int64)
    snapH = np.zeros(G, dtype=np.int64)
    totH1 = np.zeros(nblocks, dtype=np.int64)
    totH2 = np.zeros(nblocks, dtype=np.int64)
    maxabs = np.zeros(nblocks, dtype=np.int64)
    for b in nb.prange(nblocks):
        C = start[b].copy()
        C3 = np.zeros(3, dtype=np.int64)
        for r in range(Mod):
            C3[r % 3] += C[r]
        S = C3[0] + C3[1] + C3[2]
        lo = b * B
        n0 = lo if lo > 0 else 1
        r = n0 % Mod
        s1 = 0
        s2 = 0
        mx = 0
        gi = glo[b]
        ge = ghi[b]
        for n in range(n0, lo + B):
            t = tval(P, n, low, guard)
            C[r] += t
            C3[r % 3] += t
            S += t
            h = 3 * C3[2] - S
            s1 += h
            s2 += h * h
            ah = h if h >= 0 else -h
            if ah > mx:
                mx = ah
            while gi < ge and grid[gi] == n:
                for a in range(Mod):
                    snapC[gi, a] = C[a]
                snapH1[gi] = s1
                snapH2[gi] = s2
                snapH[gi] = h
                gi += 1
            r += 1
            if r == Mod:
                r = 0
        totH1[b] = s1
        totH2[b] = s2
        maxabs[b] = mx
    return snapC, snapH1, snapH2, snapH, totH1, totH2, maxabs


if __name__ == '__main__':
    t0 = time.time()
    P, patched = ext.parity_bits(LOW)
    print(f'parity table to 2^{LOW}: {time.time() - t0:.1f} s; {len(patched)} words recomputed exactly', flush=True)
    low = 1 << LOW
    rng = random.Random(4040)
    bad = 0
    for _ in range(3000):
        n = rng.randrange(low, 1 << LOGMAX)
        if tval(P, n, low, GUARD) != 1 - 2 * ext.tau_parity(n):
            bad += 1
    print(f'spot check of 3000 random n in [2^{LOW}, 2^{LOGMAX}) against Python trajectories: {bad} mismatches', flush=True)
    if bad:
        raise SystemExit('parity mismatch')
    nblocks = (1 << LOGMAX) // BLOCK
    t1 = time.time()
    cnt, mom, badb = pass_a(P, nblocks, BLOCK, SUB, M, low, GUARD)
    if badb.any():
        raise SystemExit(f'128-bit overflow at n = {int(badb[np.nonzero(badb)[0][0]])}')
    start = np.zeros_like(cnt)
    start[1:] = np.cumsum(cnt[:-1], axis=0)
    print(f'pass A: {time.time() - t1:.1f} s', flush=True)
    raw = sorted(set(int(round(2 ** (10 + j / 64))) for j in range(0, (LOGMAX - 10) * 64 + 1)))
    grid = np.array([g for g in raw if g < (1 << LOGMAX)], dtype=np.int64)
    bidx = grid // BLOCK
    glo = np.searchsorted(bidx, np.arange(nblocks), side='left').astype(np.int64)
    ghi = np.searchsorted(bidx, np.arange(nblocks), side='right').astype(np.int64)
    t2 = time.time()
    snapC, snapH1, snapH2, snapH, totH1, totH2, maxabs = pass_b(P, nblocks, BLOCK, M, start, grid, glo, ghi, low, GUARD)
    mx = int(maxabs.max())
    print(f'pass B: {time.time() - t2:.1f} s; max |H| = {mx}', flush=True)
    if BLOCK * mx * mx >= 2 ** 63:
        raise SystemExit('possible int64 overflow in block sums of H^2')
    cum1 = [0] * (nblocks + 1)
    cum2 = [0] * (nblocks + 1)
    for b in range(nblocks):
        cum1[b + 1] = cum1[b] + int(totH1[b])
        cum2[b + 1] = cum2[b] + int(totH2[b])
    pre1 = [cum1[int(bidx[i])] + int(snapH1[i]) for i in range(len(grid))]
    pre2 = [cum2[int(bidx[i])] + int(snapH2[i]) for i in range(len(grid))]
    Nmax = 1 << LOGMAX
    tN = -1 if LOGMAX % 2 else 1
    final = (start[-1] + cnt[-1]).astype(object)
    final[Nmax % M] += tN
    S_N = int(sum(final))
    C2 = int(sum(final[2::3]))
    H_N = 3 * C2 - S_N
    print(f'S(2^{LOGMAX}) = {S_N}, H(2^{LOGMAX}) = {H_N}', flush=True)
    same = None
    if OLD:
        with np.load(OLD) as Z:
            g0 = np.array(Z['grid'])
            c0 = np.array(Z['classsums'])
            h0 = np.array(Z['H'])
            q0 = [int(x) for x in Z['preH2_str']]
        idx = np.searchsorted(grid, g0)
        same = bool(np.array_equal(grid[idx], g0) and np.array_equal(snapC[idx], c0) and np.array_equal(snapH[idx], h0)
                    and all(pre2[int(i)] == q for i, q in zip(idx, q0)))
        print(f'all {len(g0)} grid points of the 2^37 extension reproduced exactly: {same}', flush=True)
    octs = {}
    for k in range(20, LOGMAX):
        b0, b1 = (1 << k) // BLOCK, (1 << (k + 1)) // BLOCK
        octs[str(k)] = (str(cum1[b1] - cum1[b0]), str(cum2[b1] - cum2[b0]))
    np.savez_compressed(OUTC, grid=grid, classsums=snapC, H=snapH,
                        preH1=np.array([float(x) for x in pre1]), preH2=np.array([float(x) for x in pre2]),
                        preH1_str=np.array([str(x) for x in pre1]), preH2_str=np.array([str(x) for x in pre2]),
                        totH1=totH1, totH2=totH2, maxabs=maxabs, block=BLOCK, logmax=LOGMAX)
    np.savez_compressed(OUTM, moments=mom, N=np.int64(Nmax), B=np.int64(SUB))
    json.dump(dict(N=Nmax, logN=LOGMAX, low=LOW, t_N=tN, S=S_N, H=H_N, max_abs_H=mx, reproduces_2p37=same,
                   octave_sums_H_H2=octs, seconds=time.time() - t0), open(OUTJ, 'w'), indent=1)
    print(f'wrote {OUTC}, {OUTM}, {OUTJ}; total {time.time() - t0:.1f} s', flush=True)
