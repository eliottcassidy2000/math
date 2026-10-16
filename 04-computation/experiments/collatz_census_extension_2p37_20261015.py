#!/usr/bin/env python3
"""Out-of-sample extension of the Collatz parity census to 2^LOGMAX (default 37): the pre-registered tests of HYP-9282
(octaves 2^34..2^36) and the window-energy predictions of HYP-9283.

T(n) = n/2 or (3n+1)/2, tau(n) = number of T-steps from n to 1, t(n) = (-1)^tau(n), H(N) = sum_(n<=N) w_n t(n) with
w_n = 3[n = 2 mod 3] - 1, so H = 3 C(2 mod 3) - S.

Pass 1 (parallel): parity bits p(n) = tau(n) mod 2 for 1 <= n < 2^LOGMAX by descent below the octave start and lookup,
         as in collatz_octave_memory_census_20261015.py, except that an orbit which would leave uint64 is not aborted:
         its 64-bit word is flagged and recomputed exactly with Python integers before the next octave starts.
Pass 2 (parallel over blocks of 2^20): class counts of t modulo 243 per block; exact prefix over blocks.
Pass 3 (parallel over blocks): walk each block from its exact starting class sums; record at every grid point
         N_j = round(2^(10 + j/64)) < 2^LOGMAX the class sums modulo 243 (inclusive of N_j), the value H(N_j) and the
         within-block prefix sums of H and H^2; block totals of H and H^2 (exact int64; the script checks
         2^20 max|H|^2 < 2^63 so no block sum can overflow).
Output NPZ: grid, class sums (G x 243), H at grid points, exact prefix sums of H and H^2 at grid points (as float64 and
as decimal strings), block totals.  Cross-checks: S(2^34) = 24953018, H(2^34) = 26869 (census of 2026-10-15),
spot check of random n against direct trajectories with Python integers.
Session opus-2026-10-15-S37.  Usage: python collatz_census_extension_2p37_20261015.py LOGMAX out.npz
"""
import sys, time, json, random
import numpy as np
import numba as nb

LOGMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 37
OUT = sys.argv[2] if len(sys.argv) > 2 else 'census_ext.npz'
GUARD = np.uint64((2**64 - 2) * 2 // 3 - 1)   # odd x <= GUARD keeps x + (x+1)//2 below 2^64
M = 243
BLOCK = 1 << 20
CHUNK_WORDS = 1 << 24


@nb.njit(parallel=True)
def parity_chunk(P, w0, w1, lo, guard, flags):
    # bits for n in words [w0, w1) of the octave starting at lo; lookups only below lo
    for w in nb.prange(w0, w1):
        word = np.uint64(0)
        for b in range(64):
            n = (w << 6) + b
            x = np.uint64(n)
            s = 0
            hit = False
            while x >= np.uint64(lo):
                if x & np.uint64(1):
                    if x > guard:
                        hit = True
                        break
                    x = x + (x + np.uint64(1)) // np.uint64(2)
                else:
                    x = x >> np.uint64(1)
                s += 1
            if hit:
                flags[w - w0] = 1
                continue
            par = (s + int((P[int(x) >> 6] >> np.uint64(int(x) & 63)) & np.uint64(1))) & 1
            if par:
                word |= np.uint64(1) << np.uint64(b)
        P[w] = word


def patch_word(P, w, lo):
    word = 0
    for b in range(64):
        n = (w << 6) + b
        x, s = n, 0
        while x >= lo:
            x = x // 2 if x % 2 == 0 else (3 * x + 1) // 2
            s += 1
        par = (s + ((int(P[x >> 6]) >> (x & 63)) & 1)) & 1
        word |= par << b
    P[w] = np.uint64(word)


def parity_bits(logmax):
    P = np.zeros((1 << logmax) // 64 + 1, dtype=np.uint64)
    word = 0
    for n in range(1, 64):
        x, s = n, 0
        while x != 1:
            x = x // 2 if x % 2 == 0 else (3 * x + 1) // 2
            s += 1
        if s & 1:
            word |= 1 << n
    P[0] = np.uint64(word)
    patched = []
    lo = 64
    while lo < (1 << logmax):
        hi = min(2 * lo, 1 << logmax)
        w0, w1 = lo >> 6, hi >> 6
        for c0 in range(w0, w1, CHUNK_WORDS):
            c1 = min(c0 + CHUNK_WORDS, w1)
            flags = np.zeros(c1 - c0, dtype=np.uint8)
            parity_chunk(P, c0, c1, lo, GUARD, flags)
            for i in np.nonzero(flags)[0]:
                patch_word(P, c0 + int(i), lo)
                patched.append(c0 + int(i))
        lo = hi
    return P, patched


@nb.njit(parallel=True)
def block_counts(P, nblocks, B, Mod):
    cnt = np.zeros((nblocks, Mod), dtype=np.int64)
    for b in nb.prange(nblocks):
        lo = b * B
        n0 = lo if lo > 0 else 1
        r = n0 % Mod
        for n in range(n0, lo + B):
            t = 1 - 2 * int((P[n >> 6] >> np.uint64(n & 63)) & np.uint64(1))
            cnt[b, r] += t
            r += 1
            if r == Mod:
                r = 0
    return cnt


@nb.njit(parallel=True)
def block_pass(P, nblocks, B, Mod, start, grid, glo, ghi):
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
            t = 1 - 2 * int((P[n >> 6] >> np.uint64(n & 63)) & np.uint64(1))
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


def tau_parity(n):
    s = 0
    while n != 1:
        n = n // 2 if n % 2 == 0 else (3 * n + 1) // 2
        s += 1
    return s & 1


if __name__ == '__main__':
    t0 = time.time()
    P, patched = parity_bits(LOGMAX)
    print(f'parity bits to 2^{LOGMAX}: {time.time() - t0:.1f} s; {len(patched)} words recomputed with Python integers',
          flush=True)
    rng = random.Random(20261015)
    bad = 0
    for _ in range(20000):
        n = rng.randrange(1, 1 << LOGMAX)
        if ((int(P[n >> 6]) >> (n & 63)) & 1) != tau_parity(n):
            bad += 1
    for w in patched[:200]:
        for b in range(0, 64, 7):
            n = (w << 6) + b
            if ((int(P[n >> 6]) >> (n & 63)) & 1) != tau_parity(n):
                bad += 1
    print(f'spot check (20000 random n + patched words): {bad} mismatches', flush=True)
    if bad:
        raise SystemExit('parity mismatch')
    nblocks = (1 << LOGMAX) // BLOCK
    t1 = time.time()
    cnt = block_counts(P, nblocks, BLOCK, M)
    start = np.zeros_like(cnt)
    start[1:] = np.cumsum(cnt[:-1], axis=0)
    print(f'block counts: {time.time() - t1:.1f} s', flush=True)
    raw = sorted(set(int(round(2 ** (10 + j / 64))) for j in range(0, (LOGMAX - 10) * 64 + 1)))
    grid = np.array([g for g in raw if g < (1 << LOGMAX)], dtype=np.int64)
    bidx = grid // BLOCK
    glo = np.searchsorted(bidx, np.arange(nblocks), side='left').astype(np.int64)
    ghi = np.searchsorted(bidx, np.arange(nblocks), side='right').astype(np.int64)
    t2 = time.time()
    snapC, snapH1, snapH2, snapH, totH1, totH2, maxabs = block_pass(P, nblocks, BLOCK, M, start, grid, glo, ghi)
    mx = int(maxabs.max())
    print(f'block pass: {time.time() - t2:.1f} s; max |H| = {mx}', flush=True)
    if BLOCK * mx * mx >= 2 ** 63:
        raise SystemExit('possible int64 overflow in block sums of H^2')
    cum1 = [0] * (nblocks + 1)
    cum2 = [0] * (nblocks + 1)
    for b in range(nblocks):
        cum1[b + 1] = cum1[b] + int(totH1[b])
        cum2[b + 1] = cum2[b] + int(totH2[b])
    pre1 = [cum1[int(bidx[i])] + int(snapH1[i]) for i in range(len(grid))]
    pre2 = [cum2[int(bidx[i])] + int(snapH2[i]) for i in range(len(grid))]
    S = snapC.sum(axis=1)
    i34 = int(np.nonzero(grid == (1 << 34))[0][0]) if (1 << 34) in set(grid.tolist()) else None
    if i34 is not None:
        H34 = 3 * int(snapC[i34, 2::3].sum()) - int(S[i34])
        print(f'cross-check: S(2^34) = {int(S[i34])} (census 24953018), H(2^34) = {H34} (census 26869), '
              f'H from block walk {int(snapH[i34])}', flush=True)
    # exact dyadic octave sums of H and H^2 from block totals (octaves k >= 20)
    octs = {}
    for k in range(20, LOGMAX):
        b0, b1 = (1 << k) // BLOCK, (1 << (k + 1)) // BLOCK
        octs[k] = (str(cum1[b1] - cum1[b0]), str(cum2[b1] - cum2[b0]))
    np.savez_compressed(OUT, grid=grid, classsums=snapC, H=snapH,
                        preH1=np.array([float(x) for x in pre1]), preH2=np.array([float(x) for x in pre2]),
                        preH1_str=np.array([str(x) for x in pre1]), preH2_str=np.array([str(x) for x in pre2]),
                        totH1=totH1, totH2=totH2, maxabs=maxabs, block=BLOCK, logmax=LOGMAX,
                        patched=np.array(patched, dtype=np.int64))
    json.dump(dict(logmax=LOGMAX, patched_words=len(patched), max_abs_H=mx,
                   octave_sums_H_H2={str(k): v for k, v in octs.items()}, seconds=time.time() - t0),
              open(OUT.replace('.npz', '.json'), 'w'), indent=1)
    print(f'wrote {OUT}; total {time.time() - t0:.1f} s', flush=True)
