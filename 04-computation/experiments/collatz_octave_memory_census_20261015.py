#!/usr/bin/env python3
"""Exact census of the Collatz mod-3 discrepancy H with octave means, Haar refinement energy and two-octave memory.

Session opus-2026-10-15-S37.  Notation of collatz_amplitude_boundary_20261013c.md / collatz_wavelet_energy_20261013d.md:
T(n) = n/2 or (3n+1)/2, tau(n) = number of T-steps from n to 1, t(n) = (-1)^tau(n),
psi_(1,r)(n) = 3[n = r mod 3] - 1,  psi_(2,r)(n) = 3[n = r mod 9] - [n = r mod 3],
H_(q,r)(N) = sum_(n <= N) psi_(q,r)(n) t(n),  H = H_(1,2),  S(N) = sum_(n <= N) t(n).

Exact identity used (from the refinement formula (2) of the wavelet note): for every integer N >= 0
    H(4N) = H(N) + H_(2,2)(3N) - H_(2,8)(6N),        so  Delta(N) := H(4N) - H(N)
is a level-2 quantity ('two-octave memory').

Pass 1 (parallel): parity bits p(n) = tau(n) mod 2 for n <= NMAX, by descent below the current octave and lookup;
         T is evaluated as x + (x+1)/2 on odd x in uint64 with an overflow guard (aborts rather than certifies).
Pass 2 (sequential, exact integers except the float64 squared sums):
  * class sums S, C3[a], C9[a], C27[a] -> all level-1/2/3 channels;
  * per dyadic octave [2^k, 2^(k+1)): sum H_c and sum H_c^2 for the 12 level-1/2 channels c;
  * exact streaming Haar decomposition of H on every octave: refinement (innovation) energy by level;
    check sum H^2 = X mean^2 + sum innovations;
  * Delta statistics per octave of N (N < 2^LEMMA): sum Delta^2, sum H(N) Delta(N);
  * samples of S, C3, C9, C27 and of prefix sums of H_c, H_c^2 at the grid N_j = round(2^(10 + j/64)).
Output: JSON (octave table) + NPZ (grid samples).  Cross-checks: H(2^34) = 26869, S(2^34) = 24953018 (amplitude note).
Usage: python collatz_octave_memory_census_20261015.py LOGMAX LEMMA out.json out.npz
"""
import sys, time, json, math
import numpy as np
import numba as nb

LOGMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 30
LEMMA = int(sys.argv[2]) if len(sys.argv) > 2 else min(LOGMAX - 3, 31)
OUTJ = sys.argv[3] if len(sys.argv) > 3 else 'octave_memory.json'
OUTN = sys.argv[4] if len(sys.argv) > 4 else 'octave_memory.npz'
NMAX = 1 << LOGMAX
GUARD = np.uint64((2**64 - 2) * 2 // 3 - 1)   # odd x <= GUARD keeps x + (x+1)//2 below 2^64


@nb.njit(cache=True)
def getbit(P, n):
    return (P[n >> 6] >> np.uint64(n & 63)) & np.uint64(1)


@nb.njit(parallel=True, cache=True)
def parity_block(P, lo, hi, guard, bad):
    # fill bits for n in [lo, hi), lo, hi multiples of 64, using bits below lo
    w0 = lo >> 6
    w1 = hi >> 6
    for w in nb.prange(w0, w1):
        word = np.uint64(0)
        for b in range(64):
            n = (w << 6) + b
            x = np.uint64(n)
            s = 0
            while x >= np.uint64(lo):
                if x & np.uint64(1):
                    if x > guard:
                        bad[0] = 1
                        break
                    x = x + (x + np.uint64(1)) // np.uint64(2)
                else:
                    x = x >> np.uint64(1)
                s += 1
            par = (s + int(getbit(P, int(x)))) & 1
            if par:
                word |= np.uint64(1) << np.uint64(b)
        P[w] = word


def parity_bits(logmax):
    P = np.zeros((1 << logmax) // 64 + 1, dtype=np.uint64)
    # small n directly: tau(1) = 0, tau(2) = 1, ...; compute n < 64 by brute force
    word = 0
    for n in range(1, 64):
        x, s = n, 0
        while x != 1:
            x = x // 2 if x % 2 == 0 else (3 * x + 1) // 2
            s += 1
        if s & 1:
            word |= 1 << n
    P[0] = np.uint64(word)
    bad = np.zeros(1, dtype=np.int64)
    lo = 64
    while lo < (1 << logmax):
        hi = min(2 * lo, 1 << logmax)
        parity_block(P, lo, hi, GUARD, bad)
        if bad[0]:
            raise SystemExit(f'overflow guard hit in block [{lo},{hi}) - aborting')
        lo = hi
    return P


@nb.njit(cache=True)
def stats_pass(P, logmax, lemma, grid, Hstore):
    NMAX = 1 << logmax
    K = logmax + 1
    S = 0
    C3 = np.zeros(3, dtype=np.int64)
    C9 = np.zeros(9, dtype=np.int64)
    C27 = np.zeros(27, dtype=np.int64)
    # per-octave sums for 12 channels: 0..2 = H_(1,r), 3..11 = H_(2,r), r = 0..8
    osum = np.zeros((K, 12), dtype=np.int64)
    osq = np.zeros((K, 12), dtype=np.float64)
    innov = np.zeros((K, K + 1), dtype=np.float64)
    # Delta stats
    d2 = np.zeros(K, dtype=np.float64)
    hd = np.zeros(K, dtype=np.float64)
    # prefix sums of channels and squares
    psum = np.zeros(12, dtype=np.int64)
    psq = np.zeros(12, dtype=np.float64)
    G = grid.shape[0]
    gS = np.zeros(G, dtype=np.int64)
    gC3 = np.zeros((G, 3), dtype=np.int64)
    gC9 = np.zeros((G, 9), dtype=np.int64)
    gC27 = np.zeros((G, 27), dtype=np.int64)
    gps = np.zeros((G, 12), dtype=np.int64)
    gpq = np.zeros((G, 12), dtype=np.float64)
    gi = 0
    stack = np.zeros(K + 2, dtype=np.int64)
    r3 = 1
    r9 = 1
    r27 = 1
    H = np.zeros(12, dtype=np.int64)
    k = 0
    octstart = 1
    for N in range(1, NMAX + 1):
        t = 1 - 2 * int(getbit(P, N))
        S += t
        C3[r3] += t
        C9[r9] += t
        C27[r27] += t
        for r in range(3):
            H[r] = 3 * C3[r] - S
        for r in range(9):
            H[3 + r] = 3 * C9[r] - C3[r % 3]
        # octave bookkeeping
        if N >= 2 * octstart:
            octstart *= 2
            k += 1
        i = N - octstart
        for c in range(12):
            osum[k, c] += H[c]
            osq[k, c] += float(H[c]) * float(H[c])
            psum[c] += H[c]
            psq[c] += float(H[c]) * float(H[c])
        # streaming Haar on H = H[2] within octave k
        carry = H[2]
        lev = 0
        ii = i
        while ii & 1:
            sl = stack[lev]
            diff = float(sl - carry)
            innov[k, lev + 1] += diff * diff / (2.0 * float(1 << lev))
            carry = sl + carry
            lev += 1
            ii >>= 1
        stack[lev] = carry
        # Delta = H(4N') - H(N') at N = 4N'
        if N < (1 << lemma):
            Hstore[N] = H[2]
        if (N & 3) == 0:
            Np = N >> 2
            if Np >= 1 and Np < (1 << lemma):
                hN = Hstore[Np]
                dd = float(H[2] - hN)
                kk = 0
                tmp = Np
                while tmp > 1:
                    tmp >>= 1
                    kk += 1
                d2[kk] += dd * dd
                hd[kk] += float(hN) * dd
        while gi < G and grid[gi] == N:
            gS[gi] = S
            for a in range(3):
                gC3[gi, a] = C3[a]
            for a in range(9):
                gC9[gi, a] = C9[a]
            for a in range(27):
                gC27[gi, a] = C27[a]
            for c in range(12):
                gps[gi, c] = psum[c]
                gpq[gi, c] = psq[c]
            gi += 1
        r3 += 1
        if r3 == 3:
            r3 = 0
        r9 += 1
        if r9 == 9:
            r9 = 0
        r27 += 1
        if r27 == 27:
            r27 = 0
    return S, H[2], osum, osq, innov, d2, hd, gS, gC3, gC9, gC27, gps, gpq


def direct_t(n):
    s = 0
    while n != 1:
        n = n // 2 if n % 2 == 0 else (3 * n + 1) // 2
        s += 1
    return 1 - 2 * (s & 1)


if __name__ == '__main__':
    t0 = time.time()
    P = parity_bits(LOGMAX)
    print(f'parity bits to 2^{LOGMAX}: {time.time() - t0:.1f} s', flush=True)
    # spot check against direct trajectories
    rng = np.random.default_rng(37)
    sample = rng.integers(1, NMAX, size=20000)
    bad = sum(1 for n in sample if direct_t(int(n)) != 1 - 2 * ((int(P[int(n) >> 6]) >> (int(n) & 63)) & 1))
    print(f'spot check of 20000 random n against direct trajectories: {bad} mismatches', flush=True)
    raw = sorted(set(int(round(2 ** (10 + j / 64))) for j in range(0, (LOGMAX - 10) * 64 + 1)))
    grid = np.array([g for g in raw if g <= NMAX], dtype=np.int64)
    Hstore = np.zeros((1 << LEMMA) if LEMMA > 0 else 1, dtype=np.int32)
    t1 = time.time()
    S, Hend, osum, osq, innov, d2, hd, gS, gC3, gC9, gC27, gps, gpq = stats_pass(P, LOGMAX, LEMMA, grid, Hstore)
    print(f'statistics pass: {time.time() - t1:.1f} s; S(2^{LOGMAX}) = {S}, H(2^{LOGMAX}) = {Hend}', flush=True)
    rows = []
    for k in range(LOGMAX):
        X = 1 << k
        mean = [int(v) / X for v in osum[k]]
        en = [float(v) / X ** 2 for v in osq[k]]
        inn = [float(v) / X ** 2 for v in innov[k, :k + 1]]
        tot = float(osq[k, 2])
        meanpart = X * (int(osum[k, 2]) / X) ** 2
        haar_check = (tot - meanpart - float(innov[k].sum())) / max(tot, 1.0)
        row = dict(k=k, X=X, mean_over_sqrtX=[m / math.sqrt(X) for m in mean], energy_over_X2=en,
                   top_mean_part=meanpart / X ** 2, top_innovation_total=float(innov[k].sum()) / X ** 2,
                   top_innovation_by_level=inn, haar_identity_rel_error=haar_check)
        if k < LEMMA:
            row['delta_energy_over_X2'] = float(d2[k]) / X ** 2
            row['H_delta_cross_over_X2'] = float(hd[k]) / X ** 2
        rows.append(row)
    out = dict(session='opus-2026-10-15-S37', logmax=LOGMAX, lemma=LEMMA, S_end=int(S), H_end=int(Hend),
               channels='0..2 = H_(1,r) r=0,1,2; 3..11 = H_(2,r) r=0..8', octaves=rows,
               spot_check_mismatches=int(bad))
    json.dump(out, open(OUTJ, 'w'), indent=1)
    np.savez_compressed(OUTN, grid=grid, S=gS, C3=gC3, C9=gC9, C27=gC27, prefix_sum=gps, prefix_sq=gpq)
    print(f'wrote {OUTJ}, {OUTN}; total {time.time() - t0:.1f} s', flush=True)
