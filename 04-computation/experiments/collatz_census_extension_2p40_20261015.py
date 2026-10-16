#!/usr/bin/env python3
"""Streaming extension of the exact Collatz parity census to 2^LOGMAX (default 40) with a stored table below 2^LOW.

T(n) = n/2 or (3n+1)/2, tau(n) = number of T-steps from n to 1, t(n) = (-1)^tau(n), H(N) = sum_(n<=N) w_n t(n) with
w_n = 3[n = 2 mod 3] - 1.  The parity bits below 2^LOW (default 37) come from collatz_census_extension_2p37_20261015.py
(uint64 descent with exact big-integer fallback, a 16 GB table at LOW = 37).  For 2^LOW <= n < 2^LOGMAX the parity is
computed on the fly: tau(n) = s + tau(T^s n) (mod 2) for every s (the terminal cycle 1 -> 2 -> 1 has even length), so
the script descends with the exact 8-step jump T^8(256h + l) = 3^(a(l)) h + T^8(l) (8 steps, even, parity unchanged)
until the value is below 2^LOW and looks up the table.  Values above 2^64/26 continue with single steps, and values
that would leave uint64 continue in exact two-limb (128-bit) arithmetic; the script aborts beyond 126 bits.

One heavy parallel pass over blocks of 2^20 records per block: class counts of t modulo 243, the signed moments
(sum t, sum (n - lo) t) of every sub-block of 2^17, and the within-block partial discrepancy hp(n) = sum_(lo<=m<=n)
w_m t(m) through sum hp, sum hp^2 and hp(end).  After exact prefix sums (Python integers) these give the exact block
sums of H and H^2: sum H = B' H0 + sum hp, sum H^2 = B' H0^2 + 2 H0 sum hp + sum hp^2 (H0 = H(lo - 1), B' = block
size).  A light pass re-walks only the blocks containing grid points N_j = round(2^(10 + j/64)) < 2^LOGMAX and records
the class sums modulo 243, H(N_j) and the prefix sums of H, H^2 there.
Outputs: census NPZ (same keys as the 2^37 extension), moments NPZ (N = 2^LOGMAX, B = 2^17), JSON (exact octave sums of
H and H^2, S(N), H(N) inclusive of n = N with tau(2^k) = k).
Checks: every grid point below 2^37 must reproduce the 2^37 extension exactly; S(2^37) = -152262364 and H(2^37) = -67763
(collatz_boundary_moments_2p37_20261015.json); 20000 random n in [2^LOW, 2^LOGMAX) against Python trajectories; the
jump table against single steps.
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
JGUARD = np.uint64(2**64 // 26)


def jump_tables():
    A = np.zeros(256, dtype=np.uint64)
    Cc = np.zeros(256, dtype=np.uint64)
    for l in range(256):
        x, a = l, 0
        for _ in range(8):
            if x & 1:
                x = (3 * x + 1) // 2
                a += 1
            else:
                x //= 2
        A[l] = 3 ** a
        Cc[l] = x
    return A, Cc


@nb.njit
def tval(P, n, low, A, Cc, jguard, guard):
    # t(n) in {+1,-1}; 0 signals a value beyond 126 bits (abort)
    if n < low:
        return 1 - 2 * int((P[n >> 6] >> np.uint64(n & 63)) & np.uint64(1))
    x = np.uint64(n)
    lo64 = np.uint64(low)
    s = 0
    while x >= lo64:
        if x <= jguard:
            l = x & np.uint64(255)
            x = A[l] * (x >> np.uint64(8)) + Cc[l]
            continue
        if x & np.uint64(1):
            if x > guard:
                hi = np.uint64(0)
                lo = x
                while hi != np.uint64(0) or lo >= lo64:
                    if lo & np.uint64(1):
                        if hi >= np.uint64(1) << np.uint64(62):
                            return 0
                        l1 = lo + np.uint64(1)
                        h1 = hi + (np.uint64(1) if l1 == np.uint64(0) else np.uint64(0))
                        yl = (l1 >> np.uint64(1)) | ((h1 & np.uint64(1)) << np.uint64(63))
                        yh = h1 >> np.uint64(1)
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
def heavy_pass(P, nblocks, B, Sb, Mod, low, A, Cc, jguard, guard):
    cnt = np.zeros((nblocks, Mod), dtype=np.int64)
    nsub = B // Sb
    mom = np.zeros((nblocks * nsub, 2), dtype=np.int64)
    hs1 = np.zeros(nblocks, dtype=np.int64)
    hs2 = np.zeros(nblocks, dtype=np.int64)
    hend = np.zeros(nblocks, dtype=np.int64)
    mx = np.zeros(nblocks, dtype=np.int64)
    bad = np.zeros(nblocks, dtype=np.int64)
    for b in nb.prange(nblocks):
        lo = b * B
        n0 = lo if lo > 0 else 1
        r = n0 % Mod
        r3 = n0 % 3
        hp = 0
        s1 = 0
        s2 = 0
        m = 0
        for n in range(n0, lo + B):
            t = tval(P, n, low, A, Cc, jguard, guard)
            if t == 0:
                bad[b] = n
                t = 1
            cnt[b, r] += t
            k = (n - lo) >> 17
            mom[b * nsub + k, 0] += t
            mom[b * nsub + k, 1] += (n - lo - (k << 17)) * t
            if r3 == 2:
                hp += 2 * t
            else:
                hp -= t
            s1 += hp
            s2 += hp * hp
            ah = hp if hp >= 0 else -hp
            if ah > m:
                m = ah
            r += 1
            if r == Mod:
                r = 0
            r3 += 1
            if r3 == 3:
                r3 = 0
        hs1[b] = s1
        hs2[b] = s2
        hend[b] = hp
        mx[b] = m
    return cnt, mom, hs1, hs2, hend, mx, bad


@nb.njit(parallel=True)
def grid_pass(P, gblocks, B, Mod, start, Hstart, grid, glo, ghi, low, A, Cc, jguard, guard):
    G = grid.shape[0]
    snapC = np.zeros((G, Mod), dtype=np.int64)
    snapH = np.zeros(G, dtype=np.int64)
    snapP1 = np.zeros(G, dtype=np.int64)
    snapP2 = np.zeros(G, dtype=np.int64)
    for i in nb.prange(gblocks.shape[0]):
        b = gblocks[i]
        C = start[i].copy()
        lo = b * B
        n0 = lo if lo > 0 else 1
        r = n0 % Mod
        r3 = n0 % 3
        hp = 0
        s1 = 0
        s2 = 0
        gi = glo[i]
        ge = ghi[i]
        n = n0
        while gi < ge:
            t = tval(P, n, low, A, Cc, jguard, guard)
            C[r] += t
            if r3 == 2:
                hp += 2 * t
            else:
                hp -= t
            s1 += hp
            s2 += hp * hp
            while gi < ge and grid[gi] == n:
                for a in range(Mod):
                    snapC[gi, a] = C[a]
                snapH[gi] = Hstart[i] + hp
                snapP1[gi] = s1
                snapP2[gi] = s2
                gi += 1
            n += 1
            r += 1
            if r == Mod:
                r = 0
            r3 += 1
            if r3 == 3:
                r3 = 0
    return snapC, snapH, snapP1, snapP2


if __name__ == '__main__':
    t0 = time.time()
    A, Cc = jump_tables()
    # jump table against single steps
    rng = random.Random(4040)
    for _ in range(2000):
        x = rng.randrange(1, 1 << 50)
        y = x
        for _ in range(8):
            y = y // 2 if y % 2 == 0 else (3 * y + 1) // 2
        assert y == int(A[x & 255]) * (x >> 8) + int(Cc[x & 255])
    P, patched = ext.parity_bits(LOW)
    print(f'parity table to 2^{LOW}: {time.time() - t0:.1f} s; {len(patched)} words recomputed exactly', flush=True)
    low = 1 << LOW
    bad = 0
    for _ in range(20000):
        n = rng.randrange(low, 1 << LOGMAX)
        if tval(P, n, low, A, Cc, JGUARD, GUARD) != 1 - 2 * ext.tau_parity(n):
            bad += 1
    print(f'spot check of 20000 random n in [2^{LOW}, 2^{LOGMAX}) against Python trajectories: {bad} mismatches', flush=True)
    if bad:
        raise SystemExit('parity mismatch')
    nblocks = (1 << LOGMAX) // BLOCK
    t1 = time.time()
    cnt, mom, hs1, hs2, hend, mxb, badb = heavy_pass(P, nblocks, BLOCK, SUB, M, low, A, Cc, JGUARD, GUARD)
    if badb.any():
        raise SystemExit(f'128-bit overflow at n = {int(badb[np.nonzero(badb)[0][0]])}')
    mx = int(mxb.max())
    if BLOCK * mx * mx >= 2 ** 63:
        raise SystemExit('possible int64 overflow in block sums of hp^2')
    print(f'heavy pass: {time.time() - t1:.1f} s; max within-block |hp| = {mx}', flush=True)
    start = np.zeros_like(cnt)
    start[1:] = np.cumsum(cnt[:-1], axis=0)
    # exact block sums of H and H^2 (Python integers)
    H0 = [0] * nblocks
    cum1 = [0] * (nblocks + 1)
    cum2 = [0] * (nblocks + 1)
    h = 0
    for b in range(nblocks):
        H0[b] = h
        size = BLOCK if b > 0 else BLOCK - 1
        sH = size * h + int(hs1[b])
        sH2 = size * h * h + 2 * h * int(hs1[b]) + int(hs2[b])
        cum1[b + 1] = cum1[b] + sH
        cum2[b + 1] = cum2[b] + sH2
        h += int(hend[b])
    raw = sorted(set(int(round(2 ** (10 + j / 64))) for j in range(0, (LOGMAX - 10) * 64 + 1)))
    grid = np.array([g for g in raw if g < (1 << LOGMAX)], dtype=np.int64)
    bidx = grid // BLOCK
    gblocks = np.unique(bidx)
    glo = np.searchsorted(bidx, gblocks, side='left').astype(np.int64)
    ghi = np.searchsorted(bidx, gblocks, side='right').astype(np.int64)
    t2 = time.time()
    snapC, snapH, snapP1, snapP2 = grid_pass(P, gblocks, BLOCK, M, start[gblocks].copy(),
                                             np.array([H0[int(b)] for b in gblocks], dtype=np.int64), grid, glo, ghi,
                                             low, A, Cc, JGUARD, GUARD)
    print(f'grid pass over {len(gblocks)} blocks: {time.time() - t2:.1f} s', flush=True)
    pre1, pre2 = [], []
    for i in range(len(grid)):
        b = int(bidx[i])
        hb = H0[b]
        lo = b * BLOCK
        cntn = int(grid[i]) - max(lo, 1) + 1
        pre1.append(cum1[b] + cntn * hb + int(snapP1[i]))
        pre2.append(cum2[b] + cntn * hb * hb + 2 * hb * int(snapP1[i]) + int(snapP2[i]))
    Nmax = 1 << LOGMAX
    tN = -1 if LOGMAX % 2 else 1
    final = (start[-1] + cnt[-1]).astype(object)
    final[Nmax % M] += tN
    S_N = int(sum(final))
    C2 = int(sum(final[2::3]))
    H_N = 3 * C2 - S_N
    S = snapC.sum(axis=1)
    i37 = int(np.nonzero(grid == (1 << 37))[0][0]) if LOGMAX > 37 else None
    if i37 is not None:
        print(f'S(2^37) = {int(S[i37])} (moments script -152262364), H(2^37) = {int(snapH[i37])} (-67763)', flush=True)
    print(f'S(2^{LOGMAX}) = {S_N}, H(2^{LOGMAX}) = {H_N}; max |H| at grid points {int(np.abs(snapH).max())}', flush=True)
    same = None
    if OLD:
        with np.load(OLD) as Z:
            g0 = np.array(Z['grid'])
            c0 = np.array(Z['classsums'])
            h0 = np.array(Z['H'])
            q0 = [int(x) for x in Z['preH2_str']]
            q1 = [int(x) for x in Z['preH1_str']]
        idx = np.searchsorted(grid, g0)
        same = bool(np.array_equal(grid[idx], g0) and np.array_equal(snapC[idx], c0) and np.array_equal(snapH[idx], h0)
                    and all(pre2[int(i)] == q for i, q in zip(idx, q0)) and all(pre1[int(i)] == q for i, q in zip(idx, q1)))
        print(f'all {len(g0)} grid points of the 2^37 extension reproduced exactly (class sums, H, prefix sums): {same}',
              flush=True)
    octs = {}
    for k in range(20, LOGMAX):
        b0, b1 = (1 << k) // BLOCK, (1 << (k + 1)) // BLOCK
        octs[str(k)] = (str(cum1[b1] - cum1[b0]), str(cum2[b1] - cum2[b0]))
    np.savez_compressed(OUTC, grid=grid, classsums=snapC, H=snapH,
                        preH1=np.array([float(x) for x in pre1]), preH2=np.array([float(x) for x in pre2]),
                        preH1_str=np.array([str(x) for x in pre1]), preH2_str=np.array([str(x) for x in pre2]),
                        block=BLOCK, logmax=LOGMAX)
    np.savez_compressed(OUTM, moments=mom, N=np.int64(Nmax), B=np.int64(SUB))
    json.dump(dict(N=Nmax, logN=LOGMAX, low=LOW, t_N=tN, S=S_N, H=H_N, reproduces_2p37=same,
                   octave_sums_H_H2=octs, seconds=time.time() - t0), open(OUTJ, 'w'), indent=1)
    print(f'wrote {OUTC}, {OUTM}, {OUTJ}; total {time.time() - t0:.1f} s', flush=True)
