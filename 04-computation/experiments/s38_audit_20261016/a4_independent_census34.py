#!/usr/bin/env python3
"""S38 audit, item 3: an independent exact Collatz parity census to 2^34 (own code, single thread) checked against the
dense S38 census (dense40.npz / dense40.json / dense40_moments.npz), plus independent checks above 2^34.

T(n) = n/2 or (3n+1)/2, tau(n) = number of T-steps to 1, t(n) = (-1)^tau(n), S(N) = sum_(n<=N) t(n),
H(N) = sum_(n<=N) w_n t(n), w_n = 3[n = 2 mod 3] - 1.
(1) Census: parity bits tau(n) mod 2 for all n <= 2^34 by stopping-time descent (T-steps until the value drops below n,
    then the stored bit; even n by halving).  Plain uint64 single steps (no jump table, no block structure shared with
    the S38 code); the largest T-value reached below 2^34 is about 1.04e19 < 2^64 and the script aborts if an odd value
    above (2^65 - 4)/3 ever occurs.  Recorded at every dense grid point N_j = round(2^(10 + j/1024)) <= 2^34: S, the class
    sums mod 3, H, and the exact prefix sums of H and H^2; per block of 2^17: sum t, sum (n - lo) t, sum H, sum H^2.
(2) Comparison with dense40.npz at every grid point <= 2^34, with the octave sums k = 20..33 of dense40.json, and with
    the first 2^17 rows of dense40_moments.npz.
(3) 96 random blocks of 2^17 in [2^34, 2^40) (16 per octave): all t(n) by descent below 2^34 with the table of (1)
    (exact Python big integers when a value would leave uint64), moments compared with dense40_moments.npz.
(4) 20000 random n in [2^40, 2^42) and the path-record starting values in [2^37, 2^42): t(n) from the full Python
    big-integer trajectory to 1 against the descent of (3).
Usage: python a4_independent_census34.py SCRATCH_S38_DIR
"""
import sys, os, json, time, math, random
import numpy as np
import numba as nb

SD = sys.argv[1]
LOG = 34
NMAX = 1 << LOG
SUB = 1 << 17
GUARD = (2 ** 64 - 2) * 2 // 3 - 1          # odd x <= GUARD: x + (x+1)//2 < 2^64
HERE = os.path.dirname(os.path.abspath(__file__))


@nb.njit
def census(P, nmax, grid, gS, gC, gH, gp1, gp2, bS1, bS2, bm0, bm1, guard):
    one = np.uint64(1)
    G = grid.shape[0]
    gi = 0
    S = 0
    C0 = 0
    C1 = 0
    C2 = 0
    H = 0
    s1 = 0
    s2 = 0
    m0 = 0
    m1 = 0
    r3 = 1
    lo = 0
    for n in range(1, nmax + 1):
        if n == 1:
            par = 0
        elif (n & 1) == 0:
            h = n >> 1
            par = 1 - int((P[h >> 6] >> np.uint64(h & 63)) & one)
        else:
            x = np.uint64(n)
            un = np.uint64(n)
            s = 0
            while x >= un:
                if x & one:
                    if x > guard:
                        return n
                    x = x + (x + one) // np.uint64(2)
                else:
                    x = x >> one
                s += 1
            xi = int(x)
            par = (s + int((P[xi >> 6] >> np.uint64(xi & 63)) & one)) & 1
        if par:
            P[n >> 6] |= one << np.uint64(n & 63)
        t = 1 - 2 * par
        S += t
        if r3 == 0:
            C0 += t
            H -= t
        elif r3 == 1:
            C1 += t
            H -= t
        else:
            C2 += t
            H += 2 * t
        s1 += H
        s2 += H * H
        m0 += t
        m1 += (n - lo) * t
        while gi < G and grid[gi] == n:
            gS[gi] = S
            gC[gi, 0] = C0
            gC[gi, 1] = C1
            gC[gi, 2] = C2
            gH[gi] = H
            gp1[gi] = s1
            gp2[gi] = s2
            gi += 1
        if (n & (SUB - 1)) == SUB - 1 or n == nmax:
            b = n >> 17
            bS1[b] = s1
            bS2[b] = s2
            bm0[b] = m0
            bm1[b] = m1
            s1 = 0
            s2 = 0
            m0 = 0
            m1 = 0
            lo = (b + 1) << 17
        r3 += 1
        if r3 == 3:
            r3 = 0
    return 0


@nb.njit
def descend_block(P, lo, B, low, guard, flags):
    # t(n) for n in [lo, lo + B) by descent below `low` (table P); flags[i] = 1 if uint64 would overflow
    one = np.uint64(1)
    out = np.zeros(B, dtype=np.int64)
    ul = np.uint64(low)
    for i in range(B):
        n = lo + i
        x = np.uint64(n)
        s = 0
        bad = False
        while x >= ul:
            if x & one:
                if x > guard:
                    bad = True
                    break
                x = x + (x + one) // np.uint64(2)
            else:
                x = x >> one
            s += 1
        if bad:
            flags[i] = 1
            continue
        xi = int(x)
        par = (s + int((P[xi >> 6] >> np.uint64(xi & 63)) & one)) & 1
        out[i] = 1 - 2 * par
    return out


def tau_parity_big(n):
    s = 0
    while n != 1:
        n = n // 2 if n % 2 == 0 else (3 * n + 1) // 2
        s += 1
    return s & 1


def t_by_descent_big(P, n, low):
    s = 0
    x = n
    while x >= low:
        x = x // 2 if x % 2 == 0 else (3 * x + 1) // 2
        s += 1
    return 1 - 2 * ((s + ((int(P[x >> 6]) >> (x & 63)) & 1)) & 1)


if __name__ == '__main__':
    t0 = time.time()
    raw = sorted(set(int(round(2 ** (10 + j / 1024))) for j in range(0, 30 * 1024 + 1)))
    grid = np.array([g for g in raw if g <= NMAX], dtype=np.int64)
    G = len(grid)
    P = np.zeros((NMAX >> 6) + 1, dtype=np.uint64)
    gS = np.zeros(G, np.int64); gC = np.zeros((G, 3), np.int64); gH = np.zeros(G, np.int64)
    gp1 = np.zeros(G, np.int64); gp2 = np.zeros(G, np.int64)
    nb_ = (NMAX >> 17) + 1
    bS1 = np.zeros(nb_, np.int64); bS2 = np.zeros(nb_, np.int64); bm0 = np.zeros(nb_, np.int64); bm1 = np.zeros(nb_, np.int64)
    # warm-up compile on a small range
    census(np.zeros(1 << 10, np.uint64), 1 << 15, grid[:5].copy(), np.zeros(5, np.int64), np.zeros((5, 3), np.int64),
           np.zeros(5, np.int64), np.zeros(5, np.int64), np.zeros(5, np.int64), np.zeros(2, np.int64), np.zeros(2, np.int64),
           np.zeros(2, np.int64), np.zeros(2, np.int64), np.uint64(GUARD))
    t1 = time.time()
    rc = census(P, NMAX, grid, gS, gC, gH, gp1, gp2, bS1, bS2, bm0, bm1, np.uint64(GUARD))
    if rc:
        raise SystemExit(f'uint64 guard hit at n = {rc}')
    print(f'census to 2^{LOG} (inclusive): {time.time() - t1:.1f} s; S(2^34) = {int(gS[-1])}, H(2^34) = {int(gH[-1])}', flush=True)
    res = dict(S_2p34=int(gS[-1]), H_2p34=int(gH[-1]), n_grid=G, census_seconds=time.time() - t1)
    # exact prefix sums at grid points
    cumS1 = [0] * (nb_ + 1)
    cumS2 = [0] * (nb_ + 1)
    for b in range(nb_):
        cumS1[b + 1] = cumS1[b] + int(bS1[b])
        cumS2[b + 1] = cumS2[b] + int(bS2[b])
    pre1 = [cumS1[int(g) >> 17] + int(gp1[i]) for i, g in enumerate(grid)]
    pre2 = [cumS2[int(g) >> 17] + int(gp2[i]) for i, g in enumerate(grid)]
    # (2) comparison with the dense census
    d = np.load(os.path.join(SD, 'dense40.npz'))
    dg = np.array(d['grid'])
    sel = dg <= NMAX
    idx = np.searchsorted(dg, grid)
    res['grid_identical'] = bool(np.array_equal(dg[sel], grid))
    dS = np.array(d['S'])[idx]; dC = np.array(d['C3'])[idx]; dH = np.array(d['H'])[idx]
    dp1 = [int(x) for x in np.array(d['preH1_str'])[idx]]
    dp2 = [int(x) for x in np.array(d['preH2_str'])[idx]]
    res['S_mismatches'] = int((dS != gS).sum())
    res['C3_mismatches'] = int((dC != gC).any(axis=1).sum())
    res['H_mismatches'] = int((dH != gH).sum())
    res['preH1_mismatches'] = sum(1 for a, b in zip(dp1, pre1) if a != b)
    res['preH2_mismatches'] = sum(1 for a, b in zip(dp2, pre2) if a != b)
    dj = json.load(open(os.path.join(SD, 'dense40.json')))
    om = []
    for k in range(20, LOG):
        b0, b1 = (1 << k) >> 17, (1 << (k + 1)) >> 17
        mine = (str(cumS1[b1] - cumS1[b0]), str(cumS2[b1] - cumS2[b0]))
        if tuple(dj['octave_sums_H_H2'][str(k)]) != mine:
            om.append(k)
    res['octave_sum_mismatches_k20_33'] = om
    mm = np.load(os.path.join(SD, 'dense40_moments.npz'))
    mom = np.array(mm['moments'])
    nrow = NMAX >> 17
    res['moment_row_mismatches_below_2p34'] = int(((mom[:nrow, 0] != bm0[:nrow]) | (mom[:nrow, 1] != bm1[:nrow])).sum())
    print({k: v for k, v in res.items()}, flush=True)
    # (3) random blocks above 2^34
    rng = random.Random(3838)
    blocks = []
    for k in range(34, 40):
        for _ in range(16):
            blocks.append(rng.randrange((1 << k) >> 17, (1 << (k + 1)) >> 17))
    descend_block(P, NMAX, 64, NMAX, np.uint64(GUARD), np.zeros(64, np.uint8))
    t2 = time.time()
    nbad_blocks = 0
    nflag = 0
    for b in blocks:
        lo = b << 17
        flags = np.zeros(SUB, np.uint8)
        tv = descend_block(P, lo, SUB, NMAX, np.uint64(GUARD), flags)
        for i in np.nonzero(flags)[0]:
            tv[i] = t_by_descent_big(P, lo + int(i), NMAX)
            nflag += 1
        m0 = int(tv.sum())
        m1 = int((np.arange(SUB, dtype=np.int64) * tv).sum())
        if m0 != int(mom[b, 0]) or m1 != int(mom[b, 1]):
            nbad_blocks += 1
            print('block mismatch', b, m0, int(mom[b, 0]), m1, int(mom[b, 1]), flush=True)
    res['random_blocks_above_2p34'] = len(blocks)
    res['random_block_moment_mismatches'] = nbad_blocks
    res['random_blocks_values_beyond_uint64'] = nflag
    print(f'(3) {len(blocks)} random blocks of 2^17 in [2^34, 2^40): {nbad_blocks} moment mismatches '
          f'({nflag} n needed big integers; {time.time() - t2:.1f} s)', flush=True)
    # (4) random n in [2^40, 2^42) and path records
    t3 = time.time()
    records = [204430613247, 231913730799, 272025660543, 446559217279, 567839862631, 871673828443, 2674309547647,
               3716509988199]
    sample = [rng.randrange(1 << 40, 1 << 42) for _ in range(20000)] + records
    bad = 0
    big = 0
    for n in sample:
        a = 1 - 2 * tau_parity_big(n)
        fl = np.zeros(1, np.uint8)
        tv = descend_block(P, n, 1, NMAX, np.uint64(GUARD), fl)
        if fl[0]:
            big += 1
            b_ = t_by_descent_big(P, n, NMAX)
        else:
            b_ = int(tv[0])
        if a != b_:
            bad += 1
    res['random_n_2p40_2p42_plus_records'] = len(sample)
    res['random_n_mismatches'] = bad
    res['random_n_needing_big_integers'] = big
    # trajectory successors of the two path records above 2^41 that stay >= 2^40: all checked
    succ_bad = 0
    succ_n = 0
    for r in (2674309547647, 3716509988199):
        x = r
        while x >= (1 << 40) and succ_n < 4000:
            ref = 1 - 2 * tau_parity_big(x)
            got = t_by_descent_big(P, x, NMAX)
            if got != ref:
                succ_bad += 1
            succ_n += 1
            x = x // 2 if x % 2 == 0 else (3 * x + 1) // 2
    res['record_trajectory_points_checked'] = succ_n
    res['record_trajectory_mismatches'] = succ_bad
    res['part4_seconds'] = time.time() - t3
    res['total_seconds'] = time.time() - t0
    print(json.dumps(res, indent=1), flush=True)
    json.dump(res, open(os.path.join(HERE, 'a4_independent_census34.json'), 'w'), indent=1)
    # keep the per-grid census values for later scripts (small)
    np.savez_compressed(os.path.join(HERE, 'a4_census34_grid.npz'), grid=grid, S=gS, C3=gC, H=gH,
                        preH1_str=np.array([str(x) for x in pre1]), preH2_str=np.array([str(x) for x in pre2]),
                        mom0=bm0, mom1=bm1)
