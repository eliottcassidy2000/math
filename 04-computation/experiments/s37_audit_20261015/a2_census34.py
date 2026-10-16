#!/usr/bin/env python3
"""S37 audit, part 2: independent census of t(n) = (-1)^tau(n) to 2^LOGMAX (default 34).

Algorithm (deliberately different from collatz_octave_memory_census_20261015.py):
  * parity bits by dyadic blocks [lo, 2lo): even n read bit(n/2); odd n are advanced by an 8-step jump table
    T^8(256a + r) = 3^c(r) a + T^8(r) until x < lo, using tau(n) = s + tau(T^s n) mod 2 (valid for every s because the
    terminal cycle 1 -> 2 -> 1 has even length); uint64 with an explicit overflow guard (abort, never certify);
  * per-64-bit-word popcount prefix tables of S and of the class 2 mod 3, giving H(n) = 3 C_2(n) - S(n) at any n in O(1);
  * exact int64 chunk sums (summed as Python ints) for octave energies E_X, octave sums of H, Delta statistics
    (Delta(N) = H(4N) - H(N), X <= N < 2X), and class sums mod 243 / prefix sums of H, H^2 at the session's grid
    round(2^(10 + j/64)).  Haar refinement totals are obtained exactly as E_X - (sum H)^2/X (no streaming Haar).
Usage: python a2_census34.py LOGMAX out.json out.npz [bitsfile] [all|bits|stats]
  (for LOGMAX = 34 run 'bits' then 'stats' so each process stays well under 3 minutes)
"""
import sys, time, json, random
import numpy as np
import numba as nb

LOGMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 34
OUTJ = sys.argv[2] if len(sys.argv) > 2 else f'a2_census{LOGMAX}.json'
OUTN = sys.argv[3] if len(sys.argv) > 3 else f'a2_census{LOGMAX}.npz'
NMAX = 1 << LOGMAX
KJ = 8
DELTA_MAX_K = min(30, LOGMAX - 3)
t00 = time.time()


def log(msg):
    print(f'[{time.time() - t00:7.1f} s] {msg}', flush=True)


# jump table
MJ = 1 << KJ
p3 = np.zeros(MJ, np.uint64)
dj = np.zeros(MJ, np.uint64)
for r in range(MJ):
    x, c = r, 0
    for _ in range(KJ):
        if x & 1:
            x = (3 * x + 1) // 2
            c += 1
        else:
            x //= 2
    p3[r] = 3 ** c
    dj[r] = x


@nb.njit(cache=False)
def popc(x):
    x = x - ((x >> np.uint64(1)) & np.uint64(0x5555555555555555))
    x = (x & np.uint64(0x3333333333333333)) + ((x >> np.uint64(2)) & np.uint64(0x3333333333333333))
    x = (x + (x >> np.uint64(4))) & np.uint64(0x0F0F0F0F0F0F0F0F)
    return np.int64((x * np.uint64(0x0101010101010101)) >> np.uint64(56))


@nb.njit(cache=False)
def bit(P, n):
    return np.int64((P[n >> 6] >> np.uint64(n & 63)) & np.uint64(1))


@nb.njit(parallel=True, cache=False)
def fill(P, lo, hi, p3, dj, K, bad):
    w0 = lo >> 6
    w1 = hi >> 6
    lou = np.uint64(lo)
    mask = np.uint64((1 << K) - 1)
    Ku = np.uint64(K)
    UM = np.uint64(0xFFFFFFFFFFFFFFFF)
    for w in nb.prange(w0, w1):
        word = np.uint64(0)
        base = w << 6
        for b in range(64):
            n = base + b
            if (n & 1) == 0:
                par = (1 + bit(P, n >> 1)) & 1
            else:
                x = np.uint64(n)
                s = 0
                while x >= lou:
                    r = np.int64(x & mask)
                    a = x >> Ku
                    pp = p3[r]
                    dd = dj[r]
                    if a > (UM - dd) // pp:
                        bad[0] = 1
                        break
                    x = a * pp + dd
                    s += K
                    if s > 4000000:
                        bad[1] = 1
                        break
                par = (s + bit(P, np.int64(x))) & 1
            if par:
                word |= np.uint64(1) << np.uint64(b)
        P[w] = word


@nb.njit(parallel=True, cache=False)
def word_sums(P, nw, M3):
    tw = np.zeros(nw, np.int64)
    cw = np.zeros(nw, np.int64)
    for w in nb.prange(nw):
        x = P[w]
        m = M3[w % 3]
        tw[w] = 64 - 2 * popc(x)
        cw[w] = popc(m) - 2 * popc(x & m)
    return tw, cw


@nb.njit(cache=False)
def H_at(P, Sw, Cw, M3, n):
    w = n >> 6
    b = n & 63
    pm = (np.uint64(2) << np.uint64(b)) - np.uint64(1)
    x = P[w]
    m = M3[w % 3]
    s = Sw[w] + popc(pm) - 2 * popc(x & pm)
    if w == 0:
        s -= 1
    c = Cw[w] + popc(pm & m) - 2 * popc(x & pm & m)
    return 3 * c - s, s


@nb.njit(cache=False)
def wt(P, n):
    t = 1 - 2 * bit(P, n)
    if n % 3 == 2:
        return 2 * t
    return -t


@nb.njit(parallel=True, cache=False)
def octave_sums(P, Sw, Cw, w0, w1, nchunk):
    out_s = np.zeros(nchunk, np.int64)
    out_q = np.zeros(nchunk, np.int64)
    span = w1 - w0
    for ci in nb.prange(nchunk):
        a = w0 + span * ci // nchunk
        e = w0 + span * (ci + 1) // nchunk
        if a >= e:
            continue
        Hc = 3 * Cw[a] - Sw[a]
        sh = 0
        sq = 0
        for w in range(a, e):
            x = P[w]
            r = w % 3
            for b in range(64):
                t = 1 - 2 * np.int64((x >> np.uint64(b)) & np.uint64(1))
                if (r + b) % 3 == 2:
                    Hc += 2 * t
                else:
                    Hc -= t
                sh += Hc
                sq += Hc * Hc
        out_s[ci] = sh
        out_q[ci] = sq
    return out_s, out_q


@nb.njit(parallel=True, cache=False)
def delta_sums(P, Sw, Cw, M3, X, nchunk):
    d2 = np.zeros(nchunk, np.int64)
    hd = np.zeros(nchunk, np.int64)
    for ci in nb.prange(nchunk):
        a = X + X * ci // nchunk
        e = X + X * (ci + 1) // nchunk
        if a >= e:
            continue
        h1, s1 = H_at(P, Sw, Cw, M3, a)
        h4, s4 = H_at(P, Sw, Cw, M3, 4 * a)
        s2 = 0
        sx = 0
        for N in range(a, e):
            if N > a:
                h1 += wt(P, N)
                for i in range(4 * N - 3, 4 * N + 1):
                    h4 += wt(P, i)
            dd = h4 - h1
            s2 += dd * dd
            sx += h1 * dd
        d2[ci] = s2
        hd[ci] = sx
    return d2, hd


@nb.njit(parallel=True, cache=False)
def grid_pieces(P, Sw, Cw, M3, pst, pen, M):
    npc = pst.shape[0]
    cnt = np.zeros((npc, M), np.int64)
    sh = np.zeros(npc, np.int64)
    sq = np.zeros(npc, np.int64)
    for k in nb.prange(npc):
        a = pst[k]
        e = pen[k]
        Hc, sdummy = H_at(P, Sw, Cw, M3, a - 1)
        r = a % M
        r3 = a % 3
        s1 = 0
        s2 = 0
        for n in range(a, e + 1):
            t = 1 - 2 * bit(P, n)
            cnt[k, r] += t
            if r3 == 2:
                Hc += 2 * t
            else:
                Hc -= t
            s1 += Hc
            s2 += Hc * Hc
            r += 1
            if r == M:
                r = 0
            r3 += 1
            if r3 == 3:
                r3 = 0
        sh[k] = s1
        sq[k] = s2
    return cnt, sh, sq


def tau_parity_direct(n):
    s = 0
    while n != 1:
        n = n // 2 if n % 2 == 0 else (3 * n + 1) // 2
        s += 1
    return s & 1


BITS = sys.argv[4] if len(sys.argv) > 4 else None
STAGE = sys.argv[5] if len(sys.argv) > 5 else 'all'

if __name__ == '__main__':
    nw = NMAX // 64 + 1
    if STAGE in ('all', 'bits'):
        P = np.zeros(nw, np.uint64)
        w0word = 0
        for n in range(1, 64):
            if tau_parity_direct(n):
                w0word |= 1 << n
        P[0] = np.uint64(w0word)
        bad = np.zeros(2, np.int64)
        lo = 64
        while lo < NMAX:
            hi = min(2 * lo, NMAX)
            fill(P, lo, hi, p3, dj, KJ, bad)
            if bad[0] or bad[1]:
                raise SystemExit(f'guard hit in block [{lo},{hi}): overflow={bad[0]} cap={bad[1]} - aborting')
            lo = hi
        if LOGMAX & 1:
            P[NMAX >> 6] |= np.uint64(1)          # tau(2^k) = k
        log(f'parity bits to 2^{LOGMAX} (jump table, K = {KJ})')
    else:
        P = np.load(BITS)
        assert P.shape[0] == nw
        log(f'loaded parity bits from {BITS}')
    rng = random.Random(3737 if STAGE != 'stats' else 7373)
    samp = [rng.randint(1, NMAX) for _ in range(1500)] + list(range(1, 3000)) + [NMAX - i for i in range(200)]
    mism = sum(1 for n in samp if tau_parity_direct(n) != ((int(P[n >> 6]) >> (n & 63)) & 1))
    log(f'validation of {len(samp)} n against big-int trajectories: {mism} mismatches')
    if mism:
        raise SystemExit('validation failed')
    if STAGE == 'bits':
        np.save(BITS, P)
        log(f'saved bits to {BITS}')
        raise SystemExit(0)
    M3 = np.zeros(3, np.uint64)
    for j in range(3):
        m = 0
        for b in range(64):
            if (j + b) % 3 == 2:
                m |= 1 << b
        M3[j] = np.uint64(m)
    tw, cw = word_sums(P, nw, M3)
    tw[0] -= 1                       # exclude n = 0 (bit 0 clear, counted as +1)
    Sw = np.zeros(nw + 1, np.int64)
    Cw = np.zeros(nw + 1, np.int64)
    np.cumsum(tw, out=Sw[1:])
    np.cumsum(cw, out=Cw[1:])
    del tw, cw
    log('prefix tables')
    # values at powers of two
    pow2 = {}
    for k in range(0, LOGMAX + 1):
        h, s = H_at(P, Sw, Cw, M3, 1 << k)
        pow2[k] = (int(s), int(h))
    log(f'S(2^{LOGMAX}) = {pow2[LOGMAX][0]}, H(2^{LOGMAX}) = {pow2[LOGMAX][1]}')
    # direct small-n H for octave sums with X < 64
    Hs = [0]
    for n in range(1, 129):
        Hs.append(Hs[-1] + (2 if n % 3 == 2 else -1) * (1 - 2 * ((int(P[n >> 6]) >> (n & 63)) & 1)))
    octs = []
    for k in range(0, LOGMAX):
        X = 1 << k
        if X < 64:
            sh = sum(Hs[n] for n in range(X, 2 * X))
            sq = sum(Hs[n] ** 2 for n in range(X, 2 * X))
        else:
            nchunk = max(1, min(4096, (X >> 6) // 256))
            a1, a2 = octave_sums(P, Sw, Cw, X >> 6, (2 * X) >> 6, nchunk)
            sh = sum(int(v) for v in a1)
            sq = sum(int(v) for v in a2)
        octs.append(dict(k=k, X=X, sumH=sh, sumH2=str(sq), e=sq / X ** 2, mean_part=sh * sh / X / X ** 2,
                         refinement=(sq - sh * sh / X) / X ** 2))
    log('octave energies')
    for k in range(0, DELTA_MAX_K + 1):
        X = 1 << k
        nchunk = max(1, min(4096, X // 4096))
        d2, hd = delta_sums(P, Sw, Cw, M3, X, nchunk)
        D = sum(int(v) for v in d2)
        C = sum(int(v) for v in hd)
        octs[k].update(delta=D / X ** 2, cross=C / X ** 2, D_exact=str(D), C_exact=str(C))
    log(f'Delta statistics to k = {DELTA_MAX_K}')
    raw = sorted(set(int(round(2 ** (10 + j / 64))) for j in range(0, (LOGMAX - 10) * 64 + 1)))
    grid = np.array([g for g in raw if g <= NMAX], dtype=np.int64)
    starts, ends, owner = [], [], []
    prev = 0
    PIECE = 1 << 22
    for i, g in enumerate(grid):
        a = prev + 1
        while a <= g:
            e = min(g, a + PIECE - 1)
            starts.append(a)
            ends.append(e)
            owner.append(i)
            a = e + 1
        prev = int(g)
    pst = np.array(starts, np.int64)
    pen = np.array(ends, np.int64)
    own = np.array(owner, np.int64)
    cnt, sh, sq = grid_pieces(P, Sw, Cw, M3, pst, pen, 243)
    G = len(grid)
    icnt = np.zeros((G, 243), np.int64)
    ish = np.zeros(G, np.int64)
    isq = np.zeros(G, np.int64)
    np.add.at(icnt, own, cnt)
    np.add.at(ish, own, sh)
    np.add.at(isq, own, sq)
    C243 = np.cumsum(icnt, axis=0)
    PS = np.cumsum(ish)
    PQ = np.cumsum(isq.astype(np.float64))
    log(f'grid samples ({G} points, {len(pst)} pieces)')
    # consistency: S at grid points from class counts equals table values
    chk = 0
    for i in range(0, G, 37):
        h, s = H_at(P, Sw, Cw, M3, int(grid[i]))
        Cl3 = np.array([C243[i, r::3].sum() for r in range(3)])
        chk += (C243[i].sum() != s) or (3 * Cl3[2] - s != h)
    log(f'grid/table consistency failures: {chk}')
    out = dict(logmax=LOGMAX, S_end=pow2[LOGMAX][0], H_end=pow2[LOGMAX][1], pow2={str(k): v for k, v in pow2.items()},
               octaves=octs, validation_mismatches=mism, grid_table_check_failures=int(chk))
    json.dump(out, open(OUTJ, 'w'), indent=1)
    np.savez_compressed(OUTN, grid=grid, C243=C243, prefix_sum=PS, prefix_sq=PQ, interval_sq=isq, interval_sum=ish)
    log(f'wrote {OUTJ}, {OUTN}')
