#!/usr/bin/env python3
"""S37 audit 2 (independent referee), script b3: an independent exact census of t(n) = (-1)^tau(n).

Algorithm (deliberately different from collatz_census_extension_2p37_20261015.py, which descends below the octave
start in single T-steps and stores a 2^37-bit table):
  * a parity table only for n < 2^TAB (default 2^28, 32 MB), built octave by octave;
  * every n >= 2^TAB is pushed below 2^TAB with 16-step 2-adic jumps T^16(2^16 h + l) = 3^c(l) h + T^16(l)
    (an even number of steps, so t(n) = t(T^16 n) because t(m) = -t(T(m)) for every m >= 1, including the
    terminal cycle 1 -> 2 -> 1), switching to single T-steps above 2^54 and to exact Python integers above 2^63;
  * class sums modulo 243, the H-sums and the B' = 2^17 block moments are accumulated per block of 2^20 as local
    quantities (no prefix needed inside the parallel pass) and assembled exactly with Python integers later
    (b3b_census_compare.py).
Usage: python b3_census_independent.py LOGLO LOGHI OUT.npz [census_ext.npz for the grid]
       processes the blocks covering [2^LOGLO, 2^LOGHI) (LOGLO = 0 means from n = 1).
"""
import sys, os, time, json, random
import numpy as np
import numba as nb

LOGLO = int(sys.argv[1]); LOGHI = int(sys.argv[2]); OUT = sys.argv[3]
GRIDNPZ = sys.argv[4] if len(sys.argv) > 4 else None
TAB = 28
B = 1 << 20
SBLOG = 17
M = 243

# ---- 16-step jump table ----
c16 = np.zeros(1 << 16, dtype=np.uint8)
d16 = np.zeros(1 << 16, dtype=np.uint64)
for l in range(1 << 16):
    y, c = l, 0
    for _ in range(16):
        if y & 1:
            y = (3 * y + 1) // 2
            c += 1
        else:
            y //= 2
    c16[l] = c
    d16[l] = y
p3 = np.array([3 ** k for k in range(17)], dtype=np.uint64)


@nb.njit(inline='always')
def par_of(n, tab, lim, c16, d16, p3):
    x = np.uint64(n)
    p = 0
    J = np.uint64(1) << np.uint64(54)
    G = np.uint64(1) << np.uint64(63)
    one = np.uint64(1)
    while x >= lim:
        if x < J:
            l = x & np.uint64(0xFFFF)
            x = p3[c16[l]] * (x >> np.uint64(16)) + d16[l]
        else:
            if x & one:
                if x >= G:
                    return -1
                x = x + ((x + one) >> one)
            else:
                x = x >> one
            p ^= 1
    xi = np.int64(x)
    return p ^ int((tab[xi >> 6] >> np.uint64(xi & 63)) & one)


@nb.njit(parallel=True)
def table_octave(tab, w0, w1, lo, c16, d16, p3, bad):
    lim = np.uint64(lo)
    for w in nb.prange(w0, w1):
        word = np.uint64(0)
        for b in range(64):
            n = (w << 6) + b
            pb = par_of(n, tab, lim, c16, d16, p3)
            if pb < 0:
                bad[0] = 1
            elif pb == 1:
                word |= np.uint64(1) << np.uint64(b)
        tab[w] = word


def tau_parity_py(n):
    s = 0
    while n != 1:
        n = n // 2 if n % 2 == 0 else (3 * n + 1) // 2
        s += 1
    return s & 1


def build_table():
    tab = np.zeros((1 << TAB) // 64, dtype=np.uint64)
    word = 0
    for n in range(1, 64):
        if tau_parity_py(n):
            word |= 1 << n
    tab[0] = np.uint64(word)
    lo = 64
    bad = np.zeros(1, dtype=np.int64)
    while lo < (1 << TAB):
        table_octave(tab, lo >> 6, (2 * lo) >> 6, lo, c16, d16, p3, bad)
        lo *= 2
    if bad[0]:
        raise SystemExit('overflow while building the table')
    return tab


@nb.njit(parallel=True)
def block_pass(tab, lim, c16, d16, p3, blist, grid, glo, ghi, ovn, ovt,
               LOC, HEND, SH1, SH2, MOM, GLOC, GH, GSH1, GSH2, FL, FLC):
    nb_ = blist.shape[0]
    for bi in nb.prange(nb_):
        b = blist[bi]
        lo = b << 20
        loc = np.zeros(243, dtype=np.int64)
        S = 0
        C2 = 0
        s1 = 0
        s2 = 0
        m0 = np.zeros(8, dtype=np.int64)
        m1 = np.zeros(8, dtype=np.int64)
        gi = glo[b]
        ge = ghi[b]
        n0 = lo if lo > 0 else 1
        r = n0 % 243
        for n in range(n0, lo + (1 << 20)):
            pb = par_of(n, tab, lim, c16, d16, p3)
            if pb < 0:
                k = np.searchsorted(ovn, n)
                if k < ovn.shape[0] and ovn[k] == n:
                    t = int(ovt[k])
                else:
                    c = FLC[bi]
                    if c < FL.shape[1]:
                        FL[bi, c] = n
                    FLC[bi] = c + 1
                    t = 0
            else:
                t = 1 - 2 * pb
            loc[r] += t
            S += t
            if r % 3 == 2:
                C2 += t
            h = 3 * C2 - S
            s1 += h
            s2 += h * h
            off = n - lo
            sb = off >> 17
            m0[sb] += t
            m1[sb] += (off - (sb << 17)) * t
            while gi < ge and grid[gi] == n:
                for a in range(243):
                    GLOC[gi, a] = loc[a]
                GH[gi] = h
                GSH1[gi] = s1
                GSH2[gi] = s2
                gi += 1
            r += 1
            if r == 243:
                r = 0
        for a in range(243):
            LOC[bi, a] = loc[a]
        HEND[bi] = 3 * C2 - S
        SH1[bi] = s1
        SH2[bi] = s2
        for k in range(8):
            MOM[bi, k, 0] = m0[k]
            MOM[bi, k, 1] = m1[k]


def run_blocks(tab, blist, grid, glo, ghi, ovn, ovt, GLOC, GH, GSH1, GSH2):
    nbk = len(blist)
    LOC = np.zeros((nbk, M), dtype=np.int32)
    HEND = np.zeros(nbk, dtype=np.int64)
    SH1 = np.zeros(nbk, dtype=np.int64)
    SH2 = np.zeros(nbk, dtype=np.int64)
    MOM = np.zeros((nbk, 8, 2), dtype=np.int64)
    FL = np.zeros((nbk, 64), dtype=np.int64)
    FLC = np.zeros(nbk, dtype=np.int64)
    block_pass(tab, np.uint64(1 << TAB), c16, d16, p3, blist, grid, glo, ghi, ovn, ovt,
               LOC, HEND, SH1, SH2, MOM, GLOC, GH, GSH1, GSH2, FL, FLC)
    return LOC, HEND, SH1, SH2, MOM, FL, FLC


if __name__ == '__main__':
    t0 = time.time()
    tab = build_table()
    print(f'table of 2^{TAB} parities: {time.time() - t0:.1f} s', flush=True)
    # self-checks of the jump descent against exact Python trajectories
    rng = random.Random(15102026)
    lim = np.uint64(1 << TAB)

    @nb.njit
    def single(n, tab, lim, c16, d16, p3):
        return par_of(n, tab, lim, c16, d16, p3)
    bad = 0
    nchk = 0
    for _ in range(4000):
        n = rng.randrange(1, 1 << 37)
        pb = single(n, tab, lim, c16, d16, p3)
        if pb >= 0:
            nchk += 1
            bad += pb != tau_parity_py(n)
    for n in range(1, 1 << 15):
        nchk += 1
        bad += single(n, tab, lim, c16, d16, p3) != tau_parity_py(n)
    print(f'self-check: {nchk} values against exact Python trajectories, {bad} mismatches', flush=True)
    if bad:
        raise SystemExit('self-check failed')
    # grid
    if GRIDNPZ:
        with np.load(GRIDNPZ) as Z:
            grid = np.array(Z['grid'], dtype=np.int64)
    else:
        raw = sorted(set(int(round(2 ** (10 + j / 64))) for j in range(0, 27 * 64 + 1)))
        grid = np.array([g for g in raw if g < (1 << 37)], dtype=np.int64)
    NB_ALL = (1 << 37) // B
    bidx = grid // B
    glo = np.searchsorted(bidx, np.arange(NB_ALL), side='left').astype(np.int64)
    ghi = np.searchsorted(bidx, np.arange(NB_ALL), side='right').astype(np.int64)
    G = len(grid)
    GLOC = np.zeros((G, M), dtype=np.int64)
    GH = np.zeros(G, dtype=np.int64)
    GSH1 = np.zeros(G, dtype=np.int64)
    GSH2 = np.zeros(G, dtype=np.int64)
    b0 = 0 if LOGLO == 0 else (1 << LOGLO) // B
    b1 = (1 << LOGHI) // B
    blist = np.arange(b0, b1, dtype=np.int64)
    empty_n = np.zeros(0, dtype=np.int64)
    empty_t = np.zeros(0, dtype=np.int8)
    t1 = time.time()
    LOC, HEND, SH1, SH2, MOM, FL, FLC = run_blocks(tab, blist, grid, glo, ghi, empty_n, empty_t, GLOC, GH, GSH1, GSH2)
    print(f'block pass over blocks [{b0}, {b1}): {time.time() - t1:.1f} s; flagged n (orbit above 2^63 before 2^{TAB}): '
          f'{int(FLC.sum())} in {int((FLC > 0).sum())} blocks', flush=True)
    if (FLC > 64).any():
        raise SystemExit('too many flagged values in a block')
    flagged = sorted(int(x) for bi in np.nonzero(FLC)[0] for x in FL[bi, :FLC[bi]])
    ovt_list = [1 - 2 * tau_parity_py(n) for n in flagged]
    ovn = np.array(flagged, dtype=np.int64)
    ovt = np.array(ovt_list, dtype=np.int8)
    dirty = np.nonzero(FLC)[0]
    if len(dirty):
        sub = blist[dirty]
        L2, HE2, S12, S22, MO2, FL2, FC2 = run_blocks(tab, sub, grid, glo, ghi, ovn, ovt, GLOC, GH, GSH1, GSH2)
        if FC2.sum():
            raise SystemExit('unresolved overflow values')
        LOC[dirty] = L2; HEND[dirty] = HE2; SH1[dirty] = S12; SH2[dirty] = S22; MOM[dirty] = MO2
    print(f'patched {len(dirty)} blocks with {len(flagged)} exact Python values; '
          f'largest flagged n = {max(flagged) if flagged else None}', flush=True)
    gsel = (bidx >= b0) & (bidx < b1)
    np.savez_compressed(OUT, b0=b0, b1=b1, LOC=LOC, HEND=HEND, SH1=SH1, SH2=SH2, MOM=MOM, grid=grid,
                        gsel=gsel, GLOC=GLOC, GH=GH, GSH1=GSH1, GSH2=GSH2, flagged=ovn, flagged_t=ovt)
    print(f'wrote {OUT}; total {time.time() - t0:.1f} s', flush=True)
