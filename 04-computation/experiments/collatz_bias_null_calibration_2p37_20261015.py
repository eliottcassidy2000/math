#!/usr/bin/env python3
"""Null calibration of HYP-9282 / HYP-9283 at the scale of the 2^37 census, with identical statistics for Collatz.

Null ('halving null'): t(n) = (-1)^v2(n) s(odd part of n), with s(m) = +-1 from a splitmix64 hash of (m, seed).  It keeps
the exact relation t(2m) = -t(m) of the Collatz parity, hence the octave inheritance C(A; N) = -C(A/2; N/2) + (odd part)
of every class sum, and has no odd-step structure.  For each seed the script computes, block-parallel and exactly as in
collatz_census_extension_2p37_20261015.py, the class sums modulo 243 and the prefix sums of H, H^2 at the grid
N_j = round(2^(10+j/64)) < 2^LOGMAX, and then the statistics below; the Collatz census NPZ goes through the same code.

Statistics (in-sample octaves 14..32, out-of-sample octaves 34..36, as pre-registered in HYP-9282):
  naive t-values of the 363 channel octave means and the counts |t| > 3, 5, 8; rms octave mean and sd by level;
  pooled octave-mean autocorrelations by level (octaves 14..36); top-5 |t| channels followed out of sample (signs kept
  per octave, majority rule, averages within +-0.6 of the in-sample means); channels with |t| > 5 followed out of sample
  (slope through the origin, correlation, sign of the out-of-sample mean); profile g = H/sqrt(2N) (mean, sd, lag-1-octave
  autocorrelation, two-octave regression slope); window energies e = E/Y^2 of H (mean, max for window starts up to 2^36);
  per-octave mean square of the normalized channels by level, and its out/in (34..36 vs 14..32) and late/early
  (30..36 vs 14..29) ratios; the number of channels whose 19 in-sample octave means share one sign.
The octave-mean matrices (octaves 14..36 x 363 channels), level energies and H profiles are saved to an NPZ next to the
JSON (row 0 = Collatz when a census NPZ is given).
Session opus-2026-10-15-S37.  Usage: python collatz_bias_null_calibration_2p37_20261015.py LOGMAX NSEEDS collatz.npz out.json
"""
import sys, json, math, time
import numpy as np
import numba as nb

LOGMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 37
NSEEDS = int(sys.argv[2]) if len(sys.argv) > 2 else 16
CNPZ = sys.argv[3] if len(sys.argv) > 3 else None
OUT = sys.argv[4] if len(sys.argv) > 4 else 'null_calibration.json'
M = 243
BLOCK = 1 << 20
INS = list(range(14, 33))
OUTK = [34, 35, 36]


@nb.njit
def tnull(n, seed):
    v = 0
    m = n
    while (m & 1) == 0:
        m >>= 1
        v += 1
    z = np.uint64(m) ^ seed
    z = z + np.uint64(0x9E3779B97F4A7C15)
    z = (z ^ (z >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
    z = (z ^ (z >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
    z = z ^ (z >> np.uint64(31))
    s = 1 - 2 * int(z & np.uint64(1))
    return s if (v & 1) == 0 else -s


@nb.njit(parallel=True)
def null_counts(seed, nblocks, B, Mod):
    cnt = np.zeros((nblocks, Mod), dtype=np.int64)
    for b in nb.prange(nblocks):
        lo = b * B
        n0 = lo if lo > 0 else 1
        r = n0 % Mod
        for n in range(n0, lo + B):
            cnt[b, r] += tnull(n, seed)
            r += 1
            if r == Mod:
                r = 0
    return cnt


@nb.njit(parallel=True)
def null_pass(seed, nblocks, B, Mod, start, grid, glo, ghi):
    G = grid.shape[0]
    snapC = np.zeros((G, Mod), dtype=np.int64)
    snapH1 = np.zeros(G, dtype=np.int64)
    snapH2 = np.zeros(G, dtype=np.int64)
    snapH = np.zeros(G, dtype=np.int64)
    totH1 = np.zeros(nblocks, dtype=np.int64)
    totH2 = np.zeros(nblocks, dtype=np.int64)
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
        gi = glo[b]
        ge = ghi[b]
        for n in range(n0, lo + B):
            t = tnull(n, seed)
            C[r] += t
            C3[r % 3] += t
            S += t
            h = 3 * C3[2] - S
            s1 += h
            s2 += h * h
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
    return snapC, snapH1, snapH2, snapH, totH1, totH2


def make_grid(logmax):
    raw = sorted(set(int(round(2 ** (10 + j / 64))) for j in range(0, (logmax - 10) * 64 + 1)))
    return np.array([g for g in raw if g < (1 << logmax)], dtype=np.int64)


def null_census(seed, logmax):
    nblocks = (1 << logmax) // BLOCK
    cnt = null_counts(np.uint64(seed), nblocks, BLOCK, M)
    start = np.zeros_like(cnt)
    start[1:] = np.cumsum(cnt[:-1], axis=0)
    grid = make_grid(logmax)
    bidx = grid // BLOCK
    glo = np.searchsorted(bidx, np.arange(nblocks), side='left').astype(np.int64)
    ghi = np.searchsorted(bidx, np.arange(nblocks), side='right').astype(np.int64)
    snapC, sH1, sH2, sH, tH1, tH2 = null_pass(np.uint64(seed), nblocks, BLOCK, M, start, grid, glo, ghi)
    c1 = np.concatenate([[0.0], np.cumsum(tH1.astype(np.float64))])
    c2 = np.concatenate([[0.0], np.cumsum(tH2.astype(np.float64))])
    pre1 = c1[bidx] + sH1
    pre2 = c2[bidx] + sH2
    return dict(grid=grid, classsums=snapC, H=sH, preH1=pre1, preH2=pre2)


def statistics(D):
    grid = D['grid'].astype(np.int64)
    C = D['classsums'].astype(np.float64)
    N = grid.astype(np.float64)
    u = np.log2(N)
    lev = {0: C.sum(axis=1)[:, None]}
    for q in range(1, 6):
        m = 3 ** q
        lev[q] = np.stack([C[:, a::m].sum(axis=1) for a in range(m)], axis=1)
    cols, levs = [], []
    for q in range(1, 6):
        m, mp = 3 ** q, 3 ** (q - 1)
        for r in range(m):
            cols.append((3 * lev[q][:, r] - lev[q - 1][:, r % mp]) / np.sqrt(6 * N / m))
            levs.append(q)
    levs = np.array(levs)
    Gm = np.stack(cols, axis=1)
    ks = list(range(14, 37))
    OM = np.stack([Gm[(u >= k) & (u < k + 1)].mean(axis=0) for k in ks])
    ii = [ks.index(k) for k in INS]
    oo = [ks.index(k) for k in OUTK]
    mi = OM[ii].mean(axis=0)
    si = OM[ii].std(axis=0, ddof=1)
    ti = mi / (si / math.sqrt(len(ii)))
    st = dict(n_t3=int((np.abs(ti) > 3).sum()), n_t5=int((np.abs(ti) > 5).sum()), n_t8=int((np.abs(ti) > 8).sum()),
              max_abs_t=float(np.abs(ti).max()),
              n_same_sign_19=int(((OM[ii] > 0).all(axis=0) | (OM[ii] < 0).all(axis=0)).sum()))
    st['rms_mean_by_level'] = [float(np.sqrt((mi[levs == q] ** 2).mean())) for q in range(1, 6)]
    # per-octave mean square of the normalized channels (all samples) by level; out/in and late/early ratios
    msq = np.stack([(Gm[(u >= k) & (u < k + 1)] ** 2).mean(axis=0) for k in ks])     # octaves x channels
    lev_ms = np.stack([msq[:, levs == q].mean(axis=1) for q in range(1, 6)], axis=1)  # octaves x levels
    st['channel_energy_by_octave_level'] = [[float(x) for x in row] for row in lev_ms]
    st['energy_ratio_out_in_by_level'] = [float(lev_ms[oo, q].mean() / lev_ms[ii, q].mean()) for q in range(5)]
    late = [ks.index(k) for k in range(30, 37)]
    early = [ks.index(k) for k in range(14, 30)]
    st['energy_ratio_late_early_by_level'] = [float(lev_ms[late, q].mean() / lev_ms[early, q].mean()) for q in range(5)]
    st['energy_ratio_out_in_levels_2_5'] = float(lev_ms[oo, 1:].mean() / lev_ms[ii, 1:].mean())
    st['sd_by_level'] = [float(np.sqrt((si[levs == q] ** 2).mean())) for q in range(1, 6)]
    X = OM - OM.mean(axis=0)
    v = (X ** 2).mean(axis=0)
    st['acf_lag1_by_level'] = [float(np.mean((X[1:, levs == q] * X[:-1, levs == q]).mean(axis=0) / v[levs == q])) for q in range(1, 6)]
    st['acf_lag2_by_level'] = [float(np.mean((X[2:, levs == q] * X[:-2, levs == q]).mean(axis=0) / v[levs == q])) for q in range(1, 6)]
    top = list(np.argsort(-np.abs(ti))[:5])
    kept = 0
    maj = 0
    band = 0
    for c in top:
        signs = [np.sign(OM[i, c]) == np.sign(mi[c]) for i in oo]
        kept += sum(signs)
        maj += int(sum(signs) >= 2)
        band += int(abs(OM[oo, c].mean() - mi[c]) <= 0.6)
    st['top5'] = dict(signs_kept_of_15=int(kept), majority_kept_of_5=int(maj), in_band_of_5=int(band),
                      passes_all=bool(maj == 5 and band == 5))
    strong = np.nonzero(np.abs(ti) > 5)[0]
    if len(strong):
        mo = OM[oo][:, strong].mean(axis=0)
        st['t5_follow'] = dict(n=int(len(strong)), sign_kept=int((np.sign(mo) == np.sign(mi[strong])).sum()),
                               slope=float((mo * mi[strong]).sum() / (mi[strong] ** 2).sum()))
    else:
        st['t5_follow'] = dict(n=0, sign_kept=0, slope=None)
    mo_all = OM[oo].mean(axis=0)
    st['all_follow'] = dict(slope=float((mo_all * mi).sum() / (mi ** 2).sum()), corr=float(np.corrcoef(mo_all, mi)[0, 1]))
    # level-1 channel H = H_(1,2) and profile
    H = np.asarray(D['H'], dtype=np.float64)
    g = H / np.sqrt(2 * N)
    sel = u >= 14
    gs = g[sel]
    x = gs - gs.mean()
    st['profile'] = dict(mean=float(gs.mean()), sd=float(gs.std()), acf_1oct=float((x[64:] * x[:-64]).mean() / x.var()),
                         slope_2oct=float(np.polyfit(gs[:-128], gs[128:], 1)[0]), mean_in=float(g[(u >= 14) & (u < 33)].mean()),
                         mean_out=float(g[(u >= 34) & (u < 37)].mean()))
    p1 = np.asarray(D['preH1'], dtype=np.float64)
    p2 = np.asarray(D['preH2'], dtype=np.float64)
    jj = np.arange(len(grid) - 64)
    Y = N[jj]
    E = (p2[jj + 64] - H[jj + 64] ** 2) - (p2[jj] - H[jj] ** 2)
    e = E / Y ** 2
    uY = np.log2(Y)
    w_in = (uY >= 14) & (uY <= 33 + 1e-9)
    w_out = uY > 33 + 1e-9
    st['window'] = dict(mean_in=float(e[w_in].mean()), max_in=float(e[w_in].max()), mean_out=float(e[w_out].mean()),
                        max_out=float(e[w_out].max()))
    return st, OM, lev_ms, g


if __name__ == '__main__':
    res = dict(logmax=LOGMAX, insample=INS, out=OUTK, null=[])
    OMs, LMs, GPs = [], [], []
    if CNPZ:
        with np.load(CNPZ) as Z:
            D = {k: np.array(Z[k]) for k in ('grid', 'classsums', 'H', 'preH1', 'preH2')}
        st, OM, LM, gp = statistics(D)
        res['collatz'] = st
        OMs.append(OM); LMs.append(LM); GPs.append(gp)
        print('collatz', json.dumps(st), flush=True)
    for s in range(NSEEDS):
        t0 = time.time()
        D = null_census(7919 * (s + 1), LOGMAX)
        st, OM, LM, gp = statistics(D)
        res['null'].append(st)
        OMs.append(OM); LMs.append(LM); GPs.append(gp)
        print(f'null {s} ({time.time() - t0:.0f} s)', json.dumps(st), flush=True)
        json.dump(res, open(OUT, 'w'), indent=1)
        np.savez_compressed(OUT.replace('.json', '.npz'), octave_means=np.array(OMs), level_energy=np.array(LMs),
                            profile=np.array(GPs), first_is_collatz=bool(CNPZ))
    json.dump(res, open(OUT, 'w'), indent=1)
