#!/usr/bin/env python3
"""S37 audit 2, script b8: HYP-9285 least-squares fits, holdouts, residue estimates, and grid aliasing.

 (1) The 44 zeros of D with 1/2 < Re < 1, 0 < Im < 3000 (Newton; counts verified by b7).
 (2) Holdouts on the session grid (64 points per octave, S from the verified 2^37 census NPZ): fit
     S(N) ~ Re sum_rho c_rho N^rho on [2^22, 2^34), relative rms on [2^34, 2^37), for the zero sets s0 only,
     Re > 0.85 (12), Re > 0.65 (28), Re > 0.5 (44); in-sample relative rms.
 (3) The 28-zero fit on [2^22, 2^37): coefficients versus collatz_mercator_waveform_prediction_20261015.json and
     |r_fit| = |c| |rho| / 2 versus collatz_other_zeros_residue_table_20261015.json.
 (4) Aliasing: the grid has 64 points per octave, so frequencies above 32 cycles per octave (Im rho > 290) are folded;
     folded frequency of every model zero.
 (5) Fine-grid control: exact S(N) at 1024 points per octave on [2^22, 2^31) from an independent census to 2^31
     (bit table, prefix popcounts); the same 28-zero model fitted on the coarse and on the fine grid of the same range;
     |r_fit| per zero compared; in-sample relative rms on both grids; coarse-grid vs fine-grid holdout [2^29, 2^31).
Usage: python b8_explicit_formula_fit.py census_ext.npz prediction.json residue_table.json
"""
import sys, json, math, time
import numpy as np
import numba as nb

CEN = np.load(sys.argv[1]); PRED = json.load(open(sys.argv[2])); RT = json.load(open(sys.argv[3]))
L2, L15 = math.log(2), math.log(1.5)


def D(s):
    return 1 + np.exp(-s * L2) + np.exp(s * L15) / 3


def Dp(s):
    return -L2 * np.exp(-s * L2) + L15 * np.exp(s * L15) / 3


zs = []
for t in np.arange(0.5, 3000.0, 0.5):
    for sg in (0.55, 0.7, 0.85, 0.95, 1.0):
        s = complex(sg, t)
        for _ in range(60):
            ds = D(s) / Dp(s)
            s -= ds
            if abs(ds) < 1e-15:
                break
        if abs(D(s)) < 1e-11 and 0 < s.imag < 3000 and 0.5 < s.real < 1.0 and all(abs(s - z) > 1e-7 for z in zs):
            zs.append(s)
zs.sort(key=lambda z: -z.real)
print(f'(1) zeros with 1/2 < Re < 1, 0 < Im < 3000: {len(zs)}')

grid = CEN['grid'].astype(np.float64)
S = CEN['classsums'].sum(axis=1).astype(np.float64)


def design(Nv, zeros):
    u = np.log(Nv)
    cols = []
    for z in zeros:
        w = np.exp(z.real * u) * np.exp(1j * z.imag * u)
        cols += [w.real, -w.imag]
    return np.stack(cols, axis=1)


def fit(Nv, Sv, zeros, mask):
    A = design(Nv, zeros)
    coef, *_ = np.linalg.lstsq(A[mask], Sv[mask], rcond=None)
    return coef, A @ coef


def relrms(res, ref):
    return float(np.sqrt((res ** 2).mean()) / np.sqrt((ref ** 2).mean()))


sets = {'s0 only': [z for z in zs if z.real > 0.998], 'Re > 0.85': [z for z in zs if z.real > 0.85],
        'Re > 0.65': [z for z in zs if z.real > 0.65], 'Re > 0.5': zs}
fitm = (grid >= 2 ** 22) & (grid < 2 ** 34)
test = (grid >= 2 ** 34) & (grid < 2 ** 37)
print('(2) holdout: fit [2^22, 2^34), test [2^34, 2^37) on the session grid')
for lab, zz in sets.items():
    coef, pr = fit(grid, S, zz, fitm)
    print(f'  {lab:10s} ({len(zz):2d} zeros, {2 * len(zz)} real parameters): in-sample rel rms {relrms(S[fitm] - pr[fitm], S[fitm]):.4f}; '
          f'holdout rel rms {relrms(S[test] - pr[test], S[test]):.4f}')
# robustness of the 28-zero holdout to the fit start
for lo in (20, 24, 26):
    m = (grid >= 2 ** lo) & (grid < 2 ** 34)
    coef, pr = fit(grid, S, sets['Re > 0.65'], m)
    print(f'  Re > 0.65, fit from 2^{lo}: holdout rel rms {relrms(S[test] - pr[test], S[test]):.4f}')

print('(3) 28-zero fit on [2^22, 2^37) versus the session prediction JSON and residue table')
z28 = sets['Re > 0.65']
m37 = (grid >= 2 ** 22) & (grid < 2 ** 37)
coef, pr = fit(grid, S, z28, m37)
pz = [complex(float(a), float(b)) for a, b in PRED['zeros']]
pc = [complex(a, b) for a, b in PRED['coefficients']]
mx = 0.0
for i, z in enumerate(z28):
    k = int(np.argmin([abs(z - w) for w in pz]))
    mx = max(mx, abs(complex(coef[2 * i], coef[2 * i + 1]) - pc[k]) / abs(pc[k]))
print(f'  in-sample rel rms {relrms(S[m37] - pr[m37], S[m37]):.4f} (session {PRED["insample_relative_rms_residual"]:.4f}); '
      f'max relative difference of coefficients from the session JSON {mx:.1e}')
rt = {(round(r['zero'][0], 4), round(r['zero'][1], 2)): r for r in RT['rows']}
rfit37 = {}
for i, z in enumerate(z28):
    c = complex(coef[2 * i], coef[2 * i + 1])
    rfit37[i] = abs(c) * abs(z) / 2

print('(4) aliasing on the 64-per-octave grid (Nyquist 32 cycles/octave = Im 290.1)')
f_of = lambda z: z.imag * L2 / (2 * math.pi)
fold = lambda f: abs(((f + 32) % 64) - 32)
above = [z for z in z28 if f_of(z) > 32]
print(f'  {len(above)} of the 28 model zeros lie above the grid Nyquist frequency; folded frequencies (cycles/octave):')
print('  ' + ', '.join(f'{f_of(z):.2f}->{fold(f_of(z)):.2f}' for z in sorted(z28, key=lambda z: z.imag)))

print('(5) fine-grid control on [2^22, 2^31)')
LOG = 31
c16 = np.zeros(1 << 16, dtype=np.uint8); d16 = np.zeros(1 << 16, dtype=np.uint64)
for l in range(1 << 16):
    y, cc = l, 0
    for _ in range(16):
        if y & 1:
            y = (3 * y + 1) // 2; cc += 1
        else:
            y //= 2
    c16[l] = cc; d16[l] = y
p3 = np.array([3 ** k for k in range(17)], dtype=np.uint64)


@nb.njit(parallel=True)
def octave(tab, w0, w1, lo, c16, d16, p3, bad):
    lim = np.uint64(lo)
    J = np.uint64(1) << np.uint64(54)
    G = np.uint64(1) << np.uint64(63)
    one = np.uint64(1)
    for w in nb.prange(w0, w1):
        word = np.uint64(0)
        for b in range(64):
            x = np.uint64((w << 6) + b)
            p = 0
            while x >= lim:
                if x < J:
                    l = x & np.uint64(0xFFFF)
                    x = p3[c16[l]] * (x >> np.uint64(16)) + d16[l]
                else:
                    if x & one:
                        if x >= G:
                            bad[0] = 1
                            break
                        x = x + ((x + one) >> one)
                    else:
                        x = x >> one
                    p ^= 1
            xi = np.int64(x)
            p ^= int((tab[xi >> 6] >> np.uint64(xi & 63)) & one)
            if p:
                word |= one << np.uint64(b)
        tab[w] = word


@nb.njit
def popc(x):
    x = x - ((x >> np.uint64(1)) & np.uint64(0x5555555555555555))
    x = (x & np.uint64(0x3333333333333333)) + ((x >> np.uint64(2)) & np.uint64(0x3333333333333333))
    x = (x + (x >> np.uint64(4))) & np.uint64(0x0F0F0F0F0F0F0F0F)
    return int((x * np.uint64(0x0101010101010101)) >> np.uint64(56))


@nb.njit
def S_at(tab, Ns):
    # S(N) = N - 2 #{1 <= n <= N : parity bit 1}; bit 0 (n = 0) is 0
    out = np.zeros(Ns.shape[0], dtype=np.int64)
    acc = 0
    w = 0
    for i in range(Ns.shape[0]):
        N = Ns[i]
        wN = N >> 6
        while w < wN:
            acc += popc(tab[w]); w += 1
        part = popc(tab[wN] & ((np.uint64(2) << np.uint64(N & 63)) - np.uint64(1))) if (N & 63) != 63 else popc(tab[wN])
        out[i] = N - 2 * (acc + part)
    return out


def tau_par(n):
    s = 0
    while n != 1:
        n = n // 2 if n % 2 == 0 else (3 * n + 1) // 2
        s += 1
    return s & 1


t0 = time.time()
tab = np.zeros((1 << LOG) // 64 + 1, dtype=np.uint64)
word = 0
for n in range(1, 64):
    if tau_par(n):
        word |= 1 << n
tab[0] = np.uint64(word)
bad = np.zeros(1, dtype=np.int64)
lo = 64
while lo < (1 << LOG):
    octave(tab, lo >> 6, (2 * lo) >> 6, lo, c16, d16, p3, bad)
    lo *= 2
assert not bad[0]
# consistency with the verified census on the session grid below 2^31
g31 = CEN['grid'][CEN['grid'] < (1 << LOG)].astype(np.int64)
Sg = S_at(tab, g31)
print(f'  census to 2^31 ({time.time() - t0:.0f} s); S on the session grid below 2^31 identical to the 2^37 census: '
      f'{np.array_equal(Sg, CEN["classsums"][CEN["grid"] < (1 << LOG)].sum(axis=1))}')
fine = np.array(sorted(set(int(round(2 ** (22 + j / 1024))) for j in range(0, 9 * 1024))), dtype=np.int64)
Sf = S_at(tab, fine).astype(np.float64)
Nf = fine.astype(np.float64)
mc = (grid >= 2 ** 22) & (grid < 2 ** 31)
coef_c, pr_c = fit(grid, S, z28, mc)
coef_f, pr_f = fit(Nf, Sf, z28, np.ones(len(fine), dtype=bool))
print(f'  28-zero fit in-sample rel rms: coarse grid {relrms(S[mc] - pr_c[mc], S[mc]):.4f}, fine grid {relrms(Sf - pr_f, Sf):.4f}; '
      f'coarse-grid model evaluated on the fine grid: {relrms(Sf - design(Nf, z28) @ coef_c, Sf):.4f}')
print('  |r_fit| per zero: [2^22,2^37) coarse (session) | [2^22,2^31) coarse | [2^22,2^31) fine | conditional interval (2,0)')
for i, z in enumerate(sorted(range(len(z28)), key=lambda i: z28[i].imag)):
    zz = z28[z]
    rc = abs(complex(coef_c[2 * z], coef_c[2 * z + 1])) * abs(zz) / 2
    rf = abs(complex(coef_f[2 * z], coef_f[2 * z + 1])) * abs(zz) / 2
    key = (round(zz.real, 4), round(zz.imag, 2))
    row = rt.get(key)
    iv = row['abs_r_interval_2_0'] if row else None
    flag = '' if f_of(zz) <= 32 else ' (aliased)'
    print(f'    {zz.real:.5f}+{zz.imag:9.3f}i{flag}: {rfit37[z]:.4f} | {rc:.4f} | {rf:.4f} | '
          f'[{iv[0]:.4f}, {iv[1]:.4f}]' + (' certified' if row and row['nonzero_certified'] else ''))
# holdout on both grids over the same ranges
mfc = (grid >= 2 ** 22) & (grid < 2 ** 29); mtc = (grid >= 2 ** 29) & (grid < 2 ** 31)
cc_, pc_ = fit(grid, S, z28, mfc)
mff = Nf < 2 ** 29; mtf = Nf >= 2 ** 29
cf_, pf_ = fit(Nf, Sf, z28, mff)
print(f'  holdout fit [2^22,2^29) test [2^29,2^31): coarse grid {relrms(S[mtc] - pc_[mtc], S[mtc]):.4f}; fine grid '
      f'{relrms(Sf[mtf] - pf_[mtf], Sf[mtf]):.4f}; coarse-fitted model tested on the fine grid {relrms(Sf[mtf] - (design(Nf, z28) @ cc_)[mtf], Sf[mtf]):.4f}')
for lab, zz in (('s0 only', sets['s0 only']), ('Re > 0.85', sets['Re > 0.85']), ('below-Nyquist zeros only (3)',
                [z for z in z28 if f_of(z) <= 32])):
    cfa, pfa = fit(Nf, Sf, zz, mff)
    cca, pca = fit(grid, S, zz, mfc)
    print(f'  {lab}: holdout [2^29,2^31) fine grid {relrms(Sf[mtf] - pfa[mtf], Sf[mtf]):.4f}; coarse grid {relrms(S[mtc] - pca[mtc], S[mtc]):.4f}')
