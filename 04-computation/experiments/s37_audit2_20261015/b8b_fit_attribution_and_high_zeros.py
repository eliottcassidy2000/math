#!/usr/bin/env python3
"""S37 audit 2, script b8b: two claims of note section 8 / HYP-9285 about the explicit-formula fit.

 (1) "The large step at Re > 0.65 comes from rho = 0.6939 + 22.25i": holdout (fit [2^22, 2^34), test [2^34, 2^37), session
     grid) of the 12-zero model plus 0.6939 + 22.25i alone, and of the 28-zero model minus that zero.
 (2) "The remaining error is consistent with zeros above Im 3000, whose oscillations alias on this grid": on an exact
     fine grid (1024 points per octave, Nyquist 512 cycles/octave = Im 4640) on [2^22, 2^31) (independent census to 2^31),
     add the zeros with 0.65 < Re < 1 and 3000 < Im < 4600 to the 28-zero model and compare in-sample and holdout errors.
Usage: python b8b_fit_attribution_and_high_zeros.py census_ext.npz
"""
import sys, math
import numpy as np
import numba as nb

CEN = np.load(sys.argv[1])
L2, L15 = math.log(2), math.log(1.5)
D = lambda s: 1 + np.exp(-s * L2) + np.exp(s * L15) / 3
Dp = lambda s: -L2 * np.exp(-s * L2) + L15 * np.exp(s * L15) / 3


def zeros(imax, rmin):
    zs = []
    for t in np.arange(0.5, imax, 0.5):
        for sg in (0.7, 0.85, 0.95, 1.0):
            s = complex(sg, t)
            for _ in range(60):
                ds = D(s) / Dp(s); s -= ds
                if abs(ds) < 1e-15:
                    break
            if abs(D(s)) < 1e-10 and 0 < s.imag < imax and rmin < s.real < 1 and all(abs(s - z) > 1e-7 for z in zs):
                zs.append(s)
    return sorted(zs, key=lambda z: -z.real)


def design(Nv, zz):
    u = np.log(Nv); cols = []
    for z in zz:
        w = np.exp(z.real * u) * np.exp(1j * z.imag * u); cols += [w.real, -w.imag]
    return np.stack(cols, axis=1)


def holdout(Nv, Sv, zz, fm, tm):
    A = design(Nv, zz)
    c, *_ = np.linalg.lstsq(A[fm], Sv[fm], rcond=None)
    p = A @ c
    r = lambda m: float(np.sqrt(((Sv[m] - p[m]) ** 2).mean()) / np.sqrt((Sv[m] ** 2).mean()))
    return r(fm), r(tm)


Z46 = zeros(4600.0, 0.65)
z28 = [z for z in Z46 if z.imag < 3000]
zhi = [z for z in Z46 if z.imag >= 3000]
print(f'zeros with 0.65 < Re < 1: {len(z28)} below Im 3000, {len(zhi)} in [3000, 4600)')
g = CEN['grid'].astype(np.float64); S = CEN['classsums'].sum(axis=1).astype(np.float64)
fm = (g >= 2 ** 22) & (g < 2 ** 34); tm = (g >= 2 ** 34) & (g < 2 ** 37)
z12 = [z for z in z28 if z.real > 0.85]
zlow = [z for z in z28 if abs(z.imag - 22.2532) < 0.01]
print('(1) holdout on the session grid (in-sample, test):')
for lab, zz in (('12 zeros', z12), ('12 zeros + 0.6939+22.25i', z12 + zlow),
                ('28 zeros', z28), ('28 zeros - 0.6939+22.25i', [z for z in z28 if z not in zlow])):
    a, b = holdout(g, S, zz, fm, tm)
    print(f'  {lab:28s}: {a:.4f}, {b:.4f}')

# fine grid census to 2^31
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
def octave(tab, w0, w1, lo, c16, d16, p3):
    lim = np.uint64(lo); J = np.uint64(1) << np.uint64(54); one = np.uint64(1)
    for w in nb.prange(w0, w1):
        word = np.uint64(0)
        for b in range(64):
            x = np.uint64((w << 6) + b); p = 0
            while x >= lim:
                if x < J:
                    l = x & np.uint64(0xFFFF)
                    x = p3[c16[l]] * (x >> np.uint64(16)) + d16[l]
                else:
                    if x & one:
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
    out = np.zeros(Ns.shape[0], dtype=np.int64); acc = 0; w = 0
    for i in range(Ns.shape[0]):
        N = Ns[i]; wN = N >> 6
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


tab = np.zeros((1 << LOG) // 64 + 1, dtype=np.uint64)
word = 0
for n in range(1, 64):
    if tau_par(n):
        word |= 1 << n
tab[0] = np.uint64(word)
lo = 64
while lo < (1 << LOG):
    octave(tab, lo >> 6, (2 * lo) >> 6, lo, c16, d16, p3); lo *= 2
g31 = CEN['grid'][CEN['grid'] < (1 << LOG)].astype(np.int64)
assert np.array_equal(S_at(tab, g31), CEN['classsums'][CEN['grid'] < (1 << LOG)].sum(axis=1))
fine = np.array(sorted(set(int(round(2 ** (22 + j / 1024))) for j in range(0, 9 * 1024))), dtype=np.int64)
Sf = S_at(tab, fine).astype(np.float64); Nf = fine.astype(np.float64)
ffm = Nf < 2 ** 29; ftm = Nf >= 2 ** 29
print('(2) fine grid (1024/octave) on [2^22, 2^31): fit [2^22, 2^29), test [2^29, 2^31) (in-sample, test)')
for lab, zz in (('28 zeros (Im < 3000)', z28), (f'{len(Z46)} zeros (Im < 4600)', Z46)):
    a, b = holdout(Nf, Sf, zz, ffm, ftm)
    a2, _ = holdout(Nf, Sf, zz, np.ones(len(Nf), bool), np.ones(len(Nf), bool))
    print(f'  {lab:22s}: in-sample {a:.4f}, holdout {b:.4f}; full-range in-sample {a2:.4f}')
