#!/usr/bin/env python3
"""S37 audit 2, script b1: THM-4658 statement 7 (finite-depth spectrum of the 3-adic transfer operator).

Independent construction (does not import the session's scan script):
 (a) density form from the measure form of statement 4: for a density f constant modulo 3^Q put
     beta(b + 3^q Z_3) = 3^(-q) f(b mod 3^Q) (q >= Q) and evaluate (P beta)(A) = -2^(-s) beta(A/2) - (3/2)^s beta((3A+1)/2)
     on every cylinder A of level q = Q, Q+1, Q+2, computing the image cylinders by brute force (images of all
     representatives); compare with 3^(-q) (P_Q f)(a mod 3^Q) for the matrix of statement 7.
 (b) spectra of P_Q(s) (all densities) and of P~_Q(s) on mean-zero densities for Q = 1..5 at several s, matched against
     {1 - D(s)} u {-2^(-s) zeta: zeta^(2 3^k) = 1, k < Q} with the block multiplicities of the proof
     (one full set of 2 3^(Q-j-1)-th roots for each j); characteristic polynomial compared at random points
     (numpy for Q <= 5, mpmath at 60 digits for Q <= 3).
 (c) Jordan structure at the repeated eigenvalues -2^(-s) (zeta = 1) and +2^(-s) (zeta = -1): geometric multiplicities.
 (d) mean-zero invariance of P~, eigenvalue 1 of P_Q at the Mercator zero and the distance of spec(P~) from 1 there.
"""
import itertools, cmath, math
import numpy as np
import mpmath as mp


def mats(Q, s, excess):
    m = 3 ** Q
    i2 = pow(2, -1, m)
    P = np.zeros((m, m), dtype=complex)
    a2 = -(2.0 ** (-s))
    b3 = -(1.5 ** s) / 3
    for a in range(m):
        P[a, (a * i2) % m] += a2
        P[a, ((3 * a + 1) * i2) % m] += b3
        if excess:
            for c in range(2, m, 3):
                P[a, c] += (1.5 ** s) / m
    return P


def Dfun(s):
    return 1 + 2.0 ** (-s) + (1.5 ** s) / 3


def predicted(Q, s, full):
    ev = []
    for j in range(Q):
        L = 2 * 3 ** (Q - j - 1)
        ev += [-(2.0 ** (-s)) * cmath.exp(2j * math.pi * r / L) for r in range(L)]
    if full:
        ev.append(1 - Dfun(s))
    return np.array(ev)


def match_err(ev, pred):
    # optimal assignment by brute greedy on sorted distance matrix (sizes <= 243); report max matched distance
    ev = list(ev); pred = list(pred)
    used = [False] * len(ev)
    worst = 0.0
    # greedy over all pairs sorted by distance
    pairs = sorted((abs(e - p), i, k) for i, p in enumerate(pred) for k, e in enumerate(ev))
    donep = [False] * len(pred)
    cnt = 0
    for d, i, k in pairs:
        if donep[i] or used[k]:
            continue
        donep[i] = True; used[k] = True
        worst = max(worst, d); cnt += 1
        if cnt == len(pred):
            break
    return worst


def meanzero_basis(m):
    # orthonormal Fourier basis of the mean-zero functions on Z/m (columns e_k / sqrt(m), k = 1..m-1)
    a = np.arange(m)
    return np.stack([np.exp(2j * np.pi * k * a / m) / math.sqrt(m) for k in range(1, m)], axis=1)


print('(a) density form from the measure form (brute-force image cylinders)')
rng = np.random.default_rng(1)
for Q in (1, 2, 3):
    s = 0.73 + 5.1j
    f = rng.normal(size=3 ** Q) + 1j * rng.normal(size=3 ** Q)
    PQ = mats(Q, s, False) @ f
    worst = 0.0
    for q in (Q, Q + 1, Q + 2):
        mod = 3 ** q
        fine = 3 ** (q + 3)          # representatives of a + 3^q Z_3 modulo 3^(q+3)
        for a in range(mod):
            reps = [a + mod * k for k in range(fine // mod)]
            # A/2: image set of x -> x/2 (2-adic inverse taken modulo 3^(q+3)); must be one level-q cylinder
            img = sorted(set((x * pow(2, -1, fine)) % mod for x in reps))
            assert len(img) == 1
            beta_half = f[img[0] % 3 ** Q] / mod
            # (3A+1)/2: image of x -> (3x+1)/2 modulo 3^(q+1) must be one level-(q+1) cylinder
            img2 = sorted(set(((3 * x + 1) * pow(2, -1, 3 * fine)) % (3 * mod) for x in reps))
            assert len(img2) == 1
            beta_odd = f[img2[0] % 3 ** Q] / (3 * mod)
            lhs = -(2.0 ** (-s)) * beta_half - (1.5 ** s) * beta_odd
            rhs = PQ[a % 3 ** Q] / mod
            worst = max(worst, abs(lhs - rhs))
    print(f'  Q = {Q}: max |(P beta)(A) - 3^-q (P_Q f)(a)| over cylinders of levels Q..Q+2: {worst:.2e}')

print('(b) spectra versus statement 7')
S_LIST = [0.7 + 3.1j, 0.3 - 11.0j, 1.2 + 0.4j, 0.5 + 14.134725j, 2.0 + 1.0j, 0.05 + 40.0j,
          0.99820922076525681 + 240.24783014069124j, 0.69391836 + 22.2532j]
for Q in range(1, 6):
    m = 3 ** Q
    U = meanzero_basis(m)
    worst_full = worst_ex = worst_inv = worst_cp = 0.0
    for s in S_LIST:
        P = mats(Q, s, False)
        Pt = mats(Q, s, True)
        evf = np.linalg.eigvals(P)
        worst_full = max(worst_full, match_err(evf, predicted(Q, s, True)))
        R = U.conj().T @ Pt @ U
        worst_inv = max(worst_inv, np.abs((np.eye(m) - U @ U.conj().T) @ Pt @ U).max())
        eve = np.linalg.eigvals(R)
        worst_ex = max(worst_ex, match_err(eve, predicted(Q, s, False)))
        # characteristic polynomial at random points (relative error of log|det|)
        for z in (0.3 + 0.2j, -0.7 + 1.1j, 1.5 - 0.4j):
            pr = predicted(Q, s, True)
            lhs = np.linalg.slogdet(z * np.eye(m) - P)
            val = np.sum(np.log(np.abs(z - pr)))
            ang = np.angle(np.prod((z - pr) / np.abs(z - pr)))
            worst_cp = max(worst_cp, abs(lhs[1] - val) / max(1.0, abs(val)), abs(np.angle(lhs[0] / np.exp(1j * ang))))
    print(f'  Q = {Q} (dim {m}): max matched eigenvalue error full P {worst_full:.1e}; P~ on mean-zero {worst_ex:.1e}; '
          f'||(I - Pi) P~ Pi|| {worst_inv:.1e}; char. poly (rel. log|det| / phase) {worst_cp:.1e}')

print('    exact characteristic polynomial check with mpmath (60 digits), Q <= 3:')
mp.mp.dps = 60
for Q in (1, 2, 3):
    m = 3 ** Q
    i2 = pow(2, -1, m)
    worst = mp.mpf(0)
    for s in (mp.mpc('0.7', '3.1'), mp.mpc('1.2', '0.4'), mp.mpc('0.99820922076525681364', '240.24783014069124432530')):
        A = mp.matrix(m, m)
        a2 = -mp.power(2, -s)
        b3 = -mp.power(mp.mpf(3) / 2, s) / 3
        for a in range(m):
            A[a, (a * i2) % m] += a2
            A[a, ((3 * a + 1) * i2) % m] += b3
        Dv = 1 + mp.power(2, -s) + mp.power(mp.mpf(3) / 2, s) / 3
        for z in (mp.mpc('0.3', '0.2'), mp.mpc('-0.7', '1.1')):
            lhs = mp.det(z * mp.eye(m) - A)
            rhs = (z - (1 - Dv))
            for k in range(Q):
                L = 2 * 3 ** k
                rhs *= (z ** L - a2 ** L)       # prod over zeta^L = 1 of (z + 2^-s zeta) = z^L - (-2^-s)^L
            worst = max(worst, abs(lhs - rhs) / abs(rhs))
    print(f'    Q = {Q}: max relative |det(zI - P) - (z - 1 + D) prod_k (z^(2 3^k) - (-2^-s)^(2 3^k))| = {mp.nstr(worst, 3)}')

print('(c) geometric multiplicities at the repeated eigenvalues (algebraic multiplicity Q each)')
for Q in range(1, 6):
    m = 3 ** Q
    U = meanzero_basis(m)
    out = []
    for s in (0.7 + 3.1j, 1.2 + 0.4j):
        Pt = mats(Q, s, True)
        R = U.conj().T @ Pt @ U
        for zeta in (1, -1):
            lam = -(2.0 ** (-s)) * zeta
            sv = np.linalg.svd(R - lam * np.eye(m - 1), compute_uv=False)
            nul = int((sv < 1e-9 * sv[0]).sum())
            out.append(f's={s}: zeta={zeta:+d} geometric {nul} (smallest sv {sv[-1]:.1e}, next {sv[-2] if m > 2 else float("nan"):.1e})')
    print(f'  Q = {Q}: ' + ' | '.join(out))

print('(d) eigenvalue 1 at the Mercator zero; mean-zero invariance')
s0 = 0.99820922076525681 + 240.24783014069124j
print(f'  D(s0) = {Dfun(s0):.2e}')
for Q in range(1, 6):
    m = 3 ** Q
    P = mats(Q, s0, False)
    Pt = mats(Q, s0, True)
    one = np.ones(m)
    U = meanzero_basis(m)
    ev = np.linalg.eigvals(U.conj().T @ Pt @ U)
    print(f'  Q = {Q}: |P 1 - 1| = {np.abs(P @ one - one).max():.1e}; |1^T P~ + 2^-s 1^T| = '
          f'{np.abs(one @ Pt + (2.0 ** (-s0)) * one).max():.1e}; min |lambda - 1| over spec(P~|mean-zero) = '
          f'{np.abs(ev - 1).min():.4f} (1 - 2^-Re s0 = {1 - 2 ** (-s0.real):.4f})')
