#!/usr/bin/env python3
"""procgen_cauchy_20260925 -- part 4: the harmonic (pointwise-covering) lower bound, which uses which bad orbits
meet which hubs, and beats every bound that depends on the law of W_L alone.

    delta_L >= H_L := E[ 1_Bad(omega) / max_{j<L} W_L(X_j) ]
where omega is a uniform parity word of length L, X_0 is Haar on Z_3 (independent: CRT), X_{j+1} = X_j/2 or
(3X_j+1)/2 according to omega_j; X_j is exactly T^j(n) mod 3^(L-1) for n in the class (omega, X_0).
Exact for L <= 13 (DFS over bad prefixes, all X_0 mod 3^(L-1)); Monte Carlo for L = 14..16.
Also: quenched vs annealed conditional hub weights along bad orbits.  Output: stdout.  Peak memory about 300 MB.
"""
import math, sys, time
import numpy as np

LEX = 13          # exact up to here
LMC = 16          # Monte Carlo up to here
LOG32 = math.log(2) / math.log(3)

def thr(k):
    a = 0
    while 3 ** a <= 2 ** k:
        a += 1
    return a

def W_table(L, dtype=np.float64, want_m2=False):
    """W_L on Z/3^(L-1) via J_L = 2^(L-1) W_L,  J_L(v) = 2^(L-1) + J_{L-1}(2v) + 3[v=2] J_{L-1}((2v-1)/3).
    J is stored in int32 (J_16 < 7.3e7).  With want_m2, also return M2(L) = E[W_L^2] (float64, chunked)."""
    assert L <= 16
    J = np.ones(1, dtype=np.int32)
    for l in range(2, L + 1):
        m = 3 ** (l - 1); mp = 3 ** (l - 2)
        v = np.arange(m, dtype=np.int32); v *= 2; v %= mp
        cur = J[v]; del v
        cur += 2 ** (l - 1)
        t = np.arange(mp, dtype=np.int32); t *= 2; t += 1; t %= mp
        cur[2::3] += 3 * J[t]; del t
        J = cur
    m2 = 0.0
    ch = 1 << 20
    for i in range(0, len(J), ch):
        x = J[i:i + ch].astype(np.float64)
        m2 += float(np.dot(x, x))
    m2 /= len(J) * 4.0 ** (L - 1)
    W = J.astype(dtype); del J
    W /= 2 ** (L - 1)
    return (W, m2) if want_m2 else W

def exact(L):
    N = 3 ** (L - 1)
    W = W_table(L)
    x = np.arange(N, dtype=np.int64)
    inv2 = pow(2, -1, N)
    ev = ((x * inv2) % N).astype(np.int32)
    od = (((3 * x + 1) % N * inv2) % N).astype(np.int32)
    T = [thr(k) for k in range(L + 1)]
    S = math.log2(3)
    tot_inv = 0.0; nbad = 0
    quenched_log_mean = 0.0; annealed_mean = 0.0; quenched_mean_log = 0.0; rises = 0.0
    def rec(j, a, X, Mx, smin, rise, s):
        nonlocal tot_inv, nbad, quenched_log_mean, annealed_mean, quenched_mean_log, rises
        if j == L - 1:
            cnt = sum(1 for up in (0, 1) if a + up >= T[L])
            if cnt:
                inv = float(np.sum(1.0 / Mx)); mean = float(Mx.mean()); ml = float(np.log2(Mx).mean())
                nbad += cnt; tot_inv += cnt * inv
                quenched_log_mean += cnt * math.log2(mean); annealed_mean += cnt * mean
                quenched_mean_log += cnt * ml; rises += cnt * rise
            return
        for up in (0, 1):
            b = a + up
            if b >= T[j + 1]:
                Xn = (od if up else ev)[X]
                s2 = s + (S - 1 if up else -1)
                rec(j + 1, b, Xn, np.maximum(Mx, W[Xn]), min(smin, s2), max(rise, s2 - smin), s2)
    rec(0, 0, x.astype(np.int32), W.copy(), 0.0, 0.0, 0.0)
    rho = nbad / 2 ** L
    H = tot_inv / (2 ** L * N)
    EW2 = float(np.dot(W, W)) / N
    return dict(rho=rho, H=H, EW2=EW2, qlm=quenched_log_mean / nbad, ann=annealed_mean / nbad,
                qml=quenched_mean_log / nbad, rise=rises / nbad, maxW=float(W.max()))

def sample_bad_words(L, n, rng):
    T = [thr(k) for k in range(L + 1)]
    cnt = [dict() for _ in range(L + 1)]
    for a in range(L + 1):
        cnt[L][a] = 1
    for j in range(L - 1, -1, -1):
        for a in range(j + 1):
            cnt[j][a] = sum(cnt[j + 1].get(a + up, 0) for up in (0, 1) if a + up >= T[j + 1])
    words = np.zeros((n, L), dtype=np.int8)
    a = np.zeros(n, dtype=np.int64)
    for j in range(L):
        c1 = np.array([cnt[j + 1].get(b + 1, 0) if b + 1 >= T[j + 1] else 0 for b in range(L + 1)], dtype=np.float64)
        c0 = np.array([cnt[j + 1].get(b, 0) if b >= T[j + 1] else 0 for b in range(L + 1)], dtype=np.float64)
        p1 = c1[a] / (c1[a] + c0[a])
        up = (rng.random(n) < p1).astype(np.int8)
        words[:, j] = up; a += up
    return words, cnt[0][0]

def hat_W(W, L, words, X0):
    N = 3 ** (L - 1); inv2 = pow(2, -1, N)
    X = X0.copy(); M = W[X].astype(np.float64)
    for j in range(L - 1):
        up = words[:, j].astype(bool)
        X = np.where(up, ((3 * X + 1) % N * inv2) % N, (X * inv2) % N)
        M = np.maximum(M, W[X])
    return M

def main():
    t0 = time.time()
    print('=' * 100)
    print('PART 4. The harmonic (pointwise covering) bound H_L; quenched vs annealed hub weights on bad orbits')
    print('=' * 100)
    print('(4a) exact values (L <= %d): H_L against rho_L and the CS bound rho^2/M2 (for delta*_L see PART 1 (1d)); hub weights on bad orbits' % LEX)
    print('   L   rho_L      H_L          H/rho     H*M2/rho   CS=rho^2/M2   E_Bad[log2 What]  log2 E_Bad[E(What|w)]  E_Bad[log2 E(What|w)]  E_Bad[max rise]')
    for L in range(4, LEX + 1):
        d = exact(L)
        print('  %2d   %-9.6f  %-11.5e  %-8.5f  %-9.3f  %-12.4e  %-16.3f  %-21.3f  %-21.3f  %.3f' % (
            L, d['rho'], d['H'], d['H'] / d['rho'], d['H'] * d['EW2'] / d['rho'], d['rho'] ** 2 / d['EW2'],
            d['qml'], math.log2(d['ann']), d['qlm'], d['rise']), flush=True)
    print('   (annealed log2 E[What | Bad] grows faster than the quenched E_Bad log2 E[What | w]: rare bad words with')
    print('    heavy prefixes dominate the annealed mean; the harmonic bound only sees the typical, quenched behaviour.)')

    print('\n(4b) Monte Carlo (L = %d..%d; L = %d, %d cross-check (4a)), 400000 (word, X_0) samples; uniform words give the unconditioned comparison:' % (LEX + 1, LMC, LEX - 1, LEX))
    rng = np.random.default_rng(20260925)
    n = 400000
    print('   L   rho_L      H_L (+- s.e.)            H/rho     H*M2/rho   E[1/What|Bad]/E[1/What]   E_Bad log2 What - E log2 What   E[What|Bad]/E[What]')
    for L in range(LEX - 1, LMC + 1):
        W, EW2 = W_table(L, dtype=np.float32, want_m2=True)
        N = 3 ** (L - 1)
        wb, tot = sample_bad_words(L, n, rng)
        rho = tot / 2 ** L
        Mb = hat_W(W, L, wb, rng.integers(0, N, n))
        wu = (rng.random((n, L)) < 0.5).astype(np.int8)
        Mu = hat_W(W, L, wu, rng.integers(0, N, n))
        Hb = float(np.mean(1 / Mb)); se = float(np.std(1 / Mb)) / math.sqrt(n)
        print('  %2d   %-9.6f  %-11.5e (+-%.1e)    %-8.5f  %-9.3f  %-24.4f  %-30.3f  %.3f' % (
            L, rho, rho * Hb, rho * se, Hb, Hb * EW2, Hb / float(np.mean(1 / Mu)), float(np.mean(np.log2(Mb)) - np.mean(np.log2(Mu))),
            float(np.mean(Mb)) / float(np.mean(Mu))), flush=True)
        del W, Mb, Mu, wb, wu
    print('\nPART 4 done in %.1fs' % (time.time() - t0))

if __name__ == '__main__':
    main()
