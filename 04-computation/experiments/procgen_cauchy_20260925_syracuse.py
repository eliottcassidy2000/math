#!/usr/bin/env python3
"""procgen_cauchy_20260925 -- part 2: the fixed point of the transfer operator (the Syracuse law pi on Z_3),
its L^2 level energies, the Poisson-kernel identity, ladder autocorrelations, and the tail of its density.

The weighted backward-tree counts are Z_k = L^k 1, where L f(v) = f(2v)/2 + (3/2)[v=2 mod 3] f((2v-1)/3) is the
Perron-Frobenius operator of the forward 3-adic walk X -> X/2 or (3X+1)/2 (prob 1/2 each).  Its adjoint is a
Wasserstein-1 contraction (constant 2/3) in the 3-adic metric; the fixed point is the Syracuse law pi.
pi mod 3^n is computed as the stationary law of the projected chain on Z/3^n, started from the lifted law mod
3^(n-1): after t steps only all-even step words carry the lifting error, so the total-variation error is <= 2^-t
(t = 64 here).  Output: stdout.  Peak memory about 400 MB (n = 15).
"""
import math, time
import numpy as np

NMAX = 15
NLAD = 14

def maps(N):
    inv2 = pow(2, -1, N)
    x = np.arange(N, dtype=np.int64)
    ev = ((x * inv2) % N).astype(np.int32)
    od = (((3 * x + 1) % N * inv2) % N).astype(np.int32)
    del x
    return ev, od

def main():
    t0 = time.time()
    print('=' * 100)
    print('PART 2. The fixed point: the Syracuse law pi (stationary law of the 3-adic walk), L^2 energies, tail')
    print('=' * 100)
    p = np.array([1.0])
    prev_col = 1.0
    E = {}
    col = {}
    rows = []
    for n in range(1, NMAX + 1):
        N = 3 ** n; Nb = 3 ** (n - 1)
        # stationary law mod 3^n: pi_n(v) = pi_n(2v)/2 + [v = 2 mod 3] pi_{n-1}(E v)/2, E v = (2v-1)/3 mod 3^(n-1);
        # the second term is fixed (the marginal mod 3^(n-1) is already stationary), so iterate the first.
        p_prev = p
        twov = np.arange(N, dtype=np.int32 if 2 * N < 2 ** 31 else np.int64); twov *= 2; twov %= N
        t = np.arange(Nb, dtype=np.int32 if 2 * N < 2 ** 31 else np.int64); t *= 2; t += 1; t %= Nb
        const = np.zeros(N); const[2::3] = 0.5 * p_prev[t]; del t
        p = np.tile(p_prev, 3) / 3.0
        for it in range(64):
            q = p[twov]; q *= 0.5; q += const
            p = q
        del twov, const, q
        l1 = 0.0
        pr = p.reshape(3, Nb)
        for row in range(3):
            l1 += float(np.abs(pr[row] - p_prev / 3.0).sum())       # || p_n - p_{n-1} ||_1 (densities)
        del pr
        c = N * float(np.dot(p, p))
        col[n] = c; E[n] = c - prev_col; prev_col = c
        m1 = N * float(np.dot(p[1::3], p[1::3])); m2 = N * float(np.dot(p[2::3], p[2::3])); m0 = N * float(np.dot(p[0::3], p[0::3]))
        imax = int(np.argmax(p)); pmax = float(p[imax])
        rows.append((n, c, E[n], l1))
        print('  n=%2d  ||pi||_n^2=%.12f  E_n=%.12f  ||p_n-p_{n-1}||_1=%.3e  L2-mass classes 1:2 = %.6f : %.6f  (class 0: %.1e)  max ball mass*2^n=%.6f at -1: %s  %.1fs'
              % (n, c, E[n], l1, m1 / c, m2 / c, m0, pmax * 2 ** n, imax == N - 1, time.time() - t0), flush=True)
        if n == NLAD:
            p_lad = p.copy()
        del p_prev
    print('\n(2a) exact low levels: E_1 = 2/3 (err %.1e), E_2 = 10/21 (err %.1e); ||pi||_2^2 = 15/7.' % (abs(E[1] - 2 / 3), abs(E[2] - 10 / 21)))
    print('     L^2 mass on classes 1 and 2 mod 3 is exactly 1/5 : 4/5 at every level (|kappa_+^|^2 : |kappa_-^|^2 = 1 : 4 pointwise).')
    print('     max_y pi(y + 3^n Z_3) = pi(-1 + 3^n Z_3) ~ c 2^-n: the L^infinity dimension is log_3 2 = 0.631 > 1/2.')
    print('     level energies E_n rise slowly towards ~0.476 (increments shrink like ~1/n or faster);  ||pi||_n^2 ~ 0.47 n.')
    print('     L^1 increments ||p_n - p_{n-1}||_1 shrink slowly (last ratio %.3f); consistent with L^1 convergence (absolute continuity of pi), not a proof.'
          % (rows[-1][3] / rows[-2][3]))

    # tail of the fixed-point density at the finest level
    N = 3 ** NMAX
    dens = p; del p
    dens *= N; dens.sort()
    print('\n(2b) tail of the fixed-point density p_%d = d pi/d mu (at resolution 3^-%d) under Haar:  x^2 P(p > x)' % (NMAX, NMAX))
    for x in (2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128):
        cnt = N - int(np.searchsorted(dens, x, side='right'))
        print('   x=%-4d  P(p>x)=%.4e   x^2 P=%.4f   x E[p;p>x]=%.4f' % (x, cnt / N, x * x * cnt / N, x * float(dens[N - cnt:].sum()) / N if cnt else 0.0))
    del dens

    # ladder autocorrelations and the Poisson-kernel identity
    print('\n(2c) ladder autocorrelations C_j^(n) = 3^n sum_y pi_n(y) pi_n(S^j y), S(y) = 4y+1 (translation by j in the')
    print('     coordinates Lambda = log_4(1+3y)), and the identity ||pi||_{n+1}^2 = sum_j 4^-|j| C_j^(n):')
    pl = p_lad
    Nn = 3 ** NLAD
    y = np.arange(Nn, dtype=np.int64)
    C = [Nn * float(np.dot(pl, pl))]
    for j in range(1, 41):
        Sj = (pow(4, j, Nn) * y + ((4 ** j - 1) // 3) % Nn) % Nn
        C.append(Nn * float(np.dot(pl, pl[Sj])))
        del Sj
    pred = C[0] + 2 * sum(4.0 ** (-j) * C[j] for j in range(1, 41))
    print('   n=%d: C_1..C_12 = %s' % (NLAD, ' '.join('%.4f' % C[j] for j in range(1, 13))))
    print('   C_27 = %.4f (v_3 = 3),  C_9 = %.4f (v_3 = 2),  C_3 = %.4f (v_3 = 1): C_j grows with v_3(j) (log singularity at 0).' % (C[27], C[9], C[3]))
    print('   sum_j 4^-|j| C_j^(%d) = %.12f   vs   ||pi||_%d^2 = %.12f   (difference %.1e; tail |j|>40 < 1e-23)' % (NLAD, pred, NLAD + 1, col[NLAD + 1], abs(pred - col[NLAD + 1])))
    print('   so E_{n+1} = 2 sum_{j>=1} 4^-j C_j^(n): the level energy is bounded iff the off-diagonal ladder correlations are.')
    del y

    # Fourier check of the Poisson-kernel identity for n <= 11 (Lambda coordinates)
    print('\n(2d) Fourier form: in Lambda = log_4(1+3y) coordinates,  ||pi||_n^2 = sum_{lev(xi) <= n-1} P(xi) |nu^(xi)|^2,')
    print('     P(xi) = 15/(17 - 8 cos 2 pi {xi}) = Poisson kernel P_{1/4}, whose average over any full level is (1+4^-N)/(1-4^-N) ~ 1:')
    pp = np.array([1.0])
    for n in range(1, 12):
        Nn = 3 ** n
        pp = np.tile(pp, 3) / 3.0
        ev, od = maps(Nn)
        for it in range(64):
            q = np.bincount(ev, weights=pp, minlength=Nn); q += np.bincount(od, weights=pp, minlength=Nn); q *= 0.5; pp = q
        del ev, od
        # nu(t) = pi(ell^{-1}(t)),  ell^{-1}(t) = (4^t - 1)/3 mod 3^n
        M = 3 ** (n + 1)
        four = np.empty(Nn, dtype=np.int64); acc = 1
        for t in range(Nn):
            four[t] = acc; acc = acc * 4 % M
        yy = ((four - 1) // 3) % Nn
        assert len(np.unique(yy)) == Nn
        nu = pp[yy]
        F = np.fft.fft(nu)
        xi = np.arange(Nn)
        P = 15.0 / (17.0 - 8.0 * np.cos(2 * np.pi * xi / Nn))
        mask = (xi % 3 == 0)
        lhs = Nn * float(np.dot(pp, pp))
        rhs = float(np.sum(P[mask] * np.abs(F[mask]) ** 2))
        # tilt of each level: sum_{lev=m} (P-1)|nu^|^2 / E_m
        print('   n=%2d  ||pi||_n^2 = %.12f   sum P|nu^|^2 = %.12f   diff %.1e' % (n, lhs, rhs, abs(lhs - rhs)))
    print('\nPART 2 done in %.1fs' % (time.time() - t0))

if __name__ == '__main__':
    main()
