#!/usr/bin/env python3
"""procgen_cauchy_20260925 -- part 1: exact second moments of the weighted backward tree, the Cauchy-Schwarz
lower bound for the provability price, the moment-method limit, SHEET and DRIFT controls.

Objects (all exact integer arithmetic unless stated):
  Z_k(v)  = sum over backward T-paths of length k from v of 3^a / 2^k   (a = number of odd steps),
            a function of v mod 3^k;  I_k = 2^k Z_k is an integer table on Z/3^k:
            I_k(v) = I_{k-1}(2v) + 3 [v = 2 mod 3] I_{k-1}((2v-1)/3).
  W_L     = sum_{k<L} Z_k, a function of v mod 3^(L-1)  (THM-4475's hub weight, before the max).
  M2(L)   = E_Haar[W_L^2] = sum_{j,k<L} <Z_j, Z_k>.
  gamma_k = int Z_k(3u+1) Z_k(u) du  (cross term);  identity  ||Z_{k+1}||^2 = ||Z_k||^2 + gamma_k / 2.
Output: plain text on stdout.  Peak memory about 420 MB (K = 15 tables).
"""
import math, sys, time
import numpy as np
from fractions import Fraction as Fr

K = 15                     # tables Z_0..Z_K, exact Gram matrix, M2(L) for L <= K+1
LOG32 = math.log(2) / math.log(3)
H = -(LOG32 * math.log2(LOG32) + (1 - LOG32) * math.log2(1 - LOG32))
ETA = 1 - H               # 0.050044

def thr(k, q=3):
    a = 0
    while q ** a <= 2 ** k:
        a += 1
    return a

def bad_count(L, q=3):
    T = [thr(k, q) for k in range(L + 1)]
    dist = {0: 1}
    for k in range(1, L + 1):
        nd = {}
        for a, c in dist.items():
            for up in (0, 1):
                b = a + up
                if b >= T[k]:
                    nd[b] = nd.get(b, 0) + c
        dist = nd
    return sum(dist.values())

def exact_dot(a, b):
    """exact integer dot product of two integer arrays (chunks cast to int64; no partial sum can overflow)."""
    ma = max(int(a.max()), -int(a.min())) if len(a) else 0
    mb = max(int(b.max()), -int(b.min())) if len(b) else 0
    ch = max(1, min(1 << 16, (9 * 10 ** 18) // max(1, ma * mb)))
    s = 0
    for i in range(0, len(a), ch):
        s += int(np.dot(a[i:i + ch].astype(np.int64), b[i:i + ch].astype(np.int64)))
    return s

def build_tables(K, q=3, legal=2, sheet=+1):
    """I_k on Z/q^k for k <= K.  sheet=+1: E(v)=(2v-1)/q legal iff v = legal mod q (3n+1 type);
    sheet=-1 (q=3 only): E(v)=(2v+1)/3 legal iff v = 1 mod 3 (the 3n-1 map)."""
    # int32 storage: I_k <= 2^k max Z_k < 2^31 for the sizes used here (q = 3: K <= 15; q = 5: K <= 9)
    assert (q == 3 and K <= 15) or (q == 5 and K <= 9)
    tabs = [np.ones(1, dtype=np.int32)]
    for k in range(1, K + 1):
        m = q ** k; mp = q ** (k - 1)
        prev = tabs[-1]
        v = np.arange(m, dtype=np.int32 if 2 * m < 2 ** 31 else np.int64)
        v *= 2; v %= mp                                   # index 2v mod q^(k-1), in place
        cur = prev[v]
        del v
        t = np.arange(mp, dtype=np.int32 if 2 * m < 2 ** 31 else np.int64)
        # sheet +1: v = q t + legal, E(v) = 2t + (2 legal - 1)/q ;  sheet -1: v = 3t + 1, (2v+1)/3 = 2t + 1
        off = (2 * legal - 1) // q if sheet == +1 else 1
        t *= 2; t += off; t %= mp
        pos = legal if sheet == +1 else 1
        cur[pos::q] += q * prev[t]
        del t
        tabs.append(cur)
    return tabs

def gram(tabs, q=3):
    Kk = len(tabs) - 1
    G = [[None] * (Kk + 1) for _ in range(Kk + 1)]
    for k in range(Kk + 1):
        f = tabs[k]
        for j in range(k, -1, -1):
            if j < k:
                f = f.reshape(q, q ** j).sum(axis=0)     # fold one level
            s = exact_dot(tabs[j], f)
            G[j][k] = G[k][j] = Fr(s, 2 ** (j + k) * q ** k)
    return G

def main():
    t0 = time.time()
    print('=' * 100)
    print('PART 1. Exact second moments of the weighted backward tree (3-adic, exact integers)')
    print('=' * 100)
    tabs = build_tables(K)
    print('tables I_0..I_%d built (%.1fs); max_v Z_k = I_k/2^k and max/1.5^k:' % (K, time.time() - t0))
    for k in (5, 10, 15):
        mx = int(tabs[k].max()); print('   k=%2d  max Z_k = %.3f   max Z_k / 1.5^k = %.4f' % (k, mx / 2 ** k, mx / 2 ** k / 1.5 ** k))
    # mean check: E[Z_k] = 1 exactly
    for k in range(K + 1):
        assert int(tabs[k].sum()) == 2 ** k * 3 ** k, 'E[Z_k] != 1'
    print('check: E_Haar[Z_k] = 1 exactly for k <= %d  (g(1) = 1: the mean-1 martingale)' % K)
    G = gram(tabs)
    print('Gram matrix <Z_j,Z_k> exact for j,k <= %d (%.1fs)' % (K, time.time() - t0))
    # cross terms and the criticality identity
    print('\n(1a) ||Z_k||^2, the cross term gamma_k = int Z_k(3u+1)Z_k(u)du, the identity ||Z_{k+1}||^2 = ||Z_k||^2 + gamma_k/2,')
    print('     and the i.i.d.-cascade value 1 + k/2 (independent siblings):')
    print('   k   ||Z_k||^2            increment    gamma_k        identity   iid 1+k/2   <Z_{k-1},Z_k>/||Z_{k-1}||^2')
    gam = []
    for k in range(K + 1):
        m = 3 ** k
        u = np.arange(m, dtype=np.int32 if 3 * m < 2 ** 31 else np.int64)
        u *= 3; u += 1; u %= m
        g = Fr(exact_dot(tabs[k][u], tabs[k]), 4 ** k * m)
        del u
        gam.append(g)
    for k in range(K + 1):
        inc = G[k][k] - G[k - 1][k - 1] if k else Fr(0)
        ok = (G[k][k] == G[k - 1][k - 1] + gam[k - 1] / 2) if k else True
        mart = float(G[k - 1][k] / G[k - 1][k - 1]) if k else 1.0
        print('  %2d   %-18.12f  %-11.6f  %-13.9f  %-9s  %-10.2f  %.6f' % (k, float(G[k][k]), float(inc), float(gam[k]), ok, 1 + k / 2, mart))
    assert all(G[k][k] == G[k - 1][k - 1] + gam[k - 1] / 2 for k in range(1, K + 1))
    print('   identity verified exactly (rationals) for k = 1..%d.' % K)
    print('   exact values: ||Z_2||^2 = %s, ||Z_3||^2 = %s, ||Z_4||^2 = %s, gamma_3 = %s' % (G[2][2], G[3][3], G[4][4], gam[3]))
    # ladder expansion of the cross term: Z_k(3x+1) = sum_{j=1}^{J} (3/4) 4^(1-j) Z_{k-2j}(S^j x) + 4^(-J) Z_{k-2J}(3 S^J x + 1),
    # S x = 4x + 1, J = floor(k/2), Z_0(.) = 1, Z_1(3y+1) = 1/2;  hence
    # gamma_k = 3 sum_j 4^-j <Z_{k-2j} o S^j, Z_k> + 4^-J c_k,  c_k = 1 (k even), 1/2 (k odd)   [E Z_k = 1]
    lad = []
    for k in range(1, 13):
        m = 3 ** k; x = np.arange(m, dtype=np.int64); J = k // 2
        tot = Fr(0)
        for j in range(1, J + 1):
            a = k - 2 * j; ma = 3 ** a
            Sj = (pow(4, j, m) * x + ((4 ** j - 1) // 3) % m) % m
            tot += 3 * Fr(1, 4 ** j) * Fr(exact_dot(tabs[a][Sj % ma], tabs[k]), 2 ** (a + k) * m)
        tot += Fr(1, 4 ** J) * (Fr(1) if k % 2 == 0 else Fr(1, 2))
        lad.append(tot == gam[k])
    print('   ladder expansion gamma_k = 3 sum_j 4^-j <Z_{k-2j} o S^j, Z_k> + 4^-J c_k verified exactly for k = 1..12: %s' % all(lad))
    assert all(lad)
    print('   increments ~0.357 and still slowly decreasing (gamma_k ~0.713); iid siblings would give 1/2 (gamma = 1).  See PART 2.')

    print('\n(1b) M2(L) = E_Haar[W_L^2] exactly (rational), L <= %d:' % (K + 1))
    rho = {}
    M2 = {}
    print('   L   M2(L)              M2/L^3    iid-model M2    rho_L        CS bound rho^2/M2   THM-4475 rho/sum(g)   CS/THM4475')
    gseq = [1.0, 2.0]
    while len(gseq) < 60:
        gseq.append(1.5 * gseq[-1] + 0.25 * gseq[-2])
    for L in range(1, K + 2):
        m2 = sum(G[j][k] for j in range(L) for k in range(L))
        M2[L] = m2
        iid = sum(1 + min(j, k) / 2 for j in range(L) for k in range(L))
        rho[L] = Fr(bad_count(L), 2 ** L)
        cs = rho[L] ** 2 / m2
        old = float(rho[L]) / sum(gseq[:L])
        print('  %2d   %-17.6f  %-8.5f  %-14.2f  %-11.6f  %-18.4e  %-20.4e  %.3f' % (L, float(m2), float(m2) / L ** 3, iid, float(rho[L]), float(cs), old, float(cs) / old))
    print('   exact: M2(8) = %s,  M2(12) = %s' % (M2[8], M2[12]))
    print('   (at small L the sup bound of THM-4475 is still better; the CS bound wins once rho_L^2 * sum(g) < rho_L * M2,')
    print('    i.e. asymptotically, since M2 grows polynomially in this range and sum(g) ~ lambda^L.)')
    # asymptotic exponents
    lam = (3 + math.sqrt(13)) / 4
    print('\n   exponents (base 2, per L):  THM-4475 sup bound: eta + log2(lambda) = %.4f;  CS with poly M2: 2 eta = %.4f;  sharp (HYP-9137): eta = %.4f'
          % (ETA + math.log2(lam), 2 * ETA, ETA))
    print('   CS/THM4475 is still < 1 for L <= %d (it rises from 0.06 at L = 8 to 0.24 at L = 16); the rigorous majorant of PART 3' % (K + 1))
    print('   puts the crossover of the proved bounds near L = 22.')

    # (1c) Holder moments
    print('\n(1c) Holder exponents: E[Z_k^p] >= g(p)^k exactly (single-path lower bound), g(p) = 2^-p (1 + 3^(p-1)).')
    def gfun(p):
        return 2.0 ** (-p) * (1 + 3.0 ** (p - 1))
    print('   p      g(p)       E[Z_%d^p]    g(p)^%d     ratio' % (K, K))
    for p in (1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0):
        z = tabs[K].astype(np.float64) / 2 ** K
        ep = float(np.mean(z ** p))
        print('   %-5.2f  %-9.5f  %-11.4f  %-10.4f  %.3f' % (p, gfun(p), ep, gfun(p) ** K, ep / gfun(p) ** K))
    print('   Holder bound delta >= rho^(p/(p-1)) / (E W^p)^(1/(p-1)); its exponent is')
    print('   e(p) = [p*eta + max(0, log2 g(p))]/(p-1)  (lower bound, using E W_L^p >= g(p)^(L-1) and E W^p >= 1):')
    for p in (1.1, 1.5, 1.9, 2.0, 2.1, 2.5, 3.0, 4.0):
        e = (p * ETA + max(0.0, math.log2(gfun(p)))) / (p - 1)
        print('   p = %-4.2f   e(p) >= %.5f' % (p, e))
    slope = math.log(27 / 16) / (4 * math.log(2))
    print('   (log2 g)\'(2) = ln(27/16)/(4 ln 2) = %.5f > eta = %.5f, so e(p) > 2 eta for every p != 2 (convexity of log g).' % (slope, ETA))
    del z
    print('   trend of E[Z_k^p] in k (bounded for p < 2 in the iid model; here FINITE-EXACT):')
    print('   k    ' + '  '.join('p=%-6.2f' % p for p in (1.25, 1.5, 1.75, 2.0)))
    for k in range(3, K + 1, 2):
        z = tabs[k].astype(np.float64) / 2 ** k
        print('   %2d   ' % k + '  '.join('%-8.4f' % float(np.mean(z ** p)) for p in (1.25, 1.5, 1.75, 2.0)))
        del z

    # (1d)-(1e): free the big tables first (keep Z_14 for the tail display)
    z14 = np.sort(tabs[14].astype(np.float64) / 2 ** 14)       # ascending, 38 MB
    del tabs, G
    print('\n(1d) the distribution-only optimum  delta*_L = min{ mu(B) : int_B W_L >= rho_L }  (best bound from the law of W_L alone)')
    print('   L   rho_L       delta*_L      rho^2/M2      delta*/(rho^2/M2)   rho/max W     max W_L    max W_L/1.5^L')
    J = np.ones(1, dtype=np.int32)          # J_L = 2^(L-1) W_L on Z/3^(L-1); W_1 = 1  (J_16 < 7.3e7 fits int32)
    tail_lines = []
    for L in range(2, K + 2):
        # W_L = 1 + L W_{L-1}:  J_L(v) = 2^(L-1) + J_{L-1}(2v) + 3[v=2] J_{L-1}(E v)
        m = 3 ** (L - 1); mp = 3 ** (L - 2)
        v = np.arange(m, dtype=np.int32); v *= 2; v %= mp
        cur = J[v]; del v
        cur += 2 ** (L - 1)
        t = np.arange(mp, dtype=np.int32); t *= 2; t += 1; t %= mp
        cur[2::3] += 3 * J[t]
        del t
        J = cur
        if L in M2:
            assert Fr(exact_dot(J, J), 4 ** (L - 1) * m) == M2[L]
        if L >= 6:
            w = J.astype(np.float64); w /= 2 ** (L - 1); w.sort(); w = w[::-1]      # descending view
            cs = np.cumsum(w); cs /= m
            idx = int(np.searchsorted(cs, float(rho[L])))
            dstar = (idx + 1) / m
            print('  %2d   %-10.6f  %-12.4e  %-12.4e  %-18.2f  %-12.4e  %-9.2f  %.3f' % (L, float(rho[L]), dstar, float(rho[L] ** 2 / M2[L]),
                  dstar / float(rho[L] ** 2 / M2[L]), float(rho[L]) / w[0], w[0], w[0] / 1.5 ** L))
            if L == K + 1:
                asc = w[::-1]                       # ascending view of the same buffer
                for x in (2, 4, 8, 16, 32, 64, 128, 256, 512, 1024):
                    iw = m - int(np.searchsorted(asc, x, side='right'))      # number of W > x
                    ew = cs[iw - 1] if iw else 0.0                              # E[W; W > x]
                    iz = len(z14) - int(np.searchsorted(z14, x, side='right'))
                    ez = float(z14[len(z14) - iz:].sum()) / len(z14) if iz else 0.0
                    pw = iw / m; pz = iz / len(z14)
                    tail_lines.append('   %-7d  %-10.3e  %-11.4f  %-11.4f |  %-10.3e  %-11.4f  %.4f' % (x, pw, x * x * pw, x * ew, pz, x * x * pz, x * ez))
            del w, cs
    print('   check: M2(L) recomputed from the W_L recursion W_L = 1 + L W_{L-1} agrees exactly.')
    del J

    print('\n(1e) tails at L = %d (W_L) and k = 14 (Z_k) under Haar: x^2 P(. > x) and x * E[. 1{. > x}]' % (K + 1))
    print('   x        P(W>x)      x^2 P(W>x)   x E[W;W>x]  |  P(Z>x)      x^2 P(Z>x)   x E[Z;Z>x]')
    for line in tail_lines:
        print(line)
    print('   (a tail ~ c x^-2 shows as a plateau of x^2 P and x E[.;.>x]; finite k, L cut the tail off near max ~ 1.7*1.5^k.)')
    del z14

    # (1f) SHEET
    print('\n(1f) SHEET: the 3n-1 map U(n) = (3n-1)/2 has backward E-child (2v+1)/3 (v = 1 mod 3).')
    tp = build_tables(12)
    tm = build_tables(12, sheet=-1)
    ok = all(np.array_equal(tm[k], tp[k][(-np.arange(3 ** k, dtype=np.int64)) % (3 ** k)]) for k in range(13))
    print('   Z^-_k(v) = Z_k(-v) for all v mod 3^k, k <= 12: %s  => M2, rho_L, the CS bound and the harmonic bound are identical on both sheets.' % ok)
    assert ok
    del tp, tm

    # (1g) DRIFT
    print('\n(1g) DRIFT: 5x+1.  Z^(5)_k on Z/5^k (E-child (2v-1)/5, legal iff v = 3 mod 5, weight 5/2).')
    t5 = build_tables(9, q=5, legal=3)
    for k in range(10):
        assert int(t5[k].sum()) == 2 ** k * 5 ** k
    print('   E[Z^(5)_k] = 1 exactly (g_5(1) = 1).')
    print('   k   ||Z^(5)_k||^2     (3/2)^k        ratio      cross term check ||Z_{k+1}||^2 - 1.5||Z_k||^2 = <Z_k o tau_5, Z_k>/2 ?')
    n5 = []
    for k in range(10):
        n5.append(Fr(exact_dot(t5[k], t5[k]), 4 ** k * 5 ** k))
    for k in range(10):
        chk = ''
        if k < 9:
            m = 5 ** k; u = np.arange(m, dtype=np.int64)
            cr = Fr(exact_dot(t5[k][(5 * u + 1) % m], t5[k]), 4 ** k * m)
            chk = str(n5[k + 1] == Fr(3, 2) * n5[k] + cr / 2)
        print('  %2d   %-15.4f  %-13.4f  %-9.4f  %s' % (k, float(n5[k]), 1.5 ** k, float(n5[k]) / 1.5 ** k, chk))
    print('   => ||Z^(5)_k||^2 >= (3/2)^k: the diagonal coefficient is g_5(2) = (1+5)/4 = 3/2 (AM of 1/2 and 5/2).')
    # 5x+1 undecided density and the CS bound
    print('   5x+1 undecided density beta_L (5^a_k > 2^k for all k <= L) and the CS bound beta^2/M2^(5) <= beta^2 (2/3)^(L-1):')
    G5 = gram(t5, q=5)
    for L in (4, 6, 8, 10):
        m25 = sum(G5[j][k] for j in range(L) for k in range(L))
        beta = bad_count(L, q=5) / 2 ** L
        print('   L=%2d beta_L=%.4f  M2^(5)(L)=%.2f  CS=%.3e  beta^2 (2/3)^(L-1)=%.3e' % (L, beta, float(m25), beta ** 2 / float(m25), beta ** 2 * (2 / 3) ** (L - 1)))
    for L in (20, 40, 100, 400):
        beta = bad_count(L, q=5) / 2 ** L
        print('   L=%3d beta_L=%.4f   CS bound <= %.3e' % (L, beta, beta ** 2 * (2 / 3) ** (L - 1)))
    print('\nPART 1 done in %.1fs' % (time.time() - t0))

if __name__ == '__main__':
    main()
