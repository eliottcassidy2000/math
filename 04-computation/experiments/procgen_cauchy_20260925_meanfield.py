#!/usr/bin/env python3
"""procgen_cauchy_20260925 -- part 5: the mean-field (i.i.d.) smoothing transform W = W'/2 + (3/2) B W''
(B ~ Bernoulli(1/3), W', W'', B independent) against the exact 3-adic cascade; the moment function
g_q(s) = 2^-s + (q/2)^s / q = E_fwd[w^(s-1)] of the forward step factor w in {1/2, q/2};
Zolotarev contraction constants; strip confinement of bad parity words (for the mean-field harmonic bound).
Output: stdout.  Peak memory < 250 MB.
"""
import math, time
import numpy as np

def g(q, s):
    return 2.0 ** (-s) + (q / 2.0) ** s / q

def root(f, a, b):
    fa = f(a)
    for _ in range(200):
        m = (a + b) / 2
        if (f(m) > 0) == (fa > 0):
            a, fa = m, f(m)
        else:
            b = m
    return (a + b) / 2

def thr(k):
    a = 0
    while 3 ** a <= 2 ** k:
        a += 1
    return a

def main():
    t0 = time.time()
    print('=' * 100)
    print('PART 5. Mean-field smoothing transform vs the exact 3-adic cascade; moment functions; strips')
    print('=' * 100)
    print('(5a) g_q(s) = 2^-s (1 + q^(s-1)) = E_fwd[w^(s-1)], w in {1/2, q/2} with probability 1/2 each:')
    print('     g_q(1) = 1 (mean), g_q(2) = (1+q)/4 = AM of the forward factors, -g_q\'(1) = -E log w = log(2/sqrt q) (> 0: contracting):')
    print('   q   g_q(2)    -g_q\'(1) = log(4/q)/2   roots of g_q(s) = 1')
    for q in (1, 3, 5, 7, 9):
        g2 = (1 + q) / 4
        drift = 0.5 * math.log(4 / q)
        f = lambda s: g(q, s) - 1
        if q == 3:
            roots = '{1, 2}'
        elif q == 1:
            roots = '{1}'
        else:
            roots = '{%.4f, 1}' % root(f, 0.01, 0.999)
        print('   %d   %-8.4f  %-23.6f  %s' % (q, g2, drift, roots))
    print('   q = 3 is the only q with g_q(2) = 1: this is THM-4470\'s AM-fairness (q+1)/4 = 1, and makes kappa = 2.')
    d1 = 0.5 * math.log(4 / 3)
    print('   q = 3: -g\'(1) = log(2/sqrt 3) = %.6f (the drift), tree-density constant R = 1/log(2/sqrt 3) = %.5f (inverse-tree note, Prop. 10);'
          % (d1, 1 / d1))
    print('          g\'(2) = ln(27/16)/4 = %.6f > 0 (Kesten-Goldie-type slope at kappa = 2); the s=2 tilt has odd-step frequency 3/4.' % (math.log(27 / 16) / 4))
    print('\n(5b) Zolotarev zeta_s contraction constants of the mean-field smoothing transform S: zeta_s(S mu, S nu) <= g(s) zeta_s(mu, nu)')
    print('     (ideal metric of order s, 1 < s <= 2, laws with equal means):')
    print('     s:    ' + '  '.join('%-6.3f' % s for s in (1.1, 1.25, 1.5, 1.75, 1.9, 2.0)))
    print('     g(s): ' + '  '.join('%-6.4f' % g(3, s) for s in (1.1, 1.25, 1.5, 1.75, 1.9, 2.0)))
    print('     strict contraction exactly on 1 < s < 2; the constant is 1 at s = 2: the Banach fixed-point argument and the')
    print('     Cauchy-Schwarz lower bound fail at the same place, g(2) = 1.')

    # mean-field second moment: E Z_k^2 = 1 + k/2 exactly; population dynamics for tails and L^1 increments
    print('\n(5c) mean-field cascade: E[Z_k^2] = 1 + k/2 exactly (siblings independent); exact 3-adic values are smaller:')
    exact3 = [1.0, 1.5, 1.75, 2.25, 2.6640625, 3.12890625, 3.53564453125, 3.963134765625, 4.35760498046875, 4.739295959472656,
              5.111232757568359, 5.47853422164917, 5.840079426765442, 6.19878226518631, 6.556690484285355, 6.913230054080]
    print('   k:        ' + ' '.join('%-6d' % k for k in range(0, 16, 3)))
    print('   1+k/2:    ' + ' '.join('%-6.3f' % (1 + k / 2) for k in range(0, 16, 3)))
    print('   3-adic:   ' + ' '.join('%-6.3f' % exact3[k] for k in range(0, 16, 3)))
    rng = np.random.default_rng(7)
    n = 2_000_000
    cur = np.ones(n); prev = np.ones(n)            # pool of pairs (Z_k, Z_{k-1}) from the same tree
    l1 = []
    for k in range(1, 41):
        i1 = rng.integers(0, n, n); i2 = rng.integers(0, n, n)
        B = rng.random(n) < 1 / 3
        ncur = 0.5 * cur[i1] + np.where(B, 1.5 * cur[i2], 0.0)
        nprev = 0.5 * prev[i1] + np.where(B, 1.5 * prev[i2], 0.0)
        cur, prev = ncur, nprev
        if k == 1:
            prev = np.ones(n)
        l1.append(float(np.mean(np.abs(cur - prev))))
        del i1, i2, B, ncur, nprev
    zs = np.sort(cur)
    print('   population dynamics (pool %d, 40 generations): mean %.4f, E|Z_{k+1}-Z_k| at k=10,20,30,39: %s' %
          (n, float(cur.mean()), ', '.join('%.4f' % l1[k] for k in (10, 20, 30, 39))))
    print('   tail of the mean-field fixed point (k = 40):  x^2 P(W > x)')
    for x in (2, 4, 8, 16, 32, 64, 128):
        c = n - int(np.searchsorted(zs, x, side='right'))
        print('     x=%-4d  P=%.4e   x^2 P=%.4f' % (x, c / n, x * x * c / n))
    print('   (compare PART 2 (2b): the 3-adic fixed point p_15 has x^2 P(p > x) ~ 0.9 on [6, 32]; both are ~ c x^-2.)')
    del zs, cur, prev

    # L^1 increments of the exact 3-adic Z_k
    print('\n(5d) exact 3-adic L^1 increments ||Z_{k+1} - Z_k||_1 (does Z_k = L^k 1 converge in L^1 to d pi/d mu?):')
    I = np.ones(1, dtype=np.int64); vals = []
    for k in range(1, 15):
        m = 3 ** k; mp = 3 ** (k - 1)
        v = np.arange(m, dtype=np.int32); v *= 2; v %= mp
        cur = I[v]; del v
        t = np.arange(mp, dtype=np.int32); t *= 2; t += 1; t %= mp
        cur[2::3] += 3 * I[t]; del t
        prevl = np.tile(I, 3).astype(np.float64) / 2 ** (k - 1)
        vals.append(float(np.mean(np.abs(cur / 2 ** k - prevl))))
        I = cur
        del prevl
    print('   ' + ', '.join('k=%d: %.4f' % (k, vals[k - 1]) for k in range(2, 15, 2)))
    print('   successive ratios: ' + ', '.join('%.3f' % (vals[k] / vals[k - 1]) for k in range(8, 14)))
    del I, cur

    # strip confinement for bad words: P(Bad_L and max_k S_k <= y) / rho_L
    print('\n(5e) strip confinement of bad words (used by the mean-field harmonic bound): S_k = a_k log2 3 - k,')
    print('     q(L,y) = P(Bad_L, max_k S_k <= y) / P(Bad_L):')
    LS = math.log2(3)
    print('     L      rho_L        y=sqrt(L)   y=2sqrt(L)   y=0.1L')
    for L in (50, 100, 200, 400, 800):
        T = [thr(k) for k in range(L + 1)]
        res = []
        for y in (math.sqrt(L), 2 * math.sqrt(L), 0.1 * L, 1e9):
            dist = {0: 1.0}
            for k in range(1, L + 1):
                nd = {}
                for a, c in dist.items():
                    for up in (0, 1):
                        b = a + up
                        if b >= T[k] and b * LS - k <= y:
                            nd[b] = nd.get(b, 0.0) + c * 0.5
                dist = nd
            res.append(sum(dist.values()))
        rho = res[-1]
        print('     %-5d  %-11.4e  %-10.4f  %-11.4f  %.4f' % (L, rho, res[0] / rho, res[1] / rho, res[2] / rho))
    print('   => a bad word typically rises only O(sqrt L): its own-path weight 3^a_j/2^j is 2^O(sqrt L).')
    print('   annealed contrast: E[2^(S_j) | Bad_L] (exact DP) grows exponentially in j (rare heavy prefixes dominate):')
    for L in (50, 100, 200):
        T = [thr(k) for k in range(L + 1)]
        out = []
        for j in (L // 4, L // 2, 3 * L // 4):
            dist = {0: (1.0, 1.0)}
            for k in range(1, L + 1):
                nd = {}
                for a, (pp, pw) in dist.items():
                    for up in (0, 1):
                        b = a + up
                        if b >= T[k]:
                            q = nd.get(b, (0.0, 0.0))
                            nd[b] = (q[0] + 0.5 * pp, q[1] + 0.5 * pw)
                if k == j:
                    nd = {a: (pp, pp * 2.0 ** (a * LS - k)) for a, (pp, pw) in nd.items()}
                dist = nd
            rho = sum(pp for pp, pw in dist.values()); ew = sum(pw for pp, pw in dist.values())
            out.append('j=%d: %.3e (log2/j = %.3f)' % (j, ew / rho, math.log2(ew / rho) / j))
        print('     L=%-4d ' % L + '   '.join(out))
    print('\nPART 5 done in %.1fs' % (time.time() - t0))

if __name__ == '__main__':
    main()
