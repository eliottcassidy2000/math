#!/usr/bin/env python3
"""procgen_cauchy_20260925 -- part 3: a certified (computer-assisted) growth bound for ||Z_k||^2, hence for M2(L),
hence a PROVED exponent for the Cauchy-Schwarz lower bound on the provability price.

Class-mass majorant at level r.  For a class c mod 3^r let u_k(c) = int_{c + 3^r Z_3} Z_k^2 dmu, and U_k the level-(r-1)
sums.  From Z_{k+1} = L Z_k and one Cauchy-Schwarz step on the cross term (proof in the note, section 3):
    u_{k+1}(c) <= Phi(u_k)(c) := u_k(2c)/4 + [c = 2 mod 3] ( 3/4 U_k(E c) + (sqrt 3)/2 sqrt( u_k(2c) U_k(E c) ) ),
E c = (2c-1)/3 mod 3^(r-1).  Phi is monotone and positively 1-homogeneous, so a vector e > 0 with Phi(e) <= theta e and
u_{k0} <= C e gives u_k <= C theta^(k-k0) e for all k >= k0.  Classes c = 0 mod 3 carry exactly u_k(c) = 4^-k 3^-r
and get e(c) = eps > 0.  The ratio max_c Phi(e)(c)/e(c) is computed in IEEE double precision (a few flops and one
correctly rounded sqrt per entry, relative error < 1e-14) and inflated by the factor 1 + 1e-9 before use.
Output: stdout.  Peak memory about 450 MB (r = 15).
"""
import math, sys, time
import numpy as np

RMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 15
K0 = 14
LOG32 = math.log(2) / math.log(3)
ETA = 1 + (LOG32 * math.log2(LOG32) + (1 - LOG32) * math.log2(1 - LOG32))
SQ = math.sqrt(3) / 2

def build_I(K):
    I = np.ones(1, dtype=np.int64)
    out = [I]
    for k in range(1, K + 1):
        m = 3 ** k; mp = 3 ** (k - 1)
        v = np.arange(m, dtype=np.int32); v *= 2; v %= mp
        cur = I[v]; del v
        t = np.arange(mp, dtype=np.int32); t *= 2; t += 1; t %= mp
        cur[2::3] += 3 * I[t]; del t
        I = cur
        out.append(I)
    return out

def Phi(e, twoc, Ec, Nb):
    U = e.reshape(3, Nb).sum(axis=0)
    new = e[twoc]; new *= 0.25
    b = U[Ec]; del U
    tmp = e[twoc[2::3]]; tmp *= b; np.sqrt(tmp, out=tmp); tmp *= SQ
    b *= 0.75; tmp += b; del b
    new[2::3] += tmp
    return new

def ratio_range(e, twoc, Ec, Nb, mask):
    ph = Phi(e, twoc, Ec, Nb)
    np.divide(ph, e, out=ph, where=mask)
    hi = float(np.max(ph, where=mask, initial=0.0)); lo = float(np.min(ph, where=mask, initial=np.inf))
    del ph
    return hi, lo

def main():
    t0 = time.time()
    print('=' * 100)
    print('PART 3. Certified growth rate of ||Z_k||^2 by the level-r class-mass majorant (computer-assisted proof)')
    print('=' * 100)
    Is = build_I(K0)
    Zsq = [float(np.dot(Is[k].astype(np.float64), Is[k].astype(np.float64))) / (4 ** k * 3 ** k) for k in range(K0 + 1)]
    Z14 = Is[K0].astype(np.float64) / 2 ** K0
    del Is
    e_prev = None
    results = []
    for r in range(1, RMAX + 1):
        N = 3 ** r; Nb = 3 ** (r - 1)
        c = np.arange(N, dtype=np.int64)
        twoc = ((2 * c) % N).astype(np.int32)
        zero = (c % 3 == 0)
        del c
        t = np.arange(Nb, dtype=np.int64)
        Ec = ((2 * t + 1) % Nb).astype(np.int32)
        del t
        if e_prev is None:
            e = np.ones(N)
        else:
            e = np.tile(e_prev, 3) / 3.0
        e[zero] = 0.0
        e /= e.sum()
        th = None
        for it in range(3000):
            new = Phi(e, twoc, Ec, Nb)
            s = new.sum()
            new /= s
            e = new; del new
            if it % 25 == 24:
                hi, lo = ratio_range(e, twoc, Ec, Nb, ~zero)
                if hi - lo < 1e-9 * hi:
                    break
        # initial condition from the exact Z_14 (constant on classes mod 3^14): u_14(c) = int_c Z_14^2
        if r <= K0:
            u0 = (Z14 ** 2).reshape(3 ** (K0 - r), N).sum(axis=0) / 3 ** K0
        else:
            u0 = np.tile(Z14 ** 2, 3 ** (r - K0)) / 3 ** r
        nz = ~zero
        q = np.divide(u0, e, out=np.zeros(N), where=nz)
        Cplus = float(q.max()); del q
        eps = (4.0 ** (-K0) * 3.0 ** (-r)) / Cplus          # so that u0 <= Cplus * e on the zero classes too
        e[zero] = eps
        allm = np.ones(N, dtype=bool)
        R, Rmin = ratio_range(e, twoc, Ec, Nb, allm)
        del allm
        theta = R * (1 + 1e-9)
        q = u0 / e
        C = max(Cplus, float(q.max())); del q
        S = float(e.sum())
        emin = float(np.min(e, where=nz, initial=np.inf)) / float(np.mean(e))
        results.append((r, theta, C, S))
        print('  r=%2d  certified theta_r = %.7f  (min ratio %.7f)  log2 theta = %.5f   C*sum(e) = %.4f  (vs ||Z_14||^2 = %.4f)  min e/mean e = %.2e  %.1fs'
              % (r, theta, Rmin, math.log2(theta), C * S, Zsq[K0], emin, time.time() - t0), flush=True)
        if r == RMAX:
            # direct iteration of the majorant from the exact u_14 (better constants for explicit L)
            NIT = 240
            u = u0; del u0; direct = [float(u.sum())]
            for it in range(NIT):
                u = Phi(u, twoc, Ec, Nb)
                direct.append(float(u.sum()) * (1 + 1e-12) ** (it + 1))
            u /= e
            Ctail = float(u.max()) * (1 + 1e-12) ** NIT
            del u
        else:
            del u0
        e_prev = e
        del e, twoc, Ec, zero
    r, theta, C, S = results[-1]
    print('\n(3a) PROVED bound: for all k >= %d,  ||Z_k||^2 <= %.4f * %.7f^(k-%d)   (level r = %d).' % (K0, C * S, theta, K0, r))
    print('     hence M2(L) = ||W_L||^2 <= (sum_{k<L} ||Z_k||)^2 <= C\' theta^L and, with rho_L >= 2^(-eta L)/poly(L) (THM-4475 B),')
    print('     delta_L >= rho_L^2 / M2(L) >= 2^(-(2 eta + log2 theta) L - O(log L)),   2 eta + log2 theta = %.5f + %.5f = %.5f'
          % (2 * ETA, math.log2(theta), 2 * ETA + math.log2(theta)))
    print('     (THM-4475: 0.7737.  Conjectured from the exact data: 2 eta = %.5f, since ||Z_k||^2 grows linearly.)' % (2 * ETA))
    print('\n(3b) explicit rigorous values of the bound  B(L) >= M2(L)  and the resulting lower bounds on delta_L:')
    print('   L     B(L)          rho_L^2/B(L)   (log2; per L)            THM-4475: rho_L/sum_{k<L} g_k')
    gs = [1.0, 2.0]
    while len(gs) < 500:
        gs.append(1.5 * gs[-1] + 0.25 * gs[-2])
    def bad_count(L):
        T = []
        for k in range(L + 1):
            a = 0
            while 3 ** a <= 2 ** k: a += 1
            T.append(a)
        dist = {0: 1}
        for k in range(1, L + 1):
            nd = {}
            for a, cc in dist.items():
                for up in (0, 1):
                    b = a + up
                    if b >= T[k]: nd[b] = nd.get(b, 0) + cc
            dist = nd
        return sum(dist.values())
    def zbound(k):     # rigorous upper bound for ||Z_k||^2, k > K0
        j = k - K0
        if j <= NIT:
            return direct[j]
        return Ctail * S * theta ** (j - NIT)
    print('   rigorous ||Z_k||^2 bounds from the direct majorant iteration: k=20: %.3f, 30: %.3f, 50: %.3f, 100: %.3f, 250: %.3f'
          % (zbound(20), zbound(30), zbound(50), zbound(100), zbound(250)))
    for L in (16, 20, 22, 24, 26, 30, 40, 60, 100, 200, 400):
        s = sum(math.sqrt(Zsq[k]) for k in range(min(L, K0 + 1)))
        s += sum(math.sqrt(zbound(k)) for k in range(K0 + 1, L))
        B = s * s
        rho = bad_count(L) / 2 ** L
        old = rho / sum(gs[:L])
        print('   %-4d  %-12.4e  %.4e     (%.2f; %.4f)         %.4e' % (L, B, rho * rho / B, math.log2(rho * rho / B), -math.log2(rho * rho / B) / L, old))
    print('\nPART 3 done in %.1fs' % (time.time() - t0))

if __name__ == '__main__':
    main()
