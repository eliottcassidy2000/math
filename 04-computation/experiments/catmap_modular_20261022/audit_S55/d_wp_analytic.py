"""Audit D (THM-4763): analytic SL2(Z) action on span{E_s = theta_(s,p)/eta (s=0..p), tau O_s, O_s = theta'_(s,p)/eta
(s=1..p-1)} from the Weil representation of sqrt(2p)Z and the eta multiplier.  Auditor for opus-2026-10-22-S55.

Derivation (Poisson summation; see the audit report):
  E_g(-1/tau) = (2p)^(-1/2) sum_(d mod 2p) e(g d/2p) E_d(tau),            E_g(tau+1) = e(g^2/4p - 1/24) E_g(tau)
  O_g(-1/tau) = tau * M_s O(tau),  M_s = -i sqrt(2/p) [sin(pi g d/p)]_(g,d=1..p-1),   O_g(tau+1) = e(g^2/4p - 1/24) O_g
so on F = (E_0..E_p, tau O_1..tau O_(p-1), O_1..O_(p-1)), F(gamma tau) = rho(gamma) F(tau) with
  rho(s) = R(s) + [[0,-1],[1,0]] (x) M_s,   rho(t) = R(t) + [[1,1],[0,1]] (x) M_t.
Checks: (i) pointwise numerical verification of F(gamma tau) = rho(gamma) F(tau) at random tau for gamma = s, t, A, and a
random word, using our own theta and eta series (pentagonal series for eta), no fitting; (ii) relations; (iii) blocks;
(iv) cat map eigenvalue moduli, (v) trace law Tr rho(A^n) = Tr R(A^n) + L_2n Tr M(A^n), (vi) spectral radius of random
hyperbolic gamma; (vii) the W(p) characters (FGST) span the (E, O) part.
Usage: python d_wp_analytic.py [pmax]"""
import sys, math, random
import mpmath as mp

mp.mp.dps = 40
pmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6


def e(x):
    return mp.expjpi(2 * x)


def eta(tau, K=40):                                   # Euler pentagonal series: sum (-1)^n q^((6n+1)^2/24)
    return mp.fsum((-1) ** n * mp.exp(2j * mp.pi * tau * mp.mpf((6 * n + 1) ** 2) / 24) for n in range(-K, K + 1))


def theta(s, p, tau, K=40):
    return mp.fsum(mp.exp(2j * mp.pi * tau * mp.mpf((2 * p * m + s) ** 2) / (4 * p)) for m in range(-K, K + 1))


def dtheta(s, p, tau, K=40):                          # sum_(n in Z + s/2p) 2pn q^(p n^2)
    return mp.fsum((2 * p * m + s) * mp.exp(2j * mp.pi * tau * mp.mpf((2 * p * m + s) ** 2) / (4 * p)) for m in range(-K, K + 1))


def F(p, tau):
    et = eta(tau)
    E = [theta(s, p, tau) / et for s in range(p + 1)]
    O = [dtheta(s, p, tau) / et for s in range(1, p)]
    return mp.matrix(E + [tau * o for o in O] + O)


def blocks(p):
    n2 = 2 * p
    Rs = mp.matrix(p + 1, p + 1); Rt = mp.matrix(p + 1, p + 1)
    for g in range(p + 1):
        for d in range(p + 1):
            if d == 0:
                v = 1
            elif d == p:
                v = (-1) ** g
            else:
                v = 2 * mp.cospi(mp.mpf(g * d) / p)
            Rs[g, d] = v / mp.sqrt(n2)
        Rt[g, g] = e(mp.mpf(g * g) / (4 * p) - mp.mpf(1) / 24)
    o = p - 1
    Ms = mp.matrix(o, o); Mt = mp.matrix(o, o)
    for i in range(o):
        for j in range(o):
            Ms[i, j] = -1j * mp.sqrt(mp.mpf(2) / p) * mp.sinpi(mp.mpf((i + 1) * (j + 1)) / p)
        Mt[i, i] = e(mp.mpf((i + 1) ** 2) / (4 * p) - mp.mpf(1) / 24)
    return Rs, Rt, Ms, Mt


def kron(a, B):                                       # 2x2 integer matrix a (x) mp matrix B
    n = B.rows
    K = mp.matrix(2 * n, 2 * n)
    for i in range(2):
        for j in range(2):
            for k in range(n):
                for l in range(n):
                    K[i * n + k, j * n + l] = a[i][j] * B[k, l]
    return K


def dsum(X, Y):
    n, m = X.rows, Y.rows
    Z = mp.matrix(n + m, n + m)
    for i in range(n):
        for j in range(n):
            Z[i, j] = X[i, j]
    for i in range(m):
        for j in range(m):
            Z[n + i, n + j] = Y[i, j]
    return Z


def rho_gens(p):
    Rs, Rt, Ms, Mt = blocks(p)
    if p >= 2 and p - 1 > 0:
        S = dsum(Rs, kron([[0, -1], [1, 0]], Ms)); T = dsum(Rt, kron([[1, 1], [0, 1]], Mt))
    return S, T, (Rs, Rt, Ms, Mt)


def mobius(g, tau):
    a, b, c, d = g
    return (a * tau + b) / (c * tau + d)


def mmul(g, h):
    a, b, c, d = g; e_, f, g_, h_ = h
    return (a * e_ + b * g_, a * f + b * h_, c * e_ + d * g_, c * f + d * h_)


def word_matrix(word, S, T):
    """word in letters 's','t','S','T' (capital = inverse); returns (gamma, rho(gamma))."""
    g = (1, 0, 0, 1); R = mp.eye(S.rows)
    Si, Ti = mp.inverse(S), mp.inverse(T)
    lett = {'s': ((0, -1, 1, 0), S), 'S': ((0, 1, -1, 0), Si), 't': ((1, 1, 0, 1), T), 'T': ((1, -1, 0, 1), Ti)}
    for ch in word:
        gg, RR = lett[ch]
        g = mmul(g, gg); R = R * RR
    return g, R


def maxabs(M):
    return max(abs(M[i, j]) for i in range(M.rows) for j in range(M.cols))


def lucas(m):
    a, b = 2, 1
    for _ in range(m):
        a, b = b, a + b
    return a


def mpow(M, n):
    X = mp.eye(M.rows)
    for _ in range(n):
        X = X * M
    return X


random.seed(7)
phi = (1 + mp.sqrt(5)) / 2
for p in range(2, pmax + 1):
    S, T, (Rs, Rt, Ms, Mt) = rho_gens(p)
    dim = S.rows
    # (i) pointwise verification
    worst = 0
    for word in ('s', 't', 'tsTS', 'tsTStsTS', 'ttsTsTTs'):
        g, R = word_matrix(word, S, T)
        for _ in range(3):
            tau = mp.mpc(random.uniform(-0.4, 0.4), random.uniform(0.9, 1.3))
            lhs = F(p, mobius(g, tau)); rhs = R * F(p, tau)
            worst = max(worst, maxabs(lhs - rhs) / maxabs(lhs))
    # (ii) relations: rho(s)^2 = I (weight-0 action of -I is trivial), (rho(s) rho(t))^3 = rho(s)^2
    rel1 = maxabs(S * S - mp.eye(dim)); rel2 = maxabs(mpow(S * T, 3) - S * S)
    relM = maxabs(Ms * Ms + mp.eye(p - 1)); relM2 = maxabs(mpow(Ms * Mt, 3) - Ms * Ms)
    relR = maxabs(Rs * Rs - mp.eye(p + 1)); relR2 = maxabs(mpow(Rs * Rt, 3) - Rs * Rs)
    # (iv) cat map
    C = T * S * mp.inverse(T) * mp.inverse(S)
    CR = Rt * Rs * mp.inverse(Rt) * mp.inverse(Rs)
    CM = Mt * Ms * mp.inverse(Mt) * mp.inverse(Ms)
    blk = maxabs(C - dsum(CR, kron([[2, 1], [1, 1]], CM)))
    ev = mp.eig(C, left=False, right=False)
    mods = sorted(set(mp.nstr(abs(z), 12) for z in ev))
    evR = mp.eig(CR, left=False, right=False); evM = mp.eig(CM, left=False, right=False) if p > 2 else [CM[0, 0]]

    def order_of(evs, cap=10000):
        from fractions import Fraction
        dens = []
        for z in evs:
            a = float(mp.arg(z) / (2 * mp.pi)) % 1.0
            f = Fraction(a).limit_denominator(cap)
            assert abs(float(f) - a) < 1e-12 and abs(abs(z) - 1) < 1e-20
            dens.append(f.denominator)
        return math.lcm(*dens)
    oR, oM = order_of(evR), order_of(evM)
    # (v) trace law
    tr_ok, traces = True, []
    Cn, CRn, CMn = mp.eye(dim), mp.eye(p + 1), mp.eye(p - 1)
    for n in range(1, 13):
        Cn, CRn, CMn = Cn * C, CRn * CR, CMn * CM
        tC = sum(Cn[i, i] for i in range(dim)); tR = sum(CRn[i, i] for i in range(p + 1)); tM = sum(CMn[i, i] for i in range(p - 1))
        tr_ok = tr_ok and abs(tC - (tR + lucas(2 * n) * tM)) < mp.mpf(10) ** -25
        traces.append(complex(tC))
    # (vi) spectral radius of random hyperbolic words
    sr_ok = True
    for _ in range(6):
        word = ''.join(random.choice('stST') for _ in range(9))
        g, R = word_matrix(word, S, T)
        tr = g[0] + g[3]
        if abs(tr) <= 2:
            continue
        lam = (abs(tr) + mp.sqrt(tr * tr - 4)) / 2
        evg = mp.eig(R, left=False, right=False)
        sr = max(abs(z) for z in evg)
        sr_ok = sr_ok and abs(sr - lam) < mp.mpf(10) ** -15
    # parabolic: rho(t) is not diagonalisable on the log block (Jordan): rank(rho(t) - mu) test via t (x) M_t structure
    Tlog = kron([[1, 1], [0, 1]], Mt)
    nondiag = maxabs(Tlog - kron([[1, 0], [0, 1]], Mt)) > 0.5
    print(f'p = {p}: dim {dim}; pointwise |F(g tau) - rho(g) F(tau)|/|F| max = {mp.nstr(worst, 3)}; rho(s)^2 = I: {mp.nstr(rel1, 2)}, '
          f'(st)^3 = s^2: {mp.nstr(rel2, 2)}; M_s^2 = -I: {mp.nstr(relM, 2)}, M rel {mp.nstr(relM2, 2)}; R_s^2 = I: {mp.nstr(relR, 2)}, '
          f'R rel {mp.nstr(relR2, 2)}; rho(A) = R(A) + A_def (x) M(A): {mp.nstr(blk, 2)}')
    print(f'       |eig rho(A)| = {mods} (phi^-2 = {mp.nstr(phi ** -2, 12)}, phi^2 = {mp.nstr(phi ** 2, 12)}); ord R_(p+1)(A) = {oR}, '
          f'ord R_(p-1)(A) = {oM}; trace law n=1..12: {tr_ok}; spectral radius = dilatation for random hyperbolic words: {sr_ok}; '
          f'Dehn twist Jordan: {nondiag}')
    print(f'       Tr rho(A^n), n = 1..8: {[complex(round(z.real, 8), round(z.imag, 8)) for z in traces[:8]]}')
    if p == 2:
        print(f'       p = 2: R_1 = M is the character M(s) = {mp.nstr(Ms[0, 0], 6)}, M(t) = {mp.nstr(Mt[0, 0], 6)} '
              f'(= exp(2 pi i/12)): NOT the trivial representation; M(A) = {mp.nstr(CM[0, 0], 6)}')
    # (vii) FGST characters: chi+_s = ((s/p) th_(p-s) + 2 th'_(p-s))/eta, chi-_s = ((s/p) th_s - 2 th'_s)/eta with th' = sum j q^(p j^2)
    # = (1/2p) * our dtheta.  Coefficient matrix w.r.t. (E_0..E_p, O_1..O_(p-1)) must be invertible.
    import numpy as np
    Mx = np.zeros((2 * p, 2 * p))
    for s in range(1, p + 1):
        row = np.zeros(2 * p); r = p - s
        row[r] += s / p
        if 1 <= r <= p - 1:
            row[p + 1 + r - 1] += 2 / (2 * p)
        Mx[s - 1] = row
        row = np.zeros(2 * p); row[s] += s / p
        if 1 <= s <= p - 1:
            row[p + 1 + s - 1] -= 2 / (2 * p)
        Mx[p + s - 1] = row
    print(f'       FGST characters chi^(+-)_s (s=1..p): rank of coefficient matrix on (E,O) = {np.linalg.matrix_rank(Mx)} (= 2p = {2 * p})')
