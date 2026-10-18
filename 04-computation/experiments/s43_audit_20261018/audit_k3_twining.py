#!/usr/bin/env python3
"""S43 referee: independent checks of section 2 of factorization_vertex_supersymmetric_bridge_20261018.md (2026-10-18).

A. Quartic LG orbifold (Kawai-Yamada-Yang, chat = 2) vs 2 phi_{0,1} at random points (mpmath theta functions).
B. Mathieu moonshine numbers by an independent route: H(tau) = Z_K3 * eta^3/theta_1^2 - 24 mu(tau,z) with the
   Appell-Lerch sum mu (Eguchi-Ooguri-Tachikawa decomposition), evaluated numerically at small q with high precision,
   then coefficients peeled off one by one.
C. Twining genera: exact q-series from T_g formulas; comparison with Cheng 2010 (Table 1 and Appendix C):
   11A: (11/6) phi2^(11) - (22/5) eta^2 eta(11)^2;
   23A: (23/12) phi2^(23) - (23/22) f_{23,1} - (161/22) f_{23,2}, f_{23,1} = f1 + f2, f_{23,2} = (f1 - f2)/sqrt5 with
        f1, f2 the Hecke eigenforms of S_2(Gamma_0(23)), f1 = q + ((-1+sqrt5)/2) q^2 - sqrt5 q^3 + ... (Cheng (B.5));
   vs the session's -(23/11)[2 eta^2 eta(23)^2 + eta eta(23) theta_{x^2+xy+6y^2}].
   Integrality and Tr(g|K_n) == dim K_n (mod p) for n <= NMAX.
Usage: python audit_k3_twining.py [A B C]
"""
import sys, random
from fractions import Fraction as Fr

# ---------------- A ----------------
def partA():
    import mpmath as mp
    mp.mp.dps = 30
    def th(i, tau, z):
        q = mp.exp(1j * mp.pi * tau)
        return mp.jtheta(i, mp.pi * z, q)
    def Zw(tau, z, w):
        return th(1, tau, (1 - w) * z) / th(1, tau, w * z)
    def phi01(tau, z):
        return 4 * sum((th(i, tau, z) / th(i, tau, 0)) ** 2 for i in (2, 3, 4))
    def orb(tau, z):
        h = 4; chat = 2
        tot = 0
        for a in range(h):
            for b in range(h):
                pre = mp.exp(2j * mp.pi * tau * chat * a * a / 2) * mp.exp(2j * mp.pi * z * chat * a) * (-1) ** (chat * (a + b + a * b))
                tot += pre * Zw(tau, z + a * tau + b, mp.mpf(1) / 4) ** 4
        return tot / h
    random.seed(11)
    worst = 0
    for _ in range(8):
        tau = mp.mpc(random.uniform(-0.5, 0.5), random.uniform(0.7, 1.5))
        z = mp.mpc(random.uniform(-0.5, 0.5), random.uniform(-0.15, 0.15))
        e = abs(orb(tau, z) - 2 * phi01(tau, z)) / abs(2 * phi01(tau, z))
        worst = max(worst, e)
    lim = [mp.nstr(Zw(mp.mpc(0.1, 1.1), mp.mpf('1e-12'), mp.mpf(1) / n), 12) for n in (3, 4, 5, 6)]
    print('[A] quartic Z_4 LG orbifold = 2 phi_{0,1}: max relative error over 8 random points %s; Z_{x^n}(tau, z->0) for n = 3..6: %s'
          % (mp.nstr(worst, 3), lim), flush=True)

# ---------------- B ----------------
def partB(NC=13):
    import mpmath as mp
    mp.mp.dps = 260
    tau = mp.mpc('0.0', '4.0')            # q = e^{-8 pi} ~ 1.2e-11
    z = mp.mpf('0.2137')
    q = mp.exp(2j * mp.pi * tau)
    y = mp.exp(2j * mp.pi * z)
    nome = mp.exp(1j * mp.pi * tau)
    th = lambda i, zz: mp.jtheta(i, mp.pi * zz, nome)
    eta = mp.exp(2j * mp.pi * tau / 24) * mp.qp(q)
    ZK3 = 8 * sum((th(i, z) / th(i, 0)) ** 2 for i in (2, 3, 4))
    # Appell-Lerch mu(tau,z) = -i y^{1/2}/theta1(tau,z) * sum_n (-1)^n y^n q^{n(n+1)/2}/(1 - y q^n)
    # (theta1 in the convention theta1 = -i q^{1/8} y^{1/2} prod(...); with mpmath's jtheta(1) = i*that, we fix the sign
    #  by requiring the q^(-1/8) coefficient of H to be -2.)
    s = 0
    for n in range(-60, 61):
        s += (-1) ** n * y ** n * q ** (mp.mpf(n * (n + 1)) / 2) / (1 - y * q ** n)
    th1 = th(1, z)
    best = None
    for sgn in (1, -1):
        mu = sgn * (-1j) * mp.sqrt(y) / th1 * s
        Hval = ZK3 * eta ** 3 / th1 ** 2 - 24 * mu
        c = Hval * q ** (mp.mpf(1) / 8) / 2
        if best is None or abs(c - (-1)) < abs(best[0] - (-1)):
            best = (c, sgn)
    c, sgn = best
    coeffs = []
    rem = c
    for n in range(NC):
        cn = rem / q ** n
        k = int(mp.nint(mp.re(cn)))
        coeffs.append(k)
        rem -= k * q ** n
    print('[B] numerical EOT decomposition at tau = 4i: H = 2 q^(-1/8) (sum A_n q^n), A_0..A_%d = %s' % (NC - 1, coeffs), flush=True)
    return coeffs

# ---------------- C ----------------
def sigma(n):
    return sum(d for d in range(1, n + 1) if n % d == 0)

def E2(N, s=1):
    e = [Fr(0)] * (N + 1)
    e[0] = Fr(1)
    for n in range(1, N // s + 1):
        e[n * s] = Fr(-24 * sigma(n))
    return e

def phi2(Nlev, N):
    a = E2(N, Nlev); b = E2(N)
    return [(Nlev * a[i] - b[i]) / (Nlev - 1) for i in range(N + 1)]

def mul(a, b, N):
    out = [Fr(0)] * (N + 1)
    for i in range(N + 1):
        if a[i]:
            for j in range(N + 1 - i):
                if b[j]:
                    out[i + j] += a[i] * b[j]
    return out

def inv(a, N):
    r = [Fr(0)] * (N + 1)
    r[0] = 1 / a[0]
    for n in range(1, N + 1):
        r[n] = -sum(a[k] * r[n - k] for k in range(1, n + 1)) / a[0]
    return r

def eta_prod(shape, N):
    """q^(sum a m /24) * returned series"""
    c = [Fr(0)] * (N + 1)
    c[0] = Fr(1)
    for a, m in shape.items():
        for _ in range(m):
            for n in range(1, N // a + 1):
                st = a * n
                for k in range(N, st - 1, -1):
                    c[k] -= c[k - st]
    return c

def shift(series, k, N):
    return ([Fr(0)] * k + series)[:N + 1]

def H_series(N):
    """independent exact H: from the theta-function form of the K3 elliptic genus at z = 1/2 ... (not used);
    here we use the classical formula H = (-2 E2 + 48 F2)/eta^3 only as a cross-check of part B."""
    e3 = [Fr(0)] * (N + 1)
    n = 0
    while n * (n + 1) // 2 <= N:
        e3[n * (n + 1) // 2] += (-1) ** n * (2 * n + 1)
        n += 1
    F2 = [Fr(0)] * (N + 1)
    for r in range(1, 2 * N + 3):
        for s in range(1, r):
            if (r - s) % 2 and r * s // 2 <= N:
                F2[r * s // 2] += (-1) ** r * s
    num = [-2 * x for x in E2(N)]
    num = [num[i] + 48 * F2[i] for i in range(N + 1)]
    return mul(num, inv(e3, N), N), inv(e3, N)

def partC(N=60):
    H, ie3 = H_series(N)
    A = [x / 2 for x in H]
    print('[C] A_n from (-2E2+48F2)/eta^3: %s' % [int(x) for x in A[:13]], flush=True)
    def Hg(chi, T):
        t = mul(T, ie3, N)
        return [Fr(chi, 24) * H[i] - t[i] for i in range(N + 1)]
    results = {}
    for p, chi in ((2, 8), (3, 6), (5, 4), (7, 3)):
        T = [Fr(2 * p, p + 1) * x for x in phi2(p, N)]
        Ag = [x / 2 for x in Hg(chi, T)]
        results[p] = Ag
        print('[C] %dA: Tr(g|K_n) = 2 A_n(g), n = 1..7: %s' % (p, [str(2 * x) for x in Ag[1:8]]), flush=True)
    # 11A (Cheng): (11/6) phi2^(11) - (22/5) eta^2 eta(11)^2
    f11 = shift(eta_prod({1: 2, 11: 2}, N), 1, N)
    T11 = [Fr(11, 6) * x for x in phi2(11, N)]
    T11 = [T11[i] - Fr(22, 5) * f11[i] for i in range(N + 1)]
    A11 = [x / 2 for x in Hg(2, T11)]
    ok11 = all((2 * A11[n]).denominator == 1 and (2 * A11[n] - 2 * A[n]) % 11 == 0 for n in range(1, N + 1))
    print('[C] 11A (Cheng App. C): A_n, n = 0..12: %s; traces integral and == dim mod 11 for n <= %d: %s'
          % ([str(x) for x in A11[:13]], N, ok11), flush=True)
    # level 23: basis u = eta^2 eta(23)^2 (starts q^2), w = eta eta(23) theta_Q1, Q1 = x^2 + xy + 6y^2
    u = shift(eta_prod({1: 2, 23: 2}, N), 2, N)
    e1 = shift(eta_prod({1: 1, 23: 1}, N), 1, N)
    th1 = [Fr(0)] * (N + 1)
    th2 = [Fr(0)] * (N + 1)
    R = int(N ** 0.5) + 3
    for x in range(-2 * R, 2 * R + 1):
        for yv in range(-R, R + 1):
            v1 = x * x + x * yv + 6 * yv * yv
            v2 = 2 * x * x + x * yv + 3 * yv * yv
            if v1 <= N:
                th1[v1] += 1
            if v2 <= N:
                th2[v2] += 1
    w = mul(e1, th1, N)
    # check eta(tau)eta(23 tau) = (theta_Q1 - theta_Q2)/2
    okq = all(e1[i] == (th1[i] - th2[i]) / 2 for i in range(N + 1))
    # Hecke T_2 on the 2-dim space spanned by u, w (weight 2, trivial character): (T2 f)_n = a_{2n} + 2 a_{n/2}
    def T2(f, M):
        return [f[2 * n] + (2 * f[n // 2] if n % 2 == 0 else 0) for n in range(M + 1)]
    M = N // 2
    tu, tw = T2(u, M), T2(w, M)
    # express T2 u, T2 w in basis (u, w) using coefficients q^1, q^2: u = (0, 1, ...), w = (1, 1, ...)
    def coords(f):
        b = f[1]            # coefficient of w
        a = f[2] - b * 1    # u has q^2 coefficient 1, w has q^2 coefficient w[2]
        a = f[2] - b * w[2]
        return a, b
    cu, cw = coords(tu), coords(tw)
    okT = all(tu[n] == cu[0] * u[n] + cu[1] * w[n] for n in range(M + 1)) and all(tw[n] == cw[0] * u[n] + cw[1] * w[n] for n in range(M + 1))
    print('[C] level 23: eta eta(23) = (theta_Q1 - theta_Q2)/2: %s; T_2 preserves span(u, w) (checked to q^%d): %s; '
          'T_2 matrix columns (u-coord, w-coord): T2u = %s, T2w = %s' % (okq, M, okT, cu, cw), flush=True)
    # eigenforms f = alpha u + w normalised (q coefficient 1); eigenvalue lambda = a_2 = alpha + w[2]
    # T2 f = lambda f  ->  alpha*cu + cw = lambda*(alpha, 1)
    # second coordinate: alpha*cu[1] + cw[1] = lambda ; first: alpha*cu[0] + cw[0] = lambda*alpha
    import sympy as sp
    al = sp.symbols('al')
    lam = al * sp.Rational(cu[1].numerator, cu[1].denominator) + sp.Rational(cw[1].numerator, cw[1].denominator)
    sols = sp.solve(sp.Eq(al * sp.Rational(cu[0].numerator, cu[0].denominator) + sp.Rational(cw[0].numerator, cw[0].denominator), lam * al), al)
    eig = []
    for s_ in sols:
        a2 = sp.nsimplify(s_ + w[2])
        eig.append((sp.simplify(s_), sp.simplify(a2)))
    print('[C] Hecke eigenforms f = alpha u + w: (alpha, a_2) = %s' % eig, flush=True)
    # Cheng's f~_{23,1}: a_2 = (-1 + sqrt5)/2
    s5 = sp.sqrt(5)
    target = (-1 + s5) / 2
    al1 = [e[0] for e in eig if sp.simplify(e[1] - target) == 0][0]
    al2 = [e[0] for e in eig if sp.simplify(e[1] - target) != 0][0]
    uS = [sp.Rational(x.numerator, x.denominator) for x in u]
    wS = [sp.Rational(x.numerator, x.denominator) for x in w]
    f1 = [sp.expand(al1 * uS[i] + wS[i]) for i in range(N + 1)]
    f2 = [sp.expand(al2 * uS[i] + wS[i]) for i in range(N + 1)]
    F1 = [sp.simplify(f1[i] + f2[i]) for i in range(N + 1)]
    F2 = [sp.simplify((f1[i] - f2[i]) / s5) for i in range(N + 1)]
    print('[C] Cheng f~_{23,1} = %s ...' % ' + '.join('(%s)q^%d' % (sp.simplify(f1[i]), i) for i in range(1, 7)), flush=True)
    print('[C] f_{23,1} = %s; f_{23,2} = %s (q^1..q^8)' % ([F1[i] for i in range(1, 9)], [F2[i] for i in range(1, 9)]), flush=True)
    cheng = [sp.Rational(-23, 22) * F1[i] - sp.Rational(161, 22) * F2[i] for i in range(N + 1)]
    sess = [sp.Rational(-23, 11) * (2 * sp.Rational(u[i].numerator, u[i].denominator) + sp.Rational(w[i].numerator, w[i].denominator)) for i in range(N + 1)]
    same = all(sp.simplify(cheng[i] - sess[i]) == 0 for i in range(N + 1))
    f232_is_u = all(sp.simplify(F2[i] - sp.Rational(u[i].numerator, u[i].denominator)) == 0 for i in range(N + 1))
    print('[C] session 23AB cusp form == Cheng -(23/22) f_{23,1} - (161/22) f_{23,2} to q^%d: %s; f_{23,2} == eta^2 eta(23)^2: %s'
          % (N, same, f232_is_u), flush=True)
    # traces for 23A
    T23 = [Fr(23, 12) * x for x in phi2(23, N)]
    T23 = [T23[i] - Fr(23, 11) * (2 * u[i] + w[i]) for i in range(N + 1)]
    A23 = [x / 2 for x in Hg(1, T23)]
    ok23 = all((2 * A23[n]).denominator == 1 and (2 * A23[n] - 2 * A[n]) % 23 == 0 for n in range(1, N + 1))
    print('[C] 23A: A_n, n = 0..12: %s; traces integral and == dim mod 23 for n <= %d: %s' % ([str(x) for x in A23[:13]], N, ok23), flush=True)

if __name__ == '__main__':
    for q in (sys.argv[1:] or ['A', 'B', 'C']):
        {'A': partA, 'B': partB, 'C': partC}[q]()
