#!/usr/bin/env python3
"""From a Landau-Ginzburg derived critical locus to Mathieu moonshine: elliptic genera computed end to end (S43).

Part 1 (numerical, complex tau, z): the elliptic genus of the LG model W = x^n (one chiral field of R-weight w = 1/n) is
Witten's Z_w(tau, z) = theta1(tau, (1-w) z)/theta1(tau, w z); check Z(tau, 0) = n - 1 = dim of the Jacobian ring
C[x]/(x^(n-1)), the cohomology of the derived critical locus (Koszul complex) of W. The LG orbifold of
W = x1^4 + x2^4 + x3^4 + x4^4 by Z_4 (Kawai-Yamada-Yang):
    Z_orb = (1/h) sum_{a,b=0}^{h-1} (-1)^{chat(a+b+ab)} q^{chat a^2/2} y^{chat a} prod_i Z_{w_i}(tau, z + a tau + b),
chat = sum(1 - 2 w_i) = 2, h = 4, is compared with the K3 elliptic genus 2 phi_{0,1} = 8 sum_{i=2,3,4}
(theta_i(tau,z)/theta_i(tau,0))^2 at random points.
Part 2 (exact q-series): Eguchi-Hikami's H(tau) = (-2 E_2(tau) + 48 F_2(tau))/eta(tau)^3, F_2 = sum_{r>s>0, r-s odd}
(-1)^r s q^(rs/2), giving the N = 4 decomposition multiplicities 2 q^(-1/8)(-1 + 45 q + 231 q^2 + ...).
Part 3: twining genera for prime-order classes g of M23 with chi_g = 24/(p+1) fixed points:
    H_g = (chi_g/24) H - Ttilde_g/eta^3,  Ttilde_g = (2p/(p+1)) (p E_2(p tau) - E_2(tau))/(p-1) [+ c_g eta(tau)^2 eta(11 tau)^2 for p = 11],
and the moonshine congruences A_n(g) == A_n (mod p) are checked; for 11A the cusp coefficient c is fitted from them.
Usage: python k3_lg_elliptic_genus_mathieu_20261018.py
Session opus-2026-10-17-S43.
"""
import os, sys, json, cmath, math, random
from fractions import Fraction as Fr

# ---------- Part 1: numerical theta functions ----------
def theta(i, tau, z, M=60):
    """Jacobi theta functions; each term computed as one exponential (no overflow for large Im z)."""
    s = 0
    for n in range(-M, M + 1):
        if i in (1, 2):
            r = n + 0.5
        else:
            r = n
        ex = 1j * math.pi * tau * r * r + 2j * math.pi * z * r
        if ex.real < -700:
            continue
        sign = (-1) ** n if i in (1, 4) else 1
        s += sign * cmath.exp(ex)
    return -1j * s if i == 1 else s


def Zw(tau, z, w):
    return theta(1, tau, (1 - w) * z) / theta(1, tau, w * z)


def phi01(tau, z):
    return 4 * sum((theta(i, tau, z) / theta(i, tau, 0)) ** 2 for i in (2, 3, 4))


def lg_orbifold(tau, z, ws, h):
    chat = sum(1 - 2 * w for w in ws)
    tot = 0
    for a in range(h):
        for b in range(h):
            sign = (-1) ** round(chat * (a + b + a * b))
            pre = cmath.exp(2j * math.pi * tau * chat * a * a / 2) * cmath.exp(2j * math.pi * z * chat * a)
            prod = 1
            for w in ws:
                prod *= Zw(tau, z + a * tau + b, w)
            tot += sign * pre * prod
    return tot / h


# ---------- Part 2/3: exact q-series ----------
def series_mul(a, b, N):
    c = [Fr(0)] * (N + 1)
    for i, x in enumerate(a[:N + 1]):
        if x:
            for j, yv in enumerate(b[:N + 1 - i]):
                if yv:
                    c[i + j] += x * yv
    return c


def series_inv(a, N):
    inv = [Fr(0)] * (N + 1)
    inv[0] = Fr(1) / a[0]
    for n in range(1, N + 1):
        s = sum(a[k] * inv[n - k] for k in range(1, n + 1) if k < len(a))
        inv[n] = -s / a[0]
    return inv


def sigma(n):
    return sum(d for d in range(1, n + 1) if n % d == 0)


def E2(N, scale=1):
    e = [Fr(0)] * (N + 1)
    e[0] = Fr(1)
    for n in range(1, N // scale + 1):
        e[n * scale] = Fr(-24 * sigma(n))
    return e


def eta3_over_q18(N):     # eta^3 / q^(1/8) = sum (-1)^n (2n+1) q^(n(n+1)/2)
    e = [Fr(0)] * (N + 1)
    n = 0
    while n * (n + 1) // 2 <= N:
        e[n * (n + 1) // 2] += (-1) ** n * (2 * n + 1)
        n += 1
    return e


def F2(N):
    f = [Fr(0)] * (N + 1)
    for r in range(1, 2 * N + 2):
        for s in range(1, r):
            if (r - s) % 2 == 1 and r * s // 2 <= N:
                f[r * s // 2] += (-1) ** r * s
    return f


def eta_prod_q(shape, N, offset_num):
    """prod eta(a tau)^n as q^(offset) * series; returns series (offset handled by caller)."""
    b = [Fr(0)] * (N + 1)
    b[0] = Fr(1)
    for a, na in shape.items():
        for _ in range(na):
            for n in range(1, N // a + 1):
                st = a * n
                for m in range(N, st - 1, -1):
                    b[m] -= b[m - st]
    return b


def main():
    out = {}
    random.seed(7)
    # Part 1
    lg_checks = []
    for n in (3, 4, 5, 6):
        tau = complex(0.13, 0.9)
        val = Zw(tau, 1e-7, 1.0 / n)
        lg_checks.append((n, val.real))
    print('Part 1a: LG W = x^n, Z(tau, z -> 0) =', [(n, round(v, 6)) for n, v in lg_checks], '(expected n - 1)')
    errs = []
    for t in range(6):
        tau = complex(random.uniform(-0.5, 0.5), random.uniform(0.6, 1.4))
        z = complex(random.uniform(-0.5, 0.5), random.uniform(-0.2, 0.2))
        a = lg_orbifold(tau, z, [0.25] * 4, 4)
        bval = 2 * phi01(tau, z)
        errs.append(abs(a - bval) / abs(bval))
    print('Part 1b: quartic LG orbifold vs 2 phi_{0,1} (K3 elliptic genus): max relative error %.2e over 6 random points' % max(errs))
    out['lg'] = dict(single=lg_checks, quartic_orbifold_rel_err=max(errs))
    # Part 2
    N = 40
    e3 = eta3_over_q18(N)
    inv_e3 = series_inv(e3, N)
    num = [-2 * x for x in E2(N)]
    f2 = F2(N)
    num = [num[i] + 48 * f2[i] for i in range(N + 1)]
    H = series_mul(num, inv_e3, N)          # H * q^(1/8)
    A = [H[n] / 2 for n in range(N + 1)]     # H = 2 q^(-1/8) (A_0 + A_1 q + ...)
    print('Part 2: H(tau) = 2 q^(-1/8) (sum A_n q^n), A_0..A_12 =', [int(x) for x in A[:13]])
    out['A'] = [str(x) for x in A]
    # Part 3: twining genera
    res = {}
    for p, cls in ((2, '2A'), (3, '3A'), (5, '5A'), (7, '7AB'), (11, '11A')):
        chi = Fr(24, p + 1)
        Tt = [Fr(0)] * (N + 1)
        eis = [p * x for x in E2(N, p)]
        e1 = E2(N)
        coef = Fr(2 * p, p + 1) / (p - 1)
        for i in range(N + 1):
            Tt[i] = coef * (eis[i] - e1[i])
        def Hg_from(T):
            t = series_mul(T, inv_e3, N)
            return [chi / 24 * H[i] - t[i] for i in range(N + 1)]
        if p == 11:
            # fit the cusp coefficient c of eta(tau)^2 eta(11 tau)^2 (= q - 2q^2 - ...) from A_n(g) == A_n mod 11, n = 1..N
            cusp = [Fr(0)] + eta_prod_q({1: 2, 11: 2}, N - 1, 1)[:N]
            Hbase = Hg_from(Tt)
            hc = series_mul(cusp, inv_e3, N)
            sols = []
            for num_c in range(-60, 61):
                for den in (1, 5, 6, 12, 10, 30, 60):
                    c = Fr(num_c, den)
                    Hg = [Hbase[i] - c * hc[i] for i in range(N + 1)]
                    Ag = [x / 2 for x in Hg]
                    if all((2 * x).denominator == 1 for x in Ag) and all((int(2 * Ag[n]) - int(2 * A[n])) % 11 == 0 for n in range(1, N + 1)):
                        sols.append(c)
            sols = sorted(set(sols))
            print('   11A: cusp coefficients c (eta^2(tau) eta^2(11 tau)) giving integral A_n(11A) == A_n mod 11 for n <= %d: %s' % (N, [str(s) for s in sols[:12]]))
            best = None
            for c in sols:
                Hg = [Hbase[i] - c * hc[i] for i in range(N + 1)]
                Ag = [x / 2 for x in Hg]
                mx = max(abs(x) for x in Ag[1:])
                if best is None or mx < best[0]:
                    best = (mx, c, Ag)
            if best:
                mx, c, Ag = best
                print('   11A: smallest-coefficient solution c = %s: A_n(11A), n = 0..12: %s' % (c, [str(x) for x in Ag[:13]]))
                res[cls] = dict(c=str(c), A=[str(x) for x in Ag])
            continue
        Hg = Hg_from(Tt)
        Ag = [x / 2 for x in Hg]
        integral = all((2 * x).denominator == 1 for x in Ag)       # traces Tr(g|K_n) = 2 A_n(g) must be integers
        cong = all((int(2 * Ag[n]) - int(2 * A[n])) % p == 0 for n in range(1, N + 1)) if integral else False
        print('   %-4s chi = %s: A_n(g), n = 0..12: %s ; integral %s ; A_n(g) == A_n mod %d for n <= %d: %s'
              % (cls, chi, [str(x) for x in Ag[:13]], integral, p, N, cong))
        res[cls] = dict(A=[str(x) for x in Ag], integral=integral, congruence=cong)

    # 23AB (Paley translations): chi = 1, Ttilde = (23/12) phi_2^(23) + c1 eta^2(tau)eta^2(23 tau) + c2 eta(tau)eta(23 tau) theta_Q0
    p = 23; chi = Fr(1)
    eis = [p * x for x in E2(N, p)]; e1 = E2(N)
    T0 = [Fr(23, 12) / (p - 1) * (eis[i] - e1[i]) for i in range(N + 1)]
    u = [Fr(0)] + [Fr(0)] + eta_prod_q({1: 2, 23: 2}, N - 2, 2)[:N - 1]            # eta^2 eta^2(23): q^2 (1 - 2q ...)
    eta1 = [Fr(0)] + eta_prod_q({1: 1, 23: 1}, N - 1, 1)[:N]                         # eta(tau) eta(23 tau) = q - q^2 - ...
    th = [Fr(0)] * (N + 1)
    for x in range(-12, 13):
        for y in range(-12, 13):
            v = x * x + x * y + 6 * y * y
            if v <= N:
                th[v] += 1
    w = series_mul(eta1, th, N)
    Hb = [chi / 24 * H[i] for i in range(N + 1)]
    tb = series_mul(T0, inv_e3, N); tu = series_mul(u, inv_e3, N); tw = series_mul(w, inv_e3, N)
    sols = []
    for d in (1, 2, 3, 4, 6, 8, 11, 12, 22, 24, 44, 66, 132):
        for n1 in range(-70, 71):
            c1 = Fr(n1, d)
            for n2 in range(-70, 71):
                c2 = Fr(n2, d)
                ok = True
                for n in range(1, N + 1):
                    val = 2 * (Hb[n] - tb[n] - c1 * tu[n] - c2 * tw[n]) / 2
                    if val.denominator != 1 or (int(val) - int(2 * A[n])) % 23:
                        ok = False
                        break
                if ok:
                    sols.append((c1, c2))
    sols = sorted(set(sols))
    print('   23AB: (c1, c2) with integral traces == dim mod 23 for n <= %d: %s' % (N, [(str(a), str(b)) for a, b in sols[:10]]))
    if sols:
        best = min(sols, key=lambda cc: max(abs(Hb[n] - tb[n] - cc[0] * tu[n] - cc[1] * tw[n]) for n in range(1, N + 1)))
        c1, c2 = best
        Ag = [(Hb[n] - tb[n] - c1 * tu[n] - c2 * tw[n]) / 2 for n in range(N + 1)]
        print('   23AB: smallest solution c1 = %s, c2 = %s: A_n(23AB), n = 0..12: %s' % (c1, c2, [str(x) for x in Ag[:13]]))
        res['23AB'] = dict(c1=str(c1), c2=str(c2), A=[str(x) for x in Ag], solutions=len(sols))
    out['twining'] = res
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(out, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'k3_lg_elliptic_genus_mathieu_20261018.json'), 'w'),
              indent=1, default=str)


if __name__ == '__main__':
    main()
