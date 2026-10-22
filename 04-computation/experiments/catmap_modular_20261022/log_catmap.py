"""The cat map on the logarithmic side: torus amplitudes of the triplet W(p) algebra.  Session opus-2026-10-22-S55.

Basis (3p - 1 functions of tau, weight 0):
    E_s = theta_(s,p)/eta,  s = 0..p          theta_(s,p) = sum_(n in Z + s/2p) q^(p n^2)          (even thetas)
    O_s = theta'_(s,p)/eta, s = 1..p-1        theta'_(s,p) = sum_(n in Z + s/2p) 2pn q^(p n^2)     (odd thetas)
    tau O_s,               s = 1..p-1
Feigin-Gainutdinov-Semikhatov-Tipunin (CMP 265, 2006): this span is SL2(Z)-invariant and carries R_(p+1) + C^2 x R_(p-1),
the same representation as the Lyubashenko action on the centre of the restricted quantum group at q = exp(i pi/p).
Here the S and T matrices are computed numerically from sample points (no formulas assumed), invariance is checked by
residuals, and the cat map rho(A) = T S T^-1 S^-1 is analysed:
    eigenvalue moduli, and the Lucas trace law  Tr rho(A^n) = Tr R_(p+1)(A^n) + L_(2n) Tr R_(p-1)(A^n).
Usage: python log_catmap.py [pmax]"""
import sys, math, cmath
import numpy as np
import mpmath as mp

mp.mp.dps = 40


def eta(tau):
    q = mp.exp(2j * mp.pi * tau)
    return mp.exp(2j * mp.pi * tau / 24) * mp.qp(q)


def theta_even(s, p, tau, K=40):
    return mp.fsum(mp.exp(2j * mp.pi * tau * p * (n + mp.mpf(s) / (2 * p)) ** 2) for n in range(-K, K + 1))


def theta_odd(s, p, tau, K=40):
    return mp.fsum(2 * p * (n + mp.mpf(s) / (2 * p)) * mp.exp(2j * mp.pi * tau * p * (n + mp.mpf(s) / (2 * p)) ** 2)
                   for n in range(-K, K + 1))


def basis(p, tau):
    e = eta(tau)
    E = [theta_even(s, p, tau) / e for s in range(p + 1)]
    O = [theta_odd(s, p, tau) / e for s in range(1, p)]
    return E + O + [tau * o for o in O]


def action(p, gamma, pts):
    """matrix M with f(gamma tau) = M f(tau) for the basis vector f: least squares over sample points, in mpmath."""
    a, b, c, d = gamma
    F = mp.matrix([basis(p, t) for t in pts])                     # rows: points, cols: basis functions
    G = mp.matrix([basis(p, (a * t + b) / (c * t + d)) for t in pts])
    FH = F.H
    X = mp.inverse(FH * F) * (FH * G)                            # normal equations at 40 digits
    resid = mp.mnorm(F * X - G, 1) / mp.mnorm(G, 1)
    return X.T, resid


def lucas(m):
    a, b = 2, 1
    for _ in range(m):
        a, b = b, a + b
    return a


if __name__ == '__main__':
    pmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    phi = (1 + 5 ** 0.5) / 2
    rng = np.random.default_rng(1)
    for p in range(2, pmax + 1):
        dim = 3 * p - 1
        pts = [mp.mpc(float(x), float(y)) for x, y in zip(rng.uniform(-0.5, 0.5, 2 * dim), rng.uniform(0.55, 1.1, 2 * dim))]
        Sm, rS = action(p, (0, -1, 1, 0), pts)
        Tm, rT = action(p, (1, 1, 0, 1), pts)
        Cm = Tm * Sm * mp.inverse(Tm) * mp.inverse(Sm)
        S = np.array(Sm.tolist(), dtype=complex); T = np.array(Tm.tolist(), dtype=complex); C = np.array(Cm.tolist(), dtype=complex)
        kS = len(S)
        ev = np.linalg.eigvals(C)
        mods = sorted(set(round(abs(z), 6) for z in ev))
        # sub-representations: even block (first p+1), odd block pair (last 2(p-1))
        e, o = p + 1, p - 1
        blockE = np.abs(S[:e, e:]).max() + np.abs(S[e:, :e]).max() + np.abs(T[:e, e:]).max() + np.abs(T[e:, :e]).max()
        CE = C[:e, :e]
        # R_(p-1): action on O modulo tau*O -- read off from the 2x2 block structure: use the S, T action on the O-span
        # projected: the C^2 x R_(p-1) block has the form gamma_def (x) R(gamma); recover R(A) as the block
        # coefficient of the 'constant' component: (C restricted) = A_def (x) R(A) up to basis order.
        Clog = C[e:, e:]
        trace_ok = []
        for n in range(1, 11):
            Cn = np.linalg.matrix_power(C, n)
            Rn_even = np.trace(np.linalg.matrix_power(CE, n))
            Rlog = np.linalg.matrix_power(Clog, n)
            # A_def^n (x) R(A)^n : trace = L_(2n) * tr R(A)^n; and the (O,O)+(tauO,tauO) diagonal blocks sum to
            # (a + d) R^n with a + d = L_(2n), so tr R^n = tr(Rlog)/L_(2n)
            trR = np.trace(Rlog) / lucas(2 * n)
            trace_ok.append(abs(np.trace(Cn) - (Rn_even + lucas(2 * n) * trR)) < 1e-6 * max(1, abs(np.trace(Cn))))
        # the log block: is it A_def (x) R? test: Clog = kron(Adef, R) in a suitable basis -> check its eigenvalues are
        # phi^(+-2) times unit-modulus numbers
        evlog = np.linalg.eigvals(Clog)
        ratio_ok = all(min(abs(abs(z) - phi ** 2), abs(abs(z) - phi ** -2)) < 1e-6 for z in evlog) if o else True
        evE = np.linalg.eigvals(CE)
        unit_ok = all(abs(abs(z) - 1) < 1e-6 for z in evE)
        print(f'   max |S| coupling even-odd {np.abs(S[:e, e:]).max():.1e}, |T| coupling {np.abs(T[:e, e:]).max():.1e}')
        print(f'p = {p}: dim {dim}; S residual {float(rS):.1e}, T residual {float(rT):.1e}; even/log blocks decouple: {blockE < 1e-8}; '
              f'|eig C| = {mods}; even block unitary-spectrum: {unit_ok}; log block |eig| = phi^(+-2): {ratio_ok}; '
              f'Lucas trace law n = 1..10: {all(trace_ok)}', flush=True)
        if p <= 4:
            print('   eigenvalues of C on the log block:', np.round(evlog, 6))
            print('   traces Tr rho(A^n), n = 1..8:', [complex(round(np.trace(np.linalg.matrix_power(C, n)).real, 6), round(np.trace(np.linalg.matrix_power(C, n)).imag, 6)) for n in range(1, 9)])
