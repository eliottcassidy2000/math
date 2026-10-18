#!/usr/bin/env python3
"""Ideal ROOT-transport energy ladder to depth 10+ in clock_duality's positive gauge, with a rigorous band tail (S41).

clock_duality_oct15 (root_positive_transport_20261017g.md, P1) showed that nu_d(r, b) = (-1)^(b+d) K_d(r, b)/(4 2^(d+b))
evolves by the positive probability kernel (r, b) -> (F_A(r), b + A - 1) with probability 2^-A, A >= 1, where
F_A(r) = (3r + c_A)/2^A mod 3q, c_A = (3 + 2^a - 2^A)/3, a = 1 (A odd), 2 (A even). In this gauge

    gamma_d = sum_(r,b) |nu_d(r, b)|,        E_d = 4 2^d sum_b 2^b sum_r nu_d(r, b)^2         (THM-4667 statement 6).

Bands only increase, so truncating at b < BMAX keeps every retained value exact. The omitted energy is bounded by
4 2^d sum_(b >= BMAX) 2^b mu_d(b)^2, where mu_d is the band marginal of the positive evolution of |nu_0| (an exact 1-D
convolution with the increment law P(A - 1 = k) = 2^-(k+1)), because sum_r nu_d(r, b)^2 <= (sum_r |nu_d(r, b)|)^2
<= mu_d(b)^2. The gamma tail is sum_(b >= BMAX) mu_d(b).
Usage: python root_transport_energy_ladder_banded_20261017.py DMAX [BMAX]
Session opus-2026-10-17-S41.
"""
import sys, os, json, time
import numpy as np
from numba import njit

K0 = np.array([[-1, -1, 2], [0, -3, 3], [-1, -1, 2], [-3, 0, 3], [-1, -1, 2], [0, 0, 0]], dtype=np.float64)   # K0[b][r]


@njit(cache=True)
def step(nu, q, BMAX):
    """nu: [BMAX, q] at depth d-1 (residues mod q = 3^d); returns [BMAX, 3q]."""
    Q = 3 * q
    out = np.zeros((BMAX, Q))
    inv2 = (Q + 1) // 2                       # 2^-1 mod Q (Q odd)
    for r in range(q):
        for b in range(BMAX):
            v = nu[b, r]
            if v == 0.0:
                continue
            p = 0.5
            i2A = 1                           # 2^-A mod Q
            pow2A = 1                         # 2^A mod 3Q-safe big? use modular arithmetic below
            for A in range(1, BMAX - b + 1):
                i2A = (i2A * inv2) % Q
                a = 1 if A % 2 == 1 else 2
                # c_A = (3 + 2^a - 2^A)/3 mod Q: compute 2^A mod 3Q, then divide by 3
                pow2A = (pow2A * 2) % (3 * Q)
                num = (3 + (1 << a) - pow2A) % (3 * Q)
                cA = num // 3
                t = ((3 * r + cA) % Q) * i2A % Q
                out[b + A - 1, t] += p * v
                p *= 0.5
    return out


def band_marginal_tail(nu0_abs_band, d, BMAX, LONG=1000):
    mu = np.zeros(LONG); mu[:len(nu0_abs_band)] = nu0_abs_band
    g = 0.5 ** (np.arange(LONG) + 1)
    for _ in range(d):
        mu = np.convolve(mu, g)[:LONG]
    b = np.arange(LONG)
    etail = 4 * 2.0 ** d * np.sum((2.0 ** b[BMAX:]) * mu[BMAX:] ** 2)
    gtail = np.sum(mu[BMAX:])
    return etail, gtail


def main():
    DMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    BMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    # depth 0: nu_0(r, b) = (-1)^b K_0(r, b)/(4 2^b), K_0 period 6 in b
    nu = np.zeros((BMAX, 3))
    for b in range(BMAX):
        nu[b] = (-1) ** b * K0[b % 6] / (4.0 * 2.0 ** b)
    nu0_abs_band = np.abs(nu).sum(axis=1)
    rows = []; q = 3; prevE = None
    t0 = time.time()
    for d in range(DMAX + 1):
        if d > 0:
            nu = step(nu, q, BMAX); q *= 3
        bw = 2.0 ** np.arange(BMAX)
        E = 4 * 2.0 ** d * float(np.sum(bw[:, None] * nu ** 2))
        gam = float(np.abs(nu).sum())
        bal = float(np.abs(nu.reshape(BMAX, 3, q // 3).sum(axis=1)).max())
        etail, gtail = band_marginal_tail(nu0_abs_band, d, BMAX)
        cross = None if prevE is None else E - prevE
        rows.append(dict(d=d, E=E, E_tail_bound=etail, gamma=gam, gamma_tail_bound=gtail, cross=cross, balance=bal))
        print('d = %2d  E_d = %.9f (+ <= %.1e)  gamma_d = %.12f (+ <= %.1e)  cross %s  ratio %s  lift-balance %.1e  [%.0f s]'
              % (d, E, etail, gam, gtail, '%+.6f' % cross if cross is not None else '   -   ',
                 '%.4f' % (cross / rows[-2]['cross']) if d >= 2 and rows[-2]['cross'] else '  -  ', bal, time.time() - t0), flush=True)
        prevE = E
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(dict(BMAX=BMAX, rows=rows), open(os.path.join(here, '..', '..', '05-knowledge', 'results',
                                                            'root_transport_energy_ladder_banded_20261017.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
