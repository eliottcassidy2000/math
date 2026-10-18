#!/usr/bin/env python3
"""The ideal finite-depth ROOT transport envelope gamma_d to high depth, with its l2 energy (S41).

clock_duality_oct15 (root_braid_logistic_20261017f.md, section 3c, R8) proved the exact integer recursion

    K_d(r, b) = sum_(j >= 0, B^j(r) != 0 (3)) (-1)^(a_j) K_(d-1)( (2^(a_j) B^j(r) - 1)/3 mod 3^d, b - (a_j - 1) - 2j ),

with B(x) = 4x + 1, a_j = 2 if B^j(r) = 1 (mod 3) and 1 if = 2, K_0(r, b) = h_r(3/2^b), residue modulus q_d = 3^(d+1),
binary-band period P_d = 2 q_d, and the l1 envelope after d odd transports A^(d)(N) = gamma_d N + O_d(log^(d+1) N),

    gamma_d = sum_(b < P_d) 2^(-b) sum_r |K_d(r, b)| / (4 2^d (1 - 2^(-P_d))).

This script recomputes gamma_d (validated against 16/7, 155896/87381, 1.4488781081, 1.2280701518) and pushes to higher
depth, together with the band-weighted l2 energy

    E_d = sum_(b < P_d) 2^(-b) sum_r K_d(r, b)^2 / (4 2^d (1 - 2^(-P_d)))      (same density normalisation),

the maximal coefficient, and the mean |K| per nonzero entry. If E_d stays bounded while the support spreads over 3^(d+1)
residues, Cauchy-Schwarz forces gamma_d <= sqrt(E_d * (density of support)) ~ 3^(-d/2): an 'inherently controlled' depth
law. Usage: python root_transport_depth_ladder_20261017.py DMAX
Session opus-2026-10-17-S41.
"""
import sys, os, json, time
import numpy as np
from numba import njit

K0 = np.array([[-1, -1, 2], [0, -3, 3], [-1, -1, 2], [-3, 0, 3], [-1, -1, 2], [0, 0, 0]], dtype=np.int64)   # K0[b][r], period 6


@njit(cache=True)
def step(prev, d):
    """prev: array [oldperiod, oldmod] of K_(d-1); returns [period, mod] of K_d."""
    mod = 3 ** (d + 1); period = 2 * mod
    oldmod = mod // 3; oldperiod = 2 * oldmod
    out = np.zeros((period, mod), dtype=np.int64)
    for r in range(mod):
        rho = r
        j = 0
        while 2 * j <= period - 1:
            if rho % 3 != 0:
                a = 2 if rho % 3 == 1 else 1
                s = ((2 ** a * rho - 1) // 3) % oldmod
                sign = 1 if a == 2 else -1
                off = (a - 1) + 2 * j
                # add sign * prev[(b - off) % oldperiod, s] for b >= off
                for b in range(off, period):
                    out[b, r] += sign * prev[(b - off) % oldperiod, s]
            rho = (4 * rho + 1) % mod
            j += 1
    return out


@njit(cache=True)
def stats(K, d):
    period, mod = K.shape
    g = 0.0; e = 0.0; w = 1.0; mx = 0; nz = 0; sabs = 0
    for b in range(period):
        sa = 0.0; sq = 0.0
        for r in range(mod):
            v = K[b, r]
            if v != 0:
                nz += 1
                av = v if v > 0 else -v
                sabs += av
                if av > mx:
                    mx = av
            sa += abs(v); sq += v * v
        g += w * sa; e += w * sq
        w *= 0.5
    norm = 4.0 * 2.0 ** d * (1.0 - 2.0 ** (-period))
    return g / norm, e / norm, mx, nz, sabs


def main():
    DMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    K = K0.copy()
    rows = []
    t0 = time.time()
    for d in range(DMAX + 1):
        if d > 0:
            K = step(K, d)
        g, e, mx, nz, sabs = stats(K, d)
        # last-digit balance check
        mod = 3 ** (d + 1); q = mod // 3
        bal = bool(np.all(K[:, :q] + K[:, q:2 * q] + K[:, 2 * q:] == 0))
        rows.append(dict(d=d, modulus=mod, period=2 * mod, gamma=g, energy=e, max_abs=int(mx), nonzero=int(nz),
                         mean_abs_nonzero=sabs / max(nz, 1), last_digit_balance=bal, seconds=round(time.time() - t0, 1)))
        print('d = %d  modulus %8d  gamma_d = %.10f  E_d = %.6f  max|K| = %6d  nonzero %10d (%.3f)  mean|K| %.3f  balance %s  [%.0f s]'
              % (d, mod, g, e, mx, nz, nz / (mod * 2 * mod), sabs / max(nz, 1), bal, time.time() - t0), flush=True)
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(rows, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'root_transport_depth_ladder_20261017.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
