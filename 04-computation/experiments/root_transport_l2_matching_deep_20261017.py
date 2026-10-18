#!/usr/bin/env python3
"""l2 matching between actual and ideal ROOT-transport coefficients to depth 10, via the banded positive gauge (S41).

Ideal coefficients K_d(r, b) = (-1)^(b+d) 4 2^(d+b) nu_d(r, b) are recovered from clock_duality's positive gauge
(root_transport_energy_ladder_banded_20261017.step, bands b < BMAX), rounded to integers (a check reports the largest
rounding residue), and compared with the actual transported coefficients at every seed of the actual support.
Usage: python root_transport_l2_matching_deep_20261017.py DMAX k1 k2 ...
Session opus-2026-10-17-S41.
"""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from root_transport_energy_profile_20261017 import initial, step
from root_transport_energy_ladder_banded_20261017 import K0, step as nstep
from root_transport_l2_matching_20261017 import bands

BMAX = 120


def main():
    DMAX = int(sys.argv[1]); ks = [int(a) for a in sys.argv[2:]]
    nus = []
    nu = np.zeros((BMAX, 3))
    for b in range(BMAX):
        nu[b] = (-1) ** b * K0[b % 6] / (4.0 * 2.0 ** b)
    q = 3; nus.append(nu)
    for d in range(1, DMAX + 1):
        nu = nstep(nu, q, BMAX); q *= 3; nus.append(nu)
    rows = []
    for k in ks:
        N = 2 ** k
        S, C = initial(N)
        for d in range(DMAX + 1):
            nu = nus[d]; mod = nu.shape[1]
            b = bands(S, N, d)
            ok = (b >= 0) & (b < BMAX)
            Kf = np.zeros(len(S))
            bb = b[ok]; rr = S[ok] % mod
            Kf[ok] = ((-1.0) ** (bb + d)) * 4.0 * (2.0 ** (d + bb)) * nu[bb, rr]
            ideal = np.rint(Kf).astype(np.int64)
            resid = float(np.max(np.abs(Kf - ideal))) if len(S) else 0.0
            err = C - ideal
            nz = np.nonzero(err)[0]
            off_root = [int(S[i]) for i in nz if S[i] != 1]
            e_err = float(np.sum(err.astype(np.float64) ** 2)) / N
            e_off = float(np.sum(err[S != 1].astype(np.float64) ** 2)) / N
            rows.append(dict(k=k, d=d, err_energy=e_err, err_energy_off_root=e_off, n_err=int(len(nz)), seeds=off_root[:20],
                             rounding_residue=resid, beyond_band=int(np.sum(b >= BMAX))))
            print('N = 2^%d d = %2d: ||c_d - K_d||^2/N = %.6f (off ROOT %.6f); seeds with error %d %s; rounding residue %.2e'
                  % (k, d, e_err, e_off, len(nz), off_root[:8], resid), flush=True)
            if d < DMAX:
                S, C = step(S, C)
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(rows, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'root_transport_l2_matching_deep_20261017.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
