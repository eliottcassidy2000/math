#!/usr/bin/env python3
"""l2 matching error between the actual ROOT transport and clock_duality's ideal coefficients (S41).

clock_duality_oct15's R8 gives ideal coefficients K_d(r, b): the actual coefficient of a seed t after d odd transports
equals K_d(t mod 3^(d+1), b(t)) up to threshold effects, where M_d/2^(b+1) < t/N <= M_d/2^b, M_d = 4 (3/2)^d, and K_d is
periodic in b >= 0 with period 2 3^(d+1). Their P3 bounds the l1 discrepancy uniformly for d <= m/16. This script
measures the l2 discrepancy ||c_d(N) - K_d||^2 / N for d <= DMAX (integer arithmetic for the bands), the item "matching
errors add in quadrature" of the merge plan in root_transport_energy_20261017.md, together with ||c_d||^2/N and the
ideal energy.
Usage: python root_transport_l2_matching_20261017.py DMAX k1 k2 ...   (N = 2^k)
Session opus-2026-10-17-S41.
"""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from root_transport_energy_profile_20261017 import initial, step
from root_transport_depth_ladder_20261017 import K0, step as kstep


def bands(S, N, d):
    """b = largest integer with 2^b 2^d S <= 4 3^d N (b may be negative: then the seed lies above the first band)."""
    num = 4 * 3 ** d * N
    out = np.empty(len(S), dtype=np.int64)
    for i, s in enumerate(S.tolist()):
        den = (2 ** d) * s
        if den > num:
            out[i] = -1; continue
        b = (num // den).bit_length() - 1
        out[i] = b
    return out


def main():
    DMAX = int(sys.argv[1]); ks = [int(a) for a in sys.argv[2:]]
    Ks = [K0.copy()]
    for d in range(1, DMAX + 1):
        Ks.append(kstep(Ks[-1], d))
    rows = []
    for k in ks:
        N = 2 ** k
        S, C = initial(N)
        for d in range(DMAX + 1):
            K = Ks[d]; period, mod = K.shape
            b = bands(S, N, d)
            ideal = np.zeros(len(S), dtype=np.int64)
            ok = b >= 0
            ideal[ok] = K[b[ok] % period, S[ok] % mod]
            # seeds absent from the actual support but with nonzero ideal coefficient are not enumerated here;
            # count them separately through the ideal energy comparison below
            err = C - ideal
            e_err = float(np.sum(err.astype(np.float64) ** 2)) / N
            e_act = float(np.sum(C.astype(np.float64) ** 2)) / N
            e_id = float(np.sum(ideal.astype(np.float64) ** 2)) / N
            l1err = float(np.sum(np.abs(err))) / N
            nz = np.nonzero(err)[0]
            off_root = [int(S[i]) for i in nz if S[i] != 1]
            e_err_off = float(np.sum(err[S != 1].astype(np.float64) ** 2)) / N
            rows.append(dict(k=k, d=d, err_energy=e_err, err_energy_off_root=e_err_off, actual_energy=e_act,
                             ideal_energy_on_support=e_id, l1_err=l1err, support=int(len(S)), n_err=int(len(nz)),
                             err_seeds_off_root=off_root[:20]))
            print('N = 2^%d d = %d: ||c_d - K_d||^2/N = %.6f (off ROOT %.6f)  seeds with error: %d %s  (||c_d||^2/N = %.4f), l1/N = %.6f'
                  % (k, d, e_err, e_err_off, len(nz), off_root[:6], e_act, l1err), flush=True)
            if d < DMAX:
                S, C = step(S, C)
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(rows, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'root_transport_l2_matching_20261017.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
