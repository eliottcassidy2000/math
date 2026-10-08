#!/usr/bin/env python3
"""audit H: bracket the optimal fixed-length margin mu* = max_Q min_blocks (tr - 2 lmax)/tr for a few (map, k).
Upper bound mu* < mu: EXACT dual certificate at mu (sum_l c_l [(1-mu) A_l - 2 (A_l z_l)(A_l z_l)^T/(z_l^T A_l z_l)] negative
definite).  Lower bound mu* > mu: a form whose float margin exceeds mu (from the cutting-plane primal).
Usage: python3 f_optimum.py"""
import sys, time
from fractions import Fraction as Fr
import numpy as np
from f_dual_certificates import blocks, psd_sqrt, cutting_plane, exact_cert, margins_at

def bracket(m, k, lo, hi, steps):
    t0 = time.time()
    U = blocks(m, k)
    Uf = U.astype(float)
    An = Uf / np.trace(Uf, axis1=1, axis2=2)[:, None, None]
    Ah, Aih = psd_sqrt(An)
    best_lo = -np.inf
    for _ in range(steps):
        mu = 0.5 * (lo + hi)
        s, best, duals, cuts_A, cuts_w = cutting_plane(An, Ah, Aih, mu=mu)
        if s < 0:
            ok, used, eig = exact_cert(U, cuts_A, cuts_w, duals, mu=Fr(mu).limit_denominator(10 ** 9))
            if ok:
                hi = mu
                continue
            print(f"   mu = {mu:+.5f}: LP says infeasible but exact certificate failed (eigs {eig})", flush=True)
            break
        else:
            mg = margins_at(best[1], U).min()
            best_lo = max(best_lo, mg)
            lo = mu
    print(f"Z_5 m = {m}, k = {k}: {len(U)} blocks; optimal margin in ({max(lo, best_lo):+.5f}, {hi:+.5f}]  "
          f"(upper end exact by dual certificate, lower end float)  [{time.time() - t0:.1f}s]", flush=True)

if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'reading'
    if which == 'reading':
        bracket([1, 2, 3, 7, 1], 2, -0.6, 0.0, 14)
        bracket([1, 2, 3, 7, 1], 3, -0.3, 0.0, 14)
    elif which == 'k4':
        bracket([1, 2, 3, 7, 1], 4, 0.0100, 0.0170, 3)
    elif which == 'pm1':
        for m, c in (([1, 4, 1, 11, 34], -0.0454), ([1, 6, 11, 11, 4], -0.0761), ([1, 1, 6, 39, 11], -0.0705),
                     ([1, 14, 4, 1, 29], -0.0899), ([1, 4, 31, 1, 6], -0.0549)):
            bracket(m, 4, c - 0.004, c + 0.004, 2)
