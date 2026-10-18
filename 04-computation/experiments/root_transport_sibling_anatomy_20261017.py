#!/usr/bin/env python3
"""Sibling anatomy of the actual ROOT-transport energy changes (S41).

Every odd-transport fibre is a union of sibling pairs (s, 2s+1), s = 1 (mod 8), at scales spaced by 64, plus, for
u = 2 (mod 3), the singleton (2u-1)/3 (THM-4667 statement 8). In a pair, s has a = 2 (sign +) and 2s+1 has a = 1
(sign -), so (Tc)(u) = -[u = 2 (3)] c((2u-1)/3) + sum_k (c(s_k) - c(2s_k+1)). Hence

    E_(d+1) - E_d = INTRA_d + INTER_d,   INTRA_d = -2 sum_(s = 1 mod 8) c_d(s) c_d(2s+1)   (sibling covariance),

and INTER_d collects the cross terms between different blocks of one fibre (the singleton and pair differences at
scales 64^k apart). This script measures both, depth by depth, for every N in chosen octaves (octave sums) and for
sampled N, and checks the pair property on the fly.
Usage: python root_transport_sibling_anatomy_20261017.py k [k ...] [--sample K S]
Session opus-2026-10-17-S41.
"""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from root_transport_energy_profile_20261017 import initial, step


def intra(S, C):
    m = S % 8 == 1
    s = S[m]; cs = C[m]
    t = 2 * s + 1
    idx = np.searchsorted(S, t)
    idx[idx >= len(S)] = 0
    hit = S[idx] == t
    return float(-2.0 * np.sum(cs[hit].astype(np.float64) * C[idx[hit]].astype(np.float64)))


def anatomy(N):
    S, C = initial(N)
    E = []; IN = []; CR = []
    while True:
        Cf = C.astype(np.float64)
        E.append(float(Cf @ Cf))
        if len(S) == 0 or (len(S) == 1 and S[0] == 1):
            break
        IN.append(intra(S, C))
        S, C = step(S, C)
    E = np.array(E)
    dE = np.diff(E)
    IN = np.array(IN)
    return E, dE, IN


def main():
    args = sys.argv[1:]
    sample = None
    if '--sample' in args:
        i = args.index('--sample'); sample = (int(args[i + 1]), int(args[i + 2])); args = args[:i]
    out = []
    for k in [int(a) for a in args]:
        X = 2 ** k
        L = 400
        Es = np.zeros(L); dEs = np.zeros(L); INs = np.zeros(L)
        for N in range(X, 2 * X):
            E, dE, IN = anatomy(N)
            n = len(dE)
            Es[:len(E)] += E; Es[len(E):] += E[-1]
            dEs[:n] += dE; INs[:n] += IN
        X2 = float(X) ** 2
        last = int(np.nonzero(dEs)[0].max()) + 1
        print('octave 2^%d (sums / X^2): depth, E_d, dE_d, INTRA (sibling), INTER' % k)
        for d in list(range(0, 12)) + list(range(12, last, 6)):
            print('  d = %3d  E %.4f  dE %+.4f  INTRA %+.4f  INTER %+.4f' % (d, Es[d] / X2, dEs[d] / X2, INs[d] / X2,
                                                                           (dEs[d] - INs[d]) / X2))
        cum_in = np.cumsum(INs[:last]) / X2; cum_inter = np.cumsum(dEs[:last] - INs[:last]) / X2
        dmin = int(np.argmin(Es[:last + 1]))
        print('  total: INTRA %.4f, INTER %.4f; up to the trace minimum (d = %d): INTRA %.4f INTER %.4f; after: INTRA %.4f INTER %.4f'
              % (cum_in[-1], cum_inter[-1], dmin, cum_in[dmin - 1], cum_inter[dmin - 1], cum_in[-1] - cum_in[dmin - 1],
                 cum_inter[-1] - cum_inter[dmin - 1]), flush=True)
        out.append(dict(k=k, E=(Es[:last + 1] / X2).tolist(), dE=(dEs[:last] / X2).tolist(), INTRA=(INs[:last] / X2).tolist()))
    if sample:
        K, SMP = sample
        rng = np.random.default_rng(41)
        for N in [int(x) for x in rng.integers(2 ** K, 2 ** (K + 1), SMP)]:
            E, dE, IN = anatomy(N)
            dmin = int(np.argmin(E))
            print('N = %d: E_0/N %.3f, min E/N %.3f at d = %d, final %.3f; INTRA before/after min %+.3f / %+.3f, INTER %+.3f / %+.3f'
                  % (N, E[0] / N, E[dmin] / N, dmin, E[-1] / N, IN[:dmin].sum() / N, IN[dmin:].sum() / N,
                     (dE[:dmin] - IN[:dmin]).sum() / N, (dE[dmin:] - IN[dmin:]).sum() / N), flush=True)
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(out, open(os.path.join(here, '..', '..', '05-knowledge', 'results',
                                     'root_transport_sibling_anatomy_%s_20261017.json' % '_'.join(args)), 'w'))


if __name__ == '__main__':
    main()
