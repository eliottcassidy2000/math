#!/usr/bin/env python3
"""Octave-summed transport energy: T_d(X) = sum_(X <= N < 2X) E_d(N), the trace of the transported energy kernel (S41).

At depth 0 it is sum_N sum_s c_0(s)^2 (clock_duality's kernel G_X on the diagonal, aggregated on seeds); at full depth it
is the current energy D_X = sum_N R(N)^2 exactly. The script computes the whole profile d -> T_d(X)/X^2 for complete
octaves, checks T_infinity = D_X against the direct current, and reports the maximum of the profile.
Usage: python root_transport_octave_trace_20261017.py k1 k2 ...
Session opus-2026-10-17-S41.
"""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from root_transport_energy_profile_20261017 import profile


def main():
    out = []
    for k in [int(a) for a in sys.argv[1:]]:
        X = 2 ** k
        prof = None; DX = 0
        for N in range(X, 2 * X):
            R, E, L1 = profile(N)
            DX += R * R
            E = np.array(E) * N                     # back to absolute energy
            if prof is None:
                prof = E.copy()
            else:
                if len(E) > len(prof):
                    prof = np.concatenate([prof, np.full(len(E) - len(prof), prof[-1])])
                elif len(E) < len(prof):
                    E = np.concatenate([E, np.full(len(prof) - len(E), E[-1])])
                prof += E
        T = prof / X ** 2
        print('X = 2^%d: T_0/X^2 = %.4f, T_6/X^2 = %.4f, max_d T_d/X^2 = %.4f at d = %d, final T/X^2 = %.4f, D_X/X^2 = %.4f'
              % (k, T[0], T[6], T.max(), int(T.argmax()), T[-1], DX / X ** 2), flush=True)
        print('   profile (every 4th depth):', ' '.join('%.2f' % t for t in T[::4]))
        out.append(dict(k=k, profile=T.tolist(), D_X=DX))
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(out, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'root_transport_octave_trace_%s_20261017.json' % '_'.join(sys.argv[1:])), 'w'))


if __name__ == '__main__':
    main()
