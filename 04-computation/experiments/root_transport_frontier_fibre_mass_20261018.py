#!/usr/bin/env python3
"""Mass versus coherence on the Solfeggio fibres of the frontier (S42).

preprint_connections_checks_20261017.out found that frontier seeds u = 2 (mod 3) carry about half their seed share of
the frontier energy E_Y(N). This script separates mass from coherence: for each fibre r = u mod 3 of the frontier seeds
it reports (i) the share of sources landing there (count of odd m <= 4N whose first seed <= Y is in the fibre), (ii) the
share of the diagonal energy sum k_N(m)^2 over those sources (the energy if no two contributions interfered), and
(iii) the share of the actual frontier energy sum c_Y(u)^2, together with the coherence ratio (iii)/(ii) per fibre.
Usage: python root_transport_frontier_fibre_mass_20261018.py "N1,N2,..." Y
Session opus-2026-10-17-S42.
"""
import sys, os, json
import numpy as np
from numba import njit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from root_transport_frontier_energy_big_20261017 import frontier_map, seed, wfun


@njit(cache=True)
def landing(N, Y, fu, fs):
    cnt = np.zeros(Y + 1, dtype=np.int64)
    diag = np.zeros(Y + 1, dtype=np.int64)
    c = np.zeros(Y + 1, dtype=np.int64)
    for m in range(1, 4 * N, 2):
        r = m % 3
        if m < 2 * N:
            k = -3 if r == 1 else (3 if r == 2 else 0)
        else:
            k = wfun(m)
        s = seed(m)
        i = (s - 1) // 2
        u = fu[i]
        cnt[u] += 1
        diag[u] += k * k
        c[u] += fs[i] * k
    return cnt, diag, c


def main():
    Ns = [int(a) for a in sys.argv[1].split(',')]; Y = int(sys.argv[2])
    M = 4 * max(Ns) + 8
    fu, fs = frontier_map(M, Y)
    seeds = np.array([u for u in range(1, Y + 1) if u % 2 == 1 and u % 8 != 5])
    out = []
    for N in Ns:
        cnt, diag, c = landing(N, Y, fu, fs)
        E = float((c[seeds].astype(np.float64) ** 2).sum()); D = float(diag[seeds].sum()); C = float(cnt[seeds].sum())
        row = dict(N=N, Y=Y, E_over_N=E / N, diag_over_N=D / N)
        line = 'N = %d, Y = %d: E_Y/N = %.4f, diagonal/N = %.4f |' % (N, Y, E / N, D / N)
        for r in range(3):
            sel = seeds[seeds % 3 == r]
            sc = cnt[sel].sum() / C; sd = diag[sel].sum() / D; se = float((c[sel].astype(np.float64) ** 2).sum()) / E
            row['fibre%d' % r] = dict(seed_share=len(sel) / len(seeds), landing_share=float(sc), diag_share=float(sd),
                                       energy_share=se, coherence=se / sd)
            line += '  fibre %d: seeds %.3f landing %.3f diag %.3f energy %.3f (ratio %.2f)' % (
                r, len(sel) / len(seeds), sc, sd, se, se / sd)
        print(line, flush=True)
        out.append(row)
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(out, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'root_transport_frontier_fibre_mass_20261018.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
