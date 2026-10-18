#!/usr/bin/env python3
"""Statistics of the frontier plateau: are the frontier coefficients c_Y(u) Gaussian-like and nearly independent? (S42)

Tests P13-P14 of HYP-9290. For Y in {100, 1000}, the incremental frontier vector is sampled at a regular grid of N in
octaves 2^24..2^(KMAX) (every STEP-th N). With s_u = mean c_Y(u)^2/N over the samples:
    model mean m = sum s_u,  model sd = sqrt(2 sum s_u^2),  F_eff = m^2/sum s_u^2,
pooled kurtosis of c_Y(u)/sqrt(N s_u) over the seeds carrying 90% of the energy, the mean absolute correlation between
the normalised coefficients of different seeds, and the observed mean, sd and maximum of E_Y/N over the samples.
Usage: python root_transport_frontier_statistics_20261017.py KMAX
Session opus-2026-10-17-S42.
"""
import sys, os, json, time
import numpy as np
from numba import njit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from root_transport_frontier_energy_big_20261017 import clock_parity8, frontier_map, seed, wfun, thr


@njit(cache=True)
def sampled(Nmax, Y, fu, fs, Nstart, step):
    c = np.zeros(Y + 1, dtype=np.int64)
    ns = (Nmax - Nstart) // step + 1
    S = np.zeros((ns, Y + 1), dtype=np.float64)
    Ns = np.zeros(ns, dtype=np.int64)
    j = 0
    for N in range(1, Nmax):
        for t in range(3):
            if t == 0:
                m = 4 * N - 3; amt = wfun(m)
            elif t == 1:
                m = 4 * N - 1; amt = wfun(m)
            else:
                m = 2 * N - 1; amt = thr(m)
            s = seed(m)
            i = (s - 1) // 2
            c[fu[i]] += fs[i] * amt
        if N >= Nstart and (N - Nstart) % step == 0 and j < ns:
            for u in range(Y + 1):
                S[j, u] = c[u] / np.sqrt(N)
            Ns[j] = N
            j += 1
    return S[:j], Ns[:j]


def main():
    KMAX = int(sys.argv[1])
    Nmax = 2 ** (KMAX + 1); M = 4 * Nmax
    par = clock_parity8(M)
    out = {}
    for Y, step in ((100, 2 ** 12), (1000, 2 ** 14)):
        fu, fs = frontier_map(M, Y)
        S, Ns = sampled(Nmax, Y, fu, fs, 2 ** 24, step)
        seeds = np.array([u for u in range(1, Y + 1) if u % 2 == 1 and u % 8 != 5])
        Z = S[:, seeds]                       # c_Y(u)/sqrt(N)
        su = (Z ** 2).mean(axis=0)
        m = su.sum(); sd_model = np.sqrt(2 * (su ** 2).sum()); Feff = m ** 2 / (su ** 2).sum()
        E = (Z ** 2).sum(axis=1)
        order = np.argsort(-su); cum = np.cumsum(su[order]) / m
        top = order[:int(np.searchsorted(cum, 0.9)) + 1]
        W = Z[:, top] / np.sqrt(su[top])
        kurt = float((W ** 4).mean() / (W ** 2).mean() ** 2)
        means = Z[:, top].mean(axis=0) / np.sqrt(su[top])
        Cm = np.corrcoef(W.T)
        off = Cm[~np.eye(len(top), dtype=bool)]
        res = dict(Y=Y, samples=int(len(Ns)), F_Y=int(len(seeds)), F_eff=float(Feff), model_mean=float(m),
                   model_sd=float(sd_model), observed_mean=float(E.mean()), observed_sd=float(E.std()),
                   observed_max=float(E.max()), top90_seeds=int(len(top)), pooled_kurtosis=kurt,
                   mean_abs_offdiag_corr=float(np.abs(off).mean()), max_abs_offdiag_corr=float(np.abs(off).max()),
                   max_abs_normalised_mean=float(np.abs(means).max()),
                   top_seeds=[[int(seeds[k]), float(su[k] / m)] for k in order[:10]])
        out[Y] = res
        print('Y = %d: %d samples (octaves 2^24..2^%d); F_Y = %d, F_eff = %.1f' % (Y, len(Ns), KMAX, len(seeds), Feff))
        print('   model mean %.4f sd %.4f | observed mean %.4f sd %.4f max %.4f' % (m, sd_model, E.mean(), E.std(), E.max()))
        print('   seeds carrying 90%% of the energy: %d; pooled kurtosis %.3f; mean |corr| %.3f (max %.3f); max |mean| %.3f'
              % (len(top), kurt, np.abs(off).mean(), np.abs(off).max(), np.abs(means).max()))
        print('   top seeds (share): ' + ', '.join('%d (%.3f)' % (seeds[k], su[k] / m) for k in order[:10]), flush=True)
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump({str(k): v for k, v in out.items()}, open(os.path.join(here, '..', '..', '05-knowledge', 'results',
              'root_transport_frontier_statistics_20261017.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
