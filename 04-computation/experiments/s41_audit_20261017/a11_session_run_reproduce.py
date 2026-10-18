#!/usr/bin/env python3
"""S41 audit A11: reproduce the committed octave statistics of root_transport_incremental_trace_22_20261017.json for
every octave up to 2^KMAX by calling the session's own numba kernels (build_tree, clock_parity, run) from this driver,
with all outputs written here (the session script would overwrite its committed receipts). This is a reproducibility
check of the session code, not an independent algorithm (A3 and A10 are the independent checks).
Usage: python a11_session_run_reproduce.py [KMAX]
"""
import os, sys, json, time
sys.dont_write_bytecode = True
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(os.environ.get('TEMP', '/tmp'), 's41_audit_numba_cache'))
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.dirname(HERE)
RES = os.path.join(EXP, '..', '..', '05-knowledge', 'results')
sys.path.insert(0, EXP)
from root_transport_incremental_trace_20261017 import build_tree, clock_parity, run


def main():
    KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 21
    t0 = time.time()
    Nmax = 2 ** (KMAX + 1)
    M = 4 * Nmax
    level0, succ, sgn, sizes = build_tree(M)
    DMAX = len(sizes) + 1
    print('tree: %d sources, %d nodes, DMAX %d [%.0f s]' % (len(level0), len(succ), DMAX, time.time() - t0), flush=True)
    clk = clock_parity(M)
    f = 1 - 2 * (clk & 1)
    n = np.arange(M + 1)
    wv = np.where(n % 3 == 2, 2, -1)
    wv[0] = 0
    H = np.cumsum(wv * f)
    del clk, n, wv, f
    T, supE, dsup, climb, E6, Rf, shellmax = run(Nmax, KMAX, level0, succ, sgn, DMAX)
    Ns = np.arange(1, Nmax)
    Rdir = H[4 * Ns] - H[Ns]
    ok = bool(np.array_equal(Rdir, Rf[1:]))
    print('run done [%.0f s]; invariant at every N < 2^%d: %s' % (time.time() - t0, KMAX + 1, ok), flush=True)
    committed = {r['k']: r for r in json.load(open(os.path.join(RES, 'root_transport_incremental_trace_22_20261017.json')))['octaves']}
    rows = []
    worst = 0.0
    for k in range(1, KMAX + 1):
        X = 2 ** k
        prof = T[k] / X ** 2
        sl = slice(X, 2 * X)
        DX = float(np.sum(Rdir[X - 1:2 * X - 1].astype(np.float64) ** 2)) / X ** 2
        row = dict(k=k, T0=float(prof[0]), Tmax=float(prof.max()), d_Tmax=int(prof.argmax()), Tfinal=float(prof[-1]), DX=DX,
                   sup_pointwise=float(supE[sl].max()), N_at_sup=int(X + np.argmax(supE[sl])),
                   frac_sup_at_0=float(np.mean(dsup[sl] == 0)), E6_mean=float(E6[sl].mean()),
                   shell_max_sum=float(np.sum(shellmax[sl] * Ns[X - 1:2 * X - 1])) / X ** 2)
        c = committed[k]
        diffs = {key: abs(row[key] - c[key]) for key in row if key != 'k'}
        worst = max(worst, max(diffs.values()))
        rows.append(row)
        print('2^%2d: T0 %.4f Tmax %.4f (d %d) final %.4f = D_X %.4f; sup %.4f at %d; peak at 0: %.6f; E6 mean %.5f; shell %.4f; '
              'max |diff vs committed| %.1e' % (k, row['T0'], row['Tmax'], row['d_Tmax'], row['Tfinal'], DX,
                                                row['sup_pointwise'], row['N_at_sup'], row['frac_sup_at_0'], row['E6_mean'],
                                                row['shell_max_sum'], max(diffs.values())), flush=True)
    g = int(np.argmax(supE[1:])) + 1
    top = np.argsort(-supE[1:])[:5] + 1
    print('global: max over 1 <= N < 2^%d of sup_d E_d(N)/N = %.6f at N = %d; top five %s; worst diff vs committed %.1e'
          % (KMAX + 1, supE[g], g, [(int(t), round(float(supE[t]), 4)) for t in top], worst))
    json.dump(dict(KMAX=KMAX, invariant_ok=ok, rows=rows, global_max=dict(N=g, value=float(supE[g])), worst_diff=worst),
              open(os.path.join(HERE, 'a11_session_run_reproduce.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
