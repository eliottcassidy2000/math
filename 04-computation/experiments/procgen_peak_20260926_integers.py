#!/usr/bin/env python3
"""procgen_peak_20260926_integers -- the peak-catch modification on actual integers (T4 of the peak lane).

For T_q and horizon L: Bad^act = {n >= 2 : T^j(n) >= n for 1 <= j <= L}; P(n) = max_{0<=j<=L-1} T^j(n)
(the value at the first maximiser); E = {P(n) : n in Bad^act}; G = 1 on E, G = T_q elsewhere.
We build E from all sources n <= N (this gives E cap [1, N] exactly, since P(n) >= n), verify that every
2 <= n <= N descends within L under G, and measure |E cap [1,N]|.  Also the exceptional set
X = Bad^act minus the bad residue classes (actual trajectory bad, word not bad), and the number of
sources whose peak is <= N (before identification of shared peak points)."""
import numpy as np


def _step(x, q):
    """one step of T_q, in place-friendly: returns new array."""
    y = x >> 1
    odd = (x & 1).astype(bool)
    y[odd] = (q * x[odd] + 1) >> 1
    return y, odd


def run_integers(q, L, N, fl, M, chunk=250_000):
    lim = (2 ** 62) // q
    peak_all = np.empty(N + 1, dtype=np.int64)
    bad_all = np.zeros(N + 1, dtype=bool)
    exc, wrong = [], 0
    for a in range(0, N + 1, chunk):
        n = np.arange(a, min(N + 1, a + chunk), dtype=np.int64)
        x = n.copy()
        peak = n.copy()
        desc = np.zeros(len(n), dtype=bool)
        ecount = np.zeros(len(n), dtype=np.int64)
        wordbad = np.ones(len(n), dtype=bool)
        for j in range(1, L + 1):
            assert int(x.max()) < lim, "int64 overflow risk"
            x, odd = _step(x, q)
            ecount += odd
            desc |= x < n
            wordbad &= ecount >= int(fl[M + j]) + 1          # exact: q^e_j > 2^j
            if j <= L - 1:
                np.maximum(peak, x, out=peak)
        bad = ~desc & (n >= 2)
        wb = wordbad & (n >= 2)
        exc += n[bad & ~wb].tolist()                          # actual-bad, word not bad
        wrong += int(np.sum(wb & ~bad))                       # must be 0: bad words have only bad lifts
        peak_all[a:a + len(n)] = peak
        bad_all[a:a + len(n)] = bad
    E = np.unique(peak_all[bad_all])
    nsrc_peak_le_N = int(np.sum(bad_all & (peak_all <= N)))
    nbad = int(bad_all.sum())
    del peak_all, bad_all
    # verification of G = 1 on E, T elsewhere
    alldesc = True
    for a in range(2, N + 1, chunk):
        start = np.arange(a, min(N + 1, a + chunk), dtype=np.int64)
        xs = start.copy()
        done = np.zeros(len(xs), dtype=bool)
        for j in range(1, L + 1):
            pos = np.searchsorted(E, xs)
            pos[pos >= len(E)] = len(E) - 1
            inE = E[pos] == xs
            xs, _ = _step(xs, q)
            xs[inE] = 1
            done |= xs < start
        alldesc &= bool(done.all())
    return {
        "nbad": nbad, "nE": int(np.sum(E <= N)), "all_descend": alldesc,
        "exc": exc, "wrong": wrong, "nsrc_peak_le_N": nsrc_peak_le_N, "maxE": int(E.max()) if len(E) else 0,
    }
