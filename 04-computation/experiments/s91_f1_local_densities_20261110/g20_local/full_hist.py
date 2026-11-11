#!/usr/bin/env python3
"""S91 helper (g20_local), part B(i): exact residue distribution of H over ALL 4^n words of the two-apex grammar G(2,0).

For each n the histogram of H modulo 720720 = 16 * 9 * 5 * 7 * 11 * 13 over all 4^n words is computed by a depth-first
walk of the 4-ary word tree with the 36-state transfer of g20_core.py (exact int64 counts; max H is printed and is far
below 2^63).  From it (full_hist_report.py): class factors  p * share(H = c mod p)  for p = 3, 5, 7, 11, 13 (and 9),
and  (M/2) * share  over the odd classes modulo M = 8, 16 (H is odd, Redei).

Two methods:
  direct   (n <= NDIRECT): all 4^n words;
  reduced  (n >= 3): only the 9 * 4^(n-2) words with first letter != 0 and last letter != 3, and
           N_n = N'_n + 2 N_(n-1) - N_(n-2)   (a leading 0 is a source and a trailing 3 is a sink: H unchanged;
           inclusion-exclusion over "first letter 0" / "last letter 3").
For 3 <= n <= NDIRECT both are computed and compared (check).

Usage: python3 full_hist.py NMAX [threads] [NDIRECT]   (writes full_hist_n{n}.json; existing files are reused)
Checks: histogram total = 4^n; no even value; n <= 6 against the brute-force subset DP; direct = reduced.
"""
import itertools
import json
import os
import sys
import time

import numpy as np

MOD = 720720


def main():
    nmax = int(sys.argv[1])
    threads = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    ndirect = int(sys.argv[3]) if len(sys.argv) > 3 else 12
    os.environ["NUMBA_NUM_THREADS"] = str(threads)
    import numba
    import g20_core as g
    numba.set_num_threads(threads)
    here = os.path.dirname(os.path.abspath(__file__))
    ok = True
    hists = {}
    for n in range(0, nmax + 1):
        path = os.path.join(here, f"full_hist_n{n}.json")
        if os.path.exists(path):
            d = json.load(open(path))
            h = np.zeros(MOD, dtype=np.int64)
            h[np.array(d["residues"], dtype=np.int64)] = np.array(d["counts"], dtype=np.int64)
            hists[n] = h
            print(f"n = {n}: loaded {path}")
            continue
        t0 = time.time()
        method = []
        hist = None
        mx = 0
        if n >= 3:
            hred, mx = g.nb_reduced_hist(n, MOD, g.T_TO, g.T_CNT, g.ACC, g.INIT)
            hist = hred + 2 * hists[n - 1] - hists[n - 2]
            method.append("reduced")
        if n <= ndirect:
            hdir, mxd = g.nb_full_hist(n, MOD, min(n, 2), g.T_TO, g.T_CNT, g.ACC, g.INIT)
            mx = max(int(mx), int(mxd))
            if hist is not None:
                same = bool((hist == hdir).all())
                ok &= same
                if not same:
                    print(f"n = {n}: direct and reduced histograms DIFFER")
            hist = hdir
            method.append("direct")
        dt = time.time() - t0
        tot = int(hist.sum())
        good = tot == 4 ** n and int(hist[0::2].sum()) == 0 and int(hist.min()) >= 0
        if n <= 6:
            ref = np.zeros(MOD, dtype=np.int64)
            for w in itertools.product(range(4), repeat=n):
                ref[int(g.nb_bruteforce(np.array(w, dtype=np.int64))) % MOD] += 1
            good &= bool((ref == hist).all())
        ok &= good
        hists[n] = hist
        nz = np.nonzero(hist)[0]
        out = {"n": n, "mod": MOD, "words": tot, "max_H_seen": int(mx), "seconds": round(dt, 2), "method": method,
               "residues": [int(r) for r in nz], "counts": [int(hist[r]) for r in nz]}
        with open(path, "w") as fh:
            json.dump(out, fh)
        print(f"n = {n}: {tot} words, max H = {int(mx)}, {dt:.1f} s, {'+'.join(method)}, "
              f"checks {'ok' if good else 'FAILED'}", flush=True)
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
