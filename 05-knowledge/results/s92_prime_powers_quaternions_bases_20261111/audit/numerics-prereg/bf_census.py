# Auditor numerics-prereg: independent BRUTE-FORCE census of F1 (no symmetry, no sharp loop bounds).
# usage: python3 bf_census.py K [KLO]
#   cumulative: all ordered triples (h1,h2,s) with N <= 2^K
#   window:     2^KLO <= N < 2^K
# Only the trivial bound N >= 2^s y1 y2 (x >= 1) is used to stop the loops; every pair is tested exactly.
import sys
import json
import time
import numpy as np
import numba as nb
from aud_common import psi_table, gcd_nb, MODULI, LCM, fold


@nb.njit
def brute(psi, A, X, L):
    hist = np.zeros((2, L), np.int64)
    typ = np.zeros((2, L + 1), np.int64)
    diag = np.zeros((2, L + 1), np.int64)
    tot = np.zeros(2, np.int64)
    s = 0
    while (1 << s) <= X + 1:
        e = 1 if s == 0 else 0
        sh = 0 if s == 0 else 1
        h1 = 0
        while True:
            y1 = 2 * h1 + 1
            if (y1 << s) > X + 1:
                break
            x1 = 2 * psi[h1] + 1
            h2 = 0
            while True:
                y2 = 2 * h2 + 1
                if (y1 << s) * y2 > X + 1:
                    break
                x2 = 2 * psi[h2] + 1
                N = x1 * x2 + (y1 << s) * y2 - e
                if N <= X and N >= A:
                    hist[sh, N % L] += 1
                    tot[sh] += 1
                    if h1 == h2:
                        diag[sh, gcd_nb(gcd_nb(x1, y1), L)] += 1
                    else:
                        hm = h1 if h1 < h2 else h2
                        xm = 2 * psi[hm] + 1
                        ym = 2 * hm + 1
                        typ[sh, gcd_nb(gcd_nb(xm, ym), L)] += 1
                h2 += 1
            h1 += 1
        s += 1
    return hist, typ, diag, tot


def main():
    K = int(sys.argv[1])
    X = 1 << K
    A = 0
    tag = "cum_%d" % K
    if len(sys.argv) > 2:
        KLO = int(sys.argv[2])
        A = 1 << KLO
        X = X - 1
        tag = "win_%d_%d" % (K, KLO)
    t0 = time.time()
    psi = psi_table(X // 2 + 2)
    hist, typ, diag, tot = brute(psi, A, X, LCM)
    out = {"K": K, "A": A, "X": X, "T": [int(tot[0]), int(tot[1])], "moduli": {}}
    for m in MODULI:
        cls, t, dg = fold(hist, typ, diag, m)
        out["moduli"][str(m)] = {"cls": cls, "typ": [{str(k): v for k, v in t[0].items()}, {str(k): v for k, v in t[1].items()}],
                                 "diag": [{str(k): v for k, v in dg[0].items()}, {str(k): v for k, v in dg[1].items()}]}
    out["diag_total"] = [int(diag[0].sum()), int(diag[1].sum())]
    with open("bf_%s.json" % tag, "w") as f:
        json.dump(out, f)
    print("brute force %s: triples sheet 0: %d, sheet 1: %d; diagonal %d / %d  (%.1f s)" % (
        tag, tot[0], tot[1], out["diag_total"][0], out["diag_total"][1], time.time() - t0))


if __name__ == "__main__":
    main()
