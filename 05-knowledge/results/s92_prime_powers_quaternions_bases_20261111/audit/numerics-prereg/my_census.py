# Auditor numerics-prereg: independent FAST census of F1 (written before reading odd_modulus_census.py).
# usage: python3 -B my_census.py K [KLO] [nthreads]
#   cumulative: ordered triples (h1,h2,s) with N <= 2^K;  window: 2^KLO <= N < 2^K.
# Method: unordered pairs h1 <= h2 (weight 2 off the diagonal), one histogram of N modulo
# LCM = 3^4 5^3 7^2 per sheet, one histogram of gcd(g(h1), LCM) for the smaller state; all moduli
# are folded from these afterwards.  Loop bounds: y/2 < x <= y, i.e. (y+1)/2 <= x <= y.
import sys
import json
import time
import threading
import numpy as np
import numba as nb
from aud_common import psi_table, gcd_nb, MODULI, LCM, fold


@nb.njit(nogil=True)
def run_batch(psi, tasks, lo_i, hi_i, A, X, L, hist, typ, diag, tot):
    for t in range(lo_i, hi_i):
        s = tasks[t, 0]
        h1 = tasks[t, 1]
        lo = tasks[t, 2]
        hi = tasks[t, 3]
        e = 1 if s == 0 else 0
        sh = 0 if s == 0 else 1
        x1 = 2 * psi[h1] + 1
        y1 = 2 * h1 + 1
        b = y1 << s
        gL = gcd_nb(gcd_nb(x1, y1), L)
        npairs = 0
        for h2 in range(lo, hi + 1):
            N = x1 * (2 * psi[h2] + 1) + b * (2 * h2 + 1) - e
            if N <= X and N >= A:
                if h2 > h1:
                    hist[sh, N % L] += 2
                    npairs += 1
                else:
                    hist[sh, N % L] += 1
                    diag[sh, gL] += 1
                    tot[sh] += 1
        typ[sh, gL] += 2 * npairs
        tot[sh] += 2 * npairs


def build_tasks(psi, A, X, chunk):
    tasks = []
    Hcap = len(psi) - 1
    s = 0
    while True:
        e = 1 if s == 0 else 0
        if 1 + (1 << s) - e > X:
            break
        h1 = 0
        while True:
            x1 = 2 * int(psi[h1]) + 1
            y1 = 2 * h1 + 1
            b = y1 << s
            # N >= x1 (y2+1)/2 + b y2 - e  =>  y2 <= (2 (X + e) - x1) / (x1 + 2 b)
            ymax = (2 * (X + e) - x1) // (x1 + 2 * b)
            hi = (ymax - 1) // 2 + 1
            if hi > Hcap:
                hi = Hcap
            # N <= y2 (x1 + b) - e      =>  y2 >= (A + e) / (x1 + b)
            ymin = (A + e) // (x1 + b)
            lo = max(h1, (ymin - 1) // 2 - 1)
            # x is not monotone in h: stop only when even the smallest possible x1 = (y1+1)/2
            # leaves no partner y2 >= y1, i.e. when y1 ((y1+1)/2 + 2 b) > 2 (X + e) - (y1+1)/2
            xmin = (y1 + 1) // 2
            if y1 * (xmin + 2 * b) > 2 * (X + e) - xmin:
                break
            if ymax >= y1 and hi >= lo:
                a = lo
                while a <= hi:
                    bnd = min(hi, a + chunk - 1)
                    tasks.append((s, h1, a, bnd))
                    a = bnd + 1
            h1 += 1
        s += 1
    return np.array(tasks, dtype=np.int64)


def main():
    K = int(sys.argv[1])
    X = 1 << K
    A = 0
    tag = "cum_%d" % K
    nthreads = 8
    if len(sys.argv) > 2 and int(sys.argv[2]) > 0:
        KLO = int(sys.argv[2])
        A = 1 << KLO
        X = X - 1
        tag = "win_%d_%d" % (K, KLO)
    if len(sys.argv) > 3:
        nthreads = int(sys.argv[3])
    t0 = time.time()
    Hcap = (2 * (X + 1)) // 6 + 4
    psi = psi_table(Hcap)
    chunk = 1 << 20
    tasks = build_tasks(psi, A, X, chunk)
    work = tasks[:, 3] - tasks[:, 2] + 1
    order = np.argsort(-work, kind="stable")
    tasks = tasks[order]
    work = work[order]
    # batches of about `chunk` work
    bounds = [0]
    acc = 0
    for i in range(len(tasks)):
        acc += int(work[i])
        if acc >= chunk:
            bounds.append(i + 1)
            acc = 0
    if bounds[-1] != len(tasks):
        bounds.append(len(tasks))
    nb_batches = len(bounds) - 1
    lock = threading.Lock()
    nxt = [0]
    L = LCM
    hists = [np.zeros((2, L), np.int64) for _ in range(nthreads)]
    typs = [np.zeros((2, L + 1), np.int64) for _ in range(nthreads)]
    diags = [np.zeros((2, L + 1), np.int64) for _ in range(nthreads)]
    tots = [np.zeros(2, np.int64) for _ in range(nthreads)]
    # compile
    run_batch(psi, tasks, 0, 0, A, X, L, hists[0], typs[0], diags[0], tots[0])

    def worker(k):
        while True:
            with lock:
                i = nxt[0]
                nxt[0] += 1
            if i >= nb_batches:
                return
            run_batch(psi, tasks, bounds[i], bounds[i + 1], A, X, L, hists[k], typs[k], diags[k], tots[k])

    ths = [threading.Thread(target=worker, args=(k,)) for k in range(nthreads)]
    for th in ths:
        th.start()
    for th in ths:
        th.join()
    hist = sum(hists)
    typ = sum(typs)
    diag = sum(diags)
    tot = sum(tots)
    assert int(hist[0].sum()) == int(tot[0]) and int(hist[1].sum()) == int(tot[1])
    out = {"K": K, "A": A, "X": X, "T": [int(tot[0]), int(tot[1])], "moduli": {}}
    for m in MODULI:
        cls, t, dg = fold(hist, typ, diag, m)
        out["moduli"][str(m)] = {"cls": cls, "typ": [{str(k): v for k, v in t[0].items()}, {str(k): v for k, v in t[1].items()}],
                                 "diag": [{str(k): v for k, v in dg[0].items()}, {str(k): v for k, v in dg[1].items()}]}
    out["diag_total"] = [int(diag[0].sum()), int(diag[1].sum())]
    out["work_tested"] = int(work.sum())
    with open("my_%s.json" % tag, "w") as f:
        json.dump(out, f)
    np.save  # (no binary files are written; everything is in the json)
    print("my census %s: triples sheet 0: %d, sheet 1: %d; diagonal %d / %d; pairs tested %d; tasks %d  (%.1f s)" % (
        tag, tot[0], tot[1], out["diag_total"][0], out["diag_total"][1], int(work.sum()), len(tasks), time.time() - t0))


if __name__ == "__main__":
    main()
