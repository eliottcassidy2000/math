# Auditor recompute-5128, tasks 2, 3, 6: exact enumeration of all F1 triples (h1, h2, s) with N < 2^B,
#   N = x1 x2 + 2^s y1 y2 - [s = 0],  ordered triples, by sheet (s = 0 / s >= 1) and dyadic bin of N.
# Everything is exact integer arithmetic (int64; N < 2^29).
# Outputs (saved to census_B.npz), all indexed [sheet or s, bin] with bin b <=> 2^(b-1) <= N < 2^b:
#   tot[sheet, b]           ordered triples
#   tots[s, b]              ordered triples with padding s
#   diag[sheet, b]          triples with h1 = h2
#   hist[sheet, b, js, c]   ordered triples with N = c (mod small_p[js])
#   dist[sheet, b, j]       ordered triples with M = N + [s = 0] = 0 (mod all_p[j])   (the distinguished class)
#   Dsm, Dlg, Dboth[sheet, b, j]  ordered triples with h1 != h2 whose smaller / larger / both states are p-divisible
#   Ddiag[sheet, b, j]      triples with h1 = h2 and the state p-divisible
#   viol[j]                 triples with a p-divisible state NOT in the distinguished class (must be 0)
#   trk[t, s, b]            ordered triples with FIRST state tracked[t] and padding s
import sys
import time
import numpy as np
from numba import njit
from concurrent.futures import ThreadPoolExecutor
from a0_defs import psi_table32

SMALL_P = np.array([3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43], dtype=np.int64)
BIG_P = np.array([61, 73, 127, 251, 257, 263, 677, 691, 1039, 2729, 2731, 2741, 8191, 65521, 65537, 65539], dtype=np.int64)
ALL_P = np.concatenate([SMALL_P, BIG_P])
TRACKED = np.array([0, 1, 2, 3, 4, 8, 16, 32, 38, 64, 128, 218, 256, 1024], dtype=np.int64)
TRK_MAX = 1024


@njit(cache=True)
def div_mask(psi, all_p):
    n = psi.shape[0]
    mask = np.zeros(n, dtype=np.uint32)
    for j in range(all_p.shape[0]):
        p = all_p[j]
        h = (p - 1) // 2  # y = 2h + 1 = 0 (mod p)
        while h < n:
            x = 2 * np.int64(psi[h]) + 1
            if x % p == 0:
                mask[h] |= np.uint32(1) << np.uint32(j)
            h += p
    return mask


@njit(cache=True, nogil=True)
def census_part(psi, mask, B, small_p, all_p, trkidx, ntrk, res, nres):
    ns = small_p.shape[0]
    na = all_p.shape[0]
    pmax = 0
    for j in range(ns):
        if small_p[j] > pmax:
            pmax = small_p[j]
    tot = np.zeros((2, B + 1), dtype=np.int64)
    tots = np.zeros((B + 1, B + 1), dtype=np.int64)
    diag = np.zeros((2, B + 1), dtype=np.int64)
    hist = np.zeros((2, B + 1, ns, pmax), dtype=np.int64)
    dist = np.zeros((2, B + 1, na), dtype=np.int64)
    Dsm = np.zeros((2, B + 1, na), dtype=np.int64)
    Dlg = np.zeros((2, B + 1, na), dtype=np.int64)
    Dboth = np.zeros((2, B + 1, na), dtype=np.int64)
    Ddiag = np.zeros((2, B + 1, na), dtype=np.int64)
    viol = np.zeros(na, dtype=np.int64)
    trk = np.zeros((ntrk, B + 1, B + 1), dtype=np.int64)
    thr = np.zeros(B + 2, dtype=np.int64)
    for i in range(B + 2):
        thr[i] = np.int64(1) << i
    hsize = psi.shape[0]
    for s in range(0, B):
        sheet = 0 if s == 0 else 1
        corr = 1 if s == 0 else 0
        Xs = (np.int64(1) << B) - 1 + corr  # M <= Xs  <=>  N <= 2^B - 1
        h1 = np.int64(0)
        while True:
            y1 = 2 * h1 + 1
            # crude monotone bound: a + 2b > y1 (1/2 + 2^(s+1))  =>  y2 < 4 Xs / ((1 + 2^(s+2)) y1)
            y2crude = (4 * Xs) // ((1 + (np.int64(1) << (s + 2))) * y1)
            if (y2crude - 1) // 2 < h1:
                break
            if h1 % nres != res:
                h1 += 1
                continue
            x1 = 2 * np.int64(psi[h1]) + 1
            a = x1
            b = y1 << s
            y2max = (2 * Xs) // (a + 2 * b)
            h2max = (y2max - 1) // 2
            if h2max >= hsize:
                h2max = hsize - 1  # cannot happen if hsize is chosen correctly; checked by the caller
            m1 = mask[h1]
            t1 = -1
            if h1 <= 1024:
                t1 = trkidx[h1]
            h2 = h1
            while h2 <= h2max:
                x2 = 2 * np.int64(psi[h2]) + 1
                y2 = 2 * h2 + 1
                M = a * x2 + b * y2
                if M <= Xs:
                    N = M - corr
                    bb = B
                    while N < thr[bb - 1]:
                        bb -= 1
                    w = 2 if h2 > h1 else 1
                    tot[sheet, bb] += w
                    tots[s, bb] += w
                    if h2 == h1:
                        diag[sheet, bb] += 1
                    for js in range(ns):
                        hist[sheet, bb, js, N % small_p[js]] += w
                    for j in range(ns, na):
                        if M % all_p[j] == 0:
                            dist[sheet, bb, j] += w
                    m2 = mask[h2]
                    if (m1 | m2) != 0:
                        for j in range(na):
                            d1 = (m1 >> np.uint32(j)) & np.uint32(1)
                            d2 = (m2 >> np.uint32(j)) & np.uint32(1)
                            if d1 != 0 or d2 != 0:
                                if M % all_p[j] != 0:
                                    viol[j] += 1
                                if h2 > h1:
                                    if d1 != 0:
                                        Dsm[sheet, bb, j] += 2
                                    if d2 != 0:
                                        Dlg[sheet, bb, j] += 2
                                    if d1 != 0 and d2 != 0:
                                        Dboth[sheet, bb, j] += 2
                                else:
                                    Ddiag[sheet, bb, j] += 1
                    if t1 >= 0:
                        trk[t1, s, bb] += 1
                    if h2 != h1 and h2 <= 1024:
                        t2 = trkidx[h2]
                        if t2 >= 0:
                            trk[t2, s, bb] += 1
                h2 += 1
            h1 += 1
    # distinguished class for the small primes from the histogram
    for sheet in range(2):
        for bb in range(B + 1):
            for js in range(ns):
                p = small_p[js]
                c = (p - 1) if sheet == 0 else 0
                dist[sheet, bb, js] = hist[sheet, bb, js, c]
    return tot, tots, diag, hist, dist, Dsm, Dlg, Dboth, Ddiag, viol, trk


def run(B, nthreads=3, verbose=True):
    t0 = time.time()
    hsize = ((1 << B) + 1) // 3 + 4
    psi = psi_table32(hsize)
    mask = div_mask(psi, ALL_P)
    trkidx = -np.ones(TRK_MAX + 1, dtype=np.int64)
    for i, h in enumerate(TRACKED):
        trkidx[h] = i
    if verbose:
        print(f"# B = {B}: tables for h < {hsize} built in {time.time() - t0:.1f} s", flush=True)
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=nthreads) as ex:
        futs = [ex.submit(census_part, psi, mask, B, SMALL_P, ALL_P, trkidx, len(TRACKED), r, nthreads) for r in range(nthreads)]
        parts = [f.result() for f in futs]
    out = [sum(p[i] for p in parts) for i in range(len(parts[0]))]
    if verbose:
        print(f"# enumeration done in {time.time() - t0:.1f} s", flush=True)
    names = ["tot", "tots", "diag", "hist", "dist", "Dsm", "Dlg", "Dboth", "Ddiag", "viol", "trk"]
    return dict(zip(names, out))


def brute(B):
    """Pure-python reference for small B: all ordered triples, independent psi recursion, no pruning shortcuts."""
    from a0_defs import psi_py
    X = (1 << B) - 1
    hmax = X  # generous
    xs = [2 * psi_py(h) + 1 for h in range(hmax + 1)]
    tot = np.zeros((2, B + 1), dtype=np.int64)
    hist = {}
    Dsm = np.zeros((2, B + 1, len(ALL_P)), dtype=np.int64)
    first = np.zeros((2, B + 1, len(ALL_P)), dtype=np.int64)
    for s in range(0, B + 1):
        sheet = 0 if s == 0 else 1
        for h1 in range(hmax + 1):
            y1 = 2 * h1 + 1
            if (y1 << s) > X + 1:
                break
            for h2 in range(hmax + 1):
                y2 = 2 * h2 + 1
                if (y1 << s) * y2 > X + 1:
                    break
                N = xs[h1] * xs[h2] + (y1 << s) * y2 - (1 if s == 0 else 0)
                if N > X:
                    continue
                bb = N.bit_length()
                tot[sheet, bb] += 1
                for js, p in enumerate(SMALL_P):
                    key = (sheet, bb, js, N % int(p))
                    hist[key] = hist.get(key, 0) + 1
                for j, p in enumerate(ALL_P):
                    p = int(p)
                    lo = min(h1, h2)
                    if h1 != h2 and xs[lo] % p == 0 and (2 * lo + 1) % p == 0:
                        Dsm[sheet, bb, j] += 1
                    if xs[h1] % p == 0 and y1 % p == 0:
                        first[sheet, bb, j] += 1
    return tot, hist, Dsm, first


def main():
    B = int(sys.argv[1]) if len(sys.argv) > 1 else 24
    nthreads = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    if len(sys.argv) > 3 and sys.argv[3] == "selftest":
        # compare with the pure-python reference
        Bs = 13
        r = run(Bs, nthreads=nthreads, verbose=False)
        tot, hist, Dsm, first = brute(Bs)
        assert np.array_equal(tot, r["tot"]), (tot, r["tot"])
        for (sheet, bb, js, c), v in hist.items():
            assert r["hist"][sheet, bb, js, c] == v
        assert int(r["hist"].sum()) == sum(hist.values())
        assert np.array_equal(Dsm, r["Dsm"])
        myfirst = (r["Dsm"] + r["Dlg"]) // 2 + r["Ddiag"]
        assert np.array_equal(first, myfirst)
        assert int(r["viol"].sum()) == 0
        print(f"selftest OK: numba census == pure-python reference for N < 2^{Bs} (totals {tot.sum(axis=1)}, histograms, D_small, first-state counts)")
        return
    r = run(B, nthreads=nthreads)
    np.savez_compressed(f"census_{B}.npz", small_p=SMALL_P, all_p=ALL_P, tracked=TRACKED, **r)
    tot = r["tot"]
    print(f"# violations (triples with a p-divisible state outside the distinguished class): {int(r['viol'].sum())}")
    print("# cumulative ordered triples with N < 2^b:  b  sheet0  sheet1")
    c0 = np.cumsum(tot[0])
    c1 = np.cumsum(tot[1])
    for b in range(10, B + 1):
        print(f"  {b:2d}  {c0[b]:12d}  {c1[b]:12d}")


if __name__ == "__main__":
    main()
