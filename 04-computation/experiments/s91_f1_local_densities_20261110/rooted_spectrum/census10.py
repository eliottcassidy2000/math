# census10.py - the rooted census of order 10, aggregated (162 million (base, cut) pairs are not stored).
#
# Bases: the list codes9.npy of census.py (316,153 labelled tournaments of order 9 containing every one of the
# 191,536 classes; verified there against the direct definition and against the class census of H).  For every base
# B and every out-set S of a new root w the rooted state (x, s, e, Yplus, Yminus, d) of (B + w, w) is formed as in
# census.py.  Aggregates written to agg10.npz:
#   sset, eset     start / end counts of order 10
#   ypset, ymset   Yplus / Yminus values of order 10
#   xset           Hamiltonian path counts of order 10 (checked against out/hist_10.txt of the mac-mini packet)
#   Z, Zd          pendant pairs (h, n) = (x, x - s + Yplus) and (x, x - e + Yminus), h < 1024, (n - h)/2 < 4096
#   A              all-out-set pairs (h, n) = (d, x) with bases of order 9
#   pend           the set of pendant values n = x - s + Yplus (n < 2^19)
#   sd             pairs (d, s) with d <= 255 (host defect of a start count)
#   mx             (max Yplus, max pendant value, least H of a strong tournament of order 10)
#   cnt            violation counters of the inequalities listed in NAMES (all must be 0) and row counts
import os
import sys
import time
import numpy as np
from numba import njit
import core

HERE = os.path.dirname(os.path.abspath(__file__))
HM, KM = 1024, 4096
NAMES = ["rows", "parity Y+ = s mod 2", "Y+ >= x", "non-source: 2s <= x - 1", "non-source: 2(Y+ - s) >= x + 1",
         "non-source, x >= 5: 2(Y+ - s) >= x + 3", "s = 0: Y+ >= 2x", "parity Y- = e mod 2", "Y- >= x",
         "non-sink: 2e <= x - 1", "non-sink: 2(Y- - e) >= x + 1", "strong rows", "Y+ - x >= indeg (indeg - 1)",
         "overflow Y", "overflow pend", "s + e <= x", "strong: s >= 1 and e >= 1", "Y+ >= d + e"]


@njit(cache=True)
def run(codes, sset, eset, ypset, ymset, xset, Z, Zd, A, pend, sd, cnt, mx):
    n = 9
    size = 1 << n
    full = size - 1
    f = np.zeros((size, n), np.int64)
    g = np.zeros((size, n), np.int64)
    hS = np.zeros(size, np.int64)
    out = np.zeros(n, np.int64)
    s = np.zeros(n, np.int64)
    e = np.zeros(n, np.int64)
    yp = np.zeros(n, np.int64)
    ym = np.zeros(n, np.int64)
    Q = np.zeros((n, n), np.int64)
    comp = np.zeros(2, np.int64)
    t = np.zeros(8, np.int64)
    YM = ypset.shape[0]
    PM = pend.shape[0]
    for k in range(codes.shape[0]):
        core.out_masks(codes[k], n, out)
        core.dp_tables(n, out, f, g, hS)
        core.base_data(n, f, g, hS, s, e, yp, ym, Q)
        core.end_components(n, out, comp)
        H = hS[full]
        for S in range(size):
            core.cut_tuple(n, S, H, s, e, yp, ym, Q, comp, t)
            x = t[0]
            ss = t[1]
            ee = t[2]
            a = t[3]
            b = t[4]
            d = t[5]
            st = t[6]
            indeg = n - t[7]
            cnt[0] += 1
            sset[ss] = 1
            eset[ee] = 1
            xset[x] = 1
            if a < YM:
                ypset[a] = 1
            else:
                cnt[13] += 1
            if b < YM:
                ymset[b] = 1
            else:
                cnt[13] += 1
            if a > mx[0]:
                mx[0] = a
            if ((a - ss) & 1) != 0:
                cnt[1] += 1
            if a < x:
                cnt[2] += 1
            if S != full:
                if 2 * ss > x - 1:
                    cnt[3] += 1
                if 2 * (a - ss) < x + 1:
                    cnt[4] += 1
                if x >= 5 and 2 * (a - ss) < x + 3:
                    cnt[5] += 1
            if ss == 0 and a < 2 * x:
                cnt[6] += 1
            if ((b - ee) & 1) != 0:
                cnt[7] += 1
            if b < x:
                cnt[8] += 1
            if S != 0:
                if 2 * ee > x - 1:
                    cnt[9] += 1
                if 2 * (b - ee) < x + 1:
                    cnt[10] += 1
            if st == 1:
                cnt[11] += 1
                if x < mx[2]:
                    mx[2] = x
                if ss < 1 or ee < 1:
                    cnt[16] += 1
            if a - x < indeg * (indeg - 1):
                cnt[12] += 1
            if ss + ee > x:
                cnt[15] += 1
            if a < d + ee:
                cnt[17] += 1
            nn = x - ss + a
            if nn < PM:
                pend[nn] = 1
            else:
                cnt[14] += 1
            if nn > mx[1]:
                mx[1] = nn
            if x < HM:
                kk = (nn - x) >> 1
                if kk < KM:
                    Z[x >> 1, kk] = 1
                kk = (b - ee) >> 1
                if kk < KM:
                    Zd[x >> 1, kk] = 1
            if d < HM:
                kk = (x - d) >> 1
                if kk < KM:
                    A[d >> 1, kk] = 1
            if d < 256:
                sd[d >> 1, ss] = 1


def main():
    t0 = time.time()
    codes = np.load(os.path.join(HERE, "codes9.npy"))
    sset = np.zeros(4096, np.uint8)
    eset = np.zeros(4096, np.uint8)
    ypset = np.zeros(1 << 19, np.uint8)
    ymset = np.zeros(1 << 19, np.uint8)
    xset = np.zeros(1 << 15, np.uint8)
    Z = np.zeros((HM // 2, KM), np.uint8)
    Zd = np.zeros((HM // 2, KM), np.uint8)
    A = np.zeros((HM // 2, KM), np.uint8)
    pend = np.zeros(1 << 20, np.uint8)
    sd = np.zeros((128, 256), np.uint8)
    cnt = np.zeros(len(NAMES), np.int64)
    mx = np.zeros(3, np.int64)
    mx[2] = 1 << 40
    run(codes, sset, eset, ypset, ymset, xset, Z, Zd, A, pend, sd, cnt, mx)
    print("bases %d, rows %d, strong rows %d, max Y+ %d, max pendant value %d, least H of a strong tournament of "
          "order 10: %d  (%.1f s)" % (len(codes), cnt[0], cnt[11], mx[0], mx[1], mx[2], time.time() - t0), flush=True)
    fail = []

    def check(name, ok, detail=""):
        print(("ok   " if ok else "FAIL ") + name + ((" : " + detail) if detail else ""), flush=True)
        if not ok:
            fail.append(name)

    check("row count = bases * 512", cnt[0] == len(codes) * 512)
    check("least H of a strong tournament of order 10 is 125", mx[2] == 125, str(mx[2]))
    for i in (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17):
        check("order 10, no violation of: " + NAMES[i], cnt[i] == 0, "violations %d" % cnt[i])
    check("converse symmetry: start set = end set, Y+ set = Y- set, Z = Zd",
          np.array_equal(sset, eset) and np.array_equal(ypset, ymset) and np.array_equal(Z, Zd))
    hist = os.path.join(HERE, "..", "..", "..", "..", "05-knowledge", "results",
                        "tilted_ensemble_local_limit_20261110", "out", "hist_10.txt")
    if os.path.exists(hist):
        hv = sorted(int(l.split()[0]) for l in open(hist) if l.strip())
        mine = [int(v) for v in np.nonzero(xset)[0]]
        check("H-values of order 10: %d values, equal to the support of out/hist_10.txt" % len(mine), hv == mine)
    else:
        check("hist_10.txt of the mac-mini packet found", False)
    np.savez_compressed(os.path.join(HERE, "agg10.npz"), sset=sset, eset=eset, ypset=ypset, ymset=ymset, xset=xset,
                        Z=Z, A=A, pend=pend, sd=sd, cnt=cnt, mx=mx)
    print("saved agg10.npz (%.1f s)" % (time.time() - t0))
    if fail:
        print("FAILED:", fail)
        sys.exit(1)
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
