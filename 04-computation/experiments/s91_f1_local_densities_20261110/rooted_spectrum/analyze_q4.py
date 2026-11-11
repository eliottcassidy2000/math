# analyze_q4.py - window statistics of the pendant move (question Q4).
#
# Pendant move (F2): w loses only to b, H(T + w) = x - s + Yplus at the rooted tournament (T, b).  By converse symmetry
# the pairs (h, n) = (H(T), H(T + w)) are those of the out-degree-one extensions.  All-out-set pairs are (d, x) of the
# rooted census ((T - v, T) with any cut).
# Input: rooted_le9.npz, agg10.npz and the bitsets out/pA_le10.bin.gz, out/pZ_le10.bin.gz of the mac-mini packet
# 05-knowledge/results/tilted_ensemble_local_limit_20261110 (bases of orders 2..10; A = all out-sets, Z = out-degree
# one), h < 1024, (n - h)/2 < 4096.
# Output: stdout (q4.out) and q4_table.txt (one line per odd h <= 400).
import gzip
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..", "..", "..")
HM, KM = 1024, 4096
FAIL = []


def check(name, ok, detail=""):
    print(("ok   " if ok else "FAIL ") + name + ((" : " + detail) if detail else ""), flush=True)
    if not ok:
        FAIL.append(name)


def load_bits(fn):
    b = np.frombuffer(gzip.open(fn, "rb").read(), dtype=np.uint8)
    b = b.reshape(HM // 2, KM // 8)
    return np.unpackbits(b, axis=1, bitorder="little").astype(np.uint8)


def mark(tab, h, n):
    ok = (h < HM) & (n >= h)
    h = h[ok]
    k = (n[ok] - h) >> 1
    ok2 = k < KM
    tab[h[ok2] >> 1, k[ok2]] = 1


def cover(tab, h, lo, hi):
    # number of odd n with lo <= n < hi and (h, n) present, and the size of the window
    tot = 0
    pres = 0
    n = lo if lo % 2 == 1 else lo + 1
    row = tab[h >> 1]
    while n < hi:
        k = (n - h) >> 1
        if 0 <= k < KM:
            tot += 1
            pres += int(row[k])
        n += 2
    return pres, tot


def main():
    z = np.load(os.path.join(HERE, "rooted_le9.npz"))
    a10 = np.load(os.path.join(HERE, "agg10.npz"))
    Z = {}
    A = {}
    zc = np.zeros((HM // 2, KM), np.uint8)
    ac = np.zeros((HM // 2, KM), np.uint8)
    for m in range(1, 10):
        R = z["m%d" % m]
        mark(zc, R[:, 0], R[:, 0] - R[:, 1] + R[:, 3])     # base order m, pendant
        Z[m] = zc.copy()
        mark(ac, R[:, 5], R[:, 0])                         # base order m - 1, all out-sets
        A[m - 1] = ac.copy()
    Z[10] = np.maximum(Z[9], a10["Z"])
    A[9] = np.maximum(A[8], a10["A"])
    pdir = os.path.join(ROOT, "05-knowledge", "results", "tilted_ensemble_local_limit_20261110", "out")
    pA = load_bits(os.path.join(pdir, "pA_le10.bin.gz"))
    pZ = load_bits(os.path.join(pdir, "pZ_le10.bin.gz"))
    check("pendant pairs with bases of order <= 10 = the mac-mini bitset pZ_le10 (out-degree-one extensions), "
          "bit for bit (%d pairs)" % int(pZ.sum()), np.array_equal(Z[10], pZ))
    check("all-out-set pairs with bases of order <= 9 are contained in the mac-mini bitset pA_le10 "
          "(%d of %d pairs)" % (int(A[9].sum()), int(pA.sum())), bool(np.all(A[9] <= pA)))
    check("pendant pairs are all-out-set pairs", bool(np.all(Z[10] <= pA)) and bool(np.all(Z[9] <= A[9])))
    A[10] = pA

    realised = [h for h in range(1, HM, 2) if A[10][h >> 1].any()]
    check("every odd h < 1024 except 7, 21 is realised", [h for h in range(1, HM, 2) if h not in realised] == [7, 21])

    def dens(tab, lo, hi, hmin, hmax):
        pres = tot = 0
        for h in realised:
            if hmin <= h <= hmax:
                p, t = cover(tab, h, -(-lo * h // 1) if isinstance(lo, int) else int(np.ceil(lo * h)), hi * h)
                pres += p
                tot += t
        return pres / tot if tot else float("nan")

    d1 = dens(Z[10], 3, 5, 237, 1023)
    d2 = dens(Z[10], 2, 3, 237, 1023)
    print("reproduction of the packet's line: out-degree-one extensions (h >= 237): density in [3h,5h): %.4f   "
          "in [2h,3h): %.4f" % (d1, d2))
    check("these are 0.7993 and 0.8641", abs(d1 - 0.7993) < 5e-5 and abs(d2 - 0.8641) < 5e-5)

    # table per h <= 400
    lines = []
    full = {8: [], 9: [], 10: []}
    for h in range(1, 401, 2):
        if h in (7, 21):
            continue
        rowtxt = "%3d" % h
        for (lo, hi) in ((1, 3), (3, 5)):
            vals = []
            for tab in (Z[8], Z[9], Z[10], A[8], A[9], A[10]):
                lo_n = h + 2 if lo == 1 else lo * h
                p, t = cover(tab, h, lo_n, hi * h)
                vals.append(p)
            rowtxt += "  | %4d: %4d %4d %4d  %4d %4d %4d" % ((t,) + tuple(vals))
            if lo == 3:
                for i, o in enumerate((8, 9, 10)):
                    if vals[i] == t:
                        full[o].append(h)
        lines.append(rowtxt)
    with open(os.path.join(HERE, "q4_table.txt"), "w") as fo:
        fo.write("# h | window (h, 3h): size, pendant with bases of order <= 8, 9, 10, all out-sets with bases of order "
                 "<= 8, 9, 10 | window [3h, 5h): the same six counts\n")
        fo.write("\n".join(lines) + "\n")
    print("pendant-only windows [3h, 5h) that are complete, h <= 400: bases <= 8: %d values of h (largest %s); "
          "<= 9: %d (largest %s); <= 10: %d (largest %s)"
          % (len(full[8]), max(full[8]) if full[8] else None, len(full[9]), max(full[9]) if full[9] else None,
             len(full[10]), max(full[10]) if full[10] else None))
    print("   complete at order <= 10: %s" % full[10])
    fullA = [h for h in realised if h <= 400 and cover(A[10], h, 3 * h, 5 * h)[0] == cover(A[10], h, 3 * h, 5 * h)[1]]
    print("all-out-set windows [3h, 5h) complete with bases <= 10, h <= 400: %d of %d values of h"
          % (len(fullA), len([h for h in realised if h <= 400])))
    print()
    print("mean coverage (covered odd n / odd n in the window, pooled over h in the band)")
    print("band of h      window        pendant<=8  pendant<=9  pendant<=10   all<=8   all<=9   all<=10")
    for (hmin, hmax) in ((3, 99), (101, 235), (237, 399), (401, 1023)):
        for (lo, hi, name) in ((1.5, 2, "[1.5h, 2h)"), (2, 3, "[2h, 3h)  "), (3, 5, "[3h, 5h)  ")):
            vals = [dens(t, lo, hi, hmin, hmax) for t in (Z[8], Z[9], Z[10], A[8], A[9], A[10])]
            print("[%4d, %4d]   %s    %.4f      %.4f      %.4f       %.4f   %.4f   %.4f"
                  % ((hmin, hmax, name) + tuple(vals)))
    print()
    # the forced gap and the first values
    gap_ok = True
    for h in realised:
        if h >= 5:
            row = Z[10][h >> 1]
            ks = np.nonzero(row)[0]
            ks = ks[ks > 0]
            if len(ks) and 2 * int(ks[0]) < (h + 3) // 2:
                gap_ok = False
    check("no pendant value in (h, (3h + 3)/2) for 5 <= h < 1024 (bases of order <= 10)", gap_ok)
    # which n in the window come only from non-pendant cuts: per-n existential form
    lo_n, hi_n = 1025, 3069
    cnt_can = cnt_any = cnt_anyA = tot = 0
    for n in range(lo_n, hi_n + 1, 2):
        tot += 1
        hc = n // 3
        if hc % 2 == 0:
            hc -= 1
        if Z[10][hc >> 1, (n - hc) >> 1]:
            cnt_can += 1
        found = False
        foundA = False
        for h in range(n // 5 + 1, n // 3 + 1):
            if h % 2 == 1 and h < HM:
                k = (n - h) >> 1
                if k < KM:
                    if Z[10][h >> 1, k]:
                        found = True
                    if A[10][h >> 1, k]:
                        foundA = True
        cnt_any += found
        cnt_anyA += foundA
    print("odd n in [%d, %d] (%d values): n is a pendant extension of a realiser of the largest odd h <= n/3 for "
          "%d of them (%.3f); of some odd h in (n/5, n/3] for %d; an all-out-set extension of some such h for %d"
          % (lo_n, hi_n, tot, cnt_can, cnt_can / tot, cnt_any, cnt_anyA))
    check("every odd n in [1025, 3069] is a pendant extension of a realiser of some odd h in (n/5, n/3]",
          cnt_any == tot)
    # ---- saturated windows: the pendant values of h in [h, 5h) at ALL orders, for h <= 123
    # T = A => C => Z with C the strong component of the root; N = hZ (hA xC - [A empty] s + pi(A) Y) and
    # pi(A) >= (|A| + 1) H(A), Y >= xC, so N < 5h forces A empty or A a single vertex (pi = 2) with Y < 2 xC.
    # Strong tournaments with H <= 123 have order <= 9, so the strong pairs below are complete.
    print()
    strong = {}
    for m in range(1, 10):
        R = z["m%d" % m]
        R = R[(R[:, 6] == 1) & (R[:, 0] <= 123)]
        for x, s, y in set(zip(R[:, 0].tolist(), R[:, 1].tolist(), R[:, 3].tolist())):
            strong.setdefault(x, set()).add((s, y))
    spec = set(h for h in range(1, 124, 2) if h not in (7, 21))
    sat = {}
    for h in sorted(spec):
        P = set()
        for xc in strong:
            if h % xc == 0 and (h // xc) in spec:
                hz = h // xc
                for (s, y) in strong[xc]:
                    n1 = hz * (xc - s + y)
                    if n1 < 5 * h:
                        P.add(n1)
                    if y < 2 * xc:
                        n2 = hz * (xc + 2 * y)
                        if n2 < 5 * h:
                            P.add(n2)
        sat[h] = P
    sub = True
    equal = 0
    covs = []
    for h in sorted(spec):
        cen = set(h + 2 * int(k) for k in np.nonzero(Z[10][h >> 1])[0] if h + 2 * int(k) < 5 * h)
        if not cen <= sat[h]:
            sub = False
        if cen == sat[h]:
            equal += 1
        if h >= 3:
            w = [n for n in range(3 * h, 5 * h, 2)]
            covs.append((len([n for n in w if n in sat[h]]) / len(w), h))
    check("h <= 123: the census pendant values below 5h (bases of order <= 10) lie in the all-order set given by "
          "the composition law; equal for %d of the %d values of h" % (equal, len(spec)), sub)
    best = max(covs)
    worst = min(covs)
    print("exact coverage of [3h, 5h) by pendant extensions of ALL orders, 3 <= h <= 123: mean %.3f, "
          "largest %.3f (h = %d), least %.3f (h = %d); no window is complete: %s"
          % (float(np.mean([c for c, _ in covs])), best[0], best[1], worst[0], worst[1], best[0] < 1))
    check("the pendant-only window law fails for every 3 <= h <= 123 (all orders)", best[0] < 1)
    for h in (13, 25, 45, 75, 105, 123):
        w = [n for n in range(3 * h, 5 * h, 2)]
        print("   h = %3d: pendant values in [3h, 5h) at all orders: %d of %d" % (h, len([n for n in w if n in sat[h]]), len(w)))
    if FAIL:
        print("FAILED:", FAIL)
        sys.exit(1)
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
