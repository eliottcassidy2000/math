# analyze_q1.py - rooted Chan-Pak (question Q1): which integers are start counts s(T, v), which are Yplus(T, v).
#
# Input: rooted_le9.npz (census.py), agg10.npz (census10.py), reps_le8.npz, codes9.npy.
# Universe: every rooted tournament of order 1..10.
# Output (stdout): value sets per order, complete initial intervals, least orders for m <= 200, the single-base
# constructions, the host defect table.
import os
import sys
import numpy as np
from numba import njit
import core

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = []
MAXH = [1, 1, 1, 3, 5, 15, 45, 189, 661, 3357, 15745]      # max H at orders 0..10


def check(name, ok, detail=""):
    print(("ok   " if ok else "FAIL ") + name + ((" : " + detail) if detail else ""), flush=True)
    if not ok:
        FAIL.append(name)


def runs(d, lo, hi):
    # d: value -> order; compact run-length print
    out = []
    v = lo
    while v <= hi:
        w = v
        while w + 1 <= hi and d.get(w + 1) == d.get(v):
            w += 1
        out.append(("%d" % v if v == w else "%d-%d" % (v, w)) + ":" + str(d.get(v)))
        v = w + 1
    return "  ".join(out)


@njit(cache=True)
def complete_scan(n, codes, res):
    # res[0] = max H(B) with complete start vector, res[1] = its code, res[2] = max number of distinct subset sums,
    # res[3] = number of bases with complete start vector
    size = 1 << n
    f = np.zeros((size, n), np.int64)
    g = np.zeros((size, n), np.int64)
    hS = np.zeros(size, np.int64)
    out = np.zeros(n, np.int64)
    sv = np.zeros(n, np.int64)
    res[0] = 0
    res[1] = -1
    res[2] = 0
    res[3] = 0
    for k in range(codes.shape[0]):
        core.out_masks(codes[k], n, out)
        core.dp_tables(n, out, f, g, hS)
        H = hS[size - 1]
        for v in range(n):
            sv[v] = g[size - 1, v]
        srt = np.sort(sv)
        pref = 0
        ok = True
        for v in range(n):
            if srt[v] > pref + 1:
                ok = False
                break
            pref += srt[v]
        if ok:
            res[3] += 1
            if H > res[0]:
                res[0] = H
                res[1] = codes[k]
        seen = np.zeros(H + 1, np.uint8)
        cntd = 0
        for S in range(size):
            t = 0
            for v in range(n):
                if (S >> v) & 1:
                    t += sv[v]
            if seen[t] == 0:
                seen[t] = 1
                cntd += 1
        if cntd > res[2]:
            res[2] = cntd


def near_transitive(n):
    # transitive on 0..n-1 with the arc 0 -> n-1 reversed
    out = np.zeros(n, np.int64)
    for i in range(n):
        for j in range(i + 1, n):
            out[i] |= (1 << j)
    out[0] &= ~(1 << (n - 1))
    out[n - 1] |= 1
    return out


def vectors(n, out):
    size = 1 << n
    f = np.zeros((size, n), np.int64)
    g = np.zeros((size, n), np.int64)
    hS = np.zeros(size, np.int64)
    core.dp_tables(n, out, f, g, hS)
    s = np.zeros(n, np.int64); e = np.zeros(n, np.int64); yp = np.zeros(n, np.int64); ym = np.zeros(n, np.int64)
    Q = np.zeros((n, n), np.int64)
    core.base_data(n, f, g, hS, s, e, yp, ym, Q)
    return int(hS[size - 1]), s, e, yp, ym


def subset_sums(vec):
    sums = {0}
    for a in vec:
        sums |= {x + int(a) for x in sums}
    return sums


def main():
    z = np.load(os.path.join(HERE, "rooted_le9.npz"))
    a10 = np.load(os.path.join(HERE, "agg10.npz"))
    S = {}
    Y = {}
    for m in range(1, 10):
        R = z["m%d" % m]
        S[m] = set(R[:, 1].tolist())
        Y[m] = set(R[:, 3].tolist())
    S[10] = set(np.nonzero(a10["sset"])[0].tolist())
    Y[10] = set(np.nonzero(a10["ypset"])[0].tolist())

    print("== start counts s(T, v) by order of T ==")
    I = {}
    for m in range(1, 11):
        mx = max(S[m])
        i = 0
        while i + 1 in S[m]:
            i += 1
        I[m] = i
        miss = [v for v in range(0, mx) if v not in S[m]]
        print("order %2d: %5d values, max %5d, every integer of [%d, %d] occurs, %4d missing below the max, first: %s"
              % (m, len(S[m]), mx, 0 if 0 in S[m] else 1, i, len(miss), miss[:14]))
        check("order %d: largest start count = max H at order %d" % (m, m - 1), mx == MAXH[m - 1])
    print("ratio (complete interval) / (max H at order n - 1), n = 6..10: "
          + ", ".join("%.3f" % (I[m] / MAXH[m - 1]) for m in range(6, 11)))
    print("binary single-base bound 2^(n-3) + 1 at order n (THM-5001 family): "
          + ", ".join("%d:%d" % (m, 2 ** (m - 3) + 1) for m in range(5, 11)))
    print()
    print("== Yplus(T, v) by order of T ==")
    J = {}
    for m in range(1, 11):
        mx = max(Y[m])
        i = 0
        while i + 1 in Y[m]:
            i += 1
        J[m] = i
        miss = [v for v in range(1, mx) if v not in Y[m]]
        print("order %2d: %6d values, max %6d, every integer of [1, %d] occurs, %5d missing below the max, first: %s"
              % (m, len(Y[m]), mx, i, len(miss), miss[:14]))
        check("order %d: Yplus is never 0 and [1, 2^(n-1)] occurs (one-apex roots)" % m,
              0 not in Y[m] and i >= 2 ** (m - 1))
    print()
    fs = {}
    fy = {}
    for v in range(0, 201):
        for m in range(1, 11):
            if v in S[m] and v not in fs:
                fs[v] = m
            if v in Y[m] and v not in fy:
                fy[v] = m
    check("every m in [0, 200] is a start count and every m in [1, 200] a Yplus value at order <= 9",
          all(v in fs and fs[v] <= 9 for v in range(0, 201)) and all(v in fy and fy[v] <= 9 for v in range(1, 201)))
    print("least order of T with a vertex starting exactly m Hamiltonian paths (m:order, runs), m = 0..200:")
    print("   " + runs(fs, 0, 200))
    print("least order of T with Yplus(T, v) = m, m = 1..200:")
    print("   " + runs(fy, 1, 200))
    # lower bounds: s(T, v) <= max H at order n - 1
    okl = all(v <= MAXH[fs[v] - 1] for v in range(1, 201))
    check("least order of a start count m is at least the least n with max H(n - 1) >= m; equality for %d of the "
          "200 values m = 1..200" % sum(1 for v in range(1, 201)
                                        if fs[v] == min(n for n in range(1, 11) if MAXH[n - 1] >= v)), okl)
    print()

    print("== single-base constructions ==")
    for n in range(4, 13):
        out = near_transitive(n)
        H, s, e, yp, ym = vectors(n, out)
        ss = subset_sums(s)
        es = subset_sums(e)
        ok = (H == 2 ** (n - 2) + 1 and sorted(e.tolist()) == sorted([1, 1, 1] + [2 ** i for i in range(1, n - 2)])
              and ss == set(range(0, H + 1)) and es == ss)
        check("near-transitive base of order %d: H = 2^(n-2) + 1 = %d, end vector (1, 1, 2, ..., 2^(n-3), 1), "
              "subset sums of the start vector = [0, %d]" % (n, H, H), ok)
    for m in range(1, 11):
        out = np.zeros(m, np.int64)
        for i in range(m):
            for j in range(i + 1, m):
                out[i] |= (1 << j)
        H, s, e, yp, ym = vectors(m, out)
        ok = (H == 1 and yp.tolist() == [2 ** i for i in range(m)] and ym.tolist() == [2 ** (m - 1 - i) for i in range(m)]
              and {1 + x for x in subset_sums(yp)} == set(range(1, 2 ** m + 1)))
        check("transitive base of order %d: Yplus(TT, v_i) = 2^i, so Yplus(TT + w, w) runs over [1, 2^%d]" % (m, m), ok)
    reps = np.load(os.path.join(HERE, "reps_le8.npz"))
    codes9 = np.load(os.path.join(HERE, "codes9.npy"))
    print("largest H(B) over bases B of order m whose start vector has subset sums [0, H(B)] (single-base optimum),")
    print("against the near-transitive value 2^(m-2) + 1 and the trivial bound 2^m - 1:")
    res = np.zeros(4, np.int64)
    best = {}
    for m in range(1, 10):
        codes = reps["n%d" % m] if m <= 8 else codes9
        complete_scan(m, codes, res)
        best[m] = int(res[0])
        print("   order %d: %4d   (2^(m-2)+1 = %3d, 2^m - 1 = %3d; most distinct subset sums of one base: %d)"
              % (m, res[0], 2 ** (m - 2) + 1 if m >= 2 else 1, 2 ** m - 1, res[2]))
    check("single-base optimum never exceeds 2^m - 1", all(best[m] <= 2 ** m - 1 for m in best))
    print()

    print("== host defect: delta(m) = min { H(T - v) - m : s(T, v) = m } ==")
    INF = 10 ** 9
    mu = {}
    for m in range(1, 10):
        R = z["m%d" % m]
        for s, d in set(zip(R[:, 1].tolist(), R[:, 5].tolist())):
            if d < mu.get(s, INF):
                mu[s] = d
    sd = a10["sd"]
    for dd in range(128):
        for s in np.nonzero(sd[dd])[0]:
            if 2 * dd + 1 < mu.get(int(s), INF):
                mu[int(s)] = 2 * dd + 1
    LIM = 122
    delta = {m: mu[m] - m for m in range(0, LIM + 1)}
    print("delta(m) for m = 0..%d: 0 for odd m, 1 for even m, except" % LIM)
    exc = {m: d for m, d in delta.items() if d != (m + 1) % 2}
    print("   " + ", ".join("delta(%d) = %d" % (m, d) for m, d in sorted(exc.items())))
    check("delta(m) <= 3 for every m <= %d" % LIM, max(delta.values()) <= 3)
    check("odd m <= %d: delta = 0 except delta(7) = delta(21) = 2" % LIM,
          {m for m in exc if m % 2 == 1} == {7, 21} and exc[7] == 2 and exc[21] == 2)
    even_exc = sorted(m for m in exc if m % 2 == 0)
    check("even m <= %d: delta = 1 except delta = 3 at %s" % (LIM, even_exc),
          even_exc == [6, 12, 20, 40, 42, 54, 66, 76, 82] and all(exc[m] == 3 for m in even_exc))
    # the even exceptions are the m with m + 1 a path count none of whose realisers has a vertex starting exactly one
    # path; realisers of m + 1 <= 123 reduce to strong tournaments of order <= 9, all inside the census
    strong_s1 = set()
    for m in range(1, 10):
        R = z["m%d" % m]
        R = R[(R[:, 6] == 1) & (R[:, 1] == 1)]
        strong_s1 |= set(R[:, 0].tolist())
    strong_x = set()
    for m in range(1, 10):
        R = z["m%d" % m]
        strong_x |= set(R[R[:, 6] == 1][:, 0].tolist())
    not_strong = sorted(x for x in range(1, 124, 2) if x not in strong_x)
    no_s1 = sorted(x for x in strong_x if x <= 123 and x not in strong_s1)
    print("odd x <= 123 that are not H of a strong tournament (order <= 9; orders >= 10 have H >= 125): %s"
          % not_strong)
    print("H-values x <= 123 of strong tournaments none of which has a vertex starting exactly one path: %s" % no_s1)
    check("even exceptions = {x - 1 : x odd <= 123, no strong realiser of x has a vertex with s = 1}",
          sorted(x - 1 for x in range(3, 124, 2) if x not in strong_s1) == even_exc)
    if FAIL:
        print("FAILED:", FAIL)
        sys.exit(1)
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
