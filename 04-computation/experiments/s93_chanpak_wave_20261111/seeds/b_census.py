#!/usr/bin/env python3
"""S93 lane "seeds", deliverable (b): census of seeds (r = +1) and anti-seeds (r = -1).

For every tournament B of order n <= 9 (one representative per isomorphism class, cls_<n>.txt of b_classes.py) and
every tournament of order 10 (all min-score one-vertex extensions of the order-9 classes; every class of order 10
occurs at least once), and EVERY vertex set O:
    r = H(B) - W_B(O),  M = U_B(O);  seed: r = 1, anti-seed: r = -1.
Recorded:
    minv[sign][c]  least 2-adic valuation v of a slope M = c 2^v (c odd < 2^21),
    the number of unit pairs by |O| (labelled subsets of the representatives; order 10: of the extensions),
    one certificate (least v, then least order) for every odd c <= 65535 and sign.
Usage:
    python3 b_census.py small            orders 1..9, writes b_census_small.npz
    python3 b_census.py part P K         order 10, representatives with index = P mod K, writes b_census_part_P_K.npz
    python3 b_census.py merge K          tables from b_census_small.npz and the K parts (if present): b_census.out
"""
import os
import sys
import time

import numpy as np
from numba import njit, prange

from seedlib import dp_tables

CN = 1 << 20            # odd c < 2^21, index c >> 1
NCERT = 1 << 15         # certificates for odd c <= 65535
NCH = 48

CTZ = np.zeros(1 << 11, np.int64)
for _i in range(11):
    CTZ[1 << _i] = _i


@njit(cache=True)
def popcount(x):
    c = 0
    while x:
        x &= x - 1
        c += 1
    return c


@njit(cache=True)
def decode(n, code, out):
    for v in range(n):
        out[v] = 0
    b = 0
    for j in range(1, n):
        for i in range(j):
            if (code >> b) & 1:
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
            b += 1


@njit(cache=True)
def dp_range(n, out, fe, gs, hS, lo):
    """tables for the subsets S in [lo, 2^n) (lo = 1: all; lo = 2^(n-1): those containing the vertex n - 1)."""
    full = (1 << n) - 1
    if lo <= 1:
        hS[0] = 1
        for v in range(n):
            fe[0, v] = 0
            gs[0, v] = 0
        lo = 1
    for S in range(lo, full + 1):
        for v in range(n):
            fe[S, v] = 0
            gs[S, v] = 0
        if (S & (S - 1)) == 0:
            v = CTZ[S]
            fe[S, v] = 1
            gs[S, v] = 1
            hS[S] = 1
            continue
        t = 0
        rem = S
        while rem:
            lb = rem & (-rem)
            rem ^= lb
            v = CTZ[lb]
            R = S ^ lb
            ov = out[v]
            a = 0
            x = R & ~ov
            while x:
                l2 = x & (-x)
                x ^= l2
                a += fe[R, CTZ[l2]]
            b = 0
            x = R & ov
            while x:
                l2 = x & (-x)
                x ^= l2
                b += gs[R, CTZ[l2]]
            fe[S, v] = a
            gs[S, v] = b
            t += a
        hS[S] = t


@njit(cache=True)
def wu_tables(n, fe, gs, hS, Wm, U):
    """Wm[e, f] = W_ef = Q(e, f) + Q(f, e) (symmetric), U[o]."""
    full = (1 << n) - 1
    top = 1 << (n - 1)
    for e in range(n):
        U[e] = 0
        for f in range(n):
            Wm[e, f] = 0
    for S in range(0, full):
        C = full ^ S
        h = hS[S]
        x = C
        while x:
            lb = x & (-x)
            x ^= lb
            f = CTZ[lb]
            g = gs[C, f]
            U[f] += h * g
        if S & top:                    # each unordered partition {S, C} once (S contains the top vertex, C nonempty)
            x = C
            while x:
                lb = x & (-x)
                x ^= lb
                f = CTZ[lb]
                g = gs[C, f]
                a = fe[C, f]
                y = S
                while y:
                    l2 = y & (-y)
                    y ^= l2
                    e = CTZ[l2]
                    w = fe[S, e] * g + gs[S, e] * a
                    Wm[e, f] += w
                    Wm[f, e] += w


@njit(cache=True)
def subset_sums_w(n, Wm, U, Wsub, Msub):
    full = (1 << n) - 1
    Wsub[0] = 0
    Msub[0] = 0
    for O in range(1, full + 1):
        low = O & (-O)
        l = CTZ[low]
        R = O ^ low
        w = Wsub[R]
        x = R
        while x:
            l2 = x & (-x)
            x ^= l2
            w += Wm[l, CTZ[l2]]
        Wsub[O] = w
        Msub[O] = Msub[R] + U[l]


@njit(cache=True)
def record(n, code, H, Wsub, Msub, minv_p, minv_m, cnt_p, cnt_m, cv, cn, cc, cO, cM, stat):
    full = (1 << n) - 1
    for O in range(1, full + 1):
        w = Wsub[O]
        if w == H - 1:
            sg = 0
        elif w == H + 1:
            sg = 1
        else:
            continue
        M = Msub[O]
        sz = popcount(O)
        v = 0
        c = M
        while (c & 1) == 0:
            c >>= 1
            v += 1
        ci = c >> 1
        if sg == 0:
            cnt_p[sz] += 1
            if ci < CN and v < minv_p[ci]:
                minv_p[ci] = v
        else:
            cnt_m[sz] += 1
            if ci < CN and v < minv_m[ci]:
                minv_m[ci] = v
        if M > stat[sg]:
            stat[sg] = M
        if ci < NCERT:
            if v < cv[sg, ci] or (v == cv[sg, ci] and n < cn[sg, ci]):
                cv[sg, ci] = v
                cn[sg, ci] = n
                cc[sg, ci] = code
                cO[sg, ci] = O
                cM[sg, ci] = M


@njit(parallel=True, cache=True)
def census(n, reps, ext, part, nparts, minv_p, minv_m, cnt_p, cnt_m, cv, cn, cc, cO, cM, stat, ncand):
    """ext = 0: reps are codes of order n.  ext = 1: reps are codes of order n - 1; all min-score extensions of the
    representatives with index = part mod nparts."""
    nrep = reps.shape[0]
    for ch in prange(NCH):
        fe = np.zeros((1 << n, n), np.int64)
        gs = np.zeros((1 << n, n), np.int64)
        hS = np.zeros(1 << n, np.int64)
        Wm = np.zeros((n, n), np.int64)
        U = np.zeros(n, np.int64)
        Wsub = np.zeros(1 << n, np.int64)
        Msub = np.zeros(1 << n, np.int64)
        out = np.zeros(n, np.int64)
        base = np.zeros(n, np.int64)
        sc = np.zeros(n, np.int64)
        full = (1 << n) - 1
        for idx in range(ch, nrep, NCH):
            if ext == 0:
                decode(n, reps[idx], out)
                dp_range(n, out, fe, gs, hS, 1)
                wu_tables(n, fe, gs, hS, Wm, U)
                subset_sums_w(n, Wm, U, Wsub, Msub)
                record(n, reps[idx], hS[full], Wsub, Msub, minv_p[ch], minv_m[ch],
                       cnt_p[ch], cnt_m[ch], cv[ch], cn[ch], cc[ch], cO[ch], cM[ch], stat[ch])
                ncand[ch] += 1
            else:
                if idx % nparts != part:
                    continue
                m = n - 1
                decode(m, reps[idx], base)
                for u in range(m):
                    out[u] = base[u] | (1 << m)
                    sc[u] = popcount(base[u])
                out[m] = 0
                dp_range(n, out, fe, gs, hS, 1)          # the subsets without the new vertex are final
                for D in range(1 << m):
                    d = popcount(D)
                    ok = True
                    for u in range(m):
                        if d > sc[u] + (0 if (D >> u) & 1 else 1):
                            ok = False
                            break
                    if not ok:
                        continue
                    code = reps[idx]
                    for u in range(m):
                        out[u] = base[u]
                        if not ((D >> u) & 1):
                            out[u] |= 1 << m
                            code |= 1 << (m * (m - 1) // 2 + u)
                    out[m] = D
                    dp_range(n, out, fe, gs, hS, 1 << m)
                    wu_tables(n, fe, gs, hS, Wm, U)
                    subset_sums_w(n, Wm, U, Wsub, Msub)
                    record(n, code, hS[full], Wsub, Msub, minv_p[ch], minv_m[ch],
                           cnt_p[ch], cnt_m[ch], cv[ch], cn[ch], cc[ch], cO[ch], cM[ch], stat[ch])
                    ncand[ch] += 1


def run(n, reps, ext, part, nparts):
    minv_p = np.full((NCH, CN), 255, np.uint8)
    minv_m = np.full((NCH, CN), 255, np.uint8)
    cnt_p = np.zeros((NCH, 17), np.int64)
    cnt_m = np.zeros((NCH, 17), np.int64)
    cv = np.full((NCH, 2, NCERT), 255, np.uint8)
    cn = np.full((NCH, 2, NCERT), 255, np.uint8)
    cc = np.zeros((NCH, 2, NCERT), np.int64)
    cO = np.zeros((NCH, 2, NCERT), np.int64)
    cM = np.zeros((NCH, 2, NCERT), np.int64)
    stat = np.zeros((NCH, 2), np.int64)
    ncand = np.zeros(NCH, np.int64)
    census(n, reps, ext, part, nparts, minv_p, minv_m, cnt_p, cnt_m, cv, cn, cc, cO, cM, stat, ncand)
    minv = np.stack([minv_p.min(axis=0), minv_m.min(axis=0)])
    cert_v = np.full((2, NCERT), 255, np.uint8)
    cert_n = np.full((2, NCERT), 255, np.uint8)
    cert_c = np.zeros((2, NCERT), np.int64)
    cert_O = np.zeros((2, NCERT), np.int64)
    cert_M = np.zeros((2, NCERT), np.int64)
    for ch in range(NCH):
        for sg in range(2):
            better = (cv[ch, sg] < cert_v[sg]) | ((cv[ch, sg] == cert_v[sg]) & (cn[ch, sg] < cert_n[sg]))
            cert_v[sg][better] = cv[ch, sg][better]
            cert_n[sg][better] = cn[ch, sg][better]
            cert_c[sg][better] = cc[ch, sg][better]
            cert_O[sg][better] = cO[ch, sg][better]
            cert_M[sg][better] = cM[ch, sg][better]
    return dict(minv=minv, cnt_p=cnt_p.sum(axis=0), cnt_m=cnt_m.sum(axis=0), cert_v=cert_v, cert_n=cert_n,
                cert_c=cert_c, cert_O=cert_O, cert_M=cert_M, maxM=stat.max(axis=0), ncand=int(ncand.sum()))


def merge_certs(dst, src):
    for sg in range(2):
        better = (src["cert_v"][sg] < dst["cert_v"][sg]) | ((src["cert_v"][sg] == dst["cert_v"][sg]) & (src["cert_n"][sg] < dst["cert_n"][sg]))
        for key in ("cert_v", "cert_n", "cert_c", "cert_O", "cert_M"):
            dst[key][sg][better] = src[key][sg][better]


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "small"
    if mode == "small":
        res = {}
        for n in range(1, 10):
            t0 = time.time()
            reps = np.atleast_1d(np.loadtxt(f"cls_{n}.txt", dtype=np.int64))
            d = run(n, reps, 0, 0, 1)
            for k, v in d.items():
                res[f"{k}_{n}"] = v
            print(f"order {n}: {d['ncand']} classes, seeds by |O| {[(k, int(x)) for k, x in enumerate(d['cnt_p']) if x]}, "
                  f"anti-seeds {[(k, int(x)) for k, x in enumerate(d['cnt_m']) if x]}, max slopes {d['maxM']} [{time.time() - t0:.1f} s]", flush=True)
        np.savez_compressed("b_census_small.npz", **res)
        return
    if mode == "part":
        part, nparts = int(sys.argv[2]), int(sys.argv[3])
        t0 = time.time()
        reps = np.atleast_1d(np.loadtxt("cls_9.txt", dtype=np.int64))
        d = run(10, reps, 1, part, nparts)
        np.savez_compressed(f"b_census_part_{part}_{nparts}.npz", **d)
        print(f"order 10, part {part} of {nparts}: {d['ncand']} extensions, seeds by |O| {[(k, int(x)) for k, x in enumerate(d['cnt_p']) if x]}, "
              f"anti-seeds {[(k, int(x)) for k, x in enumerate(d['cnt_m']) if x]}, max slopes {d['maxM']} [{time.time() - t0:.1f} s]", flush=True)
        return
    # merge
    nparts = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    lines = []

    def P(s=""):
        print(s, flush=True)
        lines.append(s)

    small = np.load("b_census_small.npz")
    per_order = {}
    cum = np.full((2, CN), 255, np.uint8)
    certs = dict(cert_v=np.full((2, NCERT), 255, np.uint8), cert_n=np.full((2, NCERT), 255, np.uint8),
                 cert_c=np.zeros((2, NCERT), np.int64), cert_O=np.zeros((2, NCERT), np.int64), cert_M=np.zeros((2, NCERT), np.int64))
    P("S93 lane seeds: census of unit pairs (B, O), every vertex set O; r = +1 seeds, r = -1 anti-seeds")
    for n in range(1, 10):
        cum = np.minimum(cum, small[f"minv_{n}"])
        per_order[n] = cum.copy()
        merge_certs(certs, {k: small[f"{k}_{n}"] for k in certs})
        P(f"order {n}: {int(small[f'ncand_{n}'])} classes; seeds by |O|: {[(k, int(x)) for k, x in enumerate(small[f'cnt_p_{n}']) if x]}; "
          f"anti-seeds by |O|: {[(k, int(x)) for k, x in enumerate(small[f'cnt_m_{n}']) if x]}; largest slopes {small[f'maxM_{n}']}")
    nmax = 9
    if nparts and all(os.path.exists(f"b_census_part_{p}_{nparts}.npz") for p in range(nparts)):
        cp = np.zeros(17, np.int64)
        cm = np.zeros(17, np.int64)
        nc = 0
        mx = np.zeros(2, np.int64)
        for p in range(nparts):
            d = np.load(f"b_census_part_{p}_{nparts}.npz")
            cum = np.minimum(cum, d["minv"])
            cp += d["cnt_p"]
            cm += d["cnt_m"]
            nc += int(d["ncand"])
            mx = np.maximum(mx, d["maxM"])
            merge_certs(certs, d)
        per_order[10] = cum.copy()
        nmax = 10
        P(f"order 10: {nc} min-score extensions of the order-9 classes (every class of order 10 at least once, classes repeated); "
          f"seeds by |O|: {[(k, int(x)) for k, x in enumerate(cp) if x]}; anti-seeds by |O|: {[(k, int(x)) for k, x in enumerate(cm) if x]}; largest slopes {mx}")
    np.savez_compressed("b_census_minv.npz", **{f"minv_{n}": per_order[n] for n in per_order})

    names = ("seeds (r = +1)", "anti-seeds (r = -1)")
    ks = list(range(3, 21))
    for sg in range(2):
        P()
        P(f"== {names[sg]} ==")
        P("T1. number of odd c <= 2^k that are the odd part of a slope, B of order <= n   (row n; columns k; last row: number of odd c <= 2^k)")
        P("      k: " + " ".join(f"{k:7d}" for k in ks))
        for n in sorted(per_order):
            mv = per_order[n][sg]
            P(f"  n={n:2d}:  " + " ".join(f"{int((mv[:(1 << k) // 2] < 255).sum()):7d}" for k in ks))
        P("  all :  " + " ".join(f"{(1 << k) // 2:7d}" for k in ks))
        P(f"T2. number of odd c <= 2^k with a slope c 2^v, v <= j, B of order <= {nmax}   (row j)")
        mv = per_order[nmax][sg]
        for j in (0, 1, 2, 3, 4, 6, 8, 12, 254):
            P(f"  j={j:3d}:  " + " ".join(f"{int((mv[:(1 << k) // 2] <= j).sum()):7d}" for k in ks))
        P("T2b. number of odd c <= 2^k that are slopes themselves (v = 0), B of order <= n")
        for n in sorted(per_order):
            mv = per_order[n][sg]
            P(f"  n={n:2d}:  " + " ".join(f"{int((mv[:(1 << k) // 2] == 0).sum()):7d}" for k in ks))
        P("T3. number of integers M in [2^k, 2^(k+1)), k = 3..19, in the closure under doubling of the slopes with B of order <= n; second line: share")
        for n in sorted(per_order):
            mv = per_order[n][sg].astype(np.int64)
            row = []
            for k in range(3, 20):
                tot = 0
                for b in range(0, k + 1):          # c in [2^b, 2^(b+1)), v = k - b
                    lo, hi = (1 << b) >> 1, (1 << (b + 1)) >> 1
                    seg = mv[lo:hi]
                    tot += int((seg <= k - b).sum())
                row.append(tot)
            P(f"  n={n:2d}:  " + " ".join(f"{x:6d}" for x in row))
            P("   share " + " ".join(f"{x / (1 << k):6.3f}" for x, k in zip(row, range(3, 20))))
        mv = per_order[nmax][sg]
        first_missing = [2 * i + 1 for i in range(1 << 15) if mv[i] == 255][:40]
        P(f"least odd c without a {'seed' if sg == 0 else 'anti-seed'} among B of order <= {nmax}: {first_missing}")
        for n in sorted(per_order):
            mvn = per_order[n][sg]
            miss = [2 * i + 1 for i in range(2049) if mvn[i] == 255]
            miss2 = int((mvn[2049:8193] == 255).sum())
            P(f"  order <= {n}: {2049 - len(miss)} of the 2049 odd c <= 4097 covered (first missing {miss[:10]}); {6144 - miss2} of the 6144 odd c in (4097, 16385]")
        # slopes below 256: explicit list (complete once the order bound of PROOFS.txt applies)
        small_slopes = sorted(c * (1 << v) for c in range(1, 256, 2) for v in range(int(mv[c >> 1]), 9) if mv[c >> 1] < 255 and c * (1 << v) < 256)
        P(f"slopes below 256 (closure under doubling), B of order <= {nmax}: {len(small_slopes)}: {small_slopes}")

    # coverage of the odd N by the unit rays: N = 1 + c 2^s with s >= v_plus(c) + 1, or N = c 2^s - 1 with s >= v_minus(c) + 1
    P()
    P(f"T4. odd N in [2^k, 2^(k+1)) on a ray of offset +1 / of offset -1 / of either, pairs of order <= {nmax} (share of the 2^(k-1) odd N)")
    mvp = per_order[nmax][0].astype(np.int64)
    mvm = per_order[nmax][1].astype(np.int64)
    for k in range(3, 21):
        N = np.arange((1 << k) + 1, 1 << (k + 1), 2, dtype=np.int64)
        x = N - 1
        s = np.zeros(N.shape[0], np.int64)
        c = x.copy()
        while True:
            ev = (c & 1) == 0
            if not ev.any():
                break
            c[ev] >>= 1
            s[ev] += 1
        onp = s >= mvp[c >> 1] + 1
        y = N + 1
        s2 = np.zeros(N.shape[0], np.int64)
        c2 = y.copy()
        while True:
            ev = (c2 & 1) == 0
            if not ev.any():
                break
            c2[ev] >>= 1
            s2[ev] += 1
        onm = s2 >= mvm[c2 >> 1] + 1
        nodd = N.shape[0]
        P(f"  k = {k:2d}: +1: {int(onp.sum()):7d} ({onp.mean():.4f})   -1: {int(onm.sum()):7d} ({onm.mean():.4f})   either: {int((onp | onm).sum()):7d} ({(onp | onm).mean():.4f})   neither: {int((~(onp | onm)).sum())}")
    small_unc = [int(n) for n in range(3, 1 << 12, 2)
                 if not ((lambda x: (lambda s: s >= mvp[(x >> s) >> 1] + 1)((x & -x).bit_length() - 1))(n - 1)
                         or (lambda y: (lambda s: s >= mvm[(y >> s) >> 1] + 1)((y & -y).bit_length() - 1))(n + 1))]
    P(f"  odd N < 4096 on no unit ray of a pair of order <= {nmax}: {len(small_unc)}: {small_unc[:60]}")

    for sg, fn in ((0, "certs_plus.txt"), (1, "certs_minus.txt")):
        with open(fn, "w") as fh:
            fh.write("# c v order colcode Omask M     (pair (B, O), offset %+d, slope M = c 2^v; B = seedlib.out_from_colcode(order, colcode))\n" % (1 if sg == 0 else -1))
            for i in range(NCERT):
                if certs["cert_v"][sg, i] < 255:
                    fh.write(f"{2 * i + 1} {certs['cert_v'][sg, i]} {certs['cert_n'][sg, i]} {certs['cert_c'][sg, i]} {certs['cert_O'][sg, i]} {certs['cert_M'][sg, i]}\n")
    with open("b_census.out", "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
