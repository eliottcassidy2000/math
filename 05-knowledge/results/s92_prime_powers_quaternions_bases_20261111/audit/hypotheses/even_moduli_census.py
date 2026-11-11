#!/usr/bin/env python3
# Auditor "hypotheses": HYP-9665 (A). Exact census of F1 by value, joint histogram of
#   (sheet, type of the smaller state w.r.t. 315 = 9*5*7, N mod 2520 = 8*9*5*7).
# From it: for a modulus Mod = 2^j m (m | 315 odd) the remainder of the product law
#   #{N = c mod 2^j m} ?= [c odd] 2^(1-j) sum_{d | gcd(m, c - c_sigma)} (d/m) D^(d)
# and the flatness of N mod 2^j inside every class mod m.
# usage: python3 even_moduli_census.py E      (X = 2^E)
import sys, time, math
import numpy as np
from numba import njit, prange, get_num_threads, get_thread_id
from f1core import build_psi, gcd

LCM = 2520
MODD = 315


@njit(parallel=True, cache=True)
def census(X, xs, ys, tyidx, nthreads, chunk):
    # hist[thread, sheet, typeindex, N mod 2520]
    nH = xs.shape[0]
    hist = np.zeros((nthreads, 2, 12, LCM), dtype=np.int64)
    nchunks = (nH + chunk - 1) // chunk
    for ci in prange(nchunks):
        tid = get_thread_id()
        lo = ci * chunk
        hi = min(nH, lo + chunk)
        for h1 in range(lo, hi):
            x1 = xs[h1]
            y1 = ys[h1]
            s = 0
            while True:
                p2 = np.int64(1) << s
                sub = 1 if s == 0 else 0
                sheet = 0 if s == 0 else 1
                # N >= x1 (h2+1) + p2 y1 (2 h2 + 1) - sub
                num = X + sub - x1 - p2 * y1
                if num < 0:
                    break
                h2max = num // (x1 + 2 * p2 * y1)
                if h2max >= nH:
                    h2max = nH - 1
                a = x1
                b = p2 * y1
                for h2 in range(h2max + 1):
                    N = a * xs[h2] + b * ys[h2] - sub
                    if N <= X:
                        hs = h1 if h1 <= h2 else h2
                        hist[tid, sheet, tyidx[hs], N % LCM] += 1
                s += 1
    return hist


def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]


def L_m(m, e):
    # L_m(e) = sum_{d | gcd(m,e)} d^-1 prod_{p | m/d} (1 - p^-2)
    g = math.gcd(m, e % m) if e % m else m
    tot = 0.0
    for d in divisors(g):
        md = m // d
        pr = 1.0
        q = md
        p = 3
        ps = set()
        while q > 1:
            if q % p == 0:
                ps.add(p)
                while q % p == 0:
                    q //= p
            p += 2
        for p in ps:
            pr *= 1 - p ** -2
        tot += pr / d
    return tot


def main():
    E = int(sys.argv[1])
    X = 1 << E
    t0 = time.time()
    nH = (X - 1) // 3 + 2
    psi = build_psi(nH)
    h = np.arange(nH, dtype=np.int64)
    xs = 2 * psi + 1
    ys = 2 * h + 1
    g = np.gcd(xs, ys)
    g315 = np.gcd(g, MODD)
    divs = divisors(MODD)  # 12 divisors
    didx = {d: i for i, d in enumerate(divs)}
    lut = np.zeros(MODD + 1, dtype=np.int64)
    for d in divs:
        lut[d] = didx[d]
    tyidx = lut[g315]
    nthreads = get_num_threads()
    hist = census(X, xs, ys, tyidx, nthreads, 64).sum(axis=0)
    print("# even_moduli_census.py  X = 2^%d  (auditor hypotheses; own code)" % E)
    print("# triples: sheet0 %d  sheet1 %d   (%.1f s)" % (hist[0].sum(), hist[1].sum(), time.time() - t0))
    csig = {0: -1, 1: 0}
    # exactness check: every triple of type d has d | N - c_sigma
    bad = 0
    for sh in (0, 1):
        for d in divs:
            for c in range(LCM):
                if hist[sh, didx[d], c] and (c - csig[sh]) % d != 0:
                    bad += hist[sh, didx[d], c]
    print("# triples of type d with d not dividing N - c_sigma (must be 0): %d" % bad)
    import json
    json.dump({"E": E, "shape": list(hist.shape), "axes": "sheet, type index (divisors of 315 in increasing order), N mod 2520", "hist": hist.tolist()}, open("even_moduli_hist_%d.json" % E, "w"))

    print()
    print("ODD moduli (THM-5134 (B) remainder, for reference):  m * max_c |#{N=c} - sum_d (d/m) D^(d)| / #T")
    for m in (3, 9, 5, 7, 15, 45, 21, 35, 63, 105, 315):
        out = []
        for sh in (0, 1):
            T = hist[sh].sum()
            cnt = np.zeros(m)
            D = {}
            for d in divisors(m):
                D[d] = 0
            for di, dd in enumerate(divs):
                dm = math.gcd(dd, m)
                D[dm] += hist[sh, di].sum()
                for c in range(LCM):
                    cnt[c % m] += hist[sh, di, c]
            worst = 0.0
            for c in range(m):
                pred = sum((d / m) * D[d] for d in divisors(math.gcd(m, (c - csig[sh]) % m) if (c - csig[sh]) % m else m))
                worst = max(worst, abs(cnt[c] - pred))
            out.append(m * worst / T)
        print("  m = %3d   %.4f / %.4f" % (m, out[0], out[1]))

    print()
    print("EVEN moduli 2^j m. Columns (sheet 0 / sheet 1):")
    print("  R_prod  = (2^j m) max_{c odd} |#{N=c} - 2^(1-j) sum_{d | gcd(m,c-c_s)} (d/m) D^(d)| / #T   (first-order product law)")
    print("  R_m     = the same remainder at the odd modulus m alone")
    print("  R_flat  = max over classes e mod m and odd a mod 2^j of |2^(j-1) #{N=a mod 2^j, N=e mod m}/#{N=e mod m} - 1|")
    print("  R_lim   = max_{c odd} |(2^j m) share(c) - 2 L_m(c - c_s)|   (distance to the LIMIT law)")
    print("  noise   = sqrt(2^(j-1) m / #T) (Poisson scale of one class share times the number of odd classes)")
    for (j, m) in ((2, 3), (3, 3), (2, 9), (3, 9), (2, 5), (3, 5), (2, 7), (3, 7), (2, 15), (3, 15), (2, 45), (3, 45), (3, 315)):
        Mod = (1 << j) * m
        res = []
        for sh in (0, 1):
            T = hist[sh].sum()
            cnt = np.zeros(Mod)
            cntm = np.zeros(m)
            D = {d: 0 for d in divisors(m)}
            for di, dd in enumerate(divs):
                dm = math.gcd(dd, m)
                D[dm] += hist[sh, di].sum()
                for c in range(LCM):
                    cnt[c % Mod] += hist[sh, di, c]
                    cntm[c % m] += hist[sh, di, c]
            Rprod = 0.0
            Rm = 0.0
            Rflat = 0.0
            Rlim = 0.0
            for c in range(m):
                e = (c - csig[sh]) % m
                gg = math.gcd(m, e) if e else m
                pred = sum((d / m) * D[d] for d in divisors(gg))
                Rm = max(Rm, m * abs(cntm[c] - pred) / T)
            for c in range(Mod):
                if c % 2 == 0:
                    assert cnt[c] == 0
                    continue
                e = (c - csig[sh]) % m
                gg = math.gcd(m, e) if e else m
                pred = 2.0 ** (1 - j) * sum((d / m) * D[d] for d in divisors(gg))
                Rprod = max(Rprod, Mod * abs(cnt[c] - pred) / T)
                Rflat = max(Rflat, abs(2 ** (j - 1) * cnt[c] / cntm[c % m] - 1))
                Rlim = max(Rlim, abs(Mod * cnt[c] / T - 2 * L_m(m, e)))
            res.append((Rprod, Rm, Rflat, Rlim, math.sqrt(2 ** (j - 1) * m / T)))
        print("  2^%d * %3d = %4d   R_prod %.4f / %.4f   R_m %.4f / %.4f   R_flat %.4f / %.4f   R_lim %.4f / %.4f   noise %.4f / %.4f" % (
            j, m, Mod, res[0][0], res[1][0], res[0][1], res[1][1], res[0][2], res[1][2], res[0][3], res[1][3], res[0][4], res[1][4]))

    print()
    print("Class factors (2^j m) * share at 12 (odd classes), measured against 2 L_3(c - c_s):")
    for sh in (0, 1):
        T = hist[sh].sum()
        cnt = np.zeros(12)
        for c in range(LCM):
            cnt[c % 12] += hist[sh, :, c].sum()
        print("  sheet %d: " % sh + "  ".join("c=%d: %.4f (lim %.4f)" % (c, 12 * cnt[c] / T, 2 * L_m(3, c - csig[sh])) for c in range(1, 12, 2)))
    print("# total time %.1f s" % (time.time() - t0))


if __name__ == "__main__":
    main()
