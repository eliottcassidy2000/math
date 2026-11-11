#!/usr/bin/env python3
"""S92 (THM-5134): the by-value class law of the join family F1 modulo an odd modulus m, prime powers included.

F1: N = M - [s = 0], M = x1 x2 + 2^s y1 y2, over ordered triples (h1, h2, s); x(h) = 2 psi(h) + 1 = y - 2 q(h),
y(h) = 2h + 1, q(2h) = q(h), q(2h + 1) = q(h) + ceil(h/2).  Sheet 0: s = 0 (c_0 = -1); sheet 1: s >= 1 (c_1 = 0).
g(h) = gcd(x(h), y(h)).  For a divisor d of m, a triple with h1 != h2 is of type d if its state of SMALLER word has
gcd(g, m) = d.

THM-5134 (B):   #{N = c mod m} = sum over d | gcd(m, c - c_sigma) of (d/m) #{triples of type d} + O(X^(1 - theta_m/2)).
THM-5134 (D):   the class factor m * share tends to L_m(e) = sum_{d | gcd(m, e)} (1/d) prod_{p | m/d} (1 - p^-2),
                e = c - c_sigma; for m = p^k and v = v_p(e) < k this is (1 - p^-2) sigma(p^v)/p^v.

Modes
  python3 odd_modulus_census.py cum K              cumulative census N <= 2^K: the exact remainder of (B)
  python3 odd_modulus_census.py win KHI KLO        window [2^KLO, 2^KHI): measured factors by gcd pattern, the
                                                   a-priori prediction (no enumeration of the window is used in it),
                                                   and the limits
  python3 odd_modulus_census.py pred KHI KLO FILE  write the a-priori predictions only (JSON)
  python3 odd_modulus_census.py eval KHI KLO FILE  enumerate the window and evaluate the registered tests of FILE

A-priori prediction.  For a smaller state (h1, s) the number of h2 > h1 with value in [A, B) is modelled by
    W(h1, s) = int_{1/2}^{1} 2 d rho * (1/2) max(0, B/f - max(A/f, y1)),     f = x1 rho + 2^s y1,
(the second state has y2 odd, y2 > y1, and x2/y2 uniform on [1/2, 1]: THM-5128 (B3)); the predicted share of the
triples of type d is sum over the smaller states of type d of W, divided by the sum over all smaller states, and the
predicted class factor of a class with gcd(m, c - c_sigma) = g0 is  Phi(g0) = sum_{d | g0} d F_d.
"""
import json
import math
import sys
import time

import numpy as np
from numba import njit, prange

MODULI = [9, 27, 81, 25, 125, 49, 45, 75, 63, 225, 15, 21]
NT = 8


def divisors(m):
    return [d for d in range(1, m + 1) if m % d == 0]


@njit(cache=True)
def build_q(size):
    q = np.zeros(size, np.int32)
    for h in range(1, size):
        half = h >> 1
        q[h] = q[half] + (((half + 1) >> 1) if (h & 1) else 0)
    return q


@njit(cache=True)
def gcd(a, b):
    while b:
        a, b = b, a % b
    return a


@njit(parallel=True, cache=True)
def census(q, A, B, K, moduli, dindex, cnt, diag):
    """cnt[t, sheet, j, di, c]: ordered triples with h1 != h2, A <= N < B, type index di for modulus j, class c.
    diag[t, sheet, j, c]: the triples with h1 = h2.  dindex[j, d] = index of the divisor d of moduli[j]."""
    nm = moduli.shape[0]
    H = q.shape[0]
    for t in prange(NT):
        for s in range(K + 1):
            sheet = 0 if s == 0 else 1
            d0 = 1 if s == 0 else 0
            h1 = np.int64(t)
            while h1 < H:
                y1 = 2 * h1 + 1
                Y = y1 << s
                x1 = y1 - 2 * np.int64(q[h1])
                # every value with h2 >= h1 exceeds (x1/2 + Y) y1 - d0; x1 > y1/2 is not monotone in h1, so the loop
                # stops only when (y1/4 + Y/2) y1 > B + d0, which then holds for every larger h1
                if (x1 + 2 * Y) * y1 > 2 * (B + d0):
                    if (y1 + 2 * Y) * y1 > 4 * (B + d0):
                        break
                    h1 += NT
                    continue
                # h2 range: N < B needs (x1/2 + Y) y2 < B + d0; N >= A needs (x1 + Y) y2 >= A + d0
                y2hi = (2 * (B + d0)) // (x1 + 2 * Y) + 1
                y2lo = (A + d0) // (x1 + Y)
                h2lo = (y2lo - 1) // 2
                if h2lo < h1:
                    h2lo = h1
                h2hi = y2hi // 2 + 1
                if h2hi >= H:
                    h2hi = H - 1
                gs = gcd(x1, y1)
                for h2 in range(h2lo, h2hi + 1):
                    y2 = 2 * h2 + 1
                    x2 = y2 - 2 * np.int64(q[h2])
                    N = x1 * x2 + Y * y2 - d0
                    if N < A or N >= B:
                        continue
                    if h2 == h1:
                        for j in range(nm):
                            diag[t, sheet, j, N % moduli[j]] += 1
                    else:
                        for j in range(nm):
                            m = moduli[j]
                            cnt[t, sheet, j, dindex[j, gcd(gs, m)], N % m] += 2
                h1 += NT


def weights(q, A, B, K, moduli, divs):
    """a-priori shares F[sheet][j][di] of the types, from the model W (Gauss-Legendre in rho)."""
    xs, ws = np.polynomial.legendre.leggauss(24)
    rho = 0.75 + 0.25 * xs
    wr = 0.25 * ws * 2.0            # density 2 d rho on [1/2, 1]
    F = [[np.zeros(len(divs[j])) for j in range(len(moduli))] for _ in (0, 1)]
    for s in range(K + 1):
        sheet = 0 if s == 0 else 1
        d0 = 1 if s == 0 else 0
        Ap, Bp = A + d0, B + d0
        h1 = 0
        while True:
            y1 = 2 * h1 + 1
            Y = y1 << s
            if (y1 / 4 + Y) * y1 > Bp:        # x1 > y1/2, so even rho = 1/2 and y2 = y1 give a value above the window
                break
            x1 = y1 - 2 * int(q[h1])
            f = x1 * rho + Y
            w = float(np.sum(wr * 0.5 * np.maximum(0.0, Bp / f - np.maximum(Ap / f, y1))))
            if w > 0:
                g = math.gcd(x1, y1)
                for j, m in enumerate(moduli):
                    F[sheet][j][divs[j].index(math.gcd(g, m))] += w
            h1 += 1
    out = []
    for sheet in (0, 1):
        out.append([(F[sheet][j] / F[sheet][j].sum()).tolist() for j in range(len(moduli))])
    return out


def limit_factor(m, g0):
    tot = 0.0
    for d in divisors(g0):
        v = 1.0 / d
        r = m // d
        p = 3
        while p <= r:
            if r % p == 0:
                v *= 1 - 1.0 / (p * p)
                while r % p == 0:
                    r //= p
            p += 2
        tot += v
    return tot


def run_census(A, B, K, moduli):
    divs = [divisors(m) for m in moduli]
    q = build_q((B >> 1) + 4)
    dindex = np.zeros((len(moduli), max(moduli) + 1), np.int64)
    for j, m in enumerate(moduli):
        for i, d in enumerate(divs[j]):
            dindex[j, d] = i
    nd = max(len(d) for d in divs)
    cnt = np.zeros((NT, 2, len(moduli), nd, max(moduli)), np.int64)
    diag = np.zeros((NT, 2, len(moduli), max(moduli)), np.int64)
    census(q, A, B, K, np.array(moduli, dtype=np.int64), dindex, cnt, diag)
    return q, divs, cnt.sum(axis=0), diag.sum(axis=0)


def table(moduli, divs, cnt, diag):
    """per sheet and modulus: total, type shares, class factors, exact remainder of (B)."""
    out = {}
    for sheet in (0, 1):
        cs = -1 if sheet == 0 else 0
        for j, m in enumerate(moduli):
            c = cnt[sheet, j, :len(divs[j]), :m].astype(np.float64)
            dg = diag[sheet, j, :m].astype(np.float64)
            tot = c.sum() + dg.sum()
            types = c.sum(axis=1)
            count = c.sum(axis=0) + dg
            pred = np.zeros(m)
            g0s = np.zeros(m, np.int64)
            for cl in range(m):
                g0 = math.gcd(m, (cl - cs) % m)
                g0s[cl] = g0
                for i, d in enumerate(divs[j]):
                    if g0 % d == 0:
                        pred[cl] += d / m * types[i]
            rem = m * np.abs(count - pred).max() / tot
            exact = True                  # a triple of type d has N - c_sigma divisible by d (exact part of (A))
            for i, d in enumerate(divs[j]):
                for cl in range(m):
                    if (cl - cs) % d != 0 and c[i, cl] != 0:
                        exact = False
            groups = {}
            for g0 in divs[j]:
                sel = g0s == g0
                fac = m * count[sel] / tot
                phi = sum(d * types[i] / tot for i, d in enumerate(divs[j]) if g0 % d == 0)
                groups[g0] = (float(fac.mean()), float(fac.min()), float(fac.max()), float(phi), limit_factor(m, g0),
                              int(sel.sum()))
            out[(sheet, m)] = dict(total=tot, types=(types / tot).tolist(), rem=rem, groups=groups,
                                   ndiag=float(dg.sum()), exact=exact)
    return out


def main():
    mode = sys.argv[1]
    t0 = time.time()
    if mode == "cum":
        K = int(sys.argv[2])
        A, B = 1, (1 << K) + 1
        q, divs, cnt, diag = run_census(A, B, K, MODULI)
        tab = table(MODULI, divs, cnt, diag)
        ok = True
        print(f"cumulative census, odd N <= 2^{K}: triples sheet 0: {tab[(0, 9)]['total']:.0f}, sheet 1: "
              f"{tab[(1, 9)]['total']:.0f}  ({time.time() - t0:.0f} s)")
        print("exact remainder of THM-5134 (B):  m max_c |#{N = c} - sum_{d | gcd(m, c - c_s)} (d/m) #type d| / #T"
              "   (noise = sqrt(m/#T))")
        print(f"  {'m':>4} {'sheet 0':>9} {'sheet 1':>9} {'noise':>8}")
        for m in MODULI:
            r0, r1 = tab[(0, m)]['rem'], tab[(1, m)]['rem']
            print(f"  {m:4d} {r0:9.5f} {r1:9.5f} {math.sqrt(m / tab[(1, m)]['total']):8.5f}")
            ok &= r0 < 0.25 and r1 < 0.25 and tab[(0, m)]['exact'] and tab[(1, m)]['exact']
        print("every triple of type d lies in a class with d | c - c_s (exact), and every remainder is below 0.25:", ok)
        print("class factors m * share by g0 = gcd(m, c - c_s): measured mean [min, max] | first-order sum_{d|g0} d F_d"
              " | limit L_m(g0)")
        for m in MODULI:
            for sheet in (0, 1):
                cells = []
                for g0, (mean, lo, hi, phi, lim, n) in sorted(tab[(sheet, m)]['groups'].items()):
                    cells.append(f"g0={g0}: {mean:.4f} [{lo:.4f}, {hi:.4f}] | {phi:.4f} | {lim:.4f}")
                print(f"  m = {m:3d} sheet {sheet}: " + ";  ".join(cells))
        print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
        sys.exit(0 if ok else 1)

    KHI, KLO = int(sys.argv[2]), int(sys.argv[3])
    A, B = 1 << KLO, 1 << KHI
    moduli = [9, 27, 81, 25, 125, 49, 45, 75]
    divs = [divisors(m) for m in moduli]
    if mode in ("pred", "win"):
        q = build_q(1 << 16)
        F = weights(q, A, B, KHI, moduli, divs)
        pred = {}
        for sheet in (0, 1):
            for j, m in enumerate(moduli):
                for g0 in divs[j]:
                    pred[f"{sheet},{m},{g0}"] = sum(d * F[sheet][j][i] for i, d in enumerate(divs[j]) if g0 % d == 0)
        if mode == "pred":
            json.dump(dict(KHI=KHI, KLO=KLO, moduli=moduli, phi=pred, F=F), open(sys.argv[4], "w"), indent=1)
            print(f"a-priori predictions for odd N in [2^{KLO}, 2^{KHI}) written to {sys.argv[4]}  ({time.time() - t0:.0f} s)")
            for m in (9, 27, 25, 49):
                j = moduli.index(m)
                print(f"  m = {m}: predicted class factors by g0 (sheet 0 / sheet 1): "
                      + ", ".join(f"{g0}: {pred[f'0,{m},{g0}']:.4f} / {pred[f'1,{m},{g0}']:.4f}" for g0 in divs[j]))
            return
    else:
        reg = json.load(open(sys.argv[4]))
        assert reg["KHI"] == KHI and reg["KLO"] == KLO and reg["moduli"] == moduli
        pred = reg["phi"]
    q, divs, cnt, diag = run_census(A, B, KHI, moduli)
    tab = table(moduli, divs, cnt, diag)
    print(f"window odd N in [2^{KLO}, 2^{KHI}): triples sheet 0: {tab[(0, 9)]['total']:.0f}, sheet 1: "
          f"{tab[(1, 9)]['total']:.0f}  ({time.time() - t0:.0f} s)")
    print("g0 = gcd(m, c - c_s): measured class factor (mean over the classes, [min, max]) | a-priori prediction | "
          "first-order value with the measured types | limit")
    dev = {}
    flat = {}
    order_ok = True
    closer = []
    for m in moduli:
        for sheet in (0, 1):
            cells = []
            means = []
            for g0, (mean, lo, hi, phi, lim, n) in sorted(tab[(sheet, m)]['groups'].items()):
                pr = pred[f"{sheet},{m},{g0}"]
                cells.append(f"g0={g0}: {mean:.4f} [{lo:.4f}, {hi:.4f}] | {pr:.4f} | {phi:.4f} | {lim:.4f}")
                dev[(m, sheet)] = max(dev.get((m, sheet), 0.0), abs(mean - pr))
                flat[(m, sheet)] = max(flat.get((m, sheet), 0.0), hi - mean, mean - lo)
                means.append((g0, mean))
                if g0 > 1 and m in (9, 27, 25):
                    closer.append(abs(mean - pr) < abs(mean - lim))
            # ordering along chains of divisors: g0 | g0' implies mean(g0) < mean(g0')
            for g0, a in means:
                for g1, b in means:
                    if g1 != g0 and g1 % g0 == 0 and not a < b:
                        order_ok = False
            print(f"  m = {m:3d} sheet {sheet}: " + ";  ".join(cells))
    print("largest |measured group mean - a-priori prediction| (sheet 0 / sheet 1) and largest deviation of a class "
          "from its group mean:")
    for m in moduli:
        print(f"  m = {m:3d}: {dev[(m, 0)]:.4f} / {dev[(m, 1)]:.4f}   flatness {flat[(m, 0)]:.4f} / {flat[(m, 1)]:.4f}")
    print(f"class factors increase along divisor chains g0 | g0': {order_ok}")
    print(f"for g0 > 1 at m = 9, 27, 25 the measured factor is closer to the prediction than to the limit: "
          f"{sum(closer)} of {len(closer)}")
    if mode == "eval":
        tol = reg["tol"]
        fl = reg["flat"]
        t1 = all(dev[(m, sh)] <= tol[str(m)] for m in moduli for sh in (0, 1))
        t2 = all(flat[(m, sh)] <= fl[str(m)] for m in moduli for sh in (0, 1))
        t4 = sum(closer) == len(closer)
        print(f"P1 (|group mean - prediction| <= registered tolerance, 16 cases): {'CONFIRMED' if t1 else 'FAILED'}")
        print(f"P2 (flatness inside a group <= registered bound, 16 cases): {'CONFIRMED' if t2 else 'FAILED'}")
        print(f"P3 (closer to the first-order prediction than to the limit, g0 > 1, m = 9, 27, 25): "
              f"{'CONFIRMED' if t4 else 'FAILED'}")
        print("(the ordering along divisor chains is NOT a registered test: it fails in sample)")
        print("PRE-REGISTERED WINDOW TESTS P1-P3 CONFIRMED" if (t1 and t2 and t4) else "SOME PRE-REGISTERED TEST FAILED")


if __name__ == "__main__":
    main()
