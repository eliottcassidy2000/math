#!/usr/bin/env python3
"""Audit proof-5128-AB: THM-5128 (B2) cells and (B4) regions, numerical tests with independent code.

(B2)  for m in (3, 5, 7, 9, 15) and n = NB: every cell (t top digits including the leading 1, j bottom digits,
      r = n - t - j free digits), every residue pair: max |count - 2^r/m^2| / (2 rate_m)^r, rate_m = max_k rho(T_k)/2
      (the infimum of the admissible rho_m).  A bounded value = an empirical C_m.
(B4)  (a) cut cells: for random integer a <= b, X and all 1 <= j <= n/2: number of cells that contain points of Omega
          and points outside (truly cut), and number of cells failing both sufficient conditions of the proof; claims:
          the proof's classification is sound; at most 4 cut top words per bottom word; at most 4 * 2^j cut cells.
      (b) statement (i): | #(B_n cap I cap Omega) - 2^n area | <= 6 * 2^(n - floor(n/2)) (the proof's constant).
      (c) statement (ii): max over residue vectors of |dev| / 2^(n (1 - theta_p)), theta_p from the limit rate.
      (d) statement (iii) with W = w_0(h) = (1/x) ln((x + y)/(x/2 + y)) (the weight of (D2), s = 0), W_0 = 2^-n.
      Also (a) and (b) are repeated with a > b (hypothesis b >= a dropped) to see whether it is needed.
"""
import math
import sys

import numpy as np

RATE = {3: 0.8660254038, 5: 0.7417603, 7: 0.9320393, 9: 0.8660254038, 15: 0.8921914, 11: 0.8587103}


def q_table(nb):
    size = 1 << nb
    q = np.zeros(size, dtype=np.int64)
    for k in range(nb):
        lo, hi = 1 << k, 1 << (k + 1)
        idx = np.arange(lo, hi, dtype=np.int64)
        half = idx >> 1
        q[lo:hi] = q[half] + (idx & 1) * ((half + 1) >> 1)
    return q


def bitrev(h, n):
    r = np.zeros_like(h)
    t = h.copy()
    for _ in range(n):
        r = (r << 1) | (t & 1)
        t >>= 1
    return r


def cells_B2(x, y, n, m, rate):
    lo, hi = 1 << (n - 1), 1 << n
    h = np.arange(lo, hi, dtype=np.int64)
    code = (x[lo:hi] % m) * m + (y[lo:hi] % m)
    worst = 0.0
    worst_at = None
    for t in range(1, n):
        for j in range(0, n - t):
            r = n - t - j
            if r < 1:
                continue
            top = (h >> (n - t)) - (1 << (t - 1))
            bottom = h & ((1 << j) - 1)
            cell = top * (1 << j) + bottom
            ncell = 1 << (t - 1 + j)
            cnt = np.bincount(cell * (m * m) + code, minlength=ncell * m * m).astype(np.float64)
            dev = float(np.max(np.abs(cnt - (1 << r) / (m * m))))
            ratio = dev / (2 * rate) ** r
            if ratio > worst:
                worst, worst_at = ratio, (t, j, r)
    return worst, worst_at


def area_region(n, a, b, X, i0, i1):
    """area{(u, v) in [1/2,1) x [0,1): i0 <= 2^n u < i1, 2^(n+1) u (a (1 - v/2) + b) <= X}."""
    ua = max(0.5, i0 / 2.0 ** n)
    ub = min(1.0, i1 / 2.0 ** n)
    if ub <= ua:
        return 0.0
    c2 = X / (a * 2.0 ** n)
    c1 = -1.0 - 2.0 * b / a
    u1 = c2 / (1.0 - c1)          # f = 1 for u <= u1
    u0 = c2 / (-c1)               # f = 0 for u >= u0
    area = 0.0
    # part with f = 1
    hi1 = min(ub, u1)
    if hi1 > ua:
        area += hi1 - ua
    # part with 0 < f < 1
    lo2, hi2 = max(ua, u1), min(ub, u0)
    if hi2 > lo2:
        area += c1 * (hi2 - lo2) + c2 * math.log(hi2 / lo2)
    return area


def main():
    nb = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    rng = np.random.default_rng(911)
    q = q_table(nb)
    hall = np.arange(1 << nb, dtype=np.int64)
    y = 2 * hall + 1
    x = y - 2 * q
    ok = True

    print("(B2) cells: max over all cells (t top, j bottom, r free digits) and residue pairs of |count - 2^r/m^2| / (2 rate)^r")
    for m in (3, 5, 7, 9, 15):
        for n in (12, 16, min(nb, 18)):
            w, at = cells_B2(x, y, n, m, RATE[m])
            print(f"   m = {m:2d}, n = {n:2d}: {w:8.4f}  at (t, j, r) = {at}")

    print("(B4)(a) cut cells, b >= a: per (n, j) the largest numbers over random (a, b, X)")
    print(f"{'n':>3} {'j':>3} {'truly cut':>10} {'proof-cut':>10} {'4*2^j':>7} {'max cut tops per bottom':>24} {'sound':>6}")
    for n in (6, 8, 10, 12, 14, 16):
        lo, hi = 1 << (n - 1), 1 << n
        h = hall[lo:hi]
        xs, ys = x[lo:hi], y[lo:hi]
        for j in range(1, n // 2 + 1):
            top = (h >> (n - j))
            bottom = h & ((1 << j) - 1)
            cell = (top - (1 << (j - 1))) * (1 << j) + bottom
            ncell = 1 << (2 * j - 1)
            percell = 1 << (n - 2 * j)
            brev = bitrev(np.arange(1 << j, dtype=np.int64), j)
            eps = n * 2.0 ** (-n)
            mx_true = mx_proof = mx_per = 0
            sound = True
            for trial in range(40):
                kind = trial % 4
                a = int(rng.integers(1, 1 << int(rng.integers(1, 20))))
                if kind == 0:
                    b = a
                elif kind == 1:
                    b = a + int(rng.integers(0, 4))
                else:
                    b = a * int(rng.integers(1, 1 << int(rng.integers(1, 10)))) + int(rng.integers(0, a))
                vals = a * xs + b * ys
                X = int(np.quantile(vals, rng.uniform(0.0, 1.0)))
                if trial % 7 == 0:
                    X = int(vals[int(rng.integers(0, len(vals)))])
                inside = vals <= X
                cin = np.bincount(cell[inside], minlength=ncell)
                truly = (cin > 0) & (cin < percell)
                # the proof's sufficient conditions
                A = np.arange(1 << (j - 1), 1 << j, dtype=np.float64)[:, None]
                Phi = (brev / 2.0 ** j)[None, :]
                R = 1.0 - Phi / 2.0
                scale = 2.0 ** (n - j + 1)
                suff_in = scale * (A + 1) * (a * (R + eps) + b) <= X
                suff_out = scale * A * (a * (R - 2.0 ** (-j - 1) - eps) + b) >= X
                proof_cut = ~(suff_in | suff_out)
                cin2 = cin.reshape(1 << (j - 1), 1 << j)
                # soundness: suff_in => all points inside, suff_out => none
                sound &= bool(np.all(cin2[suff_in] == percell) and np.all(cin2[suff_out] == 0))
                mx_true = max(mx_true, int(truly.sum()))
                mx_proof = max(mx_proof, int(proof_cut.sum()))
                mx_per = max(mx_per, int(proof_cut.sum(axis=0).max()))
            ok &= sound and mx_proof <= 4 * (1 << j) and mx_per <= 4
            if j in (1, 2, n // 4, n // 2) or n <= 8:
                print(f"{n:3d} {j:3d} {mx_true:10d} {mx_proof:10d} {4 * (1 << j):7d} {mx_per:24d} {str(sound):>6}")

    print("(B4)(b) statement (i): max over random regions of |count - 2^n area| / 2^(n - floor(n/2)) (proof: <= 6)")
    for label in ("b >= a", "a > b (hypothesis dropped)"):
        row = []
        for n in range(6, nb + 1, 2):
            lo, hi = 1 << (n - 1), 1 << n
            h = hall[lo:hi]
            xs, ys = x[lo:hi], y[lo:hi]
            worst = 0.0
            for trial in range(300 if n <= 14 else 80):
                a = int(rng.integers(1, 1 << int(rng.integers(1, 16))))
                if label == "b >= a":
                    b = a if trial % 5 == 0 else a * int(rng.integers(1, 1 << int(rng.integers(1, 8)))) + int(rng.integers(0, a))
                else:
                    b = int(rng.integers(0, a)) if trial % 3 else 0
                    if a == 1:
                        a, b = 5, int(rng.integers(0, 5))
                vals = a * xs + b * ys
                X = int(np.quantile(vals, rng.uniform(0.0, 1.0)))
                if trial % 2:
                    i0, i1 = sorted(int(v) for v in rng.integers(lo - 5, hi + 5, size=2))
                else:
                    i0, i1 = 0, 1 << 40
                cnt = int(np.sum((vals <= X) & (h >= i0) & (h < i1)))
                ar = area_region(n, a, b, X, i0, i1)
                worst = max(worst, abs(cnt - 2.0 ** n * ar) / 2.0 ** (n - n // 2))
            row.append((n, worst))
            if label == "b >= a":
                ok &= worst <= 6.0
        print(f"   {label:28s}: " + "  ".join(f"n={n}: {w:.3f}" for n, w in row))

    print("(B4)(c) statement (ii): max over random regions and residue vectors of |dev| / 2^(n (1 - theta_p))")
    for p in (3, 5, 7, 11):
        rate = RATE[p]
        L = math.log(1 / rate)
        theta = L / (math.log(2) + 2 * L)
        code = (x % p) * p + (y % p)
        row = []
        for n in range(10, nb + 1, 2):
            lo, hi = 1 << (n - 1), 1 << n
            h = hall[lo:hi]
            xs, ys, cs = x[lo:hi], y[lo:hi], code[lo:hi]
            worst = 0.0
            for trial in range(60):
                a = int(rng.integers(1, 1 << int(rng.integers(1, 16))))
                b = a if trial % 5 == 0 else a * int(rng.integers(1, 1 << int(rng.integers(1, 8)))) + int(rng.integers(0, a))
                vals = a * xs + b * ys
                X = int(np.quantile(vals, rng.uniform(0.05, 1.0)))
                if trial % 3 == 0:
                    i0, i1 = 0, 1 << 40
                else:
                    i0, i1 = sorted(int(v) for v in rng.integers(lo, hi, size=2))
                sel = (vals <= X) & (h >= i0) & (h < i1)
                cnt = np.bincount(cs[sel], minlength=p * p).astype(np.float64)
                dev = float(np.max(np.abs(cnt - cnt.sum() / (p * p))))
                worst = max(worst, dev / 2.0 ** (n * (1 - theta)))
            row.append((n, worst))
        print(f"   p = {p:2d} (theta = {theta:.4f}): " + "  ".join(f"n={n}: {w:.3f}" for n, w in row))

    print("(B4)(d) statement (iii) with W = w_0(h), W_0 = 2^-n: |sum over B_n of chi W| / (W_0 2^(n (1 - theta_p)))")
    for p in (3, 5, 7, 11):
        rate = RATE[p]
        L = math.log(1 / rate)
        theta = L / (math.log(2) + 2 * L)
        div = ((x % p) == 0) & ((y % p) == 0)
        chi = div.astype(np.float64) - 1.0 / (p * p)
        xf, yf = x.astype(np.float64), y.astype(np.float64)
        W = np.log((xf + yf) / (xf / 2 + yf)) / xf
        row = []
        for n in range(8, nb + 1, 2):
            lo, hi = 1 << (n - 1), 1 << n
            s = float(np.sum(chi[lo:hi] * W[lo:hi]))
            row.append((n, abs(s) / (2.0 ** (-n) * 2.0 ** (n * (1 - theta)))))
        print(f"   p = {p:2d}: " + "  ".join(f"n={n}: {w:.3f}" for n, w in row))
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
