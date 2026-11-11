# Auditor recompute-5128, task 6: a-priori predictions for the window [A, B) by the rule of THM-5128 (G),
# recomputed from the text of the rule only:
#   W(h1, s) = (1/x1) int_{f0}^{f1} w(f) df,  f0 = x1/2 + 2^s y1,  f1 = x1 + 2^s y1,
#   w(f) = (B - A)/f for f < A,  B/f for A <= f < B  (0 for f >= B),
#   f_p = W-weighted share of p-divisible first states (per sheet),  D = 2 f_p - 1/p^2,
#   E = p D/(1 - D),   E_nd = p (f_p^2 - 1/p^2)/(1 - D).
# No window data are used here; the measured values are read from census_28.npz (a2_census.py) for comparison.
import sys
import math
import numpy as np
from numba import njit
from a0_defs import psi_table32


@njit(cache=True)
def W_one(x1, y1, s, A, B):
    b = y1 * (2.0 ** s)
    f0 = x1 / 2.0 + b
    f1 = x1 + b
    if f0 >= B:
        return 0.0
    tot = 0.0
    if f0 < A:
        hi = f1 if f1 < A else A
        tot += (B - A) * math.log(hi / f0)
    if f1 > A:
        lo = f0 if f0 > A else A
        hi = f1 if f1 < B else B
        if hi > lo:
            tot += B * math.log(hi / lo)
    return tot / x1


@njit(cache=True)
def predict(psi, mask_p, A, B, smax, tracked):
    """returns Wtot[sheet], Wdiv[sheet, j], Wtrk[t, sheet], number of (h1, s) used"""
    npr = mask_p.shape[0]
    Wtot = np.zeros(2, dtype=np.float64)
    Wdiv = np.zeros((2, npr), dtype=np.float64)
    Wtrk = np.zeros((tracked.shape[0], 2), dtype=np.float64)
    Wsmall = np.zeros(2, dtype=np.float64)  # part with y1^2 2^s < B (first state is "the smaller one", roughly)
    nterms = 0
    hsize = psi.shape[0]
    for s in range(0, smax + 1):
        sheet = 0 if s == 0 else 1
        two_s = 2.0 ** s
        for h1 in range(hsize):
            y1 = 2.0 * h1 + 1.0
            if y1 * (0.5 + two_s) >= B * 2.0:  # f0 > y1 (1/4 + 2^s) ... generous stop: f0 >= y1 two_s
                break
            x1 = 2.0 * psi[h1] + 1.0
            w = W_one(x1, y1, s, A, B)
            if w == 0.0:
                continue
            nterms += 1
            Wtot[sheet] += w
            xi = 2 * np.int64(psi[h1]) + 1
            yi = 2 * np.int64(h1) + 1
            for j in range(npr):
                p = mask_p[j]
                if yi % p == 0:
                    if xi % p == 0:
                        Wdiv[sheet, j] += w
            for t in range(tracked.shape[0]):
                if tracked[t] == h1:
                    Wtrk[t, sheet] += w
    return Wtot, Wdiv, Wtrk, nterms


def main():
    a = int(sys.argv[1]) if len(sys.argv) > 1 else 24
    b = int(sys.argv[2]) if len(sys.argv) > 2 else 27
    A = float(2 ** a)
    B = float(2 ** b)
    hsize = (2 ** b) // 2 + 4  # f0 >= 2^s y1 * (1/4 + 1) ... y1 < B/1.25 is enough; take all h1 < B/2
    psi = psi_table32(hsize)
    primes = np.array([3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 61, 73, 127, 251, 257, 263, 677, 691, 1039, 2729, 2731, 2741, 8191, 65521, 65537, 65539], dtype=np.int64)
    tracked = np.array([0, 1, 2, 3, 4, 8, 16, 32, 38, 64, 128, 218, 256, 1024], dtype=np.int64)
    Wtot, Wdiv, Wtrk, nterms = predict(psi, primes, A, B, b + 1, tracked)
    nodd = (2 ** b - 2 ** a) // 2
    print(f"== predictions for [2^{a}, 2^{b}) by the rule of (G); {nterms} pairs (h1, s) with W > 0 ==")
    print(f"predicted triples: sheet 0 {Wtot[0]:.1f}, sheet 1 {Wtot[1]:.1f}, total {Wtot.sum():.1f}; per odd N: {Wtot.sum() / nodd:.5f}   (claimed predicted 14.1417 on [2^24, 2^27))")
    Z = None
    try:
        Z = np.load("census_28.npz")
    except Exception:
        pass
    meas = {}
    if Z is not None and b <= 28:
        wb = list(range(a + 1, b + 1))
        all_p = [int(p) for p in Z["all_p"]]
        tot = Z["tot"]
        T = [int(tot[0, wb].sum()), int(tot[1, wb].sum())]
        print(f"measured triples:  sheet 0 {T[0]}, sheet 1 {T[1]}, total {sum(T)}; per odd N: {sum(T) / nodd:.5f}")
        print(f"predicted / measured: sheet 0 {Wtot[0] / T[0]:.6f}, sheet 1 {Wtot[1] / T[1]:.6f}, total {Wtot.sum() / sum(T):.6f}")
        for p in all_p:
            j = all_p.index(p)
            vals = []
            for sheet in (0, 1):
                d = int(Z["dist"][sheet, wb, j].sum())
                dsm = int(Z["Dsm"][sheet, wb, j].sum())
                dlg = int(Z["Dlg"][sheet, wb, j].sum())
                dbo = int(Z["Dboth"][sheet, wb, j].sum())
                ddg = int(Z["Ddiag"][sheet, wb, j].sum())
                E = d / ((T[sheet] - d) / (p - 1)) - 1
                anyd = dsm + dlg - dbo + ddg
                End = (d - anyd) / ((T[sheet] - anyd - (d - anyd)) / (p - 1)) - 1
                fmeas = (dsm / 2 + dlg / 2 + ddg) / T[sheet]
                vals.append((E, End, fmeas))
            meas[p] = vals
    claimed = {5: (0.3120, 0.3284), 11: (0.1705, 0.1749), 17: (0.2081, 0.2288), 19: (0.0511, 0.0497), 29: (0.0032, 0.0007),
               43: (0.0877, 0.0941), 257: (0.1815, 0.2044), 2731: (0.061, 0.069), 65537: (0.183, 0.206), 3: (0.515, 0.532),
               7: (0.174, 0.175), 13: (0.0459, 0.0448), 31: (-0.0061, -0.0089)}
    print()
    print("  p      f_p*p^2 (pred) s0 s1   f*p^2 (measured) s0 s1 |  E pred s0   s1    | E claimed-pred s0 s1 | E measured s0  s1   | E_nd pred s0 s1  | E_nd measured s0 s1")
    for j, p in enumerate(primes):
        p = int(p)
        row = f" {p:6d}"
        f = [Wdiv[0, j] / Wtot[0], Wdiv[1, j] / Wtot[1]]
        D = [2 * f[i] - 1 / p ** 2 for i in range(2)]
        E = [p * D[i] / (1 - D[i]) for i in range(2)]
        End = [p * (f[i] ** 2 - 1 / p ** 2) / (1 - D[i]) for i in range(2)]
        row += f"   {f[0] * p * p:8.4f} {f[1] * p * p:8.4f}"
        if p in meas:
            row += f"   {meas[p][0][2] * p * p:8.4f} {meas[p][1][2] * p * p:8.4f}"
        row += f" | {E[0]:+.4f} {E[1]:+.4f}"
        if p in claimed:
            row += f" |  {claimed[p][0]:+.4f} {claimed[p][1]:+.4f} "
        else:
            row += " |      -       -   "
        if p in meas:
            row += f" | {meas[p][0][0]:+.4f} {meas[p][1][0]:+.4f}"
        row += f" | {End[0]:+.4f} {End[1]:+.4f}"
        if p in meas:
            row += f" | {meas[p][0][1]:+.4f} {meas[p][1][1]:+.4f}"
        print(row)
    print()
    print("  tracked first states: predicted W(h1) = sum_s W(h1, s) against the measured number of ordered triples with first state h1")
    if Z is not None and b <= 28:
        trk = Z["trk"]
        ztr = [int(h) for h in Z["tracked"]]
        for t, h in enumerate(tracked):
            h = int(h)
            tz = ztr.index(h)
            m0 = int(trk[tz, 0, wb].sum())
            m1 = int(trk[tz, 1:, :][:, wb].sum())
            pr = Wtrk[t, 0] + Wtrk[t, 1]
            print(f"   h1 = {h:5d}: predicted {pr:14.1f}  measured {m0 + m1:10d}   relative difference {100 * (pr / (m0 + m1) - 1):+.5f} %   (sheet 0: {100 * (Wtrk[t, 0] / m0 - 1):+.5f} %, sheet 1: {100 * (Wtrk[t, 1] / m1 - 1):+.5f} %)")


if __name__ == "__main__":
    main()
