# fixed_state_A.py -- auditor proof-5134-AE.  Direct check of THM-5134 (A): fixed smaller state (h1, s).
# For a modulus m, d = gcd(g(h1), m), G = #{h2 > h1 : M <= X}: the classes c' of M modulo m.
#   exact part : classes with d not dividing c' are empty
#   approximate: every class with d | c' has (d/m) G + O((Y/y1)^(1 - theta_m))
# Printed: G, the largest |count - (d/m) G|, the same relative to (d/m) G, and the exponent
#   e = ln(max deviation)/ln(Y/y1)   (the theorem needs e <= 1 - theta_m + o(1); a trivial bound is e = 1).
import math, sys
import numpy as np
import numba as nb


@nb.njit(cache=True)
def build_psi(H):
    psi = np.zeros(H + 1, np.int32)
    for h in range(1, H + 1):
        k = h >> 1
        if (h & 1) == 0:
            psi[h] = psi[k] + k
        else:
            psi[h] = psi[k] + (k >> 1) + 1
    return psi


@nb.njit(cache=True)
def classes(psi, h1, s, X, m):
    x1 = np.int64(2 * psi[h1] + 1)
    y1 = np.int64(2 * h1 + 1)
    p2 = np.int64(1) << np.int64(s)
    cnt = np.zeros(m, np.int64)
    den = x1 + 2 * p2 * y1
    y2max = (2 * X - 1) // den
    h2max = (y2max - 1) // 2
    A = p2 * y1
    for h2 in range(h1 + 1, h2max + 1):
        M = x1 * np.int64(2 * psi[h2] + 1) + A * np.int64(2 * h2 + 1)
        if M <= X:
            cnt[M % m] += 1
    return cnt


def main():
    LMAX = 28
    Xmax = 1 << LMAX
    H = ((2 * (Xmax + 1) - 1) // 3) // 2 + 2
    psi = build_psi(H)
    g = lambda h: math.gcd(2 * int(psi[h]) + 1, 2 * h + 1)
    # first states of each type
    first = {}
    for mm in (3, 5, 7, 9, 15, 25, 27, 45, 49, 81, 125):
        for h in range(0, 400000):
            if g(h) % mm == 0:
                first[mm] = h
                break
    print("first states h with mm | g(h):", {k: (v, 2 * int(psi[v]) + 1, 2 * v + 1) for k, v in first.items()})
    cases = [(9, [0, 1, first[9], 3]), (27, [0, 1, first[9], first[27]]), (25, [0, first[5], first[25]]),
             (45, [0, 1, first[5], first[9], first[15], first[45]]), (49, [0, first[7], first[49]]), (81, [first[27], first[81]]), (125, [first[25], first[125]])]
    rates = {9: 0.866025, 27: 0.866025, 81: 0.866025, 25: 0.815263, 125: 0.815263, 45: 0.892191, 49: 0.932039}
    for m, hs in cases:
        th = math.log(1 / rates[m]) / math.log(2 / rates[m] ** 2)
        print()
        print("modulus %d   (1 - theta_m > %.4f)" % (m, 1 - th))
        for h1 in hs:
            d = math.gcd(g(h1), m)
            for s in (0, 1, 3):
                row = []
                for L in (20, 24, 28):
                    X = (1 << L) + (1 if s == 0 else 0)
                    cnt = classes(psi, h1, s, X, m)
                    G = int(cnt.sum())
                    if G == 0:
                        row.append("X=2^%d: G=0" % L)
                        continue
                    bad = sum(int(cnt[c]) for c in range(m) if c % d != 0)
                    assert bad == 0, (m, h1, s, L)
                    exp_ = d * G / m
                    dev = max(abs(int(cnt[c]) - exp_) for c in range(0, m, d))
                    Yy = (X / 2 ** s) / (2 * h1 + 1)
                    e = math.log(max(dev, 1e-9)) / math.log(Yy) if dev > 0 else float("nan")
                    row.append("X=2^%d: G=%d maxdev=%.1f rel=%.4f e=%.3f" % (L, G, dev, dev / exp_, e))
                print("  h1=%6d (x,y)=(%d,%d) d=%3d s=%d | " % (h1, 2 * int(psi[h1]) + 1, 2 * h1 + 1, d, s) + " | ".join(row))
    print()
    print("exact part (no class with d not dividing c is ever hit): verified in every case above (assert)")


if __name__ == "__main__":
    main()
