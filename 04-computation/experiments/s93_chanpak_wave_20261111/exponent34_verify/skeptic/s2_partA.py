#!/usr/bin/env python3
"""Skeptic check 2: part (A) of THM-5141 (draft), exhaustively, with independent code.

 (A)1: no nonzero delta in {-2..2}^w with delta(2) = 0 and delta(1/2) = 0     (all 5^w strings, w <= 10)
 (A)2: eps -> (D, D*) injective on {-1,0,1}^w                                  (all 3^w strings, w <= 13)
       r_w(a, b) computed DIRECTLY from all 4^w pairs (u, u') and compared with 2^z(eps)   (w <= 11)
 (A)3: sum r^2 = 6^w                                                           (direct, w <= 11)
 (A)4: sum_{a != 0} r^2 |a|^(-1/2) versus 42 * 6^w 2^(-w/2); the constant      (direct w <= 11, via eps w <= 16)
       the intermediate claim F_j <= (5/3)(2/3)^j 6^w                          (via eps, w <= 13)
       the signed strings of value 0, +1, -1
"""
import sys
import numpy as np

out = []
def P(*a):
    s = " ".join(str(t) for t in a)
    print(s, flush=True)
    out.append(s)

def digits_table(w, lo, hi):
    """all strings in {lo..hi}^w as rows; returns (D, Dstar, zeros) arrays of int64."""
    base = hi - lo + 1
    n = base ** w
    idx = np.arange(n, dtype=np.int64)
    D = np.zeros(n, dtype=np.int64)
    Ds = np.zeros(n, dtype=np.int64)
    z = np.zeros(n, dtype=np.int64)
    t = idx.copy()
    for i in range(w):
        d = t % base + lo
        t //= base
        D += d << i if False else d * (1 << i)
        Ds += d * (1 << (w - 1 - i))
        z += (d == 0)
    return D, Ds, z

def revbits(arr, w):
    r = np.zeros_like(arr)
    for i in range(w):
        r |= ((arr >> i) & 1) << (w - 1 - i)
    return r

def main():
    # (A)1
    for w in range(1, 11):
        D, Ds, z = digits_table(w, -2, 2)
        both = np.count_nonzero((D == 0) & (Ds == 0))
        P("(A)1 w=%d: strings in {-2..2}^w with D = D* = 0: %d (must be 1, the zero string)" % (w, both))
    # (A)2 injectivity via eps
    for w in range(1, 14):
        D, Ds, z = digits_table(w, -1, 1)
        key = D * (1 << (w + 2)) + Ds   # |Ds| < 2^w; unique encoding
        nuniq = np.unique(key).size
        P("(A)2 w=%d: 3^w = %d strings, distinct (D,D*) = %d, injective: %s; D(eps)=0 only for eps=0: %s"
          % (w, 3 ** w, nuniq, nuniq == 3 ** w, np.count_nonzero(D == 0) == 1))
    # direct r_w from all pairs
    for w in range(1, 12):
        n = 1 << w
        u = np.arange(n, dtype=np.int64)
        ru = revbits(u, w)
        a = (u[:, None] - u[None, :]).ravel()
        b = (ru[:, None] - ru[None, :]).ravel()
        key = (a + n) * (2 * n) + (b + n)
        uk, cnt = np.unique(key, return_counts=True)
        energy = int(np.sum(cnt.astype(np.int64) ** 2))
        # compare with eps description
        D, Ds, z = digits_table(w, -1, 1)
        key2 = (D + n) * (2 * n) + (Ds + n)
        order = np.argsort(key2)
        same_support = np.array_equal(key2[order], uk)
        same_mult = np.array_equal((1 << z[order]), cnt) if same_support else False
        aa = uk // (2 * n) - n
        nz = aa != 0
        neg = float(np.sum(cnt[nz].astype(np.float64) ** 2 / np.sqrt(np.abs(aa[nz]).astype(np.float64))))
        bound = 42.0 * 6.0 ** w * 2.0 ** (-w / 2)
        P("direct w=%d: distinct differences %d (3^w=%d) support==image:%s r=2^z:%s  energy=%d 6^w=%d equal:%s  "
          "(A)4 sum=%.4f  constant=sum/(6^w 2^(-w/2))=%.4f (<=42: %s)"
          % (w, uk.size, 3 ** w, same_support, same_mult, energy, 6 ** w, energy == 6 ** w,
             neg, neg / (6.0 ** w * 2.0 ** (-w / 2)), neg <= bound))
    # (A)4 via eps for larger w, and F_j
    for w in range(2, 17):
        D, Ds, z = digits_table(w, -1, 1)
        wt = (4.0 ** z)
        nzm = D != 0
        neg = float(np.sum(wt[nzm] / np.sqrt(np.abs(D[nzm]).astype(np.float64))))
        const = neg / (6.0 ** w * 2.0 ** (-w / 2))
        line = "(A)4 eps w=%d: constant %.4f" % (w, const)
        if w <= 13:
            absD = np.abs(D)
            worst = 0.0
            for j in range(0, w + 1):
                Fj = float(np.sum(wt[absD < (1 << (w - j))]))
                ratio = Fj / ((5.0 / 3.0) * (2.0 / 3.0) ** j * 6.0 ** w)
                worst = max(worst, ratio)
            line += "  max_j F_j/((5/3)(2/3)^j 6^w) = %.4f (<=1: %s)" % (worst, worst <= 1.0)
        P(line)
    # signed strings of value 0, 1, -1
    for j in range(1, 13):
        D, Ds, z = digits_table(j, -1, 1)
        w0 = int(np.sum(4 ** z[D == 0]))
        w1 = int(np.sum(4 ** z[D == 1]))
        wm = int(np.sum(4 ** z[D == -1]))
        n0 = int(np.count_nonzero(D == 0)); n1 = int(np.count_nonzero(D == 1))
        P("value-0/+-1 strings j=%d: #value0=%d weight=%d (4^j=%d); #value1=%d weight=%d ((4^j-1)/3=%d); weight(-1)=%d"
          % (j, n0, w0, 4 ** j, n1, w1, (4 ** j - 1) // 3, wm))
    with open(sys.argv[1] if len(sys.argv) > 1 else "s2_partA.out", "w") as f:
        f.write("\n".join(out) + "\n")

if __name__ == "__main__":
    main()
