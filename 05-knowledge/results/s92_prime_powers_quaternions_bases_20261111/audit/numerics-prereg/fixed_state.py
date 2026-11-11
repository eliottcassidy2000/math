# Auditor numerics-prereg: direct look at THM-5134 (A) (fixed smaller state): for fixed (h1, s) and d = gcd(g(h1), m),
# the value M = x1 x2 + 2^s y1 y2 over h2 > h1 with M <= X should be uniform on the multiples of d modulo m.
# Prints, for a few small states, G = #{h2 > h1 : M <= X}, the number of h2 with d not dividing M (must be 0), and
#   dev = max over c' in d Z/m of | (m/d) #{M = c'} / G - 1 |    (relative deviation from the uniform share d/m).
import sys
import math
import numpy as np
import numba as nb
from aud_common import psi_table, gcd_nb

MODS = np.array([9, 27, 81, 25, 125, 49], dtype=np.int64)


@nb.njit
def fixed(psi, h1, s, X, mods):
    x1 = 2 * psi[h1] + 1
    y1 = 2 * h1 + 1
    b = y1 << s
    nm = mods.shape[0]
    cnt = np.zeros((nm, 128), np.int64)
    G = 0
    h2 = h1 + 1
    while True:
        y2 = 2 * h2 + 1
        if x1 * ((y2 + 1) // 2) + b * y2 > X:
            break
        M = x1 * (2 * psi[h2] + 1) + b * y2
        if M <= X:
            G += 1
            for i in range(nm):
                cnt[i, M % mods[i]] += 1
        h2 += 1
    return G, cnt


def main():
    Ks = (22, 24, 26, 28)
    psi = psi_table((1 << 28) // 3 + 8)
    states = [(0, 0), (0, 1), (1, 0), (1, 1), (2, 1), (4, 1), (3, 1), (38, 1), (256, 1), (512, 1)]
    print("fixed smaller state (h1, s): relative deviation from uniformity on the multiples of d = gcd(g(h1), m)")
    for (h1, s) in states:
        x1 = 2 * int(psi[h1]) + 1
        y1 = 2 * h1 + 1
        g = math.gcd(x1, y1)
        print("  h1 = %d (x, y) = (%d, %d), g = %d, s = %d" % (h1, x1, y1, g, s))
        for K in Ks:
            G, cnt = fixed(psi, h1, s, 1 << K, MODS)
            cells = []
            for i, m in enumerate(MODS):
                m = int(m)
                d = math.gcd(g, m)
                off = sum(int(cnt[i, c]) for c in range(m) if c % d != 0)
                if d == m:
                    cells.append("m=%d: d=m, all in class 0: %s" % (m, int(cnt[i, 0]) == G))
                    continue
                dev = max(abs((m // d) * int(cnt[i, c]) / G - 1.0) for c in range(0, m, d))
                noise = math.sqrt((m // d) / G)
                cells.append("m=%d (d=%d): dev %.4f (noise %.4f)%s" % (m, d, dev, noise, "" if off == 0 else " OFF-LATTICE %d" % off))
            print("     X = 2^%d: G = %9d;  %s" % (K, G, ";  ".join(cells)))


if __name__ == "__main__":
    main()
