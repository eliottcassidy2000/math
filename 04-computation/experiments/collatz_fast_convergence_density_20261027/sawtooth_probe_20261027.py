#!/usr/bin/env python3
"""
sawtooth_probe_20261027.py  (opus-2026-10-27-S64) -- POST-HOC probe, not pre-registered.

The pre-registered predictions P3 and P4 of fast_convergence_density_20261027.py failed on a period-two
oscillation in the dyadic index k: the floor count #{n in [2^(k-1), 2^k) : tau(n) < 3.307 ln n} grows by a factor
about 1.35 from an even k to the next odd k and by about 1.9-2.1 from an odd k to the next even k, and the block
density G_k(13) is not monotone. This probe asks what the period of the oscillation is in ln n: it recomputes the
floor count and the densities G(8), G(13) on
  (a) quarter-dyadic blocks [2^(k-1) (1 + i/4), 2^(k-1) (1 + (i+1)/4)), i = 0..3, k = 20..28 (resolution ln 2 / 4);
  (b) blocks of ratio 3, [3^j, 3^(j+1)), 3^j < 2^28;
  (c) blocks of ratio 4, [4^j, 4^(j+1)).
If the oscillation has period ln 4 in ln n (the trunk entries (4^j - 1)/3 land on 4^j), the 4-adic blocks (c) show
a smooth sequence of ratios, the quarter blocks (a) show a period of 8 quarter-steps, and the 3-adic blocks (b)
show beats. Output: counts and ratios; no verdict is pre-registered.
"""
import math
import sys

import numpy as np

sys.path.insert(0, __file__.rsplit("\\", 1)[0] if "\\" in __file__ else (__file__.rsplit("/", 1)[0] or "."))
from fast_convergence_density_20261027 import tau_all, FLOOR_C  # noqa: E402

KMAX = 28
CHUNK = 1 << 22


def stats(tau, lo, hi):
    """(count, floor count, #tau <= 8 ln n, #tau <= 13 ln n) on [lo, hi)."""
    fl = c8 = c13 = 0
    for a in range(lo, hi, CHUNK):
        b = min(hi, a + CHUNK)
        n = np.arange(a, b, dtype=np.float64)
        ln = np.log(n)
        t = tau[a:b].astype(np.float64)
        fl += int(np.count_nonzero(t < FLOOR_C * ln))
        c8 += int(np.count_nonzero(t <= 8.0 * ln))
        c13 += int(np.count_nonzero(t <= 13.0 * ln))
    return hi - lo, fl, c8, c13


def main():
    X = 1 << KMAX
    tau, _ = tau_all(X)
    print("(a) quarter-dyadic blocks: k, i, lower bound, size, floor#, G(8), G(13), floor density x 1e4")
    prev = None
    for k in range(20, KMAX + 1):
        base = 1 << (k - 1)
        for i in range(4):
            lo = base + (base * i) // 4
            hi = base + (base * (i + 1)) // 4
            cnt, fl, c8, c13 = stats(tau, lo, hi)
            ratio = (fl / cnt) / prev if prev else float("nan")
            print("%3d %d %12d %10d %8d  %.4f  %.4f  %.3f  (ratio to previous quarter %.3f)"
                  % (k, i, lo, cnt, fl, c8 / cnt, c13 / cnt, 1e4 * fl / cnt, ratio))
            prev = fl / cnt
    print()
    print("(b) blocks of ratio 3: j, [3^j, 3^(j+1)), size, floor#, floor density x 1e4, G(8), G(13), count ratio")
    prevfl = None
    j = 10
    while 3 ** (j + 1) <= X:
        lo, hi = 3 ** j, 3 ** (j + 1)
        cnt, fl, c8, c13 = stats(tau, lo, hi)
        print("%3d %12d %12d %10d %8d  %.3f  %.4f  %.4f  %s"
              % (j, lo, hi, cnt, fl, 1e4 * fl / cnt, c8 / cnt, c13 / cnt,
                 ("%.3f" % (fl / prevfl)) if prevfl else "nan"))
        prevfl = fl
        j += 1
    print()
    print("(c) blocks of ratio 4: j, [4^j, 4^(j+1)), size, floor#, floor density x 1e4, G(8), G(13), count ratio")
    prevfl = None
    j = 8
    while 4 ** (j + 1) <= X:
        lo, hi = 4 ** j, 4 ** (j + 1)
        cnt, fl, c8, c13 = stats(tau, lo, hi)
        print("%3d %12d %12d %10d %8d  %.3f  %.4f  %.4f  %s"
              % (j, lo, hi, cnt, fl, 1e4 * fl / cnt, c8 / cnt, c13 / cnt,
                 ("%.3f" % (fl / prevfl)) if prevfl else "nan"))
        prevfl = fl
        j += 1
    print()
    print("(d) dyadic blocks once more, with the two-block ratio (k vs k-2) and the implied exponent log_4(ratio):")
    fls = {}
    for k in range(16, KMAX + 1):
        cnt, fl, c8, c13 = stats(tau, 1 << (k - 1), 1 << k)
        fls[k] = fl
        if k >= 18:
            r2 = fl / fls[k - 2]
            print("%3d floor# %8d  two-block ratio %.3f  exponent %.3f" % (k, fl, r2, math.log(r2, 4)))
        else:
            print("%3d floor# %8d" % (k, fl))


if __name__ == "__main__":
    main()
