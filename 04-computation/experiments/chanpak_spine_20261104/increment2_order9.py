#!/usr/bin/env python3
"""increment2_order9.py -- over all tournaments X of order 9 (191,536): does any cut give H(X + v) = H(X) + 2 with H(X) >= 5?
Records the smallest positive increment per X and the h-values where increment 2 occurs."""
import time
from ext_spectrum_true import tournaments, ext_spectrum
t0 = time.time(); inc2 = {}; mininc = {}; cnt = 0
for adj in tournaments(9):
    H = ext_spectrum(9, adj); h = int(H[0]); E = set(int(x) for x in H)
    pos = sorted(x for x in E if x > h); mp = (pos[0] - h) // 2 if pos else None
    cnt += 1
    if mp == 1: inc2[h] = inc2.get(h, 0) + 1
    mininc[h] = min(mininc.get(h, 10**9), mp if mp is not None else 10**9)
    if cnt % 20000 == 0: print("  ", cnt, "tournaments,", round(time.time() - t0), "s", flush=True)
print("order 9 (X): %d tournaments; increment-2 cases by h: %s" % (cnt, inc2))
print("smallest positive increment by h (h <= 60):", {h: mininc[h] for h in sorted(mininc) if h <= 60})
print("min positive increment over all X with h >= 5:", min(v for h, v in mininc.items() if h >= 5), "time", round(time.time() - t0), "s")
