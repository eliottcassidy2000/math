#!/usr/bin/env python3
"""S93 lane seeds: two-apex bases, seeds (r = +1): for each spine length n, the number of distinct slopes that are
new at length n (not a slope of a shorter word), by dyadic window [2^k, 2^(k+1)); cumulative share of each window
(doubling closure is automatic: a sink doubles the slope and lengthens the word by one letter)."""
import numpy as np

recs = {}
d9 = np.load("d_two_9.npz")
for n in range(1, 10):
    recs[n] = d9[f"rec_{n}"]
recs[10] = np.load("d_two_10_10.npz")["rec_10"]
for sgn, name in ((1, "seeds"), (-1, "anti-seeds")):
    seen = set()
    print(name)
    print("   n   pairs  new slopes | new slopes by window k = 3..19")
    cum = np.zeros(21, np.int64)
    for n in range(1, 11):
        R = recs[n]
        sub = R[R[:, 2] == sgn]
        sl = set(sub[:, 3].tolist())
        new = sl - seen
        seen |= sl
        hist = np.zeros(21, np.int64)
        for M in new:
            hist[M.bit_length() - 1] += 1
        cum += hist
        print(f"  {n:2d} {len(sub):7d} {len(new):7d} | " + " ".join(f"{int(hist[k]):5d}" for k in range(3, 20)))
    print("  cumulative share of the window: " + " ".join(f"{cum[k] / (1 << k):5.3f}" for k in range(3, 20)))
    # median and quartiles of log2(slope) of the pairs of length n
    for n in range(4, 11):
        R = recs[n]
        sub = R[R[:, 2] == sgn][:, 3]
        lg = np.log2(sub.astype(np.float64))
        print(f"  n = {n}: log2(slope) of the pairs: quartiles {np.percentile(lg, 25):.2f} {np.percentile(lg, 50):.2f} {np.percentile(lg, 75):.2f}, max {lg.max():.2f}")
