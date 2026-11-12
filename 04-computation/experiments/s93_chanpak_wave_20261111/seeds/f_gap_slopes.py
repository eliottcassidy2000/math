#!/usr/bin/env python3
"""S93 lane seeds, FINISHER: the integers below 2^16 that are not seed (anti-seed) slopes of any pair of order <= 10
(b_small_slopes.out / f_recheck.out F3): are they slopes of the two-apex bases of spine length <= 10 (orders up to 12,
d_two_9.npz, d_two_10_10.npz) or of the one-apex bases of length <= 14 (d_one_14.npz)?  (Closure under doubling: an
integer M is a slope iff M / 2^j is a slope of some pair for some j >= 0.)"""
import numpy as np

mv = np.load("b_census_minv.npz")["minv_10"].astype(np.int64)
two = {}
d9 = np.load("d_two_9.npz")
for n in range(1, 10):
    two[n] = d9[f"rec_{n}"]
two[10] = np.load("d_two_10_10.npz")["rec_10"]
one = np.load("d_one_14.npz")
for sg, sgn, name in ((0, 1, "seed"), (1, -1, "anti-seed")):
    def is_slope(M):
        v = (M & -M).bit_length() - 1
        return mv[sg][(M >> v) >> 1] <= v
    S2 = set()
    for n in two:
        R = two[n]
        S2 |= set(R[R[:, 2] == sgn][:, 3].tolist())
    S1 = set()
    for n in range(1, 15):
        R = one[f"rec_{n}"]
        S1 |= set(R[R[:, 2] == sgn][:, 3].tolist())
    def in_closure(M, S):
        while True:
            if M in S:
                return True
            if M & 1:
                return False
            M >>= 1
    for lo, hi in ((1 << 9, 1 << 13), (1 << 13, 1 << 14), (1 << 14, 1 << 15), (1 << 15, 1 << 16)):
        non = [M for M in range(lo, hi) if not is_slope(M)]
        f2 = [M for M in non if in_closure(M, S2)]
        f1 = [M for M in non if in_closure(M, S1)]
        print(f"{name}: [{lo}, {hi}): {len(non)} integers are not slopes of a pair of order <= 10; of these {len(f2)} are slopes of two-apex bases (length <= 10), {len(f1)} of one-apex bases (length <= 14)"
              + (f"; two-apex: {f2}" if len(non) <= 30 else f"; first two-apex ones: {f2[:10]}"))
