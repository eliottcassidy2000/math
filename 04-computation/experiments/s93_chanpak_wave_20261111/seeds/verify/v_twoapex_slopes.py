#!/usr/bin/env python3
"""Verifier: seed / anti-seed slopes (< 2^17) of the two-apex bases of spine length <= 10, from my own enumeration
(vk.exe twoapex), against PROOFS.txt section 7 and the finisher's f_gap_slopes.out."""
import os
import numpy as np

SL = 1 << 17
sl = np.zeros((2, SL), np.uint8)
files = [f"two_slopes_{n}.bin" for n in range(1, 10)] + ["two_slopes_10a.bin", "two_slopes_10b.bin"]
missing = [f for f in files if not os.path.exists(f)]
print("files missing:", missing)
for f in files:
    if os.path.exists(f):
        sl |= np.fromfile(f, dtype=np.uint8).reshape(2, SL)
raw = sl.copy()
clo = sl.copy()
for x in range(1, SL // 2):
    for sg in range(2):
        if clo[sg, x]:
            clo[sg, 2 * x] = 1
for name, a in (("raw", raw), ("doubling closure", clo)):
    print(f"two-apex seed slopes, length <= 10 ({name}): share of [2^k, 2^(k+1)), k = 6..15:",
          [round(float(a[0, 1 << k: 1 << (k + 1)].mean()), 3) for k in range(6, 16)])
    print(f"   counts in [64, 128), [128, 256), [256, 512): {[int(a[0, lo:hi].sum()) for lo, hi in ((64, 128), (128, 256), (256, 512))]}  (PROOFS: 29, 53, 135)")

lane = np.load("../b_census_minv.npz")["minv_10"].astype(np.int64)


def isslope(sg, x):
    a = (x & -x).bit_length() - 1
    return lane[sg][(x >> a) >> 1] <= a


for sg, nm in ((0, "seed"), (1, "anti-seed")):
    for lo, hi in ((512, 8192), (8192, 16384), (16384, 32768), (32768, 65536)):
        gaps = [x for x in range(lo, hi) if not isslope(sg, x)]
        hit_raw = [x for x in gaps if raw[sg, x]]
        hit_clo = [x for x in gaps if clo[sg, x]]
        print(f"{nm}: [{lo}, {hi}): {len(gaps)} integers are not slopes of a pair of order <= 10; two-apex (length <= 10) slopes among them: raw {len(hit_raw)}, closure {len(hit_clo)}; first {hit_raw[:6]}")
# every census seed slope below 512 that is not a two-apex slope
print("census seed slopes in [256, 512) that are not two-apex slopes (length <= 10):", [x for x in range(256, 512) if isslope(0, x) and not clo[0, x]])
