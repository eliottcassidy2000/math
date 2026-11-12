#!/usr/bin/env python3
"""Verifier: compare the verifier's own census tables (vk.exe, C) with the lane's tables, and recompute the
statistics quoted in PROOFS.txt section 5 from them."""
import numpy as np

CN = 1 << 20


def loadbin(fn):
    a = np.fromfile(fn, dtype=np.uint8)
    assert a.size == 2 * CN
    return a.reshape(2, CN)


lane = np.load("../b_census_minv.npz")
# 1. orders <= 9: cumulative minimum of my per-order tables against minv_9 of the lane
cum = np.full((2, CN), 255, np.uint8)
ok = True
for n in range(1, 10):
    cum = np.minimum(cum, loadbin(f"minv_small_{n}.bin"))
    same = bool((cum == lane[f"minv_{n}"]).all())
    ok &= same
print("V1  orders 1..9: my cumulative table minv (2 x 2^20 entries) = lane's minv_n for every n = 1..9:", ok)
print("     odd c <= 4097 with a seed / anti-seed of order <= 9 (mine):", int((cum[0][:2049] < 255).sum()), int((cum[1][:2049] < 255).sum()), "(lane, THM-5130 (E): 2016, 2013)")

# 2. order-10 slices
for P in (4, 13, 22, 31):
    mine = loadbin(f"minv_part_{P}_32.bin")
    d = np.load(f"../b_census_part_{P}_32.npz")
    print(f"V2  order 10, part {P} of 32: my table = lane's table: {bool((mine == d['minv']).all())}; entries set: {int((mine < 255).sum())}")

# 3. Theorem 5 (d): seed slopes below 450 from MY tables only (orders <= 9 complete + order 10 with H <= 224)
low = loadbin("minv_lowH_224.bin")
mine450 = np.minimum(cum, low)[0]
full10 = lane["minv_10"]


def slopes_below(mv, X):
    s = set()
    for c in range(1, X, 2):
        v = int(mv[c >> 1])
        if v == 255:
            continue
        x = c << v
        while x < X:
            s.add(x)
            x <<= 1
    return s


S_mine = slopes_below(mine450, 450)
S_lane = slopes_below(full10[0], 450)
print("V3  seed slopes below 450: mine", len(S_mine), "lane", len(S_lane), "equal:", S_mine == S_lane)
wins = [(1, 8), (8, 16), (16, 32), (32, 64), (64, 128), (128, 256), (256, 450)]
print("     windows:", [(a, b, sum(1 for x in S_mine if a <= x < b), b - a) for a, b in wins])
odd = sorted(x for x in S_mine if x & 1)
print("     odd seed slopes below 450:", len(odd), odd[:14], "...")
print("     small non-slopes (seeds):", [x for x in range(1, 40) if x not in S_mine])
# order needed: is every slope < 450 already present at order <= 10 of the lane's per-order tables (sanity)?
A_mine = slopes_below(np.minimum(cum, low)[1], 450)
A_lane = slopes_below(full10[1], 450)
print("     anti-seed slopes below 450 from bases of order <= 10 with H <= 224 (mine)", len(A_mine), "lane (order <= 10)", len(A_lane), "equal:", A_mine == A_lane,
      "[not a complete list of anti-seed slopes: see FINISHER R2]")

# 4. statistics of PROOFS section 5 from the lane's full table
for sg, name in ((0, "seeds"), (1, "anti-seeds")):
    mv = full10[sg]
    first = next(2 * i + 1 for i in range(1 << 16) if mv[i] == 255)
    mx = int(mv[: first >> 1].max())
    arg = [2 * i + 1 for i in range(first >> 1) if mv[i] == mx][:5]
    print(f"V4  {name}: first odd c without one (order <= 10): {first}; largest least-v below it: {mx} at c = {arg}")
    miss = [2 * i + 1 for i in range(1 << 15) if mv[i] == 255]
    print(f"     odd c <= 65535 without one: {len(miss)}: {miss[:6]} ... {miss[-3:]}")
    print("     odd M <= 2^k that are slopes themselves (v = 0), k = 13, 14, 15:", [int((mv[: (1 << k) // 2] == 0).sum()) for k in (13, 14, 15)])
    # integers in windows
    mvl = mv.astype(np.int64)
    shares = []
    for k in range(3, 16):
        tot = 0
        for b in range(0, k + 1):
            lo, hi = (1 << b) >> 1, (1 << (b + 1)) >> 1
            tot += int((mvl[lo:hi] <= k - b).sum())
        shares.append(round(tot / (1 << k), 4))
    print("     share of the integers of [2^k, 2^(k+1)) that are slopes, k = 3..15:", shares)

    def isslope(x):
        a = (x & -x).bit_length() - 1
        return mvl[(x >> a) >> 1] <= a
    hi = first - 1
    lo = hi
    while lo >= 1 and isslope(lo):
        lo -= 1
    print(f"     every integer of [{lo + 1}, {hi}] is a slope; non-slopes in [2^14, 2^15):", [x for x in range(1 << 14, 1 << 15) if not isslope(x)])

# 5. odd N on a unit ray
mvp = full10[0].astype(np.int64)
mvm = full10[1].astype(np.int64)


def onray(N):
    x = N - 1
    s = (x & -x).bit_length() - 1
    y = N + 1
    s2 = (y & -y).bit_length() - 1
    return s >= mvp[(x >> s) >> 1] + 1 or s2 >= mvm[(y >> s2) >> 1] + 1


unc = [N for N in range(3, 4096, 2) if not onray(N)]
print("V5  odd N < 4096 on no unit ray of a pair of order <= 10:", len(unc), unc[:8])
for k in (12, 13, 14, 15, 16):
    c = sum(1 for N in range((1 << k) + 1, 1 << (k + 1), 2) if onray(N))
    print(f"     odd N in [2^{k}, 2^{k+1}) on a unit ray: {c} of {1 << (k - 1)}")
