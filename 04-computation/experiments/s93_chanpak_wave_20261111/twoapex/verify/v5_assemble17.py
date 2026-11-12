#!/usr/bin/env python3
"""Adversarial verifier, lane twoapex: assemble the independent census at length 17 from the verifier's state after
length 16 and the chunks of v3_chunk17.py; statistics; array-level comparison with the lane's state after length 17.

usage: python v5_assemble17.py VERIFIER_STATE_DIR NCHUNKS LANE_STATE_DIR
"""
import os
import pickle
import sys

import numpy as np

BEXP = 26
HALF = 1 << (BEXP - 1)


def main():
    vdir, nch, ldir = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    small = pickle.load(open(os.path.join(vdir, "small.pkl"), "rb"))
    assert small["n"] == 16
    P = np.zeros(HALF, np.uint32)
    words = nb = 0
    pc = np.zeros(24, np.int64)
    for ch in range(nch):
        P += np.load(os.path.join(vdir, f"P17_part{ch}of{nch}.npy"))
        s = np.load(os.path.join(vdir, f"P17_part{ch}of{nch}_small.npz"))
        words += int(s["words"])
        nb += int(s["nbig"])
        pc += s["cls"]
    n = 17
    assert words == 9 * 4 ** (n - 2), words
    assert int(P.sum(dtype=np.int64)) + nb == words
    cumprim = np.load(os.path.join(vdir, "cumprim.npy"))
    c_n = np.load(os.path.join(vdir, "c_n.npy"))
    cum = np.load(os.path.join(vdir, "cum.npy"))
    minlen = np.load(os.path.join(vdir, "minlen.npy"))
    cumprim += P
    c_n += cumprim
    newly = (cum == 0) & (c_n != 0)
    minlen[newly] = n
    cum += c_n
    cumcls = small["cumcls"] + pc
    R = small["R"] + cumcls
    cumbig = small["cumbig"] + nb
    cbig = small["cbig"] + cumbig
    print(f"n = 17: primitive words {words}, with H >= 2^26: {nb} (expected from the lane's figures 301885533 - 2 * 667 = "
          f"{301885533 - 2 * 667}); all words of length 17 with H >= 2^26: {cbig} (lane: 301885533)")
    print(f"  all 4^17 words accounted for (value histogram + above the bound): {int(c_n.sum(dtype=np.int64)) + cbig == 4 ** n}; "
          f"class recursion gives 4^17 words: {int(R.sum()) == 4 ** n}")
    fac3 = [3 * float(R.reshape(3, 8)[c].sum()) / 4 ** n for c in range(3)]
    r8 = R.reshape(3, 8).sum(axis=0)
    print("  by-index class factors of H modulo 3 over all 4^17 words: " + " ".join(f"{t:.6f}" for t in fac3)
          + "  (lane: 1.003721 1.002709 0.993570); odd classes modulo 8: "
          + " ".join(f"{4 * float(r8[c]) / 4 ** n:.6f}" for c in (1, 3, 5, 7)))
    c64 = c_n.astype(np.int64)
    E = int((c64 * c64).sum())
    print(f"  distinct odd values < 2^26 at exact length 17: {int(np.count_nonzero(c64))}; cumulative: {int(np.count_nonzero(cum))} "
          f"(2^25 - 2 = {HALF - 2})")
    print(f"  E_17 >= {E} (values >= 2^26 not counted); E_17 (3/16)^17 >= {E * (3 / 16) ** n:.5f}; "
          f"max_(N < 2^26) c_17(N) = {int(c64.max())} at N = {2 * int(np.argmax(c64)) + 1}; max/(4/3)^17 = {c64.max() / (4 / 3) ** n:.2f}")
    del c64
    miss = 2 * np.flatnonzero(cum == 0) + 1
    print(f"  odd N < 2^26 not attained by words of length <= 17: {len(miss)} -> {miss[:10].tolist()}")
    print(f"  total multiplicity below the bound {int(cum.sum(dtype=np.int64))} (lane: 22604606045); primitive words recorded "
          f"{int(cumprim.sum(dtype=np.int64))} (lane: 12583017022)")
    print("      window          zeros      min   median       mean | prim: min  median     mean | max minlen")
    for k in range(18, BEXP):
        lo, hi = (1 << k) >> 1, (1 << (k + 1)) >> 1
        seg, ps, ml = cum[lo:hi], cumprim[lo:hi], minlen[lo:hi]
        print(f"      [2^{k}, 2^{k + 1}) {int((seg == 0).sum()):9d} {int(seg.min()):8d} {int(np.median(seg)):8d} "
              f"{seg.mean():10.1f} | {int(ps.min()):9d} {int(np.median(ps)):7d} {ps.mean():8.1f} | {int(ml.max())}")
    print("      number of odd N of minimal length m:", [(m, int((minlen == m).sum())) for m in range(13, 18)])
    print("  --- dispersion of the primitive multiplicities (length <= 17), blocks of 256 consecutive odd N")
    for k in range(18, BEXP):
        lo, hi = (1 << k) >> 1, (1 << (k + 1)) >> 1
        seg = cumprim[lo:hi].astype(np.float64)
        blk = seg.reshape(-1, 256)
        res = blk - blk.mean(axis=1, keepdims=True)
        var = float((res ** 2).sum() / (blk.shape[0] * 255))
        lag = float((res[:, 1:] * res[:, :-1]).sum() / (res ** 2).sum())
        odd = float((cumprim[lo:hi] & 1).mean())
        blk2 = seg.reshape(-1, 4096)
        var2 = float(((blk2 - blk2.mean(axis=1, keepdims=True)) ** 2).sum() / (blk2.shape[0] * 4095))
        gvar = float(seg.var())
        N = 2 * np.arange(lo, hi, dtype=np.int64) + 1
        m3 = [float(seg[N % 3 == c].mean()) / float(seg.mean()) for c in range(3)]
        print(f"      [2^{k}, 2^{k + 1}): mean {seg.mean():9.2f}  var/mean {var / seg.mean():.4f}  of prim/2 {var / seg.mean() / 2:.4f}"
              f"  (blocks of 4096: {var2 / seg.mean() / 2:.4f}; whole window: {gvar / seg.mean() / 2:.3f})  lag-1 {lag:+.4f}  "
              f"odd share {odd:.6f}  min {int(seg.min())}  class factors of N mod 3: " + " ".join(f"{t:.5f}" for t in m3))
    # comparison with the lane's state after length 17
    st = np.load(os.path.join(ldir, "census_state_b26_L17.npz"))
    print("lane state: upto =", int(st["upto"]))
    print("  cum_17 identical:", bool((cum == st["cum"]).all()))
    print("  c_17 identical:", bool((c_n == st["last"]).all()))
    print("  primitive counts: verifier's direct cumprim_17 = lane's c_17 - c_16:",
          bool((cumprim.astype(np.int64) == st["last"].astype(np.int64) - st["prev"].astype(np.int64)).all()))
    att = minlen != 255
    print("  minimal lengths identical on the attained N:", bool((minlen[att] == st["minlen"][att]).all()),
          "; attained:", int(att.sum()))


if __name__ == "__main__":
    main()
