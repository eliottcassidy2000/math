#!/usr/bin/env python3
"""Adversarial verifier, lane twoapex: array-level comparison of the verifier's independent census (primitive words,
own transfer; state after length 16) with the lane's census state after length 16 (all words, lane's transfer).

usage: python v4_compare16.py VERIFIER_STATE_DIR LANE_STATE_DIR
"""
import os
import pickle
import sys

import numpy as np


def main():
    vdir, ldir = sys.argv[1], sys.argv[2]
    small = pickle.load(open(os.path.join(vdir, "small.pkl"), "rb"))
    assert small["n"] == 16
    st = np.load(os.path.join(ldir, "census_state_b26_L16.npz"))
    print("lane state: upto =", int(st["upto"]))
    for mine, theirs, what in (("cum.npy", "cum", "cum_16 (all words of length <= 16)"),
                               ("c_n.npy", "last", "c_16 (all words of length exactly 16)"),
                               ("minlen.npy", "minlen", "minimal length")):
        a = np.load(os.path.join(vdir, mine))
        b = st[theirs]
        if theirs == "minlen":
            att = a != 255
            same = bool((a[att] == b[att]).all()) and bool(((b == 0) | att)[1:].all())
            print(f"{what}: equal on the attained N: {same}; attained {int(att.sum())}")
        else:
            print(f"{what}: arrays identical: {bool((a.astype(np.int64) == b.astype(np.int64)).all())}; "
                  f"sum {int(a.sum(dtype=np.int64))} vs {int(b.sum(dtype=np.int64))}")
        del a, b
    cp = np.load(os.path.join(vdir, "cumprim.npy")).astype(np.int64)
    prim_lane = st["last"].astype(np.int64) - st["prev"].astype(np.int64)
    print(f"cumprim_16 (primitive words of length <= 16, enumerated directly) = lane's c_16 - c_15: "
          f"{bool((cp == prim_lane).all())}; sum {int(cp.sum())} vs {int(prim_lane.sum())}")


if __name__ == "__main__":
    main()
