#!/usr/bin/env python3
"""Lane twoapex, deliverable (b): exhaustive census of the two-apex words by length (meet in the middle).

usage: python b1_census.py LSTART LEND BOUND_EXP STATE_DIR
For n = LSTART..LEND: all 4^n words w = w1 w2 (|w1| = ceil(n/2)), H = <alpha(w1), beta(w2)> (24-state transfer),
  c_n(N) = number of words of length exactly n with H = N, for odd N < 2^BOUND_EXP;
  cum(N) = sum_{m <= n} c_m(N);  minlen(N) = least n with c_n(N) > 0;
  histogram of H modulo 720720 over ALL words of length n; histogram of H / (2 * 3^n) in 4096 bins.
The state (cum, last, prev, minlen) is kept in STATE_DIR so that the lengths can be run in separate calls.
Per-length summaries are appended to census_by_length.jsonl; the residue histograms go to STATE_DIR.
"""
import json
import os
import sys
import time

import numpy as np

from ta_fast import FIN24, INIT24, M24, all_alpha, all_beta, census_pairs, word_of_index

MOD = 720720            # 16 * 9 * 5 * 7 * 11 * 13
NBIN = 4096


def main():
    lstart, lend, bexp, sdir = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    bound = 1 << bexp
    half = bound >> 1
    nchunk = int(os.environ.get("TA_CHUNKS", "12"))
    spath = os.path.join(sdir, f"census_state_b{bexp}.npz")
    if lstart > 1:
        st = np.load(spath)
        assert int(st["upto"]) == lstart - 1, (int(st["upto"]), lstart)
        cum, last, minlen = st["cum"], st["last"], st["minlen"]
    else:
        cum = np.zeros(half, np.uint32)
        last = np.zeros(half, np.uint32)
        minlen = np.zeros(half, np.uint8)
        # the empty word (n = 0): H = 1
        cum[0] = 1
        last[0] = 1
        minlen[0] = 0
    for n in range(lstart, lend + 1):
        t0 = time.time()
        n1, n2 = (n + 1) // 2, n // 2
        al = all_alpha(n1, M24, INIT24)
        be = all_beta(n2, M24, FIN24)
        keep = (al != 0).any(axis=0) & (be != 0).any(axis=0)
        al = np.ascontiguousarray(al[:, keep])
        be = np.ascontiguousarray(be[:, keep])
        hist, modhist, bins, big, mx, a1, a2 = census_pairs(al, be, bound, nchunk, MOD, NBIN, 1, 2 * 3 ** n)
        prev = last
        last = hist.sum(axis=0, dtype=np.uint32)
        del hist
        mh = modhist.sum(axis=0)
        bn = bins.sum(axis=0)
        newly = (cum == 0) & (last != 0)
        minlen[newly] = n
        cum = cum + last
        k = int(np.argmax(mx))
        word = word_of_index(int(a1[k]), n1) + word_of_index(int(a2[k]), n2)
        rec = {
            "n": n, "words": 4 ** n, "ge_bound": int(big.sum()), "bound_exp": bexp,
            "max_H": int(mx.max()), "max_word": "".join(map(str, word)),
            "max_over_2_3n": float(mx.max()) / (2 * 3 ** n),
            "distinct_odd_values_below_bound_exact_length": int(np.count_nonzero(last)),
            "distinct_odd_values_below_bound_cumulative": int(np.count_nonzero(cum)),
            "mean_H": float((np.arange(NBIN + 1) + 0.5) @ bn) / NBIN * 2 * 3 ** n / 4 ** n,
            "kept_coordinates": int(keep.sum()),
            "seconds": round(time.time() - t0, 1),
        }
        assert int(mh.sum()) == 4 ** n and int(last.sum(dtype=np.int64)) + int(big.sum()) == 4 ** n
        np.save(os.path.join(sdir, f"modhist_n{n:02d}.npy"), mh)
        np.save(os.path.join(sdir, f"valbins_n{n:02d}.npy"), bn)
        with open("census_by_length.jsonl", "a") as f:
            f.write(json.dumps(rec) + "\n")
        print(json.dumps(rec), flush=True)
        np.savez(spath, cum=cum, last=last, prev=prev, minlen=minlen, upto=n)
    print("state saved:", spath)


if __name__ == "__main__":
    main()
