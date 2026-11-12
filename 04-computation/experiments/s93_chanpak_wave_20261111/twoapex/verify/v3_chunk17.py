#!/usr/bin/env python3
"""Adversarial verifier, lane twoapex: the primitive words of length 17, in chunks of the prefix range (each run short).

usage: python v3_chunk17.py STATE_DIR CHUNK NCHUNKS THREADS
Writes STATE_DIR/P17_part{CHUNK}of{NCHUNKS}.npy (uint32 histogram of H >> 1 for H < 2^26) and a small .npz with the
number of words, the number of words with H >= 2^26 and the classes of H modulo 3 and 8.
"""
import os
import sys
import time

import numpy as np
from numba import set_num_threads

from v2_census import HALF, census, prefixes, suffixes


def main():
    sdir, ch, nch, T = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    set_num_threads(T)
    n = 17
    n2 = 6
    n1 = n - n2
    t0 = time.time()
    al = prefixes(n1)
    be = suffixes(n2)
    lo, hi = (al.shape[0] * ch) // nch, (al.shape[0] * (ch + 1)) // nch
    sub = np.ascontiguousarray(al[lo:hi])
    del al
    hist = np.zeros((T, HALF), np.uint16)
    cls = np.zeros((T, 24), np.int64)
    bigbuf = np.zeros((T, 16), np.int64)
    tc = time.time()
    nbig = census(sub, be, HALF, hist, cls, bigbuf)
    tc = time.time() - tc
    P = hist.sum(axis=0, dtype=np.uint32)
    words = sub.shape[0] * be.shape[0]
    nb = int(nbig.sum())
    assert int(P.sum(dtype=np.int64)) + nb == words
    assert int(cls.sum()) == words
    np.save(os.path.join(sdir, f"P17_part{ch}of{nch}.npy"), P)
    np.savez(os.path.join(sdir, f"P17_part{ch}of{nch}_small.npz"), words=words, nbig=nb, cls=cls.sum(axis=0))
    print(f"chunk {ch} of {nch}: prefixes [{lo}, {hi}), {words} primitive words of length 17, {nb} with H >= 2^26, "
          f"census {tc:.1f} s, total {time.time() - t0:.1f} s", flush=True)


if __name__ == "__main__":
    main()
