#!/usr/bin/env python3
"""S93 coordinator probe: the thin sub-family of two-apex words in which the second apex v has ONE out-arc to the spine.

Spine s_0 -> ... -> s_(m-1) transitive, apexes u -> v, a_i = [u -> s_i] arbitrary, b = indicator of one position j.
These are the closed tournaments B => (trailing sinks) => v with v -> y over one-apex bases B rooted at a spine vertex
(closing law N = D + 2^t U).  Which odd N below 2^K are attained with m <= MMAX?
Usage: python3 single_arc_second_apex.py MMAX K
"""
import sys

import numpy as np
from numba import njit, prange

sys.path.insert(0, '.')
from eval_w1_w2_two_apex import step, final  # noqa: E402  (36-state transfer, checked against brute force)


@njit(parallel=True, cache=True)
def census(m, bound, mult):
    size = 1 << m
    for j in prange(m):
        vec = np.zeros(36, np.int64)
        out = np.zeros(36, np.int64)
        for A in range(size):
            for k in range(36):
                vec[k] = 0
            vec[0] = 1
            vec[18] = 1
            for i in range(m):
                a = (A >> (m - 1 - i)) & 1
                b = 1 if i == j else 0
                step(vec, a, b, out)
                for k in range(36):
                    vec[k] = out[k]
            h = final(vec)
            if h < bound:
                mult[j, h >> 1] += 1


def main():
    mmax = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    K = int(sys.argv[2]) if len(sys.argv) > 2 else 16
    bound = 1 << K
    total = np.zeros(bound >> 1, np.int64)
    for m in range(1, mmax + 1):
        mult = np.zeros((m, bound >> 1), np.int32)
        census(m, bound, mult)
        total += mult.sum(axis=0, dtype=np.int64)
        miss = np.flatnonzero(total == 0)
        missing = [int(2 * i + 1) for i in miss]
        first_big = [x for x in missing if x > 21][:12]
        covered_upto = missing[2] - 2 if len(missing) > 2 and missing[:2] == [7, 21] else None
        print(f"m <= {m:2d}: odd N < 2^{K} not attained: {len(missing)}; first misses: {missing[:8]}; "
              f"first misses above 21: {first_big[:6]}", flush=True)
    for k in range(4, K):
        seg = total[(1 << k) >> 1:(1 << (k + 1)) >> 1]
        print(f"  odd N in [2^{k}, 2^{k+1}): zeros {int((seg == 0).sum())} of {seg.shape[0]}, min multiplicity {int(seg.min())}, "
              f"mean {seg.mean():.1f}")


if __name__ == "__main__":
    main()
