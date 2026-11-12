#!/usr/bin/env python3
"""Verifier: the table entry 'two-apex spine G(2,0): every odd < 2^21 but 7, 21 at length <= 14 (seed)' rests on a
scratch result that the seed lane did not store.  Recompute it with the lane's own 36 x 36 letters (two_apex_algebra.
step_matrix, checked there against direct counts): all words of length <= NMAX, pruned when H >= LIM (H is monotone
under extension of the word: the prefix tournament is a vertex-deleted subtournament and H(T) >= H(T - v)).
Usage: v_twoapex_cover.py NMAX LOG2LIM
"""
import sys
import numpy as np
from numba import njit

sys.path.insert(0, "C:/tmp/wt-s93/04-computation/experiments/s93_chanpak_wave_20261111/thinorbit")
from two_apex_algebra import step_matrix, final_vec, brute_two_apex


@njit(cache=True)
def run(dst, cnt, f, v0, nmax, lim, hit):
    # dst[l, src, j] (j < cnt[l, src]) = target rows of column src of letter l (with multiplicity)
    vec = np.zeros((nmax + 1, 36), np.int64)
    for i in range(36):
        vec[0, i] = v0[i]
    letter = np.zeros(nmax + 1, np.int64)
    depth = 0
    letter[0] = 0
    nodes = 0
    while True:
        if depth < nmax and letter[depth] < 4:
            l = letter[depth]
            letter[depth] += 1
            for i in range(36):
                vec[depth + 1, i] = 0
            for src in range(36):
                c = vec[depth, src]
                if c != 0:
                    for j in range(cnt[l, src]):
                        vec[depth + 1, dst[l, src, j]] += c
            H = 0
            for i in range(36):
                H += f[i] * vec[depth + 1, i]
            nodes += 1
            if H < lim:
                hit[H] = 1
                depth += 1
                letter[depth] = 0
        else:
            if depth == 0:
                break
            depth -= 1
    return nodes


def main():
    nmax = int(sys.argv[1])
    lim = 1 << int(sys.argv[2])
    letters = [(0, 0), (1, 0), (0, 1), (1, 1)]
    Ls = [step_matrix(a, b) for (a, b) in letters]
    dst = np.zeros((4, 36, 8), np.int64)
    cnt = np.zeros((4, 36), np.int64)
    for l in range(4):
        for src in range(36):
            for d in range(36):
                for _ in range(int(Ls[l][d, src])):
                    dst[l, src, cnt[l, src]] = d
                    cnt[l, src] += 1
    f = final_vec()
    v0 = np.zeros(36, np.int64)
    v0[0] = 1
    v0[18] = 1
    # spot check of the transfer against a direct count (random words of length 7, 8)
    rng = np.random.default_rng(11)
    bad = 0
    for n in (7, 7, 8, 8, 8):
        wd = [int(x) for x in rng.integers(0, 4, n)]
        v = v0.copy()
        for l in wd:
            v = Ls[l].dot(v)
        if int(f.dot(v)) != brute_two_apex(n, [letters[l][0] for l in wd], [letters[l][1] for l in wd]):
            bad += 1
    print("spot check transfer = direct count on 5 random words of length 7, 8: failures", bad)
    hit = np.zeros(lim, np.uint8)
    nodes = run(dst, cnt, f, v0, nmax, lim, hit)
    miss = [v for v in range(1, lim, 2) if not hit[v]]
    print(f"words of length <= {nmax} (pruned at H >= 2^{sys.argv[2]}): {nodes} nodes; odd v < 2^{sys.argv[2]} not attained: "
          f"{len(miss)}; first {miss[:20]}; even values attained: {int(hit[0::2].sum())}")
    for k in range(3, int(sys.argv[2])):
        lo, hi = 1 << k, 1 << (k + 1)
        z = sum(1 for v in range(lo + 1, hi, 2) if not hit[v])
        if z:
            print(f"   window [2^{k}, 2^{k+1}): {z} odd values not attained")


if __name__ == "__main__":
    main()
