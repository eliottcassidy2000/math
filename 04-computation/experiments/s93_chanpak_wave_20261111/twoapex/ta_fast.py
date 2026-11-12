#!/usr/bin/env python3
"""Lane twoapex: fast (numba) machinery on the 24-state aggregated transfer.

State of a prefix, per arrangement arr in {0, 1}: (c1, r2, r3) with
  c1 = 1 iff R1 is nonempty and its last vertex does not beat x (R1 needs one more vertex), else 0;
  r2 in {0 empty, 1 last ok, 2 last bad}; r3 in {0 empty, 1 nonempty}.
index = arr * 12 + c1 * 6 + r2 * 2 + r3.  The 36-state transfer of ta_core.py projects onto it (r1 in {0, 1} merge).
alpha(w)[s] = number of 3-colourings of the prefix w in state s;  beta(w)[s] = number of completions of state s by
the suffix w;  H(w1 w2) = sum_s alpha(w1)[s] beta(w2)[s]   (all terms nonnegative: no overflow below 2^63).
"""
import numpy as np
from numba import njit, prange

NS = 24


def build24():
    M = np.zeros((4, NS, NS), np.int64)
    for l in range(4):
        a, b = l & 1, (l >> 1) & 1
        for arr in range(2):
            al, be = (a, b) if arr == 0 else (b, a)
            for c1 in range(2):
                for r2 in range(3):
                    for r3 in range(2):
                        s = arr * 12 + c1 * 6 + r2 * 2 + r3
                        M[l, s, arr * 12 + (1 if al == 1 else 0) * 6 + r2 * 2 + r3] += 1
                        if r2 != 0 or al == 1:
                            M[l, s, arr * 12 + c1 * 6 + (1 if be == 0 else 2) * 2 + r3] += 1
                        if r3 != 0 or be == 1:
                            M[l, s, arr * 12 + c1 * 6 + r2 * 2 + 1] += 1
    init = np.zeros(NS, np.int64)
    init[0] = 1
    init[12] = 1
    fin = np.zeros(NS, np.int64)
    for arr in range(2):
        for r2 in (0, 1):
            if r2 == 0 and arr == 1:
                continue
            for r3 in range(2):
                fin[arr * 12 + 0 * 6 + r2 * 2 + r3] = 1
    return M, init, fin


M24, INIT24, FIN24 = build24()


def H24(word):
    v = INIT24.copy()
    for l in word:
        v = v @ M24[l]
    return int(v @ FIN24)


@njit(cache=True)
def all_alpha(n, M, init):
    """alpha vectors of all words of length n; index of w = l_1 ... l_n is sum l_i 4^(n-i) (first letter = top digit)."""
    cur = np.zeros((1, NS), np.int64)
    cur[0, :] = init
    for step in range(n):
        m = cur.shape[0]
        nxt = np.zeros((m * 4, NS), np.int64)
        for i in range(m):
            for l in range(4):
                for s in range(NS):
                    c = cur[i, s]
                    if c != 0:
                        for t in range(NS):
                            if M[l, s, t] != 0:
                                nxt[i * 4 + l, t] += c * M[l, s, t]
        cur = nxt
    return cur


@njit(cache=True)
def all_beta(n, M, fin):
    """beta vectors of all words of length n; same indexing (first letter = top digit)."""
    cur = np.zeros((1, NS), np.int64)
    cur[0, :] = fin
    size = 1
    for step in range(n):
        m = cur.shape[0]
        nxt = np.zeros((m * 4, NS), np.int64)
        # prepend letter l: beta(l w)[s] = sum_t M[l, s, t] beta(w)[t]; index of l w = l * size + idx(w)
        for l in range(4):
            for i in range(m):
                for s in range(NS):
                    acc = 0
                    for t in range(NS):
                        if M[l, s, t] != 0:
                            acc += M[l, s, t] * cur[i, t]
                    nxt[l * size + i, s] = acc
        cur = nxt
        size *= 4
    return cur


@njit(cache=True)
def compress(vecs):
    """Drop the coordinates that vanish on every vector; returns the list of kept coordinates."""
    keep = np.zeros(NS, np.bool_)
    for i in range(vecs.shape[0]):
        for s in range(NS):
            if vecs[i, s] != 0:
                keep[s] = True
    return keep


@njit(parallel=True, cache=True)
def census_pairs(al, be, bound, nchunk, modulus, nbin, scale_num, scale_den):
    """All words w1 w2 (w1 over the rows of al, w2 over the rows of be).
    Returns per-chunk arrays: hist[c, N >> 1] for odd N < bound, modhist[c, N % modulus], bins[c, .] of
    floor(nbin * N * scale_num / scale_den) (value distribution, clipped), and big[c], mx[c], arg1[c], arg2[c]."""
    n1 = al.shape[0]
    n2 = be.shape[0]
    d = al.shape[1]
    hist = np.zeros((nchunk, bound >> 1), np.uint32)
    modhist = np.zeros((nchunk, modulus), np.int64)
    bins = np.zeros((nchunk, nbin + 1), np.int64)
    big = np.zeros(nchunk, np.int64)
    mx = np.zeros(nchunk, np.int64)
    arg1 = np.zeros(nchunk, np.int64)
    arg2 = np.zeros(nchunk, np.int64)
    for c in prange(nchunk):
        lo = (n1 * c) // nchunk
        hi = (n1 * (c + 1)) // nchunk
        for i in range(lo, hi):
            for j in range(n2):
                h = 0
                for s in range(d):
                    h += al[i, s] * be[j, s]
                if h < bound:
                    hist[c, h >> 1] += 1
                else:
                    big[c] += 1
                modhist[c, h % modulus] += 1
                b = (h * scale_num * nbin) // scale_den
                if b > nbin:
                    b = nbin
                bins[c, b] += 1
                if h > mx[c]:
                    mx[c] = h
                    arg1[c] = i
                    arg2[c] = j
    return hist, modhist, bins, big, mx, arg1, arg2


def word_of_index(idx, n):
    return [(idx >> (2 * (n - 1 - i))) & 3 for i in range(n)]
