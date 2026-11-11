# core.py - numba kernels for the rooted spectrum (helper of session opus-2026-11-10-S91).
#
# Conventions.  A tournament on vertices 0..n-1 is coded by an integer: for i < j the bit with index
# j(j-1)/2 + i is 1 iff i -> j (colex order of pairs, so that adding the vertex n appends n bits).
# out[v] is the bit mask of the out-neighbours of v.
# f[S, v] = number of Hamiltonian paths of T[S] ending at v, g[S, v] = number starting at v,
# hS[S] = number of Hamiltonian paths of T[S], hS[0] = 1.
# Rooted invariants of (T, v):
#   H      = hS[full]
#   s      = g[full, v]                      paths starting at v
#   e      = f[full, v]                      paths ending at v
#   Yplus  = sum_{S contains v} f[S, v] hS[full - S]   ordered pairs (P1, P2), P1 ends at v, P2 possibly empty
#   Yminus = sum_{S contains v} g[S, v] hS[full - S]   ordered pairs (P1, P2), P2 starts at v, P1 possibly empty
#   d      = hS[full - v]                    H(T - v)
#   P2T    = ordered pairs of non-empty disjoint paths covering T
#   P2B    = the same for T - v
import itertools
import numpy as np
from numba import njit


def perm_table(n):
    return np.array(list(itertools.permutations(range(n))), dtype=np.int64)


@njit(cache=True)
def pidx(i, j):
    if i < j:
        return j * (j - 1) // 2 + i
    return i * (i - 1) // 2 + j


@njit(cache=True)
def out_masks(code, n, out):
    for i in range(n):
        out[i] = 0
    for j in range(1, n):
        base = j * (j - 1) // 2
        for i in range(j):
            if (code >> (base + i)) & 1:
                out[i] |= (1 << j)
            else:
                out[j] |= (1 << i)


@njit(cache=True)
def code_from_out(n, out):
    code = 0
    for j in range(1, n):
        base = j * (j - 1) // 2
        for i in range(j):
            if (out[i] >> j) & 1:
                code |= (1 << (base + i))
    return code


@njit(cache=True)
def class_reps(n, perms):
    # orbit marking: the least code of every isomorphism class of tournaments of order n (n <= 8)
    npairs = n * (n - 1) // 2
    total = 1 << npairs
    seen = np.zeros((total >> 3) + 1, np.uint8)
    nperm = perms.shape[0]
    tgt = np.zeros((nperm, max(npairs, 1)), np.int64)
    flip = np.zeros((nperm, max(npairs, 1)), np.int64)
    for p in range(nperm):
        for j in range(1, n):
            for i in range(j):
                k = j * (j - 1) // 2 + i
                a = perms[p, i]
                b = perms[p, j]
                tgt[p, k] = pidx(a, b)
                flip[p, k] = 0 if a < b else 1
    reps = np.zeros(8192, np.int64)
    nr = 0
    for code in range(total):
        if (seen[code >> 3] >> (code & 7)) & 1:
            continue
        reps[nr] = code
        nr += 1
        for p in range(nperm):
            c2 = 0
            for k in range(npairs):
                bit = ((code >> k) & 1) ^ flip[p, k]
                c2 |= (bit << tgt[p, k])
            seen[c2 >> 3] |= (1 << (c2 & 7))
    return reps[:nr].copy()


@njit(cache=True)
def dp_tables(n, out, f, g, hS):
    full = (1 << n) - 1
    hS[0] = 1
    for S in range(1, full + 1):
        tot = 0
        for v in range(n):
            if ((S >> v) & 1) == 0:
                f[S, v] = 0
                g[S, v] = 0
                continue
            R = S ^ (1 << v)
            if R == 0:
                f[S, v] = 1
                g[S, v] = 1
                tot += 1
                continue
            a = 0
            b = 0
            for w in range(n):
                if (R >> w) & 1:
                    if (out[w] >> v) & 1:
                        a += f[R, w]
                    else:
                        b += g[R, w]
            f[S, v] = a
            g[S, v] = b
            tot += a
        hS[S] = tot


@njit(cache=True)
def rooted_direct(n, f, g, hS, res):
    # res[v] = (H, s, e, Yplus, Yminus, d, P2T, P2B) for every root v, straight from the definitions
    full = (1 << n) - 1
    H = hS[full]
    P2T = 0
    for S in range(1, full):
        P2T += hS[S] * hS[full ^ S]
    for v in range(n):
        yp = 0
        ym = 0
        for S in range(1, full + 1):
            if (S >> v) & 1:
                yp += f[S, v] * hS[full ^ S]
                ym += g[S, v] * hS[full ^ S]
        rest = full ^ (1 << v)
        p2b = 0
        S = (rest - 1) & rest
        while S > 0:
            p2b += hS[S] * hS[rest ^ S]
            S = (S - 1) & rest
        res[v, 0] = H
        res[v, 1] = g[full, v]
        res[v, 2] = f[full, v]
        res[v, 3] = yp
        res[v, 4] = ym
        res[v, 5] = hS[rest]
        res[v, 6] = P2T
        res[v, 7] = p2b


@njit(cache=True)
def base_data(n, f, g, hS, s, e, yp, ym, Q):
    # start / end vectors, Yplus / Yminus vectors and the exposed-slot matrix Q of THM-4097 (5)
    full = (1 << n) - 1
    for v in range(n):
        s[v] = g[full, v]
        e[v] = f[full, v]
        a = 0
        b = 0
        for S in range(1, full + 1):
            if (S >> v) & 1:
                a += f[S, v] * hS[full ^ S]
                b += g[S, v] * hS[full ^ S]
        yp[v] = a
        ym[v] = b
    for a in range(n):
        for b in range(n):
            Q[a, b] = 0
            if a == b:
                continue
            rest = full ^ (1 << a) ^ (1 << b)
            S0 = rest
            while True:
                A = S0 | (1 << a)
                Q[a, b] += f[A, a] * g[full ^ A, b]
                if S0 == 0:
                    break
                S0 = (S0 - 1) & rest


@njit(cache=True)
def end_components(n, out, comp):
    # comp[0] = vertex mask of the initial strong component, comp[1] = mask of the terminal one
    full = (1 << n) - 1
    first = 0
    last = 0
    for v in range(n):
        R = 1 << v
        while True:
            N = R
            for u in range(n):
                if (R >> u) & 1:
                    N |= out[u]
            if N == R:
                break
            R = N
        if R == full:
            first |= (1 << v)
        R = 1 << v
        while True:
            N = R
            for u in range(n):
                if ((R >> u) & 1) == 0 and (out[u] & R) != 0:
                    N |= (1 << u)
            if N == R:
                break
            R = N
        if R == full:
            last |= (1 << v)
    comp[0] = first
    comp[1] = last


@njit(cache=True)
def cut_tuple(n, S, H, s, e, yp, ym, Q, comp, t):
    # rooted invariants of (B + w, w) where S is the out-set of the new vertex w:
    # t = (x, s, e, Yplus, Yminus, d, strong, out-degree of the root)
    x = 0
    ss = 0
    ee = 0
    a = H
    b = H
    deg = 0
    for u in range(n):
        if (S >> u) & 1:
            ss += s[u]
            b += ym[u]
            deg += 1
        else:
            ee += e[u]
            a += yp[u]
            for v in range(n):
                if (S >> v) & 1:
                    x += Q[u, v]
    t[0] = x + ss + ee
    t[1] = ss
    t[2] = ee
    t[3] = a
    t[4] = b
    t[5] = H
    full = (1 << n) - 1
    # B + w is strong iff w beats a vertex of the initial component and loses to one of the terminal component
    if (S & comp[0]) != 0 and (((~S) & full) & comp[1]) != 0:
        t[6] = 1
    else:
        t[6] = 0
    t[7] = deg


@njit(cache=True)
def rooted_rows_from_bases(n, codes, rows):
    # rows[(k * 2^n + S)] = cut_tuple for base codes[k] and out-set S; returns number of rows
    nb = codes.shape[0]
    size = 1 << n
    f = np.zeros((size, n), np.int64)
    g = np.zeros((size, n), np.int64)
    hS = np.zeros(size, np.int64)
    out = np.zeros(n, np.int64)
    s = np.zeros(n, np.int64)
    e = np.zeros(n, np.int64)
    yp = np.zeros(n, np.int64)
    ym = np.zeros(n, np.int64)
    Q = np.zeros((n, n), np.int64)
    comp = np.zeros(2, np.int64)
    t = np.zeros(8, np.int64)
    r = 0
    for k in range(nb):
        out_masks(codes[k], n, out)
        dp_tables(n, out, f, g, hS)
        base_data(n, f, g, hS, s, e, yp, ym, Q)
        end_components(n, out, comp)
        H = hS[size - 1]
        for S in range(size):
            cut_tuple(n, S, H, s, e, yp, ym, Q, comp, t)
            for c in range(8):
                rows[r, c] = t[c]
            r += 1
    return r


@njit(cache=True)
def rooted_rows_direct(n, codes, rows):
    # rows = (H, s, e, Yplus, Yminus, d, strong, out-degree, P2T, P2B) for every root of every tournament in codes
    nb = codes.shape[0]
    size = 1 << n
    f = np.zeros((size, n), np.int64)
    g = np.zeros((size, n), np.int64)
    hS = np.zeros(size, np.int64)
    out = np.zeros(n, np.int64)
    res = np.zeros((n, 8), np.int64)
    comp = np.zeros(2, np.int64)
    r = 0
    for k in range(nb):
        out_masks(codes[k], n, out)
        dp_tables(n, out, f, g, hS)
        rooted_direct(n, f, g, hS, res)
        end_components(n, out, comp)
        strong = 1 if comp[0] == size - 1 else 0
        for v in range(n):
            for c in range(6):
                rows[r, c] = res[v, c]
            rows[r, 6] = strong
            deg = 0
            m = out[v]
            while m:
                deg += 1
                m &= m - 1
            rows[r, 7] = deg
            rows[r, 8] = res[v, 6]
            rows[r, 9] = res[v, 7]
            r += 1
    return r


@njit(cache=True)
def count_h(n, code):
    size = 1 << n
    f = np.zeros((size, n), np.int64)
    g = np.zeros((size, n), np.int64)
    hS = np.zeros(size, np.int64)
    out = np.zeros(n, np.int64)
    out_masks(code, n, out)
    dp_tables(n, out, f, g, hS)
    return hS[size - 1]


def extend_code(code, n, S):
    # code of B + w (w = vertex n) where S is the out-set of w
    full = (1 << n) - 1
    return code | (((~S) & full) << (n * (n - 1) // 2))
