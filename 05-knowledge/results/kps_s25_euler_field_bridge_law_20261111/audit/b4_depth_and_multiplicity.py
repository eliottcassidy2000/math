"""b4_depth_and_multiplicity.py - audit of C2, C3, C5 (scripts/chain_semigroup.py: depth).

Auditor's own enumeration of the chains of depth 1, 2, 3 in FULLY STRIPPED form:

    blocks  p_0, p_1, p_2 in E = {0} + {even p >= 2}     (no leading letter 0, no trailing letter 1),
    chains  c_1, c_2 >= 1  (all vertices that can be moved into a chain are in the chain).

Every chain Y_0 => C_1 => Y_1 (=> C_2 => Y_2) of the note is one of these, up to sources and sinks of the whole
tournament (leading letters 0 of Y_0, trailing letters 1 of the last block), which do not change H:
  - leading letters 0 of Y_j (j >= 1) are chain vertices               (README, Corollary B2)
  - trailing letters 1 of Y_j (j < last) are chain vertices as well    (audit: M(p', i') M(2p+1, i) = M(p', i'+1) M(p, i))
So a stripped chain with chain lengths (c_1) resp. (c_1, c_2) is counted
  - c_1 (c_1 + 1)/2 times resp. c_1 * c_2 (c_2 + 1)/2 times by depth_cover() of the author
        (first block: any p, t trailing 1's; middle block: value p_1, any trailing 1's, exponent i_1 >= 1;
         last block: all q, i.e. k leading 0's, exponent i >= 1)
  - c_1 resp. c_1 c_2 times as a tuple of matrices M(p, i) with the last block stripped of trailing 1's
  - once as a tournament presentation.
This program prints the three counts, the hit sets, and compares with runs/depth_18.txt.

Also: completeness of the depth-1 hit set used by the author in depth(): G = g_list(X) only covers p <= X, while
even cores up to 2X - 2 can have g(p) <= X (g(p) >= p/2 + 1): values of g on (X, 2X) lie in the top dyadic block.

usage: python b4_depth_and_multiplicity.py 18
"""
import sys

import numpy as np
from numba import njit


def g_list(n):
    g = [1] * (n + 1)
    for v in range(1, n + 1):
        if v & 1:
            g[v] = g[v >> 1]
        else:
            h = v >> 1
            g[v] = g[h] + h + (h & 1)
    return g


def R_of(p):
    L = p.bit_length()
    comp = ((1 << L) - 1) ^ p
    return int(format(comp, "0%db" % L)[::-1], 2) if L else 0


@njit(cache=False)
def enumerate_stripped(X, G, Lb, Rb, hit2, hit3, cmax):
    # counts: [0] stripped presentations, [1] matrix tuples (last block stripped), [2] author's tuples
    c2 = np.zeros(3, dtype=np.int64)
    c3 = np.zeros(3, dtype=np.int64)
    hist2 = np.zeros(cmax, dtype=np.int64)          # stripped depth-2 chains by chain length c_1
    # depth 2: N = g(q) g(p0) + 2^c (q + 1)(p0 + 1), p0, q in E
    p0 = 0
    while 2 * (p0 + 1) <= X:
        q = 0
        while 2 * (q + 1) * (p0 + 1) <= X:
            base = G[q] * G[p0]
            t = 2 * (q + 1) * (p0 + 1)
            c = 1
            while base + t <= X:
                hit2[base + t] = 1
                c2[0] += 1
                c2[1] += c
                c2[2] += c * (c + 1) // 2
                hist2[c] += 1
                t *= 2
                c += 1
            q = 2 if q == 0 else q + 2
        p0 = 2 if p0 == 0 else p0 + 2
    # depth 3
    p1 = 0
    while (1 << (Lb[p1] + 2)) <= X:
        L1 = Lb[p1]
        g1 = G[p1]
        r1 = Rb[p1] + 1
        c1 = 1
        while (1 << (c1 + L1 + 1)) <= X:
            e1 = 1 << c1
            f1 = 1 << (c1 + L1)
            p0 = 0
            while 2 * f1 * (p0 + 1) <= X:
                g0 = G[p0]
                A1 = g1 * g0 + e1 * r1 * (p0 + 1)
                W1 = (p1 + 1) * g0 + f1 * (p0 + 1)
                if A1 + 2 * W1 <= X:
                    t = 2 * W1
                    cc = 1
                    while t <= X:
                        q = 0
                        while (q + 1) * t <= X:
                            N = G[q] * A1 + (q + 1) * t
                            if N <= X:
                                hit3[N] = 1
                                c3[0] += 1
                                c3[1] += c1 * cc
                                c3[2] += c1 * cc * (cc + 1) // 2
                            q = 2 if q == 0 else q + 2
                        t *= 2
                        cc += 1
                p0 = 2 if p0 == 0 else p0 + 2
            c1 += 1
        p1 = 2 if p1 == 0 else p1 + 2
    return c2, c3, hist2


def main():
    K = int(sys.argv[1])
    X = 1 << K
    G = np.array(g_list(2 * X + 2), dtype=np.int64)
    Lb = np.array([p.bit_length() for p in range(X + 2)], dtype=np.int64)
    Rb = np.array([R_of(p) for p in range(X + 2)], dtype=np.int64)
    # converse symmetry used by the q-parametrisation: g(R(p)) = g(p), R an involution on even numbers of L bits
    bad = sum(1 for p in range(2, 1 << 12, 2) if G[R_of(p)] != G[p] or R_of(R_of(p)) != p or R_of(p) % 2 or R_of(p).bit_length() != p.bit_length())
    print("converse symmetry on even p < 2^12 (g(R(p)) = g(p), R(R(p)) = p, R(p) even of the same length): %d failures" % bad)

    hit2 = np.zeros(X + 1, dtype=np.uint8)
    hit3 = np.zeros(X + 1, dtype=np.uint8)
    c2, c3, hist2 = enumerate_stripped(X, G, Lb, Rb, hit2, hit3, 64)

    # depth 1: complete (cores up to 2X) and as in the author's depth() (p <= X only)
    hit1_full = np.zeros(X + 1, dtype=np.uint8)
    vals = G[:2 * X + 1]
    hit1_full[vals[vals <= X]] = 1
    hit1_auth = np.zeros(X + 1, dtype=np.uint8)
    va = G[:X + 1]
    hit1_auth[va[va <= X]] = 1
    n1_cores = int(np.sum(G[0:2 * X + 1:2] <= X))            # stripped depth-1 presentations: p in E with g(p) <= X
    extra = np.nonzero((hit1_full == 1) & (hit1_auth == 0))[0]
    print("depth 1: one-apex values <= 2^%d: %d (author's table p <= X: %d); values of g attained only for p in (X, 2X]: %d, all in the top dyadic block: %s"
          % (K, int(hit1_full.sum()), int(hit1_auth.sum()), len(extra), bool(len(extra) == 0 or extra.min() > X // 2)))

    odd = np.arange(1, X + 1, 2)
    lo = X // 2
    for label, h1 in (("author's depth-1 table (p <= X)", hit1_auth), ("complete depth-1 table (p <= 2X)", hit1_full)):
        h12 = h1 | hit2
        h123 = h12 | hit3
        m2 = odd[h12[odd] == 0]
        m3 = odd[h123[odd] == 0]
        t2 = int(np.sum(m2 >= lo))
        t3 = int(np.sum(m3 >= lo))
        print("%s: depth <= 2 misses %d below 2^%d, %d in the top block (%.4f); depth <= 3 misses %d, %d in the top block (%.5f)"
              % (label, len(m2), K, t2, t2 / (X / 4), len(m3), t3, t3 / (X / 4)))
        if K >= 14:
            print("    depth <= 3 misses below 2^14: %d" % int(np.sum(m3 < (1 << 14))))
    lost2 = [int(v) for v in extra if hit2[v] == 0]
    lost3 = [int(v) for v in extra if hit2[v] == 0 and hit3[v] == 0]
    print("one-apex values missing from the author's depth-1 table and not reached at depth 2: %d; not reached at depth 2 or 3: %d %s"
          % (len(lost2), len(lost3), lost3[:20]))

    half_odds = X / 2
    print("counts of representations with N <= 2^%d:" % K)
    print("  depth 2: stripped presentations %d; matrix tuples (last block stripped) %d; author's tuples %d" % (c2[0], c2[1], c2[2]))
    print("  depth 3: stripped presentations %d; matrix tuples (last block stripped) %d; author's tuples %d" % (c3[0], c3[1], c3[2]))
    print("  depth 1: stripped presentations (cores p in E with g(p) <= X): %d" % n1_cores)
    for name, idx in (("stripped presentations", 0), ("matrix tuples, last block stripped", 1), ("author's tuples", 2)):
        tot = int(c2[idx] + c3[idx])
        print("  mean multiplicity per odd number, depth 2 + depth 3, %-36s: %.1f   (depth 2 alone: %.1f)   threshold of C5: 40"
              % (name, tot / half_odds, c2[idx] / half_odds))
    print("  stripped depth-2 chains by chain length c_1 = 1..12:", hist2[1:13].tolist())


if __name__ == "__main__":
    main()
