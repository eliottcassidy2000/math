"""direct_joins.py - direct joins of one-apex blocks rooted at their apexes, without a bank
(kind-pasteur-2026-11-11-S25; pre-registered addendum D).

State of a rooted tournament: (A, o, e), A = H, o (e) = sum of H(X - P) over the paths P from the root with an odd
(even) number of vertices.  Joining a one-apex block (apex word of length n, P = value of the letters 1 with the
first spine vertex as top bit, Q = value of the letters 0 with the first spine vertex as unit bit), the arc between
the old root and the new apex reversed:

    A' = A g(P) + 2 (o E + e O),   o' = A (1 + P//2) + o t_o + e t_e,   e' = A ((P+1)//2) + o t_e + e t_o,
    O = 1 + Q//2, E = (Q+1)//2, t_e = t_o = 2^(n-1) (n >= 1), t_e = 1, t_o = 0 (n = 0).

A chain of i >= 1 vertices before the block is the same as i leading letters 0 of its word.

usage: python direct_joins.py d1          brute force
       python direct_joins.py closure K   full closure below 2^K
"""
import sys
from itertools import product

import numpy as np
from numba import njit

from chain_semigroup import g_list, ham, build_chain, word_data, L_of, R_of


def d1():
    G = g_list(1 << 6)
    tot = bad = 0
    for depth in (2, 3):
        for lens in product(range(0, 3), repeat=depth):
            for words in product(*[list(product((0, 1), repeat=l)) for l in lens]):
                for chains in product((0, 1), repeat=depth - 1):
                    n, adj = build_chain(words, chains)
                    if n > 13:
                        continue
                    h = ham(n, adj)
                    g0, P0, Q0, n0 = word_data(words[0], G)
                    A, o, e = g0, 1 + P0 // 2, (P0 + 1) // 2
                    for j in range(1, depth):
                        gj, Pj, Qj, nj = word_data(words[j], G)
                        i = chains[j - 1]
                        if i >= 1:
                            w = (o + e) << (i - 1)
                            o, e = w, w
                        O, E = 1 + Qj // 2, (Qj + 1) // 2
                        te, to = (1 << (nj - 1), 1 << (nj - 1)) if nj >= 1 else (1, 0)
                        A, o, e = (A * gj + 2 * (o * E + e * O),
                                   A * (1 + Pj // 2) + o * to + e * te,
                                   A * ((Pj + 1) // 2) + o * te + e * to)
                    tot += 1
                    bad += (A != h)
    print("D1  direct joins, depth 2 and 3, words of length <= 2, chain lengths 0, 1: %d cases, %d mismatches" % (tot, bad))


@njit(cache=True)
def insert(table, key):
    mask = table.shape[0] - 1
    h = (key * 0x9E3779B97F4A7C15) & 0x7FFFFFFFFFFFFFFF
    idx = h & mask
    while True:
        cur = table[idx]
        if cur == key:
            return False
        if cur == -1:
            table[idx] = key
            return True
        idx = (idx + 1) & mask


@njit(cache=True)
def closure(X, coreP, coreG, coreQ, coreN, hit, table, stK, ns):
    half = X // 2
    ncore = coreP.shape[0]
    states = 0
    cap = stK.shape[0]
    while ns > 0:
        ns -= 1
        kk = stK[ns]
        A = kk >> 40
        o = (kk >> 20) & 0xFFFFF
        e = kk & 0xFFFFF
        states += 1
        for c in range(ncore):
            g = coreG[c]
            if g * A + 2 * e > X:
                if coreP[c] > 2 * (X // A) + 4:
                    break
                continue
            q1 = coreQ[c] + 1
            while True:                              # leading letters 0: Q + 1 = 2^i (Q' + 1)
                Q = q1 - 1
                O = 1 + Q // 2
                E = (Q + 1) // 2
                A2 = A * g + 2 * (o * E + e * O)
                if A2 > X:
                    break
                hit[A2] = 1
                p1 = coreP[c] + 1
                n = coreN[c]
                lead = 0
                tmp = q1 // (coreQ[c] + 1)
                while tmp > 1:
                    lead += 1
                    tmp //= 2
                n = n + lead
                while True:                          # trailing letters 1: P + 1 = 2^t (P' + 1)
                    P = p1 - 1
                    if n >= 1:
                        te = 1 << (n - 1)
                        to = te
                    else:
                        te = 1
                        to = 0
                    o2 = A * (1 + P // 2) + o * to + e * te
                    e2 = A * ((P + 1) // 2) + o * te + e * to
                    if e2 > half:
                        break
                    if o2 > half:
                        if A2 + 2 * e2 <= X:
                            hit[A2 + 2 * e2] = 1
                    else:
                        key = (A2 << 40) | (o2 << 20) | e2
                        if insert(table, key):
                            if ns >= cap:
                                return -1
                            stK[ns] = key
                            ns += 1
                    p1 *= 2
                    n += 1
                q1 *= 2
    return states


def run_closure(K):
    X = 1 << K
    half = X // 2
    G = g_list(4 * X + 8)
    cores = [0] + [p for p in range(2, 2 * X + 6, 2)]
    coreP = np.array(cores, dtype=np.int64)
    coreG = np.array([G[p] for p in cores], dtype=np.int64)
    coreQ = np.array([R_of(p) for p in cores], dtype=np.int64)
    coreN = np.array([L_of(p) for p in cores], dtype=np.int64)
    hit = np.zeros(X + 1, dtype=np.uint8)
    tbits = 25 if K <= 16 else 27
    table = np.full(1 << tbits, -1, dtype=np.int64)
    cap = 1 << (tbits - 1)
    stK = np.zeros(cap, dtype=np.int64)
    ns = 0
    for P in range(0, 4 * X + 8):
        A, o, e = G[P], 1 + P // 2, (P + 1) // 2
        if A <= X:
            hit[A] = 1
            if o <= half and e <= half:
                key = (A << 40) | (o << 20) | e
                if insert(table, key):
                    stK[ns] = key
                    ns += 1
    states = closure(X, coreP, coreG, coreQ, coreN, hit, table, stK, ns)
    odd = np.arange(1, X + 1, 2)
    miss = odd[hit[odd] == 0]
    print("direct joins, full closure below 2^%d: %d live states; odd numbers missed: %d" % (K, states, len(miss)))
    print("   misses:", miss[:120].tolist())


if __name__ == "__main__":
    if sys.argv[1] == "d1":
        d1()
    else:
        run_closure(int(sys.argv[2]))
