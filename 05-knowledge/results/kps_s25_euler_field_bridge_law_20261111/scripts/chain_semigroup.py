"""chain_semigroup.py - the one-apex chain semigroup (kind-pasteur-2026-11-11-S25; pre-registered C1-C5).

For p >= 0: L(p) = bit length, R(p) = reversal of the complement of p on L(p) bits, g = the one-apex function
(g(0) = 1, g(2v+1) = g(v), g(2v) = g(v) + v + (v mod 2)).  For i >= 1:

    M(p, i) = [[ g(p), 2^i (R(p) + 1) ], [ p + 1, 2^(i + L(p)) ]].

Claim: the chain Y_0 => C_1 => Y_1 => ... => C_k => Y_k of one-apex tournaments (apex words of values p_j), with
transitive chains C_j of i_j >= 1 vertices between them and the arcs between consecutive apexes reversed, has
H = top-left entry of M(p_k, i_k) ... M(p_1, i_1) M(p_0, 1).

usage: python chain_semigroup.py c1            brute force (dynamic programming over subsets)
       python chain_semigroup.py depth K       coverage of depth <= 2 and <= 3 below 2^K, tuple count
       python chain_semigroup.py closure K     full closure below 2^K
"""
import sys
from itertools import product

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


def L_of(p):
    return p.bit_length()


def R_of(p):
    L = p.bit_length()
    r = 0
    for j in range(L):
        if not (p >> j) & 1:                 # complement bit at position j
            r |= 1 << (L - 1 - j)
    return r


# ---------------------------------------------------------------- C1: brute force

def ham(n, adj):
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1, 1 << n):
        row = dp[mask]
        for last in range(n):
            c = row[last]
            if c:
                nxt = adj[last] & ~mask
                while nxt:
                    b = nxt & -nxt
                    dp[mask | b][b.bit_length() - 1] += c
                    nxt ^= b
    return sum(dp[full])


def build_chain(words, chains):
    """blocks Y_0, C_1, Y_1, ..., C_k, Y_k in this order (ordinal sum), apex arcs reversed.
    words[j] = tuple of letters (1 = apex beats the spine vertex); chains[j-1] = length of C_j."""
    blocks = []                      # list of vertex lists; apex is the first vertex of a Y-block
    n = 0
    apexes = []
    kinds = []
    for j, w in enumerate(words):
        if j > 0:
            blocks.append(list(range(n, n + chains[j - 1])))
            kinds.append(("C", None))
            n += chains[j - 1]
        blk = list(range(n, n + 1 + len(w)))
        blocks.append(blk)
        kinds.append(("Y", w))
        apexes.append(blk[0])
        n += 1 + len(w)
    adj = [0] * n

    def arc(a, b):
        adj[a] |= 1 << b

    for bi, blk in enumerate(blocks):
        kind, w = kinds[bi]
        if kind == "C":
            for x in range(len(blk)):
                for y in range(x + 1, len(blk)):
                    arc(blk[x], blk[y])
        else:
            apex, spine = blk[0], blk[1:]
            for x in range(len(spine)):
                for y in range(x + 1, len(spine)):
                    arc(spine[x], spine[y])
                if w[x]:
                    arc(apex, spine[x])
                else:
                    arc(spine[x], apex)
        for bj in range(bi + 1, len(blocks)):
            for a in blk:
                for b in blocks[bj]:
                    arc(a, b)
    for j in range(1, len(apexes)):         # reverse the arc between consecutive apexes
        a, b = apexes[j - 1], apexes[j]
        adj[a] &= ~(1 << b)
        adj[b] |= 1 << a
    return n, adj


def word_data(w, G):
    n = len(w)
    P = 0
    for a in w:
        P = 2 * P + a
    Q = sum(1 << j for j, a in enumerate(w) if not a)
    return G[P], P, Q, n


def c1():
    G = g_list(1 << 6)
    tot = bad = 0
    for depth in (2, 3):
        for lens in product(range(0, 3), repeat=depth):
            for words in product(*[list(product((0, 1), repeat=l)) for l in lens]):
                for chains in product((1, 2), repeat=depth - 1):
                    n, adj = build_chain(words, chains)
                    if n > 13:
                        continue
                    h = ham(n, adj)
                    g0, P0, Q0, n0 = word_data(words[0], G)
                    A, W = g0, P0 + 1
                    for j in range(1, depth):
                        gj, Pj, Qj, nj = word_data(words[j], G)
                        i = chains[j - 1]
                        A, W = gj * A + (1 << i) * (Qj + 1) * W, (Pj + 1) * A + (1 << (i + nj)) * W
                    tot += 1
                    bad += (A != h)
    print("C1  chains of depth 2 and 3, words of length <= 2, chain lengths 1, 2: %d cases, %d mismatches" % (tot, bad))
    # normal form: leading B's are absorbed in the exponent
    bad2 = 0
    for w in product((0, 1), repeat=4):
        g0, P, Q, n = word_data(w, G)
        k = 0
        while k < n and w[k] == 0:
            k += 1
        p = P
        if (Q + 1) != (1 << k) * (R_of(p) + 1) or n != k + L_of(p):
            bad2 += 1
    print("    normal form M(B^k w; i) = M(w; i + k) on all words of length 4: %d failures" % bad2)


# ---------------------------------------------------------------- C2, C3, C5: depth <= 3

@njit(cache=True)
def depth_cover(X, G, Lb, Rb, hit2, hit3):
    cnt2 = 0
    cnt3 = 0
    # depth 2: N = g(q) g(p) + 2^i (q + 1)(p + 1)
    for p in range(0, X):
        if 2 * (p + 1) > X:
            break
        for q in range(0, X):
            if 2 * (q + 1) * (p + 1) > X:
                break
            base = G[q] * G[p]
            t = 2 * (q + 1) * (p + 1)
            while base + t <= X:
                hit2[base + t] = 1
                cnt2 += 1
                t *= 2
    # depth 3
    for p1 in range(0, X):
        L1 = Lb[p1]
        if (1 << (L1 + 2)) > X:
            break
        g1 = G[p1]
        r1 = Rb[p1] + 1
        i1 = 1
        while (1 << (i1 + L1 + 1)) <= X:
            e1 = 1 << i1
            f1 = 1 << (i1 + L1)
            for p0 in range(0, X):
                if 2 * f1 * (p0 + 1) > X:
                    break
                g0 = G[p0]
                A1 = g1 * g0 + e1 * r1 * (p0 + 1)
                W1 = (p1 + 1) * g0 + f1 * (p0 + 1)
                if A1 + 2 * W1 > X:
                    continue
                t = 2 * W1
                while t <= X:
                    for q in range(0, X):
                        if (q + 1) * t > X:
                            break
                        N = G[q] * A1 + (q + 1) * t
                        if N <= X:
                            hit3[N] = 1
                            cnt3 += 1
                    t *= 2
            i1 += 1
    return cnt2, cnt3


def depth(K):
    X = 1 << K
    G = np.array(g_list(X), dtype=np.int64)
    Lb = np.array([L_of(p) for p in range(X + 1)], dtype=np.int64)
    Rb = np.array([R_of(p) for p in range(X + 1)], dtype=np.int64)
    hit1 = np.zeros(X + 1, dtype=np.uint8)
    hit1[G[G <= X]] = 1
    hit2 = np.zeros(X + 1, dtype=np.uint8)
    hit3 = np.zeros(X + 1, dtype=np.uint8)
    cnt2, cnt3 = depth_cover(X, G, Lb, Rb, hit2, hit3)
    h12 = hit1 | hit2
    h123 = h12 | hit3
    odd = np.arange(1, X + 1, 2)
    m2 = odd[h12[odd] == 0]
    m3 = odd[h123[odd] == 0]
    lo, hi = X // 2, X
    f2 = np.sum((m2 >= lo) & (m2 < hi)) / ((hi - lo) / 2)
    f3 = np.sum((m3 >= lo) & (m3 < hi)) / ((hi - lo) / 2)
    print("depth <= 2 below 2^%d: %d odd numbers missed; fraction missed in the top dyadic block: %.4f" % (K, len(m2), f2))
    print("depth <= 3 below 2^%d: %d odd numbers missed; fraction missed in the top dyadic block: %.5f" % (K, len(m3), f3))
    print("   first misses at depth <= 3:", m3[:40].tolist())
    print("   last misses at depth <= 3:", m3[-12:].tolist())
    print("tuples: depth 2: %d, depth 3: %d; mean multiplicity per odd number below 2^%d: depth 2: %.1f, depth <= 3: %.1f"
          % (cnt2, cnt3, K, cnt2 / (X / 2), (cnt2 + cnt3) / (X / 2)))
    k = 6
    while (1 << (k + 1)) <= X:
        a, b = 1 << k, 1 << (k + 1)
        print("   [2^%d, 2^%d): depth<=2 misses %d, depth<=3 misses %d of %d" % (k, k + 1, np.sum((m2 >= a) & (m2 < b)), np.sum((m3 >= a) & (m3 < b)), (b - a) // 2))
        k += 1


# ---------------------------------------------------------------- C4: full closure

@njit(cache=True)
def closure(X, G, Lb, Rb, hit, seen, stackA, stackW, nstack):
    half = X // 2
    width = half + 1
    states = 0
    while nstack > 0:
        nstack -= 1
        A = stackA[nstack]
        W = stackW[nstack]
        states += 1
        p = 0
        while True:
            g = G[p]
            if g * A + 2 * W > X:
                if p > 0 and p > 2 * X // A + 4:
                    break
                p = 2 if p == 0 else p + 2
                if p > X:
                    break
                continue
            r = Rb[p] + 1
            L = Lb[p]
            e = 2
            f = 1 << (L + 1)
            while g * A + e * r * W <= X:
                A2 = g * A + e * r * W
                hit[A2] = 1
                W2 = (p + 1) * A + f * W
                while W2 <= half:                      # (A2, 2^t W2): trailing A's of the word
                    idx = ((A2 - 1) // 2) * width + W2
                    if seen[idx] == 0:
                        seen[idx] = 1
                        stackA[nstack] = A2
                        stackW[nstack] = W2
                        nstack += 1
                    W2 *= 2
                e *= 2
                f *= 2
            p = 2 if p == 0 else p + 2
            if p > X:
                break
    return states


def run_closure(K):
    X = 1 << K
    half = X // 2
    G = np.array(g_list(2 * X + 8), dtype=np.int64)
    Lb = np.array([L_of(p) for p in range(2 * X + 9)], dtype=np.int64)
    Rb = np.array([R_of(p) for p in range(2 * X + 9)], dtype=np.int64)
    hit = np.zeros(X + 1, dtype=np.uint8)
    width = half + 1
    seen = np.zeros((X // 2 + 1) * width, dtype=np.uint8)
    cap = 1 << 26
    stackA = np.zeros(cap, dtype=np.int64)
    stackW = np.zeros(cap, dtype=np.int64)
    ns = 0
    # initial states (g(p0), p0 + 1), all p0 (the first block); one-apex numbers are hits
    for p0 in range(0, 2 * X + 8):
        a, w = int(G[p0]), p0 + 1
        if a <= X:
            hit[a] = 1
            if w <= half:
                idx = ((a - 1) // 2) * width + w
                if not seen[idx]:
                    seen[idx] = 1
                    stackA[ns] = a
                    stackW[ns] = w
                    ns += 1
    states = closure(X, G, Lb, Rb, hit, seen, stackA, stackW, ns)
    odd = np.arange(1, X + 1, 2)
    miss = odd[hit[odd] == 0]
    print("full closure below 2^%d: %d live states; odd numbers missed: %d" % (K, states, len(miss)))
    print("   misses:", miss[:80].tolist())


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "c1":
        c1()
    elif mode == "depth":
        depth(int(sys.argv[2]))
    elif mode == "closure":
        run_closure(int(sys.argv[2]))
