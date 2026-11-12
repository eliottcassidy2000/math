"""b3_independent_closure.py - INDEPENDENT closure of the one-apex chain semigroup (audit of C4 and of
scripts/chain_semigroup.py: closure).

No normal form: every binary word w of length 0..NMAX is a block, with its own data
    g(w) = number of Hamiltonian paths of the one-apex tournament (closed formula below, checked against the
           subset DP of tt.py for all words of length <= 8 and against the recursion for all words used),
    P(w) + 1 = W+,  Q(w) + 1 = W-,  n = length,  tau = 2^n,
and every chain length i >= 1.  Transition (Proposition B1):  A' = g A + 2^i (Q+1) W,  W' = (P+1) A + 2^(i+n) W.
All states (A, W) with A <= X are generated; a state is expanded iff A + 2 W <= X (otherwise no successor has
A' <= X, because A' >= A + 2 W).  The loop over i stops when A + 2^i W > X (same reason).

Word length: a word with A' <= X has Q + 1 <= X/2, so its letters at positions >= K - 1 are all 1 (X = 2^K) and it
has the same (g, Q) as its prefix of length K - 1; the first block needs even cores P <= 2X - 2, i.e. length K + 1.
NMAX = K + 1 is used for all blocks (more than needed).

Closed formula for g (auditor's): a Hamiltonian path is (increasing set S_1 of spine vertices, apex, increasing
set S_2), S_1 + S_2 = spine, last(S_1) -> apex (letter 0), apex -> first(S_2) (letter 1):
    g(w) = [n = 0 or w_0 = 1] + [n >= 1 and w_(n-1) = 0] - [n = 0] ... (see code) + #{i: w_i = 0, w_(i+1) = 1}
           + sum over j < i with w_j = 1, w_i = 0 of 2^(i-j-1).

usage: python b3_independent_closure.py 12          (no cut at all)
       python b3_independent_closure.py 14 cut      (words sorted by Q+1; only the prefix with 2^i W (Q+1) <= X - A
                                                     is scanned: an exact cut, not a normal form)
"""
import os
import sys
from itertools import product

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def g_closed(w):
    n = len(w)
    if n == 0:
        return 1
    tot = (1 if w[0] == 1 else 0) + (1 if w[-1] == 0 else 0)
    for i in range(n - 1):
        if w[i] == 0 and w[i + 1] == 1:
            tot += 1
    for j in range(n):
        if w[j] == 1:
            for i in range(j + 1, n):
                if w[i] == 0:
                    tot += 1 << (i - j - 1)
    return tot


def word_arrays(nmax):
    g, P1, Q1, nn = [], [], [], []
    for n in range(0, nmax + 1):
        for w in product((0, 1), repeat=n):
            P = 0
            for a in w:
                P = 2 * P + a
            Q = sum(1 << j for j, a in enumerate(w) if not a)
            g.append(g_closed(w))
            P1.append(P + 1)
            Q1.append(Q + 1)
            nn.append(n)
    return (np.array(g, dtype=np.int64), np.array(P1, dtype=np.int64), np.array(Q1, dtype=np.int64),
            np.array(nn, dtype=np.int64))


def self_test():
    from tt import one_apex, H, g_fun
    bad = 0
    for n in range(0, 9):
        for w in product((0, 1), repeat=n):
            m, adj = one_apex(w)
            bad += (H(m, adj) != g_closed(w))
    bad2 = 0
    for n in range(0, 14):
        for w in product((0, 1), repeat=n):
            P = 0
            for a in w:
                P = 2 * P + a
            bad2 += (g_fun(P) != g_closed(w))
    print("self-test: closed formula for g against the subset DP (all words of length <= 8): %d failures; against the recursion (length <= 13): %d failures" % (bad, bad2))


def closure(K, cut):
    X = 1 << K
    nmax = K + 1
    g, P1, Q1, nn = word_arrays(nmax)
    if cut:
        order = np.argsort(Q1, kind="stable")
        g, P1, Q1, nn = g[order], P1[order], Q1[order], nn[order]
    pow_n = (np.int64(1) << nn)
    hit = np.zeros(X + 1, dtype=bool)
    seen = set()
    stack = []
    # first block: all words
    sel = g <= X
    hit[g[sel]] = True
    for a, w in set(zip(g[sel].tolist(), P1[sel].tolist())):
        if a + 2 * w <= X:
            seen.add((a, w))
            stack.append((a, w))
    expanded = 0
    transitions = 0
    while stack:
        A, W = stack.pop()
        expanded += 1
        t = 2 * W
        while A + t <= X:
            if cut:
                m = int(np.searchsorted(Q1, (X - A) // t, side="right"))
                vals = g[:m] * A + t * Q1[:m]
                ok = vals <= X
                A2 = vals[ok]
                W2 = P1[:m][ok] * A + t * pow_n[:m][ok]
            else:
                vals = g * A + t * Q1
                ok = vals <= X
                A2 = vals[ok]
                W2 = P1[ok] * A + t * pow_n[ok]
            transitions += int(A2.shape[0])
            hit[A2] = True
            live = (A2 + 2 * W2) <= X
            for a, w in set(zip(A2[live].tolist(), W2[live].tolist())):
                if (a, w) not in seen:
                    seen.add((a, w))
                    stack.append((a, w))
            t *= 2
    odd = np.arange(1, X + 1, 2)
    miss = odd[~hit[odd]].tolist()
    return miss, expanded, transitions, int(g.shape[0])


def author_closure(K):
    """the author's closure (scripts/chain_semigroup.py), full list of misses"""
    sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "scripts")))
    import chain_semigroup as cs
    X = 1 << K
    half = X // 2
    G = np.array(cs.g_list(2 * X + 8), dtype=np.int64)
    Lb = np.array([cs.L_of(p) for p in range(2 * X + 9)], dtype=np.int64)
    Rb = np.array([cs.R_of(p) for p in range(2 * X + 9)], dtype=np.int64)
    hit = np.zeros(X + 1, dtype=np.uint8)
    width = half + 1
    seen = np.zeros((X // 2 + 1) * width, dtype=np.uint8)
    cap = 1 << 24
    stackA = np.zeros(cap, dtype=np.int64)
    stackW = np.zeros(cap, dtype=np.int64)
    ns = 0
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
    states = cs.closure(X, G, Lb, Rb, hit, seen, stackA, stackW, ns)
    odd = np.arange(1, X + 1, 2)
    return odd[hit[odd] == 0].tolist(), int(states)


def stored_first80(K):
    path = os.path.join(HERE, "..", "runs", "closure_%d.txt" % K)
    with open(path) as f:
        lines = f.read().splitlines()
    count = int(lines[0].split("odd numbers missed:")[1])
    lst = [int(s) for s in lines[1].split("[")[1].split("]")[0].split(",")]
    return count, lst


def main():
    K = int(sys.argv[1])
    cut = len(sys.argv) > 2 and sys.argv[2] == "cut"
    self_test()
    miss, expanded, transitions, nwords = closure(K, cut)
    print("independent closure below 2^%d (%s; %d words of length <= %d as blocks): %d states expanded, %d transitions with A' <= X; odd numbers missed: %d"
          % (K, "exact prefix cut" if cut else "no cut", nwords, K + 1, expanded, transitions, len(miss)))
    print("   first 30 misses:", miss[:30])
    print("   last 10 misses:", miss[-10:])
    ref, states = author_closure(K)
    print("author's closure below 2^%d: %d live states, %d misses; full lists identical: %s" % (K, states, len(ref), ref == miss))
    if os.path.exists(os.path.join(HERE, "..", "runs", "closure_%d.txt" % K)):
        cnt, first = stored_first80(K)
        print("stored runs/closure_%d.txt: count %d (equal: %s); its first %d misses equal mine: %s" % (K, cnt, cnt == len(miss), len(first), first == miss[:len(first)]))
    if K >= 14:
        sub = [m for m in miss if m <= 4096]
        ref12, _ = author_closure(12)
        print("restriction of the misses to [1, 2^12]: %d; equal to the closure below 2^12: %s" % (len(sub), sub == ref12))
    with open(os.path.join(HERE, "b3_misses_closure_%d.txt" % K), "w") as f:
        f.write("misses of the full closure of the one-apex chain semigroup below 2^%d (auditor's independent computation): %d\n" % (K, len(miss)))
        f.write(" ".join(str(m) for m in miss) + "\n")


if __name__ == "__main__":
    main()
