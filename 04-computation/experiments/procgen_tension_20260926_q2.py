#!/usr/bin/env python3
"""
procgen_tension_20260926_q2.py -- Q2: the maximum-density cycles of G_(sigma_k) (sigma_k = flip Bad_k).
They are exactly the periodic concatenations of first-descent blocks of density F_k (plus flip blocks '10'
when F_k = 1/2).  The upper Christoffel word of slope F_k is one such block, the unique balanced one; for
k >= 4 there are others, so 'every maximizer is Christoffel' is refuted.  Run through the runner.
"""
from fractions import Fraction
from math import comb, gcd

import numpy as np

from procgen_tension_20260926_lib import (
    as_cycle, below_c, best_lower, block_decomposition, canon, check, claim, cycle_word, density,
    is_fd_word, karp_single, necklace_count, periodic_point, potential, potential_np, potential_ok,
    residue_of, rho_all, sigma_k_mask, simple_cycles_adj, strictly_above_words, succ, tight_edges,
    upper_christoffel)


def primitive_root(w):
    """the primitive root of a word (the shortest u with w = u^m)"""
    L = len(w)
    for d in range(1, L + 1):
        if L % d == 0 and w[:d] * (L // d) == w:
            return w[:d]
    return w


def fd_blocks(F, k):
    """all first-descent words of density exactly F and length <= k (exact enumeration)"""
    a, d = F.numerator, F.denominator
    out = []
    m = 1
    while m * d <= k:
        for w in strictly_above_words(a, d, m):
            check(is_fd_word(w), "strictly-above word of density F is first-descent")
            out.append(w)
        m += 1
    return out


def balanced_rows(W, chunk=20000):
    """vectorized cyclic balance test for the rows of a 0/1 matrix (all rows of equal length L), in chunks"""
    n, L = W.shape
    ok = np.ones(n, dtype=bool)
    for st in range(0, n, chunk):
        Wc = W[st:st + chunk].astype(np.int16)
        WW = np.concatenate([Wc, Wc], axis=1)
        P = np.concatenate([np.zeros((len(Wc), 1), dtype=np.int16), np.cumsum(WW, axis=1, dtype=np.int16)], axis=1)
        okc = np.ones(len(Wc), dtype=bool)
        for m in range(1, L):
            win = P[:, m:m + L] - P[:, 0:L]
            okc &= (win.max(axis=1) - win.min(axis=1)) <= 1
        ok[st:st + chunk] = okc
    return ok


def realize_block_cycle(k, mask, word):
    """the Collatz periodic orbit of a first-descent word, read in G_(sigma_k): returns the node cycle"""
    x = periodic_point(word)          # Collatz signs
    nodes = []
    y = x
    for b in word:
        r = residue_of(y, k)
        check((r & 1) == b, "parity of the periodic orbit matches the word")
        nodes.append(r)
        y = (3 * y + 1) / 2 if b else y / 2
    check(y == x, "orbit closes")
    for i in range(len(nodes)):
        check(nodes[(i + 1) % len(nodes)] in succ(k, mask, nodes[i]), "consecutive residues form an edge of G_(sigma_k)")
    return nodes


def run():
    print("Q2. maximum-density cycles of G_(sigma_k)")
    # (1) rho_max(sigma_k) = F_k: potential at F_k on every edge (k <= 18) + the Christoffel block cycle; exact
    #     Karp cross-check for k <= 13
    rows = []
    for k in range(2, 19):
        m = sigma_k_mask(k)
        F = best_lower(k)
        psi = potential_np(k, m, F) if k > 10 else potential(k, m, F)
        check(psi is not None and potential_ok(k, m, F, psi), "potential at F_k for sigma_k, k=%d" % k)
        a, d = F.numerator, F.denominator
        cw = upper_christoffel(a, d)
        check(is_fd_word(cw), "upper Christoffel word of slope F_k is a first-descent word")
        nodes = realize_block_cycle(k, m, cw)
        check(density(nodes) == F, "Christoffel cycle has density F_k")
        if k <= 13:
            check(karp_single(k, m) == F, "exact Karp rho_max(sigma_k) = F_k, k=%d" % k)
        rows.append((k, str(F)))
    claim(True, "rho_max(sigma_k) = F_k for k = 2..18 (edge-checked integer potential at F_k plus the realized "
          "upper-Christoffel cycle; exact Karp agrees for k <= 13)")
    # (2) the block catalogue: first-descent words of density F_k and length <= k
    print("    k : F_k : #blocks by length : Christoffel block : #balanced among length-d blocks")
    cache = {}
    for k in range(2, 31):
        F = best_lower(k)
        a, d = F.numerator, F.denominator
        key = (F, k // d)                      # the block set depends only on F_k and floor(k/d)
        if key not in cache:
            cache.clear()
            cache[key] = fd_blocks(F, k)
        blocks = cache[key]
        bylen = {}
        for w in blocks:
            bylen[len(w)] = bylen.get(len(w), 0) + 1
        # cycle lemma: exactly one strictly-above rotation per necklace of (d, a)-words (gcd(a,d) = 1)
        check(gcd(a, d) == 1 and bylen.get(d, 0) == comb(d, a) // d == necklace_count(d, a),
              "cycle lemma count, k=%d" % k)
        cw = upper_christoffel(a, d)
        check(cw in blocks, "Christoffel block present")
        Wd = np.array([w for w in blocks if len(w) == d], dtype=np.int8)
        bal_mask = balanced_rows(Wd)
        bal = [list(map(int, Wd[i])) for i in np.nonzero(bal_mask)[0]]
        check(bal == [cw], "the upper Christoffel word is the unique balanced block of length d, k=%d" % k)
        # pointwise minimality of the Christoffel lattice path among length-d blocks, and (hence) its
        # minimal maximal excursion max_j 3^(a_j)/2^j
        pc = np.cumsum(np.array(cw, dtype=np.int16))
        check(bool(np.all(np.cumsum(Wd, axis=1, dtype=np.int16) >= pc[None, :])), "Christoffel path is lowest, k=%d" % k)
        del Wd
        if k in (2, 3, 4, 5, 7, 8, 16, 24, 26, 27, 30):
            print("    %2d : %s : %s : %s : %d" % (k, F, dict(sorted(bylen.items())), ''.join(map(str, cw)), len(bal)))
    cache.clear()
    claim(True, "k = 2..30: the density-F_k first-descent blocks of length d (F_k = a/d) are one per necklace "
          "(C(d,a)/d of them, cycle lemma); the upper Christoffel word is among them, is the unique balanced one, "
          "and has the pointwise lowest lattice path; longer blocks (length m d <= k, m >= 2) exist from k = 16 on")
    # (2') independent brute force: for k <= 18, the first-descent words of length <= k and density F_k, found by
    #      testing every word with the right number of ones, are exactly the strictly-above words
    from itertools import combinations
    for k in range(2, 19):
        F = best_lower(k)
        a, d = F.numerator, F.denominator
        brute = set()
        for L in range(1, k + 1):
            if (L * a) % d:
                continue
            ones = L * a // d
            for pos in combinations(range(L), ones):
                w = [0] * L
                for i in pos:
                    w[i] = 1
                if is_fd_word(w):
                    brute.add(tuple(w))
        check(brute == set(tuple(w) for w in fd_blocks(F, k)), "brute-force block set, k=%d" % k)
    claim(True, "k = 2..18: brute force over all words of each length confirms that the density-F_k first-descent words "
          "of length <= k are exactly the strictly-above words (Theorem M(ii))")
    # (3) every block of length <= k is realized as a cycle of G_(sigma_k) (k <= 14)
    for k in range(2, 15):
        m = sigma_k_mask(k)
        F = best_lower(k)
        for w in fd_blocks(F, k):
            nodes = realize_block_cycle(k, m, w)
            check(density(nodes) == F, "block cycle density")
    claim(True, "k = 2..14: every density-F_k first-descent block of length <= k is the word of a cycle of G_(sigma_k) "
          "(its Collatz periodic orbit, read mod 2^k, is a closed walk of G_(sigma_k))")
    # (4) exhaustive: all simple max-density cycles of G_(sigma_k), k <= 7, and their block decompositions
    print("    k : #max-density simple cycles : their canonical words (count)")
    for k in range(2, 8):
        m = sigma_k_mask(k)
        F = best_lower(k)
        psi = potential(k, m, F)
        T = tight_edges(k, m, F, psi)
        adj = {}
        for s, t in T:
            adj.setdefault(s, []).append(t)
        for s in range(1 << k):
            adj.setdefault(s, [])
        cycles = simple_cycles_adj(adj)
        words = {}
        for cyc in cycles:
            cyc = as_cycle(k, m, cyc)
            check(density(cyc) == F, "tight cycles have density F_k")
            bl = block_decomposition(k, m, cyc)
            check(bl is not None, "block decomposition")
            for (i0, L, kind) in bl:
                seg = [cyc[(i0 + j) % len(cyc)] & 1 for j in range(L)]
                if kind == 'C':
                    check(is_fd_word(seg) and Fraction(sum(seg), L) == F, "C-block of density F_k")
                else:
                    check(seg == [1, 0] and F == Fraction(1, 2), "F-block only when F_k = 1/2")
            cw = canon(cycle_word(cyc))
            words[cw] = words.get(cw, 0) + 1
        chris = canon(''.join(map(str, upper_christoffel(F.numerator, F.denominator))))
        non_chris = [w for w in words if canon(primitive_root(w)) != chris]
        if k <= 3:
            check(len(non_chris) == 0, "k <= 3: all maximizers carry (powers of) the Christoffel word")
        else:
            check(len(non_chris) > 0, "k >= 4: non-Christoffel maximizers exist")
            check(('1100' in words) if k == 4 else ('11100' in words), "the named non-Christoffel maximizer, k=%d" % k)
        print("    %d : %d : %s" % (k, len(cycles), sorted(words.items(), key=lambda x: (len(x[0]), x[0]))[:12]))
    claim(True, "k = 2..7: every simple cycle of the tight subgraph has density F_k and factors into first-descent blocks "
          "of density F_k (and flip blocks 10, only when F_k = 1/2); for k = 2, 3 all of them carry the Christoffel word "
          "10, while for k = 4..7 non-Christoffel maximizers occur (1100 at k = 4, 11100 at k = 5..7)")
    # (5) entropy of the maximizing set: Perron root of the tight subgraph vs the block equation sum z^|w| = 1
    print("    k : Perron root of the tight subgraph : root of sum_blocks z^-|w| = 1 (k >= 5)")
    for k in range(4, 12):
        m = sigma_k_mask(k)
        F = best_lower(k)
        psi = potential(k, m, F)
        A = np.zeros((1 << k, 1 << k))
        for s, t in tight_edges(k, m, F, psi):
            A[s, t] = 1
        lam = max(abs(np.linalg.eigvals(A)))
        blocks = fd_blocks(F, k)
        lens = [len(w) for w in blocks]
        if k >= 5:
            lo, hi = 1.0, 2.0
            for _ in range(200):
                mid = (lo + hi) / 2
                if sum(mid ** (-L) for L in lens) > 1:
                    lo = mid
                else:
                    hi = mid
            check(abs(lam - lo) < 1e-6, "Perron root = block growth rate, k=%d" % k)
            print("    %2d : %.6f : %.6f" % (k, lam, lo))
        else:
            check(lam > 1 + 1e-6, "positive entropy at k = 4")
            print("    %2d : %.6f : (flip blocks 10 also tight)" % (k, lam))
    claim(True, "k = 4..11: the maximizing set has positive entropy; for k >= 5 the Perron root of the tight subgraph "
          "equals the growth rate of block concatenations (numerical, 1e-6): e.g. 2^(1/5) = 1.1487 for k = 5..7, "
          "7^(1/8) = 1.2754 for k = 8..11")
    # (6) the global maximizers over class (i) at levels 4, 5 are not Christoffel either
    for k, expect in ((4, {'11100'}), (5, {'11100110', '11111000', '11101100', '11110010'})):
        R = rho_all(k)
        Mk = max(f for f in R if below_c(f))
        opt = [m for m in range(len(R)) if R[m] == Mk]
        words = set()
        for m in opt:
            adj = {s: list(succ(k, m, s)) for s in range(1 << k)}
            for cyc in simple_cycles_adj(adj):
                if density(cyc) == Mk:
                    words.add(canon(cycle_word(as_cycle(k, m, cyc))))
        chris = canon(''.join(map(str, upper_christoffel(Mk.numerator, Mk.denominator))))
        check(words == expect and chris not in words, "global maximizers' words, k=%d" % k)
        print("    level %d: M_k = %s, %d optimal strategies, max-cycle words %s, Christoffel %s absent"
              % (k, Mk, len(opt), sorted(words), chris))
    claim(True, "the class-(i) strategies of largest rho_max at levels 4 and 5 (M_4 = 3/5, M_5 = 5/8) have no "
          "Christoffel maximizer at all: their maximal cycles carry 11100 (k=4) and 11100110, 11111000, 11101100, "
          "11110010 (k=5)")


if __name__ == "__main__":
    run()
