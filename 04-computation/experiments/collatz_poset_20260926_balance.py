#!/usr/bin/env python3
"""collatz_poset_20260926_balance.py -- the no-descent words as linear extensions of a width-2 poset, and
the 1/3-2/3 balance constant of these "Collatz posets" (session collatz-crossings-20260926, opus, 2026-09-26).

A parity word of length k with o odd letters (O_1 < ... < O_o in order of occurrence) and e = k - o halvings
(H_1 < ... < H_e) has no descent (all prefix multipliers 3^(#odd)/2^j > 1) iff for every prefix the number of
halvings e' and odd letters o' satisfy 2^(e') ... precisely 3^(o') > 2^(o'+e'), i.e. e' < o' log_2(3/2)... in
integers: the j-th halving H_j must be preceded by at least f(j) odd letters, where f(j) = least o' with
3^(o') > 2^(o'+j), i.e. H_j > O_(f(j)). So

    Bad_k(o) = { no-descent words with o odd letters } = linear extensions of the poset P(o, e):
               two chains O_1 < ... < O_o and H_1 < ... < H_e, plus the cross relations O_(f(j)) < H_j.

P(o, e) has width 2, so Linial's theorem (the 1/3-2/3 conjecture for width 2) applies: some incomparable pair
(x, y) has P(x < y) in [1/3, 2/3] over uniformly random linear extensions. This script enumerates the linear
extensions (= the no-descent words) for k <= KMAX, computes for every incomparable pair (O_i, H_j) the
probability that O_i precedes H_j, and reports the balance constant delta(P) = max over incomparable pairs of
min(p, 1 - p), the most balanced pair, and the number of linear extensions e(P) (= C(k,o) restricted to no
descent). It also reports delta for the union over o (the whole Bad_k, which is not a poset's extension set:
different o give different posets) as the best single "comparison" O_i vs H_j across Bad_k.
Usage: python3 collatz_poset_20260926_balance.py [KMAX=22]
"""
import math, sys
from itertools import combinations


def f_of(j):
    o = 0
    while not (3 ** o > 2 ** (o + j)):
        o += 1
    return o


def nodescent_words(k):
    """all no-descent words of length k as tuples of letters (1 = odd, 0 = halving); exact test 3^o > 2^j."""
    out = []
    def rec(prefix, o, j):
        if j == k:
            out.append(tuple(prefix)); return
        for step in (1, 0):
            o2 = o + step
            if 3 ** o2 > 2 ** (j + 1):
                prefix.append(step); rec(prefix, o2, j + 1); prefix.pop()
    rec([], 0, 0)
    return out


def main():
    KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 22
    print("f(j) = least number of odd letters that must precede the j-th halving:", [f_of(j) for j in range(1, 13)])
    print("k   o    e(P)=#words   delta(P)   most balanced pair (O_i before H_j: p)   pairs>=1/3   width-2 check")
    for k in range(4, KMAX + 1):
        words = nodescent_words(k)
        by_o = {}
        for w in words:
            by_o.setdefault(sum(w), []).append(w)
        best_overall = (0.0, None)
        for o in sorted(by_o):
            ws = by_o[o]
            e = k - o
            n = len(ws)
            # position lists
            best = (0.0, None); cnt13 = 0; npairs = 0
            for i in range(1, o + 1):
                for j in range(1, e + 1):
                    # comparable if O_i < H_j forced (i <= f(j)) or ... H_j < O_i is never forced; O_i < H_j forced iff i <= f(j)
                    if i <= f_of(j):
                        continue
                    npairs += 1
                    c = 0
                    for w in ws:
                        # position of the i-th odd and the j-th halving
                        oi = hj = None; co = ch = 0
                        for pos, letter in enumerate(w):
                            if letter == 1:
                                co += 1
                                if co == i: oi = pos
                            else:
                                ch += 1
                                if ch == j: hj = pos
                            if oi is not None and hj is not None: break
                        if oi < hj: c += 1
                    p = c / n
                    bal = min(p, 1 - p)
                    if bal >= 1 / 3 - 1e-12: cnt13 += 1
                    if bal > best[0]: best = (bal, (i, j, p))
            if best[1] is not None and best[0] > best_overall[0]:
                best_overall = (best[0], (o,) + best[1])
            if n >= 2:
                print("%2d  %2d   %10d   %.4f     O_%d before H_%d: p = %.4f            %3d/%-3d   %s" % (
                    k, o, n, best[0], best[1][0], best[1][1], best[1][2], cnt13, npairs, "ok (>=1/3)" if best[0] >= 1/3 - 1e-12 else "VIOLATION"))
        print("    k=%d total |Bad_k| = %d; best balanced comparison over o: %s" % (k, len(words), best_overall))


if __name__ == '__main__':
    main()
