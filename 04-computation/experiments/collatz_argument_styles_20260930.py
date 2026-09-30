#!/usr/bin/env python3
"""collatz_argument_styles_20260930.py -- exact reformulations behind the argument-style proposals of the thirteenth note
(session collatz-posets-zeta5-20260927, opus, 2026-09-30, part C).

 (1) The Hadamard identity: with F_n(z) = sum 2^(d_k) z^k (valuation series) and G_n(z) = sum m_k z^k (Syracuse orbit series),
     (1 - 3z) (F_n * G_n)(z) = n + z F_n(z), where * is the Hadamard (termwise) product.  Checked to order 40.
 (2) The weighted-mediant theorem: x_(uv) = (alpha x_u + beta x_v)/(alpha + beta) with alpha = 3^(p_v) D_u, beta = 2^(A_u) D_v,
     and the Farey determinant S_u D_v - S_v D_u = -beta(u,v) of the twelfth note; the letter points 1/(2^v - 3) are the
     one-step rational cycles; words without 1-steps have fixed points in (0, 1] (convexity); census of the rational cycles
     (denominators of x_w) over all words with A <= 18, and the positive near-misses (x_w > 1 closest to an integer).
 (3) The never-descending Cantor set K: K ∩ [-10^5, -1] = {-1, -5, -17} and K ∩ [1, 10^5] = {} (every positive n descends;
     the record stopping times), so Collatz over Z is "K ∩ Z = {-1, -5, -17}".
Usage: python3 collatz_argument_styles_20260930.py
"""
import math, random
from fractions import Fraction
from math import gcd
from collections import Counter


def syracuse_word(n, steps):
    """valuations of the Syracuse orbit of odd n (n != 0), and the orbit."""
    x = n; word = []; orbit = [n]
    for _ in range(steps):
        y = 3 * x + 1; v = 0
        while y % 2 == 0:
            y //= 2; v += 1
        word.append(v); x = y; orbit.append(x)
    return word, orbit


def part1():
    print("== (1) the Hadamard identity (1 - 3z)(F * G) = n + z F ==")
    ok = True
    for n in (7, 27, 97, 871, 6171, 77031):
        word, orbit = syracuse_word(n, 40)
        d = [0]
        for v in word:
            d.append(d[-1] + v)
        F = [2 ** d[k] for k in range(41)]
        H = [F[k] * orbit[k] for k in range(41)]        # Hadamard product coefficients
        lhs = [H[k] - 3 * (H[k - 1] if k else 0) for k in range(41)]
        rhs = [n] + [F[k - 1] for k in range(1, 41)]
        ok &= lhs == rhs
    print(" identity holds to order 40 for n = 7, 27, 97, 871, 6171, 77031: %s" % ok)
    print(" reading: the orbit series is the termwise quotient of the rational transform (n + zF)/(1 - 3z) by F; Collatz over N is 'G_n rational for every n', and F_n rational iff G_n rational (Polya / Prop 10 of the first note)")


def carry(w):
    S = 0; d = 0
    for v in w:
        S = 3 * S + 2 ** d; d += v
    return S


def clock(w):
    return 2 ** sum(w) - 3 ** len(w)


def x_of(w):
    return Fraction(carry(w), clock(w))


def part2():
    print("== (2) the weighted-mediant theorem and the rational cycles ==")
    random.seed(11); ok = True; okdet = True
    for _ in range(500):
        u = tuple(random.randint(1, 4) for _ in range(random.randint(1, 5)))
        v = tuple(random.randint(1, 4) for _ in range(random.randint(1, 5)))
        alpha = 3 ** len(v) * clock(u); beta = 2 ** sum(u) * clock(v)
        ok &= x_of(u + v) == (alpha * x_of(u) + beta * x_of(v)) / (alpha + beta)
        okdet &= carry(u) * clock(v) - carry(v) * clock(u) == -(carry(u + v) - carry(v + u))
    print(" x_(uv) = (3^(p_v) D_u x_u + 2^(A_u) D_v x_v)/D_(uv): %s; Farey determinant S_u D_v - S_v D_u = -beta(u,v): %s" % (ok, okdet))
    u, v = (1, 1, 1), (2, 1, 1, 4)
    alpha = 3 ** 4 * clock(u); beta = 2 ** 3 * clock(v)
    print(" the -17 word = (1,1,1)(2,1,1,4): x_u = %s, x_v = %s, weights alpha = %d, beta = %d, mediant = %s" % (x_of(u), x_of(v), alpha, beta, (alpha * x_of(u) + beta * x_of(v)) / (alpha + beta)))
    print(" letter points x_(v) = 1/(2^v - 3) (the one-step rational cycles):", [str(x_of((v,))) for v in range(1, 8)])
    # census over all words with A <= 18
    dens = Counter(); words_by_den = {}; near = []; noone_ok = True; total = 0
    for A in range(1, 19):
        for cuts in range(1 << (A - 1)):
            # composition of A from the bit pattern: a cut after position i iff bit i set
            w = []; run = 1
            for i in range(A - 1):
                if (cuts >> i) & 1:
                    w.append(run); run = 1
                else:
                    run += 1
            w.append(run); w = tuple(w); total += 1
            S, D = carry(w), clock(w); g = gcd(S, abs(D)); den = abs(D) // g
            dens[den] += 1
            if den <= 15 and len(words_by_den.get(den, [])) < 6:
                words_by_den.setdefault(den, []).append((w, str(Fraction(S, D))))
            if 1 not in w:
                noone_ok &= 0 < Fraction(S, D) <= 1
            x = Fraction(S, D)
            if x > 1 and D > 0:
                N = round(x)
                if N >= 2:
                    near.append((abs(x - N), w, str(x), N))
    print(" words with A <= 18: %d; words without 1-steps have x_w in (0, 1]: %s" % (total, noone_ok))
    print(" denominators of x_w (rational cycles), smallest 20 with counts:", sorted(dens.items())[:20])
    for den in sorted(words_by_den)[:9]:
        print("  den %2d: %s" % (den, words_by_den[den][:4]))
    near.sort()
    print(" closest positive near-misses x_w > 1 to an integer N >= 2 (gap, word, x_w, N):", [(float(g), w, x, N) for g, w, x, N in near[:6]])


def part3():
    print("== (3) the never-descending set K and the integers ==")
    def first_word_descent(n, steps=3000):
        """K is defined by the word: no descent means 3^j >= 2^(d_j) for every prefix (the sibling ladder's Bad(3)).
        Returns the first j with 2^(d_j) > 3^j, or None if none within `steps` odd steps (or the word became periodic)."""
        x = n; d = 0; p3 = 1
        for j in range(1, steps + 1):
            y = 3 * x + 1; v = 0
            while y % 2 == 0:
                y //= 2; v += 1
            d += v; p3 *= 3
            if 2 ** d > p3:
                return j
            x = y
            if x == n:
                return None
        return None
    neg = [n for n in range(-1, -100001, -2) if first_word_descent(n) is None]
    print(" negative odd integers n >= -10^5 in K (word never descends): %s" % neg)
    worst = (0, 0); stuck = []
    for n in range(1, 100001, 2):
        j = first_word_descent(n)
        if j is None:
            stuck.append(n)
        elif j > worst[0]:
            worst = (j, n)
    print(" positive odd n <= 10^5 in K: %s (1 is not in K: its word (2)^inf descends at once); record first word-descent: n = %d at odd step %d" % (stuck, worst[1], worst[0]))
    print(" reading: Collatz over Z is the statement K cap Z = {-1, -5, -17}; K is a closed 2-adic set of Hausdorff dimension h(log_3 2) = 0.95 (sibling ladder Thm 1a) that contains three negative integers and, conjecturally, no positive one: the distinction is the sign, invisible 2-adically")


if __name__ == "__main__":
    part1(); part2(); part3()


def part4():
    print("== (4) local equidistribution of the carries, and the absence of unimodular word pairs ==")
    # carries of all words of shape (A, p) modulo small primes l: max relative deviation from uniform
    A = 18
    from itertools import combinations
    for p in (9, 11, 12):
        cnt = {l: [0] * l for l in (5, 7, 11, 13)}
        total = 0
        base = pow(3, p - 1)
        for cuts in combinations(range(1, A), p - 1):
            S = base
            for i, d in enumerate(cuts, start=1):
                S += pow(3, p - 1 - i) * (1 << d)
            total += 1
            for l in cnt:
                cnt[l][S % l] += 1
        devs = {l: round(max(abs(c * l / total - 1) for c in cnt[l]), 4) for l in cnt}
        print(" shape (%d,%d): %d words; max relative deviation of S_w mod l from uniform: %s" % (A, p, total, devs))
    # unimodular pairs: min |beta(u,v)| = |S_u D_v - S_v D_u| over words with A <= 8
    words = []
    for A in range(1, 9):
        for cuts in range(1 << (A - 1)):
            w = []; run = 1
            for i in range(A - 1):
                if (cuts >> i) & 1:
                    w.append(run); run = 1
                else:
                    run += 1
            w.append(run); words.append(tuple(w))
    best = None; parity_even = True
    for u in words:
        Su, Du = carry(u), clock(u)
        for v in words:
            if v <= u:
                continue
            b = Su * clock(v) - carry(v) * Du
            parity_even &= (b % 2 == 0)
            if b != 0 and (best is None or abs(b) < best[0]):
                best = (abs(b), u, v)
    print(" over %d words with A <= 8: the Farey determinant S_u D_v - S_v D_u is always even: %s; minimum nonzero |det| = %d at %s, %s (no unimodular pairs: the twisted mediant tree is not Stern-Brocot)" % (len(words), parity_even, best[0], best[1], best[2]))


if __name__ == "__main__":
    part4()
