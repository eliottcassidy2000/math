#!/usr/bin/env python3
"""collatz_landing_20260926_probe.py -- the landing multiplicity of THM-4476's recursion, probed
(session collatz-landing-20260926, opus, 2026-09-26).

(A) Strip entropy. THM-4476/4499 count words by one-sided barriers (entropy h = h(log_3 2) = 0.94996 per
    letter). A segment that HOVERS in a band of width W bits is a word with all partial sums in [-W, 0] (up
    to a shift); such two-sided words are exponentially rarer: count ~ 2^(h_W m) with h_W < h. We estimate
    h_W for W = 1, 2, 3, 4, 6, 8 bits from exact DP counts (words of length m whose partial sums S_j =
    o_j log_2 3 - j all lie in [-W, 0]) by the growth rate between m = M1 and m = M2.
(B) Landing multiplicity on actual orbit segments. For a segment y_0..y_M and X = 2^L, theta: a (D) index
    i (y_i <= X) has first s <= L with y_(i+s) < y_i X^(-theta); the landing point is j = i + s. We tabulate
    the multiplicity distribution of landing points, the average multiplicity, and the fraction of (D) among
    indices <= X, for (i) the 5n+1 orbit of 7 (genuinely divergent? unknown; it grows), (ii) long 3n+1
    segments (orbit of 27, and of 2^L - 1 which climbs L steps), (iii) a synthetic hover-then-drop segment
    built from a residue class (Terras), to show that multiplicity ~ L is realized by integers.
Usage: python3 collatz_landing_20260926_probe.py [M2=300]
"""
import math, sys

ALPHA = math.log2(3)


def strip_count(m, W):
    """words of length m with all partial sums in [-W, 0] (float comparison; no ties for W integer except o=0,
    where S_j = -j >= -W iff j <= W: compare exactly)."""
    cur = {0: 1}
    for j in range(1, m + 1):
        nxt = {}
        for o, c in cur.items():
            for st in (0, 1):
                o2 = o + st
                if o2 == 0:
                    ok = j <= W
                else:
                    v = o2 * ALPHA - j
                    ok = (v <= 0) and (v >= -W)
                if ok:
                    nxt[o2] = nxt.get(o2, 0) + c
        cur = nxt
    return sum(cur.values())


def T(n, q=3, b=1):
    return n // 2 if n % 2 == 0 else (q * n + b) // 2


def orbit(n, steps, q=3, b=1):
    ys = [n]
    for _ in range(steps):
        n = T(n, q, b)
        ys.append(n)
    return ys


def landing_stats(ys, L, theta):
    X = 2 ** L
    thr = 2.0 ** (-theta * L)
    mult = {}
    nd = 0
    total = 0
    for i, y in enumerate(ys):
        if y > X or i + L >= len(ys):
            continue
        total += 1
        land = None
        for s in range(1, L + 1):
            if ys[i + s] < y * thr:
                land = i + s
                break
        if land is None:
            nd += 1
        else:
            mult[land] = mult.get(land, 0) + 1
    ms = sorted(mult.values(), reverse=True)
    D = sum(ms)
    return total, nd, D, len(ms), (ms[:8] if ms else []), (D / len(ms) if ms else 0.0)


def synthetic_hover(L, K):
    """an integer whose word keeps its partial sums in [-1.49, 0] for K steps (a strict [-1, 0] hover of length >= 3 is impossible) and then has L/4 zeros and L/4 ones:
    built by choosing the parity word and solving for the residue (Terras)."""
    word = []
    S = 0.0
    for j in range(K):
        # choose the letter keeping S in [-1, 0]
        if S + (ALPHA - 1) <= 0:
            word.append(1); S += ALPHA - 1
        else:
            word.append(0); S -= 1
    word += [0] * (L // 4) + [1] * (L // 4)
    k = len(word)
    # residue r mod 2^k with the given parity word: iterate
    r = 0
    for j in range(k):
        # find bit: T^j(r + 2^j c) parity must equal word[j]; parity of T^j depends on the j-th bit
        x = r
        for _ in range(j):
            x = T(x)
        if x % 2 != word[j]:
            r += 2 ** j
    return r + 2 ** k  # a positive representative


def main():
    M2 = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    M1 = M2 // 2
    h = -(math.log(2) / math.log(3)) * math.log2(math.log(2) / math.log(3)) - (1 - math.log(2) / math.log(3)) * math.log2(1 - math.log(2) / math.log(3))
    print("(A) strip entropy h_W (growth rate of #words with all partial sums in [-W, 0], m = %d -> %d), h = %.5f:" % (M1, M2, h))
    print("   W=1: the band [-1, 0] admits no word of length >= 3 (steps -1 and +0.585 span 1.585 > 1): count(m=150) = %d" % strip_count(150, 1))
    for W in (2, 3, 4, 6, 8, 12, 16):
        c1, c2 = strip_count(M1, W), strip_count(M2, W)
        if c1 == 0 or c2 == 0:
            print("   W=%2d bits: no words (count %d, %d)" % (W, c1, c2)); continue
        hW = (math.log2(c2) - math.log2(c1)) / (M2 - M1)
        print("   W=%2d bits: h_W = %.4f   (1 - h_W = %.4f vs 1 - h = %.4f); count(m=%d) = 2^%.2f" % (W, hW, 1 - hW, 1 - h, M2, math.log2(c2)))
    print("(B) landing multiplicities: total indices <= X, #ND, #D, #landing points, largest multiplicities, mean multiplicity")
    L = 20
    for theta in (0.05, 0.1, 0.2):
        print("  L = %d, theta = %.2f (X^theta = %.1f):" % (L, theta, 2 ** (theta * L)))
        cases = [
            ("5n+1 orbit of 7 (2000 steps)", orbit(7, 2000, q=5)),
            ("3n+1 orbit of 27 (until 1, padded)", orbit(27, 200)),
            ("3n+1 orbit of 2^18-1 (climbs 18 steps, then 400)", orbit(2 ** 18 - 1, 400)),
            ("3n+1 orbit of 2^19-1 (until 1-ish)", orbit(2 ** 19 - 1, 600)),
            ("synthetic hover(K=15)+drop", orbit(synthetic_hover(L, 15), 200)),
        ]
        for name, ys in cases:
            tot, nd, D, lp, top, mean = landing_stats(ys, L, theta)
            print("     %-46s total=%4d ND=%4d D=%4d landing=%3d top=%s mean=%.2f" % (name, tot, nd, D, lp, top, mean))


if __name__ == '__main__':
    main()
