#!/usr/bin/env python3
"""collatz_oscillation_20260926_shells.py -- the oscillation lemma as a shell-visit statement, tested
(session collatz-oscillation-20260926, opus, 2026-09-26).

THM-4506: all dippers of a landing point j lie in one dyadic shell S_j = (2^D y_j, 2^(D+1) y_j] and consecutive
dippers are separated by an odd letter; so the multiplicity m(j) is the number of (odd-separated) visits of the
orbit to the shell S_j during the k steps before j that are not preceded by a drop. HYP-9161 in its averaged
form (#Dip(X,D) <= C L^mu N(X 2^-D) + C L^2) is therefore a statement about how often an orbit revisits one
dyadic shell within a window. This script measures, on long 3n+1 and 5n+1 segments and for the bootstrap's own
depth D = ceil(1.05 log_2 L):
  (i)  the multiplicity distribution and its mean over landing points, and the number of shell visits (returns
       to S_j inside the window) for the heaviest landing points;
  (ii) the scale regularity R(X) = N(X) / N(X 2^-D) and the fraction of elements <= X that are dippers, which is
       what the averaged hypothesis constrains: #Dip(X,D) <= C L^mu N(X 2^-D) forces N(X) <= (C L^mu + 1) N(X 2^-D) + #ND + C L^2.
Usage: python3 collatz_oscillation_20260926_shells.py
"""
import math, random


def T(n, q=3):
    return n // 2 if n % 2 == 0 else (q * n + 1) // 2


def orbit(n, steps, q=3):
    ys = [n]
    for _ in range(steps):
        n = T(n, q); ys.append(n)
        if n == 1: break
    return ys


def stats(ys, L, D):
    X = 2 ** L; thr = 2.0 ** (-D)
    mult = {}; nd = 0; tot = 0; below_lower = 0
    for i, y in enumerate(ys):
        if y <= X * thr and i + L < len(ys):
            below_lower += 1
        if y > X or i + L >= len(ys):
            continue
        tot += 1
        land = None
        for s in range(1, L + 1):
            if ys[i + s] < y * thr:
                land = i + s; break
        if land is None:
            nd += 1
        else:
            mult.setdefault(land, []).append(i)
    ms = sorted(((len(v), j, v) for j, v in mult.items()), reverse=True)
    # shell visits for the heaviest landing points: count maximal runs of window indices inside S_j
    heavy = []
    for m, j, dippers in ms[:3]:
        yj = ys[j]; lo, hi = (2 ** D) * yj, (2 ** (D + 1)) * yj
        runs = 0; inside = False
        for t in range(max(0, j - L), j):
            now = lo < ys[t] <= hi
            if now and not inside: runs += 1
            inside = now
        heavy.append((m, runs, j))
    D_total = sum(m for m, _, _ in ms)
    return tot, nd, D_total, len(ms), (D_total / len(ms) if ms else 0.0), heavy, below_lower


def main():
    random.seed(20260926)
    print("L    D  segment                       total   ND     D  landing  mean  heaviest (mult, shell runs)   N(X)/N(X 2^-D)")
    for L in (24, 32, 48, 64):
        D = math.ceil(1.05 * math.log2(L))
        segs = [("3n+1 record 63728127", orbit(63728127, 60 * L)), ("3n+1 from 2^L-1", orbit(2 ** L - 1, 60 * L)),
                ("3n+1 random odd ~2^L", orbit(random.randrange(2 ** (L - 1), 2 ** L) | 1, 60 * L)),
                ("5n+1 from 7", orbit(7, 60 * L, q=5)), ("5n+1 random odd ~2^(L/2)", orbit(random.randrange(2 ** (L // 2 - 1), 2 ** (L // 2)) | 1, 60 * L, q=5))]
        for name, ys in segs:
            tot, nd, Dt, lp, mean, heavy, bl = stats(ys, L, D)
            ratio = (tot / bl) if bl else float('inf')
            print("%2d  %2d  %-28s %5d %5d %5d %7d  %5.2f  %-28s %6.2f" % (L, D, name, tot, nd, Dt, lp, mean, " ".join("(%d,%d)" % (m, r) for m, r, _ in heavy), ratio))


if __name__ == '__main__':
    main()
