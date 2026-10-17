#!/usr/bin/env python3
"""The Legendre walk of a Paley tournament: area = class number, positivity = local optimality of the natural order (S40).

For p = 3 (mod 4) prime let S(t) = sum_(j <= t) (j/p). Moving vertex i of the natural order of P_p forward or backward
over t vertices changes the forward-minus-backward excess by exactly -2 S(t); so the natural order admits no improving
single-vertex insertion iff S(t) >= 0 for all t (Kalmynin's set L+, arXiv:2105.06910), and its excess is
sum_t S(t) = p h(-p). By Fekete's identity f_p(x) = (1 - x) sum_t S(t) x^t, a nonnegative walk makes the Fekete
polynomial positive on (0, 1), so L(s, chi_(-p)) has no real zero s in (0, 1) (Fekete 1912; Polya).
This script tabulates, for every p = 3 (mod 4) prime up to PMAX: membership in L+, min_t S(t), the area check
sum_t S(t) = p h(-p), and L(1, chi) = pi h / sqrt(p); and summarises L+ against L(1, chi).
Usage: python paley_legendre_walk_positivity_20261017.py PMAX
Session opus-2026-10-17-S40.
"""
import sys, os, json, math
import numpy as np
from numba import njit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paley_class_number_ranking_20261017 import class_number, is_prime


@njit(cache=True)
def walk(p):
    chi = np.full(p, -1, dtype=np.int64)
    chi[0] = 0
    for x in range(1, (p - 1) // 2 + 1):
        chi[(x * x) % p] = 1
    s = 0; area = 0; mn = 10 ** 9; first_neg = -1
    for t in range(1, p):
        s += chi[t]
        area += s
        if s < mn:
            mn = s
        if s < 0 and first_neg < 0:
            first_neg = t
    return mn, area, first_neg


def main():
    PMAX = int(sys.argv[1])
    rows = []
    bad_area = 0
    for p in range(7, PMAX + 1, 4):
        if not is_prime(p):
            continue
        mn, area, fn = walk(p)
        h = class_number(-p)
        if area != p * h:
            bad_area += 1
        rows.append((p, h, math.pi * h / math.sqrt(p), int(mn), int(fn)))
    Lplus = [r for r in rows if r[3] >= 0]
    print('p = 3 (mod 4) primes 7..%d: %d; area = p h(-p) fails at %d' % (PMAX, len(rows), bad_area))
    print('L+ (natural order insertion-locally optimal): %d primes (%.4f); first ones: %s'
          % (len(Lplus), len(Lplus) / len(rows), [r[0] for r in Lplus[:30]]))
    for lo, hi in [(0, 1), (1, 1.5), (1.5, 2), (2, 2.5), (2.5, 3), (3, 10)]:
        sel = [r for r in rows if lo <= r[2] < hi]
        inn = [r for r in sel if r[3] >= 0]
        print('  L(1,chi) in [%.1f, %.1f): %6d primes, in L+ %6d (%.3f)' % (lo, hi, len(sel), len(inn),
                                                                        len(inn) / max(1, len(sel))))
    notL = [r for r in rows if r[3] < 0]
    print('largest L(1,chi) outside L+: %s' % sorted(notL, key=lambda r: -r[2])[:5])
    print('smallest L(1,chi) inside L+: %s' % sorted(Lplus, key=lambda r: r[2])[:5])
    # density decay of L+ by dyadic ranges
    for k in range(3, int(math.log2(PMAX)) + 1):
        sel = [r for r in rows if 2 ** (k - 1) < r[0] <= 2 ** k]
        if sel:
            print('  p in (2^%d, 2^%d]: %5d primes, L+ fraction %.4f' % (k - 1, k, len(sel),
                                                                        sum(r[3] >= 0 for r in sel) / len(sel)))
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(dict(PMAX=PMAX, n=len(rows), area_failures=bad_area, Lplus=[r[0] for r in Lplus]),
              open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'paley_legendre_walk_positivity_20261017.json'), 'w'))


if __name__ == '__main__':
    main()
