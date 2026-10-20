#!/usr/bin/env python3
"""Small exact combinatorial controls; numerical controls of analytic formulas.

No repository mathematics code is imported. The sibling verifier supplies
independent path DP and the earlier proved scalar word recurrence.
"""
from bisect import bisect_right
from collections import Counter
from fractions import Fraction
from itertools import product
from math import log

from verify import h_cut, h_dp, h_word, need, tournament


K = 4 * log(2)


def F(z):
    if z <= 0.25:
        return 0.0
    if z <= 0.5:
        return 4*z*log(4*z)-4*z+1
    if z <= 1:
        return 4*z-3-4*z*log(z)
    return 1.0


def error_bound(L):
    return K * (2**(L//2) + 2**((L+1)//2) + (2*L+7)/8)


def binary(bits):
    n = 0
    for bit in bits:
        n = 2*n+bit
    return n


def main():
    by_length = {}
    rows = []
    exact = dp = 0
    for m in range(2, 13):
        s = (m-2)//2
        prefixes = Counter()
        counts = Counter()
        fixed = 0
        for middle in product((0, 1), repeat=m-2):
            b = (1,)+middle+(0,)
            dual = tuple(1-z for z in reversed(b))
            h = h_word(b)
            need(h == h_cut(b) == h_word(dual), ('word/duality', b))
            N = binary(b)
            R = binary(dual)
            A = sum(b[i] == 1 and b[i+1] == 0 for i in range(m-1))
            Bscaled = sum(2**(m+j-i-1) for i in range(m) for j in range(i)
                          if b[i] and not b[j])
            difference = h*2**m-N*R
            need(difference == (1+A)*2**m-Bscaled, ('identity', b))
            need(2*A <= m and 2*Bscaled < m*2**m, ('error ingredients', b))
            need(2*abs(difference) <= (m+2)*2**m, ('error bound', b))
            need(2**(m-2) < h < 2**m, ('cutoff', b, h))
            a_prefix = binary(b[1:s+1])
            c_prefix = binary(dual[1:s+1])
            prefixes[a_prefix, c_prefix] += 1
            # The actual coordinates lie in the cells used for the coupling.
            for numerator, prefix in [(N, a_prefix), (R, c_prefix)]:
                bottom = 2**(m-1)+prefix*2**(m-s-1)
                need(bottom <= numerator < bottom+2**(m-s-1),
                     ('prefix cell', b, numerator, prefix))
            fixed += b == dual
            counts[h] += 1
            exact += 1
            if m <= 8:
                need(h == h_dp(tournament(b)), ('independent path DP', b))
                dp += 1
        need(len(prefixes) == 2**(2*s), ('prefix support', m))
        need(set(prefixes.values()) == {2**(m-2-2*s)}, ('prefix independence', m))
        need(fixed == (0 if m % 2 else 2**(m//2-1)), ('fixed words', m))
        total = 2**(m-2)
        if m >= 3:
            exact_mean = Fraction(sum(h*n for h, n in counts.items()), total)
            need(exact_mean == Fraction(9*2**m, 16)+Fraction(1, 2),
                 ('exact mean', m, exact_mean))
        seen = 0
        distance = 0.0
        for h, multiplicity in sorted(counts.items()):
            cdf = F(h/2**m)
            distance = max(distance, abs(seen/total-cdf))
            seen += multiplicity
            distance = max(distance, abs(seen/total-cdf))
        delta = 2**(-s)+(1+m/2)*2**(-m)
        need(distance <= K*delta+1e-12, ('uniform CDF bound', m, distance))
        by_length[m] = sorted(counts.elements())
        rows.append((m, round(distance, 9)))
    print('exact reduced-word identities, cutoff, duality and prefix cells:', exact)
    print('independent subset-DP path checks:', dp)
    print('prefix independence and fixed-word counts: all lengths 2..12 PASS')
    print('exact reduced-word mean: all lengths 3..12 PASS')
    print('finite CDF distances by length:', rows)

    phase_checks = 0
    for j in range(65):
        t = 1+j/64
        need(abs(1+F(t/2)+2*F(t/4)-2*t*log(2)) < 3e-14,
             ('phase cancellation', t))
        phase_checks += 1
    print('non-dyadic phase cancellation controls:', phase_checks)

    bank = sorted(h for values in by_length.values() for h in values)
    spectrum = sorted(set(bank) | {1})
    cutoffs = [Fraction(x) for x in range(4, 2048)]
    cutoffs += [Fraction((2**L)*j, 8) for L in range(2, 11) for j in (9, 11, 13, 15)]
    for X in cutoffs:
        L = (X.numerator//X.denominator).bit_length()-1
        W = bisect_right(bank, X)
        exact_two_length = (2**(L-1)-1
                            +bisect_right(by_length[L+1], X)
                            +bisect_right(by_length[L+2], X))
        need(W == exact_two_length, ('two lengths', X, W, exact_two_length))
        E = W-float(X)*log(2)+1
        need(abs(E) <= error_bound(L)+1e-10, ('all-cutoff error', X, E))
        upper = float(X)*log(2)/2+error_bound(L)/2+2**((L+2)//2-1)
        need(bisect_right(spectrum, X) <= upper+1e-10, ('spectrum upper bound', X))
    print('exact two-length/all-cutoff error and spectrum-bound controls:', len(cutoffs))
    print('word-count ratios W(X)/X at X=16,32,...,1024:',
          [(2**L, round(bisect_right(bank, 2**L)/2**L, 9)) for L in range(4, 11)])
    print('ALL CHECKS PASSED')


if __name__ == '__main__':
    main()
