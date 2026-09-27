#!/usr/bin/env python3
"""collatz_shadow_error_20260927_itineraries.py -- the bounded-shadow regime B = 2: itineraries of
f(eta) = (3 eta - 1)/2 on [1, 5/3), (3 eta - 1)/4 on [5/3, 2), their Bernstein self-consistency, and the 2-adic
points they define (session collatz-shadow-flp-20260927, opus, 2026-09-27).

For rational eta_0 = a/b in [1, 2) (b <= BMAX) the itinerary w = (v_1, v_2, ...) is computed exactly to length L.
 (a) growth: log_2(3^L / 2^(d_L)) > 0 means the word is a growth word over the window;
 (b) self-consistency: the partial Bernstein sums sum_(k<K) 2^(d_k)/3^(k+1) -> eta_0 (they must, by the uniqueness
     of bounded solutions of 2^v eta' = 3 eta - 1 along a growth word);
 (c) the 2-adic point x(w): the residue rho_K = (2^(d_K) - S_K) 3^(-K) mod 2^(d_K + 1) of the exact K-prefix
     class (THM-4512's exact cylinder); x(w) is a positive integer iff rho_K stabilises to a fixed value for all
     large K; we report whether rho_K is constant over K in [L/2, L] (it never is for the samples, as expected:
     the itineraries' 2-adic points are not integers).
 (d) the frequency of v = 2 along the itineraries (must exceed 0 and stay below 0.585 for growth).
Usage: python3 collatz_shadow_error_20260927_itineraries.py [BMAX=40] [L=80]
"""
import math, sys
from fractions import Fraction

LOG23 = math.log2(3)


def itinerary(eta0, L):
    eta = eta0; w = []; etas = [eta]
    for _ in range(L):
        if eta < Fraction(5, 3):
            eta = (3 * eta - 1) / 2; w.append(1)
        else:
            eta = (3 * eta - 1) / 4; w.append(2)
        etas.append(eta)
    return w, etas


def analyse(w, eta0):
    L = len(w)
    d = [0]
    for v in w:
        d.append(d[-1] + v)
    partial = sum(Fraction(2 ** d[k], 3 ** (k + 1)) for k in range(L))
    growth = L * LOG23 - d[L]
    rhos = []
    S = 0
    for K in range(1, L + 1):
        S = 3 * S + 2 ** d[K - 1]
        mod = 2 ** (d[K] + 1)
        rho = ((2 ** d[K] - S) * pow(3, -K, mod)) % mod
        rhos.append(rho)
    stable = len(set(rhos[L // 2:])) == 1
    return growth, partial, rhos, stable, d[L]


def main():
    BMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    L = int(sys.argv[2]) if len(sys.argv) > 2 else 80
    samples = []
    for b in range(1, BMAX + 1):
        for a in range(b, 2 * b):
            q = Fraction(a, b)
            if q.denominator == b and 1 <= q < 2:
                samples.append(q)
    samples = sorted(set(samples))
    print("== %d rational starting errors in [1,2) with denominators <= %d, itineraries of length %d ==" % (len(samples), BMAX, L))
    n_growth = 0; n_stable = 0; worst_dev = 0.0; freqs = []
    worst_dev_at = None
    for q in samples:
        w, etas = itinerary(q, L)
        growth, partial, rhos, stable, dL = analyse(w, q)
        assert all(1 <= e < 2 for e in etas)
        if growth > 0:
            n_growth += 1
        dev = abs(float(partial - q))
        if dev > worst_dev:
            worst_dev = dev; worst_dev_at = q
        n_stable += stable
        freqs.append(w.count(2) / L)
    print(" all itineraries stay in [1,2): True; growth words over the window: %d of %d; max |partial Bernstein sum - eta_0| = %.3e (at eta_0 = %s; the tail is (2^(d_L)/3^L) eta_L < 2 * 2^(-0.335 L) ~ %.1e by the three-ones rule)" % (n_growth, len(samples), worst_dev, worst_dev_at, 2 * 2 ** (-0.335 * L)))
    print(" frequency of v = 2 along itineraries: min %.3f, mean %.3f, max %.3f (growth needs < 0.585)" % (min(freqs), sum(freqs) / len(freqs), max(freqs)))
    print(" itineraries whose exact 2-adic residue rho_K stabilises over K in [%d, %d] (would signal an integer point): %d of %d" % (L // 2, L, n_stable, len(samples)))
    # show two examples
    for q in (Fraction(1, 1), Fraction(3, 2), Fraction(211, 179), Fraction(7, 5)):
        w, etas = itinerary(q, 40)
        growth, partial, rhos, stable, dL = analyse(w, q)
        print(" eta_0 = %-8s itinerary %s ... growth %.2f, partial sum %.6f, rho_K for K = 36..40: %s" % (q, ''.join(str(v) for v in w[:24]), growth, float(partial), rhos[35:40]))
    # the periodic itineraries: eta_0 = R(w) reproduces w
    for w in ((1, 1, 1, 1, 2), (1, 1, 1, 1, 1, 2), (1, 1, 1, 1, 1, 2, 1, 1, 1, 2)):
        p = len(w); A = sum(w); head = Fraction(0); d = 0
        for t in range(p):
            head += Fraction(2 ** d, 3 ** (t + 1)); d += w[t]
        R = head / (1 - Fraction(2 ** A, 3 ** p))
        it, _ = itinerary(R, 3 * p)
        print(" periodic: R(%s) = %s, itinerary reproduces the word: %s" % (w, R, tuple(it[:p]) == w and it == list(w) * 3))


if __name__ == '__main__':
    main()
