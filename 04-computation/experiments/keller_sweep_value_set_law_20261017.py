#!/usr/bin/env python3
"""The tangent-sweep Keller maps (THM-4625; THM-1300 is d = 2) modulo p: image density and fibre-size law (S39).

THM-4625 (1): over F_3 = c != 0 and F_2 = b, the fibre of F over (a, b, c) is the set of simple roots w of
H_(b,c)(w) = c^2 a, where H_(b,c)(w) = q(w) + w (cb - p(w))/2, H' = (cb - p(w))/2 (critical points = punctures).
So over F_p the fibre size N(a, b, c) is the number of w in F_p with H(w) = c^2 a and p(w) != cb, and the image of
{F_3 != 0} is the union over (b, c) of the value sets of H_(b,c) on its non-critical points.  For a polynomial with
Galois group S_(d+1) the proportion of values with exactly k simple preimages tends to the proportion of permutations
of d+1 letters with k fixed points (Chebotarev for the polynomial's splitting field; for the value set this is the
Birch-Swinnerton-Dyer / S. D. Cohen theorem, CITED), so the image density is c_(d+1) = sum_(k=1)^(d+1) (-1)^(k-1)/k!
(2/3, 5/8, 19/30, 91/144, ... -> 1 - 1/e), and the fibre sizes avoid d (no permutation fixes exactly d of d+1 letters:
the finite-field shadow of THM-4625 (6)).
The d = 2 and d = 3 rows are the actual Keller maps (THM-1300, Gao's map); the d = 4, 5 polynomials p(w) fail Gao's side
conditions (S39 audit, MISTAKE-624), so those rows illustrate only the one-variable family H_(b,c), not Keller maps.
Usage: python keller_sweep_value_set_law_20261017.py
Session opus-2026-10-17-S39.
"""
from fractions import Fraction
from math import factorial
import numpy as np
from numba import njit


def poly_coeffs_q(pc):
    """q(w) = int_0^w (s/2) p'(s) ds for p = sum pc[i] w^i: coefficient of w^(i+1) is i pc[i] / (2(i+1))."""
    qc = [Fraction(0)] * (len(pc) + 1)
    for i, a in enumerate(pc):
        if i >= 1:
            qc[i + 1] += Fraction(i * a, 2 * (i + 1))
    return qc


@njit(cache=True)
def fibre_hist(p, hc_num, pc_mod, dmax):
    """hc_num[k][i]: coefficients of H's pieces mod p; returns histogram of N over (a, b, c), c != 0."""
    hist = np.zeros(dmax + 2, dtype=np.int64)
    cnt = np.zeros(p, dtype=np.int64)
    for b in range(p):
        for c in range(1, p):
            cb = (c * b) % p
            cnt[:] = 0
            for w in range(p):
                # p(w), q(w) mod p
                pw = 0; qw = 0; wp = 1
                for i in range(pc_mod.shape[0]):
                    pw = (pw + pc_mod[i] * wp) % p
                    qw = (qw + hc_num[i] * wp) % p
                    wp = (wp * w) % p
                qw = (qw + hc_num[pc_mod.shape[0]] * wp) % p
                if pw == cb:
                    continue
                inv2 = (p + 1) // 2
                H = (qw + w * ((cb - pw) % p) % p * inv2) % p
                cnt[H] += 1
            # value H = c^2 a ranges over all of F_p as a does (c != 0)
            for v in range(p):
                hist[cnt[v]] += 1
    return hist


def law(n):
    """proportion of permutations of n letters with exactly k fixed points, k = 0..n."""
    return [Fraction(1, factorial(k)) * sum(Fraction((-1) ** j, factorial(j)) for j in range(n - k + 1)) for k in range(n + 1)]


if __name__ == '__main__':
    cases = {2: [0, 4, -3], 3: [0, 6, -6, 1], 4: [0, 1, 3, -2, 1], 5: [0, -1, 1, -3, 1, 2]}   # p(w) coefficients, p(0) = 0
    for d, pc in cases.items():
        qc = poly_coeffs_q([Fraction(a) for a in pc])
        n = d + 1
        L = law(n)
        cdens = 1 - L[0]
        for p in [61, 101]:
            qmod = np.array([(x.numerator * pow(x.denominator, -1, p)) % p for x in qc], dtype=np.int64)
            pmod = np.array([a % p for a in pc], dtype=np.int64)
            hist = fibre_hist(p, qmod, pmod, n)
            tot = hist.sum()
            emp = [hist[k] / tot for k in range(n + 1)]
            print('d = %d, p = %3d: fibre-size frequencies %s' % (d, p, ' '.join('%d:%.4f' % (k, emp[k]) for k in range(n + 1))))
            print('              S_%d fixed-point law     %s' % (n, ' '.join('%d:%.4f' % (k, float(L[k])) for k in range(n + 1))))
            print('              image density %.4f vs c_%d = %s = %.4f; fibres of size %d: %d' % (1 - emp[0], n, cdens, float(cdens), d, hist[d]))
            assert hist[d] == 0
            assert abs((1 - emp[0]) - float(cdens)) < 0.02
    print('ALL CHECKS PASSED')
