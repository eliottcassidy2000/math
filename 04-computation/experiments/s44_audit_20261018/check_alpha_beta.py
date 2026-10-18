#!/usr/bin/env python3
"""S44 audit, check A(b): length-(N-1) part of the Virasoro singular vector.
Claim (THM-4695(b)): v = L_-2^N + alpha L_-4 L_-2^(N-2) + beta L_-3^2 L_-2^(N-3) + (length <= N-2), with
    alpha = -2N(N-1)/3 - Nc/12,   beta = 11N(N-1)/48 + 5Nc/96.
Exact singular vectors over Q (independent engine, ascending PBW order; the length-(N-1) coefficients are
ordering-independent because reordering only produces lower-length terms).  Also checks that the length-(N-1)
part has no other monomials (there are none of weight 2N) and records the coefficient count."""
import sys, time, json
from fractions import Fraction
from vcore import central_charge
from singvec import singular_vector_exact, singular_vector_modp

P = 2147483647


def formulas(N, c):
    return (Fraction(-2 * N * (N - 1), 3) - N * c / 12, Fraction(11 * N * (N - 1), 48) + 5 * N * c / 96)


def main():
    models = [(3, 4), (2, 7), (3, 5), (2, 9), (3, 7), (4, 5), (2, 13), (3, 8), (5, 6), (3, 11), (2, 21), (4, 7), (3, 10),
              (2, 19), (5, 7), (4, 9)]
    res = []
    allok = True
    for (p, q) in models:
        t0 = time.time()
        N = (p - 1) * (q - 1) // 2
        c = central_charge(p, q)
        a_pred, b_pred = formulas(N, c)
        if 2 * N <= 24:
            v = singular_vector_exact(p, q)
            a = v.get((2,) * (N - 2) + (4,), Fraction(0))
            b = v.get((2,) * (N - 3) + (3, 3), Fraction(0))
            ok = (a == a_pred and b == b_pred)
            mode = 'exact Q'
        else:
            v = singular_vector_modp(p, q, P)
            red = lambda x: x.numerator % P * pow(x.denominator % P, -1, P) % P
            a = v.get((2,) * (N - 2) + (4,), 0)
            b = v.get((2,) * (N - 3) + (3, 3), 0)
            ok = (a == red(a_pred) and b == red(b_pred))
            mode = 'mod 2^31-1'
        lenNm1 = sorted(k for k in v if len(k) == N - 1)
        allok &= ok
        print('(%d,%d) N=%2d c=%-8s alpha=%-12s beta=%-12s formula=(%s, %s)  match=%s  len-(N-1) monomials=%s  [%s, %.1fs]'
              % (p, q, N, c, a, b, a_pred, b_pred, ok, lenNm1, mode, time.time() - t0), flush=True)
        res.append(dict(p=p, q=q, N=N, c=str(c), alpha=str(a), beta=str(b), alpha_pred=str(a_pred), beta_pred=str(b_pred),
                        match=ok, mode=mode))
    print('ALL MATCH:', allok)
    json.dump(res, open('check_alpha_beta.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
