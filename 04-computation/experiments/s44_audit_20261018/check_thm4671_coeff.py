#!/usr/bin/env python3
"""S44 audit: exact (over Q) value of T_1 on the Euler syzygy s_1 = x_0 e_1 - N x_1 e_0, i.e. the length-N part of
Phi(s_1) = L_-2 L_-1 v - N L_-3 v, written as A b_N + mu f_3 with b_N = x0^(N-2) x1 x2 + ((2N-7)/12) x0^(N-3) x1^3 and
f_3 = d^3 x_0^N.  Compares A with THM-4695/THM-4671's N(2N+1)(c - c_(2,2N+1))/8.  Own engine (vcore, singvec)."""
from fractions import Fraction
from vcore import apply_neg_vec, central_charge, dpoly
from singvec import singular_vector_exact


def main():
    rows = []
    for (p, q) in [(3, 4), (3, 5), (4, 5), (3, 7), (3, 8), (5, 6), (4, 7), (3, 10), (2, 7), (2, 9), (2, 13)]:
        N = (p - 1) * (q - 1) // 2
        c = central_charge(p, q)
        c0 = central_charge(2, 2 * N + 1)
        v = singular_vector_exact(p, q)
        u1 = apply_neg_vec(1, v)
        phi = apply_neg_vec(2, u1)
        for k, x in apply_neg_vec(3, v).items():
            phi[k] = phi.get(k, 0) - N * x
        phi = {k: x for k, x in phi.items() if x}
        assert all(len(k) <= N for k in phi), 'top length did not cancel'
        T = {k: x for k, x in phi.items() if len(k) == N}
        m3 = (2,) * (N - 1) + (5,)          # x0^(N-1) x3
        m12 = (2,) * (N - 2) + (3, 4)       # x0^(N-2) x1 x2
        m111 = (2,) * (N - 3) + (3, 3, 3)   # x0^(N-3) x1^3
        assert set(T) <= {m3, m12, m111}
        f3 = {(2,) * N: 1}
        for _ in range(3):
            f3 = dpoly(f3)
        mu = Fraction(T.get(m3, 0)) / f3[m3]
        A = Fraction(T.get(m12, 0)) - mu * f3[m12]
        resid = Fraction(T.get(m111, 0)) - mu * f3.get(m111, 0) - A * Fraction(2 * N - 7, 12)
        pred = Fraction(N * (2 * N + 1), 8) * (c - c0)
        print('(%d,%d) N=%d c=%s: T_1(s_1) = A b_N + mu f_3 with A = %s, mu = %s, residual %s; predicted A = %s; match %s'
              % (p, q, N, c, A, mu, resid, pred, A == pred and resid == 0))
        rows.append(A == pred and resid == 0)
    print('ALL MATCH:', all(rows))


if __name__ == '__main__':
    main()
